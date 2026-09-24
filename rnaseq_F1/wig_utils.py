#!/usr/bin/env python3
"""
wig_utils.py
====================================================================
Funciones compartidas por pipeline_single_end.py y pipeline_paired_end.py
para separar un BAM ordenado en hebras (+ / -) y convertir la cobertura
resultante a formato WIG (variableStep).

Este archivo debe estar en la MISMA carpeta que esos dos scripts,
ya que ambos hacen `from wig_utils import ...`.

No añade dependencias nuevas: todo se hace con `samtools`, que ya
tenías instalado para el resto del pipeline.

--------------------------------------------------------------------
¿QUÉ ES "STRANDED COVERAGE" Y POR QUÉ SE HACE ASÍ AQUÍ?
--------------------------------------------------------------------
Cuando secuencias RNA-seq con un kit "stranded", cada read conserva
la información de qué hebra del ADN se transcribió originalmente.
Bowtie2 marca en el FLAG del SAM/BAM si un read mapeó en la hebra +
o en la hebra - (bit 0x10, valor 16 del flag).

PERO el bit 16 del flag NO es directamente "la hebra transcrita" en
todos los protocolos:

- Kits tipo dUTP / Illumina TruSeq Stranded mRNA (LOS MÁS HABITUALES
  hoy en día comercialmente): el read que secuencias es
  COMPLEMENTARIO a la hebra molde original -> hay que INVERTIR la
  hebra de mapeo para saber la hebra transcrita real.
  Esto es lo que aquí llamamos 'fr-firststrand'.

- Kits antiguos tipo "Ligation" (poco habituales ya): el read SÍ
  coincide con la hebra transcrita tal cual mapea.
  Esto es lo que aquí llamamos 'fr-secondstrand'.

- Si tu librería no es "stranded" (protocolo genómico estándar, o
  RNA-seq no direccional), no hay hebra que recuperar: usa
  'unstranded' y se genera un único WIG con la cobertura total.

SI NO SABES CUÁL DE LOS DOS USASTE (dUTP vs Ligation), la forma más
fiable de comprobarlo es correr UNA VEZ la herramienta
`infer_experiment.py` del paquete RSeQC sobre uno de tus BAM ya
alineados, usando una anotación de tu genoma en formato BED:

    pip install RSeQC
    infer_experiment.py -i SRRxxxx_sorted.bam -r anotacion.bed

Te devuelve el % de reads consistentes con cada convención: el que
tenga el porcentaje alto (normalmente >80-90%) es tu protocolo real.
Para bacterias, puedes generar el BED a partir del GFF de Prokka/NCBI
(hay conversores como `gff2bed` de BEDOPS).

¿Por qué te importa esto para tu tesis? Porque en genomas bacterianos
los genes están muy empaquetados en ambas hebras, con solapamientos
y antisentido frecuentes -- justo el caso de islas de patogenicidad
como la PAI de E. faecalis, donde está esp. Sin separar por hebra no
puedes saber si la señal de cobertura pertenece al gen de interés o
a algo solapante en la hebra opuesta.
--------------------------------------------------------------------
"""

import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Utilidades de bajo nivel sobre samtools
# ---------------------------------------------------------------------------

def _samtools_view_filter(input_bam, flag_f, flag_F, output_bam, threads):
    """
    Escribe en output_bam los reads de input_bam que cumplan:
      - flag_f: bits que DEBEN estar presentes (samtools -f), o None
      - flag_F: bits que NO deben estar presentes (samtools -F), o None
    """
    cmd = ["samtools", "view", "-@", str(threads), "-b"]
    if flag_f is not None:
        cmd += ["-f", str(flag_f)]
    if flag_F is not None:
        cmd += ["-F", str(flag_F)]
    cmd += [str(input_bam)]
    with open(output_bam, "wb") as out:
        subprocess.run(cmd, stdout=out, check=True)


def _merge_bams(bam_list, output_bam, threads):
    """
    Fusiona 2+ BAMs ya ordenados por coordenadas en uno solo,
    conservando el orden (todos proceden del mismo BAM ordenado
    original, así que el resultado sigue ordenado sin re-sortear).
    """
    cmd = ["samtools", "merge", "-f", "-@", str(threads), str(output_bam)]
    cmd += [str(b) for b in bam_list]
    subprocess.run(cmd, check=True)


def compute_scale_factor(bam_path, threads=4):
    """
    Devuelve el factor de escala para normalizar a CPM
    (Counts Per Million reads mapeados): 1e6 / total_reads_mapeados.
    Se calcula sobre el BAM completo (no por hebra), para que las
    hebras + y - de una misma muestra queden en la misma escala.
    """
    result = subprocess.run(
        ["samtools", "view", "-@", str(threads), "-c", "-F", "4", str(bam_path)],
        stdout=subprocess.PIPE, text=True, check=True,
    )
    total_mapped = int(result.stdout.strip() or 0)
    if total_mapped == 0:
        return 1.0
    return 1_000_000 / total_mapped


def bam_to_wig(bam_path, wig_path, threads=4, scale_factor=1.0, track_name=None):
    """
    Convierte un BAM (ya filtrado/ordenado) a un WIG (variableStep,
    span=1) con la profundidad de cobertura por posición.

    Solo se escriben posiciones con profundidad > 0 (comportamiento
    por defecto de 'samtools depth'), lo cual es válido en formato
    WIG: las posiciones no listadas se interpretan como sin datos.

    Si scale_factor != 1.0, los valores se normalizan (p.ej. a CPM)
    en vez de escribir la profundidad cruda.
    """
    name = track_name or Path(wig_path).stem
    proc = subprocess.Popen(
        ["samtools", "depth", "-@", str(threads), str(bam_path)],
        stdout=subprocess.PIPE,
        text=True,
    )

    current_chrom = None
    with open(wig_path, "w") as wig:
        wig.write(f'track type=wiggle_0 name="{name}"\n')
        for line in proc.stdout:
            chrom, pos, depth = line.rstrip("\n").split("\t")
            if chrom != current_chrom:
                wig.write(f"variableStep chrom={chrom} span=1\n")
                current_chrom = chrom
            if scale_factor != 1.0:
                wig.write(f"{pos} {float(depth) * scale_factor:.3f}\n")
            else:
                wig.write(f"{pos} {depth}\n")

    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools depth falló para {bam_path}")


# ---------------------------------------------------------------------------
# Separación por hebra: SINGLE-END
# ---------------------------------------------------------------------------

def split_strands_single_end(sorted_bam, prefix, threads, library_strand,
                              normalize_cpm=False):
    """
    Genera el/los WIG de cobertura a partir de un BAM single-end
    ya ordenado e indexado.

    Para un read individual, la hebra de mapeo se lee directamente
    del flag 16 (0x10 = alineado en hebra reversa):
      - flag 16 NO presente -> read mapeado en hebra +
      - flag 16 presente    -> read mapeado en hebra -

    library_strand:
      'fr-firststrand'  -> kit dUTP / TruSeq Stranded (el más común):
                            la hebra de mapeo está INVERTIDA respecto
                            a la hebra transcrita real.
      'fr-secondstrand' -> kit tipo Ligation: la hebra de mapeo
                            coincide con la hebra transcrita.
      'unstranded'      -> no se separa; un único WIG con cobertura total.

    Devuelve la lista de Path a los WIG generados.
    """
    sorted_bam = Path(sorted_bam)
    outdir = sorted_bam.parent
    scale = compute_scale_factor(sorted_bam, threads) if normalize_cpm else 1.0

    if library_strand == "unstranded":
        wig = outdir / f"{prefix}.wig"
        bam_to_wig(sorted_bam, wig, threads, scale, track_name=prefix)
        return [wig]

    fwd_bam = outdir / f"{prefix}_fwdmapped.tmp.bam"  # reads mapeados en hebra +
    rev_bam = outdir / f"{prefix}_revmapped.tmp.bam"  # reads mapeados en hebra -

    # -F 20 = excluye unmapped(4) y reverse(16)  -> deja solo mapeados en +
    _samtools_view_filter(sorted_bam, None, 20, fwd_bam, threads)
    # -f 16 = exige reverse(16)                  -> deja solo mapeados en -
    _samtools_view_filter(sorted_bam, 16, None, rev_bam, threads)

    if library_strand == "fr-firststrand":
        # el mapeo está invertido respecto a la hebra transcrita real
        plus_source, minus_source = rev_bam, fwd_bam
    else:  # fr-secondstrand
        plus_source, minus_source = fwd_bam, rev_bam

    plus_wig = outdir / f"{prefix}_plus.wig"
    minus_wig = outdir / f"{prefix}_minus.wig"
    bam_to_wig(plus_source, plus_wig, threads, scale, track_name=f"{prefix}_plus")
    bam_to_wig(minus_source, minus_wig, threads, scale, track_name=f"{prefix}_minus")

    fwd_bam.unlink(missing_ok=True)
    rev_bam.unlink(missing_ok=True)

    return [plus_wig, minus_wig]


# ---------------------------------------------------------------------------
# Separación por hebra: PAIRED-END
# ---------------------------------------------------------------------------

def split_strands_paired_end(sorted_bam, prefix, threads, library_strand,
                              normalize_cpm=False):
    """
    Igual que split_strands_single_end pero para paired-end.

    En paired-end la hebra de origen del FRAGMENTO no se lee de un
    solo read: hay que combinar el flag de los dos mates (esto es lo
    que en herramientas como RSeQC se conoce como la convención
    '1+-,1-+,2++,2--' / '1++,1--,2+-,2-+'):

      reads_A = (mate2 mapeado en +) O (mate1 mapeado en -)
      reads_B = (mate2 mapeado en -) O (mate1 mapeado en +)

    Con 'fr-firststrand' (dUTP, el caso más común), reads_A
    corresponde a la hebra + transcrita real; con 'fr-secondstrand'
    es al revés.

    Flags SAM usados:
      64  = mate1 (first in pair)
      128 = mate2 (second in pair)
      16  = alineado en hebra reversa
    """
    sorted_bam = Path(sorted_bam)
    outdir = sorted_bam.parent
    scale = compute_scale_factor(sorted_bam, threads) if normalize_cpm else 1.0

    if library_strand == "unstranded":
        wig = outdir / f"{prefix}.wig"
        bam_to_wig(sorted_bam, wig, threads, scale, track_name=prefix)
        return [wig]

    a1 = outdir / f"{prefix}_a1.tmp.bam"  # mate2 forward  (-f 128 -F 16)
    a2 = outdir / f"{prefix}_a2.tmp.bam"  # mate1 reverse  (-f 80,  80=64+16)
    b1 = outdir / f"{prefix}_b1.tmp.bam"  # mate2 reverse  (-f 144, 144=128+16)
    b2 = outdir / f"{prefix}_b2.tmp.bam"  # mate1 forward  (-f 64  -F 16)

    _samtools_view_filter(sorted_bam, 128, 16, a1, threads)
    _samtools_view_filter(sorted_bam, 80, None, a2, threads)
    _samtools_view_filter(sorted_bam, 144, None, b1, threads)
    _samtools_view_filter(sorted_bam, 64, 16, b2, threads)

    reads_A = outdir / f"{prefix}_readsA.tmp.bam"
    reads_B = outdir / f"{prefix}_readsB.tmp.bam"
    _merge_bams([a1, a2], reads_A, threads)
    _merge_bams([b1, b2], reads_B, threads)

    if library_strand == "fr-firststrand":
        plus_source, minus_source = reads_A, reads_B
    else:  # fr-secondstrand
        plus_source, minus_source = reads_B, reads_A

    plus_wig = outdir / f"{prefix}_plus.wig"
    minus_wig = outdir / f"{prefix}_minus.wig"
    bam_to_wig(plus_source, plus_wig, threads, scale, track_name=f"{prefix}_plus")
    bam_to_wig(minus_source, minus_wig, threads, scale, track_name=f"{prefix}_minus")

    for tmp in (a1, a2, b1, b2, reads_A, reads_B):
        tmp.unlink(missing_ok=True)

    return [plus_wig, minus_wig]
