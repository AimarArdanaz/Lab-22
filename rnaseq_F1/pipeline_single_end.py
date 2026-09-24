#!/usr/bin/env python3
"""
Pipeline SRA -> FASTQ -> Alineamiento (SINGLE-END) -> BAM ordenado -> WIG
====================================================================
Uso:
    python pipeline_single_end.py lista_srr.txt

Requisitos (deben estar en el PATH):
    - sra-tools (prefetch, fasterq-dump) >= 3.x
    - bowtie2
    - samtools

IMPORTANTE: este script importa wig_utils.py, que debe estar en la
misma carpeta (contiene las funciones de separación por hebra y
generación de WIG; no añade dependencias nuevas, todo usa samtools).

El archivo lista_srr.txt debe tener un accession SRR por línea, ej:
    SRR13200561
    SRR13200562
    SRR13200563
"""

import subprocess
import sys
import shutil
from pathlib import Path
from datetime import datetime

from wig_utils import split_strands_single_end

# ---------------------- CONFIGURACIÓN (editar aquí) ------------------------
REF_INDEX = "../F1_index"   # ruta al índice de bowtie2 (sin extensión .bt2)
THREADS = 4                 # nº de hilos para fasterq-dump / bowtie2 / samtools
OUTDIR = Path(".")          # carpeta de salida de los FASTQ/SAM/BAM/WIG
KEEP_FASTQ = False          # True = conservar los .fastq tras alinear
KEEP_SRA = False            # True = conservar la carpeta con el .sra tras extraer

# --- Config de cobertura por hebra (stranded coverage) ----------------------
# 'fr-firststrand'  -> kit dUTP / Illumina TruSeq Stranded (EL MÁS HABITUAL
#                       comercialmente): el read mapea en la hebra OPUESTA
#                       a la hebra realmente transcrita.
# 'fr-secondstrand' -> kit tipo Ligation (menos habitual hoy): el read
#                       mapea en la MISMA hebra que la transcrita.
# 'unstranded'      -> no separar por hebra; un único WIG con cobertura total.
#
# Si no sabes cuál usaste, revisa la ficha del kit de tu proveedor de
# secuenciación, o comprueba empíricamente con RSeQC (ver wig_utils.py,
# función infer_experiment.py sobre uno de tus BAM ya generados).
LIBRARY_STRAND = "fr-firststrand"

# True = normaliza la cobertura a CPM (counts per million reads mapeados)
# en vez de profundidad cruda. Útil si vas a comparar el WIG de una
# muestra con el de otra con distinta profundidad de secuenciación.
NORMALIZE_CPM = False
# -----------------------------------------------------------------------------


def log(msg, logfile):
    """Imprime en pantalla y guarda en el archivo de log."""
    print(msg)
    with open(logfile, "a") as f:
        f.write(msg + "\n")


def run(cmd, logfile, capture_to=None):
    """
    Ejecuta un comando externo. Si capture_to se indica, redirige stdout+stderr
    a ese archivo (equivalente a '> archivo 2>&1' o '2> archivo').
    Devuelve el código de salida (returncode).
    """
    log(f"$ {' '.join(cmd)}", logfile)
    if capture_to:
        with open(capture_to, "w") as out:
            result = subprocess.run(cmd, stdout=out, stderr=out)
    else:
        with open(logfile, "a") as out:
            result = subprocess.run(cmd, stdout=out, stderr=out)
    return result.returncode


def main():
    if len(sys.argv) != 2:
        print(f"Uso: python {sys.argv[0]} <archivo_lista_srr.txt>")
        sys.exit(1)

    lista_path = Path(sys.argv[1])
    if not lista_path.is_file():
        print(f"ERROR: no se encuentra el archivo '{lista_path}'")
        sys.exit(1)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = OUTDIR / f"pipeline_log_{timestamp}.txt"

    log(f"=== Pipeline single-end iniciado: {datetime.now()} ===", logfile)
    log(f"    Modo de hebra (LIBRARY_STRAND): {LIBRARY_STRAND}", logfile)
    log(f"    Normalización CPM: {NORMALIZE_CPM}", logfile)

    # Leemos la lista, quitando \r (por si el .txt viene de Windows) y líneas vacías
    with open(lista_path, "r") as f:
        sra_list = [line.strip().replace("\r", "") for line in f if line.strip()]

    total = len(sra_list)
    ok = 0
    fail = 0

    for i, srr in enumerate(sra_list, start=1):
        log("", logfile)
        log("=" * 60, logfile)
        log(f">>> [{i}/{total}] Procesando {srr}", logfile)
        log("=" * 60, logfile)

        # ---------- 1. Descarga (prefetch) ----------
        run(["prefetch", srr, "-O", str(OUTDIR)], logfile)
        sra_file = OUTDIR / srr / f"{srr}.sra"
        if not sra_file.is_file():
            log(f"!!! ERROR: prefetch falló para {srr}, se omite.", logfile)
            fail += 1
            continue

        # ---------- 2. SRA -> FASTQ ----------
        run(
            ["fasterq-dump", str(sra_file), "-e", str(THREADS), "-O", str(OUTDIR)],
            logfile,
        )
        fastq = OUTDIR / f"{srr}.fastq"
        if not fastq.is_file():
            log(
                f"!!! ERROR: fasterq-dump no generó {fastq} "
                f"(¿es realmente single-end?), se omite {srr}.",
                logfile,
            )
            fail += 1
            continue

        # ---------- 3. Alineamiento con bowtie2 (single-end) ----------
        sam = OUTDIR / f"{srr}.sam"
        summary = OUTDIR / f"{srr}_bowtie2_summary.txt"

        rc = run(
            [
                "bowtie2", "-p", str(THREADS), "-x", REF_INDEX,
                "-U", str(fastq),
                "-S", str(sam),
            ],
            logfile,
            capture_to=summary,
        )

        if not sam.is_file() or sam.stat().st_size == 0:
            log(f"!!! ERROR: bowtie2 no generó SAM para {srr}, se omite.", logfile)
            fail += 1
            continue

        log(f"Resumen bowtie2 ({srr}):", logfile)
        with open(summary) as f:
            log(f.read(), logfile)

        # ---------- 4. SAM -> BAM ordenado e indexado ----------
        bam = OUTDIR / f"{srr}_sorted.bam"

        # samtools view -F 4 (quita no-mapeados) | samtools sort
        view_proc = subprocess.Popen(
            ["samtools", "view", "-@", str(THREADS), "-b", "-F", "4", str(sam)],
            stdout=subprocess.PIPE,
        )
        with open(logfile, "a"):
            subprocess.run(
                ["samtools", "sort", "-@", str(THREADS), "-o", str(bam)],
                stdin=view_proc.stdout,
            )
        view_proc.stdout.close()
        view_proc.wait()

        run(["samtools", "index", str(bam)], logfile)

        # ---------- 5. BAM -> WIG (cobertura, separada por hebra) ----------
        try:
            wig_files = split_strands_single_end(
                bam, srr, THREADS, LIBRARY_STRAND, normalize_cpm=NORMALIZE_CPM
            )
            for w in wig_files:
                log(f"    WIG generado -> {w}", logfile)
        except Exception as e:
            log(f"!!! ERROR generando WIG para {srr}: {e}", logfile)
            # No se marca como fallo total: el BAM sí se generó correctamente.

        # ---------- 6. Limpieza ----------
        sam.unlink(missing_ok=True)

        if not KEEP_FASTQ:
            fastq.unlink(missing_ok=True)

        if not KEEP_SRA:
            shutil.rmtree(OUTDIR / srr, ignore_errors=True)

        log(f">>> {srr} completado -> {bam}", logfile)
        ok += 1

    log("", logfile)
    log(f"=== Pipeline finalizado: {datetime.now()} ===", logfile)
    log(f"Completados correctamente: {ok} / {total}", logfile)
    log(f"Fallidos: {fail} / {total}", logfile)


if __name__ == "__main__":
    main()
