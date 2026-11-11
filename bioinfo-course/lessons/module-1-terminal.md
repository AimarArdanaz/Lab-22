# Module 1: The Terminal

The Terminal is where most bioinformatics happens. In this module you learn to move around your files, inspect sequence files, chain commands together and write scripts that process hundreds of genomes at once.

[Back to the course map](../README.md)

- [Lesson 1: Navigating the file system](#lesson-1-navigating-the-file-system)
- [Lesson 2: Creating and reading files](#lesson-2-creating-and-reading-files)
- [Lesson 3: Copying, moving and deleting](#lesson-3-copying-moving-and-deleting)
- [Lesson 4: Redirection and pipes](#lesson-4-redirection-and-pipes)
- [Lesson 5: Your first script](#lesson-5-your-first-script)
- [Lesson 6: For loops](#lesson-6-for-loops)
- [Lesson 7: Conditionals and robust scripts](#lesson-7-conditionals-and-robust-scripts)
- [Lesson 8: awk and sed](#lesson-8-awk-and-sed)

---

## Lesson 1: Navigating the file system

**Goal:** know where you are, see what is there, and move between folders.

### Three basic commands

| Command | Stands for | What it does |
|---|---|---|
| `pwd` | print working directory | Shows the folder you are in |
| `ls` | list | Lists the contents of the folder |
| `cd` | change directory | Moves you to another folder |

Try them:

```bash
pwd
```

```bash
ls
```

```bash
cd Desktop
```

```bash
pwd
```

Go up one level:

```bash
cd ..
```

Go back to your home folder from anywhere:

```bash
cd ~
```

### Create folders

`mkdir` (make directory) creates a new folder:

```bash
cd ~/learn
```

```bash
mkdir data
```

### Useful tricks

- **Tab** autocompletes names. Type `cd ~/lea` and press Tab.
- **Up arrow** brings back previous commands.
- `clear` cleans the screen.
- **Avoid spaces and accents in file names.** Use `my_genomes`, not `my genomes`. Spaces break many commands.

### Exercise 1

Inside `~/learn`, create two more folders: `scripts` and `results`. Then run `ls` to check that `data`, `scripts` and `results` are there.

---

## Lesson 2: Creating and reading files

**Goal:** create a FASTA file and inspect it without opening any program.

### Create a FASTA file

Copy the whole block and paste it into the Terminal:

```bash
cat > data/example.fasta << 'EOF'
>strain01_bap
ATGAAAGCGTTACCGATTGCGCTGGCAGTGACCGGTATG
>strain02_esp
ATGAGTAACGTTCTGAAAGCCGGTCTGATCGCAGCGTTA
>strain03_bap
ATGAAAGCGTTACCGATTGCGCTGGCGGTGACCGGTATG
>strain04_control
ATGACCGATTTAGCGGCAAATCTGCGTGAAGCGCTGTAA
EOF
```

These sequences are invented, for practice only.

### How does `cat > file << 'EOF'` work?

`cat` does not really "show files". It **copies its input to its output**:

```
INPUT  ──►  cat  ──►  OUTPUT
```

- By default the input is a file and the output is the screen, so `cat file` displays it.
- Here, `<< 'EOF'` makes the input the text you type next, until a line that says `EOF`. This is called a **heredoc**.
- `> data/example.fasta` sends the output to a file instead of the screen.

So the text goes into `cat` and comes out into the file. Some details:

- `EOF` is just a label. `<< 'END'` followed by a final `END` works the same way.
- The closing `EOF` must be **alone on its line**.
- The quotes in `'EOF'` save the text exactly as written. Without them, symbols like `$` would be interpreted. Always use them.

Other ways to create files:

| Method | When to use it |
|---|---|
| `touch file.txt` | Create an empty file |
| `echo "text" > file.txt` | Write one line |
| `cat > file << 'EOF' ... EOF` | Paste several lines at once |
| `nano file.txt` | Edit by hand inside the Terminal (save: Ctrl+O, Enter; exit: Ctrl+X) |

### Look inside a file

| Command | What it does |
|---|---|
| `cat file` | Prints the whole file |
| `head -n 2 file` | Prints the first 2 lines |
| `tail -n 2 file` | Prints the last 2 lines |
| `wc -l file` | Counts lines |
| `less file` | Opens the file to scroll through it (quit with `q`) |

```bash
cat data/example.fasta
```

```bash
head -n 2 data/example.fasta
```

```bash
wc -l data/example.fasta
```

### `grep`: search inside files

`grep` prints the lines that contain a pattern. In a FASTA file, every header starts with `>`:

```bash
grep ">" data/example.fasta
```

Count the sequences with `-c` (count):

```bash
grep -c ">" data/example.fasta
```

Find the headers that mention *bap*:

```bash
grep "bap" data/example.fasta
```

Only use `grep -c ">"` on FASTA files. In other formats `>` means nothing special.

To see the folder in Finder:

```bash
open .
```

The `.` means "the folder I am in".

### Exercise 2

1. How many sequences are in `data/example.fasta`?
2. How many headers contain `bap`?
3. Run `grep -i "BAP" data/example.fasta`. What does `-i` do?

---

## Lesson 3: Copying, moving and deleting

**Goal:** organise files safely.

### `cp`: copy

```bash
cp data/example.fasta data/example_backup.fasta
```

The general form is `cp source destination`. To copy a whole folder, add `-r` (recursive):

```bash
cp -r data data_copy
```

### `mv`: move and rename

`mv` moves files. If the destination is a new name in the same folder, it renames:

```bash
mv data/example_backup.fasta scripts/
```

```bash
mv data_copy data_old
```

`mv` **overwrites the destination without asking** if it already exists. Add `-i` to be asked first:

```bash
mv -i old_name new_name
```

To rename a folder you are inside, step out of it first (`cd ..`).

### Paths

| Symbol | Meaning |
|---|---|
| `~` | Your home folder (`/Users/yourname`) |
| `.` | The current folder |
| `..` | The folder one level up |
| `/` at the start | An absolute path, from the root of the disk |

From `~/learn/scripts`, the path `../data/example.fasta` points to `~/learn/data/example.fasta`.

### Wildcards: many files at once

`touch` creates empty files:

```bash
touch data/strain01.fasta data/strain02.fasta data/strain03.fasta data/notes.txt
```

`*` means "anything":

```bash
ls data/*.fasta
```

```bash
ls data/strain*
```

### `rm`: delete (carefully)

**`rm` does not use the Bin. Deleted files are gone for good.**

```bash
rm -i data/notes.txt
```

```bash
rm -ri data_old
```

`-i` asks for confirmation; answer `y`. `-r` is needed for folders.

Golden rules:

1. **Before deleting with wildcards, run the same pattern with `ls`.** If the list is what you expect, replace `ls` with `rm`.
2. **Never run `rm -rf` with `~` or `/` unless you are completely sure.** `rm -rf ~` deletes your whole home folder without asking.
3. **Watch the spaces:** `rm -r data /strain*` (extra space) is very different from `rm -r data/strain*`.

### Exercise 3

1. Copy `data/example.fasta` into `results` with the name `bap_sequences.fasta`.
2. Delete `strain01.fasta`, `strain02.fasta` and `strain03.fasta` in a single command, checking with `ls` first.
3. Run `ls -R` to list every folder and its contents.

---

## Lesson 4: Redirection and pipes

**Goal:** save results to files and chain commands. This is where the Terminal becomes powerful.

### Redirection: `>` and `>>`

| Symbol | What it does |
|---|---|
| `>` | Saves the output to a file, **overwriting** it |
| `>>` | **Appends** the output to the end of a file |
| `<` | Uses a file as the input of a command |
| `<<` | Uses the text that follows as input (heredoc) |

Memory trick: the arrow points where the text goes. In `cat > file << 'EOF'`, `<<` points into `cat` (text goes in) and `>` points to the file (text comes out).

```bash
grep ">" data/example.fasta > results/headers.txt
```

```bash
echo "Analysis run today" > results/log.txt
```

```bash
echo "Second line" >> results/log.txt
```

```bash
cat results/log.txt
```

A common mistake is typing `>>` where `<<` was meant. `cat > file >> 'EOF'` does not start a heredoc: `cat` waits for keyboard input and the Terminal seems frozen. **If the Terminal ever stops returning the prompt, press Ctrl+C.**

### Pipes: `|`

`|` sends the output of one command into the next. On a Mac keyboard with a Spanish layout it is Option+1.

```bash
grep ">" data/example.fasta | wc -l
```

Read it as: "take the headers, **then** count them".

### Extract the *bap* sequences

`grep -A1` prints the matching line and the line after it (the sequence). `grep -v` prints lines that do **not** match; `^` means "start of line":

```bash
grep -A1 "bap" data/example.fasta | grep -v "^--" > results/bap_only.fasta
```

```bash
cat results/bap_only.fasta
```

This trick only works when every sequence is on a single line. Real FASTA files usually split sequences across lines; Lesson 10 shows the proper tool.

### Your first prevalence table

Create a strain table (made-up data):

```bash
cat > data/strains.csv << 'EOF'
strain,origin,bap,esp
C01,clinical,1,0
C02,clinical,1,1
C03,commensal,0,0
C04,clinical,0,1
C05,commensal,1,0
C06,clinical,1,1
C07,commensal,0,0
C08,clinical,1,0
C09,commensal,0,1
C10,clinical,1,1
EOF
```

New tools:

| Command | What it does |
|---|---|
| `tail -n +2` | Prints from line 2 onwards (drops the header) |
| `cut -d',' -f3` | Keeps column 3, using `,` as the separator |
| `sort` | Sorts lines |
| `uniq -c` | Counts consecutive identical lines (so sort first) |

How many strains carry *bap*?

```bash
tail -n +2 data/strains.csv | cut -d',' -f3 | sort | uniq -c
```

```
   4 0
   6 1
```

6 of 10 strains carry *bap*: a prevalence of 60%.

By origin:

```bash
tail -n +2 data/strains.csv | cut -d',' -f2,3 | sort | uniq -c
```

This is the basis of a contingency table, which you will test statistically in Lesson 16.

### Exercise 4

1. Compute the prevalence of *esp* (column 4) and save it to `results/esp_prevalence.txt`.
2. Compute *esp* by origin and **append** it to the same file.
3. How many strains carry **both** genes? Hint: `grep ",1,1$" data/strains.csv | wc -l`. `$` means "end of line".

---

## Lesson 5: Your first script

**Goal:** save commands in a file and run them all at once, on any input.

### Write the script

```bash
cat > scripts/prevalence.sh << 'EOF'
#!/bin/bash
# Computes the prevalence of bap and esp in a strain table
# Usage: bash scripts/prevalence.sh data/strains.csv

table=$1

echo "=== Analysis of: $table ==="

# Count strains (without the header). tr -d ' ' removes the spaces wc adds on macOS
total=$(tail -n +2 "$table" | wc -l | tr -d ' ')
echo "Number of strains: $total"

echo "--- bap (0 = absent, 1 = present) ---"
tail -n +2 "$table" | cut -d',' -f3 | sort | uniq -c

echo "--- esp ---"
tail -n +2 "$table" | cut -d',' -f4 | sort | uniq -c

echo "--- Strains with both bap and esp ---"
grep -c ",1,1$" "$table"
EOF
```

The quotes in `'EOF'` matter here: the script contains `$`, and without them the symbols would be interpreted while creating the file.

### The new pieces

| Piece | Meaning |
|---|---|
| `#!/bin/bash` | First line (*shebang*): which program runs the script |
| `# text` | A comment. Not executed. Use them a lot |
| `table=$1` | Creates a variable `table` holding the first argument |
| `$table` | Uses the value of the variable |
| `$( ... )` | Runs the command inside and captures its output |
| `"$table"` | Quotes protect paths that contain spaces |

**No spaces around `=`** when creating a variable: `table=$1` works, `table = $1` fails.

### Run it

```bash
bash scripts/prevalence.sh data/strains.csv
```

`data/strains.csv` is the **argument**; inside the script it becomes `$1`.

### Make it executable

```bash
chmod +x scripts/prevalence.sh
```

```bash
./scripts/prevalence.sh data/strains.csv
```

`ls -l scripts/` now shows an `x` (executable) in `-rwxr-xr-x`.

### Reuse it

Create a second table (made-up data):

```bash
cat > data/strains_hospital2.csv << 'EOF'
strain,origin,bap,esp
H01,clinical,1,1
H02,clinical,1,1
H03,clinical,0,1
H04,commensal,0,0
H05,clinical,1,0
EOF
```

```bash
./scripts/prevalence.sh data/strains_hospital2.csv
```

Save a report:

```bash
./scripts/prevalence.sh data/strains.csv > results/report_strains.txt
```

### Exercise 5

1. Edit the script with `nano scripts/prevalence.sh` and add a section that counts clinical strains (`grep -c "clinical" "$table"`).
2. Run it on both tables.
3. Run it **without** an argument: `./scripts/prevalence.sh`. What happens, and why?

---

## Lesson 6: For loops

**Goal:** repeat the same command on many files. One command, 300 genomes.

### The structure

```bash
for gene in bap esp agg; do
  echo "Searching for $gene"
done
```

- `gene` takes each value in turn.
- Everything between `do` and `done` runs once per value.
- You can type it on one line: `for gene in bap esp agg; do echo "$gene"; done`

Ranges:

```bash
for i in {1..5}; do echo "Sample $i"; done
```

### Loop over files

Create three small mock genome files with a loop. `printf` writes text; `\n` is a new line and `%s` is replaced by the value after the text:

```bash
for i in 1 2 3; do
  printf ">strain%s_contig1\nATGAAAGCGTTACCG\n>strain%s_contig2\nTTGACCGATTTAGCG\n" "$i" "$i" > data/genome_$i.fasta
done
```

```bash
ls data/genome_*.fasta
```

Count the contigs in each:

```bash
for f in data/genome_*.fasta; do
  n=$(grep -c ">" "$f")
  printf "%s\t%s\n" "$f" "$n"
done
```

`\t` is a tab, so the output is a two-column table.

### Name output files after input files

`basename` strips the folder and, optionally, the extension:

```bash
basename data/genome_1.fasta .fasta
```

Run the prevalence script on every table and save one report per table:

```bash
for f in data/strains*.csv; do
  name=$(basename "$f" .csv)
  echo "Processing $name"
  ./scripts/prevalence.sh "$f" > "results/${name}_report.txt"
done
```

`${name}_report` uses braces so the shell knows where the variable name ends. Without them, it would look for a variable called `name_report`.

### Read a list line by line

When you have a list of IDs (for example genome accessions), use `while read`:

```bash
printf "GCF_0001\nGCF_0002\nGCF_0003\n" > data/ids.txt
```

```bash
while read -r id; do
  echo "Would download $id"
done < data/ids.txt
```

### Pitfalls

- **Always quote variables:** `"$f"`, not `$f`. Unquoted names with spaces break into pieces.
- In zsh, a pattern with no matches (`for f in data/*.xyz`) stops with `no matches found`. In bash it loops once over the literal text `data/*.xyz`. Check with `ls` first.

### Exercise 6

Write a loop that saves a table `results/contig_counts.tsv` with two columns: the genome name (without folder or extension) and its number of contigs, for every `data/genome_*.fasta`.

---

## Lesson 7: Conditionals and robust scripts

**Goal:** make scripts check their inputs and stop on errors instead of producing wrong results silently.

### `if`

```bash
file="data/strains.csv"

if [ -f "$file" ]; then
  echo "$file exists"
else
  echo "$file not found"
fi
```

The spaces inside `[ ` and ` ]` are mandatory.

Common tests:

| Test | True when |
|---|---|
| `-f file` | The file exists |
| `-d folder` | The folder exists |
| `-s file` | The file exists and is not empty |
| `-z "$var"` | The variable is empty |
| `"$a" == "$b"` | Two texts are equal |
| `$n -eq 5` | Numbers are equal (also `-ne`, `-gt`, `-lt`, `-ge`, `-le`) |
| `! test` | Negates the test |

### Exit codes

Every command returns a number when it finishes: `0` means success, anything else means failure. `$?` holds the last one:

```bash
grep -q "bap" data/example.fasta
echo $?
```

```bash
grep -q "xyz" data/example.fasta
echo $?
```

`-q` makes `grep` silent. Exit codes let you chain commands:

```bash
mkdir -p results && echo "Folder ready"
```

`&&` runs the second command only if the first succeeded. `||` runs it only if the first failed.

### A robust version of the prevalence script

```bash
cat > scripts/prevalence.sh << 'EOF'
#!/bin/bash
# Computes the prevalence of bap and esp in a strain table
# Usage: ./scripts/prevalence.sh data/strains.csv

# Stop on errors (-e), on undefined variables (-u) and on failures inside pipes (pipefail)
set -euo pipefail

# Check that exactly one argument was given
if [ $# -ne 1 ]; then
  echo "Usage: $0 <table.csv>" >&2
  exit 1
fi

table=$1

if [ ! -f "$table" ]; then
  echo "Error: file '$table' not found" >&2
  exit 1
fi

echo "=== Analysis of: $table ==="

total=$(tail -n +2 "$table" | wc -l | tr -d ' ')
echo "Number of strains: $total"

echo "--- bap (0 = absent, 1 = present) ---"
tail -n +2 "$table" | cut -d',' -f3 | sort | uniq -c

echo "--- esp ---"
tail -n +2 "$table" | cut -d',' -f4 | sort | uniq -c

echo "--- Strains with both bap and esp ---"
# grep exits with 1 when it finds nothing; '|| true' stops set -e from killing the script
grep -c ",1,1$" "$table" || true
EOF
```

New pieces:

| Piece | Meaning |
|---|---|
| `set -euo pipefail` | Safety switch. Put it at the top of every script |
| `$#` | Number of arguments given |
| `$0` | The name of the script |
| `>&2` | Sends the message to the error channel (stderr) instead of normal output |
| `exit 1` | Stops the script with a failure code |
| `|| true` | Accepts a "failure" that is actually fine |

The last point is a classic trap: with `set -e`, a `grep` that finds no matches counts as an error and stops the script.

Test it:

```bash
./scripts/prevalence.sh
```

```bash
./scripts/prevalence.sh data/does_not_exist.csv
```

```bash
./scripts/prevalence.sh data/strains.csv
```

### Exercise 7

Add a check that warns (without stopping) when the table has fewer than 5 strains: "Warning: fewer than 5 strains, prevalence estimates will be very imprecise". Test it on `data/strains_hospital2.csv` and on a new table with 3 strains.

---

## Lesson 8: awk and sed

**Goal:** process tables and FASTA files in one line. `awk` works on columns; `sed` edits text.

### awk: a calculator for columns

`awk` reads a file line by line and splits each line into fields: `$1`, `$2`... `$0` is the whole line. `-F','` sets the separator.

Print columns 1 and 3:

```bash
awk -F',' '{print $1, $3}' data/strains.csv
```

Built-in variables:

| Variable | Meaning |
|---|---|
| `NR` | Number of the current line (record) |
| `NF` | Number of fields in the current line |
| `FS`, `OFS` | Input and output field separators |

Skip the header and keep only *bap*-positive strains:

```bash
awk -F',' 'NR > 1 && $3 == 1' data/strains.csv
```

Compute the prevalence of *bap*. The `END` block runs after the last line:

```bash
awk -F',' 'NR > 1 {n++; pos += $3} END {printf "bap: %d/%d (%.1f%%)\n", pos, n, 100 * pos / n}' data/strains.csv
```

`printf` formats: `%d` integer, `%.1f` one decimal, `%%` a literal `%`.

Prevalence by origin, using arrays indexed by the origin name:

```bash
awk -F',' 'NR > 1 {n[$2]++; pos[$2] += $3}
  END {for (o in n) printf "%s\t%d/%d\t%.1f%%\n", o, pos[o], n[o], 100 * pos[o] / n[o]}' data/strains.csv
```

Prevalence of **every** gene column, whatever the number of genes. The gene names are read from the header:

```bash
awk -F',' 'NR == 1 {ncol = NF; for (i = 3; i <= NF; i++) gene[i] = $i; next}
  {n++; for (i = 3; i <= ncol; i++) pos[i] += $i}
  END {for (i = 3; i <= ncol; i++) printf "%s\t%d/%d\t%.1f%%\n", gene[i], pos[i], n, 100 * pos[i] / n}' data/strains.csv
```

`next` skips to the next line, so the header is not counted.

Sequence lengths from a FASTA file (works with multi-line sequences too):

```bash
awk '/^>/ {if (name) print name "\t" len; name = substr($0, 2); len = 0; next}
  {len += length($0)}
  END {print name "\t" len}' data/example.fasta
```

### sed: find and replace

```bash
sed 's/clinical/hospital/' data/strains.csv
```

`s/old/new/` replaces the first match on each line; add `g` to replace all: `s/old/new/g`. `sed` prints the result and leaves the file untouched.

Add a prefix to every FASTA header (useful before merging genomes, so you know where each contig came from):

```bash
sed 's/^>/>hospital2_/' data/example.fasta
```

Other uses:

```bash
sed '1d' data/strains.csv
```

```bash
sed -n '2,4p' data/strains.csv
```

`1d` deletes line 1; `-n '2,4p'` prints only lines 2 to 4.

**Editing a file in place differs between macOS and Linux:**

| System | Command |
|---|---|
| macOS | `sed -i '' 's/old/new/' file` |
| Linux | `sed -i 's/old/new/' file` |

Keep a backup or work on copies before editing in place.

**Another macOS difference:** `sed` on a Mac does not understand `\t` in the replacement. To convert a CSV into a tab-separated file, use `tr`:

```bash
tr ',' '\t' < data/strains.csv > results/strains.tsv
```

### Exercise 8

1. With `awk`, compute the prevalence of *esp* by origin.
2. With `sed`, print `data/example.fasta` with the prefix `cohortA_` added to every header and save it to `results/example_renamed.fasta`.
3. With `awk`, list the strains that carry *esp* but **not** *bap*.

Next: [Module 2, bioinformatics tools](module-2-tools.md)
