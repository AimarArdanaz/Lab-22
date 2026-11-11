# Module 0: Prepare your Mac

Before writing any code you need three things: the Terminal, a package manager for scientific software (conda), and a good text editor.

[Back to the course map](../README.md)

---

## 1. Open the Terminal

- Press **Cmd + Space**, type **Terminal** and press Enter.
- You will see a line like this:

  ```
  yourname@MacBook ~ %
  ```

  This is the **prompt**. The Terminal is waiting for a command. The `~` means you are in your home folder. The `%` means your shell is **zsh**, the default on macOS. zsh is almost identical to **bash**, the shell used in most bioinformatics tutorials. Everything in this course works in both.

Tip: drag the Terminal icon to your Dock. You will use it every day.

## 2. Install Apple's command line tools

They include `git` and the compilers some tools need:

```bash
xcode-select --install
```

A window will open. Accept and wait for it to finish. If it says they are already installed, you are done.

## 3. Install conda (Miniforge)

**conda** installs bioinformatics programs (BLAST, seqkit, abricate...) and keeps each project's software in a separate **environment**, so versions never clash. We use **Miniforge**, a lightweight conda distribution configured to use the community `conda-forge` channel.

> If you already have conda (you see `(base)` or another name in brackets at the start of your prompt), skip to step 4.

Check which processor your Mac has:

```bash
uname -m
```

`arm64` means Apple Silicon (M1, M2, M3...). `x86_64` means Intel. The installer below picks the right version automatically.

Download and run the installer:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

```bash
bash Miniforge3-$(uname)-$(uname -m).sh
```

- Press Enter to read the licence, type `yes` to accept, accept the default location.
- When asked whether to initialise conda, answer `yes`.
- **Close the Terminal and open it again.** You should now see `(base)` at the start of the prompt.

## 4. Configure the bioinformatics channels

**Bioconda** is the channel where almost all bioinformatics tools live. Run these three commands once:

```bash
conda config --add channels bioconda
```

```bash
conda config --add channels conda-forge
```

```bash
conda config --set channel_priority strict
```

## 5. Install a code editor

Install **Visual Studio Code** from https://code.visualstudio.com. It highlights code in bash, R and Python and has a built-in terminal.

To open folders from the Terminal with `code .`, open VS Code, press **Cmd + Shift + P**, type `shell command` and choose **Install 'code' command in PATH**.

Useful extensions (search for them in the Extensions panel): *Python*, *R*, *Rainbow CSV*.

## 6. Create your working folder

All lessons happen inside `~/learn`:

```bash
mkdir ~/learn
```

```bash
cd ~/learn
```

## Checklist

- [ ] The Terminal opens and shows a prompt.
- [ ] `git --version` prints a version number.
- [ ] `conda --version` prints a version number.
- [ ] VS Code is installed.
- [ ] The folder `~/learn` exists.

Next: [Module 1, the Terminal](module-1-terminal.md)
