# Strand-Passage Explorer

Strand-Passage Explorer is a research tool for studying **strand passages**
(crossing changes) on knots and links given by a signed
**Dowker-Thistlethwaite (DT) code**.

The main feature is component-colour preservation: the identity of each original
link component is tracked through chained passages and simplifications, so the
same component can be followed across a branching exploration.

## What It Does

- **Interactive GUI:** load a DT code, click a crossing, and open a new diagram
  window for the after-passage result. Each result window is clickable, so
  passages can branch and chain.
- **Batch mode (`--nongui`):** run all one- and two-step strand passages, writing
  an Excel spreadsheet and an overview SVG.
- **Demo mode (`--demo`):** render a headless before/after cascade image.
- **SnapPy/Sage path:** use SnapPy simplification, Jones polynomials, linking
  invariants, and component-colour recovery.
- **Continuity fallback:** the scripts still run without Sage/SnapPy, but Jones
  polynomials and SnapPy-based colour matching are unavailable.

Development notes and version history are in [DEVELOPMENT_LOG.md](DEVELOPMENT_LOG.md).

## Repository Layout

```text
strand_passage_guiV3_4.py        Main entry point: GUI, --nongui, and --demo
link_engine_v3_4.py              Diagram engine and SnapPy bridge
draw_dt_original_labelsV3_11.py  DT parser, layout, renderer, and standalone GUI
check_two_dt.py                  Standalone SnapPy/Sage DT-comparison utility
assets/strand_passage_icon.png   Optional window/task-menu icon
bin/strand-passage               Convenience launcher
requirements.txt                 Python package notes
DEVELOPMENT_LOG.md               Development notes and version history
LICENSE                          MIT license
```

Import chain:

```text
strand_passage_guiV3_4 -> link_engine_v3_4 -> draw_dt_original_labelsV3_11
```

## Install

Clone the repository:

```bash
git clone https://github.com/DiLiuLab/dt_strand_passage_explorer.git
cd dt_strand_passage_explorer
```

`git clone` downloads the project into a new folder and sets `origin` to the
GitHub repository, so future updates know where to come from.

Install the regular Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

For full functionality, run with **Sage + SnapPy**:

```bash
sage -python strand_passage_guiV3_4.py --help
```

SnapPy is intentionally not pinned in `requirements.txt`, because this project
is intended to use the SnapPy/Sage installation on the research machine.

## Run

Interactive GUI:

```bash
sage -python strand_passage_guiV3_4.py
sage -python strand_passage_guiV3_4.py --dt "DT: [(4,6,2)]"
```

If the TkAgg backend is not available:

```bash
python3 strand_passage_guiV3_4.py --gui-backend agg
```

Batch spreadsheet and overview SVG:

```bash
sage -python strand_passage_guiV3_4.py --nongui \
  --dt "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]" \
  --out strand_passage_results.xlsx
```

In `--nongui` mode, `--out` controls both output files. If you set:

```text
--out /path/to/name.xlsx
```

the script writes:

```text
/path/to/name.xlsx
/path/to/name_overview.svg
```

If the path does not end in `.xlsx`, the script adds it first. For example,
`--out results/run1` writes `results/run1.xlsx` and
`results/run1_overview.svg`. To choose the SVG directory and basename, choose the
directory and basename of the spreadsheet path.

Headless cascade figure:

```bash
python3 strand_passage_guiV3_4.py --dt "DT: [(4,6,2)]" --demo 2 1 --out chain.png
```

Standalone DT comparison utility:

```bash
sage -python check_two_dt.py
```

## Pull Updates

To update an existing clone:

```bash
git pull
```

`git pull` asks GitHub for new commits on the current branch and applies them to
your local copy. If you have local edits, commit or save them first so Git can
merge cleanly.

To see what changed after pulling:

```bash
git log --oneline -5
```

## Executable Use

The repository includes a launcher:

```bash
./bin/strand-passage --help
./bin/strand-passage --dt "DT: [(4,6,2)]"
```

The launcher uses `sage -python` when Sage is available, and falls back to
`python3` otherwise.

To make the Python scripts directly executable on macOS/Linux:

```bash
chmod +x strand_passage_guiV3_4.py
chmod +x draw_dt_original_labelsV3_11.py
chmod +x check_two_dt.py
```

Then they can be run as:

```bash
./strand_passage_guiV3_4.py --gui-backend agg
```

For SnapPy/Jones functionality, prefer:

```bash
sage -python ./strand_passage_guiV3_4.py
```

## Notes

- Backtrack-assisted simplification is on by default: 200 rounds x 30 steps.
  Use `--no-backtrack`, `--backtrack-rounds N`, or `--backtrack-steps K` to
  tune it.
- A negative even DT label means the even visit is the over strand by default.
  Use `--negative-even under` for the opposite convention.
- Generated outputs such as `*.xlsx`, `*_overview*.svg`, and cascade PNGs are
  ignored by Git.
- The optional icon lives in `assets/`. If it is missing, the GUI scripts still
  run normally.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
