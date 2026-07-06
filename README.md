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
strand_passage_guiV3_8.py        Main entry point: GUI, --nongui, and --demo
link_engine_v3_8.py              Diagram engine and SnapPy bridge
draw_dt_original_labelsV3_12.py  DT parser, layout, renderer, and standalone GUI
check_two_dt.py                  Standalone SnapPy/Sage DT-comparison utility
find_link_in_snappy.py           Search SnapPy link databases for DT matches
assets/strand_passage_icon.png   Optional window/task-menu icon
bin/strand-passage               Convenience launcher
requirements.txt                 Python package notes
DEVELOPMENT_LOG.md               Development notes and version history
LICENSE                          MIT license
```

Import chain:

```text
strand_passage_guiV3_8 -> link_engine_v3_8 -> draw_dt_original_labelsV3_12
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
sage -python strand_passage_guiV3_8.py --help
```

SnapPy is intentionally not pinned in `requirements.txt`, because this project
is intended to use the SnapPy/Sage installation on the research machine.

## Run

Interactive GUI:

```bash
sage -python strand_passage_guiV3_8.py
sage -python strand_passage_guiV3_8.py --dt "DT: [(4,6,2)]"
```

If the TkAgg backend is not available:

```bash
python3 strand_passage_guiV3_8.py --gui-backend agg
```

Batch spreadsheet and overview SVG:

```bash
sage -python strand_passage_guiV3_8.py --nongui \
  --dt "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]" \
  --out strand_passage_results.xlsx
```

For one-crossing components, use Python tuple syntax with a trailing comma:
`DT: [(4,), (2,)]`, not `DT: [(4), (2)]`.  The latter is ambiguous because
Python reads `(4)` as the integer `4`.

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

Custom displayed crossing IDs use the same syntax as
`draw_dt_original_labelsV3_12.py`:

```bash
sage -python strand_passage_guiV3_8.py \
  --dt "DT: [(4,6,2)]" \
  --crossing-order "c1 c3 c2"
```

`--crossing-order` lists displayed crossing IDs in odd-label order
`1,3,5,...`.  Alternatively, use `--crossing-map "c1=1,c3=3,c2=5"`.
In the GUI, the same controls are available as `Crossing order` and
`Crossing map`.  These labels affect the drawing and passage notes only; the
internal strand-passage calculation still uses the underlying DT crossing ids.

Spreadsheet columns to know:

- `DT_code_chosen` is the visible DT code used for the drawn structure and for
  any next passage step.
- The second pass is run once for each merged first-step structure,
  not once for every raw first-step crossing. A first-step result is continued
  when `new_components > 2` and it has a usable `DT_code_chosen`; the CLI prints
  this criterion and the resulting continuation counts. Second-pass sheets are
  named `merged_<node>_<labels...>`, and the columns `first_step_passages` and
  `first_step_representative` show which first-step crossings merged together
  and which representative was continued.
- The `run_info` sheet records the software version, engine/drawing modules,
  runtime, command/arguments, output paths, key parameters, second-pass
  continuation criterion, and continuation counts for that workbook.
- Before a pass is enumerated, the source DT is simplified to the representative
  used for that pass's crossing-count baseline. This keeps continuation rows
  from being generated on a higher-crossing stale parent and prevents invalid
  `increase` outcomes.
- `Jones_polynomial` is computed from that exact visible DT code, so pasting
  `DT_code_chosen` into Sage/SnapPy should give the same polynomial.
- `topological_Jones_polynomial` records the SnapPy simplified link object
  before its visible DT export. When `hidden_split_unknots` is greater than
  zero, this polynomial includes the extra split-circle factor, while
  `Jones_polynomial` describes the crossing-bearing DT code that can be drawn.
- `visible_components` is the component count in the exported/drawn DT code.
  `topological_components` is the component count before any zero-crossing
  split unknots are omitted from the DT export.

Headless cascade figure:

```bash
python3 strand_passage_guiV3_8.py --dt "DT: [(4,6,2)]" --demo 2 1 --out chain.png
```

Standalone DT comparison utility:

```bash
sage -python check_two_dt.py
```

SnapPy database search utility:

```bash
sage -python find_link_in_snappy.py --dt "DT: [(4,6,2)]" --output matches.tsv
sage -python find_link_in_snappy.py --input links.tsv --output matches.tsv
```

`find_link_in_snappy.py` searches SnapPy's `HTLinkExteriors` and
`LinkExteriors` databases for links matching one or more DT codes. For each
query it compares, in order: extended alphabetic DT code with flip data, exact
numeric DT code, optional loose numeric DT code, and then meridian-preserving
exterior identification for hyperbolic links. Use `--strict` to disable the
loose numeric fallback. Input files may contain one DT code per line, or
`label<TAB>DT_code`; output is a TSV table.

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
chmod +x strand_passage_guiV3_8.py
chmod +x draw_dt_original_labelsV3_12.py
chmod +x check_two_dt.py
chmod +x find_link_in_snappy.py
```

Then they can be run as:

```bash
./strand_passage_guiV3_8.py --gui-backend agg
```

For SnapPy/Jones functionality, prefer:

```bash
sage -python ./strand_passage_guiV3_8.py
```

## Notes

- Backtrack-assisted simplification is on by default: 200 rounds x 30 steps.
  Use `--no-backtrack`, `--backtrack-rounds N`, or `--backtrack-steps K` to
  tune it.
- A negative even DT label means the even visit is the over strand by default.
  Use `--negative-even under` for the opposite convention.
- Generated outputs such as `*.xlsx`, `*_overview*.svg`, and cascade PNGs are
  ignored by Git.
- The `--nongui` overview SVG keeps labels and captions as editable text using
  Arial, so text can be selected and edited in Illustrator/Inkscape. V3.8 also
  enlarges the surrounding label boxes/circles so the exported SVG better
  matches the live Matplotlib view in Illustrator.
- Standalone SVGs from `draw_dt_original_labelsV3_12.py` use the same Arial
  editable-text policy and roomier DT-label/crossing-ID boxes.
- In the drawing helper GUI, saved SVG font sizes follow the live GUI fields:
  `DT label font` maps to `--font-size`, and `crossing ID font` maps to
  `--crossing-id-font-size`.
- The optional icon lives in `assets/`. If it is missing, the GUI scripts still
  run normally.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
