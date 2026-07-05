Strand-Passage Explorer Development Log
=======================================

This file records version-specific development notes and implementation history.
For installation and everyday usage, see README.md.

Interactive strand-passage explorer with component-colour preservation, plus a
non-interactive spreadsheet mode.  Every crossing you click opens a NEW WINDOW
showing the chosen after-passage topology; each new window is itself clickable,
so passages chain (and may branch) across as many windows as you like.

Files
-----
    strand_passage_guiV3_4.py         (V3.4) interactive GUI + --nongui + --demo
    link_engine_v3_4.py               (V3.4) passage / DT-choice / SnapPy engine
    draw_dt_original_labelsV3_11.py   (V3.11) drawing + model layer
    check_two_dt.py                   standalone SnapPy/Sage utility: compare two
                                      DT codes (topology + Jones + backtrack test)
    assets/strand_passage_icon.png    optional Tk window/task-menu icon
    bin/strand-passage                convenience launcher
    DEVELOPMENT_LOG.md

--nongui outputs
----------------
    <name>.xlsx           the two-pass spreadsheet (first_pass sheet + one sheet
                          per first-flip crossing)
    <name>_overview.svg   a large, tidy overview of every resulting structure

Run
---
    sage -python strand_passage_guiV3_4.py                     # SnapPy enabled
    sage -python strand_passage_guiV3_4.py --dt "DT: [(4,6,2)]"
    python3 strand_passage_guiV3_4.py --gui-backend agg        # if TkAgg won't load
    ./bin/strand-passage --help                                # launcher

Backtrack-assisted simplification (ON by default: 200 rounds, 30 steps):
    sage -python strand_passage_guiV3_4.py                     # backtrack ON
    sage -python strand_passage_guiV3_4.py --backtrack-rounds 400
    sage -python strand_passage_guiV3_4.py --no-backtrack      # turn it OFF
    (in the GUI, the "Backtrack simplify" checkbox + rounds/steps fields start
     ON at 200/30 and can be toggled live between clicks.)

Non-interactive two-pass spreadsheet:
    sage -python strand_passage_guiV3_4.py --nongui \
        --dt "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]" \
        --out strand_passage_results.xlsx        # backtrack ON by default

Headless cascade figure (no display needed):
    python3 strand_passage_guiV3_4.py --dt "DT: [(4,6,2)]" --demo 2 1 --out chain.png


What is new in V3.4 (GitHub packaging + icon)
---------------------------------------------
  * The public README now includes clone, pull, install, run, and executable-use
    notes for a GitHub checkout.
  * The repository includes an MIT license and a small bin/strand-passage
    launcher. The launcher uses sage -python when Sage is on PATH and otherwise
    falls back to python3.
  * draw_dt_original_labelsV3_11.py adds optional Tk window/task-menu icon
    support using assets/strand_passage_icon.png. If the icon file is absent or
    the local Tk build cannot load it, the GUI simply continues without an icon.
  * Bug fix: the GUI and --nongui outputs now report the Jones polynomial of the
    exact DT code shown/drawn and used for the next passage step. This avoids
    mismatches when a transient SnapPy simplified object and its exported
    DT_code do not round-trip with identical Jones values.


What is new in V3.3 (backtrack-assisted simplification)
------------------------------------------------------
SnapPy's simplify('global') is a greedy heuristic and can stop at a NON-minimal
diagram (a real example we hit: a 4-component link whose minimum is 10 crossings
was reported as 12 by simplify('global'); backtracking down to 10 confirmed the
minimum).  This affects every diagram the tool draws and every DT-code choice,
not just the overview merge.  V3.3 adds an OPTIONAL backtrack-assisted
simplification that repeatedly complicates the diagram, then re-simplifies, and
keeps the fewest-crossing result:

  * link_engine_v3_4.backtrack_simplify(snappy, link, mode, rounds, steps)
    is the core routine; rounds<=0 is a plain single simplify (unchanged).
  * It is threaded through the SnapPy path used by BOTH the GUI (advance ->
    snappy_global_simplification) and --nongui (every passage result and the
    original link).
  * Controls -- ON BY DEFAULT (200 rounds, 30 steps):
      CLI:  --no-backtrack to disable; --backtrack-rounds N / --backtrack-steps K
            to tune (--backtrack is kept as a no-op for compatibility).
      GUI:  the "Backtrack simplify" checkbox + "rounds"/"steps" fields in the
            top bar start ON at 200/30 and can be toggled live between clicks.
  * Cost: each round is one backtrack + one simplify, so N rounds ~ N extra
    simplifications per structure; lower the rounds if a large --nongui sweep is
    slow.
  * Because backtrack is randomised, the exact rounds needed vary run to run.

The standalone check_two_dt.py utility demonstrates the effect (plain loop vs.
backtrack, and reports the round at which a lower crossing count is reached).


Overview SVG refinements (V3.3)
-------------------------------
  * Targeted per-step reconciliation (reconcile_steps): within a step, any
    structures that share the same Jones (up to mirror q->1/q and framing q^n)
    and component count but show DIFFERENT crossing counts are the same link left
    at different simplify('global') plateaus.  The higher-crossing ones are
    driven DOWN to the step's minimum with a *targeted* backtrack (it stops the
    moment it reaches that count -- far cheaper than a big global round count),
    then the now-identical structures are merged and their arrows combined.  This
    removes the redundant near-duplicate cards.
  * Structures are deduplicated WITHIN each step only, never across steps.
  * All cards share one fixed size; each column is centred vertically so the
    figure is symmetric.  Diagrams are enlarged (bigger artboard) so DT and
    crossing labels are not crowded; the small filled component-colour dots were
    removed.  Card text (DT / crossings+components / Jones) is blank-line
    separated and set in a larger font.
  * Each arrow has its own colour; its label (the flipped DT crossing) is drawn
    in that colour and placed ON its own arrow (spread along the arc so labels
    never stack), and long labels are wrapped so they don't overrun the panels.
  * Illustrator-friendly export: all text is outlined to vector paths (no more
    "Font Problems" dialog about DejaVu/Computer Modern/etc.) and artist clipping
    is turned off (no more "Clipping will be lost on roundtrip to Tiny").


What changed from V3.1 (the five requested items)
-------------------------------------------------
1. Drawing helper is now draw_dt_original_labelsV3_11.py, imported through
   link_engine_v3_4.py.  The strand-passage views draw 2-D links with the
   helper's OWN default settings: default layout (DEFAULT_LAYOUT = "tutte"),
   top-to-bottom orientation, and the helper's false-crossing audit with a
   planar fallback.  Those defaults now live in one place
   (draw_dt_original_labelsV3_11.DEFAULT_LAYOUT / DEFAULT_Y_DIRECTION /
   DEFAULT_ROTATE) so the engine and the standalone helper stay in sync.

2. --nongui reproduces the old two-pass spreadsheet (sheet "first_pass" plus one
   sheet per first-flip crossing) using SnapPy directly.  DT-code choice rule:
   after a passage, if the DIRECT after-passage DT code has MORE crossings than
   the SnapPy simplify('global') code, the SnapPy-simplified code is used;
   otherwise (fewer-or-equal, i.e. a tie) the direct after-passage code is kept.
   New columns record both codes, the chosen code, the chosen source, and the
   direct / SnapPy / chosen crossing counts.  The second pass continues from the
   CHOSEN code.

3. The same DT-choice rule drives the interactive GUI (advance()).  DT traversal
   labels stay visible in the view even after a SnapPy simplification, because
   the chosen diagram is always built from its own DT code and carries a
   _dt_labels_valid flag; render() shows labels for any such pristine drawing
   (original, direct after-passage, or SnapPy-simplified) and hides them the
   moment the drawing is further passage-changed or reduced.

4. The mouse cursor becomes a pointing hand ("hand2") whenever it is near a
   crossing where a strand passage can be performed (same hit radius as the
   click test), in both the TkAgg and the Agg-in-Tk backends.

5. SVG-vs-live-view difference (diagnosed and FIXED in V3.10):
   * Cause: the knot geometry (deterministic Tutte layout) is IDENTICAL in both
     the live GUI preview and the saved file.  The visible difference was purely
     framing.  The live preview reframes the diagram to the wide Tk canvas
     aspect ratio and applies the interactive zoom (_apply_preview_zoom), while
     the saved image was written from a SQUARE figure cropped with
     bbox_inches='tight'.  Different aspect + zoom + crop margins made the two
     look different even though the drawing was the same.
   * Fix: a single shared routine, apply_content_framing(), now frames both the
     live preview AND the file save with the same content bounds, padding,
     preview zoom, and panel aspect ratio.  draw()/save_image()/save_all() take
     match_view={"aspect", "zoom"} and write the whole figure (no tight crop),
     so the saved SVG reproduces exactly what is on screen.  (The plain CLI
     pipeline keeps its legacy square/tight output, since it has no live view to
     match.)
   GUI layout was also optimised: a resizable diagram/properties split, a
   scrollable properties panel, a hint bar, a larger diagram, and Save PNG / Save
   SVG buttons that save exactly the on-screen view.

Earlier live files are kept only in the local maintenance archive and are not
part of the public GitHub package.


V3.2 refinements (second round)
-------------------------------
a. The strand-passage 2-D views now apply the helper's full default look,
   including crossing IDs coloured by the over-strand's component colour
   (link_engine_v3_2 tags each crossing with its odd/even visit and passes
   color_crossing_ids_by_overstrand=True).
b. draw_dt_original_labelsV3_11: the GUI window title reads "V3_11",
   and the final XYZ smoothing defaults changed to smooth window = 10 and smooth
   passes = 5 (argparse, GUI fields, and function defaults).
c. The pointing-hand cursor now uses an ADAPTIVE radius (just under half the
   nearest-crossing gap, capped small), so it only appears when the cursor is
   genuinely on a crossing instead of far away.
d. Running under plain python (not `sage -python`) prints a warning that Jones
   polynomials -- and the SnapPy colour-matching that relies on them -- cannot
   be computed, and suggests Sage.
e. --nongui now also writes <name>_overview.svg: one large figure of every
   resulting structure across the two passage steps.  Three columns (original /
   after 1 passage / after 2 passages), colour-coded cards (green = fewer
   crossings than the original, gray = same, orange = more, gold = unknot/
   unlink), each card showing the helper-default 2-D view, DT code, crossing and
   component counts, Jones polynomial, a component colour key, an outcome tag,
   and a badge when several passages merge into it.  Curved arrows labelled with
   the flipped DT crossing show the operation order.  Cards never overlap and the
   diagrams are drawn on top of the arrows so nothing is hidden.

   Merging rule: two structures are treated as the SAME card when they share the
   SnapPy-simplified crossing count, the component count, and the Jones
   polynomial UP TO (a) the mirror symmetry q -> 1/q and (b) an overall monomial
   factor q^n.  The mirror image of a link has Jones V(1/q), so a chiral link and
   its mirror merge; an overall q^n factor is only a framing/writhe (Reidemeister
   I) normalization artifact, so those merge too.  The crossing/component/Jones
   values shown
   on each card are exactly the values written to the spreadsheet (the overview
   reuses them rather than recomputing), so the two outputs always agree -- this
   also fixes the earlier case where the overview showed "n/a" Jones for
   structures that had a Jones value in the spreadsheet.


Colour tracking (unchanged from V3.1)
-------------------------------------
When SnapPy's simplified diagram is chosen, component colours are recovered by
matching the linking matrix + per-component Jones signatures over all component
permutations and orientation reversals.  A unique match => original colours are
kept and the log says "colours tracked".  If the outcome has topologically
interchangeable components (e.g. two split unknots after unlinking a Hopf link)
no invariant can label them; the diagram is still shown in default colours and
the log says "COLOURS NOT TRACKED" with the reason.  colours_tracked propagates
to all descendants.


Dependencies / testing notes
----------------------------
* Python 3 with numpy, networkx, matplotlib; tkinter for the GUI.
* SnapPy is optional but is the intended path -- run under `sage -python`.
* --nongui additionally needs pandas + openpyxl (.xlsx writer) and matplotlib
  (overview SVG).
* In environments WITHOUT networkx / SnapPy / pandas:
    - all scripts were syntax-checked (py_compile) and pass;
    - the --nongui DT-choice logic was unit-tested with a mock SnapPy (tie keeps
      the direct code; strictly-fewer crossings picks the SnapPy code);
    - the _dt_labels_valid flag and apply_content_framing() math were unit-tested;
    - the engine block-model odd/even + over-strand colouring was unit-tested;
    - the adaptive cursor radius and the Sage warning were unit-tested;
    - the overview-SVG layout/arrows/merge/badge rendering was exercised with a
      stubbed engine + mock SnapPy and produced a clean, non-overlapping figure.
  A live `sage -python` smoke test on your machine remains the first
  confirmation step: the GUI + SnapPy path (chained multi-window passages,
  colour-tracking messages, Save PNG/SVG matching the view) and a real --nongui
  run (spreadsheet + overview SVG on your 14-crossing 4-component link).
