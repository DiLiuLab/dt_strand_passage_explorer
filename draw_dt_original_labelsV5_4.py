#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_dt_original_labelsV5_4.py
==============================

Draw a smooth planar oriented link diagram from a signed Dowker-Thistlethwaite
(DT) code while preserving the original traversal labels supplied by the user.

V5.4 changes
------------
* 3D XYZ CORRECTNESS: dense diagrams (e.g. the mirror chain and heavy chain
  families) could render as the WRONG LINK: the kamada sphere layouts pinch
  strands through each other and thread foreign crossings, invisibly to the
  2D audits.  Three additions:
  - '--sphere-layout stereo-safe': correct-by-construction 3D curve -- the
    audited planar layout lifted stereographically with radial over/under
    bumps and no position-moving smoothing (strands disjoint in 2D can never
    intersect in 3D).  Use this whenever the audit flags a failure.
  - automatic TOPOLOGY AUDIT after every xyz build (audit_xyz.py, same
    folder; needs spherogram+snappy): rebuilds the link from the 3D curve
    and verifies it is isometric to the DT link; prints [ok]/[FAIL] + hint,
    while the GUI retains the result in its status log and 3D windows.
  - clearance repair for the kamada layouts (--xyz-clearance /
    --no-xyz-repair): separates non-crossing close approaches into extra
    radial layers, pinch-coherently.
* Previous V5.3 notes follow.
* 3D projection surface wireframe is now configurable: 'grid density' (number
  of meridians; parallels follow at about half; 0 hides the mesh -- needs
  'Redraw 3D projection' since the mesh is built with the curve), and 'grid
  color' / 'grid line width' (any matplotlib colour; applied live from the
  cache).  CLI: --proj-grid-density / --proj-grid-color / --proj-grid-lw.
* The neutral framework chords connecting the crossing anchors can be shown or
  hidden independently of the wireframe ('framework chords' checkbox /
  --proj-hide-chords), so the mapped shape can be shown as surface + anchors
  only.
* Mouse roll in the live projection window: RIGHT-drag (or Shift+left-drag)
  rotates the image about the view axis, pivoting on the canvas centre --
  complementing the left-drag azimuth/elevation rotation and wheel zoom.
* Projection crossing anchor dots are now optional and hidden by default
  ('crossing anchor dots' checkbox / --proj-crossing-dots).
* The mapped-surface grid can be locked to the last projection rebuild view
  while live orbit/roll rotates the strands ('fixed grid while rotating'
  checkbox / --proj-fixed-grid).

V5.2 changes
------------
* Mapped skeleton REDONE to match score_diagram's sphere panels: 'show mapped
  skeleton' now renders the WHOLE mapped surface as a light depth-shaded
  wireframe (the actual sphere / ellipsoid / cylinder / torus the strands ride
  on, produced by passing a sphere grid through the same surface warp), plus a
  NEUTRAL dark-gray framework of strand chords and black depth-faded crossing
  anchors.  Component colours are reserved for the smooth strands, so the
  overlay no longer doubles the colours (that was messy).
* New 2D layout 'sphere-stereo' -- the creative route to a FLAT diagram for
  links whose rings live in orthogonal planes (Edwards-Venn AM_n) and defeat
  boundary-pinned layouts: the crossing graph is laid out on the unit sphere
  (3D Kamada-Kawai directions + gadget compaction, as in the 3D pipeline), a
  projection pole is searched that is farthest from every strand arc, and the
  whole spherical diagram is STEREOGRAPHICALLY projected from that pole.  A
  stereographic image of a spherical diagram is a genuine planar diagram, so
  near-orthogonal rings land as clean overlapping ovals instead of collapsing
  onto a pinned boundary.  Available in the GUI layout menu and --layout
  sphere-stereo; combines with 'relax passes' / 'min separation'.

V5.1 changes
------------
* 3D projection, for links whose rings lie in near-orthogonal planes (e.g. the
  Edwards-Venn AM_n family, where every axis view leaves one ring edge-on):
  - 'Auto view' (--proj-auto-view): searches the view sphere for the direction
    that maximizes the WORST component's projected roundness (isoperimetric
    ratio of its projected loop), i.e. the view where no ring degenerates to a
    sliver.  Button in the 3D-view tab and the projection window.
  - 'depth fade' (--proj-depth-fade, 0..1): score_diagram-style depth cueing --
    nearer strand chunks are opaque, farther ones fade.  Applied to strands
    and skeleton alike; makes front/back of an edge-on ring distinguishable.
  - 'perspective' (--proj-perspective, 0 = orthographic, else camera distance
    in scene radii, e.g. 2-4): foreshortening turns an edge-on ring into an
    open lens instead of a flat line, like a real 3D viewer.
* Mapped skeleton is now the WHOLE framework in the style of
  score_diagramV2_0's sphere panels: every strand chord drawn as a smooth arc
  through its connector anchor, in its component colour, with per-segment
  depth transparency, plus depth-faded crossing anchor dots.  New 'skeleton
  only' mode (--proj-skeleton-only) shows just this framework, which makes the
  global SHAPE of a difficult link clear at a glance.
* Projection window: mouse DRAG ROTATES the view (Chimera-like: horizontal
  drag = azimuth, vertical = elevation; live fast redraw while dragging, full
  halo render on release), and the scroll wheel zooms about the view centre.

V5.0 changes
------------
* 2D, holed-tutte: new 'ring equalize' (--ring-equalize, 0..1).  The harmonic
  (Tutte) solve compresses interior detail toward the pinned rims on large
  diagrams (e.g. the 80-crossing closed fishtail), which is why strands crowd
  at the boundary.  Ring equalize histogram-equalizes the crossings' radial
  positions across the ring width (angles kept), spreading the strands evenly
  between the rims.  1 = fully even.  Applied before invert-ring / ring-tilt.
* 2D, all layouts: planarity-guarded relaxation ('relax passes' / 'relax
  strength', --relax-passes / --relax-strength).  Each pass moves non-pinned
  nodes toward their neighbours' centroid while equalizing incident edge
  lengths, and is REJECTED (step halved) if it would introduce any
  straight-edge crossing, so the layout stays planar by construction.  Evens
  out crowded/stretched regions on difficult structures.
* New '3D view' tab + persistent '3D projection' window: an illustrative 2D
  orthographic projection of the 3D Sphere-XYZ curve, drawn depth-ordered with
  white halos so nearer strands read as crossing over farther ones -- a
  presentation method for links with no readable flat diagram (e.g. the
  Edwards-Venn Brunnian links).  The expensive 3D construction runs only on
  'Redraw 3D projection'; the quick-view buttons re-project the cached curve
  instantly.
* 3D view: optional 'mapped skeleton' overlay (crossing anchors of the sphere
  layout + strand chords, warped like the curve, with optional crossing IDs)
  to relate the 3D picture crossing-by-crossing to the 2D diagram.
* 'Save projection(s)' writes selected views (current/top/front/side) as
  SVG+PNG next to the XYZ file; 'Save XYZ' can include them automatically.
  CLI: --save-projections, --proj-views, --proj-elev/azim/roll,
  --proj-skeleton, --proj-skeleton-ids, --proj-line-width.
* UI: 'View XYZ' moved to the '3D view' tab; save-button row decluttered;
  sessions include all new parameters and still load V4.x session files.

V4.5 changes
------------
* GUI: the right-hand parameter panel is now split into two tabs -- "2D diagram"
  (layout / tutte / holed-tutte / drawing options / crossing order & map) and
  "3D XYZ" (sphere layout and surface / smoothing parameters) -- so only the
  relevant fields are shown at a time instead of one long scrolling list.  Each
  tab scrolls independently.

V4.4 changes
------------
* holed-tutte now builds the balanced 'closed principal curve' from a 3D torus
  layout: a 3D Kamada-Kawai layout naturally arranges a symmetric multi-component
  link on a torus, and its azimuth about the torus axis (the longitudinal 'long way
  around', recovered by PCA) is spread far more evenly than a 2D layout -- so links
  like the 4BL / C5BL / L6BL / 6-component come out as near-perfectly symmetric
  wreaths.  The top-down projection of that torus is the annulus (outer equator ->
  outer rim, inner equator -> inner rim); its two boundary cycles are pinned to the
  holed shape and the interior is Tutte-solved, then refined by a few fixed-point
  rounds.  When a link is not annular this way (or SciPy is missing) it falls back
  to the circular-Tutte probe, which still forces a clean annulus.
* GUI: the Output-image and Sphere-XYZ path fields were removed (the Save buttons
  prompt for the directory and filename), and the 'Save all' button was removed.
* holed-tutte: the 3D torus layout is now deterministic (seeded) and cached per
  crossing graph, so it is only recomputed when the DT code changes -- tweaking
  rotate / swap / hole ratio / ring tilt no longer re-rolls the 3D layout, and the
  diagram stays put instead of jumping.
* Orientation arrowheads now scale with the strand line width (in both the live
  preview and the saved image), instead of a fixed size.
* holed-tutte: new 'invert ring (inside-out)' option (--invert-ring) that turns the
  ring inside-out by reflecting each crossing's radius about the ring mid-line, so
  the inner boundary ends up outside and vice versa.  This is distinct from 'swap
  inner/outer face', which instead re-solves with the two boundary faces exchanged.
* GUI: the Save image / Save table / Save XYZ buttons now always open a Save-As
  dialog for the directory and filename (like 'Save as session').
* GUI: 'Refresh 2D' is now 'Redraw 2D' -- it clears the cached 3D-torus layout and
  rerolls its seed, so holed-tutte is recomputed from scratch (a fresh layout) even
  when no parameter changed.  Normal parameter tweaks still reuse the cached layout.
* GUI: the DT code box is kept on a single line -- line breaks from a multi-line
  paste are stripped automatically.
* The DT parser tolerates look-alike Unicode dashes (math minus, en/em dash, ...)
  and non-ASCII spaces pasted from documents.
* holed-tutte 'ring tilt' (0-90 deg) is now a bucket view: the flat annulus is the
  top-down view of a bucket whose wall carries the crossings.  At 90 (default) you
  look straight down -> the flat top-view wreath (closed principal curve around the
  whole circumference).  Lower tilts rotate the bucket about the horizontal axis so
  the wall opens toward a side view.  The tilt lifts each crossing onto the 3D
  bucket wall and projects it, so it is intrinsic to the mapping, not a post-hoc
  rotation of the finished picture.
* The shaped-tutte / holed-tutte shape menu drops the redundant 'circle' and
  'rectangle' choices: an ellipse at aspect 1 already is a circle, and a rounded
  rectangle at corner radius 0 already is a sharp rectangle.  Only 'ellipse' and
  'rounded-rectangle' are offered now (old saved sessions naming circle/rectangle
  still load and render).
* GUI: 'Save session' is now 'Save as session', and 'Quit' moved to a second row
  of the save-button area.

V4.3 changes
------------
* New standalone 2D layout "holed-tutte": instead of pinning a single outer face
  to a convex boundary, it pins TWO faces -- the outer face and an auto-picked
  central 'hole' face -- to the outer and inner outlines of a holed shape
  (annulus / elliptical or rectangular ring) and harmonically solves the ring, so
  the diagram wraps around a central hole along its natural closed principal
  curved axis.  Controls: shape / aspect / orient (as in shaped-tutte), a 'hole
  ratio' (inner/outer size), and a 'swap inner/outer' checkbox.  The shape-outline
  overlay shows both ring outlines; the PCA-axis overlay shows the mid-ring curved
  axis.
* GUI: when 'tutte auto aspect' is on, the panel now shows the computed aspect
  value next to the checkbox (for shaped-tutte and holed-tutte).
* Crossing-ID circles are now sized to snugly enclose the ID text (measured text
  extent plus a small margin) instead of a fixed, over-large data fraction.
* New 'min separation' control: a post-layout relaxation pushes apart non-incident
  strand pieces that sit closer than the given fraction of the diagram span (with
  a spring back to the original layout), opening up shallow / near-parallel runs.
  0 = off.

V4.2 changes
------------
* GUI: new "Save session" / "Load session" buttons write and restore every
  setting/parameter on the panel to a JSON file, so a full working session can be
  reproduced later.
* GUI: the signed DT code (?) help now states that the default code is the 4BL
  diagram and lists ready-to-copy example codes (TK, HL, BR, C5BL, L6BL); the help
  popup is sized larger to hold them.
* Crossing-ID circles are now drawn at the same data-space size in both the live
  preview and the saved image, and remain true circle objects in the SVG (an
  unfilled Circle patch) rather than being flattened into an expanded path.
* GUI: editing any Sphere XYZ (3D-only) setting no longer triggers an automatic
  2D preview redraw, so those fields respond smoothly without lag.  The 2D preview
  still refreshes for every setting that actually affects it.

V4.1 changes
------------
* Auto orientation and auto aspect are independent toggles for both the 2D
  shaped-tutte layout and the 3D shaped-kamada surface.  Either can be auto while
  the other is set manually.  (Auto aspect measures the diagram's own elongation
  via PCA -- the ratio of its spread along its long vs. short principal axis --
  and sizes the boundary / surface proportions to match it; it does not equalize
  strand lengths.)
* 2D shaped-tutte orientation is always relative to the diagram's PCA elongation
  axis and works for every boundary shape.  Auto orient puts the PCA axis along
  the view x-axis; when it is off, a manual 'tutte orient' angle rotates from the
  PCA axis.  The global 'rotate degrees' still applies on top.
* 3D shaped-kamada has two manual orientation controls, both relative to the PCA
  axis and used when auto orient is off: 'surface orient' spins the mapping about
  the surface's primary axis, and 'surface tilt' rotates that primary axis away
  from the PCA axis (about the perpendicular secondary axis).
* 2D decompression now has two independent controls: 'tutte decompress' (the
  original boundary-depth weighting, referenced to the boundary/shape) and a new
  'tutte COM expand' that radially expands interior structure about the crossing
  center of mass so crowded crossings get more room and sit near the center.
* Help (?) text for the new shaped-layout parameters states the suggested range
  and whether the value is an absolute distance or a dimensionless ratio; the
  short hints that used to sit under some fields were folded into their ? help.
* The figure-size help now reports the current live-preview panel size.
* The DPI field was removed from the GUI (vector SVG/PDF output is the common
  case); PNG raster export still uses the command-line --dpi default.

V4.0 changes
------------
* 2D: new "shaped-tutte" layout.  The Tutte (barycentric) boundary can now be a
  circle, an ellipse, a rectangle, or a rounded rectangle with an adjustable
  aspect ratio and corner radius, instead of only a circle.  With auto mode on
  (default), the boundary aspect is derived from the diagram's own elongation
  and the finished diagram is rotated so its principal (PCA) axis lies
  horizontally, i.e. it is automatically oriented along the elongation axis.
* 2D: both "tutte" and "shaped-tutte" gain an internal-decompression weight.
  A single strength value (0 = classic Tutte) applies boundary-depth edge
  weights in a convex-combination Tutte solve, which pushes interior structure
  outward so it is not over-compressed toward the center.  Positive weights keep
  the barycentric map valid, so no new crossings are introduced by the weighting.
* 3D: new "shaped-kamada" sphere layout.  The sphere-native spherical-kamada
  construction is warped onto other closed/oriented surfaces -- an ellipsoid, a
  cylinder, or a torus -- with the over/under crossing offset applied along the
  local surface normal.  Auto mode (default) aligns the surface's principal axis
  to the diagram (3D PCA) and derives the surface aspect / tube ratio from the
  diagram's own shape; manual mode exposes the shape parameters directly.
* GUI: parameter fields that do not apply to the current 2D layout / sphere
  layout / shape are now dynamically greyed out (disabled), so only relevant
  controls are active.

V3.14 changes
-------------
* The chosen 2D layout is now kept as requested even when it introduces false
  crossings -- there is no automatic fallback to 'planar'.  This preserves, for
  example, a "failed" kamada layout for inspection.
* False crossings (accidental strand overlaps a layout adds away from true DT
  crossings) are now drawn with the same over/under gap style as real crossings.
  The over/under choice there is arbitrary (chosen deterministically), since a
  false crossing is not a real DT crossing.
* A red warning listing the number of false crossings is shown on the live
  preview and saved into the exported image/SVG.
* The diagram metadata (script name/version and essential facts) is now also
  drawn as a visible caption at the bottom of the saved image, in addition to
  being written as file metadata.

V3.13 changes
-------------
* Over/under gaps now render correctly on self-crossings (e.g. the trefoil
  DT: [(4,6,2)]).  Local crossing pieces are selected by arc length along the
  curve parameter instead of by spatial distance, so a component that passes
  through one crossing twice no longer re-covers its own under-strand gap.
* Saved diagrams embed metadata: the generating script name and version plus
  essential facts (component and crossing counts and the key drawing
  parameters).  For SVG these appear in the file's <metadata> block.
* GUI: the Signed DT code box (now two lines) and the output-image / Sphere-XYZ
  path fields sit at the top of the left column, directly above the live
  preview, where long values have full width.  The Sphere XYZ parameter section
  moved to the bottom of the right panel.  The mouse wheel now scrolls the right
  parameter panel.  The CSV table field was removed -- "Save table" writes the
  table beside the output image with the same name and a .csv extension.  The
  redundant "Refresh preview" button was removed ("Refresh 2D" remains).

V3.12 changes
-------------
* Editable SVG/PDF text now uses Arial consistently during preview rendering and
  file export, reducing font-metric differences when the result is opened in
  Illustrator.
* DT traversal-label boxes and crossing-ID circles have larger default padding,
  so their outlines do not look too tight around editable text in the final SVG.

V3.11 changes
-------------
* Adds optional Tk window/task-menu icon support.  GUI mode tries to load
  ``assets/strand_passage_icon.png`` when it is present, but silently continues
  without an icon if the asset is missing or the local Tk build cannot read it.

V3.10 changes
-------------
* The saved 2-D image (SVG/PDF/PNG) now reproduces the live GUI preview.  Both
  paths share one content-framing routine (``apply_content_framing``): the live
  preview and the file-save use the same content bounds, padding, preview zoom,
  and panel aspect ratio, so "what you see is what you save".  Previously the
  preview reframed the diagram to the wide Tk canvas (aspect ratio + interactive
  zoom) while the file was written from a square figure cropped with
  ``bbox_inches='tight'``, so the two looked different even though the knot
  geometry (Tutte layout) was identical.
* Adds module-level defaults ``DEFAULT_LAYOUT`` / ``DEFAULT_Y_DIRECTION`` /
  ``DEFAULT_ROTATE`` as a single source of truth for the 2-D drawing settings,
  which external consumers (e.g. the strand-passage engine) reuse so their
  diagrams match this helper's default output.

V3.9 exports a spherical 3D polyline coordinate file. By default, the
spherical coordinates are produced from a 3D Kamada-Kawai-style graph layout and
then projected to a sphere. Direct connecting remains on by default for
spherical-kamada, so connector segments travel smoothly between the inner and
outer crossing layers instead of returning to the middle sphere layer. V3.9
keeps the in-GUI 3D XYZ viewer and 2D preview zoom/rotation controls.
V3.9 also improves SVG post-editability: text stays as SVG text, component
strands are exported as closed stroked paths with local white gap masks, and
clipping masks are avoided.  Crossing IDs stay fixed at their true crossing
centers, while DT labels are moved outward when needed so the original traversal
labels remain visible. The older stereographic mapping from the 2D drawing is
still available with --sphere-layout stereographic.

Inputs:
  --dt       signed DT code, for example
             'DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]'
  --output   output image path (.svg, .pdf, .png, etc.)
  --table    optional CSV mapping table

Outputs:
  A link diagram image with original odd/even DT traversal labels, optional
  crossing IDs, orientation arrows, over/under gaps, an optional CSV table,
  and a spherical XYZ coordinate file with blank lines between components.

Example:
  python draw_dt_original_labelsV5_4.py \
    --dt 'DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]' \
    --output example_v3.svg \
    --table example_v3.csv \
    --xyz-output example_v3.xyz \
    --show-crossing-ids \
    --crossing-order 'c1 c7 c14 c12 c3 c6 c9 c5 c11 c13 c4 c2 c10 c8' \
    --color-crossing-ids-by-overstrand

If no command-line arguments are supplied, or if --gui is supplied, the script
opens a Tkinter GUI with a live diagram preview panel.

Python 3.9 compatible. Dependencies: numpy, networkx, matplotlib.
Tkinter is needed only for GUI mode.
"""

import argparse
import ast
import csv
import io
import time
import json
import math
import os
import re
import sys

import numpy as np
import networkx as nx
import matplotlib
if not (len(sys.argv) == 1 or "--gui" in sys.argv):
    matplotlib.use("Agg")
# Keep SVG/PDF text as editable font text instead of converting labels to
# outline paths.  This makes labels easier to adjust in Illustrator/Inkscape.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
DIAGRAM_FONT_FAMILY = "Arial"
SCRIPT_VERSION = "V5.4"
VERSION = SCRIPT_VERSION
DT_LABEL_BOX_PAD = 0.22
CROSSING_ID_BOX_PAD = 0.28
matplotlib.rcParams["font.family"] = DIAGRAM_FONT_FAMILY
matplotlib.rcParams["font.sans-serif"] = [DIAGRAM_FONT_FAMILY]
matplotlib.rcParams["font.monospace"] = [DIAGRAM_FONT_FAMILY]
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle


EXAMPLE_DT = (
    "DT: [(-8,-12,16),(-24,-22,-28,-26),"
    "(-10,-14,-2),(-20,-6,-18,-4)]"
)


GUI_HELP_TEXT = {
    "dt": (
        "Signed Dowker-Thistlethwaite code. Multi-component links use one "
        "tuple/list per component.\n\n"
        "The default code loaded in the box is the 4BL diagram:\n"
        "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]\n\n"
        "Other example codes (copy and paste any line into the DT box):\n\n"
        "TK:   [(4,6,2)]\n\n"
        "HL:   [(4,), (2,)]\n\n"
        "BR:   [(6, 8), (12, 10), (2, 4)]\n\n"
        "C5BL: [(18,-28,-66,74,-20,26,-72,80),(34,4,-12,42,-36,10,-2,-44),"
        "(50,30,-22,58,-64,24,-32,-56),(46,68,-78,-48,40,-70,76,-38),"
        "(-6,16,54,-62,8,-14,60,-52)]\n\n"
        "L6BL: [(14,18,-24,-16,-32,26),(2,-8),(6,34,-40,-4,10,-48,42,-12),"
        "(22,50,-58,-20,28,-52,56,-30),(62,38,-46,-64,44,-36),(60,-54)]"
    ),
    "output": "Path for the 2D diagram image. The extension controls the format: .svg, .pdf, .png, etc.",
    "table": "Optional CSV table listing each crossing, original odd/even DT labels, components, over/under label, and 2D coordinates.",
    "xyz": "Path for the spherical x y z coordinate file. The file is plain three-column coordinates. Blank lines separate link components.",
    "negative_even": "DT sign convention. Default: a negative even DT label means the even-labeled visit is over, matching the common Sage/KnotTheory convention. Choose 'under' for the opposite convention.",
    "layout": "2D preview/image layout. tutte is usually smooth and clean; shaped-tutte pins the Tutte boundary to a chosen shape (see tutte shape); holed-tutte pins two faces (outer + auto-picked central hole) to the outer/inner outlines of a holed shape (ring/annulus) so the diagram wraps around a central hole; sphere-stereo lays the diagram out on the unit sphere (3D Kamada-Kawai, as the 3D XYZ pipeline) and stereographically projects it from the pole farthest from all strands -- the method of choice for links whose rings lie in near-orthogonal planes (e.g. Edwards-Venn AM5), which collapse in boundary-pinned layouts (the pole's face becomes the outer region; combine with relax passes / min separation); planar is safest; spring and kamada are aesthetic force-directed alternatives that are audited for false crossings.",
    "tutte_shape": "Boundary shape for the shaped-tutte and holed-tutte layouts: circle, ellipse, rectangle, or rounded-rectangle. All shapes are convex, so the Tutte solve stays valid. For holed-tutte this is the shape of both the outer and inner ring outlines.",
    "hole_ratio": "Holed-tutte inner (hole) outline size as a fraction of the outer outline, in (0,1). Smaller values give a bigger hole/thinner ring wall; larger values give a small hole. Only used when layout is holed-tutte. Default: 0.4.",
    "hole_swap": "Holed-tutte: swap which of the two chosen boundary FACES is pinned to the OUTER outline and which to the INNER (hole) outline, then re-solve. This is a structural change (a different face ends up on the rim), NOT a simple inside-out flip -- for that, use 'invert ring'. Only used when layout is holed-tutte.",
    "invert_ring": "Holed-tutte: turn the ring inside-out (invert inner/outer outlines).  It reflects every crossing's radius about the ring's mid-line, so whatever is currently near the inner hole ends up on the outside and vice versa -- the same diagram, flipped radially (the mid-line curve is unchanged).  Unlike 'swap inner/outer face' (which re-solves with the two boundary faces exchanged), this keeps the layout and just flips it. Only used when layout is holed-tutte.",
    "ring_tilt": "Holed-tutte ring tilt in degrees (0-90): the viewing angle of the ring seen as the wall of a bucket. The flat annulus is the top-down view of a bucket whose wall carries the crossings. At 90 (default) you look straight down the bucket -> the flat top-view wreath, with the closed principal curve running around the whole circumference. As you lower it, the bucket rotates about the horizontal axis so its wall opens toward a side view (0 = fully side-on). The tilt lifts each crossing onto the 3D bucket wall and projects it, so it redistributes the crossings along the tilted ring rather than rigidly rotating the finished picture. Only used when layout is holed-tutte.",
    "min_sep": "Minimum separation between non-incident strand pieces, as a fraction of the diagram span (equivalently discourages very shallow near-parallel crossings). A post-layout relaxation pushes closer pieces apart while a spring keeps the overall layout; larger values separate more but can distort. 0 = off. Try 0.02-0.05. Applies to all 2D layouts.",
    "ring_equalize": "holed-tutte only. Strength 0..1 of a radial histogram equalization applied after the harmonic ring solve: each crossing keeps its angle but its radial position across the ring width is remapped so the radial distribution becomes uniform (1 = fully even, 0 = off, values between blend). This directly decongests diagrams whose strands crowd at the rims (typical for large closed structures such as the fishtail closure, because the barycentric solve compresses interior detail exponentially toward the pinned boundaries). May introduce flagged false crossings at extreme values; combine with 'relax passes' and 'min separation'. Default: 0.0.",
    "relax_passes": "Number of planarity-guarded relaxation passes applied after ANY 2D layout (0 = off). Each pass moves every non-pinned node toward the centroid of its graph neighbours while equalizing incident edge lengths toward the median; a pass that would create a straight-edge crossing is rejected and retried at half strength, so the layout can never lose planarity. Use 5-30 for congested diagrams. Runtime grows with passes x edges^2 (bounding-box pruned). Default: 0.",
    "relax_strength": "Step size 0..1 of each relaxation pass (see 'relax passes'). Larger is faster but rejected more often near tight spots. Default: 0.5.",
    "proj_elev": "3D projection viewpoint elevation in degrees (matplotlib convention: 90 = looking straight down the +z axis, 0 = side view). Used by the 'current' view of the 3D projection window and saved projections. Default: 25.",
    "proj_azim": "3D projection viewpoint azimuth in degrees (rotation about the z axis). Default: -60.",
    "proj_roll": "3D projection roll in degrees about the viewing direction (rotates the projected image). Default: 0.",
    "proj_line_width": "Line width of the projected strands. The white occlusion halo scales with it (2.8x), so thicker lines also produce wider under-strand gaps at apparent crossings. Default: 2.2.",
    "proj_skeleton": "Overlay the mapped FRAMEWORK on the projection: a light depth-shaded wireframe of the whole mapped surface (the sphere / ellipsoid / cylinder / torus the strands ride on), plus optional neutral framework chords and crossing anchor dots. Component colours stay reserved for the smooth strands. Kamada sphere layouts only (stereographic has no skeleton).",
    "proj_skeleton_ids": "Label each skeleton anchor with its crossing ID (same IDs as the 2D diagram / crossing order box).",
    "proj_views": "Comma-separated list of projections written by 'Save projection(s)' and --save-projections. Choices: current (the elev/azim/roll fields), top (down +z), front (down +y), side (down +x). Example: current,top,front,side. Default: current,top.",
    "proj_save_with_xyz": "When checked, 'Save XYZ' also writes the selected projection views (SVG+PNG) next to the chosen .xyz file, so the 2D projections are saved alongside the 3D coordinates in one step.",
    "proj_depth_fade": "Depth cueing strength 0..1 for the 3D projection: 0 = uniform opacity; larger values make farther strand chunks (and skeleton parts) more transparent, nearer ones opaque -- the score_diagram sphere-panel effect. Helps disentangle rings seen edge-on. Default: 0.55.",
    "proj_perspective": "0 = orthographic projection. A positive value switches to a perspective camera at that distance in scene radii (try 2-4): foreshortening opens up rings that are edge-on in orthographic view, like a real 3D viewer. Smaller = stronger perspective. Default: 0.",
    "proj_auto_view": "Search the view sphere (10-degree grid + 2-degree refinement) for the direction that maximizes the WORST component's projected roundness (isoperimetric ratio of its projected loop). This avoids viewpoints where any ring degenerates into a sliver -- the failure mode of links whose rings lie in near-orthogonal planes. Sets the elevation/azimuth fields.",
    "proj_skeleton_only": "Hide the smooth strands and show ONLY the mapped framework (skeleton arcs + anchors). With depth fade this gives the clearest picture of the global shape of a difficult structure.",
    "proj_grid_density": "Number of meridian lines of the mapped-surface wireframe (parallels follow at about half). 0 hides the mesh entirely. The mesh is generated together with the 3D curve, so changing this takes effect on the next 'Redraw 3D projection'. Default: 24.",
    "proj_grid_color": "Colour of the surface wireframe. Any matplotlib colour: gray levels as '0.82' (0=black, 1=white), names like 'lightsteelblue', or hex '#d0e0f0'. Applied live (cached curve is re-rendered). Default: 0.82.",
    "proj_grid_lw": "Line width of the surface wireframe. Applied live. Default: 0.4.",
    "proj_chords": "Show the neutral dark-gray framework chords connecting the crossing anchors. Uncheck to show only the surface wireframe plus any optional anchor dots/IDs (+ strands). Applied live. CLI: --proj-hide-chords.",
    "proj_crossing_dots": "Show small black dots at the crossing anchors in the mapped framework. Hidden by default so dense projections do not get peppered with dots. Applied live. CLI: --proj-crossing-dots.",
    "proj_fixed_grid": "Keep the mapped-surface grid fixed to the view used by the last 'Redraw 3D projection' while live orbit/roll rotates the strands and framework. This gives a steady reference mesh during rotation. Applied live. CLI: --proj-fixed-grid.",
    "xyz_clearance": "Minimum 3D distance between non-neighboring strand pieces. 0 selects the automatic value (0.8 x crossing offset). The V5.4 repair pass uses extra radial layers to separate crowded non-crossing approaches. CLI: --xyz-clearance.",
    "xyz_repair": "Run the V5.4 clearance repair before the topology audit. Keep this enabled for Kamada sphere layouts. Stereo-safe disables repair and smoothing internally because moving its angular positions would void its correct-by-construction guarantee. CLI: --no-xyz-repair disables it.",
    "redraw_projection": "Recompute and topology-audit the 3D Sphere-XYZ curve with the CURRENT '3D XYZ' tab parameters, then redraw the projection window. The audit result appears in both the main status log and the projection window. Quick-view buttons only re-project the audited cached curve and are instant.",
    "tutte_aspect": "Dimensionless ratio (long axis / short axis) of the shaped-tutte ellipse or rectangle boundary. Suggested range 1.0-4.0 (1.0 = round/square). Used only when layout is shaped-tutte, the shape is not a circle, and auto aspect is off. Default: 1.8.",
    "tutte_corner_radius": "Dimensionless ratio 0.0-1.0 giving the rounded-rectangle corner radius as a fraction of the short half-extent. 0 = sharp rectangle, 1 = fully rounded (stadium). Only used for the rounded-rectangle shape. Default: 0.25.",
    "tutte_decompress": "Dimensionless internal decompression strength for tutte and shaped-tutte, measured RELATIVE TO THE BOUNDARY (not a geometric shape-center push). 0 reproduces the classic Tutte layout; suggested range 0.0-1.0 (values much above ~1 can push interior nodes far enough out to create flagged false crossings). It applies boundary-depth edge weights (graph distance from the outer face) that push deep interior structure outward toward the boundary. It does not target where crossings are crowded -- for that, use 'tutte COM expand'. Default: 0.0.",
    "tutte_com_expand": "Dimensionless strength of a SEPARATE radial expansion about the crossing center of mass, for tutte and shaped-tutte. 0 = off. Unlike 'tutte decompress' (which references the boundary), this expands interior structure outward from the density-weighted crowded-crossing centroid, tapering to zero at the pinned boundary, so crowded crossings get more room and sit near the center. It is layered on top of 'tutte decompress'. Suggested range 0.0-1.0; larger values can create flagged false crossings. Default: 0.0.",
    "tutte_auto_aspect": "Shaped-tutte auto aspect. When on, the boundary long/short ratio is measured from the diagram's own elongation (the PCA spread ratio of a circular Tutte solve) so a naturally long-thin knot gets a long-thin boundary. It does NOT equalize strand lengths. When off, the manual tutte aspect is used. Independent of orientation.",
    "tutte_auto_orient": "Shaped-tutte auto orient (independent of auto aspect). Controls the on-screen framing only. When on, the whole diagram is rotated so its intrinsic (circular-Tutte) PCA elongation axis is horizontal; the tutte shape outline then appears tilted by 'tutte orient deg'. When off, the shape's long axis is left horizontal instead (the PCA axis then appears tilted). The global 'rotate degrees' field is applied on top of either mode.",
    "tutte_orient": "Shaped-tutte shape tilt in degrees: the angle by which the tutte shape's aspect (long) axis is tilted AWAY FROM the diagram's intrinsic PCA elongation axis. 0 = shape long axis aligned with the PCA axis. Unlike a final rotation, this actually re-stretches the diagram into a boundary pointed in a new direction relative to the diagram's natural elongation. Only meaningful for non-circular shapes (a circle has no aspect axis). Turn on 'show shape outline' and 'show PCA axis' to see the angle between them. Any value; wraps at 360.",
    "show_tutte_outline": "Overlay the layout's shape outline(s), in both the live preview and any saved image. For shaped-tutte: the boundary outline and its aspect (long) axis. For holed-tutte: both the outer and inner ring outlines. Only shown for the shaped-tutte / holed-tutte layouts.",
    "show_tutte_pca": "Overlay the layout's principal axis, in both the live preview and any saved image. For shaped-tutte: the diagram's intrinsic (circular-Tutte) PCA elongation axis (with 'show shape outline', the angle between the two equals 'tutte orient deg'). For holed-tutte: the mid-ring closed principal curved axis. Only shown for the shaped-tutte / holed-tutte layouts.",
    "surface_shape": "Target surface for the shaped-kamada sphere layout: ellipsoid, cylinder, or torus. The spherical construction is warped onto this surface and the over/under crossing offset follows the local surface normal. Only used when sphere layout is shaped-kamada.",
    "surface_auto_orient": "Shaped-kamada auto orient. When on, the surface's primary axis (ellipsoid major, cylinder length, torus symmetry) is aligned to the diagram's principal (3D PCA) axis. When off, 'surface orient deg' spins the mapping about that primary axis, measured from the PCA axis. Independent of auto aspect.",
    "surface_auto_aspect": "Shaped-kamada auto aspect. When on, the surface proportions (ellipsoid axis ratios, cylinder length/radius, or torus tube ratio) are derived from the diagram's own 3D PCA magnitudes so the surface is as elongated as the diagram. It does NOT equalize strand lengths. When off, the manual surface aspect / tube values are used. Independent of orientation.",
    "surface_aspect": "Dimensionless ratio for shaped-kamada when auto aspect is off: ellipsoid major/minor, or cylinder length/radius. Suggested range 1.0-3.0. Ignored by the torus. Default: 1.6.",
    "surface_tube": "Dimensionless torus tube ratio (minor radius / major radius) for shaped-kamada when auto aspect is off. Suggested range 0.15-0.60. Only used by the torus surface. Default: 0.35.",
    "surface_orient": "Manual shaped-kamada spin in degrees (any value; wraps at 360). Rotates the mapping ABOUT the surface's primary axis, measured from the diagram's PCA axis (a spin around the PCA axis). Used only when surface auto orient is off. Default: 0.",
    "surface_tilt": "Manual shaped-kamada tilt in degrees (any value; wraps at 360). Rotates the surface's primary axis AWAY from the diagram's PCA axis, about the perpendicular secondary axis, so the mapping is no longer aligned with the PCA axis. Complements 'surface orient deg' (which spins about it). Used only when surface auto orient is off. Default: 0.",
    "y_direction": "Controls the final 2D coordinate convention. top-to-bottom is the default drawing orientation; bottom-to-top flips it.",
    "rotate": "Rotates the final 2D scheme by this many degrees, applied on top of the shaped-tutte orientation (see 'tutte auto orient' / 'tutte orient deg'). Affects the 2D drawing and the stereographic sphere layout, but not the native spherical-kamada / shaped-kamada XYZ layout.",
    "dpi": "Raster resolution for PNG output. SVG/PDF are vector formats and are mostly unaffected.",
    "figsize": "Square Matplotlib figure size in inches for the saved 2D image (absolute size, not a ratio).",
    "font_size": "Font size for the original DT traversal labels such as 1, -8, 16. Default: 7.",
    "crossing_id_font_size": "Font size for displayed crossing IDs such as c1, c7, c14. Default: 6.",
    "line_width": "2D strand line width in Matplotlib points. Default: 2.0.",
    "gap_frac": "Under-strand gap size in the 2D image. This is a ratio of the overall 2D diagram span, not an absolute coordinate distance. Default: 0.025.",
    "sphere_layout": "XYZ sphere layout. spherical-kamada distributes the graph directly over the sphere and is best for symmetric spherical models. shaped-kamada warps that construction onto a shaped surface. stereographic maps the current 2D drawing onto a sphere. stereo-safe uses the audited planar layout, stereographically lifts it, and adds only radial crossing bumps; it is the correct-by-construction choice when a dense Kamada layout fails the topology audit.",
    "sphere_radius": "Base radius of the sphere in XYZ coordinate units. Default: 50.0. With crossing offset = 0, all XYZ points lie at this radius.",
    "sphere_extent": "Only used by stereographic. Dimensionless. After centering the 2D drawing, the farthest planar point is scaled to this radius before inverse stereographic projection. 1.0 reaches the equator; larger values use more of the southern hemisphere.",
    "crossing_offset": "Absolute radial offset in XYZ coordinate units, not a ratio. Over-layer radius is R + offset; under-layer radius is R - offset. Default: 5.0. Use 0 for a perfect sphere with no height separation.",
    "sphere_bump_frac": "Dimensionless ratio, not an absolute distance. It does not change the 2D preview. Without direct connecting, it controls local crossing bump width: 1.0 uses the whole crossing patch and smaller values narrow the bump. With direct connecting, local crossing arcs stay on their inner/outer layer; for dip-to-bump or bump-to-dip connectors this fraction controls how much of the connector is used for the smooth radial transition. Default: 1.0.",
    "sphere_crossing_angle": "Angular half-size of each local crossing patch, in degrees. Larger values make crossings more open on the sphere. Default: 15 degrees.",
    "direct_connecting": "Only for spherical-kamada XYZ export. Default: on. If off, connectors return to the middle sphere layer between localized crossing bumps/dips. If on, each segment from crossing to crossing directly connects its endpoint layers: dip-to-dip stays on the inner layer, bump-to-bump stays on the outer layer, and dip-to-bump or bump-to-dip transitions smoothly between layers.",
    "xyz_spacing": "Target separation between adjacent XYZ points, in the same coordinate units as the sphere radius. Default: 1.8. The script resamples each closed component so points are distributed evenly along the curve; the exact spacing is adjusted slightly so each component closes cleanly.",
    "xyz_final_smooth": "Final smoothing for the written 3D XYZ curve. This is applied after the spherical construction and even point redistribution. It smooths unit directions and radial offsets separately, so crossing layers remain recognizable while small stitched-arc kinks are reduced. Default: on.",
    "xyz_smooth_window": "Final smoothing window, measured in number of adjacent XYZ points on each side. Larger values make smoother curves but may soften local features. Default: 10.",
    "xyz_smooth_passes": "Number of repeated final smoothing passes. Larger values smooth more strongly. Default: 5. Use 0 or uncheck final smooth to disable.",
    "xyz_decimals": "Number of decimal places written for each x y z coordinate.",
    "xyz_close_components": "If checked, repeat the first coordinate at the end of each component block so plotting tools display closed loops. Default: unchecked/off.",
    "view_xyz": "Open an interactive 3D preview window for the current Sphere XYZ parameters. The viewer computes the XYZ curves in memory and does not require saving the XYZ file first.",
    "show_crossing_ids": "Draw displayed crossing IDs such as c1, c2, etc. These IDs stay fixed at the true crossing centers; if they overlap DT labels, the DT labels move outward. IDs can be remapped with crossing order or explicit map.",
    "color_crossing_ids_by_overstrand": "Color each crossing ID by the component color of the over-strand. This helps identify the top strand when the label sits near the crossing.",
    "hide_labels": "Hide original DT traversal labels in the 2D image. The internal mapping and CSV still preserve them.",
    "no_arrows": "Hide orientation arrows on components in the 2D image.",
    "title": "Optional title printed above the 2D diagram.",
    "crossing_order": "Optional list of displayed crossing IDs ordered by odd labels 1,3,5,... Example:\n\nc1 c7 c14 c12 c3 c6 c9 c5 c11 c13 c4 c2 c10 c8",
    "crossing_map": "Alternative explicit mapping from displayed crossing ID to odd traversal label. Example:\n\nc1=1,c7=3,c14=5,c12=7\n\nUse either crossing order or explicit map, not both.",
    "preview_zoom": "2D preview controls only. Zoom + and Zoom - change the live preview scale without cropping the full diagram while unused margin remains. Rotate buttons update the rotate-degrees field, so saved 2D images use the same rotation. These controls do not change the DT mapping.",
}


# --------------------------------------------------------------------------- #
#  1. Parsing and combinatorial model
# --------------------------------------------------------------------------- #
def parse_dt(text):
    """Parse a string like 'DT: [(-8,-12,16),(-24,...)]'."""
    # Copy-pasting from documents often turns the ASCII hyphen-minus into a
    # look-alike Unicode dash (math minus, en/em dash, fullwidth minus) or a
    # non-ASCII space; normalize those so ast.literal_eval accepts the code.
    text = str(text)
    for ch in ("−", "–", "—", "‐", "‑", "－", "―"):
        text = text.replace(ch, "-")
    for ch in (" ", " ", " ", " ", " ", "　"):
        text = text.replace(ch, " ")
    m = re.search(r"\[.*\]", text.strip(), re.DOTALL)
    if not m:
        raise ValueError("Could not find a '[...]' list in the DT input.")
    try:
        raw = ast.literal_eval(m.group(0))
    except Exception as exc:
        raise ValueError("Could not parse the DT list: %s" % exc)

    if not isinstance(raw, (list, tuple)) or len(raw) == 0:
        raise ValueError("Empty or invalid DT code.")

    # Accept both multi-component syntax [(...), (...)] and a single-component
    # shorthand [4, 6, -2].
    if all(isinstance(x, int) for x in raw):
        raw = [tuple(raw)]

    comps = []
    for ci, comp in enumerate(raw, start=1):
        if not isinstance(comp, (list, tuple)):
            raise ValueError(
                "Component %d is not a list/tuple. Use syntax like "
                "DT: [(-8,-12,16), (...)]" % ci
            )
        clean = []
        for x in comp:
            if not isinstance(x, int):
                raise ValueError("DT entries must be integers; got %r." % (x,))
            if x == 0 or abs(x) % 2 != 0:
                raise ValueError(
                    "DT entries must be nonzero even integers; got %r." % x
                )
            clean.append(int(x))
        comps.append(tuple(clean))

    if not comps:
        raise ValueError("Empty DT code.")
    return comps


def build_model(comps, negative_even="over"):
    """
    Build the full combinatorial model from signed DT components.

    Returned keys include:
      comp_positions : list[list[int]] positions 1..2n per component
      comp_of        : dict position -> component index, zero-based
      nextpos        : dict position -> next position along its component
      prevpos        : dict position -> previous position along its component
      label_of       : callable position -> visible signed label
      pos_cross      : dict position -> internal crossing index, zero-based
      pos_role       : dict position -> 'o' for odd visit or 'e' for even visit
      over_at        : dict position -> True if that visit is the over-strand
      crossings      : list of {odd, even, even_signed}, ordered by odd labels
      twon           : total number of traversal labels
    """
    pos = 1
    odd_partner = {}
    even_sign = {}
    comp_positions = []
    comp_of = {}

    for ci, comp in enumerate(comps):
        cp = []
        for signed_even in comp:
            odd_pos = pos
            traversal_even_pos = pos + 1
            dt_even = abs(signed_even)

            odd_partner[odd_pos] = signed_even
            even_sign[dt_even] = 1 if signed_even > 0 else -1

            cp.append(odd_pos)
            cp.append(traversal_even_pos)
            comp_of[odd_pos] = ci
            comp_of[traversal_even_pos] = ci
            pos += 2
        comp_positions.append(cp)

    twon = pos - 1

    # Validation requested by the user: the absolute even DT labels must be
    # exactly 2,4,...,2n with no repeats.
    evens = sorted(even_sign.keys())
    expected = list(range(2, twon + 1, 2))
    if evens != expected:
        raise ValueError(
            "Invalid DT code: absolute even labels are %s but must be exactly "
            "%s, with no repeats." % (evens, expected)
        )

    def label_of(p):
        if p % 2 == 1:
            return p
        return even_sign[p] * p

    nextpos = {}
    prevpos = {}
    for cp in comp_positions:
        if len(cp) == 0:
            continue
        L = len(cp)
        for i, p in enumerate(cp):
            nextpos[p] = cp[(i + 1) % L]
            prevpos[p] = cp[(i - 1) % L]

    crossings = []
    pos_cross = {}
    pos_role = {}
    for k, odd_pos in enumerate(sorted(odd_partner)):
        signed = odd_partner[odd_pos]
        even_pos = abs(signed)
        crossings.append({"odd": odd_pos, "even": even_pos, "even_signed": signed})
        pos_cross[odd_pos] = k
        pos_role[odd_pos] = "o"
        pos_cross[even_pos] = k
        pos_role[even_pos] = "e"

    # Default convention requested by the user:
    #   negative even DT label => the even-labeled visit is over.
    # Opposite convention is available with --negative-even under.
    over_at = {}
    for c in crossings:
        neg = c["even_signed"] < 0
        if negative_even == "over":
            even_over = neg
        else:
            even_over = not neg
        over_at[c["even"]] = even_over
        over_at[c["odd"]] = not even_over

    return {
        "comp_positions": comp_positions,
        "comp_of": comp_of,
        "nextpos": nextpos,
        "prevpos": prevpos,
        "label_of": label_of,
        "pos_cross": pos_cross,
        "pos_role": pos_role,
        "over_at": over_at,
        "crossings": crossings,
        "twon": twon,
    }


# --------------------------------------------------------------------------- #
#  2. Crossing-ID remapping
# --------------------------------------------------------------------------- #
def default_crossing_ids(model):
    """Default displayed crossing IDs: c1, c2, ..., in increasing odd-label order."""
    return ["c%d" % (k + 1) for k in range(len(model["crossings"]))]


def _token_to_crossing_id(token):
    m = re.fullmatch(r"[cC]?(\d+)", str(token).strip())
    if not m:
        raise ValueError("Invalid crossing ID token %r. Use values like c7 or 7." % token)
    n = int(m.group(1))
    if n <= 0:
        raise ValueError("Crossing IDs must be positive; got %r." % token)
    return "c%d" % n


def _token_to_int(token):
    m = re.fullmatch(r"[cC]?(\d+)", str(token).strip())
    if not m:
        raise ValueError("Invalid numeric token %r." % token)
    return int(m.group(1))


def _validate_crossing_ids(ids, expected_count):
    if len(ids) != expected_count:
        raise ValueError(
            "Expected %d crossing IDs, but received %d." % (expected_count, len(ids))
        )
    seen = set()
    for cid in ids:
        if cid in seen:
            raise ValueError("Duplicate displayed crossing ID %s." % cid)
        seen.add(cid)
    return ids


def parse_crossing_order(text, model):
    """
    Parse a list of displayed crossing IDs corresponding to odd labels in order.

    Example for 14 crossings:
      c1 c7 c14 c12 c3 c6 c9 c5 c11 c13 c4 c2 c10 c8

    This means:
      odd 1  -> c1
      odd 3  -> c7
      odd 5  -> c14
      ...
      odd 27 -> c8

    The parser also accepts a pasted two-column/list block that contains both
    the crossing IDs and the odd labels, as long as the odd-label sequence is
    1,3,5,... in the same order.
    """
    if not text or not text.strip():
        return default_crossing_ids(model)

    ncross = len(model["crossings"])
    expected_odds = [c["odd"] for c in model["crossings"]]
    tokens = re.findall(r"[cC]?\d+", text)
    if not tokens:
        raise ValueError("No crossing IDs found in --crossing-order.")

    chosen = None
    if len(tokens) == ncross:
        chosen = tokens
    elif len(tokens) == 2 * ncross:
        first = tokens[:ncross]
        second = tokens[ncross:]
        if [_token_to_int(x) for x in second] == expected_odds:
            chosen = first
        elif [_token_to_int(x) for x in first] == expected_odds:
            chosen = second
        elif [_token_to_int(x) for x in tokens[1::2]] == expected_odds:
            chosen = tokens[0::2]
        else:
            raise ValueError(
                "Could not infer crossing IDs from the supplied --crossing-order. "
                "Provide exactly %d IDs, ordered by odd labels %s." %
                (ncross, expected_odds)
            )
    else:
        raise ValueError(
            "--crossing-order should contain exactly %d crossing IDs ordered by "
            "odd labels %s. Found %d numeric tokens." %
            (ncross, expected_odds, len(tokens))
        )

    ids = [_token_to_crossing_id(tok) for tok in chosen]
    return _validate_crossing_ids(ids, ncross)


def parse_crossing_map(text, model):
    """
    Parse explicit crossing-to-odd-label pairs, for example:
      c1=1, c7=3, c14=5

    This is an alternative to --crossing-order. The left side is the displayed
    crossing ID and the right side is the odd traversal label.
    """
    if not text or not text.strip():
        return None

    ncross = len(model["crossings"])
    odd_to_internal = {c["odd"]: k for k, c in enumerate(model["crossings"])}
    ids = [None] * ncross
    pairs = re.findall(r"([cC]?\d+)\s*[:=]\s*(\d+)", text)
    if len(pairs) != ncross:
        raise ValueError(
            "--crossing-map should contain exactly %d pairs such as c7=3. "
            "Found %d pair(s)." % (ncross, len(pairs))
        )

    for cid_token, odd_token in pairs:
        cid = _token_to_crossing_id(cid_token)
        odd = int(odd_token)
        if odd not in odd_to_internal:
            raise ValueError(
                "Odd label %d in --crossing-map is not one of %s." %
                (odd, sorted(odd_to_internal))
            )
        k = odd_to_internal[odd]
        if ids[k] is not None:
            raise ValueError("Odd label %d was assigned more than once." % odd)
        ids[k] = cid

    if any(x is None for x in ids):
        raise ValueError("Incomplete --crossing-map.")
    return _validate_crossing_ids(ids, ncross)


def resolve_crossing_ids(model, crossing_order=None, crossing_map=None):
    """Return displayed crossing IDs indexed by internal odd-label crossing order."""
    if crossing_order and crossing_order.strip() and crossing_map and crossing_map.strip():
        raise ValueError("Use either --crossing-order or --crossing-map, not both.")
    if crossing_map and crossing_map.strip():
        return parse_crossing_map(crossing_map, model)
    return parse_crossing_order(crossing_order or "", model)


# --------------------------------------------------------------------------- #
#  3. Planar graph with crossing gadgets
# --------------------------------------------------------------------------- #
def build_gadget_graph(model):
    """
    Every crossing k becomes a 4-cycle with corners, cyclically ordered as:
      in_o -- in_e -- out_o -- out_e -- in_o

    Thus the odd strand uses opposite corners in_o/out_o, and the even strand
    uses opposite corners in_e/out_e. Traversal segments p -> next(p) are added
    through degree-2 midpoint nodes ('seg', p).
    """
    G = nx.Graph()
    ncross = len(model["crossings"])

    for k in range(ncross):
        a = (k, "in_o")
        b = (k, "in_e")
        c = (k, "out_o")
        d = (k, "out_e")
        G.add_edge(a, b)
        G.add_edge(b, c)
        G.add_edge(c, d)
        G.add_edge(d, a)

    pos_cross = model["pos_cross"]
    pos_role = model["pos_role"]
    for cp in model["comp_positions"]:
        for p in cp:
            q = model["nextpos"][p]
            out_corner = (pos_cross[p], "out_" + pos_role[p])
            in_corner = (pos_cross[q], "in_" + pos_role[q])
            seg = ("seg", p)
            G.add_edge(out_corner, seg)
            G.add_edge(seg, in_corner)
    return G


def planar_faces(emb):
    seen = set()
    faces = []
    for u in emb:
        for v in emb[u]:
            if (u, v) in seen:
                continue
            faces.append(emb.traverse_face(u, v, mark_half_edges=seen))
    return faces


TUTTE_SHAPES = ("circle", "ellipse", "rectangle", "rounded-rectangle")
# Shapes offered in the GUI / CLI.  'circle' and 'rectangle' are dropped as
# redundant (an ellipse at aspect 1 is a circle; a rounded rectangle at corner
# radius 0 is a sharp rectangle).  The full TUTTE_SHAPES set is still accepted by
# the geometry / normalization so older saved sessions keep working.
TUTTE_SHAPE_CHOICES = ("ellipse", "rounded-rectangle")


def _normalize_tutte_shape(shape):
    s = str(shape or "circle").strip().lower().replace("_", "-")
    if s in ("round-rectangle", "roundrect", "rounded-rect"):
        s = "rounded-rectangle"
    if s in ("rect",):
        s = "rectangle"
    return s if s in TUTTE_SHAPES else "circle"


def _resample_polyline_by_arclength_closed_2d(poly, n):
    """Return n points evenly spaced by arc length around a closed 2D polyline."""
    poly = np.asarray(poly, float)
    n = int(n)
    if len(poly) == 0 or n <= 0:
        return np.zeros((n, 2), float)
    seg = np.linalg.norm(np.roll(poly, -1, axis=0) - poly, axis=1)
    total = float(np.sum(seg))
    if total <= 1.0e-12:
        return np.repeat(poly[:1], n, axis=0)
    cum = np.concatenate(([0.0], np.cumsum(seg)))
    out = np.zeros((n, 2), float)
    for i in range(n):
        target = total * i / float(n)
        j = int(np.searchsorted(cum, target, side="right") - 1)
        j = max(0, min(j, len(seg) - 1))
        t = 0.0 if seg[j] <= 1.0e-12 else (target - cum[j]) / seg[j]
        j2 = (j + 1) % len(poly)
        out[i] = (1.0 - t) * poly[j] + t * poly[j2]
    return out


def _rounded_rectangle_dense(half_w, half_h, radius, per_side=200):
    """Dense CCW boundary polyline of a (possibly rounded) rectangle centered at origin."""
    W = float(half_w)
    H = float(half_h)
    r = max(0.0, min(float(radius), min(W, H)))
    pts = []
    # Right edge (bottom -> top), top-right arc, top edge, top-left arc,
    # left edge, bottom-left arc, bottom edge, bottom-right arc.  CCW order.
    pts.append((W, -(H - r)))
    pts.append((W, (H - r)))
    if r > 0:
        for a in np.linspace(0.0, np.pi / 2, per_side // 4 + 2)[1:]:
            pts.append((W - r + r * np.cos(a), H - r + r * np.sin(a)))
    pts.append((-(W - r), H))
    if r > 0:
        for a in np.linspace(np.pi / 2, np.pi, per_side // 4 + 2)[1:]:
            pts.append((-(W - r) + r * np.cos(a), H - r + r * np.sin(a)))
    pts.append((-W, -(H - r)))
    if r > 0:
        for a in np.linspace(np.pi, 1.5 * np.pi, per_side // 4 + 2)[1:]:
            pts.append((-(W - r) + r * np.cos(a), -(H - r) + r * np.sin(a)))
    pts.append(((W - r), -H))
    if r > 0:
        for a in np.linspace(1.5 * np.pi, 2.0 * np.pi, per_side // 4 + 2)[1:]:
            pts.append((W - r + r * np.cos(a), -(H - r) + r * np.sin(a)))
    return np.asarray(pts, float)


def _boundary_shape_points(n, shape="circle", aspect=1.0, corner_radius=0.0):
    """
    Return n points placed in CCW cyclic order on a convex boundary shape,
    centered at the origin.

    For ``circle``/``ellipse`` the points are evenly spaced in angle; for the
    rectangles they are evenly spaced by arc length around the perimeter.  All
    shapes are convex, so the fixed boundary keeps the Tutte solve valid.  With
    ``shape='circle'`` and ``aspect=1`` this reproduces the classic unit-circle
    boundary exactly.
    """
    n = int(n)
    if n <= 0:
        return np.zeros((0, 2), float)
    shape = _normalize_tutte_shape(shape)
    aspect = float(aspect)
    if not np.isfinite(aspect) or aspect <= 0.0:
        aspect = 1.0

    if shape in ("circle", "ellipse"):
        # a/b = aspect, a*b = 1 so a round diagram (aspect=1) stays a unit circle.
        a = math.sqrt(aspect)
        b = 1.0 / math.sqrt(aspect)
        th = 2.0 * np.pi * np.arange(n) / float(max(n, 1))
        return np.column_stack([a * np.cos(th), b * np.sin(th)])

    # Rectangle / rounded rectangle: half-extents with W/H = aspect, W*H = 1.
    W = math.sqrt(aspect)
    H = 1.0 / math.sqrt(aspect)
    radius = 0.0
    if shape == "rounded-rectangle":
        frac = max(0.0, min(float(corner_radius), 0.999))
        radius = frac * min(W, H)
    dense = _rounded_rectangle_dense(W, H, radius)
    return _resample_polyline_by_arclength_closed_2d(dense, n)


def _boundary_depth(G, boundary):
    """Multi-source BFS graph distance of every node from the boundary node set."""
    from collections import deque

    depth = {}
    dq = deque()
    for node in boundary:
        depth[node] = 0
        dq.append(node)
    while dq:
        u = dq.popleft()
        for v in G[u]:
            if v not in depth:
                depth[v] = depth[u] + 1
                dq.append(v)
    return depth


def _pca_axes_2d(coords):
    """Return (major_axis, minor_axis, eigenvalues_desc) for a 2D point cloud."""
    pts = np.asarray(coords, float)
    if len(pts) < 2:
        return np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([1.0, 1.0])
    q = pts - pts.mean(axis=0)
    cov = q.T @ q
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]
    return evecs[:, 0], evecs[:, 1], np.maximum(evals, 0.0)


def _pca_axis_ratio_2d(coords):
    """Elongation ratio (major/minor standard deviation) of a 2D point cloud."""
    _major, _minor, evals = _pca_axes_2d(coords)
    lo = float(evals[1])
    hi = float(evals[0])
    if lo <= 1.0e-12:
        return 1.0
    return math.sqrt(hi / lo)


def auto_orient_positions(P, mode="off"):
    """
    Rotate a 2D position map so its principal (PCA) axis is aligned with an axis.

    mode='horizontal' aligns the elongation axis with the x-axis; 'vertical'
    with the y-axis; 'off'/None leaves the map unchanged.
    """
    m = str(mode or "off").strip().lower()
    if not P or m in ("off", "none", ""):
        return P
    coords = np.array(list(P.values()), float)
    major, _minor, _evals = _pca_axes_2d(coords)
    angle = -math.atan2(float(major[1]), float(major[0]))
    if m == "vertical":
        angle += math.pi / 2.0
    c = math.cos(angle)
    s = math.sin(angle)
    R = np.array([[c, -s], [s, c]])
    center = coords.mean(axis=0)
    return {k: R.dot(np.asarray(xy, float) - center) + center for k, xy in P.items()}


def rotate_positions_2d(P, degrees):
    """Rotate a 2D position map about its centroid by a fixed angle in degrees."""
    if not P:
        return P
    theta = math.radians(float(degrees or 0.0))
    if abs(theta) < 1.0e-12:
        return P
    c = math.cos(theta)
    s = math.sin(theta)
    R = np.array([[c, -s], [s, c]])
    coords = np.array(list(P.values()), float)
    center = coords.mean(axis=0)
    return {k: R.dot(np.asarray(xy, float) - center) + center for k, xy in P.items()}


def orient_positions_2d(P, auto_orient=True, orient_degrees=0.0):
    """
    Orient a shaped 2D layout relative to its principal (PCA) elongation axis.

    The diagram is always first aligned so its PCA major axis is horizontal
    (along the view x-axis).  With auto_orient=True that PCA alignment is the
    final orientation.  With auto_orient=False an additional rotation of
    orient_degrees is applied, so the manual angle is measured relative to the
    PCA axis (0 = along the PCA axis / horizontal).  This works for every
    boundary shape, since it is a post-layout rotation of the whole diagram.
    """
    P2 = auto_orient_positions(P, mode="horizontal")
    if auto_orient:
        return P2
    return rotate_positions_2d(P2, orient_degrees)


def _expand_about_center_of_mass(pos, G, outer_set, com_expand):
    """
    Radially expand interior nodes outward from the density-weighted crossing
    center of mass, tapered to zero at the boundary so pinned boundary nodes do
    not move.

    ``pos`` maps every node to a 2D coordinate.  The center of mass is the mean
    of the crossing-corner node positions -- a region with more (crowded)
    crossings contributes more corners, so the COM sits in the crowded area.
    Deep interior nodes (far from the boundary in graph distance) are expanded
    the most, giving crowded crossings more room and making that region the
    focal center.  Strength 0 leaves ``pos`` unchanged.
    """
    strength = float(com_expand or 0.0)
    if strength <= 0.0:
        return pos
    corners = [n for n in pos if not (isinstance(n, tuple) and len(n) == 2 and n[0] == "seg")]
    if not corners:
        return pos
    com = np.mean([pos[n] for n in corners], axis=0)
    depth = _boundary_depth(G, outer_set)
    max_depth = max(depth.values()) if depth else 0
    if max_depth <= 0:
        return pos
    out = {}
    for n, xy in pos.items():
        f = float(depth.get(n, 0)) / float(max_depth)   # 0 at boundary, 1 deepest
        out[n] = xy + strength * f * (xy - com)
    return out


def tutte_layout_connected(G, emb, shape="circle", aspect=1.0, corner_radius=0.0,
                           decompress=0.0, auto_aspect=False, com_expand=0.0):
    """
    Barycentric/Tutte embedding for one connected planar graph.

    The outer face is pinned to a convex boundary ``shape`` (circle, ellipse,
    rectangle, or rounded rectangle).  ``decompress`` >= 0 applies boundary-depth
    edge weights in a convex-combination solve so interior structure is pushed
    outward instead of being over-compressed toward the center; 0 reproduces the
    classic uniform Tutte solve.  ``auto_aspect`` derives the boundary aspect
    ratio from the diagram's own elongation (measured from a circular solve).
    ``com_expand`` >= 0 is a separate post-solve radial expansion about the
    crossing center of mass (see _expand_about_center_of_mass).

    With shape='circle', aspect=1, decompress=0, com_expand=0 this reproduces the
    previous behavior exactly.
    """
    faces = planar_faces(emb)
    if not faces:
        return {n: np.array([0.0, 0.0]) for n in G.nodes()}

    outer = max(faces, key=len)
    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    outer_set = set(outer)
    m = len(outer)

    beta = float(decompress or 0.0)
    depth = _boundary_depth(G, outer_set) if beta > 0.0 else None

    # The left-hand side depends only on the interior weighting and which nodes
    # are pinned; the boundary shape only changes the right-hand side, so we can
    # factor A once and reuse it for both a probe (circle) and the final solve.
    A = np.zeros((N, N))
    for n in nodes:
        i = idx[n]
        if n in outer_set:
            A[i, i] = 1.0
        else:
            nbrs = list(G[n])
            if beta > 0.0 and depth is not None:
                di = depth.get(n, 0)
                w = np.array([math.exp(beta * (di - depth.get(nb, 0))) for nb in nbrs], float)
                A[i, i] = float(np.sum(w))
                for nb, wij in zip(nbrs, w):
                    A[i, idx[nb]] -= wij
            else:
                A[i, i] = float(len(nbrs))
                for nb in nbrs:
                    A[i, idx[nb]] -= 1.0

    def _solve_with(bpts):
        bx = np.zeros(N)
        by = np.zeros(N)
        for j, node in enumerate(outer):
            bx[idx[node]], by[idx[node]] = bpts[j]
        X = np.linalg.solve(A, bx)
        Y = np.linalg.solve(A, by)
        return X, Y

    shp = _normalize_tutte_shape(shape)
    if auto_aspect and shp in ("ellipse", "rectangle", "rounded-rectangle"):
        Xc, Yc = _solve_with(_boundary_shape_points(m, "circle", 1.0, 0.0))
        ratio = _pca_axis_ratio_2d(np.column_stack([Xc, Yc]))
        aspect = float(np.clip(ratio, 1.0, 6.0))

    bpts = _boundary_shape_points(m, shp, aspect, corner_radius)
    X, Y = _solve_with(bpts)
    pos = {n: np.array([X[idx[n]], Y[idx[n]]]) for n in nodes}
    if float(com_expand or 0.0) > 0.0:
        pos = _expand_about_center_of_mass(pos, G, outer_set, com_expand)
    return pos


def _shape_radius_interpolator(dense):
    """
    Build a periodic polar-radius lookup r(angle) for a convex outline centered
    at the origin.  ``dense`` is an Nx2 CCW polyline; returns a callable mapping
    an angle (radians) to the boundary radius along that ray.
    """
    dense = np.asarray(dense, float)
    th = np.arctan2(dense[:, 1], dense[:, 0])
    r = np.hypot(dense[:, 0], dense[:, 1])
    order = np.argsort(th)
    th_s = th[order]
    r_s = r[order]
    # Extend one period on each side so np.interp wraps smoothly across +-pi.
    th_ext = np.concatenate([th_s - 2.0 * np.pi, th_s, th_s + 2.0 * np.pi])
    r_ext = np.concatenate([r_s, r_s, r_s])

    def rad(angle):
        a = (float(angle) + np.pi) % (2.0 * np.pi) - np.pi
        return float(np.interp(a, th_ext, r_ext))

    return rad


def shaped_tutte_layout(G, emb, shape="ellipse", aspect=1.0, corner_radius=0.0,
                        decompress=0.0, auto_aspect=False, com_expand=0.0,
                        orient_degrees=0.0, auto_orient=True, meta_out=None):
    """
    Shaped-tutte layout in which the boundary shape's aspect (long) axis is tilted
    relative to the diagram's *intrinsic* elongation axis.

    The intrinsic axis is the principal (PCA) direction of a plain circular Tutte
    solve -- it depends only on the graph, not on the imposed shape.  The boundary
    is then pinned so its long axis makes an angle ``orient_degrees`` with that
    intrinsic PCA axis (0 = shape long axis along the PCA axis).  This genuinely
    re-stretches the diagram in a new direction rather than merely spinning the
    finished picture.

    With ``auto_orient`` the final frame is rotated so the intrinsic PCA axis is
    horizontal (the shape outline then appears tilted by ``orient_degrees``);
    otherwise the shape long axis is left horizontal.

    When ``meta_out`` is a dict it is filled with guide geometry (in the returned
    coordinate frame): ``boundary`` (Nx2 shape outline), ``shape_axis`` and
    ``pca_axis`` (each a 2x2 pair of endpoints), for optional overlay drawing.
    """
    faces = planar_faces(emb)
    if not faces:
        return {n: np.array([0.0, 0.0]) for n in G.nodes()}

    outer = max(faces, key=len)
    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    outer_set = set(outer)
    m = len(outer)

    beta = float(decompress or 0.0)
    depth = _boundary_depth(G, outer_set) if beta > 0.0 else None

    A = np.zeros((N, N))
    for n in nodes:
        i = idx[n]
        if n in outer_set:
            A[i, i] = 1.0
        else:
            nbrs = list(G[n])
            if beta > 0.0 and depth is not None:
                di = depth.get(n, 0)
                w = np.array([math.exp(beta * (di - depth.get(nb, 0))) for nb in nbrs], float)
                A[i, i] = float(np.sum(w))
                for nb, wij in zip(nbrs, w):
                    A[i, idx[nb]] -= wij
            else:
                A[i, i] = float(len(nbrs))
                for nb in nbrs:
                    A[i, idx[nb]] -= 1.0

    def _solve_with(bpts):
        bx = np.zeros(N)
        by = np.zeros(N)
        for j, node in enumerate(outer):
            bx[idx[node]], by[idx[node]] = bpts[j]
        X = np.linalg.solve(A, bx)
        Y = np.linalg.solve(A, by)
        return X, Y

    # 1) Circular probe -> intrinsic PCA elongation axis (graph-only property).
    Xc, Yc = _solve_with(_boundary_shape_points(m, "circle", 1.0, 0.0))
    coords_c = np.column_stack([Xc, Yc])
    major, _minor, _evals = _pca_axes_2d(coords_c)
    alpha = math.atan2(float(major[1]), float(major[0]))

    shp = _normalize_tutte_shape(shape)
    asp = float(aspect)
    if not np.isfinite(asp) or asp <= 0.0:
        asp = 1.0
    if auto_aspect and shp in ("ellipse", "rectangle", "rounded-rectangle"):
        asp = float(np.clip(_pca_axis_ratio_2d(coords_c), 1.0, 6.0))

    orient = math.radians(float(orient_degrees or 0.0))

    # 2) Pin the outer face onto the shape by polar angle, offsetting the whole
    #    angular assignment by -(alpha + orient).  This places the intrinsic
    #    PCA-tip node at boundary polar angle -orient, so the shape long axis
    #    (polar 0) ends up 'orient' away from the intrinsic PCA axis.
    dense = _boundary_shape_points(720, shp, asp, corner_radius)
    rad = _shape_radius_interpolator(dense)
    bpts = np.zeros((m, 2))
    for i in range(m):
        theta_i = 2.0 * np.pi * i / float(max(m, 1))
        b = theta_i - alpha - orient
        rr = rad(b)
        bpts[i] = [rr * math.cos(b), rr * math.sin(b)]

    X, Y = _solve_with(bpts)
    pos = {n: np.array([X[idx[n]], Y[idx[n]]]) for n in nodes}
    if float(com_expand or 0.0) > 0.0:
        pos = _expand_about_center_of_mass(pos, G, outer_set, com_expand)

    # 3) Final framing rotation about the content centroid.
    #    In the canonical frame the shape long axis is along +x and the intrinsic
    #    PCA axis is at angle -orient.  auto_orient rotates the PCA axis to
    #    horizontal (shape then tilts by +orient); otherwise leave as is.
    phi = orient if auto_orient else 0.0
    c = math.cos(phi)
    s = math.sin(phi)
    R = np.array([[c, -s], [s, c]])
    ctr = np.mean(list(pos.values()), axis=0) if pos else np.zeros(2)
    pos = {n: R.dot(xy - ctr) + ctr for n, xy in pos.items()}

    if meta_out is not None:
        def _frame(arr):
            return np.array([R.dot(np.asarray(p, float) - ctr) + ctr for p in arr])

        L_shape = rad(0.0)
        ang_pca = -orient
        L_pca = rad(ang_pca)
        shape_axis = np.array([[-L_shape, 0.0], [L_shape, 0.0]])
        pca_axis = np.array([
            [-L_pca * math.cos(ang_pca), -L_pca * math.sin(ang_pca)],
            [L_pca * math.cos(ang_pca), L_pca * math.sin(ang_pca)],
        ])
        meta_out.clear()
        meta_out["kind"] = "shaped"
        meta_out["aspect_value"] = float(asp)
        meta_out["boundary"] = _frame(dense)
        meta_out["shape_axis"] = _frame(shape_axis)
        meta_out["pca_axis"] = _frame(pca_axis)

    return pos


def _face_angular_spread(face, coords, idx, center):
    """Angular coverage (0..2pi) of a face's vertices about ``center``.

    Returns 2*pi minus the largest angular gap between consecutive vertex angles;
    a value near 2*pi means the face encircles the center (good annulus hole).
    """
    angs = []
    for v in face:
        p = np.asarray(coords[idx[v]], float) - center
        angs.append(math.atan2(float(p[1]), float(p[0])))
    if len(angs) < 2:
        return 0.0
    angs = sorted(angs)
    gaps = [angs[i + 1] - angs[i] for i in range(len(angs) - 1)]
    gaps.append(2.0 * math.pi - (angs[-1] - angs[0]))
    return 2.0 * math.pi - max(gaps)


def _pick_hole_face(faces, outer_set, coords, idx, G, center):
    """
    Choose an inner 'hole' face for the holed-tutte layout.

    Candidates are faces vertex-disjoint from the outer face.  We prefer faces
    that both encircle the diagram center (large angular spread in the probe
    layout, so the ring is not twisted) and sit deep/central in the graph.
    Returns the chosen face (list of nodes) or None.
    """
    depth = _boundary_depth(G, outer_set)
    best = None
    best_score = None
    for f in faces:
        if len(f) < 3:
            continue
        fs = set(f)
        if fs & outer_set:
            continue  # must be vertex-disjoint from the outer boundary
        spread = _face_angular_spread(f, coords, idx, center)
        mean_depth = float(np.mean([depth.get(v, 0) for v in f]))
        # Encircling first (rounded so near-ties defer to depth), then depth, size.
        score = (round(spread, 2), mean_depth, len(f))
        if best_score is None or score > best_score:
            best_score = score
            best = f
    return best


# Cache of the (expensive, otherwise non-deterministic) 3D-torus longitudinal
# projection used by holed-tutte, keyed on the crossing graph.  This keeps the
# holed layout stable while the user tweaks parameters (rotate, swap, hole ratio,
# ring tilt, ...): the 3D layout is only recomputed when the graph itself changes
# (i.e. a new DT code), not on every small adjustment.
_HOLED_TORUS_CACHE = {}
_HOLED_TORUS_CACHE_ORDER = []
# Seed for the 3D layout.  It is fixed (so the diagram is stable while tweaking
# parameters) but can be bumped to force a fresh recomputation -- see the GUI's
# "Redraw 2D" button, which clears the cache and rerolls this.
_HOLED_TORUS_SEED = [1]


def clear_holed_cache(reroll=True):
    """Clear the cached 3D-torus layouts so holed-tutte recomputes from scratch."""
    _HOLED_TORUS_CACHE.clear()
    _HOLED_TORUS_CACHE_ORDER.clear()
    if reroll:
        _HOLED_TORUS_SEED[0] += 1


def _holed_torus_coords(G):
    """
    Deterministic 3D-torus longitudinal projection of the crossing graph, cached
    per graph.  A 3D Kamada-Kawai layout (seeded so it is reproducible) is
    projected onto the plane of its two largest principal axes; the azimuth there
    is the torus longitudinal angle.  Returns a dict node -> (x, y), or None if
    kamada is unavailable (e.g. SciPy missing).
    """
    try:
        sig = (G.number_of_nodes(),
               frozenset(frozenset((u, v)) for u, v in G.edges()))
    except Exception:
        sig = None
    if sig is not None and sig in _HOLED_TORUS_CACHE:
        return _HOLED_TORUS_CACHE[sig]

    nodes = list(G.nodes())
    coords = None
    try:
        # A seeded initial layout makes the 3D kamada reproducible (networkx uses a
        # random init for dim>=3, which is what made the diagram jump around).
        init = nx.random_layout(G, dim=3, seed=_HOLED_TORUS_SEED[0])
        kpos3 = nx.kamada_kawai_layout(G, pos=init, dim=3)
        c3 = np.array([kpos3[n] for n in nodes], float)
        q3 = c3 - c3.mean(axis=0)
        _ev, _evec = np.linalg.eigh(q3.T @ q3)
        # [:, 0] = torus axis (thinnest spread); project onto the two largest.
        e1 = _evec[:, 2].copy()
        e2 = _evec[:, 1].copy()
        # Deterministic sign so the projected annulus has a stable orientation.
        if e1[int(np.argmax(np.abs(e1)))] < 0.0:
            e1 = -e1
        if e2[int(np.argmax(np.abs(e2)))] < 0.0:
            e2 = -e2
        proj = np.column_stack([q3 @ e1, q3 @ e2])
        coords = {n: proj[i] for i, n in enumerate(nodes)}
    except Exception:
        coords = None

    if sig is not None:
        _HOLED_TORUS_CACHE[sig] = coords
        _HOLED_TORUS_CACHE_ORDER.append(sig)
        if len(_HOLED_TORUS_CACHE_ORDER) > 16:
            _HOLED_TORUS_CACHE.pop(_HOLED_TORUS_CACHE_ORDER.pop(0), None)
    return coords


def holed_tutte_layout(G, emb, shape="circle", aspect=1.0, corner_radius=0.0,
                       hole_ratio=0.4, swap=False, auto_aspect=False,
                       orient_degrees=0.0, auto_orient=True, ring_tilt=90.0,
                       invert_ring=False, ring_equalize=0.0, meta_out=None):
    """
    'Holed' Tutte layout: two boundary cycles of the planar embedding are pinned,
    one to an outer outline and one to an inner (hole) outline of a holed shape
    (annulus / elliptical or rectangular ring).  The remaining vertices solve the
    harmonic (barycentric) system -- a discrete harmonic map onto the ring -- so
    the diagram wraps around a central hole along its natural closed 'principal
    curved axis'.

    The outer face is the largest planar face; the hole is an auto-picked central,
    encircling face (see _pick_hole_face).  ``swap`` exchanges which cycle is
    pinned to the outer vs. inner outline.  ``hole_ratio`` in (0,1) sets the inner
    outline size as a fraction of the outer.  ``shape``/``aspect``/
    ``corner_radius``/``orient_degrees`` behave as in shaped_tutte_layout and apply
    to both outlines.

    ``meta_out``, when given, is filled with ``boundary_outer``, ``boundary_inner``
    and ``medial`` (the mid-ring closed curved axis) in the returned frame.

    The angular distribution of the crossings around the ring (the 'closed
    principal curve') is taken from a Kamada-Kawai layout of the crossing graph,
    which spreads the components evenly around the diagram, so the ring is balanced
    instead of lopsided.  ``ring_tilt`` then views the resulting ring as the wall of
    a bucket: 90 deg looks straight down (flat top-view annulus), lower values tilt
    the bucket so its wall opens toward a side view.  The tilt lifts each crossing
    onto the 3D bucket wall and projects it, so it is intrinsic to the mapping (not
    a post-hoc rotation of the finished picture).
    """
    faces = planar_faces(emb)
    if not faces:
        return {n: np.array([0.0, 0.0]) for n in G.nodes()}
    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)

    outer = max(faces, key=len)
    outer_set = set(outer)
    m_out = len(outer)

    def _build_system(pinned):
        A = np.zeros((N, N))
        for n in nodes:
            i = idx[n]
            if n in pinned:
                A[i, i] = 1.0
            else:
                nbrs = list(G[n])
                A[i, i] = float(len(nbrs))
                for nb in nbrs:
                    A[i, idx[nb]] -= 1.0
        return A

    def _solve(A, bmap):
        bx = np.zeros(N)
        by = np.zeros(N)
        for node, (px, py) in bmap.items():
            bx[idx[node]] = px
            by[idx[node]] = py
        return np.linalg.solve(A, bx), np.linalg.solve(A, by)

    # 1) 'Closed principal curve' from a 3D torus layout.  A 3D Kamada-Kawai layout
    #    of the crossing graph naturally arranges a symmetric multi-component link
    #    on a torus; its azimuth about the torus axis (the longitudinal 'long way
    #    around') is spread far more evenly than a 2D kamada, which is what makes the
    #    wreath symmetric.  We take that torus and view it top-down: the projection
    #    onto the plane of its two largest principal axes is a clean annulus whose
    #    angle is the longitudinal coordinate and whose radius is the meridional
    #    coordinate (outer equator -> outer rim, inner equator -> inner rim).  When
    #    a link is not annular this way (or SciPy is missing) we fall back to the
    #    circular-Tutte probe, which forces a clean annulus.
    def _mean_radius(coords, ctr, face):
        return float(np.mean([np.linalg.norm(coords[idx[v]] - ctr) for v in face]))

    coords_p = None
    outer = None
    hole = None
    use_kamada = False

    coords_k = None
    _torus = _holed_torus_coords(G)  # cached per graph; only recomputed on DT change
    if _torus is not None:
        coords_k = np.array([_torus[n] for n in nodes], float)

    if coords_k is not None:
        ctr_k = coords_k.mean(axis=0)
        encircling = [
            f for f in faces
            if len(f) >= 3
            and _face_angular_spread(f, coords_k, idx, ctr_k) > 1.3 * np.pi
        ]
        if len(encircling) >= 2:
            outer_c = max(encircling, key=lambda f: _mean_radius(coords_k, ctr_k, f))
            outer_cs = set(outer_c)
            hole_c = None
            best_r = None
            for f in encircling:
                if set(f) & outer_cs:
                    continue
                r = _mean_radius(coords_k, ctr_k, f)
                if best_r is None or r < best_r:
                    best_r = r
                    hole_c = f
            if hole_c is not None:
                coords_p, outer, hole = coords_k, outer_c, hole_c
                use_kamada = True

    if coords_p is None:
        # Circular-probe fallback: pin the largest face to a unit circle (this
        # forces an annulus even when the link is not naturally wreath-like), and
        # pick a central encircling hole face from that probe layout.
        outer = max(faces, key=len)
        outer_set = set(outer)
        A_probe = _build_system(outer_set)
        probe_map = {}
        for j, node in enumerate(outer):
            th = 2.0 * np.pi * j / float(max(m_out, 1))
            probe_map[node] = (math.cos(th), math.sin(th))
        Xp, Yp = _solve(A_probe, probe_map)
        coords_p = np.column_stack([Xp, Yp])
        ctr_probe = coords_p.mean(axis=0)
        outer = max(faces, key=len)
        outer_set = set(outer)
        hole = _pick_hole_face(faces, outer_set, coords_p, idx, G, ctr_probe)
        if hole is None:
            return shaped_tutte_layout(
                G, emb, shape=shape, aspect=aspect, corner_radius=corner_radius,
                auto_aspect=auto_aspect, orient_degrees=orient_degrees,
                auto_orient=auto_orient, meta_out=meta_out,
            )

    shp = _normalize_tutte_shape(shape)
    if auto_aspect and shp in ("ellipse", "rectangle", "rounded-rectangle"):
        aspect = float(np.clip(_pca_axis_ratio_2d(coords_p), 1.0, 6.0))
    orient = math.radians(float(orient_degrees or 0.0))
    hr = float(hole_ratio)
    if not np.isfinite(hr) or hr <= 0.0:
        hr = 0.4
    hr = min(max(hr, 0.05), 0.95)

    dense = _boundary_shape_points(720, shp, float(aspect), corner_radius)
    rad = _shape_radius_interpolator(dense)

    def _cycle_angles(cycle, coords, ctr, alpha):
        # Angular positions of a cycle's vertices from ``coords`` (monotonically
        # unwrapped so the pinned polygon does not self-intersect).  The raw
        # spacing is kept -- forcing even angles twists the interior into a mess,
        # because the two rings then lose their natural radial correspondence.
        n = len(cycle)
        if n == 0:
            return []
        raw = []
        for v in cycle:
            p = coords[idx[v]] - ctr
            raw.append((math.atan2(float(p[1]), float(p[0])) - alpha - orient) % (2.0 * np.pi))
        start = int(np.argmin(raw))
        order = [(start + k) % n for k in range(n)]
        unwrapped = {order[0]: raw[order[0]]}
        acc = raw[order[0]]
        for k in range(1, n):
            oi = order[k]
            a = raw[oi]
            while a <= acc + 1.0e-6:
                a += 2.0 * np.pi
            acc = a
            unwrapped[oi] = a
        return [unwrapped[i] for i in range(n)]

    def _pick_ring(coords, ctr):
        enc = [
            f for f in faces
            if len(f) >= 3 and _face_angular_spread(f, coords, idx, ctr) > 1.3 * np.pi
        ]
        if len(enc) < 2:
            return None
        o = max(enc, key=lambda f: _mean_radius(coords, ctr, f))
        os = set(o)
        h = None
        br = None
        for f in enc:
            if set(f) & os:
                continue
            r = _mean_radius(coords, ctr, f)
            if br is None or r < br:
                br = r
                h = f
        return (o, h) if h is not None else None

    A_cache = {}

    def _one_round(coords, outer_face, hole_face):
        ctr = coords.mean(axis=0)
        maj, _mn, _e = _pca_axes_2d(coords)
        a0 = math.atan2(float(maj[1]), float(maj[0]))
        # 'swap' just exchanges which face is pinned to the outer vs inner outline.
        face_out = hole_face if swap else outer_face
        face_in = outer_face if swap else hole_face
        key = frozenset(set(outer_face) | set(hole_face))
        A = A_cache.get(key)
        if A is None:
            A = _build_system(set(outer_face) | set(hole_face))
            A_cache[key] = A
        pin = {}
        for nd, a in zip(face_out, _cycle_angles(face_out, coords, ctr, a0)):
            r = rad(a)
            pin[nd] = (r * math.cos(a), r * math.sin(a))
        for nd, a in zip(face_in, _cycle_angles(face_in, coords, ctr, a0)):
            r = hr * rad(a)
            pin[nd] = (r * math.cos(a), r * math.sin(a))
        Xf, Yf = _solve(A, pin)
        return np.column_stack([Xf, Yf])

    # 2) Pin the outer + hole cycles onto the outlines at their (balanced) angles
    #    and Tutte-solve the interior.  For the kamada ring we then re-run the same
    #    step a few times, each round re-deriving the ring vertices' angles from the
    #    previous result: this fixed-point refinement gradually evens out the
    #    residual kamada jitter and improves the diagram's symmetry.  The two ring
    #    faces are chosen ONCE (above) and held fixed across rounds -- re-picking
    #    them by geometry each round would silently undo the 'swap' option, because
    #    after a swap the inner face has become the geometrically-outer one.
    current = coords_p
    n_refine = 3 if use_kamada else 0
    for _it in range(n_refine + 1):
        current = _one_round(current, outer, hole)
    pos = {n: current[idx[n]] for n in nodes}

    # 2a) V5.0 'ring equalize': the harmonic solve compresses interior detail
    #     toward the pinned rims on large rings.  Histogram-equalize the free
    #     nodes' radial fractions across the ring width (angles unchanged) so the
    #     radial distribution becomes uniform; blend by the requested strength.
    s_eq = min(1.0, max(0.0, float(ring_equalize or 0.0)))
    if s_eq > 0.0:
        pinned_ring = set(outer) | set(hole)
        free_nodes = [n for n in nodes if n not in pinned_ring]
        if len(free_nodes) >= 2:
            tvals = {}
            for n in free_nodes:
                px, py = pos[n]
                a = math.atan2(py, px)
                r = math.hypot(px, py)
                ro = rad(a)
                ri = hr * ro
                t = (r - ri) / max(ro - ri, 1.0e-12)
                tvals[n] = min(1.0, max(0.0, t))
            order_eq = sorted(free_nodes, key=lambda n: tvals[n])
            m_eq = float(len(order_eq))
            for rank, n in enumerate(order_eq):
                t_new = (1.0 - s_eq) * tvals[n] + s_eq * ((rank + 1.0) / (m_eq + 1.0))
                px, py = pos[n]
                a = math.atan2(py, px)
                r = math.hypot(px, py)
                ro = rad(a)
                ri = hr * ro
                r_new = ri + t_new * (ro - ri)
                sc = r_new / max(r, 1.0e-12)
                pos[n] = np.array([px * sc, py * sc])

    # 2b) Optional 'invert ring (inside-out)': reflect every crossing's radius
    #     about the ring's mid-line so the inner boundary ends up outside and the
    #     outer boundary inside, turning the same diagram inside-out (the medial
    #     curve is unchanged).  Unlike 'swap inner/outer face' (which re-solves with
    #     the two boundary faces exchanged), this keeps the layout and just flips it
    #     radially.
    def _invert_pt(px, py):
        r = math.hypot(px, py)
        if r <= 1.0e-12:
            return px, py
        a = math.atan2(py, px)
        r2 = (1.0 + hr) * rad(a) - r   # reflect about the medial radius
        s = r2 / r
        return px * s, py * s

    if invert_ring:
        pos = {n: np.array(_invert_pt(xy[0], xy[1])) for n, xy in pos.items()}

    # 3) 'Ring tilt' as a bucket view.  The flat annulus is the top-down view of a
    #    bucket whose wall carries the crossings; the radial distance across the
    #    wall is the height up the wall.  Lift each crossing onto that 3D wall and
    #    view the bucket from an angle: ring_tilt=90 looks straight down (the flat
    #    annulus), lower tilts rotate the bucket about the horizontal axis so the
    #    wall opens toward a side view.  This redistributes the crossings along the
    #    tilted principal curve rather than rigidly rotating the finished diagram.
    beta = math.radians(90.0) - math.radians(float(ring_tilt if ring_tilt is not None else 90.0))
    cb = math.cos(beta)
    sb = math.sin(beta)

    def _bucket(px, py):
        if abs(sb) <= 1.0e-9:
            return px, py
        r = math.hypot(px, py)
        if r <= 1.0e-12:
            return px, py
        a = math.atan2(py, px)
        r_bot = hr * rad(a)          # inner rim radius along this ray
        h = r - r_bot                # height up the bucket wall (0 at inner rim)
        return px, py * cb - h * sb

    if abs(sb) > 1.0e-9:
        pos = {n: np.array(_bucket(xy[0], xy[1])) for n, xy in pos.items()}

    # 4) Framing (auto_orient -> PCA axis horizontal).
    phi = orient if auto_orient else 0.0
    c = math.cos(phi)
    s = math.sin(phi)
    R = np.array([[c, -s], [s, c]])
    ctr = np.mean(list(pos.values()), axis=0) if pos else np.zeros(2)
    pos = {n: R.dot(xy - ctr) + ctr for n, xy in pos.items()}

    if meta_out is not None:
        def _bucket_loop(arr):
            arr = np.asarray(arr, float)
            if abs(sb) <= 1.0e-9:
                return arr
            return np.array([_bucket(px, py) for px, py in arr])

        def _frame(arr):
            return np.array([R.dot(np.asarray(p, float) - ctr) + ctr for p in arr])

        def _invert_loop(arr):
            if not invert_ring:
                return arr
            return np.array([_invert_pt(px, py) for px, py in np.asarray(arr, float)])

        # When inverted, the drawn 'outer' ring is the (reflected) inner outline and
        # vice versa, so the big ring stays the outer boundary in the overlay.
        outer_loop = (hr * dense) if invert_ring else dense
        inner_loop = dense if invert_ring else (hr * dense)
        meta_out.clear()
        meta_out["kind"] = "holed"
        meta_out["pinned_nodes"] = list(set(outer) | set(hole))
        meta_out["aspect_value"] = float(aspect)
        meta_out["boundary_outer"] = _frame(_bucket_loop(_invert_loop(outer_loop)))
        meta_out["boundary_inner"] = _frame(_bucket_loop(_invert_loop(inner_loop)))
        meta_out["medial"] = _frame(_bucket_loop(_invert_loop(0.5 * (1.0 + hr) * dense)))
    return pos


def _transform_points_like(arr, center, y_direction="top-to-bottom", rotate_degrees=0.0):
    """Apply the same drawing transform as ``transform_positions`` to a point array."""
    arr = np.asarray(arr, float)
    if arr.size == 0:
        return arr
    center = np.asarray(center, float)
    theta = math.radians(float(rotate_degrees or 0.0))
    c = math.cos(theta)
    s = math.sin(theta)
    R = np.array([[c, -s], [s, c]])
    out = np.zeros_like(arr)
    for i, xy in enumerate(arr):
        q = np.asarray(xy, float) - center
        if y_direction == "top-to-bottom":
            q = np.array([q[0], -q[1]])
        elif y_direction == "bottom-to-top":
            pass
        else:
            raise ValueError("Unknown y_direction %r" % y_direction)
        out[i] = R.dot(q)
    return out


def _transform_tutte_guides(meta, center, y_direction="top-to-bottom", rotate_degrees=0.0):
    """Transform every guide-geometry array in ``meta`` like ``transform_positions``."""
    if not meta:
        return meta
    out = {}
    for k, v in meta.items():
        # Non-geometry entries (e.g. the 'kind' tag) pass through untouched.
        if not isinstance(v, np.ndarray):
            out[k] = v
            continue
        out[k] = _transform_points_like(v, center, y_direction, rotate_degrees)
    return out


def _closest_points_segments(p1, p2, q1, q2):
    """Closest points between 2D segments p1p2 and q1q2 (clamped)."""
    d1 = p2 - p1
    d2 = q2 - q1
    r = p1 - q1
    a = float(d1 @ d1)
    e = float(d2 @ d2)
    f = float(d2 @ r)
    if a <= 1.0e-12 and e <= 1.0e-12:
        return p1, q1
    if a <= 1.0e-12:
        s = 0.0
        t = min(max(f / e, 0.0), 1.0)
    else:
        c = float(d1 @ r)
        if e <= 1.0e-12:
            t = 0.0
            s = min(max(-c / a, 0.0), 1.0)
        else:
            b = float(d1 @ d2)
            denom = a * e - b * b
            s = min(max((b * f - c * e) / denom, 0.0), 1.0) if denom > 1.0e-12 else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = min(max(-c / a, 0.0), 1.0)
            elif t > 1.0:
                t = 1.0
                s = min(max((b - c) / a, 0.0), 1.0)
    return p1 + d1 * s, q1 + d2 * t


def nudge_min_separation(P, G, min_sep, iterations=8, damping=0.5, anchor=0.12):
    """
    Post-layout relaxation that pushes apart non-incident graph edges (strand
    pieces) that come closer than ``min_sep`` (a fraction of the diagram span),
    while a weak spring to each node's original position preserves the overall
    layout.  This opens up shallow/near-parallel strand runs that would otherwise
    read as a single thick line.  ``min_sep`` <= 0 disables it.

    The whole layout is treated uniformly (no pinned boundary), so it applies to
    every 2D layout; the spring keeps boundary shapes close to their target.
    """
    if not P or float(min_sep or 0.0) <= 0.0:
        return P
    nodes = list(P.keys())
    coords = {n: np.asarray(P[n], float).copy() for n in nodes}
    orig = {n: coords[n].copy() for n in nodes}
    allc = np.array([orig[n] for n in nodes])
    span = float(np.linalg.norm(allc.max(axis=0) - allc.min(axis=0))) or 1.0
    thr = float(min_sep) * span
    if thr <= 0.0:
        return P
    edges = [(u, v) for u, v in G.edges() if u != v]
    ne = len(edges)
    if ne < 2:
        return P

    for _ in range(int(iterations)):
        # Per-edge bounding boxes (expanded by thr) for fast rejection.
        boxes = []
        for (u, v) in edges:
            a = coords[u]
            b = coords[v]
            boxes.append((
                min(a[0], b[0]) - thr, max(a[0], b[0]) + thr,
                min(a[1], b[1]) - thr, max(a[1], b[1]) + thr,
            ))
        disp = {n: np.zeros(2) for n in nodes}
        active = False
        for i in range(ne):
            u1, v1 = edges[i]
            bx = boxes[i]
            for j in range(i + 1, ne):
                u2, v2 = edges[j]
                if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                    continue  # incident edges share a node -> skip
                bx2 = boxes[j]
                if bx[1] < bx2[0] or bx2[1] < bx[0] or bx[3] < bx2[2] or bx2[3] < bx[2]:
                    continue  # bounding boxes farther than thr apart
                cp1, cp2 = _closest_points_segments(
                    coords[u1], coords[v1], coords[u2], coords[v2]
                )
                d = cp2 - cp1
                dist = float(np.linalg.norm(d))
                if dist >= thr:
                    continue
                if dist < 1.0e-9:
                    seg = coords[v1] - coords[u1]
                    n_hat = np.array([-seg[1], seg[0]], float)
                    nn = float(np.linalg.norm(n_hat))
                    n_hat = n_hat / nn if nn > 1.0e-9 else np.array([1.0, 0.0])
                else:
                    n_hat = d / dist
                push = 0.5 * (thr - dist)
                disp[u1] -= 0.5 * push * n_hat
                disp[v1] -= 0.5 * push * n_hat
                disp[u2] += 0.5 * push * n_hat
                disp[v2] += 0.5 * push * n_hat
                active = True
        if not active:
            break
        for n in nodes:
            coords[n] = coords[n] + damping * disp[n]
            coords[n] = coords[n] + anchor * (orig[n] - coords[n])
    return {n: coords[n] for n in nodes}


def sphere_stereo_layout(model, G, crossing_angle_degrees=8.0, equalize=0.0):
    """
    V5.2 2D layout for links that defeat boundary-pinned layouts (rings in
    near-orthogonal planes, e.g. Edwards-Venn AM_n).

    The crossing graph is laid out on the unit sphere exactly as the 3D XYZ
    pipeline does (3D Kamada-Kawai unit directions + local crossing-gadget
    compaction), a projection pole is chosen as the direction farthest (in
    angle) from every node and strand-arc sample, and the whole spherical
    diagram is stereographically projected from that pole.  Stereographic
    projection is a homeomorphism of the sphere minus the pole, so the image
    is a genuine planar diagram of the same link; the face containing the pole
    becomes the outer region.
    """
    raw_dirs = _kamada_3d_unit_directions(G)
    P3, centers3 = _compact_crossing_gadgets_on_sphere(
        model, raw_dirs, crossing_angle=math.radians(float(crossing_angle_degrees))
    )
    nodes = list(G.nodes())
    dirs = np.array([P3[n] for n in nodes], float)

    # Sample strand arcs (out-corner -> seg-mid -> in-corner great arcs) so the
    # pole avoids strands, not just graph nodes.
    samples = [dirs]
    for cp in model["comp_positions"]:
        for p in cp:
            k = model["pos_cross"][p]
            role = model["pos_role"][p]
            q = model["nextpos"][p]
            k2 = model["pos_cross"][q]
            role2 = model["pos_role"][q]
            a = P3[(k, "out_" + role)]
            m = P3[("seg", p)]
            b = P3[(k2, "in_" + role2)]
            samples.append(_spherical_quadratic_bezier(a, m, b, 9))
    cloud = _normalize_rows(np.vstack(samples))

    # Pole = Fibonacci-sphere candidate with the largest angular clearance.
    cands = _fibonacci_sphere_directions(1500)
    dots = np.clip(cands @ cloud.T, -1.0, 1.0)
    clearance = np.arccos(dots).min(axis=1)
    pole = cands[int(np.argmax(clearance))]

    # Stereographic projection from the pole onto the plane tangent at -pole.
    t1, t2 = _tangent_basis_at(pole)
    pos = {}
    for n, d in zip(nodes, dirs):
        w = 1.0 - float(np.dot(d, pole))
        if w < 1.0e-9:
            w = 1.0e-9
        pos[n] = np.array([float(np.dot(d, t1)) / w, float(np.dot(d, t2)) / w])

    # Optional radial histogram equalization about the centroid (the 'ring
    # equalize' knob): stereographic projection compresses the anti-pole region
    # into the centre; equalizing the radial distribution (angles kept) opens
    # it up.  The false-crossing audit still applies downstream.
    s_eq = min(1.0, max(0.0, float(equalize or 0.0)))
    if s_eq > 0.0 and len(pos) >= 3:
        # Smooth power-law radial expansion about the centroid (gamma < 1
        # opens the compressed anti-pole centre).  A smooth monotone radial
        # map is a homeomorphism of the plane, so it distorts the (straight)
        # edges far less than a per-node rank remap would.
        ctr = np.mean(list(pos.values()), axis=0)
        rr = {n: float(np.linalg.norm(pos[n] - ctr)) for n in pos}
        r_max = max(rr.values()) or 1.0
        gamma = 1.0 - 0.6 * s_eq
        for n in pos:
            if rr[n] > 1.0e-12:
                r_new = r_max * (rr[n] / r_max) ** gamma
                pos[n] = ctr + (pos[n] - ctr) * (r_new / rr[n])
    return pos


def relax_planar_layout(P, G, passes=0, strength=0.5, pinned=None):
    """
    V5.0 planarity-guarded relaxation, applied after any 2D layout.

    Each pass moves every non-pinned node toward the centroid of its graph
    neighbours while nudging incident edge lengths toward the median edge
    length.  After the whole sweep the layout is checked for straight-edge
    crossings between non-incident edges; a sweep that would introduce one is
    rejected and retried at half strength, so planarity is preserved by
    construction.  ``pinned`` is an iterable of node keys to hold fixed (the
    layout's pinned boundary faces when known); when None, the convex hull of
    the input positions is held fixed so the outer shape survives.
    """
    passes = int(passes or 0)
    strength = float(strength or 0.0)
    if not P or passes <= 0 or strength <= 0.0:
        return P
    nodes = list(P.keys())
    coords = {n: np.asarray(P[n], float).copy() for n in nodes}
    edges = [(u, v) for u, v in G.edges() if u != v]
    if len(edges) < 2:
        return P

    if pinned is None:
        pts = np.array([coords[n] for n in nodes])
        try:
            hull_idx = set()
            c0 = pts.mean(axis=0)
            rel = pts - c0
            ang = np.arctan2(rel[:, 1], rel[:, 0])
            rr = np.linalg.norm(rel, axis=1)
            # Nodes that are radially extremal in their angular bin approximate
            # the hull cheaply and robustly for our near-convex boundaries.
            nbin = max(16, int(len(nodes) ** 0.5))
            for b in range(nbin):
                lo = -np.pi + 2.0 * np.pi * b / nbin
                hi = lo + 2.0 * np.pi / nbin
                sel = np.where((ang >= lo) & (ang < hi))[0]
                if len(sel):
                    hull_idx.add(int(sel[np.argmax(rr[sel])]))
            pinned_set = {nodes[i] for i in hull_idx}
        except Exception:
            pinned_set = set()
    else:
        pinned_set = set(pinned)

    node_index = {n: i for i, n in enumerate(nodes)}
    e_u = np.array([node_index[u] for (u, v) in edges])
    e_v = np.array([node_index[v] for (u, v) in edges])
    ii, jj = np.triu_indices(len(edges), k=1)
    # Exclude incident edge pairs (sharing an endpoint) once, up front.
    share = ((e_u[ii] == e_u[jj]) | (e_u[ii] == e_v[jj]) |
             (e_v[ii] == e_u[jj]) | (e_v[ii] == e_v[jj]))
    ii, jj = ii[~share], jj[~share]

    def _has_crossing(cc):
        pts = np.array([cc[n] for n in nodes])
        A = pts[e_u]; B = pts[e_v]
        # Bounding-box rejection, vectorized over the precomputed pairs.
        lo = np.minimum(A, B); hi = np.maximum(A, B)
        ok = ~((hi[ii, 0] < lo[jj, 0]) | (hi[jj, 0] < lo[ii, 0]) |
               (hi[ii, 1] < lo[jj, 1]) | (hi[jj, 1] < lo[ii, 1]))
        if not ok.any():
            return False
        i2, j2 = ii[ok], jj[ok]
        p, q = A[i2], B[i2]
        r, s = A[j2], B[j2]
        d1 = np.cross(q - p, r - p)
        d2 = np.cross(q - p, s - p)
        d3 = np.cross(s - r, p - r)
        d4 = np.cross(s - r, q - r)
        return bool(np.any((d1 * d2 < 0.0) & (d3 * d4 < 0.0)))

    step = strength
    for _p in range(passes):
        if step < 1.0e-3:
            break
        lens = [float(np.linalg.norm(coords[u] - coords[v])) for (u, v) in edges]
        L_med = float(np.median([x for x in lens if x > 0.0]) or 0.0)
        if L_med <= 0.0:
            break
        proposal = {n: coords[n].copy() for n in nodes}
        for n in nodes:
            if n in pinned_set:
                continue
            nbrs = list(G[n])
            if not nbrs:
                continue
            p0 = coords[n]
            cen = np.mean([coords[b] for b in nbrs], axis=0)
            pull = np.zeros(2)
            for b in nbrs:
                d = coords[b] - p0
                dist = float(np.linalg.norm(d))
                if dist > 1.0e-12:
                    pull += (d / dist) * (dist - L_med)
            move = 0.55 * (cen - p0) + 0.45 * (pull / max(1, len(nbrs)))
            proposal[n] = p0 + step * move
        if _has_crossing(proposal):
            step *= 0.5
            continue
        coords = proposal
    return {n: coords[n] for n in nodes}


def compute_positions_connected(G, layout, tutte_opts=None, meta_out=None):
    ok, emb = nx.check_planarity(G)
    if not ok:
        raise RuntimeError(
            "The crossing graph is not planar; the DT code may be non-realizable."
        )

    opts = tutte_opts or {}
    if layout == "holed-tutte":
        # 'holed-tutte' pins two faces (outer + auto-picked central hole) to the
        # outer/inner outlines of a holed shape and harmonically solves the ring.
        shape = _normalize_tutte_shape(opts.get("shape", "circle"))
        aspect = float(opts.get("aspect", 1.0))
        corner = float(opts.get("corner_radius", 0.0))
        auto_aspect = bool(opts.get("auto_aspect", opts.get("auto", False)))
        try:
            return holed_tutte_layout(
                G,
                emb,
                shape=shape,
                aspect=aspect,
                corner_radius=corner,
                hole_ratio=float(opts.get("hole_ratio", 0.4)),
                swap=bool(opts.get("hole_swap", False)),
                auto_aspect=auto_aspect,
                orient_degrees=float(opts.get("orient", 0.0)),
                auto_orient=bool(opts.get("auto_orient", True)),
                ring_tilt=float(opts.get("ring_tilt", 90.0)),
                invert_ring=bool(opts.get("invert_ring", False)),
                ring_equalize=float(opts.get("ring_equalize", 0.0)),
                meta_out=meta_out,
            )
        except np.linalg.LinAlgError:
            return {n: np.asarray(xy, float) for n, xy in nx.planar_layout(G).items()}
    if layout == "shaped-tutte":
        # 'shaped-tutte' pins the boundary to a chosen convex shape and tilts the
        # shape's aspect axis relative to the diagram's intrinsic PCA axis.
        shape = _normalize_tutte_shape(opts.get("shape", "circle"))
        aspect = float(opts.get("aspect", 1.0))
        corner = float(opts.get("corner_radius", 0.0))
        auto_aspect = bool(opts.get("auto_aspect", opts.get("auto", False)))
        try:
            return shaped_tutte_layout(
                G,
                emb,
                shape=shape,
                aspect=aspect,
                corner_radius=corner,
                decompress=float(opts.get("decompress", 0.0)),
                auto_aspect=auto_aspect,
                com_expand=float(opts.get("com_expand", 0.0)),
                orient_degrees=float(opts.get("orient", 0.0)),
                auto_orient=bool(opts.get("auto_orient", True)),
                meta_out=meta_out,
            )
        except np.linalg.LinAlgError:
            # Some connected planar graphs are singular for this simple
            # barycentric solve. NetworkX planar_layout is a safe fallback.
            return {n: np.asarray(xy, float) for n, xy in nx.planar_layout(G).items()}
    if layout == "tutte":
        # Plain 'tutte' always uses the classic unit-circle boundary.
        try:
            return tutte_layout_connected(
                G,
                emb,
                shape="circle",
                aspect=1.0,
                corner_radius=0.0,
                decompress=float(opts.get("decompress", 0.0)),
                auto_aspect=False,
                com_expand=float(opts.get("com_expand", 0.0)),
            )
        except np.linalg.LinAlgError:
            return {n: np.asarray(xy, float) for n, xy in nx.planar_layout(G).items()}
    if layout == "planar":
        return {n: np.asarray(xy, float) for n, xy in nx.planar_layout(G).items()}
    if layout == "spring":
        p0 = nx.planar_layout(G)
        return {
            n: np.asarray(xy, float)
            for n, xy in nx.spring_layout(G, pos=p0, seed=1, iterations=200).items()
        }
    if layout == "kamada":
        return {n: np.asarray(xy, float) for n, xy in nx.kamada_kawai_layout(G).items()}
    if layout == "sphere-stereo":
        model = opts.get("model")
        if model is None:
            raise ValueError("sphere-stereo layout needs the DT model (internal).")
        pos = sphere_stereo_layout(
            model, G, equalize=float(opts.get("ring_equalize", 0.0))
        )
        if meta_out is not None:
            meta_out.clear()
            meta_out["kind"] = "sphere-stereo"
        return pos
    raise ValueError("Unknown layout %r" % layout)


def compute_positions(G, layout, tutte_opts=None, meta_out=None):
    """Compute graph coordinates; disconnected graph pieces are packed side-by-side."""
    if G.number_of_nodes() == 0:
        return {}

    components = [list(nodes) for nodes in nx.connected_components(G)]
    if len(components) == 1:
        return compute_positions_connected(
            G, layout, tutte_opts=tutte_opts, meta_out=meta_out
        )

    # Guide overlays are only defined for a single connected shaped-tutte diagram.
    if meta_out is not None:
        meta_out.clear()
    packed = {}
    x_offset = 0.0
    gap = 0.75
    for nodes in components:
        H = G.subgraph(nodes).copy()
        P = compute_positions_connected(H, layout, tutte_opts=tutte_opts)
        coords = np.array(list(P.values()))
        mn = coords.min(axis=0)
        mx = coords.max(axis=0)
        span = max(float(mx[0] - mn[0]), float(mx[1] - mn[1]), 1.0e-9)
        width = float(mx[0] - mn[0]) / span
        for n, xy in P.items():
            q = (xy - mn) / span
            q[0] += x_offset
            packed[n] = q
        x_offset += width + gap
    return packed


# --------------------------------------------------------------------------- #
#  4. Coordinate transforms and geometry helpers
# --------------------------------------------------------------------------- #
def transform_positions(P, y_direction="top-to-bottom", rotate_degrees=0.0):
    """
    Apply final drawing-coordinate convention.

    y_direction='top-to-bottom' reflects the y coordinate and is the default.
    Rotation is then applied around the diagram center in degrees.
    """
    if not P:
        return {}
    coords = np.array(list(P.values()))
    center = coords.mean(axis=0)
    theta = math.radians(float(rotate_degrees or 0.0))
    c = math.cos(theta)
    s = math.sin(theta)
    R = np.array([[c, -s], [s, c]])

    out = {}
    for key, xy in P.items():
        q = np.asarray(xy, float) - center
        if y_direction == "top-to-bottom":
            q[1] = -q[1]
        elif y_direction == "bottom-to-top":
            pass
        else:
            raise ValueError("Unknown y_direction %r" % y_direction)
        out[key] = R.dot(q)
    return out


def catmull_rom(points, samples=26):
    """Closed Catmull-Rom spline. Returns dense points and span start indices."""
    pts = np.asarray(points, float)
    n = len(pts)
    if n == 0:
        return np.zeros((0, 2)), []
    if n == 1:
        return pts.copy(), [0]

    out = []
    starts = []
    for i in range(n):
        p0 = pts[(i - 1) % n]
        p1 = pts[i]
        p2 = pts[(i + 1) % n]
        p3 = pts[(i + 2) % n]
        starts.append(len(out))
        for s in range(samples):
            t = s / float(samples)
            q = 0.5 * (
                (2.0 * p1)
                + (-p0 + p2) * t
                + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t * t
                + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t * t * t
            )
            out.append(q)
    return np.asarray(out), starts


def crossing_centers(model, P):
    centers = {}
    for k in range(len(model["crossings"])):
        centers[k] = np.mean(
            [P[(k, r)] for r in ("in_o", "in_e", "out_o", "out_e")], axis=0
        )
    return centers


def component_controls(model, P, centers, ci):
    """
    Ordered control points for component ci.

    Returns (controls, center_marks), where center_marks contains
    (control_index, traversal_position, is_over).
    """
    pos_cross = model["pos_cross"]
    pos_role = model["pos_role"]
    over_at = model["over_at"]
    ctrl = []
    center_marks = []
    for p in model["comp_positions"][ci]:
        k = pos_cross[p]
        ctrl.append(P[(k, "in_" + pos_role[p])])
        center_marks.append((len(ctrl), p, over_at[p]))
        ctrl.append(centers[k])
        ctrl.append(P[(k, "out_" + pos_role[p])])
        ctrl.append(P[("seg", p)])
    return ctrl, center_marks


def over_position_for_crossing(model, k):
    c = model["crossings"][k]
    odd = c["odd"]
    even = c["even"]
    if model["over_at"][even]:
        return even
    return odd


# --------------------------------------------------------------------------- #
#  5. False-crossing audit
# --------------------------------------------------------------------------- #
def _seg_intersect(p1, p2, p3, p4):
    r = p2 - p1
    s = p4 - p3
    rxs = r[0] * s[1] - r[1] * s[0]
    if abs(rxs) < 1.0e-12:
        return False
    qp = p3 - p1
    t = (qp[0] * s[1] - qp[1] * s[0]) / rxs
    u = (qp[0] * r[1] - qp[1] * r[0]) / rxs
    return (1.0e-6 < t < 1.0 - 1.0e-6) and (1.0e-6 < u < 1.0 - 1.0e-6)


def audit_false_crossings(model, P, centers, samples=26):
    """Return the number of segment intersections away from true crossing disks."""
    if not centers:
        return 0

    cvals = np.array(list(centers.values()))
    scale = float(np.linalg.norm(cvals.max(axis=0) - cvals.min(axis=0))) or 1.0
    R = 0.06 * scale

    curves = []
    for ci in range(len(model["comp_positions"])):
        ctrl, _ = component_controls(model, P, centers, ci)
        dense, _ = catmull_rom(ctrl, samples)
        curves.append(dense)

    segs = []
    for ci, dense in enumerate(curves):
        n = len(dense)
        if n < 2:
            continue
        mask = np.ones(n, bool)
        for k in centers:
            mask[np.linalg.norm(dense - centers[k], axis=1) < R] = False
        for i in range(n):
            j = (i + 1) % n
            if mask[i] and mask[j]:
                segs.append((ci, i, n, dense[i], dense[j]))

    false = 0
    for a in range(len(segs)):
        cia, ia, na, pa, pb = segs[a]
        for b in range(a + 1, len(segs)):
            cib, ib, nb, pc, pd = segs[b]
            if cia == cib:
                circular_gap = min(abs(ia - ib), na - abs(ia - ib))
                if circular_gap <= 1:
                    continue
            if _seg_intersect(pa, pb, pc, pd):
                false += 1
    return false


def _seg_intersect_point(p1, p2, p3, p4):
    """Return the intersection point of segments p1p2 and p3p4, or None."""
    p1 = np.asarray(p1, float)
    p2 = np.asarray(p2, float)
    p3 = np.asarray(p3, float)
    p4 = np.asarray(p4, float)
    r = p2 - p1
    s = p4 - p3
    rxs = r[0] * s[1] - r[1] * s[0]
    if abs(rxs) < 1.0e-12:
        return None
    qp = p3 - p1
    t = (qp[0] * s[1] - qp[1] * s[0]) / rxs
    u = (qp[0] * r[1] - qp[1] * r[0]) / rxs
    if (1.0e-6 < t < 1.0 - 1.0e-6) and (1.0e-6 < u < 1.0 - 1.0e-6):
        return p1 + t * r
    return None


def find_false_crossings_in_curves(dense_curves, true_centers_xy, exclude_radius):
    """Locate false crossings among already-rendered curves.

    ``dense_curves`` is a list of Nx2 arrays (one per component, in final drawing
    coordinates).  Segment intersections that fall within ``exclude_radius`` of a
    true crossing center are ignored, so only the accidental overlaps introduced
    by the layout are reported.  Returns a list of dicts, each describing one
    false crossing with the two participating strands::

        {"point": (x, y),
         "a": {"curve": ci_a, "index": i_a},
         "b": {"curve": ci_b, "index": i_b}}
    """
    centers = [np.asarray(c, float) for c in (true_centers_xy or [])]
    excl2 = float(exclude_radius) ** 2

    def _near_true_center(pt):
        for c in centers:
            d = pt - c
            if float(d[0] * d[0] + d[1] * d[1]) < excl2:
                return True
        return False

    segs = []
    for ci, dense in enumerate(dense_curves):
        pts = np.asarray(dense, float)
        n = len(pts)
        if n < 3:
            continue
        for i in range(n):
            a = pts[i]
            b = pts[(i + 1) % n]
            mid = 0.5 * (a + b)
            if _near_true_center(a) or _near_true_center(b) or _near_true_center(mid):
                continue
            segs.append((ci, i, n, a, b))

    found = []
    for x in range(len(segs)):
        cia, ia, na, pa, pb = segs[x]
        for y in range(x + 1, len(segs)):
            cib, ib, nb, pc, pd = segs[y]
            if cia == cib:
                circular_gap = min(abs(ia - ib), na - abs(ia - ib))
                if circular_gap <= 1:
                    continue
            pt = _seg_intersect_point(pa, pb, pc, pd)
            if pt is None:
                continue
            found.append({
                "point": (float(pt[0]), float(pt[1])),
                "a": {"curve": cia, "index": ia},
                "b": {"curve": cib, "index": ib},
            })
    return found


# --------------------------------------------------------------------------- #
#  6. Drawing
# --------------------------------------------------------------------------- #
DEFAULT_PALETTE = plt.cm.tab10(np.linspace(0, 1, 10))

# Default 2-D drawing settings (single source of truth).  The CLI argparse
# defaults below match these, and external consumers such as the strand-passage
# engine import them so they can draw links exactly the way this helper does by
# default.
DEFAULT_LAYOUT = "tutte"
DEFAULT_Y_DIRECTION = "top-to-bottom"
DEFAULT_ROTATE = 0.0
ASSET_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
APP_ICON_PNG = os.path.join(ASSET_DIR, "strand_passage_icon.png")


def apply_tk_window_icon(window, tk_module=None):
    """Apply the project icon to a Tk/Toplevel window when the PNG asset exists.

    The icon is intentionally optional: source checkouts, packaged copies, and
    Sage environments should keep running even if the asset is omitted or the
    active Tk build cannot decode PNG files.
    """
    if not os.path.exists(APP_ICON_PNG):
        return False
    try:
        if tk_module is None:
            import tkinter as tk_module  # type: ignore
        image = tk_module.PhotoImage(file=APP_ICON_PNG)
        window.iconphoto(True, image)
        window._strand_passage_icon = image  # keep Tk image alive
        return True
    except Exception:  # noqa: BLE001
        return False


def _default_color_of(ci):
    return DEFAULT_PALETTE[ci % len(DEFAULT_PALETTE)]


def _add_arrows(ax, dense, centers, color, scale, n_arrows=2, lw=2.0):
    if len(dense) < 4 or not centers:
        return
    # Arrowhead size scales with the strand line width (mutation_scale is in
    # points, as is lw, so this scales identically in the preview and saved image).
    try:
        head = max(6.0, 9.0 * float(lw))
    except Exception:
        head = 18.0
    n = len(dense)
    cvals = np.array(list(centers.values()))
    dmin = np.min(
        np.linalg.norm(dense[:, None, :] - cvals[None, :, :], axis=2), axis=1
    )
    order = np.argsort(-dmin)
    chosen = []
    min_sep = max(1, n // max(1, 2 * n_arrows))
    for i in order:
        if all(min(abs(int(i) - j), n - abs(int(i) - j)) > min_sep for j in chosen):
            chosen.append(int(i))
        if len(chosen) >= n_arrows:
            break
    for i in chosen:
        a = dense[i]
        b = dense[(i + 2) % n]
        ann = ax.annotate(
            "",
            xy=b,
            xytext=a,
            arrowprops=dict(arrowstyle="-|>", color=color, lw=0, mutation_scale=head),
            zorder=5,
        )
        try:
            ann.set_clip_on(False)
            if ann.arrow_patch is not None:
                ann.arrow_patch.set_clip_on(False)
        except Exception:
            pass



def _safe_canvas_draw(fig):
    """Draw a Matplotlib figure just enough to obtain text bounding boxes."""
    try:
        fig.canvas.draw()
        return fig.canvas.get_renderer()
    except Exception:
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        return canvas.get_renderer()


def _move_dt_labels_away_from_crossing_ids(ax, dt_label_entries, crossing_id_entries,
                                           max_passes=10):
    """
    Resolve text-box collisions by moving DT labels, not crossing IDs.

    Crossing IDs stay fixed at the true crossing centers. If a displayed crossing
    ID box overlaps a DT traversal-label box, the DT label is nudged radially
    away from the crossing ID in display coordinates. The final text positions
    are copied back into the returned CSV coordinate table.
    """
    if not dt_label_entries or not crossing_id_entries:
        return

    fig = ax.figure
    renderer = _safe_canvas_draw(fig)
    inv = ax.transData.inverted()

    def _artist_bbox(artist, expand_x=1.0, expand_y=1.0):
        return artist.get_window_extent(renderer=renderer).expanded(expand_x, expand_y)

    for pass_index in range(int(max_passes)):
        id_infos = []
        for entry in crossing_id_entries:
            artist = entry.get("artist")
            if artist is None or not artist.get_visible():
                continue
            bbox = _artist_bbox(artist, 1.10, 1.15)
            center = entry.get("center")
            if center is None:
                center = getattr(artist, "center", (0.0, 0.0))
            center_data = np.asarray(center, float)
            id_infos.append((bbox, center_data))
        if not id_infos:
            break

        changed = False
        for entry in dt_label_entries:
            artist = entry.get("artist")
            if artist is None or not artist.get_visible():
                continue
            dt_bbox = _artist_bbox(artist, 1.03, 1.08)
            overlaps = [(bbox, center) for bbox, center in id_infos if dt_bbox.overlaps(bbox)]
            if not overlaps:
                continue

            current = np.asarray(artist.get_position(), float)
            nearest_center = min(
                (np.asarray(center, float) for _bbox, center in overlaps),
                key=lambda center: float(np.linalg.norm(current - center)),
            )

            current_disp = ax.transData.transform(current)
            center_disp = ax.transData.transform(nearest_center)
            direction_disp = current_disp - center_disp
            norm = float(np.linalg.norm(direction_disp))
            if norm <= 1.0e-9:
                fallback = np.asarray(entry.get("direction", [1.0, 0.0]), float)
                if float(np.linalg.norm(fallback)) <= 1.0e-12:
                    fallback = np.array([1.0, 0.0], float)
                direction_disp = ax.transData.transform(current + fallback) - current_disp
                norm = float(np.linalg.norm(direction_disp))
            if norm <= 1.0e-9:
                direction_disp = np.array([1.0, 0.0], float)
                norm = 1.0
            direction_disp = direction_disp / norm

            # Small, repeated pixel nudges keep the label near its DT visit while
            # ensuring the fixed crossing-ID disk no longer covers the label.
            step_pixels = 6.0 + 2.5 * float(pass_index)
            new_disp = current_disp + step_pixels * direction_disp
            new_data = inv.transform(new_disp)
            artist.set_position((float(new_data[0]), float(new_data[1])))
            entry["position"] = np.asarray(new_data, float)
            changed = True

        if not changed:
            break
        renderer = _safe_canvas_draw(fig)

def _axis_content_bounds(ax):
    """Return (xmin, xmax, ymin, ymax) for visible 2D preview content."""
    xs = []
    ys = []
    for patch in ax.patches:
        try:
            if not patch.get_visible():
                continue
            path = patch.get_path()
            verts = np.asarray(path.vertices, float)
        except Exception:
            continue
        if verts.size:
            mask = np.isfinite(verts[:, 0]) & np.isfinite(verts[:, 1])
            if np.any(mask):
                xs.extend(verts[mask, 0].tolist())
                ys.extend(verts[mask, 1].tolist())
    for line in ax.lines:
        try:
            x = np.asarray(line.get_xdata(), float)
            y = np.asarray(line.get_ydata(), float)
        except Exception:
            continue
        mask = np.isfinite(x) & np.isfinite(y)
        if np.any(mask):
            xs.extend(x[mask].tolist())
            ys.extend(y[mask].tolist())
    for txt in ax.texts:
        try:
            if not txt.get_visible():
                continue
            tr = txt.get_transform()
            if hasattr(tr, "contains_branch") and not tr.contains_branch(ax.transData):
                continue
            x, y = txt.get_position()
            x = float(x)
            y = float(y)
        except Exception:
            continue
        if math.isfinite(x) and math.isfinite(y):
            xs.append(x)
            ys.append(y)
    if not xs or not ys:
        return None
    xmin = float(min(xs))
    xmax = float(max(xs))
    ymin = float(min(ys))
    ymax = float(max(ys))
    if not (math.isfinite(xmin) and math.isfinite(xmax) and math.isfinite(ymin) and math.isfinite(ymax)):
        return None
    if abs(xmax - xmin) <= 1.0e-12:
        xmin -= 0.5
        xmax += 0.5
    if abs(ymax - ymin) <= 1.0e-12:
        ymin -= 0.5
        ymax += 0.5
    return (xmin, xmax, ymin, ymax)


def apply_content_framing(ax, content_bounds, aspect=1.0, zoom=1.0,
                          base_pad_frac=0.08, min_pad_frac=0.010):
    """Frame an axis around its content the same way the live preview does.

    This is the single source of truth shared by the live Tk preview and by the
    file-save path, so a saved image reproduces exactly what is shown on screen.

    ``content_bounds`` is ``(xmin, xmax, ymin, ymax)`` in data coordinates.
    ``aspect`` is the target width/height ratio of the drawing region (the
    preview canvas, or the saved figure).  ``zoom`` is the preview zoom factor:
    ``zoom = 1`` gives a comfortable padding; increasing it only trims padding
    down to a small safety margin, it never crops the diagram.
    Returns True when limits were applied.
    """
    if content_bounds is None:
        return False
    try:
        xmin, xmax, ymin, ymax = [float(v) for v in content_bounds]
    except Exception:
        return False
    try:
        zoom = float(zoom)
    except Exception:
        zoom = 1.0
    zoom = max(0.05, min(50.0, zoom if zoom else 1.0))

    dx = max(xmax - xmin, 1.0e-9)
    dy = max(ymax - ymin, 1.0e-9)
    span = max(dx, dy)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)

    pad_frac = min_pad_frac + (base_pad_frac - min_pad_frac) / zoom
    hx = 0.5 * dx + pad_frac * span
    hy = 0.5 * dy + pad_frac * span

    try:
        aspect = float(aspect)
    except Exception:
        aspect = 1.0
    aspect = max(aspect, 1.0e-9)
    if hx / hy < aspect:
        hx = hy * aspect
    else:
        hy = hx / aspect

    ax.set_xlim(cx - hx, cx + hx)
    ax.set_ylim(cy - hy, cy + hy)
    return True


def _closed_curve_path_patch(dense, color, lw, zorder=2):
    """Return a closed stroked PathPatch for one component curve."""
    pts = np.asarray(dense, float)
    if pts.ndim != 2 or pts.shape[0] < 2:
        return None
    verts = np.vstack([pts, pts[0]])
    codes = [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1) + [Path.CLOSEPOLY]
    patch = PathPatch(
        Path(verts, codes),
        facecolor="none",
        edgecolor=color,
        lw=lw,
        capstyle="round",
        joinstyle="round",
        zorder=zorder,
        clip_on=False,
    )
    patch.set_gid("closed_component_strand")
    return patch


def _closed_mask_runs(mask):
    """Return index arrays for True-runs in a circular boolean mask."""
    m = np.asarray(mask, bool)
    n = len(m)
    if n == 0 or not np.any(m):
        return []
    if np.all(m):
        return [np.arange(n, dtype=int)]
    false_idx = np.where(~m)[0]
    start = int((false_idx[-1] + 1) % n)
    ordered = np.concatenate([np.arange(start, n, dtype=int), np.arange(0, start, dtype=int)])
    vals = m[ordered]
    runs = []
    run_start = None
    for i, val in enumerate(vals):
        if val and run_start is None:
            run_start = i
        end_now = (i == len(vals) - 1)
        if run_start is not None and ((not val) or end_now):
            run_end = i if (val and end_now) else i - 1
            if run_end >= run_start:
                runs.append(ordered[run_start:run_end + 1])
            run_start = None
    return runs


def _plot_local_curve_piece(ax, dense, center, radius, color, lw, zorder,
                            capstyle="round"):
    """Plot the portions of a closed curve within radius of center."""
    pts = np.asarray(dense, float)
    if pts.ndim != 2 or pts.shape[0] < 3:
        return
    c = np.asarray(center, float)
    mask = np.linalg.norm(pts - c, axis=1) <= float(radius)
    for run in _closed_mask_runs(mask):
        if len(run) == 0:
            continue
        n = len(pts)
        padded = np.concatenate([
            np.array([(int(run[0]) - 1) % n], dtype=int),
            run.astype(int),
            np.array([(int(run[-1]) + 1) % n], dtype=int),
        ])
        seg = pts[padded]
        ax.plot(
            seg[:, 0],
            seg[:, 1],
            color=color,
            lw=lw,
            solid_capstyle=capstyle,
            solid_joinstyle="round",
            zorder=zorder,
            clip_on=False,
        )


def _local_arclen_run(dense, center_index, radius):
    """Indices within ``radius`` arc length of ``center_index`` along a closed curve.

    Unlike a spatial-distance test, this follows the curve *parameter*, so on a
    self-crossing -- where the same closed component passes through one point
    twice -- it selects only the single pass that owns ``center_index``.  A
    spatial radius would grab both passes, which is why over/under gaps were
    invisible on knots with self-crossings (e.g. the trefoil): the under gap was
    immediately re-covered when the over pass was redrawn.
    """
    pts = np.asarray(dense, float)
    n = len(pts)
    if n < 3:
        return np.arange(n, dtype=int)
    radius = float(radius)
    i0 = int(center_index) % n

    back = 0
    dist = 0.0
    cur = i0
    while back < n - 1:
        nxt = (cur - 1) % n
        dist += float(np.linalg.norm(pts[cur] - pts[nxt]))
        if dist > radius:
            break
        cur = nxt
        back += 1

    fwd = 0
    dist = 0.0
    cur = i0
    while fwd < n - 1:
        nxt = (cur + 1) % n
        dist += float(np.linalg.norm(pts[cur] - pts[nxt]))
        if dist > radius:
            break
        cur = nxt
        fwd += 1

    total = back + fwd + 1
    if total >= n:
        return np.arange(n, dtype=int)
    start = (i0 - back) % n
    return np.array([(start + t) % n for t in range(total)], dtype=int)


def _plot_local_curve_piece_by_index(ax, dense, center_index, radius, color, lw,
                                     zorder, capstyle="round"):
    """Plot the single curve pass around ``center_index`` spanning ~``radius`` arc length.

    Index-based counterpart to :func:`_plot_local_curve_piece`.  Selecting the
    piece by curve parameter (rather than by spatial distance) is what makes
    over/under gaps render correctly at self-crossings, where the component
    passes through the crossing twice at nearly the same point.
    """
    pts = np.asarray(dense, float)
    if pts.ndim != 2 or pts.shape[0] < 3:
        return
    run = _local_arclen_run(pts, center_index, radius)
    if len(run) < 2:
        return
    seg = pts[run]
    ax.plot(
        seg[:, 0],
        seg[:, 1],
        color=color,
        lw=lw,
        solid_capstyle=capstyle,
        solid_joinstyle="round",
        zorder=zorder,
        clip_on=False,
    )


def _disable_figure_clipping(fig):
    """Turn off artist clipping so SVG output avoids unnecessary clipPath masks."""
    try:
        artists = fig.findobj()
    except Exception:
        artists = []
    for artist in artists:
        if hasattr(artist, "set_clip_on"):
            try:
                artist.set_clip_on(False)
            except Exception:
                pass


def _tighten_axis_to_content(ax, pad_frac=0.08, content_bounds=None, canvas_aspect=None):
    """Set axis limits to content bounds plus padding while preserving aspect."""
    bounds = content_bounds if content_bounds is not None else _axis_content_bounds(ax)
    if bounds is None:
        return None
    xmin, xmax, ymin, ymax = [float(v) for v in bounds]
    dx = max(xmax - xmin, 1.0e-9)
    dy = max(ymax - ymin, 1.0e-9)
    pad = float(pad_frac) * max(dx, dy)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    hx = 0.5 * dx + pad
    hy = 0.5 * dy + pad
    if canvas_aspect is None:
        try:
            bbox = ax.get_window_extent()
            canvas_aspect = float(bbox.width) / max(float(bbox.height), 1.0e-9)
        except Exception:
            canvas_aspect = 1.0
    canvas_aspect = max(float(canvas_aspect), 1.0e-9)
    if hx / hy < canvas_aspect:
        hx = hy * canvas_aspect
    else:
        hy = hx / canvas_aspect
    ax.set_xlim(cx - hx, cx + hx)
    ax.set_ylim(cy - hy, cy + hy)
    return (cx, cy, hx, hy)


def _maximize_axis_in_figure(ax, has_title=False):
    """Use almost the whole figure/canvas for the axis in preview/vector export."""
    try:
        if has_title:
            ax.figure.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.94)
        else:
            ax.set_position([0.005, 0.005, 0.99, 0.99])
    except Exception:
        pass


def render_diagram(
    ax,
    model,
    P,
    centers,
    crossing_ids=None,
    color_of=None,
    show_labels=True,
    show_crossing_ids=False,
    color_crossing_ids_by_overstrand=False,
    gap_frac=0.025,
    lw=2.0,
    label_fontsize=7.0,
    crossing_id_fontsize=6.0,
    label_box_pad=DT_LABEL_BOX_PAD,
    crossing_id_box_pad=CROSSING_ID_BOX_PAD,
    font_family=DIAGRAM_FONT_FAMILY,
    arrows=True,
    origin=(0.0, 0.0),
    scale_to=None,
):
    """
    Render the diagram onto a Matplotlib axis.

    Crossing IDs are kept fixed at the crossing centers. If an ID box overlaps a
    DT traversal label, only the DT label is moved slightly outward. Returns
    (label_coords, crossing_xy), both in final drawing coordinates.
    """
    comp_positions = model["comp_positions"]
    label_of = model.get("label_of", lambda p: p)
    if color_of is None:
        color_of = _default_color_of
    if crossing_ids is None:
        crossing_ids = default_crossing_ids(model)

    cvals = np.array(list(centers.values())) if centers else np.zeros((1, 2))
    span = float(np.linalg.norm(cvals.max(axis=0) - cvals.min(axis=0))) or 1.0
    sc = 1.0 if scale_to is None else (float(scale_to) / span)
    off = np.asarray(origin, float)

    def T(pt):
        return off + sc * np.asarray(pt, float)

    gap = float(gap_frac) * span * sc
    label_coords = {}
    crossing_xy = {k: T(centers[k]) for k in centers}
    dt_label_specs = []
    curve_infos = []

    # First compute/draw each component as a single closed stroked path.  The
    # over/under gaps are added afterwards as white masks plus re-drawn over
    # pieces.  This keeps the colored component curves closed in SVG editors,
    # instead of exporting each component as many broken NaN-separated paths.
    for ci, cp in enumerate(comp_positions):
        if not cp:
            continue
        color = color_of(ci)
        ctrl, marks = component_controls(model, P, centers, ci)
        ctrl = [T(c) for c in ctrl]
        dense, starts = catmull_rom(ctrl, samples=30)
        curve_infos.append({
            "ci": ci,
            "cp": cp,
            "color": color,
            "dense": dense,
            "marks": marks,
            "starts": starts,
        })

        patch = _closed_curve_path_patch(dense, color=color, lw=lw, zorder=2)
        if patch is not None:
            ax.add_patch(patch)

        if show_labels:
            pos_cross = model["pos_cross"]
            pos_role = model["pos_role"]
            for ppos in cp:
                k = pos_cross[ppos]
                corner = P[(k, "in_" + pos_role[ppos])]
                base = T(centers[k] + 1.25 * (corner - centers[k]))
                cxy = crossing_xy[k]
                direction = np.asarray(base, float) - np.asarray(cxy, float)
                if float(np.linalg.norm(direction)) <= 1.0e-12:
                    direction = np.array([1.0, 0.0], float)
                dt_label_specs.append({
                    "position": np.asarray(base, float),
                    "direction": direction,
                    "text": str(label_of(ppos)),
                    "color": color,
                    "p": ppos,
                    "crossing_index": k,
                    "crossing_center": cxy,
                })

    # Draw under-strand gaps as local white masks.  Then redraw the over-strand
    # local pieces on top so the over strand remains visually continuous.
    gap_mask_color = "white"
    for info in curve_infos:
        dense = info["dense"]
        starts = info["starts"]
        for cidx, _ppos, is_over in info["marks"]:
            if is_over or cidx >= len(starts):
                continue
            _plot_local_curve_piece_by_index(
                ax,
                dense,
                starts[cidx],
                radius=gap,
                color=gap_mask_color,
                lw=max(float(lw) * 2.35, float(lw) + 2.0),
                zorder=3,
            )

    for info in curve_infos:
        dense = info["dense"]
        starts = info["starts"]
        for cidx, _ppos, is_over in info["marks"]:
            if (not is_over) or cidx >= len(starts):
                continue
            _plot_local_curve_piece_by_index(
                ax,
                dense,
                starts[cidx],
                radius=gap * 1.25,
                color=info["color"],
                lw=lw,
                zorder=4,
            )

    # Add the same style of over/under gap at *false* crossings -- the accidental
    # strand overlaps some layouts (e.g. kamada) introduce away from true DT
    # crossings.  These are not real crossings, so the over/under choice is
    # arbitrary; we pick it deterministically (higher (curve, index) goes over)
    # so each false crossing simply reads as a clean crossing.
    dense_curves = [info["dense"] for info in curve_infos]
    false_exclude_radius = 0.06 * span * sc
    false_crossings = find_false_crossings_in_curves(
        dense_curves, list(crossing_xy.values()), false_exclude_radius
    )
    for fc in false_crossings:
        a = fc["a"]
        b = fc["b"]
        over, under = (a, b) if (a["curve"], a["index"]) >= (b["curve"], b["index"]) else (b, a)
        _plot_local_curve_piece_by_index(
            ax,
            dense_curves[under["curve"]],
            under["index"],
            radius=gap,
            color=gap_mask_color,
            lw=max(float(lw) * 2.35, float(lw) + 2.0),
            zorder=3,
        )
        _plot_local_curve_piece_by_index(
            ax,
            dense_curves[over["curve"]],
            over["index"],
            radius=gap * 1.25,
            color=curve_infos[over["curve"]]["color"],
            lw=lw,
            zorder=4,
        )
    n_false_crossings = len(false_crossings)
    # Expose the count so callers (file save, GUI log) can report it consistently
    # with what is actually drawn.
    setattr(ax, "_false_crossing_count", n_false_crossings)

    if arrows:
        for info in curve_infos:
            _add_arrows(ax, info["dense"], crossing_xy, info["color"], span * sc, lw=lw)

    dt_label_entries = []
    if show_labels:
        for spec in dt_label_specs:
            base = spec["position"]
            artist = ax.text(
                base[0],
                base[1],
                spec["text"],
                fontsize=label_fontsize,
                color=spec["color"],
                ha="center",
                va="center",
                zorder=7,
                fontweight="bold",
                fontfamily=font_family,
                clip_on=False,
                bbox=dict(boxstyle="round,pad=%s" % label_box_pad,
                          fc="white", ec="none", alpha=0.70),
            )
            entry = dict(spec)
            entry["artist"] = artist
            dt_label_entries.append(entry)
            label_coords[spec["p"]] = np.asarray(base, float)

    crossing_id_entries = []
    if show_crossing_ids:
        # Draw each ID text, then wrap it in a real Circle patch sized to snugly
        # enclose the measured text extent (a small margin beyond the glyph box)
        # rather than an over-large fixed data fraction.  The disk stays a clean
        # circle object in the SVG.
        id_texts = []
        for k in sorted(centers):
            cxy = crossing_xy[k]
            txt_color = "0.25"
            edge_color = "0.60"
            if color_crossing_ids_by_overstrand:
                over_pos = over_position_for_crossing(model, k)
                over_comp = model["comp_of"][over_pos]
                txt_color = color_of(over_comp)
                edge_color = txt_color
            t = ax.text(
                cxy[0],
                cxy[1],
                crossing_ids[k],
                fontsize=crossing_id_fontsize,
                color=txt_color,
                ha="center",
                va="center",
                zorder=5,
                fontweight="bold" if color_crossing_ids_by_overstrand else "normal",
                fontfamily=font_family,
                clip_on=False,
            )
            id_texts.append((k, cxy, t, edge_color))

        # Fallback radius (used if the text extent cannot be measured).
        fallback_r = span * sc * 0.02
        try:
            renderer = _safe_canvas_draw(ax.figure)
            inv = ax.transData.inverted()
        except Exception:
            renderer = None
            inv = None
        for k, cxy, t, edge_color in id_texts:
            r = fallback_r
            if renderer is not None and inv is not None:
                try:
                    ext = t.get_window_extent(renderer=renderer)
                    p0 = inv.transform((ext.x0, ext.y0))
                    p1 = inv.transform((ext.x1, ext.y1))
                    half_w = 0.5 * abs(float(p1[0] - p0[0]))
                    half_h = 0.5 * abs(float(p1[1] - p0[1]))
                    # Circle through the text-box corners, plus a small margin, is
                    # just enough to cover the text without being oversized.
                    r = math.hypot(half_w, half_h) * 1.12
                except Exception:
                    r = fallback_r
            disk = Circle(
                (float(cxy[0]), float(cxy[1])),
                radius=r,
                facecolor="white",
                edgecolor=edge_color,
                linewidth=0.8,
                alpha=0.78,
                zorder=4.9,
                clip_on=False,
            )
            ax.add_patch(disk)
            crossing_id_entries.append({"artist": disk, "center": cxy, "crossing_index": k})

    if show_labels and show_crossing_ids:
        _move_dt_labels_away_from_crossing_ids(ax, dt_label_entries, crossing_id_entries)
        for entry in dt_label_entries:
            label_coords[entry["p"]] = np.asarray(entry["artist"].get_position(), float)

    # Red warning about false crossings, shown on the live preview and saved into
    # the image/SVG.  Anchored in axes fraction so it stays put regardless of the
    # data limits and does not affect content framing.
    if n_false_crossings > 0:
        ax.text(
            0.5,
            0.985,
            "Warning: %d false crossing(s) present (layout artifact, not real DT crossings)"
            % n_false_crossings,
            transform=ax.transAxes,
            ha="center",
            va="top",
            color="red",
            fontsize=max(7.0, float(label_fontsize)),
            fontweight="bold",
            fontfamily=font_family,
            zorder=20,
            clip_on=False,
        )

    return label_coords, crossing_xy

def build_diagram_description(model, title=None, **info):
    """Return ``(creator, title_text, description, stamp)`` describing the diagram.

    Shared by the file metadata and by the visible caption drawn into the image,
    so both carry identical text: the generating script/version plus essential
    facts (component and crossing counts and the key drawing parameters).
    """
    try:
        script_name = os.path.basename(os.path.abspath(__file__))
    except Exception:
        script_name = "draw_dt_original_labels.py"
    n_comp = len(model.get("comp_positions", []))
    n_cross = len(model.get("crossings", []))
    creator = "%s (%s)" % (script_name, SCRIPT_VERSION)
    from datetime import datetime
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    parts = [
        "Signed DT link diagram",
        "components=%d" % n_comp,
        "crossings=%d" % n_cross,
    ]
    for key in ("dt", "layout", "crossing_order", "y_direction", "rotate",
                "negative_even", "gap_frac", "line_width", "label_fontsize",
                "false_crossings"):
        val = info.get(key)
        if val is not None:
            parts.append("%s=%s" % (key, val))
    parts.append("script=%s" % creator)
    parts.append("generated=%s" % stamp)
    description = "; ".join(parts)
    title_text = title or ("DT link diagram (%d components, %d crossings)"
                           % (n_comp, n_cross))
    return creator, title_text, description, stamp


def build_diagram_metadata(out_path, model, title=None, **info):
    """Assemble descriptive metadata embedded in the saved diagram.

    Records the generating script and version plus essential facts about the
    diagram (component and crossing counts and the key drawing parameters) so an
    exported file -- an SVG in particular -- is self-documenting when reopened in
    Illustrator/Inkscape or inspected as text.  The returned key set is tailored
    to the output format so Matplotlib's SVG/PDF/raster backends accept it.
    """
    creator, title_text, description, stamp = build_diagram_description(
        model, title=title, **info
    )

    ext = os.path.splitext(str(out_path))[1].lower()
    if ext == ".svg":
        return {
            "Creator": creator,
            "Title": title_text,
            "Description": description,
            "Date": stamp,
        }
    if ext == ".pdf":
        return {
            "Creator": creator,
            "Title": title_text,
            "Subject": description,
            "Keywords": "DT code, knot, link diagram",
        }
    # PNG / other raster formats: the Agg backend stores arbitrary text chunks.
    return {
        "Software": creator,
        "Title": title_text,
        "Description": description,
    }


def draw(
    model,
    P,
    centers,
    out_path,
    dpi=200,
    show_crossing_ids=False,
    color_crossing_ids_by_overstrand=False,
    crossing_ids=None,
    title=None,
    label_fontsize=7.0,
    crossing_id_fontsize=6.0,
    line_width=2.0,
    gap_frac=0.025,
    figsize=12.0,
    show_labels=True,
    arrows=True,
    match_view=None,
    dt_code=None,
    layout=None,
    tutte_guides=None,
    show_tutte_outline=False,
    show_tutte_pca=False,
):
    """Render and save a 2-D diagram.

    ``match_view`` (V3.10): when given as ``{"aspect": <w/h>, "zoom": <factor>}``
    the saved image reproduces the live GUI preview.  The figure is sized to the
    preview panel's aspect ratio, framed with the shared ``apply_content_framing``
    routine (same padding/zoom as the preview), and written as the whole figure
    (no ``bbox_inches='tight'`` crop).  When ``match_view`` is None the legacy
    square-figure, tight-cropped behaviour is preserved for CLI compatibility.
    """
    if match_view:
        try:
            aspect = max(float(match_view.get("aspect", 1.0)) or 1.0, 1.0e-6)
        except Exception:
            aspect = 1.0
        try:
            zoom = float(match_view.get("zoom", 1.0)) or 1.0
        except Exception:
            zoom = 1.0
        base = float(figsize)
        fig_wh = (base * aspect, base) if aspect >= 1.0 else (base, base / aspect)
        fig, ax = plt.subplots(figsize=fig_wh)
    else:
        aspect = None
        zoom = 1.0
        fig, ax = plt.subplots(figsize=(float(figsize), float(figsize)))

    label_coords, crossing_xy = render_diagram(
        ax,
        model,
        P,
        centers,
        crossing_ids=crossing_ids,
        show_crossing_ids=show_crossing_ids,
        color_crossing_ids_by_overstrand=color_crossing_ids_by_overstrand,
        label_fontsize=label_fontsize,
        crossing_id_fontsize=crossing_id_fontsize,
        lw=line_width,
        gap_frac=gap_frac,
        show_labels=show_labels,
        arrows=arrows,
    )
    _draw_tutte_guides(
        ax, tutte_guides,
        show_outline=bool(show_tutte_outline),
        show_pca=bool(show_tutte_pca),
    )
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=11, fontfamily=DIAGRAM_FONT_FAMILY)
    _maximize_axis_in_figure(ax, has_title=bool(title))
    _tighten_axis_to_content(ax, pad_frac=0.08)
    _disable_figure_clipping(fig)
    ensure_parent_dir(out_path)

    n_false = int(getattr(ax, "_false_crossing_count", 0) or 0)
    meta_info = dict(
        gap_frac=gap_frac,
        line_width=line_width,
        label_fontsize=label_fontsize,
    )
    if dt_code:
        # Collapse any line breaks/extra spaces from the GUI text box.
        meta_info["dt"] = " ".join(str(dt_code).split())
    if layout:
        meta_info["layout"] = layout
    if crossing_ids:
        # The effective displayed crossing-ID order (odd-label order), which
        # reflects a custom "crossing order" when one was supplied.
        meta_info["crossing_order"] = " ".join(str(c) for c in crossing_ids)
    if n_false > 0:
        meta_info["false_crossings"] = n_false
    metadata = build_diagram_metadata(out_path, model, title=title, **meta_info)

    # Also display the metadata as visible text at the bottom of the saved image
    # (so an SVG carries it both as file metadata and as an on-canvas caption).
    _creator, _title_text, _description, _stamp = build_diagram_description(
        model, title=title, **meta_info
    )
    ax.text(
        0.5,
        0.004,
        _description,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        color="0.35",
        fontsize=5.0,
        fontfamily=DIAGRAM_FONT_FAMILY,
        zorder=20,
        clip_on=False,
    )

    if match_view:
        # Reproduce the preview framing exactly, then save the whole figure so
        # the file matches what the live preview shows.
        apply_content_framing(ax, _axis_content_bounds(ax), aspect=aspect, zoom=zoom)
        fig.savefig(out_path, dpi=dpi, pad_inches=0.0, metadata=metadata)
    else:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pad_inches=0.03,
                    metadata=metadata)
    plt.close(fig)
    return label_coords, crossing_xy



# --------------------------------------------------------------------------- #
#  7. Spherical XYZ export
# --------------------------------------------------------------------------- #
def _normalize_vector(vec, fallback=None):
    """Return a unit vector, using fallback if vec is nearly zero."""
    arr = np.asarray(vec, float)
    nrm = float(np.linalg.norm(arr))
    if nrm > 1.0e-12:
        return arr / nrm
    if fallback is not None:
        return _normalize_vector(fallback)
    return np.array([0.0, 0.0, 1.0], float)


def _normalize_rows(arr):
    """Normalize an Nx3 array row-by-row."""
    pts = np.asarray(arr, float)
    if pts.size == 0:
        return np.zeros((0, 3), float)
    out = np.zeros((pts.shape[0], 3), float)
    last = np.array([0.0, 0.0, 1.0], float)
    for i, row in enumerate(pts):
        out[i] = _normalize_vector(row, fallback=last)
        last = out[i]
    return out


def _inverse_stereographic_to_unit_sphere(uv):
    """Map 2D plane coordinates to the unit sphere by inverse stereographic projection."""
    arr = np.asarray(uv, float)
    if arr.size == 0:
        return np.zeros((0, 3), float)
    u = arr[:, 0]
    v = arr[:, 1]
    r2 = u * u + v * v
    denom = 1.0 + r2
    xyz = np.empty((arr.shape[0], 3), float)
    xyz[:, 0] = 2.0 * u / denom
    xyz[:, 1] = 2.0 * v / denom
    xyz[:, 2] = (1.0 - r2) / denom
    return _normalize_rows(xyz)


def _sample_planar_components_for_sphere(model, P, centers, samples_per_span):
    """Sample the same smooth closed component curves used for the 2D drawing."""
    samples = []
    nspan = max(3, int(samples_per_span))
    for ci in range(len(model["comp_positions"])):
        ctrl, marks = component_controls(model, P, centers, ci)
        dense, starts = catmull_rom(ctrl, samples=nspan)
        samples.append({
            "component_index": ci,
            "xy": np.asarray(dense, float),
            "marks": marks,
            "starts": starts,
        })
    return samples


def _normalization_for_sphere(component_samples, sphere_extent):
    """
    Return center and scale for stereographic mode.

    sphere_extent is used only by --sphere-layout stereographic. It means: after
    centering the 2D diagram, the farthest sampled planar point is rescaled to
    radius sphere_extent before inverse stereographic projection.
    """
    arrays = [item["xy"] for item in component_samples if len(item["xy"])]
    if not arrays:
        return np.array([0.0, 0.0]), 1.0
    pts = np.vstack(arrays)
    mn = pts.min(axis=0)
    mx = pts.max(axis=0)
    center = 0.5 * (mn + mx)
    max_norm = float(np.max(np.linalg.norm(pts - center, axis=1)))
    if max_norm <= 1.0e-12:
        max_norm = 1.0
    scale = float(sphere_extent) / max_norm
    return center, scale


def _diagram_span_from_samples(component_samples):
    arrays = [item["xy"] for item in component_samples if len(item["xy"])]
    if not arrays:
        return 1.0
    pts = np.vstack(arrays)
    span = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
    return span if span > 1.0e-12 else 1.0


def _closed_curve_arclength_values(points):
    """Cumulative arclength values for a sampled closed curve."""
    pts = np.asarray(points, float)
    n = len(pts)
    if n == 0:
        return np.zeros(0, float), 0.0
    if n == 1:
        return np.zeros(1, float), 0.0
    seg = np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1)
    svals = np.zeros(n, float)
    svals[1:] = np.cumsum(seg[:-1])
    return svals, float(np.sum(seg))


def _crossing_radial_offsets(points, marks, starts, bump_radius, sphere_offset):
    """
    Compute radial offsets for one sampled component in stereographic mode.

    Away from crossing neighborhoods the offset is zero. Near a crossing visit,
    a smooth cosine bump is positive for over-strands and negative for
    under-strands. The neighborhood is measured along the component's sampled
    arclength, not by coordinate distance alone.
    """
    pts = np.asarray(points, float)
    offsets = np.zeros(len(pts), float)
    strengths = np.zeros(len(pts), float)
    if len(pts) == 0 or bump_radius <= 0.0 or sphere_offset == 0.0:
        return offsets

    svals, total = _closed_curve_arclength_values(pts)
    if total <= 1.0e-12:
        return offsets
    width = min(float(bump_radius), 0.25 * total)

    for control_index, _pos, is_over in marks:
        if control_index >= len(starts):
            continue
        center_index = int(starts[control_index])
        if center_index >= len(pts):
            continue
        s0 = svals[center_index]
        d = np.abs(svals - s0)
        d = np.minimum(d, total - d)
        mask = d < width
        if not np.any(mask):
            continue
        weight = np.zeros(len(pts), float)
        weight[mask] = 0.5 * (1.0 + np.cos(np.pi * d[mask] / width))
        sign = 1.0 if is_over else -1.0
        contribution = sign * float(sphere_offset) * weight
        replace = weight > strengths
        offsets[replace] = contribution[replace]
        strengths[replace] = weight[replace]
    return offsets


def _slerp_unit(a, b, t):
    """Spherical linear interpolation between two unit-vector directions."""
    a = _normalize_vector(a)
    b = _normalize_vector(b, fallback=a)
    tt = float(t)
    dot = max(-1.0, min(1.0, float(np.dot(a, b))))
    if dot > 0.9995:
        return _normalize_vector((1.0 - tt) * a + tt * b, fallback=a)
    if dot < -0.9995:
        ref = np.array([0.0, 0.0, 1.0], float)
        tangent = np.cross(a, ref)
        if float(np.linalg.norm(tangent)) < 1.0e-9:
            ref = np.array([1.0, 0.0, 0.0], float)
            tangent = np.cross(a, ref)
        tangent = _normalize_vector(tangent)
        angle = math.pi * tt
        return _normalize_vector(math.cos(angle) * a + math.sin(angle) * tangent, fallback=a)
    theta = math.acos(dot)
    sin_theta = math.sin(theta)
    q = (math.sin((1.0 - tt) * theta) / sin_theta) * a
    q += (math.sin(tt * theta) / sin_theta) * b
    return _normalize_vector(q, fallback=a)


def _resample_closed_unit_offset_curve(unit_dirs, radial_offsets, sphere_radius, xyz_spacing):
    """
    Resample a closed spherical/radial curve to near-uniform XYZ spacing.

    The curve is represented by unit directions and radial offsets. Distances are
    measured in final XYZ coordinates. Directions are interpolated by spherical
    linear interpolation and offsets are interpolated linearly; therefore, when
    crossing_offset = 0, all output points remain exactly on the base sphere.
    """
    dirs = _normalize_rows(unit_dirs)
    offsets = np.asarray(radial_offsets, float)
    if len(dirs) == 0:
        return np.zeros((0, 3), float)
    if len(offsets) != len(dirs):
        raise ValueError("Internal error: direction/offset arrays have different lengths.")
    spacing = float(xyz_spacing)
    if spacing <= 0.0:
        raise ValueError("--xyz-spacing must be positive.")

    if len(dirs) > 2:
        same_dir = float(np.linalg.norm(dirs[-1] - dirs[0])) < 1.0e-9
        same_off = abs(float(offsets[-1] - offsets[0])) < 1.0e-9
        if same_dir and same_off:
            dirs = dirs[:-1]
            offsets = offsets[:-1]

    if len(dirs) == 1:
        return dirs * (float(sphere_radius) + offsets[:1])[:, None]

    radii = float(sphere_radius) + offsets
    xyz = dirs * radii[:, None]
    seg_lengths = np.linalg.norm(np.roll(xyz, -1, axis=0) - xyz, axis=1)
    total = float(np.sum(seg_lengths))
    if total <= 1.0e-12:
        return xyz[:1].copy()

    # Use ceil so actual spacing is no larger than requested. The first point is
    # not repeated; write_spherical_xyz handles optional explicit closure.
    n_out = max(4, int(math.ceil(total / spacing)))
    n_out = min(n_out, 200000)
    actual = total / float(n_out)
    cumulative = np.concatenate(([0.0], np.cumsum(seg_lengths)))
    out = np.zeros((n_out, 3), float)

    for i in range(n_out):
        target = i * actual
        j = int(np.searchsorted(cumulative, target, side="right") - 1)
        j = max(0, min(j, len(seg_lengths) - 1))
        if seg_lengths[j] <= 1.0e-12:
            t = 0.0
        else:
            t = (target - cumulative[j]) / seg_lengths[j]
        j2 = (j + 1) % len(dirs)
        direction = _slerp_unit(dirs[j], dirs[j2], t)
        offset = (1.0 - t) * offsets[j] + t * offsets[j2]
        out[i] = direction * (float(sphere_radius) + offset)
    return out



def _smooth_periodic_array(values, window=4, passes=2):
    """Circular Gaussian smoothing for a closed sequence."""
    arr = np.asarray(values, float)
    if len(arr) == 0:
        return arr.copy()
    radius = int(round(float(window)))
    npasses = int(round(float(passes)))
    if radius <= 0 or npasses <= 0 or len(arr) < 5:
        return arr.copy()

    radius = min(radius, max(1, len(arr) // 3))
    shifts = np.arange(-radius, radius + 1)
    sigma = max(1.0e-6, 0.5 * float(radius))
    weights = np.exp(-0.5 * (shifts.astype(float) / sigma) ** 2)
    weights = weights / float(np.sum(weights))

    out = arr.copy()
    for _ in range(npasses):
        smoothed = np.zeros_like(out, dtype=float)
        for shift, weight in zip(shifts, weights):
            smoothed += float(weight) * np.roll(out, int(shift), axis=0)
        out = smoothed
    return out


def _preserve_signed_offset_extrema(original_offsets, smoothed_offsets):
    """Rescale positive/negative smoothed offsets so layer heights stay recognizable."""
    original = np.asarray(original_offsets, float)
    smoothed = np.asarray(smoothed_offsets, float).copy()
    if len(original) == 0 or len(smoothed) == 0:
        return smoothed

    old_pos = max(0.0, float(np.max(original)))
    new_pos = max(0.0, float(np.max(smoothed)))
    if old_pos > 1.0e-12 and new_pos > 1.0e-12:
        mask = smoothed > 0.0
        smoothed[mask] *= old_pos / new_pos

    old_neg = min(0.0, float(np.min(original)))
    new_neg = min(0.0, float(np.min(smoothed)))
    if old_neg < -1.0e-12 and new_neg < -1.0e-12:
        mask = smoothed < 0.0
        smoothed[mask] *= old_neg / new_neg

    # Never create larger excursions than the unsmoothed curve had.
    lo = min(float(np.min(original)), 0.0)
    hi = max(float(np.max(original)), 0.0)
    return np.clip(smoothed, lo, hi)


def _final_smooth_xyz_component(points, sphere_radius, xyz_spacing,
                                smooth_window=10, smooth_passes=5):
    """
    Smooth a closed XYZ component after the spherical construction.

    The smoothing is done in two channels: unit direction on the sphere and
    radial offset from the base sphere.  This removes visible kink artifacts at
    stitched local arcs while keeping crossing over/under layers recognizable.
    With zero crossing offset, all returned points remain on the base sphere.
    """
    pts = np.asarray(points, float)
    if len(pts) < 5:
        return pts.copy()

    radius = float(sphere_radius)
    if radius <= 0.0:
        return pts.copy()

    r = np.linalg.norm(pts, axis=1)
    safe = r > 1.0e-12
    dirs = np.zeros_like(pts, float)
    dirs[safe] = pts[safe] / r[safe, None]
    if not np.all(safe):
        dirs = _normalize_rows(dirs)
    offsets = r - radius

    smooth_dirs = _smooth_periodic_array(dirs, window=smooth_window, passes=smooth_passes)
    smooth_dirs = _normalize_rows(smooth_dirs)

    smooth_offsets = _smooth_periodic_array(
        offsets[:, None], window=smooth_window, passes=smooth_passes
    ).reshape(-1)
    smooth_offsets = _preserve_signed_offset_extrema(offsets, smooth_offsets)

    # Re-distribute points by the requested chord spacing after smoothing.  This
    # keeps the output visually even instead of merely smoothing the old vertices.
    return _resample_closed_unit_offset_curve(
        smooth_dirs,
        smooth_offsets,
        sphere_radius=radius,
        xyz_spacing=xyz_spacing,
    )


def _apply_final_xyz_smoothing(xyz_components, sphere_radius, xyz_spacing,
                               enabled=True, smooth_window=10, smooth_passes=5):
    """Apply the optional final 3D smoothing pass to every component."""
    if not enabled:
        return xyz_components
    window = int(round(float(smooth_window)))
    passes = int(round(float(smooth_passes)))
    if window <= 0 or passes <= 0:
        return xyz_components
    return [
        _final_smooth_xyz_component(
            arr,
            sphere_radius=sphere_radius,
            xyz_spacing=xyz_spacing,
            smooth_window=window,
            smooth_passes=passes,
        )
        for arr in xyz_components
    ]

def _samples_for_length(length, xyz_spacing, oversample=5.0, min_samples=5, max_samples=2001):
    """Internal smooth-curve sample count for a span of approximate XYZ length."""
    spacing = float(xyz_spacing)
    if spacing <= 0.0:
        raise ValueError("--xyz-spacing must be positive.")
    target = spacing / max(1.0, float(oversample))
    nsamp = int(math.ceil(max(0.0, float(length)) / target)) + 1
    return max(int(min_samples), min(int(max_samples), nsamp))


def _internal_samples_per_span_for_xyz(sphere_radius, sphere_offset, xyz_spacing):
    """Conservative hidden 2D sampling density before final even XYZ resampling."""
    max_radius = float(sphere_radius) + abs(float(sphere_offset))
    approx = math.pi * max(1.0, max_radius)
    return _samples_for_length(approx, xyz_spacing, oversample=4.0, min_samples=32, max_samples=700)


def _build_stereo_xyz_components(
    model,
    P,
    centers,
    sphere_radius,
    sphere_extent,
    sphere_offset,
    sphere_bump_frac,
    xyz_spacing,
):
    """Map the 2D diagram to the sphere, then evenly resample each XYZ component."""
    internal_samples = _internal_samples_per_span_for_xyz(sphere_radius, sphere_offset, xyz_spacing)
    component_samples = _sample_planar_components_for_sphere(
        model, P, centers, internal_samples
    )
    center, scale = _normalization_for_sphere(component_samples, sphere_extent)
    span = _diagram_span_from_samples(component_samples)
    bump_radius = float(sphere_bump_frac) * 0.04 * span

    xyz_components = []
    for item in component_samples:
        xy = item["xy"]
        if len(xy) == 0:
            xyz_components.append(np.zeros((0, 3), float))
            continue
        uv = (xy - center) * scale
        unit = _inverse_stereographic_to_unit_sphere(uv)
        radial_offsets = _crossing_radial_offsets(
            xy,
            item["marks"],
            item["starts"],
            bump_radius=bump_radius,
            sphere_offset=sphere_offset,
        )
        xyz_components.append(
            _resample_closed_unit_offset_curve(
                unit,
                radial_offsets,
                sphere_radius=sphere_radius,
                xyz_spacing=xyz_spacing,
            )
        )
    return xyz_components


def _fibonacci_sphere_directions(n):
    """Deterministic, well-spaced unit vectors used as a 3D layout seed."""
    n = int(n)
    if n <= 0:
        return np.zeros((0, 3), float)
    out = np.zeros((n, 3), float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / float(n)
        r = math.sqrt(max(0.0, 1.0 - z * z))
        th = i * golden
        out[i] = [math.cos(th) * r, math.sin(th) * r, z]
    return out


def _kamada_3d_unit_directions(G):
    """
    Kamada-Kawai-style 3D graph layout followed by radial projection to S^2.
    """
    nodes = sorted(list(G.nodes()), key=repr)
    if not nodes:
        return {}
    seed = _fibonacci_sphere_directions(len(nodes))
    pos0 = {node: seed[i] for i, node in enumerate(nodes)}
    try:
        raw = nx.kamada_kawai_layout(G, pos=pos0, dim=3)
    except Exception:
        raw = nx.spring_layout(G, pos=pos0, dim=3, seed=1, iterations=400)

    coords = np.array([np.asarray(raw[node], float) for node in nodes], float)
    coords = coords - coords.mean(axis=0)
    if float(np.max(np.linalg.norm(coords, axis=1))) <= 1.0e-10:
        coords = seed.copy()
    dirs = _normalize_rows(coords)
    return {node: dirs[i] for i, node in enumerate(nodes)}


def _tangent_basis_at(center):
    c = _normalize_vector(center)
    ref = np.array([0.0, 0.0, 1.0], float)
    if abs(float(np.dot(c, ref))) > 0.90:
        ref = np.array([1.0, 0.0, 0.0], float)
    e1 = ref - float(np.dot(ref, c)) * c
    e1 = _normalize_vector(e1, fallback=np.array([1.0, 0.0, 0.0], float))
    e2 = _normalize_vector(np.cross(c, e1), fallback=np.array([0.0, 1.0, 0.0], float))
    return e1, e2


def _project_tangent(v, center, fallback):
    c = _normalize_vector(center)
    t = np.asarray(v, float) - float(np.dot(v, c)) * c
    return _normalize_vector(t, fallback=fallback)


def _unit_from_tangent(center, tangent, angle):
    c = _normalize_vector(center)
    t = _project_tangent(tangent, c, _tangent_basis_at(c)[0])
    return _normalize_vector(math.cos(float(angle)) * c + math.sin(float(angle)) * t, fallback=c)


def _angle_between_unit_vectors(a, b):
    """Return the unsigned angular distance between two unit-vector directions."""
    aa = _normalize_vector(a)
    bb = _normalize_vector(b, fallback=aa)
    dot = float(np.clip(np.dot(aa, bb), -1.0, 1.0))
    return math.acos(dot)


def _compact_crossing_gadgets_on_sphere(model, raw_dirs, crossing_angle):
    """
    Create small local crossing patches on the sphere.

    For each true crossing, the four gadget corners are compacted around a single
    crossing center. Odd and even visits are placed on opposite tangent axes.
    """
    angle = float(crossing_angle)
    if angle <= 0.0:
        raise ValueError("Internal crossing angle must be positive.")
    if angle >= 0.50:
        raise ValueError("Internal crossing angle should be less than 0.5 radians.")

    roles = ("in_o", "in_e", "out_o", "out_e")
    fallback = _fibonacci_sphere_directions(max(1, len(model["crossings"])))
    P3 = dict(raw_dirs)
    centers3 = {}

    for k in range(len(model["crossings"])):
        raw_corners = [raw_dirs.get((k, r), fallback[k % len(fallback)]) for r in roles]
        center = _normalize_vector(np.mean(raw_corners, axis=0), fallback=fallback[k % len(fallback)])
        centers3[k] = center
        b1, b2 = _tangent_basis_at(center)

        e1_guess = raw_dirs.get((k, "out_o"), center) - raw_dirs.get((k, "in_o"), center)
        e1 = _project_tangent(e1_guess, center, b1)
        e2 = _normalize_vector(np.cross(center, e1), fallback=b2)
        e2_guess = _project_tangent(
            raw_dirs.get((k, "out_e"), center) - raw_dirs.get((k, "in_e"), center),
            center,
            e2,
        )
        if float(np.dot(e2, e2_guess)) < 0.0:
            e2 = -e2

        P3[(k, "in_o")] = _unit_from_tangent(center, -e1, angle)
        P3[(k, "out_o")] = _unit_from_tangent(center, e1, angle)
        P3[(k, "in_e")] = _unit_from_tangent(center, -e2, angle)
        P3[(k, "out_e")] = _unit_from_tangent(center, e2, angle)

    for cp in model["comp_positions"]:
        for p in cp:
            key = ("seg", p)
            if key not in P3:
                P3[key] = _normalize_vector(raw_dirs.get(key, [0.0, 0.0, 1.0]))
    return P3, centers3


def _spherical_quadratic_bezier(a, b, c, samples):
    """Sample a normalized quadratic Bezier arc on S^2."""
    nsamp = max(2, int(samples))
    a = _normalize_vector(a)
    b = _normalize_vector(b, fallback=a)
    c = _normalize_vector(c, fallback=b)
    out = []
    for i in range(nsamp):
        t = i / float(nsamp - 1)
        q = (1.0 - t) * (1.0 - t) * a + 2.0 * (1.0 - t) * t * b + t * t * c
        fallback = a if t < 0.5 else c
        out.append(_normalize_vector(q, fallback=fallback))
    return np.asarray(out, float)


def _sample_local_crossing_arc(center, tangent_axis, crossing_angle, bump_half_angle,
                               is_over, crossing_offset, samples):
    """
    Analytic local crossing arc on S^2 with a radial over/under bump.

    This avoids the small U-shaped wiggle that can occur when a global spline is
    radially displaced after the fact.
    """
    nsamp = max(5, int(samples))
    if nsamp % 2 == 0:
        nsamp += 1
    center = _normalize_vector(center)
    basis0 = _tangent_basis_at(center)[0]
    axis = _project_tangent(tangent_axis, center, basis0)
    angle = float(crossing_angle)
    half_width = max(1.0e-12, min(float(bump_half_angle), angle))
    sign = 1.0 if is_over else -1.0

    dirs = []
    offsets = []
    for i in range(nsamp):
        t = -1.0 + 2.0 * i / float(nsamp - 1)
        phi = t * angle
        q = math.cos(phi) * center + math.sin(phi) * axis
        dirs.append(_normalize_vector(q, fallback=center))
        ad = abs(phi)
        if crossing_offset != 0.0 and ad <= half_width:
            w = 0.5 * (1.0 + math.cos(math.pi * ad / half_width))
            offsets.append(sign * float(crossing_offset) * w)
        else:
            offsets.append(0.0)
    return np.asarray(dirs, float), np.asarray(offsets, float)


def _append_direction_samples(direction_list, offset_list, dirs, offsets, skip_first=True):
    """Append sampled unit directions and matching radial offsets."""
    if len(dirs) == 0:
        return
    start = 1 if skip_first and len(direction_list) > 0 else 0
    for i in range(start, len(dirs)):
        direction_list.append(_normalize_vector(dirs[i]))
        offset_list.append(float(offsets[i]))


def _smooth_layer_offsets(start_offset, end_offset, samples, transition_fraction=1.0):
    """
    Radial offsets for direct connecting between two crossing layers.

    If the two endpoint layers match, the whole connector stays on that layer.
    If they differ, the offset changes with a C1 smoothstep profile.  A
    transition_fraction of 1 uses the full connector; smaller values keep the
    first/last layers longer and transition only in the middle.
    """
    nsamp = max(1, int(samples))
    start = float(start_offset)
    end = float(end_offset)
    if nsamp == 1:
        return np.asarray([start], float)
    if abs(start - end) <= 1.0e-12:
        return np.full(nsamp, start, float)

    frac = max(1.0e-6, min(1.0, float(transition_fraction)))
    lo = 0.5 * (1.0 - frac)
    hi = 0.5 * (1.0 + frac)
    vals = np.zeros(nsamp, float)
    for i in range(nsamp):
        t = i / float(nsamp - 1)
        if t <= lo:
            u = 0.0
        elif t >= hi:
            u = 1.0
        else:
            u = (t - lo) / max(1.0e-12, hi - lo)
            u = u * u * (3.0 - 2.0 * u)
        vals[i] = (1.0 - u) * start + u * end
    return vals


def _layer_offset_for_visit(model, pos, sphere_offset):
    """Return +offset for an over/bump visit and -offset for an under/dip visit."""
    return (1.0 if model["over_at"][pos] else -1.0) * float(sphere_offset)


def _build_spherical_kamada_xyz_components(
    model,
    G,
    sphere_radius,
    sphere_offset,
    sphere_bump_frac,
    xyz_spacing,
    crossing_angle,
    direct_connecting=True,
    skeleton_out=None,
    skeleton_grid=24,
):
    """Sphere-native method: distribute the diagram directly on S^2."""
    if G is None:
        G = build_gadget_graph(model)
    raw_dirs = _kamada_3d_unit_directions(G)
    P3, centers3 = _compact_crossing_gadgets_on_sphere(
        model, raw_dirs, crossing_angle=crossing_angle
    )

    # V5.0 'mapped skeleton': the per-crossing anchors and the strand chords
    # between them (through the connector midpoints), as UNIT directions; the
    # caller scales them to the sphere radius and applies the surface warp.
    if skeleton_out is not None:
        skeleton_out.clear()
        _ncross = len(model["crossings"])
        skeleton_out["points"] = np.array(
            [centers3[k] for k in range(_ncross)], float
        )
        _segs = []
        _seg_comps = []
        for _ci, cp in enumerate(model["comp_positions"]):
            for p in cp:
                k = model["pos_cross"][p]
                q = model["nextpos"][p]
                k2 = model["pos_cross"][q]
                _segs.append(np.array(
                    [centers3[k], P3[("seg", p)], centers3[k2]], float
                ))
                _seg_comps.append(_ci)
        skeleton_out["segments"] = _segs
        skeleton_out["segment_comps"] = _seg_comps
        # Whole-surface wireframe (unit sphere grid; the caller scales it and
        # puts it through the same surface warp as the strands, so the drawn
        # mesh IS the mapped surface -- sphere, ellipsoid, cylinder or torus).
        _wires = []
        _ng = max(0, int(skeleton_grid or 0))
        if _ng > 0:
            _npar = max(3, _ng // 2 - 1)
            _vv = np.linspace(0.02 * math.pi, 0.98 * math.pi, 49)
            for _u in np.linspace(0.0, 2.0 * math.pi, _ng + 1)[:-1]:
                _wires.append(np.column_stack([
                    np.cos(_u) * np.sin(_vv), np.sin(_u) * np.sin(_vv), np.cos(_vv)
                ]))
            _uu = np.linspace(0.0, 2.0 * math.pi, 73)
            for _v in np.linspace(0.0, math.pi, _npar + 2)[1:-1]:
                _wires.append(np.column_stack([
                    np.cos(_uu) * math.sin(_v), np.sin(_uu) * math.sin(_v),
                    np.full_like(_uu, math.cos(_v))
                ]))
        skeleton_out["surface_wires"] = _wires
        skeleton_out["unit"] = True

    dense_step = float(xyz_spacing) / 5.0
    if dense_step <= 0.0:
        raise ValueError("--xyz-spacing must be positive.")
    effective_radius = float(sphere_radius) + abs(float(sphere_offset))
    cross_len = max(1.0e-12, 2.0 * float(crossing_angle) * effective_radius)
    nspan_cross = max(11, int(math.ceil(cross_len / dense_step)) + 1)
    if nspan_cross % 2 == 0:
        nspan_cross += 1

    bump_half_angle = float(crossing_angle) * float(sphere_bump_frac)

    xyz_components = []
    for ci, cp in enumerate(model["comp_positions"]):
        directions = []
        radial_offsets = []
        for p in cp:
            k = model["pos_cross"][p]
            role = model["pos_role"][p]
            center = centers3[k]
            in_dir = P3[(k, "in_" + role)]
            out_dir = P3[(k, "out_" + role)]

            axis = _project_tangent(out_dir - in_dir, center, _tangent_basis_at(center)[0])
            in_proj = _project_tangent(in_dir, center, -axis)
            if float(np.dot(in_proj, axis)) > 0.0:
                axis = -axis

            cross_dirs, cross_offsets = _sample_local_crossing_arc(
                center=center,
                tangent_axis=axis,
                crossing_angle=crossing_angle,
                bump_half_angle=bump_half_angle,
                is_over=model["over_at"][p],
                crossing_offset=sphere_offset,
                samples=nspan_cross,
            )
            if direct_connecting:
                # In direct-connecting mode the local crossing visit itself lies
                # on the over/under layer from its incoming endpoint through its
                # outgoing endpoint.  This lets the following connector start
                # directly from the same layer rather than dropping back to the
                # middle sphere.
                cross_offsets = np.full(
                    len(cross_dirs),
                    _layer_offset_for_visit(model, p, sphere_offset),
                    float,
                )
            _append_direction_samples(directions, radial_offsets, cross_dirs, cross_offsets)

            q = model["nextpos"][p]
            next_k = model["pos_cross"][q]
            next_role = model["pos_role"][q]
            next_in = P3[(next_k, "in_" + next_role)]
            mid = P3[("seg", p)]

            conn_angle = _angle_between_unit_vectors(out_dir, mid) + _angle_between_unit_vectors(mid, next_in)
            conn_len = max(1.0e-12, conn_angle * effective_radius)
            nspan_conn = max(5, int(math.ceil(conn_len / dense_step)) + 1)
            conn_dirs = _spherical_quadratic_bezier(out_dir, mid, next_in, nspan_conn)
            if direct_connecting:
                conn_offsets = _smooth_layer_offsets(
                    _layer_offset_for_visit(model, p, sphere_offset),
                    _layer_offset_for_visit(model, q, sphere_offset),
                    len(conn_dirs),
                    transition_fraction=sphere_bump_frac,
                )
            else:
                conn_offsets = np.zeros(len(conn_dirs), float)
            _append_direction_samples(directions, radial_offsets, conn_dirs, conn_offsets)

        if not directions:
            xyz_components.append(np.zeros((0, 3), float))
            continue
        xyz_components.append(
            _resample_closed_unit_offset_curve(
                np.asarray(directions, float),
                np.asarray(radial_offsets, float),
                sphere_radius=sphere_radius,
                xyz_spacing=xyz_spacing,
            )
        )
    return xyz_components


SURFACE_SHAPES = ("ellipsoid", "cylinder", "torus")


def _normalize_surface_shape(shape):
    s = str(shape or "ellipsoid").strip().lower()
    return s if s in SURFACE_SHAPES else "ellipsoid"


def _pca_frame_3d(dirs):
    """
    Principal axes of a 3D direction cloud.

    Returns (evecs, evals) with columns ordered by descending eigenvalue, so
    evecs[:, 0] is the major axis and evecs[:, 2] is the minor axis.
    """
    pts = np.asarray(dirs, float)
    if len(pts) < 3:
        return np.eye(3), np.ones(3)
    q = pts - pts.mean(axis=0)
    cov = q.T @ q
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    evals = np.maximum(evals[order], 0.0)
    evecs = evecs[:, order]
    if np.linalg.det(evecs) < 0.0:
        evecs[:, 2] = -evecs[:, 2]
    return evecs, evals


def _rot_about_canonical_axis(axis, angle):
    """3x3 rotation about a canonical axis ('x', 'y', or 'z')."""
    c = math.cos(float(angle))
    s = math.sin(float(angle))
    if axis == "x":
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], float)
    if axis == "y":
        return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], float)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], float)


def _surface_frame_and_params(kind, dirs_all, aspect, tube,
                              orient_auto=True, aspect_auto=True,
                              orient_degrees=0.0, tilt_degrees=0.0):
    """
    Resolve the surface's orientation frame ``M`` and its shape parameters.

    Orientation and aspect are handled independently.

    Orientation: the surface's primary axis (ellipsoid major -> +x; cylinder
    length and torus symmetry -> +z) is always aligned to the diagram's principal
    (PCA) axis.  When ``orient_auto`` is off, two manual rotations are applied
    relative to that PCA axis: ``orient_degrees`` spins about the primary axis,
    and ``tilt_degrees`` tilts the primary axis away from the PCA axis (a
    rotation about the perpendicular secondary axis).

    Aspect: with ``aspect_auto`` on, the ellipsoid axis ratios / cylinder height
    / torus tube ratio are derived from the diagram's 3D PCA magnitudes so the
    surface proportions match the diagram's own elongation.  Otherwise the
    supplied ``aspect`` / ``tube`` are used.
    """
    kind = _normalize_surface_shape(kind)
    aspect = float(aspect)
    if not np.isfinite(aspect) or aspect <= 0.0:
        aspect = 1.0
    tube = float(tube)
    if not np.isfinite(tube) or tube <= 0.0:
        tube = 0.35

    pts = np.asarray(dirs_all, float)
    have_pca = pts.shape[0] >= 3
    if have_pca:
        evecs, evals = _pca_frame_3d(pts)
        e0, e1, e2 = evecs[:, 0], evecs[:, 1], evecs[:, 2]
        ev = np.maximum(evals, 1.0e-9)
        gm = float(np.prod(ev) ** (1.0 / 3.0))

    # Base PCA alignment of the surface primary axis, plus the canonical primary
    # axis (for the manual spin) and a perpendicular secondary axis (for the tilt
    # away from the PCA axis).
    if kind == "cylinder":
        primary, secondary = "z", "x"
        M0 = np.column_stack([e2, e1, e0]) if have_pca else np.eye(3)
    elif kind == "torus":
        primary, secondary = "z", "x"
        M0 = np.column_stack([e0, e1, e2]) if have_pca else np.eye(3)
    else:  # ellipsoid
        primary, secondary = "x", "y"
        M0 = np.column_stack([e0, e1, e2]) if have_pca else np.eye(3)

    if orient_auto:
        M = M0
    else:
        M = (M0
             @ _rot_about_canonical_axis(primary, math.radians(orient_degrees))
             @ _rot_about_canonical_axis(secondary, math.radians(tilt_degrees)))

    if aspect_auto and have_pca:
        if kind == "ellipsoid":
            resolved = {"axes": np.clip(np.sqrt(ev / gm), 0.4, 2.6)}
        elif kind == "cylinder":
            height = float(np.clip(math.sqrt(ev[0] / (0.5 * (ev[1] + ev[2]))), 1.0, 6.0))
            resolved = {"radius": 1.0, "height": height}
        else:  # torus
            tr = float(np.clip(math.sqrt(ev[2] / (0.5 * (ev[0] + ev[1]))), 0.12, 0.55))
            resolved = {"tube": tr}
    else:
        if kind == "ellipsoid":
            axes = np.array([aspect, 1.0, 1.0 / max(aspect, 1.0e-6)], float)
            gm2 = float(np.prod(axes) ** (1.0 / 3.0))
            resolved = {"axes": axes / (gm2 if gm2 > 1.0e-12 else 1.0)}
        elif kind == "cylinder":
            resolved = {"radius": 1.0, "height": aspect}
        else:  # torus
            resolved = {"tube": tube}
    return M, resolved


def _surface_warp_fn(kind, sphere_radius, resolved, M):
    """
    Build warp(dirs, offsets) -> Nx3 mapping unit-sphere directions and signed
    normal offsets onto the chosen surface.  Over/under offsets are applied along
    the local outward surface normal, so with offset 0 all points lie exactly on
    the base surface.  For the sphere this is the identity map used elsewhere.
    """
    kind = _normalize_surface_shape(kind)
    R = float(sphere_radius)
    Mrot = np.asarray(M, float)

    def warp(dirs, offsets):
        d = np.asarray(dirs, float)
        if d.size == 0:
            return np.zeros((0, 3), float)
        d = d @ Mrot
        x = d[:, 0]
        y = d[:, 1]
        z = np.clip(d[:, 2], -1.0, 1.0)
        o = np.asarray(offsets, float)
        if kind == "ellipsoid":
            ax, ay, az = resolved["axes"]
            base = np.column_stack([R * ax * x, R * ay * y, R * az * z])
            nn = np.column_stack([x / ax, y / ay, z / az])
        elif kind == "cylinder":
            rad = R * float(resolved["radius"])
            H = R * float(resolved["height"])
            phi = np.arctan2(y, x)
            base = np.column_stack([rad * np.cos(phi), rad * np.sin(phi), H * z])
            nn = np.column_stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)])
        elif kind == "torus":
            r_major = R
            r_minor = float(resolved["tube"]) * R
            u = np.arctan2(y, x)
            theta = np.arccos(z)               # 0..pi from north to south pole
            v = math.pi - 2.0 * theta          # full poloidal sweep, poles -> inner ridge
            ring = r_major + r_minor * np.cos(v)
            base = np.column_stack([ring * np.cos(u), ring * np.sin(u), r_minor * np.sin(v)])
            nn = np.column_stack([np.cos(v) * np.cos(u), np.cos(v) * np.sin(u), np.sin(v)])
        else:  # sphere
            base = R * d
            nn = d
        norm = np.linalg.norm(nn, axis=1)
        norm[norm < 1.0e-12] = 1.0
        nn = nn / norm[:, None]
        return base + o[:, None] * nn

    return warp


def _geom_resample_closed_3d(points, spacing):
    """Even arc-length resample of a closed 3D polyline (surface-agnostic)."""
    pts = np.asarray(points, float)
    spacing = float(spacing)
    if len(pts) < 4 or spacing <= 0.0:
        return pts
    seg = np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1)
    total = float(np.sum(seg))
    if total <= 1.0e-12:
        return pts[:1].copy()
    n_out = max(4, int(math.ceil(total / spacing)))
    n_out = min(n_out, 200000)
    actual = total / float(n_out)
    cum = np.concatenate(([0.0], np.cumsum(seg)))
    out = np.zeros((n_out, 3), float)
    for i in range(n_out):
        target = i * actual
        j = int(np.searchsorted(cum, target, side="right") - 1)
        j = max(0, min(j, len(seg) - 1))
        t = 0.0 if seg[j] <= 1.0e-12 else (target - cum[j]) / seg[j]
        j2 = (j + 1) % len(pts)
        out[i] = (1.0 - t) * pts[j] + t * pts[j2]
    return out


def _warp_components_to_surface(xyz_components, sphere_radius, surface_shape,
                                surface_aspect, surface_tube, xyz_spacing,
                                surface_auto_orient=True, surface_auto_aspect=True,
                                surface_orient=0.0, surface_tilt=0.0,
                                return_warp=False):
    """
    Warp spherical components (built on the base sphere of radius R) onto a
    shaped surface.  Each spherical point is decomposed into a unit direction and
    a signed radial offset, a common orientation frame + shape parameters are
    resolved once from all components' directions, and every point is remapped
    and then re-spaced evenly along the surface.
    """
    R = float(sphere_radius)
    decomposed = []
    all_dirs = []
    for arr in xyz_components:
        arr = np.asarray(arr, float)
        if len(arr) == 0:
            decomposed.append(None)
            continue
        r = np.linalg.norm(arr, axis=1)
        safe = r > 1.0e-12
        d = np.zeros_like(arr)
        d[safe] = arr[safe] / r[safe, None]
        if not np.all(safe):
            d = _normalize_rows(d)
        o = r - R
        decomposed.append((d, o))
        all_dirs.append(d)

    dirs_all = np.vstack(all_dirs) if all_dirs else np.zeros((0, 3), float)
    M, resolved = _surface_frame_and_params(
        surface_shape, dirs_all, surface_aspect, surface_tube,
        orient_auto=surface_auto_orient,
        aspect_auto=surface_auto_aspect,
        orient_degrees=surface_orient,
        tilt_degrees=surface_tilt,
    )
    warp = _surface_warp_fn(surface_shape, R, resolved, M)

    out = []
    for item in decomposed:
        if item is None:
            out.append(np.zeros((0, 3), float))
            continue
        d, o = item
        warped = warp(d, o)
        out.append(_geom_resample_closed_3d(warped, xyz_spacing))

    if return_warp:
        def _warp_pointset(arr):
            arr = np.asarray(arr, float)
            if len(arr) == 0:
                return arr
            r = np.linalg.norm(arr, axis=1)
            safe = r > 1.0e-12
            d = np.zeros_like(arr)
            d[safe] = arr[safe] / r[safe, None]
            if not np.all(safe):
                d = _normalize_rows(d)
            return warp(d, r - R)
        return out, _warp_pointset
    return out




# --------------------------------------------------------------------------
# V5.4: 3D clearance repair + topology audit.
#
# Dense diagrams (chains, heavy chains) cannot always be laid out on a single
# surface layer: the spherical layout can pinch two strands to near-contact
# away from any true crossing, and smoothing/warping can collapse the radial
# over/under separation at true crossings.  Both corrupt the 3D curve's link
# type even though the 2D diagram is audit-clean.  The repair pass below
# detects every close approach between non-neighboring arcs and pushes the
# arcs apart -- coherently per pinch cluster, preserving the current local
# order -- effectively stacking crowded regions into extra radial layers.
# It runs after the raw build, after smoothing, and after the surface warp.
# When spherogram/snappy are available (audit_xyz.py), a final topology
# audit reconstructs the link from the 3D curve and verifies it is
# isometric to the DT link, so any residual corruption is reported loudly.
# --------------------------------------------------------------------------

def _xyz_conflict_pairs(comps, clearance, skip=4, anchors=None,
                        exclude_radius=0.0):
    """Close approaches between non-neighboring points of closed polylines.

    Returns a list of (ci, i, cj, j, dist).  Point sampling is dense
    (xyz_spacing), so point-point distance is an adequate proxy for
    segment-segment distance at these clearances."""
    cell = max(float(clearance), 1.0e-9)
    grid = {}
    for ci, pts in enumerate(comps):
        for i, p in enumerate(pts):
            key = (int(math.floor(p[0] / cell)),
                   int(math.floor(p[1] / cell)),
                   int(math.floor(p[2] / cell)))
            grid.setdefault(key, []).append((ci, i))
    out = []
    for ci, pts in enumerate(comps):
        n = len(pts)
        for i, p in enumerate(pts):
            base = (int(math.floor(p[0] / cell)),
                    int(math.floor(p[1] / cell)),
                    int(math.floor(p[2] / cell)))
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for cj, j in grid.get((base[0] + dx, base[1] + dy,
                                               base[2] + dz), ()):
                            if cj < ci or (cj == ci and j <= i):
                                continue
                            if cj == ci:
                                d_idx = min(abs(i - j), n - abs(i - j))
                                if d_idx <= skip:
                                    continue
                            q = comps[cj][j]
                            dist = math.sqrt((p[0] - q[0]) ** 2 +
                                             (p[1] - q[1]) ** 2 +
                                             (p[2] - q[2]) ** 2)
                            if dist < clearance:
                                if anchors is not None and exclude_radius > 0.0:
                                    mid = 0.5 * (np.asarray(p) + np.asarray(q))
                                    dmin = np.min(np.linalg.norm(
                                        anchors - mid[None, :], axis=1))
                                    if dmin < exclude_radius:
                                        continue   # inside a true crossing
                                out.append((ci, i, cj, j, dist))
    return out


def _repair_xyz_clearance(comps, clearance, max_iters=80, halfwidth=8,
                          push_frac=0.6, anchors=None, exclude_radius=0.0):
    """Separate all conflicting arcs to >= clearance, coherently per pinch.

    Conflicts are clustered by component pair and index adjacency; each
    cluster is pushed apart along its mean separation direction (preserving
    the current order, so true crossings -- already separated by the
    over/under bumps -- are never reordered, and each false pinch resolves
    to a trivial R2 pair on one consistent side).  Returns (comps, report).
    """
    comps = [np.array(c, float) for c in comps]
    report = {"initial": None, "iters": 0, "resolved": True,
              "max_push": 0.0, "layers": 1}
    total_push = [np.zeros(len(c)) for c in comps]
    for it in range(max_iters):
        conflicts = _xyz_conflict_pairs(comps, clearance, anchors=anchors,
                                        exclude_radius=exclude_radius)
        if report["initial"] is None:
            report["initial"] = len(conflicts)
        if not conflicts:
            report["iters"] = it
            break
        # cluster by (ci, cj) and adjacency of i
        conflicts.sort(key=lambda r: (r[0], r[2], r[1], r[3]))
        clusters = []
        for rec in conflicts:
            placed = False
            for cl in clusters:
                if cl[0][0] == rec[0] and cl[0][2] == rec[2] and                    abs(cl[-1][1] - rec[1]) <= 2 * halfwidth:
                    cl.append(rec)
                    placed = True
                    break
            if not placed:
                clusters.append([rec])
        disp = [np.zeros_like(c) for c in comps]
        wsum = [np.zeros(len(c)) for c in comps]
        gain = push_frac * (1.0 + 0.15 * it)      # escalate if stuck
        for cl in clusters:
            ci, cj = cl[0][0], cl[0][2]
            # coherent mean separation direction for the whole pinch
            vecs = [comps[ci][i] - comps[cj][j] for (_, i, _, j, _) in cl]
            mean = np.sum(vecs, axis=0)
            norm = np.linalg.norm(mean)
            if norm < 1.0e-9:
                mean = comps[ci][cl[0][1]].copy()
                norm = np.linalg.norm(mean)
                if norm < 1.0e-9:
                    mean, norm = np.array([0.0, 0.0, 1.0]), 1.0
            u = mean / norm
            need = max(clearance - d for (_, _, _, _, d) in cl)
            step = 0.5 * gain * need
            for (_, i, _, j, _) in cl:
                for comp_idx, idx, sgn in ((ci, i, 1.0), (cj, j, -1.0)):
                    n = len(comps[comp_idx])
                    for k in range(-halfwidth, halfwidth + 1):
                        w = math.exp(-(k * k) / (0.5 * halfwidth ** 2))
                        m = (idx + k) % n
                        disp[comp_idx][m] += sgn * step * w * u
                        wsum[comp_idx][m] += w
        for c_idx in range(len(comps)):
            mask = wsum[c_idx] > 0
            if mask.any():
                comps[c_idx][mask] += (disp[c_idx][mask] /
                                       wsum[c_idx][mask][:, None])
                total_push[c_idx][mask] += np.linalg.norm(
                    disp[c_idx][mask] / wsum[c_idx][mask][:, None], axis=1)
    else:
        report["resolved"] = False
        report["iters"] = max_iters
    mp = max((float(t.max()) if len(t) else 0.0) for t in total_push)
    report["max_push"] = mp
    report["layers"] = 1 + int(math.ceil(mp / max(2.0 * clearance, 1e-9)))
    return comps, report






def _match_visit_sequence(runs, want_ids):
    """Order-aware matching of a component's clustered Voronoi runs
    [(id, i0, j0, dmin), ...] (cyclic) against the diagram's expected
    crossing-visit sequence want_ids (cyclic), minimizing the summed
    closest-approach distance.  Tries both orientations and all cyclic
    starts.  Returns the matched runs in want order, or None."""
    M, m = len(runs), len(want_ids)
    if M < m:
        return None
    best = None
    for orient in (1, -1):
        seq = runs if orient == 1 else runs[::-1]
        starts = [s for s in range(M) if seq[s][0] == want_ids[0]]
        for s in starts:
            INF = float("inf")
            cost = [[INF] * M for _ in range(m)]
            parent = [[-1] * M for _ in range(m)]
            cost[0][0] = seq[s][3]
            for t in range(1, m):
                run_min, run_arg = INF, -1
                for p in range(t - 1, M - (m - t)):
                    if cost[t - 1][p] < run_min:
                        run_min, run_arg = cost[t - 1][p], p
                    q = p + 1
                    if q < M and seq[(s + q) % M][0] == want_ids[t]:
                        c = run_min + seq[(s + q) % M][3]
                        if c < cost[t][q]:
                            cost[t][q] = c
                            parent[t][q] = run_arg
                # also allow q beyond consecutive: fill remaining columns
                for q in range(t, M):
                    if seq[(s + q) % M][0] != want_ids[t]:
                        continue
                    run_min2 = min(cost[t - 1][p] for p in range(0, q))
                    if run_min2 < INF:
                        c = run_min2 + seq[(s + q) % M][3]
                        if c < cost[t][q]:
                            cost[t][q] = c
                            parent[t][q] = int(np.argmin(
                                [cost[t - 1][p] for p in range(0, q)]))
            endc = min(cost[m - 1])
            if endc < INF and (best is None or endc < best[0]):
                q = int(np.argmin(cost[m - 1]))
                idxs = [0] * m
                idxs[m - 1] = q
                for t in range(m - 1, 0, -1):
                    q = parent[t][q]
                    idxs[t - 1] = q
                matched = [seq[(s + q) % M] for q in idxs]
                best = (endc, matched)
    return None if best is None else best[1]






def _match_all_visits(comps, model, anchors, ball_radius, messages):
    """Recover, for every component, its ordered crossing visits on the 3D
    curve: [(anchor_k, passage_point_index, over_flag), ...] per component,
    aligned with the diagram via order-aware minimum-distance matching.
    Returns None (with a warning) if alignment fails."""
    anchors = np.asarray(anchors, float)
    visits_by_comp = []
    for pts in comps:
        n = len(pts)
        d = np.linalg.norm(pts[:, None, :] - anchors[None, :, :], axis=2)
        lab = np.argmin(d, axis=1)
        runs = []
        i = 0
        while i < n:
            j = i
            while j + 1 < n and lab[j + 1] == lab[i]:
                j += 1
            runs.append([int(lab[i]), i, j])
            i = j + 1
        merged = []
        for r in runs:
            if merged and merged[-1][0] == r[0]:
                merged[-1][2] = r[2]
            else:
                merged.append(list(r))
        if len(merged) > 1 and merged[0][0] == merged[-1][0]:
            merged[0] = [merged[0][0], merged[-1][1] - n, merged[0][2]]
            merged.pop()
        gap_pts = max(6, int(round(ball_radius / max(
            1.0e-9, float(np.median(np.linalg.norm(
                np.diff(pts, axis=0), axis=1)))))))
        clustered = []
        last_of = {}
        for kk, i0, j0 in merged:
            if kk in last_of and i0 - clustered[last_of[kk]][2] <= gap_pts:
                clustered[last_of[kk]][2] = j0
            else:
                clustered.append([kk, i0, j0])
                last_of[kk] = len(clustered) - 1
        if len(clustered) > 1 and clustered[0][0] == clustered[-1][0] and \
                (clustered[0][1] + n - clustered[-1][2]) <= gap_pts:
            clustered[0] = [clustered[0][0],
                            clustered[-1][1] - n, clustered[0][2]]
            clustered.pop()
        out_runs = []
        for kk, i0, j0 in clustered:
            idxs = [t % n for t in range(i0, j0 + 1)]
            dmin = float(np.min(d[idxs, kk]))
            out_runs.append((kk, i0, j0, dmin))
        visits_by_comp.append(out_runs)
    result = []
    for ci, runs in enumerate(visits_by_comp):
        pts = comps[ci]
        n = len(pts)
        want = [model["pos_cross"][p] for p in model["comp_positions"][ci]]
        overs = [bool(model["over_at"][p])
                 for p in model["comp_positions"][ci]]
        matched = _match_visit_sequence(runs, want)
        if matched is None:
            from collections import Counter as _C
            have = _C(r[0] for r in runs)
            need = _C(want)
            new_runs = []
            for (kk, i0, j0, dmin) in runs:
                if have[kk] < need.get(kk, 0):
                    idxs = [t % n for t in range(i0, j0 + 1)]
                    prof = np.linalg.norm(
                        pts[idxs] - anchors[kk][None, :], axis=1)
                    if len(prof) >= 5:
                        a_ = int(np.argmin(prof))
                        prof2 = prof.copy()
                        prof2[max(0, a_ - 2):a_ + 3] = np.inf
                        b_ = int(np.argmin(prof2))
                        lo_, hi_ = sorted((a_, b_))
                        cut = lo_ + int(np.argmax(prof[lo_:hi_ + 1]))
                        new_runs.append((kk, i0, i0 + cut,
                                         float(prof[:cut + 1].min())))
                        new_runs.append((kk, i0 + cut + 1, j0,
                                         float(prof[cut + 1:].min())))
                        continue
                new_runs.append((kk, i0, j0, dmin))
            matched = _match_visit_sequence(new_runs, want)
        if matched is None:
            messages.append("[WARN] xyz visit matching failed for "
                            "component %d; layered rebuild skipped." % ci)
            return None
        vis = []
        for t, k in enumerate(want):
            kk, i0, j0 = matched[t][:3]
            idxs = [tt % n for tt in range(i0, j0 + 1)]
            dd = np.linalg.norm(pts[idxs] - anchors[k][None, :], axis=1)
            vis.append((k, idxs[int(np.argmin(dd))], overs[t]))
        result.append(vis)
    return result


def _layered_radial_rebuild(comps, model, anchors, sphere_radius,
                            sphere_offset, ball_radius, messages):
    """V5.4 principled fix: rebuild the entire radial structure from the
    diagram.  Complex diagrams cannot live on ONE surface layer: first
    compute how many layers are necessary, then render.

    Every inter-crossing arc gets a radial LAYER such that (i) arcs whose
    2D routes come close sit on different layers and (ii) arcs passing
    through a foreign crossing's gadget zone avoid the center layer (no
    silent threading between that crossing's over/under strands).  At each
    crossing passage the radius is exactly anchor_radius +/- sphere_offset
    per the DT over/under, blended smoothly into the arc layers.
    Returns (comps, n_layers) or (comps, -1) if visit matching failed."""
    comps = [np.array(c, float) for c in comps]
    anchors = np.asarray(anchors, float)
    ncross = len(anchors)
    visits = _match_all_visits(comps, model, anchors, ball_radius, messages)
    if visits is None:
        return comps, -1
    # ---- arcs: (ci, passage index v) spans passages v -> v+1 (cyclic)
    arcs = []            # (ci, t_start, t_end)  point-index span (cyclic)
    arc_of = {}
    for ci, vis in enumerate(visits):
        n = len(comps[ci])
        m = len(vis)
        for v in range(m):
            t0 = vis[v][1]
            t1 = vis[(v + 1) % m][1]
            arcs.append((ci, t0, t1 if t1 > t0 else t1 + n))
            arc_of[(ci, v)] = len(arcs) - 1
    # ---- conflict graph over arcs (angular proximity away from gadgets)
    R = float(sphere_radius)
    thr = 0.85 * float(ball_radius)
    pts_dir = [c / np.maximum(np.linalg.norm(c, axis=1)[:, None], 1e-12)
               for c in comps]
    a_dir = anchors / np.maximum(
        np.linalg.norm(anchors, axis=1)[:, None], 1e-12)
    samples = []
    for ai, (ci, t0, t1) in enumerate(arcs):
        n = len(comps[ci])
        idxs = [t % n for t in range(t0, t1 + 1, 3)]
        D = pts_dir[ci][idxs] * R
        da = np.linalg.norm(D[:, None, :] - a_dir[None, :, :] * R, axis=2)
        far = np.min(da, axis=1) > 1.1 * float(ball_radius)
        samples.append((D, far, np.min(da, axis=1)))
    conflict = [set() for _ in arcs]
    threader = [False] * len(arcs)
    for ai in range(len(arcs)):
        Da, fara, dmina = samples[ai]
        # does this arc cut through a FOREIGN gadget core?  (an arc may
        # legitimately enter only its own two endpoint gadgets)
        ci, t0, t1 = arcs[ai]
        vlist = visits[ci]
        vv = [v for v in range(len(vlist)) if arc_of[(ci, v)] == ai]
        own = set()
        for v in vv:
            own.add(vlist[v][0])
            own.add(vlist[(v + 1) % len(vlist)][0])
        da_all = np.linalg.norm(Da[:, None, :] - a_dir[None, :, :] * R,
                                axis=2)
        # ignore samples close to the arc's endpoints (inside own gadgets)
        for k in range(ncross):
            if k in own:
                continue
            if np.any(da_all[:, k] < 0.7 * float(ball_radius)):
                threader[ai] = True
                break
        for bj in range(ai + 1, len(arcs)):
            cj = arcs[bj][0]
            Db, farb, _ = samples[bj]
            sel_a = fara
            if not sel_a.any():
                continue
            dd = np.linalg.norm(Da[sel_a][:, None, :] - Db[None, :, :],
                                axis=2)
            mask_b = farb[None, :] if farb.any() else None
            if mask_b is None:
                continue
            dd = np.where(np.broadcast_to(farb[None, :], dd.shape), dd,
                          np.inf)
            if np.min(dd) < thr:
                conflict[ai].add(bj)
                conflict[bj].add(ai)
    # ---- greedy coloring (highest degree first)
    order = sorted(range(len(arcs)), key=lambda a: -len(conflict[a]))
    layer = [-1] * len(arcs)
    for a in order:
        used = {layer[b] for b in conflict[a] if layer[b] >= 0}
        l = 0
        while l in used:
            l += 1
        layer[a] = l
    L = max(layer) + 1
    # threading arcs must not sit on the center layer (radius ~ R)
    if L == 1 and any(threader):
        L = 2
    center = (L - 1) / 2.0
    for a in range(len(arcs)):
        if threader[a] and abs(layer[a] - center) < 0.25:
            used = {layer[b] for b in conflict[a] if layer[b] >= 0}
            l = 0
            while l in used or abs(l - center) < 0.25:
                l += 1
            layer[a] = l
            L = max(L, l + 1)
            center = (L - 1) / 2.0
    lam = 2.6 * float(sphere_offset)
    # ---- rebuild radii: exact passage radii at crossings, layer radii
    # mid-arc, cosine blends along each arc's OWN arclength; short arcs
    # interpolate passage-to-passage directly.
    blend = max(1.0e-9, float(ball_radius))
    for ci, vis in enumerate(visits):
        pts = comps[ci]
        n = len(pts)
        m = len(vis)
        new_r = np.zeros(n)
        assigned = np.zeros(n, bool)
        for v in range(m):
            a_idx = arc_of[(ci, v)]
            t0, t1 = arcs[a_idx][1], arcs[a_idx][2]
            r_layer = R + (layer[a_idx] - (L - 1) / 2.0) * lam
            k0, ov0 = vis[v][0], vis[v][2]
            k1, ov1 = vis[(v + 1) % m][0], vis[(v + 1) % m][2]
            r_p0 = np.linalg.norm(anchors[k0]) + \
                (1.0 if ov0 else -1.0) * float(sphere_offset)
            r_p1 = np.linalg.norm(anchors[k1]) + \
                (1.0 if ov1 else -1.0) * float(sphere_offset)
            idxs = [t % n for t in range(t0, t1 + 1)]
            seg = pts[idxs]
            ds = np.linalg.norm(np.diff(seg, axis=0), axis=1)
            s = np.concatenate([[0.0], np.cumsum(ds)])
            S = max(s[-1], 1.0e-9)
            for t_, m_ in enumerate(idxs):
                d0, d1 = s[t_], S - s[t_]
                if S < 2.0 * blend:
                    w = 0.5 * (1 + math.cos(math.pi * min(1.0, d0 / S)))
                    r = w * r_p0 + (1 - w) * r_p1
                elif d0 < blend:
                    w = 0.5 * (1 + math.cos(math.pi * d0 / blend))
                    r = w * r_p0 + (1 - w) * r_layer
                elif d1 < blend:
                    w = 0.5 * (1 + math.cos(math.pi * d1 / blend))
                    r = w * r_p1 + (1 - w) * r_layer
                else:
                    r = r_layer
                if not assigned[m_]:
                    new_r[m_] = r
                    assigned[m_] = True
        for v in range(m):                       # exact passage radii last
            k0, p_idx, ov0 = vis[v]
            new_r[p_idx % n] = np.linalg.norm(anchors[k0]) + \
                (1.0 if ov0 else -1.0) * float(sphere_offset)
        rr = np.linalg.norm(pts, axis=1)
        scale = np.where(rr > 1e-9, new_r / np.maximum(rr, 1e-9), 1.0)
        comps[ci] = pts * scale[:, None]
    messages.append("[fix] xyz layered rebuild: %d arc conflict pair(s), "
                    "%d gadget-threading arc(s); diagram needs %d radial "
                    "layer(s) (spacing %.1f)."
                    % (sum(len(c) for c in conflict) // 2,
                       sum(1 for t in threader if t), L, lam))
    return comps, L


def _enforce_crossing_orders(comps, model, anchors, sphere_offset,
                             ball_radius, messages, halfwidth=6):
    """V5.4: verify and enforce the DT over/under RADIAL order at every
    true crossing of the raw spherical curve.

    Overlapping bump/dip profiles of nearby crossings can partially cancel,
    silently flipping a crossing and changing the link type.  For each
    crossing anchor, the two strand visits are recovered by nearest-anchor
    segmentation, aligned with the model's per-component position lists,
    and the radial order at their closest approach is checked against
    model['over_at']; violations (or margins below sphere_offset) are
    repaired by pushing the two visits radially apart around the crossing.
    Returns (comps, n_fixed)."""
    comps = [np.array(c, float) for c in comps]
    anchors = np.asarray(anchors, float)
    ncross = len(anchors)
    # ---- segment each component into nearest-anchor (Voronoi) visit runs,
    # in curve order; every point belongs to its nearest crossing, so the
    # number of runs equals the number of crossing visits exactly.
    visits_by_comp = []
    for pts in comps:
        d = np.linalg.norm(pts[:, None, :] - anchors[None, :, :], axis=2)
        lab = np.argmin(d, axis=1)
        runs = []
        i = 0
        n = len(pts)
        while i < n:
            j = i
            while j + 1 < n and lab[j + 1] == lab[i]:
                j += 1
            runs.append([int(lab[i]), i, j])
            i = j + 1
        # merge fragments (wobble across a Voronoi boundary) and wraparound
        merged = []
        for r in runs:
            if merged and merged[-1][0] == r[0]:
                merged[-1][2] = r[2]
            else:
                merged.append(list(r))
        if len(merged) > 1 and merged[0][0] == merged[-1][0]:
            merged[0] = [merged[0][0],
                         merged[-1][1] - len(pts), merged[0][2]]
            merged.pop()
        # gap-tolerant clustering: Voronoi boundaries oscillate inside
        # tight gadget clusters, splitting one true visit into several
        # interleaved runs; re-merge same-anchor runs whose index gap is
        # smaller than a gadget span.
        gap_pts = max(6, int(round(ball_radius / max(
            1.0e-9, float(np.median(np.linalg.norm(
                np.diff(pts, axis=0), axis=1)))))))
        clustered = []
        last_of = {}
        for kk, i0, j0 in merged:
            if kk in last_of and i0 - clustered[last_of[kk]][2] <= gap_pts:
                clustered[last_of[kk]][2] = j0
            else:
                clustered.append([kk, i0, j0])
                last_of[kk] = len(clustered) - 1
        # wraparound cluster merge
        if len(clustered) > 1 and clustered[0][0] == clustered[-1][0] and \
                (clustered[0][1] + n - clustered[-1][2]) <= gap_pts:
            clustered[0] = [clustered[0][0],
                            clustered[-1][1] - n, clustered[0][2]]
            clustered.pop()
        # attach each run's closest approach to its anchor
        out_runs = []
        for kk, i0, j0 in clustered:
            idxs = [t % n for t in range(i0, j0 + 1)]
            dmin = float(np.min(d[idxs, kk]))
            out_runs.append((kk, i0, j0, dmin))
        visits_by_comp.append(out_runs)
    # ---- expected visit sequences from the model
    exp_cross, exp_over = [], []
    for cp in model["comp_positions"]:
        exp_cross.append([model["pos_cross"][p] for p in cp])
        exp_over.append([bool(model["over_at"][p]) for p in cp])
    # ---- order-aware matching of clustered runs to expected visits
    per_cross = {}
    ok_align = True
    for ci, runs in enumerate(visits_by_comp):
        want = exp_cross[ci]
        matched = _match_visit_sequence(runs, want)
        if matched is None:
            # fused U-turn visits: split deficient ids at the distance
            # profile's two minima, then retry once
            from collections import Counter as _C
            have = _C(r[0] for r in runs)
            need = _C(want)
            pts = comps[ci]
            n = len(pts)
            new_runs = []
            for (kk, i0, j0, dmin) in runs:
                if have[kk] < need.get(kk, 0):
                    idxs = [t % n for t in range(i0, j0 + 1)]
                    prof = np.linalg.norm(
                        pts[idxs] - anchors[kk][None, :], axis=1)
                    if len(prof) >= 5:
                        a_ = int(np.argmin(prof))
                        prof2 = prof.copy()
                        prof2[max(0, a_ - 2):a_ + 3] = np.inf
                        b_ = int(np.argmin(prof2))
                        lo_, hi_ = sorted((a_, b_))
                        cut = lo_ + int(np.argmax(prof[lo_:hi_ + 1]))
                        new_runs.append((kk, i0, i0 + cut,
                                         float(prof[:cut + 1].min())))
                        new_runs.append((kk, i0 + cut + 1, j0,
                                         float(prof[cut + 1:].min())))
                        continue
                new_runs.append((kk, i0, j0, dmin))
            matched = _match_visit_sequence(new_runs, want)
        if matched is None:
            import os as _os
            if _os.environ.get("V54_DEBUG"):
                messages.append("[dbg] comp %d want=%r" % (ci, want))
                messages.append("[dbg] comp %d runs=%r"
                                % (ci, [r[0] for r in runs]))
            ok_align = False
            break
        for t, k in enumerate(want):
            per_cross.setdefault(k, []).append(
                (ci, tuple(matched[t][:3]), exp_over[ci][t]))
    if not ok_align:
        messages.append("[WARN] xyz crossing-order check: could not align "
                        "curve visits with the diagram; skipped.")
        return comps, -1
    n_fixed = 0
    for _sweep in range(4):
        fixed_this = 0
        for k, pair in per_cross.items():
            if len(pair) != 2:
                continue
            (c1, (a1, i1, j1), ov1), (c2, (a2, i2, j2), ov2) = pair
            if ov1 == ov2:
                continue
            n1, n2 = len(comps[c1]), len(comps[c2])
            idx1 = [t % n1 for t in range(i1, j1 + 1)]
            idx2 = [t % n2 for t in range(i2, j2 + 1)]
            # the crossing passage = each run's closest point to the anchor
            d1 = np.linalg.norm(comps[c1][idx1] - anchors[k][None, :], axis=1)
            d2 = np.linalg.norm(comps[c2][idx2] - anchors[k][None, :], axis=1)
            a = int(np.argmin(d1)); b = int(np.argmin(d2))
            p1, p2 = comps[c1][idx1[a]], comps[c2][idx2[b]]
            r1, r2 = np.linalg.norm(p1), np.linalg.norm(p2)
            margin = (r1 - r2) if ov1 else (r2 - r1)
            if margin >= 0.9 * float(sphere_offset):
                continue
            fixed_this += 1
            deficit = 0.9 * float(sphere_offset) - margin
            for (cc, center, sgn_over) in ((c1, idx1[a], ov1),
                                           (c2, idx2[b], ov2)):
                sgn = 1.0 if sgn_over else -1.0
                nn = len(comps[cc])
                for t in range(-halfwidth, halfwidth + 1):
                    w = math.exp(-(t * t) / (0.5 * halfwidth ** 2))
                    m_ = (center + t) % nn
                    pt = comps[cc][m_]
                    rr = np.linalg.norm(pt)
                    if rr > 1.0e-9:
                        comps[cc][m_] = pt * (
                            1.0 + sgn * 0.55 * deficit * w / rr)
        n_fixed += fixed_this
        if fixed_this == 0:
            break
    if n_fixed:
        messages.append("[fix] xyz crossing-order repair: re-enforced the "
                        "DT over/under radial order at %d crossing(s)."
                        % n_fixed)
    # ---- evict gadget threaders: a strand that merely passes near a
    # crossing it does not belong to must NOT run radially BETWEEN that
    # crossing's over and under strands (that silently threads the
    # crossing and changes the link).  Push such runs fully above or
    # below the gadget, whichever side they already favor.
    member_pts = {}
    for k, pair in per_cross.items():
        s = set()
        for (cc, (aa, ii, jj), ov) in pair:
            nn = len(comps[cc])
            for t in range(ii - 2, jj + 3):
                s.add((cc, t % nn))
        member_pts[k] = s
    n_evict = 0
    core = 0.65 * float(ball_radius)
    for k in range(ncross):
        if k not in per_cross or len(per_cross[k]) != 2:
            continue
        Ra = float(np.linalg.norm(anchors[k]))
        lo_r = Ra - float(sphere_offset)
        hi_r = Ra + float(sphere_offset)
        pad = 0.6 * float(sphere_offset)
        for cc in range(len(comps)):
            pts = comps[cc]
            nn = len(pts)
            d = np.linalg.norm(pts - anchors[k][None, :], axis=1)
            inside = d < core
            if not inside.any():
                continue
            # contiguous runs of inside points, skipping member points
            i = 0
            while i < nn:
                if not inside[i] or (cc, i) in member_pts[k]:
                    i += 1
                    continue
                j = i
                while j + 1 < nn and inside[j + 1] and \
                        (cc, j + 1) not in member_pts[k]:
                    j += 1
                idxs = list(range(i, j + 1))
                rads = np.linalg.norm(pts[idxs], axis=1)
                if np.all(rads > hi_r + 0.2 * pad) or \
                        np.all(rads < lo_r - 0.2 * pad):
                    i = j + 1
                    continue
                n_evict += 1
                go_out = float(np.mean(rads)) >= Ra
                target = (hi_r + pad) if go_out else (lo_r - pad)
                w_edge = 4
                for t_, m_ in enumerate(idxs):
                    # smooth blend at the run edges
                    wl = min(1.0, (t_ + 1) / float(w_edge))
                    wr = min(1.0, (len(idxs) - t_) / float(w_edge))
                    w = min(wl, wr)
                    rr = np.linalg.norm(pts[m_])
                    if rr > 1.0e-9:
                        want_r = (1 - w) * rr + w * target
                        if go_out:
                            want_r = max(rr, want_r)
                        else:
                            want_r = min(rr, want_r)
                        comps[cc][m_] = pts[m_] * (want_r / rr)
                i = j + 1
    if n_evict:
        messages.append("[fix] xyz gadget-threading repair: pushed %d "
                        "passing strand run(s) clear of foreign "
                        "crossings." % n_evict)
    return comps, n_fixed + n_evict


def _run_xyz_repair_stage(xyz_components, clearance, stage, messages,
                          anchors=None, exclude_radius=0.0):
    comps, rep = _repair_xyz_clearance(xyz_components, clearance,
                                       anchors=anchors,
                                       exclude_radius=exclude_radius)
    if rep["initial"]:
        messages.append(
            "[fix] xyz %s: separated %d close approaches (< %.2f) in %d "
            "passes; %s; est. %d radial layer(s), max push %.2f"
            % (stage, rep["initial"], clearance, rep["iters"],
               "all resolved" if rep["resolved"] else
               "WARNING: NOT fully resolved", rep["layers"],
               rep["max_push"]))
    return comps, rep


def audit_xyz_components(xyz_components, dt_string):
    """Deep topology audit via audit_xyz.py (needs spherogram + snappy)."""
    try:
        import audit_xyz as _ax
    except Exception as e:
        return {"status": "unavailable", "detail": str(e)[:80]}
    return _ax.audit_components_against_dt(xyz_components, dt_string)


def build_spherical_xyz_components(
    model,
    P,
    centers,
    G=None,
    sphere_radius=50.0,
    sphere_extent=1.55,
    sphere_offset=5.0,
    sphere_bump_frac=1.0,
    xyz_spacing=1.8,
    sphere_layout="spherical-kamada",
    sphere_crossing_angle=15.0,
    direct_connecting=True,
    xyz_final_smooth=True,
    xyz_smooth_window=10,
    xyz_smooth_passes=5,
    surface_shape="ellipsoid",
    surface_auto_orient=True,
    surface_auto_aspect=True,
    surface_aspect=1.6,
    surface_tube=0.35,
    surface_orient=0.0,
    surface_tilt=0.0,
    skeleton_out=None,
    skeleton_grid=24,
    xyz_repair=True,
    xyz_clearance=None,
    audit_dt=None,
    status_messages=None,
):
    """
    Return one Nx3 array per component for the spherical diagram.

    V5.4: when ``xyz_repair`` is on (default), a clearance-repair pass runs
    after the raw build, after smoothing, and after the surface warp,
    separating any two non-neighboring arcs closer than ``xyz_clearance``
    (default 0.8 * sphere_offset) into distinct radial layers.  When
    ``audit_dt`` is given (the DT string), a final topology audit checks
    that the finished 3D curve realizes that link (requires spherogram +
    snappy; see audit_xyz.py).  Messages are appended to
    ``status_messages`` (a list) if provided, else printed.

    ``skeleton_out``, when given a dict, is filled with the 'mapped skeleton':
    ``points`` (one 3D anchor per crossing, same indexing as the crossing IDs),
    ``segments`` (one polyline per strand step between crossings), warped onto
    the same surface as the curve.  Kamada sphere layouts only.

    Parameter units:
      sphere_radius and crossing_offset are absolute XYZ coordinate distances.
      xyz_spacing is an absolute target adjacent-point separation in XYZ units.
      sphere_crossing_angle is in degrees.
      sphere_bump_frac is dimensionless, with 0 < fraction <= 1.
      xyz_final_smooth is a final closed-curve smoothing pass on the written 3D
      polyline. xyz_smooth_window is a number of neighboring points, and
      xyz_smooth_passes is a small integer count of repeated smoothing passes.
    """
    sphere_radius = float(sphere_radius)
    sphere_extent = float(sphere_extent)
    sphere_offset = float(sphere_offset)
    sphere_bump_frac = float(sphere_bump_frac)
    xyz_spacing = float(xyz_spacing)
    sphere_crossing_angle_degrees = float(sphere_crossing_angle)
    sphere_crossing_angle = math.radians(sphere_crossing_angle_degrees)
    xyz_smooth_window = int(round(float(xyz_smooth_window)))
    xyz_smooth_passes = int(round(float(xyz_smooth_passes)))

    if sphere_radius <= 0.0:
        raise ValueError("--sphere-radius must be positive.")
    if sphere_extent <= 0.0:
        raise ValueError("--sphere-extent must be positive.")
    if sphere_offset < 0.0:
        raise ValueError("--crossing-offset must be non-negative.")
    if sphere_offset >= sphere_radius:
        raise ValueError("--crossing-offset must be smaller than --sphere-radius.")
    if sphere_bump_frac <= 0.0 or sphere_bump_frac > 1.0:
        raise ValueError("--sphere-bump-frac must be in the range 0 < value <= 1.")
    if xyz_spacing <= 0.0:
        raise ValueError("--xyz-spacing must be positive.")
    if sphere_crossing_angle_degrees <= 0.0:
        raise ValueError("--sphere-crossing-angle must be positive degrees.")
    if sphere_crossing_angle_degrees >= 28.5:
        raise ValueError("--sphere-crossing-angle should be less than 28.5 degrees.")
    if xyz_smooth_window < 0:
        raise ValueError("--xyz-smooth-window must be non-negative.")
    if xyz_smooth_passes < 0:
        raise ValueError("--xyz-smooth-passes must be non-negative.")

    layout = str(sphere_layout or "spherical-kamada").strip().lower()
    if skeleton_out is None and layout in ("spherical-kamada", "shaped-kamada"):
        skeleton_out = {}          # V5.4: anchors needed for clearance repair
    _stereo_safe = (layout == "stereo-safe")
    if _stereo_safe:
        # V5.4 correct-by-construction mode: lift the AUDITED planar layout
        # stereographically with radial crossing bumps and do not move
        # angular positions afterwards (strands disjoint in the 2D layout
        # can then never intersect in 3D).  Smoothing and clearance repair
        # are disabled because both move positions and would void the
        # guarantee; the topology audit still runs as a check.
        layout = "stereographic"
        xyz_final_smooth = False
        xyz_repair = False
    if layout == "stereographic":
        xyz_components = _build_stereo_xyz_components(
            model,
            P,
            centers,
            sphere_radius=sphere_radius,
            sphere_extent=sphere_extent,
            sphere_offset=sphere_offset,
            sphere_bump_frac=sphere_bump_frac,
            xyz_spacing=(0.6 * xyz_spacing if _stereo_safe else xyz_spacing),
        )
    elif layout in ("spherical-kamada", "shaped-kamada"):
        # Both use the sphere-native construction; shaped-kamada additionally
        # warps the finished spherical curve onto the chosen surface below.
        xyz_components = _build_spherical_kamada_xyz_components(
            model,
            G,
            sphere_radius=sphere_radius,
            sphere_offset=sphere_offset,
            sphere_bump_frac=sphere_bump_frac,
            xyz_spacing=xyz_spacing,
            crossing_angle=sphere_crossing_angle,
            direct_connecting=direct_connecting,
            skeleton_out=skeleton_out,
            skeleton_grid=skeleton_grid,
        )
    else:
        raise ValueError("Unknown --sphere-layout %r." % sphere_layout)

    _msgs = status_messages if status_messages is not None else []
    _clr = float(xyz_clearance) if xyz_clearance else 0.8 * float(sphere_offset)
    _excl = max(2.4 * float(sphere_offset),
                float(sphere_radius) * math.sin(sphere_crossing_angle))
    _anch = None
    if skeleton_out is not None and skeleton_out.get("points") is not None \
            and len(skeleton_out.get("points", [])):
        _anch = np.asarray(skeleton_out["points"], float) * (
            float(sphere_radius) if skeleton_out.get("unit") else 1.0)
    if xyz_repair and _clr > 0.0:
        xyz_components, _ = _run_xyz_repair_stage(
            xyz_components, _clr, "raw build", _msgs,
            anchors=_anch, exclude_radius=_excl)

    xyz_components = _apply_final_xyz_smoothing(
        xyz_components,
        sphere_radius=sphere_radius,
        xyz_spacing=xyz_spacing,
        enabled=bool(xyz_final_smooth),
        smooth_window=xyz_smooth_window,
        smooth_passes=xyz_smooth_passes,
    )

    if xyz_repair and _clr > 0.0:
        xyz_components, _ = _run_xyz_repair_stage(
            xyz_components, _clr, "post-smoothing", _msgs,
            anchors=_anch, exclude_radius=_excl)

    _warp_pts = None
    if layout == "shaped-kamada":
        xyz_components, _warp_pts = _warp_components_to_surface(
            xyz_components,
            sphere_radius=sphere_radius,
            surface_shape=surface_shape,
            surface_aspect=surface_aspect,
            surface_tube=surface_tube,
            xyz_spacing=xyz_spacing,
            surface_auto_orient=surface_auto_orient,
            surface_auto_aspect=surface_auto_aspect,
            surface_orient=surface_orient,
            surface_tilt=surface_tilt,
            return_warp=True,
        )

    if xyz_repair and _clr > 0.0 and layout == "shaped-kamada":
        _anch_w = _anch
        if _anch is not None and _warp_pts is not None:
            _anch_w = _warp_pts(np.asarray(_anch, float))
        xyz_components, _ = _run_xyz_repair_stage(
            xyz_components, _clr, "post-warp", _msgs,
            anchors=_anch_w, exclude_radius=_excl)

    if audit_dt:
        _res = audit_xyz_components(xyz_components, audit_dt)
        if _res["status"] == "ok":
            _msgs.append("[ok] xyz topology audit: %s" % _res["detail"])
        elif _res["status"] == "unavailable":
            _msgs.append("[info] xyz topology audit skipped: %s"
                         % _res["detail"])
        else:
            _msgs.append("[%s] xyz topology audit: %s"
                         % ("WARN" if _res["status"] == "inconclusive"
                            else "FAIL", _res["detail"]))
            if _res["status"] == "FAIL":
                if _stereo_safe:
                    _msgs.append(
                        "[hint] stereo-safe failed unexpectedly; do not use "
                        "this XYZ output. Recheck the DT code and audit "
                        "dependencies before investigating the curve."
                    )
                else:
                    _msgs.append("[hint] dense diagrams overwhelm the kamada "
                                 "sphere layouts; re-run with --sphere-layout "
                                 "stereo-safe (correct by construction, "
                                 "audited).")

    if status_messages is None:
        for _m in _msgs:
            print(_m)

    # V5.0: scale the unit skeleton to the sphere and put it through the same
    # surface warp as the curve, so the overlay matches the drawn strands.
    if skeleton_out is not None and skeleton_out.get("unit"):
        _sk_pts = np.asarray(skeleton_out.get("points", np.zeros((0, 3))), float) * sphere_radius
        _sk_segs = [np.asarray(s, float) * sphere_radius
                    for s in skeleton_out.get("segments", [])]
        _sk_wires = [np.asarray(w, float) * sphere_radius
                     for w in skeleton_out.get("surface_wires", [])]
        if _warp_pts is not None:
            _sk_pts = _warp_pts(_sk_pts)
            _sk_segs = [_warp_pts(s) for s in _sk_segs]
            _sk_wires = [_warp_pts(w) for w in _sk_wires]
        skeleton_out["points"] = _sk_pts
        skeleton_out["segments"] = _sk_segs
        skeleton_out["surface_wires"] = _sk_wires
        skeleton_out["unit"] = False

    return xyz_components


def write_spherical_xyz(
    model,
    P,
    centers,
    path,
    G=None,
    sphere_radius=50.0,
    sphere_extent=1.55,
    sphere_offset=5.0,
    sphere_bump_frac=1.0,
    xyz_spacing=1.8,
    close_components=False,
    decimals=9,
    sphere_layout="spherical-kamada",
    sphere_crossing_angle=15.0,
    direct_connecting=True,
    xyz_final_smooth=True,
    xyz_smooth_window=10,
    xyz_smooth_passes=5,
    surface_shape="ellipsoid",
    surface_auto_orient=True,
    surface_auto_aspect=True,
    surface_aspect=1.6,
    surface_tube=0.35,
    surface_orient=0.0,
    surface_tilt=0.0,
    xyz_repair=True,
    xyz_clearance=None,
    audit_dt=None,
    status_stream=None,
):
    """
    Write spherical component polylines as plain x y z coordinates.

    This is not chemical XYZ format with atom counts. Every line has three
    floating-point values, and a blank line separates components. By default the
    first point is NOT repeated at the end of each component block; use
    --xyz-close-components if your viewer wants explicit closure.
    """
    xyz_components = build_spherical_xyz_components(
        model,
        P,
        centers,
        G=G,
        sphere_radius=sphere_radius,
        sphere_extent=sphere_extent,
        sphere_offset=sphere_offset,
        sphere_bump_frac=sphere_bump_frac,
        xyz_spacing=xyz_spacing,
        sphere_layout=sphere_layout,
        sphere_crossing_angle=sphere_crossing_angle,
        direct_connecting=direct_connecting,
        xyz_final_smooth=xyz_final_smooth,
        xyz_smooth_window=xyz_smooth_window,
        xyz_smooth_passes=xyz_smooth_passes,
        surface_shape=surface_shape,
        surface_auto_orient=surface_auto_orient,
        surface_auto_aspect=surface_auto_aspect,
        surface_aspect=surface_aspect,
        surface_tube=surface_tube,
        surface_orient=surface_orient,
        surface_tilt=surface_tilt,
        xyz_repair=xyz_repair,
        xyz_clearance=xyz_clearance,
        audit_dt=audit_dt,
        status_messages=(_v54_msgs := []),
    )
    for _m in _v54_msgs:
        if status_stream is not None:
            status_stream.write(_m + "\n")
        else:
            print(_m)
    ensure_parent_dir(path)
    nd = max(3, int(decimals))
    fmt = "%%.%df %%.%df %%.%df\n" % (nd, nd, nd)
    n_written = 0
    with open(path, "w") as fh:
        for ci, arr in enumerate(xyz_components):
            pts = arr
            if close_components and len(pts) > 0:
                pts = np.vstack([pts, pts[0]])
            for x, y, z in pts:
                fh.write(fmt % (x, y, z))
                n_written += 1
            if ci != len(xyz_components) - 1:
                fh.write("\n")
    return n_written




# --------------------------------------------------------------------------- #
#  7b. Illustrative 2D projections of the 3D curve (V5.0)
# --------------------------------------------------------------------------- #
PROJECTION_VIEW_PRESETS = {
    # name: (elev, azim, roll) in degrees; matplotlib view_init conventions.
    "top": (90.0, -90.0, 0.0),
    "front": (0.0, -90.0, 0.0),
    "side": (0.0, 0.0, 0.0),
}
PROJECTION_VIEW_CHOICES = ("current", "top", "front", "side")


def _parse_projection_views(text):
    """Parse a comma-separated view list; unknown names raise ValueError."""
    views = []
    for tok in str(text or "").replace(";", ",").split(","):
        name = tok.strip().lower()
        if not name:
            continue
        if name not in PROJECTION_VIEW_CHOICES:
            raise ValueError(
                "Unknown projection view %r (choices: %s)."
                % (name, ", ".join(PROJECTION_VIEW_CHOICES))
            )
        if name not in views:
            views.append(name)
    return views or ["current"]


def _projection_basis(elev_deg, azim_deg, roll_deg=0.0):
    """
    Orthonormal (right, up, toward-viewer) basis for an orthographic projection
    with matplotlib-like elev/azim angles plus a roll about the view direction.
    """
    el = math.radians(float(elev_deg))
    az = math.radians(float(azim_deg))
    f = np.array([
        math.cos(el) * math.cos(az),
        math.cos(el) * math.sin(az),
        math.sin(el),
    ], float)                              # unit vector pointing at the viewer
    zaxis = np.array([0.0, 0.0, 1.0])
    r = np.cross(zaxis, f)
    nr = float(np.linalg.norm(r))
    if nr < 1.0e-9:                        # looking along +-z: pick x as right
        r = np.array([1.0, 0.0, 0.0])
    else:
        r = r / nr
    u = np.cross(f, r)
    ro = math.radians(float(roll_deg or 0.0))
    if abs(ro) > 1.0e-12:
        cr, sr = math.cos(ro), math.sin(ro)
        r, u = cr * r + sr * u, -sr * r + cr * u
    return r, u, f


def project_points_3d(arr, basis):
    """Project Nx3 points; returns (Nx2 screen coords, N depths toward viewer)."""
    r, u, f = basis
    arr = np.asarray(arr, float)
    if arr.size == 0:
        return np.zeros((0, 2)), np.zeros(0)
    xy = np.column_stack([arr.dot(r), arr.dot(u)])
    return xy, arr.dot(f)


def auto_projection_view(xyz_components, coarse_step=10.0, refine_step=2.0):
    """
    V5.1: search the view sphere for the direction that maximizes the WORST
    component's projected roundness.

    For a candidate view, each component's closed curve is projected and scored
    by its isoperimetric ratio 4*pi*area/perimeter^2 (1 for a circle, ~0 for a
    sliver).  The score of the view is the MINIMUM over components, so the
    optimum is the view in which no ring degenerates to an edge-on line --
    the failure mode of links whose rings lie in near-orthogonal planes.
    Returns (elev, azim).
    """
    curves = []
    for arr in xyz_components:
        arr = np.asarray(arr, float)
        if len(arr) < 8:
            continue
        k = max(1, len(arr) // 120)
        curves.append(arr[::k])
    if not curves:
        return 25.0, -60.0

    def _score(elev, azim):
        basis = _projection_basis(elev, azim, 0.0)
        worst = 1.0e9
        for cur in curves:
            xy, _d = project_points_3d(cur, basis)
            x = xy[:, 0]; y = xy[:, 1]
            x2 = np.roll(x, -1); y2 = np.roll(y, -1)
            area = abs(0.5 * float(np.sum(x * y2 - x2 * y)))
            per = float(np.sum(np.hypot(x2 - x, y2 - y)))
            if per <= 1.0e-9:
                return 0.0
            worst = min(worst, 4.0 * math.pi * area / (per * per))
        return worst

    best = (25.0, -60.0)
    best_s = -1.0
    e = -80.0
    while e <= 80.0 + 1.0e-9:
        a = -180.0
        while a < 180.0 - 1.0e-9:
            s = _score(e, a)
            if s > best_s:
                best_s, best = s, (e, a)
            a += float(coarse_step)
        e += float(coarse_step)
    e0, a0 = best
    e = e0 - float(coarse_step)
    while e <= e0 + float(coarse_step) + 1.0e-9:
        a = a0 - float(coarse_step)
        while a <= a0 + float(coarse_step) + 1.0e-9:
            s = _score(e, a)
            if s > best_s:
                best_s, best = s, (e, a)
            a += float(refine_step)
        e += float(refine_step)
    return best


def render_projection_on_axis(ax, xyz_components, elev=25.0, azim=-60.0,
                              roll=0.0, line_width=2.2, skeleton=None,
                              show_skeleton=False, skeleton_ids=False,
                              skeleton_fontsize=6.0, chunk=8, title=None,
                              depth_fade=0.55, perspective=0.0,
                              skeleton_only=False, fast=False,
                              grid_color="0.82", grid_lw=0.4,
                              show_chords=True, show_crossing_dots=False,
                              grid_basis=None):
    """
    Draw an ILLUSTRATIVE 2D orthographic projection of the 3D curve.

    The curve of each component is cut into short chunks that are depth-sorted
    (matplotlib zorder) and stroked over a slightly wider white halo, so at
    every apparent crossing of the projection the nearer strand visually passes
    over the farther one with a clean under-strand gap -- the same reading as a
    knot diagram, but produced from the 3D geometry.  This is the V5.0 route
    for links whose flat 2D diagrams are unreadable.
    """
    import matplotlib.patheffects as _pe

    basis = _projection_basis(elev, azim, roll)
    grid_basis = basis if grid_basis is None else grid_basis
    ax.clear()
    ax.axis("off")

    # Depth normalization across everything drawn (strands AND skeleton).
    all_dep = []
    proj = []
    for arr in xyz_components:
        arr = np.asarray(arr, float)
        if len(arr) == 0:
            proj.append((np.zeros((0, 2)), np.zeros(0)))
            continue
        closed = np.vstack([arr, arr[:1]])
        xy, dep = project_points_3d(closed, basis)
        proj.append((xy, dep))
        all_dep.append(dep)
    if not all_dep:
        return
    dmin = min(float(d.min()) for d in all_dep)
    dmax = max(float(d.max()) for d in all_dep)
    dspan = max(dmax - dmin, 1.0e-9)

    # Optional perspective: a camera on the view axis at 'perspective' scene
    # radii; nearer points are magnified.  Edge-on rings open into a lens.
    persp = float(perspective or 0.0)
    if persp > 0.0:
        persp = max(persp, 1.2)
        allxy = np.vstack([xy for xy, d in proj if len(xy)])
        c2 = allxy.mean(axis=0)
        r_scene = max(float(np.linalg.norm(allxy - c2, axis=1).max()),
                      0.5 * dspan, 1.0e-9)
        dcen = 0.5 * (dmin + dmax)
        D = persp * r_scene

        def _persp_apply(xy, dep):
            fac = D / np.maximum(D - (dep - dcen), 1.0e-3 * r_scene)
            return c2 + (xy - c2) * fac[:, None]
    else:
        def _persp_apply(xy, dep):
            return xy

    fade = min(1.0, max(0.0, float(depth_fade or 0.0)))

    def _alpha_of(depth_mean):
        if fade <= 0.0:
            return 1.0
        t = (float(depth_mean) - dmin) / dspan
        return max(0.08, (1.0 - fade) + fade * (0.15 + 0.85 * t))

    # Butt-capped halo: a round/projecting cap would overhang the chunk end and
    # erase the strand's own continuation (drawn at a lower depth), producing a
    # dashed artifact.  With a butt cap the halo ends exactly at the joint,
    # which the overlapping neighbour chunk repaints.
    halo = [_pe.Stroke(linewidth=float(line_width) * 2.8, foreground="white",
                       capstyle="butt"),
            _pe.Normal()]
    if not skeleton_only:
        for ci, (xy, dep) in enumerate(proj):
            npts = len(xy)
            if npts < 2:
                continue
            color = _default_color_of(ci)
            xyp = _persp_apply(xy, dep)
            if fast:
                # Fast mode (live dragging): one plain line per component.
                ax.plot(xyp[:, 0], xyp[:, 1], color=color,
                        lw=max(1.0, 0.6 * float(line_width)), alpha=0.9)
                continue
            step = max(2, int(chunk))
            # Chunks overlap by a full step on each side: the colour of a chunk
            # then always covers the white-halo bleed of its neighbours at the
            # joints, so the strand reads as continuous (gaps appear ONLY where
            # a genuinely nearer strand crosses).  Closed curve -> wrap.
            xy_ext = np.vstack([xyp, xyp[1:min(step + 1, npts)]])
            next_len = len(xy_ext)
            for i in range(0, npts - 1, step):
                core = slice(i, min(i + step + 1, npts))
                z = 10.0 + 5.0 * ((float(dep[core].mean()) - dmin) / dspan)
                lo = max(0, i - step)
                hi = min(i + 2 * step + 1, next_len)
                ax.plot(xy_ext[lo:hi, 0], xy_ext[lo:hi, 1], color=color,
                        lw=float(line_width), solid_capstyle="round",
                        zorder=z, alpha=_alpha_of(dep[core].mean()),
                        path_effects=halo)

    if (show_skeleton or skeleton_only) and skeleton and len(skeleton.get("points", [])):
        # V5.2: score_diagram-style mapped framework -- a light depth-shaded
        # wireframe of the WHOLE mapped surface, NEUTRAL dark-gray strand
        # chords (component colours stay reserved for the smooth strands), and
        # black depth-faded crossing anchors.
        pts = np.asarray(skeleton["points"], float)
        pxy, pdep = project_points_3d(pts, basis)
        pxy = _persp_apply(pxy, pdep)
        # -- surface wireframe (drawn beneath everything, subtle depth shade).
        for wirearr in skeleton.get("surface_wires", []):
            warr = np.asarray(wirearr, float)
            if len(warr) < 2:
                continue
            wxy, wdep = project_points_3d(warr, grid_basis)
            wxy = _persp_apply(wxy, wdep)
            half = max(4, len(wxy) // 10)
            for j0 in range(0, len(wxy) - 1, half):
                j1 = min(j0 + half + 1, len(wxy))
                zt = ((float(wdep[j0:j1].mean()) - dmin) / dspan)
                ax.plot(wxy[j0:j1, 0], wxy[j0:j1, 1], color=grid_color,
                        lw=float(grid_lw), zorder=1.0 + zt,
                        alpha=0.25 + 0.45 * max(0.0, min(1.0, zt)))
        ss = np.linspace(0.0, 1.0, 13)
        w0 = ((1.0 - ss) ** 2)[:, None]
        w1 = (2.0 * (1.0 - ss) * ss)[:, None]
        w2 = (ss ** 2)[:, None]
        base_z = 30.0 if skeleton_only else 20.0
        lw_fw = 1.6 if skeleton_only else 1.0
        for si, segarr in enumerate(
                skeleton.get("segments", []) if show_chords else []):
            seg = np.asarray(segarr, float)
            if len(seg) == 3:
                arc = w0 * seg[0] + w1 * seg[1] + w2 * seg[2]
            else:
                arc = seg
            sxy, sdep = project_points_3d(arc, basis)
            sxy = _persp_apply(sxy, sdep)
            for j in range(len(sxy) - 1):
                dmid = 0.5 * (sdep[j] + sdep[j + 1])
                zt = (float(dmid) - dmin) / dspan
                ax.plot(sxy[j:j + 2, 0], sxy[j:j + 2, 1], color="0.30",
                        lw=lw_fw, solid_capstyle="round",
                        zorder=base_z + 5.0 * zt,
                        alpha=_alpha_of(dmid) if fade > 0.0
                        else (0.12 + 0.83 * zt))
        if show_crossing_dots:
            pa = np.array([(0.2 + 0.8 * ((d - dmin) / dspan)) for d in pdep])
            ax.scatter(pxy[:, 0], pxy[:, 1], s=11.0, c="#222222",
                       alpha=None, zorder=base_z + 5.5,
                       edgecolors="none")
            # Per-point alpha via colours (matplotlib scatter alpha is scalar).
            try:
                from matplotlib.colors import to_rgba
                cols = [to_rgba("#222222", float(a)) for a in pa]
                ax.collections[-1].set_facecolor(cols)
            except Exception:
                pass
        if skeleton_ids:
            for k in range(len(pxy)):
                ax.annotate(
                    str(k + 1), (pxy[k, 0], pxy[k, 1]),
                    textcoords="offset points", xytext=(3, 3),
                    fontsize=float(skeleton_fontsize), zorder=base_z + 6.0,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white",
                              ec="0.6", lw=0.4, alpha=0.85),
                )

    ax.set_aspect("equal")
    ax.relim()
    ax.autoscale_view()
    ax.margins(0.04)
    if title:
        ax.set_title(title, fontsize=9)


def save_projection_images(xyz_components, base_path, views, elev, azim,
                           roll=0.0, line_width=2.2, skeleton=None,
                           show_skeleton=False, skeleton_ids=False,
                           figsize=8.0, formats=("svg", "png"), dpi=200,
                           depth_fade=0.55, perspective=0.0,
                           skeleton_only=False, grid_color="0.82",
                           grid_lw=0.4, show_chords=True,
                           show_crossing_dots=False, fixed_grid=False,
                           fixed_grid_view=None):
    """
    Write the requested projection views next to ``base_path`` (typically the
    .xyz output).  File pattern: <base>_proj_<view>.<ext>.  Returns the list of
    written paths.
    """
    import matplotlib.pyplot as _plt

    stem = os.path.splitext(str(base_path))[0]
    written = []
    if fixed_grid_view is None:
        fixed_grid_view = (float(elev), float(azim), float(roll))
    grid_basis = (_projection_basis(*fixed_grid_view)
                  if bool(fixed_grid) else None)
    for name in views:
        if name == "current":
            e, a, r = float(elev), float(azim), float(roll)
        else:
            e, a, r = PROJECTION_VIEW_PRESETS[name]
        fig, ax = _plt.subplots(figsize=(float(figsize), float(figsize)))
        try:
            render_projection_on_axis(
                ax, xyz_components, elev=e, azim=a, roll=r,
                line_width=line_width, skeleton=skeleton,
                show_skeleton=show_skeleton, skeleton_ids=skeleton_ids,
                depth_fade=depth_fade, perspective=perspective,
                skeleton_only=skeleton_only,
                grid_color=grid_color, grid_lw=grid_lw, show_chords=show_chords,
                show_crossing_dots=show_crossing_dots,
                grid_basis=grid_basis,
                title="projection '%s' (elev=%g, azim=%g, roll=%g)" % (name, e, a, r),
            )
            for ext in formats:
                path = "%s_proj_%s.%s" % (stem, name, ext)
                ensure_parent_dir(path)
                fig.savefig(path, dpi=dpi, bbox_inches="tight")
                written.append(path)
        finally:
            _plt.close(fig)
    return written


# --------------------------------------------------------------------------- #
#  8. Lightweight Tk 3D XYZ viewer
# --------------------------------------------------------------------------- #
def component_label(index):
    """Return spreadsheet-style component labels: A, B, ..., Z, AA, AB, ..."""
    if index < 0:
        raise ValueError("Component index must be non-negative.")
    label = ""
    n = int(index)
    while True:
        label = chr(ord("A") + (n % 26)) + label
        n = n // 26 - 1
        if n < 0:
            break
    return label


def set_equal_aspect_3d(ax, points):
    """Set a 3D axis to equal aspect ratio for the supplied points."""
    pts = np.asarray(points, float)
    if pts.ndim != 2 or pts.shape[1] != 3 or len(pts) == 0:
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-1.0, 1.0)
        ax.set_zlim(-1.0, 1.0)
        return
    x = pts[:, 0]
    y = pts[:, 1]
    z = pts[:, 2]
    max_range = max(float(x.max() - x.min()), float(y.max() - y.min()), float(z.max() - z.min()))
    if max_range <= 1.0e-12:
        max_range = 1.0
    mid_x = 0.5 * float(x.max() + x.min())
    mid_y = 0.5 * float(y.max() + y.min())
    mid_z = 0.5 * float(z.max() + z.min())
    half = 0.5 * max_range
    ax.set_xlim(mid_x - half, mid_x + half)
    ax.set_ylim(mid_y - half, mid_y + half)
    ax.set_zlim(mid_z - half, mid_z + half)


def _xyz_audit_display(messages):
    """Return (text, level) for the topology-audit lines in build messages."""
    lines = [str(msg) for msg in (messages or ())
             if "xyz topology audit" in str(msg) or str(msg).startswith("[hint]")]
    text = "\n".join(lines)
    if any(line.startswith("[FAIL]") for line in lines):
        return text, "fail"
    if any(line.startswith("[WARN]") for line in lines):
        return text, "warn"
    if any(line.startswith("[ok]") for line in lines):
        return text, "ok"
    if lines:
        return text, "info"
    return "[info] xyz topology audit was not run.", "info"


def _audit_foreground(level):
    return {
        "ok": "#176b2c",
        "fail": "#a40000",
        "warn": "#8a5200",
        "info": "#4b5563",
    }.get(level, "#4b5563")


def _open_xyz_viewer_window(parent, xyz_components, closed=True,
                            title="Sphere XYZ preview", audit_messages=None):
    """Open a Tk popup with an interactive Matplotlib 3D view of components.

    This follows the same basic idea as view_xyzV3.py in the curve_it project:
    blank-line-separated components are shown with distinct colors, start points
    are marked, and component visibility can be toggled.
    """
    import tkinter as tk
    from tkinter import ttk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    clean_components = []
    for arr in xyz_components:
        pts = np.asarray(arr, float)
        if pts.ndim == 2 and pts.shape[1] == 3 and len(pts) > 0:
            clean_components.append(pts)
    if not clean_components:
        raise ValueError("No valid XYZ components are available to view.")

    # Ensure the 3D projection is registered even in older Matplotlib builds.
    try:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except Exception:
        pass

    top = tk.Toplevel(parent)
    top.title(title)
    top.geometry("980x720")
    top.minsize(760, 560)
    top.columnconfigure(0, weight=1)
    top.columnconfigure(1, weight=0)
    top.rowconfigure(0, weight=1)
    top.rowconfigure(1, weight=0)
    top.rowconfigure(2, weight=0)

    fig = Figure(figsize=(7.6, 6.2), dpi=100)
    ax3 = fig.add_subplot(111, projection="3d")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.04)

    color_map = plt.get_cmap("tab10")
    artists_by_component = []
    for i, pts in enumerate(clean_components):
        color = color_map(i % 10)
        label = "%s (%d pts)" % (component_label(i), len(pts))
        artists = []
        if len(pts) >= 2:
            line, = ax3.plot(
                pts[:, 0], pts[:, 1], pts[:, 2],
                color=color, linewidth=1.4, label=label,
            )
            artists.append(line)
            if closed and len(pts) > 2:
                close_line, = ax3.plot(
                    [pts[-1, 0], pts[0, 0]],
                    [pts[-1, 1], pts[0, 1]],
                    [pts[-1, 2], pts[0, 2]],
                    color=color, linestyle=":", linewidth=1.0,
                )
                artists.append(close_line)
        if len(pts) <= 6000:
            scatter = ax3.scatter(
                pts[:, 0], pts[:, 1], pts[:, 2],
                color=[color], s=8, depthshade=True,
            )
            artists.append(scatter)
        start = ax3.scatter(
            pts[0, 0], pts[0, 1], pts[0, 2],
            color="red", s=34, edgecolor="k", depthshade=True,
        )
        artists.append(start)
        artists_by_component.append(artists)

    all_points = np.vstack(clean_components)
    set_equal_aspect_3d(ax3, all_points)
    ax3.set_xlabel("X")
    ax3.set_ylabel("Y")
    ax3.set_zlabel("Z")
    ax3.set_title(title, fontfamily=DIAGRAM_FONT_FAMILY)
    ax3.legend(loc="best")

    plot_frame = ttk.Frame(top)
    plot_frame.grid(row=0, column=0, sticky="nsew")
    plot_frame.columnconfigure(0, weight=1)
    plot_frame.rowconfigure(0, weight=1)
    canvas3 = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas3.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    toolbar_frame = ttk.Frame(top)
    toolbar_frame.grid(row=1, column=0, sticky="ew")
    toolbar = NavigationToolbar2Tk(canvas3, toolbar_frame)
    toolbar.update()

    controls = ttk.Frame(top, padding=(8, 8, 8, 8))
    controls.grid(row=0, column=1, rowspan=2, sticky="ns")
    ttk.Label(controls, text="Show components", font=("TkDefaultFont", 10, "bold")).grid(
        row=0, column=0, sticky="w", pady=(0, 6)
    )
    visibility = []

    def update_visibility():
        for idx, var in enumerate(visibility):
            visible = bool(var.get())
            for artist in artists_by_component[idx]:
                artist.set_visible(visible)
        canvas3.draw_idle()

    for i in range(len(clean_components)):
        var = tk.BooleanVar(value=True)
        visibility.append(var)
        ttk.Checkbutton(
            controls,
            text="%s (%d pts)" % (component_label(i), len(clean_components[i])),
            variable=var,
            command=update_visibility,
        ).grid(row=i + 1, column=0, sticky="w", pady=1)

    def set_all(value):
        for var in visibility:
            var.set(bool(value))
        update_visibility()

    btn_row = len(clean_components) + 2
    ttk.Button(controls, text="All", command=lambda: set_all(True)).grid(
        row=btn_row, column=0, sticky="ew", pady=(8, 2)
    )
    ttk.Button(controls, text="None", command=lambda: set_all(False)).grid(
        row=btn_row + 1, column=0, sticky="ew", pady=2
    )
    ttk.Button(controls, text="Close", command=top.destroy).grid(
        row=btn_row + 2, column=0, sticky="ew", pady=(10, 2)
    )

    audit_text, audit_level = _xyz_audit_display(audit_messages)
    audit_label = tk.Label(
        top,
        text=audit_text,
        anchor="w",
        justify="left",
        wraplength=900,
        padx=8,
        pady=5,
        relief="groove",
        foreground=_audit_foreground(audit_level),
    )
    audit_label.grid(row=2, column=0, columnspan=2, sticky="ew", padx=4, pady=4)

    canvas3.draw_idle()
    return top


# --------------------------------------------------------------------------- #
#  8. CSV mapping table
# --------------------------------------------------------------------------- #
def write_table(model, label_coords, crossing_xy, path, crossing_ids=None):
    label_of = model["label_of"]
    comp_of = model["comp_of"]
    over_at = model["over_at"]
    if crossing_ids is None:
        crossing_ids = default_crossing_ids(model)

    ensure_parent_dir(path)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "crossing_index",
            "odd_label",
            "even_label_signed",
            "odd_component",
            "even_component",
            "over_label",
            "under_label",
            "crossing_x",
            "crossing_y",
            "odd_x",
            "odd_y",
            "even_x",
            "even_y",
            "odd_order_crossing_index",
        ])
        for k, c in enumerate(model["crossings"]):
            odd = c["odd"]
            even = c["even"]
            odd_label = label_of(odd)
            even_label = label_of(even)
            over_label = even_label if over_at[even] else odd_label
            under_label = odd_label if over_at[even] else even_label

            cx = cy = ox = oy = ex = ey = ""
            if k in crossing_xy:
                cx, cy = ("%.6f" % crossing_xy[k][0], "%.6f" % crossing_xy[k][1])
            if odd in label_coords:
                ox, oy = ("%.6f" % label_coords[odd][0], "%.6f" % label_coords[odd][1])
            if even in label_coords:
                ex, ey = ("%.6f" % label_coords[even][0], "%.6f" % label_coords[even][1])

            w.writerow([
                crossing_ids[k],
                odd_label,
                even_label,
                comp_of[odd] + 1,
                comp_of[even] + 1,
                over_label,
                under_label,
                cx,
                cy,
                ox,
                oy,
                ex,
                ey,
                "c%d" % (k + 1),
            ])


# --------------------------------------------------------------------------- #
#  9. Pipeline, CLI, and GUI
# --------------------------------------------------------------------------- #

def ensure_parent_dir(path):
    if path is None or not str(path).strip():
        return
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent, exist_ok=True)



def prepare_diagram(args, status_stream=None):
    """
    Parse, validate, lay out, transform, and audit the DT diagram.

    This function does not save files and does not draw on a figure.  It is used
    by both the command-line pipeline and the live GUI preview.
    """
    out = status_stream if status_stream is not None else io.StringIO()

    comps = parse_dt(args.dt)
    model = build_model(comps, negative_even=args.negative_even)
    crossing_ids = resolve_crossing_ids(
        model,
        crossing_order=getattr(args, "crossing_order", None),
        crossing_map=getattr(args, "crossing_map", None),
    )

    G = build_gadget_graph(model)
    tutte_opts = {
        "shape": getattr(args, "tutte_shape", "ellipse"),
        "aspect": getattr(args, "tutte_aspect", 1.8),
        "corner_radius": getattr(args, "tutte_corner_radius", 0.25),
        "decompress": getattr(args, "tutte_decompress", 0.0),
        "com_expand": getattr(args, "tutte_com_expand", 0.0),
        "auto_aspect": bool(getattr(args, "tutte_auto_aspect", True)),
        # Shaped-tutte orientation is now handled inside the layout solve: the
        # boundary shape's aspect axis is tilted 'tutte orient' degrees from the
        # diagram's intrinsic (circular-Tutte) PCA axis.  'auto orient' frames the
        # PCA axis horizontally so the shape appears tilted by that angle.
        "orient": getattr(args, "tutte_orient", 0.0),
        "auto_orient": bool(getattr(args, "tutte_auto_orient", True)),
        # holed-tutte controls.
        "hole_ratio": getattr(args, "hole_ratio", 0.4),
        "hole_swap": bool(getattr(args, "hole_swap", False)),
        "ring_tilt": getattr(args, "ring_tilt", 90.0),
        "invert_ring": bool(getattr(args, "invert_ring", False)),
        "ring_equalize": getattr(args, "ring_equalize", 0.0),
        # sphere-stereo needs the parsed model to trace the strand arcs.
        "model": model,
    }
    tutte_guides = {}
    P = compute_positions(G, args.layout, tutte_opts=tutte_opts, meta_out=tutte_guides)
    # V5.0: optional planarity-guarded relaxation that evens out the node/edge
    # distribution (0 passes = off), then the minimum-separation nudge.
    _relax_passes = int(getattr(args, "relax_passes", 0) or 0)
    if _relax_passes > 0:
        P = relax_planar_layout(
            P, G, passes=_relax_passes,
            strength=float(getattr(args, "relax_strength", 0.5) or 0.0),
            pinned=(tutte_guides or {}).get("pinned_nodes"),
        )
    P = nudge_min_separation(P, G, getattr(args, "min_sep", 0.0))
    # Guide geometry (shape outline / PCA axis) must ride through the same drawing
    # transform as the node positions, using the pre-transform node centroid.
    guide_center = None
    if P:
        guide_center = np.mean(list(P.values()), axis=0)
    P = transform_positions(P, args.y_direction, args.rotate)
    if tutte_guides and guide_center is not None:
        tutte_guides = _transform_tutte_guides(
            tutte_guides, guide_center, args.y_direction, args.rotate
        )
    else:
        tutte_guides = {}
    centers = crossing_centers(model, P)

    false = audit_false_crossings(model, P, centers)
    used_layout = args.layout
    # The chosen layout is kept as-is even when it introduces false crossings
    # (there is no automatic fallback to 'planar').  False crossings are drawn
    # with over/under gaps and flagged in red on the diagram, so the requested
    # layout -- including a "failed" kamada -- is preserved for inspection.
    if false > 0:
        out.write(
            "[warn] layout '%s' has %d false crossing(s); kept as requested. "
            "They are drawn with gaps and flagged in red on the diagram.\n"
            % (args.layout, false)
        )
    else:
        out.write("[ok] audit: no false crossings.\n")

    out.write(
        "[info] %d components, %d crossings, 2n = %d.\n" %
        (len(comps), len(model["crossings"]), model["twon"])
    )

    return {
        "comps": comps,
        "model": model,
        "crossing_ids": crossing_ids,
        "graph": G,
        "P": P,
        "centers": centers,
        "false_crossings": false,
        "used_layout": used_layout,
        "tutte_guides": tutte_guides,
        "aspect_value": tutte_guides.get("aspect_value") if tutte_guides else None,
    }


def _draw_tutte_guides(ax, guides, show_outline=False, show_pca=False):
    """
    Overlay layout guides on an axis.  ``guides`` is the metadata dict produced by
    ``shaped_tutte_layout`` or ``holed_tutte_layout`` (already in drawing
    coordinates).  For shaped-tutte: the boundary outline + aspect axis (outline
    toggle) and the intrinsic PCA axis (pca toggle).  For holed-tutte: the outer
    and inner outlines (outline toggle) and the mid-ring closed 'principal curved
    axis' (pca toggle).
    """
    if not guides:
        return

    def _plot_loop(arr, **kw):
        arr = np.asarray(arr, float)
        if len(arr) >= 2:
            loop = np.vstack([arr, arr[:1]])
            ax.plot(loop[:, 0], loop[:, 1], clip_on=False, **kw)

    if guides.get("kind") == "holed":
        if show_outline:
            for key in ("boundary_outer", "boundary_inner"):
                b = guides.get(key)
                if b is not None:
                    _plot_loop(b, color="#1f77b4", lw=1.2, ls=(0, (6, 4)),
                               alpha=0.85, zorder=15)
        if show_pca:
            med = guides.get("medial")
            if med is not None:
                _plot_loop(med, color="#ff7f0e", lw=1.5, ls=(0, (2, 2)),
                           alpha=0.95, zorder=16)
        return

    if show_outline:
        b = guides.get("boundary")
        if b is not None:
            _plot_loop(b, color="#1f77b4", lw=1.2, ls=(0, (6, 4)),
                       alpha=0.85, zorder=15)
        sa = guides.get("shape_axis")
        if sa is not None and len(sa) >= 2:
            sa = np.asarray(sa, float)
            ax.plot(sa[:, 0], sa[:, 1], color="#1f77b4", lw=1.3,
                    ls="-", alpha=0.85, zorder=15, clip_on=False)
    if show_pca:
        pa = guides.get("pca_axis")
        if pa is not None and len(pa) >= 2:
            pa = np.asarray(pa, float)
            ax.plot(pa[:, 0], pa[:, 1], color="#ff7f0e", lw=1.5,
                    ls=(0, (2, 2)), alpha=0.95, zorder=16, clip_on=False)


def render_prepared_diagram(ax, data, args):
    """Draw already prepared geometry on an existing Matplotlib axis."""
    ax.clear()
    label_coords, crossing_xy = render_diagram(
        ax,
        data["model"],
        data["P"],
        data["centers"],
        crossing_ids=data["crossing_ids"],
        show_crossing_ids=args.show_crossing_ids,
        color_crossing_ids_by_overstrand=args.color_crossing_ids_by_overstrand,
        label_fontsize=args.font_size,
        crossing_id_fontsize=args.crossing_id_font_size,
        lw=args.line_width,
        gap_frac=args.gap_frac,
        show_labels=not args.hide_labels,
        arrows=not args.no_arrows,
    )
    _draw_tutte_guides(
        ax,
        data.get("tutte_guides"),
        show_outline=bool(getattr(args, "show_tutte_outline", False)),
        show_pca=bool(getattr(args, "show_tutte_pca", False)),
    )
    ax.set_aspect("equal")
    ax.axis("off")
    if args.title:
        ax.set_title(args.title, fontsize=11, fontfamily=DIAGRAM_FONT_FAMILY)
    _maximize_axis_in_figure(ax, has_title=bool(args.title))
    _tighten_axis_to_content(ax, pad_frac=0.08)
    _disable_figure_clipping(ax.figure)
    return label_coords, crossing_xy


def compute_diagram_state(args, status_stream=None):
    """Backward-compatible name used by GUI event handlers."""
    return prepare_diagram(args, status_stream=status_stream)


def compute_diagram_objects(args, status_stream=None):
    """Name used by the live GUI; returns the prepared diagram state dict."""
    return prepare_diagram(args, status_stream=status_stream)


def render_state_on_axis(ax, state, args):
    """Name used by the live GUI; render a prepared state on an axis."""
    return render_prepared_diagram(ax, state, args)


def render_figure(
    fig,
    model,
    P,
    centers,
    show_crossing_ids=True,
    color_crossing_ids_by_overstrand=True,
    crossing_ids=None,
    title=None,
    label_fontsize=7.0,
    crossing_id_fontsize=6.0,
    line_width=2.0,
    gap_frac=0.025,
    show_labels=True,
    arrows=True,
):
    """Render a diagram into an existing Matplotlib Figure.

    Used by the live Tk preview and by table-only saves, where we need the
    label/crossing coordinates without necessarily writing a new image.
    """
    fig.clear()
    ax = fig.add_subplot(111)
    label_coords, crossing_xy = render_diagram(
        ax,
        model,
        P,
        centers,
        crossing_ids=crossing_ids,
        show_crossing_ids=show_crossing_ids,
        color_crossing_ids_by_overstrand=color_crossing_ids_by_overstrand,
        label_fontsize=label_fontsize,
        crossing_id_fontsize=crossing_id_fontsize,
        lw=line_width,
        gap_frac=gap_frac,
        show_labels=show_labels,
        arrows=arrows,
    )
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=11, fontfamily=DIAGRAM_FONT_FAMILY)
    _maximize_axis_in_figure(ax, has_title=bool(title))
    _tighten_axis_to_content(ax, pad_frac=0.08)
    _disable_figure_clipping(fig)
    return label_coords, crossing_xy


def _render_table_coordinates_for_pipeline(data, args):
    """Render off-screen only to recover 2D label/crossing coordinates for CSV."""
    tmp_fig, tmp_ax = plt.subplots(figsize=(float(args.figsize), float(args.figsize)))
    try:
        label_coords, crossing_xy = render_prepared_diagram(tmp_ax, data, args)
    finally:
        plt.close(tmp_fig)
    return label_coords, crossing_xy


def run_pipeline(args, status_stream=None):
    out = status_stream if status_stream is not None else sys.stderr
    data = prepare_diagram(args, status_stream=out)

    label_coords = None
    crossing_xy = None

    if not getattr(args, "no_image", False):
        label_coords, crossing_xy = draw(
            data["model"],
            data["P"],
            data["centers"],
            args.output,
            dpi=args.dpi,
            show_crossing_ids=args.show_crossing_ids,
            color_crossing_ids_by_overstrand=args.color_crossing_ids_by_overstrand,
            crossing_ids=data["crossing_ids"],
            title=args.title,
            label_fontsize=args.font_size,
            crossing_id_fontsize=args.crossing_id_font_size,
            line_width=args.line_width,
            gap_frac=args.gap_frac,
            figsize=args.figsize,
            show_labels=not args.hide_labels,
            arrows=not args.no_arrows,
            dt_code=args.dt,
            layout=args.layout,
            tutte_guides=data.get("tutte_guides"),
            show_tutte_outline=bool(getattr(args, "show_tutte_outline", False)),
            show_tutte_pca=bool(getattr(args, "show_tutte_pca", False)),
        )
        out.write("[ok] wrote %s\n" % args.output)

    if args.table:
        if label_coords is None or crossing_xy is None:
            label_coords, crossing_xy = _render_table_coordinates_for_pipeline(data, args)
        write_table(
            data["model"],
            label_coords,
            crossing_xy,
            args.table,
            crossing_ids=data["crossing_ids"],
        )
        out.write("[ok] wrote %s\n" % args.table)

    if getattr(args, "write_xyz", True) and args.xyz_output:
        n_points = write_spherical_xyz(
            data["model"],
            data["P"],
            data["centers"],
            args.xyz_output,
            G=data["graph"],
            sphere_radius=args.sphere_radius,
            sphere_extent=args.sphere_extent,
            sphere_offset=args.sphere_offset,
            sphere_bump_frac=args.sphere_bump_frac,
            xyz_spacing=args.xyz_spacing,
            close_components=args.xyz_close_components,
            decimals=args.xyz_decimals,
            sphere_layout=args.sphere_layout,
            sphere_crossing_angle=args.sphere_crossing_angle,
            direct_connecting=args.direct_connecting,
            xyz_final_smooth=args.xyz_final_smooth,
            xyz_smooth_window=args.xyz_smooth_window,
            xyz_smooth_passes=args.xyz_smooth_passes,
            surface_shape=getattr(args, "surface_shape", "ellipsoid"),
            surface_auto_orient=bool(getattr(args, "surface_auto_orient", True)),
            surface_auto_aspect=bool(getattr(args, "surface_auto_aspect", True)),
            surface_aspect=getattr(args, "surface_aspect", 1.6),
            surface_tube=getattr(args, "surface_tube", 0.35),
            surface_orient=getattr(args, "surface_orient", 0.0),
            surface_tilt=getattr(args, "surface_tilt", 0.0),
            xyz_repair=not getattr(args, "no_xyz_repair", False),
            xyz_clearance=(getattr(args, "xyz_clearance", 0.0) or None),
            audit_dt=getattr(args, "dt", None),
            status_stream=out,
        )
        out.write("[ok] wrote %s (%d xyz rows)\n" % (args.xyz_output, n_points))

        if getattr(args, "save_projections", False):
            _skel = {}
            _comps3 = build_spherical_xyz_components(
                data["model"], data["P"], data["centers"], G=data["graph"],
                sphere_radius=args.sphere_radius,
                sphere_extent=args.sphere_extent,
                sphere_offset=args.sphere_offset,
                sphere_bump_frac=args.sphere_bump_frac,
                xyz_spacing=args.xyz_spacing,
                sphere_layout=args.sphere_layout,
                sphere_crossing_angle=args.sphere_crossing_angle,
                direct_connecting=args.direct_connecting,
                xyz_final_smooth=args.xyz_final_smooth,
                xyz_smooth_window=args.xyz_smooth_window,
                xyz_smooth_passes=args.xyz_smooth_passes,
                surface_shape=getattr(args, "surface_shape", "ellipsoid"),
                surface_auto_orient=bool(getattr(args, "surface_auto_orient", True)),
                surface_auto_aspect=bool(getattr(args, "surface_auto_aspect", True)),
                surface_aspect=getattr(args, "surface_aspect", 1.6),
                surface_tube=getattr(args, "surface_tube", 0.35),
                surface_orient=getattr(args, "surface_orient", 0.0),
                surface_tilt=getattr(args, "surface_tilt", 0.0),
                skeleton_out=_skel,
                skeleton_grid=getattr(args, "proj_grid_density", 24),
            )
            _views = _parse_projection_views(getattr(args, "proj_views", "current,top"))
            _pe_ = getattr(args, "proj_elev", 25.0)
            _pa_ = getattr(args, "proj_azim", -60.0)
            if getattr(args, "proj_auto_view", False):
                _pe_, _pa_ = auto_projection_view(_comps3)
                out.write("[proj] auto view: elev=%.1f azim=%.1f\n" % (_pe_, _pa_))
            _written = save_projection_images(
                _comps3, args.xyz_output, _views,
                elev=_pe_,
                azim=_pa_,
                roll=getattr(args, "proj_roll", 0.0),
                line_width=getattr(args, "proj_line_width", 2.2),
                skeleton=_skel,
                show_skeleton=bool(getattr(args, "proj_skeleton", False)),
                skeleton_ids=bool(getattr(args, "proj_skeleton_ids", False)),
                figsize=getattr(args, "figsize", 8.0),
                dpi=getattr(args, "dpi", 200),
                depth_fade=getattr(args, "proj_depth_fade", 0.55),
                perspective=getattr(args, "proj_perspective", 0.0),
                skeleton_only=bool(getattr(args, "proj_skeleton_only", False)),
                grid_color=getattr(args, "proj_grid_color", "0.82"),
                grid_lw=getattr(args, "proj_grid_lw", 0.4),
                show_chords=not bool(getattr(args, "proj_hide_chords", False)),
                show_crossing_dots=bool(getattr(args, "proj_crossing_dots", False)),
                fixed_grid=bool(getattr(args, "proj_fixed_grid", False)),
                fixed_grid_view=(
                    getattr(args, "proj_elev", 25.0),
                    getattr(args, "proj_azim", -60.0),
                    getattr(args, "proj_roll", 0.0),
                ),
            )
            for _pth in _written:
                out.write("[ok] wrote %s\n" % _pth)

    if getattr(args, "no_image", False) and not args.table and not (getattr(args, "write_xyz", True) and args.xyz_output):
        out.write("[warn] nothing was written; enable an image, table, or XYZ output.\n")

    return 0

def build_arg_parser():
    ap = argparse.ArgumentParser(
        description=(
            "Draw a planar link diagram from a signed DT code while preserving "
            "the original traversal labels. With no arguments, GUI mode opens."
        )
    )
    ap.add_argument("--gui", action="store_true", help="Open GUI mode.")
    ap.add_argument(
        "--dt",
        default=None,
        help=(
            "Signed DT code, e.g. "
            "'DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]'"
        ),
    )
    ap.add_argument(
        "--output",
        default="link_diagram.svg",
        help="Output image path (.svg / .pdf / .png). Default: link_diagram.svg",
    )
    ap.add_argument(
        "--table",
        default=None,
        help="Optional CSV path for the crossing/label mapping table.",
    )
    ap.add_argument(
        "--negative-even",
        choices=["over", "under"],
        default="over",
        help=(
            "What a negative even DT label means for the even visit. "
            "'over' is the default Sage/KnotTheory-style convention."
        ),
    )
    ap.add_argument(
        "--crossing-order",
        default=None,
        help=(
            "Displayed crossing IDs ordered by odd labels. Example: "
            "'c1 c7 c14 c12 c3 c6 c9 c5 c11 c13 c4 c2 c10 c8' "
            "means odd labels 1,3,5,... map to those crossing IDs."
        ),
    )
    ap.add_argument(
        "--crossing-map",
        default=None,
        help=(
            "Alternative explicit map: displayed crossing ID = odd label, e.g. "
            "'c1=1,c7=3,c14=5'. Do not combine with --crossing-order."
        ),
    )
    ap.add_argument(
        "--show-crossing-ids",
        dest="show_crossing_ids",
        action="store_true",
        default=True,
        help="Draw crossing IDs at each crossing. Default: on.",
    )
    ap.add_argument(
        "--hide-crossing-ids",
        dest="show_crossing_ids",
        action="store_false",
        help="Do not draw crossing IDs.",
    )
    ap.add_argument(
        "--color-crossing-ids-by-overstrand",
        dest="color_crossing_ids_by_overstrand",
        action="store_true",
        default=True,
        help="Color each crossing ID using the color of the over-strand component. Default: on.",
    )
    ap.add_argument(
        "--no-color-crossing-ids-by-overstrand",
        dest="color_crossing_ids_by_overstrand",
        action="store_false",
        help="Use a neutral color for crossing IDs instead of the over-strand color.",
    )
    ap.add_argument(
        "--y-direction",
        choices=["top-to-bottom", "bottom-to-top"],
        default="top-to-bottom",
        help="Final y-coordinate direction. Default: top-to-bottom.",
    )
    ap.add_argument(
        "--rotate",
        type=float,
        default=0.0,
        help="Rotate the final scheme by this many degrees after y-direction handling.",
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI for raster output. Default: 200.",
    )
    ap.add_argument(
        "--layout",
        choices=["tutte", "shaped-tutte", "holed-tutte", "sphere-stereo", "planar", "spring", "kamada"],
        default="tutte",
        help=(
            "Layout engine. 'tutte' is default. 'shaped-tutte' pins the Tutte "
            "boundary to a chosen convex shape (see --tutte-shape). 'holed-tutte' "
            "pins two faces (outer + auto-picked central hole) to the outer/inner "
            "outlines of a holed shape and solves the ring. 'sphere-stereo' lays "
            "the diagram on the unit sphere and stereographically projects it "
            "from the pole farthest from all strands (best for links whose rings "
            "lie in near-orthogonal planes, e.g. Edwards-Venn AM5). If false crossings are "
            "detected, the chosen layout is kept and the artifacts are highlighted."
        ),
    )
    ap.add_argument(
        "--tutte-shape",
        choices=list(TUTTE_SHAPES),
        default="ellipse",
        metavar="{%s}" % ",".join(TUTTE_SHAPE_CHOICES),
        help=(
            "Boundary shape for --layout shaped-tutte / holed-tutte: ellipse or "
            "rounded-rectangle (use aspect 1 for a circle, corner radius 0 for a "
            "sharp rectangle). Default: ellipse."
        ),
    )
    ap.add_argument(
        "--tutte-aspect", type=float, default=1.8,
        help=(
            "Boundary aspect ratio (long/short) for shaped-tutte ellipse/"
            "rectangle shapes when auto is off. Default: 1.8."
        ),
    )
    ap.add_argument(
        "--tutte-corner-radius", type=float, default=0.25,
        help=(
            "Corner radius as a fraction (0..1) of the short half-extent for the "
            "shaped-tutte rounded-rectangle shape. Default: 0.25."
        ),
    )
    ap.add_argument(
        "--tutte-decompress", type=float, default=0.0,
        help=(
            "Internal decompression strength for 'tutte' and 'shaped-tutte'. "
            "0 = classic Tutte. Larger values apply boundary-depth edge weights "
            "that push interior structure outward from the boundary so it is not "
            "over-compressed toward the center. Default: 0.0."
        ),
    )
    ap.add_argument(
        "--tutte-com-expand", type=float, default=0.0,
        help=(
            "Separate radial expansion strength about the crossing center of "
            "mass for 'tutte' and 'shaped-tutte'. 0 = off. Expands interior "
            "structure outward from the (density-weighted) crowded-crossing "
            "centroid, tapering to zero at the boundary, so crowded crossings get "
            "more room and sit near the center. Layered on top of "
            "--tutte-decompress. Suggested range 0.0-1.0. Default: 0.0."
        ),
    )
    ap.add_argument(
        "--tutte-auto-orient",
        dest="tutte_auto_orient",
        action="store_true",
        default=True,
        help=(
            "For shaped-tutte: frame the diagram so its intrinsic (circular-Tutte) "
            "PCA elongation axis lies along the view x-axis. Default: on. The shape "
            "outline then appears tilted by --tutte-orient. When off, the shape's "
            "long axis is left horizontal instead. Works for every boundary shape."
        ),
    )
    ap.add_argument(
        "--no-tutte-auto-orient",
        dest="tutte_auto_orient",
        action="store_false",
        help="Frame the shaped-tutte shape long axis horizontally instead of the PCA axis.",
    )
    ap.add_argument(
        "--tutte-orient", type=float, default=0.0,
        help=(
            "Shaped-tutte shape tilt in degrees: the angle by which the boundary "
            "shape's aspect (long) axis is tilted away from the diagram's intrinsic "
            "PCA elongation axis. 0 = shape long axis along the PCA axis. This "
            "re-stretches the layout in a new direction (it does not merely spin the "
            "finished picture). Only meaningful for non-circular shapes. Default: 0."
        ),
    )
    ap.add_argument(
        "--show-tutte-outline",
        dest="show_tutte_outline",
        action="store_true",
        default=False,
        help="Overlay the shaped-tutte boundary outline and its aspect axis.",
    )
    ap.add_argument(
        "--show-tutte-pca",
        dest="show_tutte_pca",
        action="store_true",
        default=False,
        help="Overlay the diagram's intrinsic PCA elongation axis.",
    )
    ap.add_argument(
        "--tutte-auto-aspect",
        dest="tutte_auto_aspect",
        action="store_true",
        default=True,
        help=(
            "For shaped-tutte: derive the boundary aspect ratio from the "
            "diagram's own elongation (PCA). Default: on. When off, "
            "--tutte-aspect is used."
        ),
    )
    ap.add_argument(
        "--no-tutte-auto-aspect",
        dest="tutte_auto_aspect",
        action="store_false",
        help="Disable shaped-tutte auto aspect; use --tutte-aspect instead.",
    )
    ap.add_argument(
        "--hole-ratio", type=float, default=0.4,
        help=(
            "For --layout holed-tutte: inner (hole) outline size as a fraction of "
            "the outer outline, in (0,1). Smaller = bigger ring. Default: 0.4."
        ),
    )
    ap.add_argument(
        "--hole-swap",
        dest="hole_swap",
        action="store_true",
        default=False,
        help="For holed-tutte: swap which face is pinned to the outer vs inner outline.",
    )
    ap.add_argument(
        "--invert-ring",
        dest="invert_ring",
        action="store_true",
        default=False,
        help=(
            "For holed-tutte: turn the ring inside-out (reflect each crossing's "
            "radius about the ring mid-line) so the inner boundary ends up outside "
            "and the outer boundary inside."
        ),
    )
    ap.add_argument(
        "--ring-tilt", type=float, default=90.0,
        help=(
            "For holed-tutte: how far the closed principal ring is tilted UP to face "
            "the viewer, 0-90 degrees. 90 = ring faces you -> the flat annulus / "
            "donut 'top view' (curved axis a loop in the canvas plane). 0 = ring "
            "edge-on -> 'side view' band with the curved axis perpendicular to the "
            "canvas. Default: 90."
        ),
    )
    ap.add_argument(
        "--min-sep", type=float, default=0.0,
        help=(
            "Minimum separation between non-incident strand pieces, as a fraction "
            "of the diagram span. A post-layout relaxation pushes closer pieces "
            "apart (with a spring back to the layout). 0 = off. Try 0.02-0.05."
        ),
    )
    ap.add_argument(
        "--ring-equalize", type=float, default=0.0,
        help=(
            "holed-tutte only: strength 0..1 of a radial histogram equalization "
            "across the ring width (angles kept), spreading crowded strands "
            "evenly between the rims. 0 = off (classic harmonic ring)."
        ),
    )
    ap.add_argument(
        "--relax-passes", type=int, default=0,
        help=(
            "Planarity-guarded relaxation passes applied after any 2D layout "
            "(rejected if a pass would introduce a straight-edge crossing). "
            "0 = off; try 5-30 for congested diagrams."
        ),
    )
    ap.add_argument(
        "--relax-strength", type=float, default=0.5,
        help="Step size 0..1 of each relaxation pass. Default: 0.5.",
    )
    ap.add_argument(
        "--save-projections", action="store_true",
        help=(
            "Also write illustrative 2D orthographic projections of the 3D "
            "curve (SVG+PNG per view in --proj-views) next to the XYZ output."
        ),
    )
    ap.add_argument(
        "--proj-views", default="current,top",
        help=(
            "Comma-separated projection views for --save-projections: "
            "current, top, front, side. Default: current,top."
        ),
    )
    ap.add_argument("--proj-elev", type=float, default=25.0,
                    help="Elevation (deg) of the 'current' projection view. Default: 25.")
    ap.add_argument("--proj-azim", type=float, default=-60.0,
                    help="Azimuth (deg) of the 'current' projection view. Default: -60.")
    ap.add_argument("--proj-roll", type=float, default=0.0,
                    help="Roll (deg) about the view direction. Default: 0.")
    ap.add_argument("--proj-line-width", type=float, default=2.2,
                    help="Line width of projected strands (halo scales with it). Default: 2.2.")
    ap.add_argument("--proj-skeleton", action="store_true",
                    help="Overlay the mapped skeleton (crossing anchors + strand chords).")
    ap.add_argument("--proj-skeleton-ids", action="store_true",
                    help="Label skeleton anchors with crossing IDs.")
    ap.add_argument("--proj-skeleton-only", action="store_true",
                    help="Projections show ONLY the mapped framework (no smooth strands).")
    ap.add_argument("--proj-depth-fade", type=float, default=0.55,
                    help="Depth-cue transparency strength 0..1 for projections. Default: 0.55.")
    ap.add_argument("--proj-perspective", type=float, default=0.0,
                    help="0 = orthographic; else perspective camera distance in scene radii (try 2-4).")
    ap.add_argument("--proj-grid-density", type=int, default=24,
                    help="Meridian count of the surface wireframe (0 = no mesh). Default: 24.")
    ap.add_argument("--proj-grid-color", default="0.82",
                    help="Matplotlib colour of the surface wireframe. Default: 0.82.")
    ap.add_argument("--proj-grid-lw", type=float, default=0.4,
                    help="Line width of the surface wireframe. Default: 0.4.")
    ap.add_argument("--proj-hide-chords", action="store_true",
                    help="Hide the neutral framework chords between crossing anchors.")
    ap.add_argument("--proj-crossing-dots", action="store_true",
                    help="Show black crossing-anchor dots in the mapped framework (default: hidden).")
    ap.add_argument("--proj-fixed-grid", action="store_true",
                    help="Keep the mapped-surface grid fixed to the current view while strands rotate.")
    ap.add_argument("--proj-auto-view", action="store_true",
                    help="Pick the 'current' view automatically (maximizes the worst ring's projected roundness).")
    ap.add_argument("--title", default=None, help="Optional title for the figure.")
    ap.add_argument("--font-size", type=float, default=7.0, help="DT label font size. Default: 7.")
    ap.add_argument(
        "--crossing-id-font-size",
        type=float,
        default=6.0,
        help="Crossing ID font size. Default: 6.",
    )
    ap.add_argument("--line-width", type=float, default=2.0, help="Strand line width. Default: 2.0.")
    ap.add_argument(
        "--gap-frac",
        type=float,
        default=0.025,
        help="Under-strand gap size as a fraction of diagram span. Default: 0.025.",
    )
    ap.add_argument(
        "--figsize",
        type=float,
        default=12.0,
        help="Square Matplotlib figure size in inches. Default: 12.",
    )
    ap.add_argument("--hide-labels", action="store_true", help="Do not draw DT labels.")
    ap.add_argument("--no-arrows", action="store_true", help="Do not draw orientation arrows.")
    ap.add_argument(
        "--xyz-output",
        "--xyz",
        dest="xyz_output",
        default="link_sphere.xyz",
        help="Output path for the spherical x y z coordinate file. Default: link_sphere.xyz.",
    )
    ap.add_argument(
        "--no-xyz",
        dest="write_xyz",
        action="store_false",
        default=True,
        help="Do not write the spherical XYZ coordinate file.",
    )
    ap.add_argument(
        "--no-image",
        action="store_true",
        help="Do not save the 2D image; useful when only the XYZ file is needed.",
    )
    ap.add_argument(
        "--sphere-layout",
        choices=["spherical-kamada", "shaped-kamada", "stereographic", "stereo-safe"],
        default="spherical-kamada",
        help=(
            "XYZ sphere layout. 'spherical-kamada' spreads the graph directly "
            "on S^2 and is the default; 'shaped-kamada' warps that spherical "
            "construction onto a shaped surface (see --surface-shape); "
            "'stereographic' uses the 2D diagram and inverse stereographic "
            "projection; 'stereo-safe' preserves the audited 2D angular "
            "positions and is the correct-by-construction choice after a "
            "Kamada topology-audit failure."
        ),
    )
    ap.add_argument(
        "--surface-shape",
        choices=list(SURFACE_SHAPES),
        default="ellipsoid",
        help=(
            "Target surface for --sphere-layout shaped-kamada: ellipsoid, "
            "cylinder, or torus. Default: ellipsoid."
        ),
    )
    ap.add_argument(
        "--surface-auto-orient",
        dest="surface_auto_orient",
        action="store_true",
        default=True,
        help=(
            "For shaped-kamada: auto-align the surface's primary axis to the "
            "diagram (3D PCA). Default: on. When off, --surface-orient spins the "
            "mapping about that axis relative to the PCA axis."
        ),
    )
    ap.add_argument(
        "--no-surface-auto-orient",
        dest="surface_auto_orient",
        action="store_false",
        help="Disable shaped-kamada auto orientation; use --surface-orient instead.",
    )
    ap.add_argument(
        "--surface-orient", type=float, default=0.0,
        help=(
            "Manual shaped-kamada spin in degrees about the surface's primary "
            "axis, measured from the diagram's PCA axis. Used only when "
            "--no-surface-auto-orient. Default: 0."
        ),
    )
    ap.add_argument(
        "--surface-tilt", type=float, default=0.0,
        help=(
            "Manual shaped-kamada tilt in degrees that rotates the surface's "
            "primary axis away from the diagram's PCA axis (about the "
            "perpendicular secondary axis). Used only when "
            "--no-surface-auto-orient. Default: 0."
        ),
    )
    ap.add_argument(
        "--surface-auto-aspect",
        dest="surface_auto_aspect",
        action="store_true",
        default=True,
        help=(
            "For shaped-kamada: derive the surface aspect / tube ratio from the "
            "diagram's own 3D shape (PCA). Default: on. When off, "
            "--surface-aspect / --surface-tube are used."
        ),
    )
    ap.add_argument(
        "--no-surface-auto-aspect",
        dest="surface_auto_aspect",
        action="store_false",
        help="Disable shaped-kamada auto aspect; use manual surface values.",
    )
    ap.add_argument(
        "--surface-aspect", type=float, default=1.6,
        help=(
            "Manual surface aspect for shaped-kamada when auto aspect is off: "
            "ellipsoid major/minor, or cylinder length/radius. Ignored by torus. "
            "Default: 1.6."
        ),
    )
    ap.add_argument(
        "--surface-tube", type=float, default=0.35,
        help=(
            "Manual torus tube ratio (minor/major radius) for shaped-kamada when "
            "auto aspect is off. Default: 0.35."
        ),
    )
    ap.add_argument(
        "--direct-connecting",
        dest="direct_connecting",
        action="store_true",
        default=True,
        help=(
            "Only for spherical-kamada XYZ export. Connect each segment directly "
            "between its endpoint crossing layers instead of returning to the "
            "middle sphere layer. Default: on."
        ),
    )
    ap.add_argument(
        "--no-direct-connecting",
        dest="direct_connecting",
        action="store_false",
        help="Disable direct connecting and return connector segments to the middle sphere layer.",
    )
    ap.add_argument(
        "--sphere-radius", type=float, default=50.0,
        help="Base sphere radius in XYZ coordinate units. Default: 50.0.",
    )
    ap.add_argument(
        "--sphere-extent", type=float, default=1.55,
        help=(
            "Only for --sphere-layout stereographic: dimensionless farthest "
            "normalized 2D radius before inverse stereographic projection. "
            "1.0 reaches the equator; larger values use more southern hemisphere. "
            "Default: 1.55."
        ),
    )
    ap.add_argument(
        "--crossing-offset", dest="sphere_offset", type=float, default=5.0,
        help=(
            "Absolute radial over/under displacement near crossings, in XYZ "
            "coordinate units. Over = R + offset; under = R - offset. "
            "Default: 5.0. Use 0 for all points exactly on the base sphere."
        ),
    )
    ap.add_argument("--sphere-offset", dest="sphere_offset", type=float, help=argparse.SUPPRESS)
    ap.add_argument(
        "--sphere-crossing-angle", "--crossing-angle",
        dest="sphere_crossing_angle", type=float, default=15.0,
        help=(
            "Angular half-size of each compact crossing gadget for spherical-kamada "
            "XYZ export, in degrees. Default: 15."
        ),
    )
    ap.add_argument(
        "--sphere-bump-frac", "--bump-fraction",
        dest="sphere_bump_frac", type=float, default=1.0,
        help=(
            "Dimensionless fraction in the range 0 < value <= 1. Without "
            "--direct-connecting, it controls local crossing bump width. With "
            "--direct-connecting, it controls how much of a dip-to-bump or "
            "bump-to-dip connector is used for the smooth layer transition. "
            "Default: 1.0."
        ),
    )
    ap.add_argument(
        "--xyz-spacing", "--xyz-point-spacing", "--point-spacing",
        dest="xyz_spacing", type=float, default=1.8,
        help=(
            "Target separation between adjacent XYZ points, in XYZ coordinate units. "
            "Points are redistributed approximately evenly along each component. "
            "Default: 1.8."
        ),
    )
    ap.add_argument(
        "--xyz-final-smooth",
        dest="xyz_final_smooth",
        action="store_true",
        default=True,
        help=(
            "Apply a final smoothing pass to each written 3D XYZ component. "
            "Default: on. This removes small kink artifacts from stitched "
            "spherical arcs."
        ),
    )
    ap.add_argument(
        "--no-xyz-final-smooth",
        dest="xyz_final_smooth",
        action="store_false",
        help="Disable the final 3D XYZ smoothing pass.",
    )
    ap.add_argument(
        "--xyz-smooth-window",
        type=int,
        default=10,
        help=(
            "Final XYZ smoothing window: number of neighboring points on each "
            "side used by the circular smoothing kernel. Default: 10."
        ),
    )
    ap.add_argument(
        "--xyz-smooth-passes",
        type=int,
        default=5,
        help=(
            "Number of repeated final XYZ smoothing passes. Default: 5. "
            "Use 0 or --no-xyz-final-smooth for no smoothing."
        ),
    )
    ap.add_argument(
        "--xyz-samples", dest="_deprecated_xyz_samples", type=int,
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--xyz-clearance", type=float, default=0.0,
        help="V5.4: minimum 3D distance enforced between non-neighboring "
             "strands (0 = auto: 0.8 * crossing offset). The clearance "
             "repair separates crowded regions into extra radial layers.")
    ap.add_argument(
        "--no-xyz-repair", action="store_true",
        help="V5.4: disable the 3D clearance repair (not recommended; dense "
             "diagrams can otherwise render as the WRONG link).")
    ap.add_argument(
        "--xyz-decimals", type=int, default=9,
        help="Decimal places written in the XYZ file. Default: 9.",
    )
    ap.add_argument(
        "--xyz-close-components",
        dest="xyz_close_components",
        action="store_true",
        default=False,
        help="Repeat the first point at the end of each component block. Default: off.",
    )
    ap.add_argument(
        "--no-xyz-close-components",
        dest="xyz_close_components",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    return ap



def run_gui(initial_args):
    try:
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure
    except Exception as exc:
        sys.stderr.write("[error] GUI mode requires Tkinter: %s\n" % exc)
        return 1

    try:
        root = tk.Tk()
    except Exception as exc:
        sys.stderr.write("[error] Could not start GUI mode: %s\n" % exc)
        return 1

    apply_tk_window_icon(root, tk)
    root.title("draw_dt_original_labelsV5_4")
    root.geometry("1320x860")
    root.minsize(1050, 680)

    main = ttk.Frame(root, padding=8)
    main.grid(row=0, column=0, sticky="nsew")
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    main.columnconfigure(0, weight=3)
    main.columnconfigure(1, weight=2)
    main.rowconfigure(0, weight=1)

    # ------------------------------------------------------------------
    # Shared helpers and the file-path inputs.  These are defined before the
    # left panel because the DT code box and the output-path fields now live at
    # the top of the left column (they can be long), directly above the preview.
    # ------------------------------------------------------------------
    def show_arg_help(key):
        title = key.replace("_", " ").strip().title()
        message = GUI_HELP_TEXT.get(key, "No help text is available for this parameter.")
        if key == "figsize":
            try:
                w = int(canvas_widget.winfo_width())
                h = int(canvas_widget.winfo_height())
                if w > 1 and h > 1:
                    message += (
                        "\n\nCurrent live preview panel: %d x %d px "
                        "(about %.1f x %.1f in at the preview's 100 dpi)."
                        % (w, h, w / 100.0, h / 100.0)
                    )
            except Exception:
                pass
        # Show help in a resizable popup with a read-only, selectable Text widget
        # so long entries (e.g. the DT example codes) fit and can be copied.
        win = tk.Toplevel(root)
        win.title(title)
        try:
            apply_tk_window_icon(win)
        except Exception:
            pass
        win.transient(root)
        # Size the popup to the content; the DT entry gets a wide, tall window.
        lines = message.splitlines() or [""]
        longest = max((len(ln) for ln in lines), default=40)
        width_chars = max(48, min(100, longest + 2))
        height_lines = max(8, min(34, len(lines) + 2))
        if key == "dt":
            width_chars = max(width_chars, 92)
            height_lines = max(height_lines, 26)
        frame = ttk.Frame(win, padding=8)
        frame.pack(fill="both", expand=True)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        txt = tk.Text(frame, wrap="word", width=width_chars, height=height_lines,
                      font=("TkFixedFont",))
        yscroll = ttk.Scrollbar(frame, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=yscroll.set)
        txt.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        txt.insert("1.0", message)
        txt.configure(state="disabled")
        btns = ttk.Frame(win, padding=(8, 0, 8, 8))
        btns.pack(fill="x")
        ttk.Button(btns, text="Close", command=win.destroy).pack(side="right")
        win.bind("<Escape>", lambda _e: win.destroy())
        try:
            win.update_idletasks()
            win.minsize(win.winfo_reqwidth(), win.winfo_reqheight())
        except Exception:
            pass

    def _mk_help_button(parent, key):
        return tk.Button(
            parent,
            text="?",
            width=2,
            bg="#cfeeff",
            activebackground="#aee3ff",
            relief="raised",
            command=lambda k=key: show_arg_help(k),
        )

    output_var = tk.StringVar(value=initial_args.output or "link_diagram.svg")
    xyz_var = tk.StringVar(value=initial_args.xyz_output or "link_sphere.xyz")

    def _derived_table_path():
        """CSV table path derived from the output image path (same name, .csv)."""
        base = output_var.get().strip() or "link_diagram.svg"
        return os.path.splitext(base)[0] + ".csv"

    def browse_output():
        path = filedialog.asksaveasfilename(
            title="Save diagram as",
            defaultextension=".svg",
            filetypes=[
                ("SVG", "*.svg"),
                ("PDF", "*.pdf"),
                ("PNG", "*.png"),
                ("All files", "*.*"),
            ],
        )
        if path:
            output_var.set(path)

    def browse_xyz():
        path = filedialog.asksaveasfilename(
            title="Save spherical XYZ coordinates as",
            defaultextension=".xyz",
            filetypes=[("XYZ coordinate file", "*.xyz"), ("Text", "*.txt"), ("All files", "*.*")],
        )
        if path:
            xyz_var.set(path)

    # ------------------------------------------------------------------
    # Left side: long inputs + live preview + save buttons + status log.
    # ------------------------------------------------------------------
    left = ttk.Frame(main)
    left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))
    left.columnconfigure(0, weight=1)
    left.rowconfigure(1, weight=1)

    # DT code and output-path fields.  Placed above the preview because these
    # values can be very long.
    io_frame = ttk.Frame(left)
    io_frame.grid(row=0, column=0, sticky="ew", pady=(0, 6))
    io_frame.columnconfigure(1, weight=1)

    # Only the DT code is entered here; the Save buttons prompt for the output
    # directory and filename, so the output-image / sphere-XYZ path fields are gone.
    ttk.Label(io_frame, text="Signed DT code").grid(row=0, column=0, sticky="nw", pady=3)
    dt_text = tk.Text(io_frame, height=2, width=48, wrap="word")
    dt_text.grid(row=0, column=1, sticky="ew", pady=3)
    dt_text.insert("1.0", initial_args.dt if initial_args.dt else EXAMPLE_DT)
    _mk_help_button(io_frame, "dt").grid(row=0, column=3, sticky="nw", padx=(5, 0), pady=3)

    # The "Live preview" caption now lives inside the preview canvas (as a
    # figure-level text that survives axis clears) to save vertical space, so the
    # dedicated label row above the canvas is gone.
    preview_frame = ttk.Frame(left, relief="sunken", padding=2)
    preview_frame.grid(row=1, column=0, sticky="nsew", pady=(0, 6))
    preview_frame.columnconfigure(0, weight=1)
    preview_frame.rowconfigure(0, weight=1)

    fig = Figure(figsize=(6.4, 6.4), dpi=100)
    fig.text(
        0.012, 0.985, "Live preview",
        ha="left", va="top",
        fontsize=11, fontweight="bold", color="0.35",
    )
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        "Enter or edit parameters on the right.\nThe diagram updates automatically.",
        ha="center",
        va="center",
        transform=ax.transAxes,
    )
    canvas = FigureCanvasTkAgg(fig, master=preview_frame)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(row=0, column=0, sticky="nsew")

    save_bar = ttk.Frame(left)
    save_bar.grid(row=2, column=0, sticky="ew", pady=(0, 4))
    for j in range(7):
        save_bar.columnconfigure(j, weight=1)

    preview_control_bar = ttk.Frame(left)
    preview_control_bar.grid(row=3, column=0, sticky="ew", pady=(0, 6))
    for j in range(7):
        preview_control_bar.columnconfigure(j, weight=1)

    log_text = tk.Text(left, height=8, state="disabled", wrap="word")
    log_text.grid(row=4, column=0, sticky="nsew")

    def set_log(text):
        log_text.configure(state="normal")
        log_text.delete("1.0", "end")
        log_text.insert("1.0", text)
        log_text.see("end")
        log_text.configure(state="disabled")

    # ------------------------------------------------------------------
    # Right side: scrollable parameter panel.
    # ------------------------------------------------------------------
    right_outer = ttk.Frame(main)
    right_outer.grid(row=0, column=1, sticky="nsew")
    right_outer.columnconfigure(0, weight=1)
    right_outer.rowconfigure(0, weight=1)

    # The parameters are split into two tabs: 2D diagram settings and 3D XYZ
    # settings, so only the relevant fields are shown at a time.  Each tab is an
    # independently scrollable panel.
    param_nb = ttk.Notebook(right_outer)
    param_nb.grid(row=0, column=0, sticky="nsew")

    def _make_scroll_tab(title):
        outer = ttk.Frame(param_nb)
        param_nb.add(outer, text=title)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(0, weight=1)
        cv = tk.Canvas(outer, borderwidth=0, highlightthickness=0)
        sb = ttk.Scrollbar(outer, orient="vertical", command=cv.yview)
        cv.configure(yscrollcommand=sb.set)
        cv.grid(row=0, column=0, sticky="nsew")
        sb.grid(row=0, column=1, sticky="ns")
        inner = ttk.Frame(cv, padding=(0, 4, 6, 0))
        win = cv.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>",
                   lambda _e, c=cv: c.configure(scrollregion=c.bbox("all")))
        cv.bind("<Configure>",
                lambda e, c=cv, w=win: c.itemconfigure(w, width=e.width))

        # Mouse-wheel scrolling only while the pointer is over this panel, so it
        # does not fight the 2D preview canvas or the other tab.
        def _wheel(event, c=cv):
            if getattr(event, "num", None) == 4:
                c.yview_scroll(-1, "units")
            elif getattr(event, "num", None) == 5:
                c.yview_scroll(1, "units")
            else:
                d = getattr(event, "delta", 0)
                if d:
                    c.yview_scroll(-1 if d > 0 else 1, "units")
            return "break"

        def _bind(_e=None, c=cv):
            c.bind_all("<MouseWheel>", _wheel)
            c.bind_all("<Button-4>", _wheel)
            c.bind_all("<Button-5>", _wheel)

        def _unbind(_e=None, c=cv):
            c.unbind_all("<MouseWheel>")
            c.unbind_all("<Button-4>")
            c.unbind_all("<Button-5>")

        for _w in (cv, inner, sb):
            _w.bind("<Enter>", _bind)
            _w.bind("<Leave>", _unbind)

        inner.columnconfigure(1, weight=1)
        inner.columnconfigure(2, weight=0)
        inner.columnconfigure(3, weight=0)
        return inner

    tab_2d = _make_scroll_tab("2D diagram")
    tab_3d = _make_scroll_tab("3D XYZ")
    tab_view3d = _make_scroll_tab("3D view")

    # 'settings' points at the tab that new fields are added to; it is switched to
    # the 3D tab before the sphere-XYZ block below.  add_help_button / add_entry
    # read it lazily, so they always target the current tab.
    settings = tab_2d

    def add_help_button(grid_row, key):
        btn = _mk_help_button(settings, key)
        btn.grid(row=grid_row, column=3, sticky="w", padx=(5, 0), pady=2)
        return btn

    # The Signed DT code box and the output-path fields now live at the top of
    # the left column (above the preview), so the parameter panel starts with the
    # drawing options.
    row = 1

    neg_var = tk.StringVar(value=initial_args.negative_even)
    layout_var = tk.StringVar(value=initial_args.layout)
    ydir_var = tk.StringVar(value=initial_args.y_direction)
    rotate_var = tk.StringVar(value=str(initial_args.rotate))
    font_var = tk.StringVar(value=str(initial_args.font_size))
    cid_font_var = tk.StringVar(value=str(initial_args.crossing_id_font_size))
    lw_var = tk.StringVar(value=str(initial_args.line_width))
    gap_var = tk.StringVar(value=str(initial_args.gap_frac))
    figsize_var = tk.StringVar(value=str(initial_args.figsize))
    title_var = tk.StringVar(value=initial_args.title or "")
    sphere_layout_var = tk.StringVar(value=getattr(initial_args, "sphere_layout", "spherical-kamada"))
    direct_connecting_var = tk.BooleanVar(value=bool(getattr(initial_args, "direct_connecting", True)))
    sphere_radius_var = tk.StringVar(value=str(initial_args.sphere_radius))
    sphere_extent_var = tk.StringVar(value=str(initial_args.sphere_extent))
    sphere_offset_var = tk.StringVar(value=str(initial_args.sphere_offset))
    sphere_bump_var = tk.StringVar(value=str(initial_args.sphere_bump_frac))
    sphere_angle_var = tk.StringVar(value=str(getattr(initial_args, "sphere_crossing_angle", 15.0)))
    xyz_spacing_var = tk.StringVar(value=str(initial_args.xyz_spacing))
    xyz_final_smooth_var = tk.BooleanVar(value=bool(getattr(initial_args, "xyz_final_smooth", True)))
    xyz_smooth_window_var = tk.StringVar(value=str(getattr(initial_args, "xyz_smooth_window", 10)))
    xyz_smooth_passes_var = tk.StringVar(value=str(getattr(initial_args, "xyz_smooth_passes", 5)))
    xyz_clearance_var = tk.StringVar(value=str(getattr(initial_args, "xyz_clearance", 0.0)))
    xyz_repair_var = tk.BooleanVar(value=not bool(getattr(initial_args, "no_xyz_repair", False)))
    xyz_decimals_var = tk.StringVar(value=str(initial_args.xyz_decimals))
    xyz_close_var = tk.BooleanVar(value=bool(initial_args.xyz_close_components))

    show_ids_var = tk.BooleanVar(value=bool(initial_args.show_crossing_ids))
    color_ids_var = tk.BooleanVar(
        value=bool(initial_args.color_crossing_ids_by_overstrand)
    )
    hide_labels_var = tk.BooleanVar(value=bool(initial_args.hide_labels))
    no_arrows_var = tk.BooleanVar(value=bool(initial_args.no_arrows))

    # V4.0/V4.1 shaped-tutte (2D) and shaped-kamada (3D) controls.  Auto orient
    # and auto aspect are independent.
    tutte_shape_var = tk.StringVar(value=getattr(initial_args, "tutte_shape", "ellipse"))
    tutte_aspect_var = tk.StringVar(value=str(getattr(initial_args, "tutte_aspect", 1.8)))
    tutte_corner_var = tk.StringVar(value=str(getattr(initial_args, "tutte_corner_radius", 0.25)))
    tutte_decompress_var = tk.StringVar(value=str(getattr(initial_args, "tutte_decompress", 0.0)))
    tutte_com_expand_var = tk.StringVar(value=str(getattr(initial_args, "tutte_com_expand", 0.0)))
    tutte_auto_aspect_var = tk.BooleanVar(value=bool(getattr(initial_args, "tutte_auto_aspect", True)))
    tutte_auto_orient_var = tk.BooleanVar(value=bool(getattr(initial_args, "tutte_auto_orient", True)))
    tutte_orient_var = tk.StringVar(value=str(getattr(initial_args, "tutte_orient", 0.0)))
    show_tutte_outline_var = tk.BooleanVar(value=bool(getattr(initial_args, "show_tutte_outline", False)))
    show_tutte_pca_var = tk.BooleanVar(value=bool(getattr(initial_args, "show_tutte_pca", False)))
    hole_ratio_var = tk.StringVar(value=str(getattr(initial_args, "hole_ratio", 0.4)))
    hole_swap_var = tk.BooleanVar(value=bool(getattr(initial_args, "hole_swap", False)))
    invert_ring_var = tk.BooleanVar(value=bool(getattr(initial_args, "invert_ring", False)))
    ring_tilt_var = tk.StringVar(value=str(getattr(initial_args, "ring_tilt", 90.0)))
    min_sep_var = tk.StringVar(value=str(getattr(initial_args, "min_sep", 0.0)))
    ring_equalize_var = tk.StringVar(value=str(getattr(initial_args, "ring_equalize", 0.0)))
    relax_passes_var = tk.StringVar(value=str(getattr(initial_args, "relax_passes", 0)))
    relax_strength_var = tk.StringVar(value=str(getattr(initial_args, "relax_strength", 0.5)))
    surface_shape_var = tk.StringVar(value=getattr(initial_args, "surface_shape", "ellipsoid"))
    surface_auto_orient_var = tk.BooleanVar(value=bool(getattr(initial_args, "surface_auto_orient", True)))
    surface_auto_aspect_var = tk.BooleanVar(value=bool(getattr(initial_args, "surface_auto_aspect", True)))
    surface_aspect_var = tk.StringVar(value=str(getattr(initial_args, "surface_aspect", 1.6)))
    surface_tube_var = tk.StringVar(value=str(getattr(initial_args, "surface_tube", 0.35)))
    surface_orient_var = tk.StringVar(value=str(getattr(initial_args, "surface_orient", 0.0)))
    surface_tilt_var = tk.StringVar(value=str(getattr(initial_args, "surface_tilt", 0.0)))
    proj_elev_var = tk.StringVar(value=str(getattr(initial_args, "proj_elev", 25.0)))
    proj_azim_var = tk.StringVar(value=str(getattr(initial_args, "proj_azim", -60.0)))
    proj_roll_var = tk.StringVar(value=str(getattr(initial_args, "proj_roll", 0.0)))
    proj_lw_var = tk.StringVar(value=str(getattr(initial_args, "proj_line_width", 2.2)))
    proj_views_var = tk.StringVar(value=str(getattr(initial_args, "proj_views", "current,top")))
    proj_skeleton_var = tk.BooleanVar(value=bool(getattr(initial_args, "proj_skeleton", False)))
    proj_skeleton_ids_var = tk.BooleanVar(value=bool(getattr(initial_args, "proj_skeleton_ids", False)))
    proj_save_with_xyz_var = tk.BooleanVar(value=bool(getattr(initial_args, "proj_save_with_xyz", False)))
    proj_depth_fade_var = tk.StringVar(value=str(getattr(initial_args, "proj_depth_fade", 0.55)))
    proj_perspective_var = tk.StringVar(value=str(getattr(initial_args, "proj_perspective", 0.0)))
    proj_skeleton_only_var = tk.BooleanVar(value=bool(getattr(initial_args, "proj_skeleton_only", False)))
    proj_grid_density_var = tk.StringVar(value=str(getattr(initial_args, "proj_grid_density", 24)))
    proj_grid_color_var = tk.StringVar(value=str(getattr(initial_args, "proj_grid_color", "0.82")))
    proj_grid_lw_var = tk.StringVar(value=str(getattr(initial_args, "proj_grid_lw", 0.4)))
    proj_chords_var = tk.BooleanVar(value=not bool(getattr(initial_args, "proj_hide_chords", False)))
    proj_crossing_dots_var = tk.BooleanVar(value=bool(getattr(initial_args, "proj_crossing_dots", False)))
    proj_fixed_grid_var = tk.BooleanVar(value=bool(getattr(initial_args, "proj_fixed_grid", False)))

    # Registry of parameter widgets that are dynamically greyed out when they do
    # not apply to the current layout / sphere-layout / shape selection.  Each
    # entry maps a key to a list of (widget, kind) pairs.
    dynamic_widgets = {}

    def add_entry(label, var, width=12, help_key=None, key=None):
        nonlocal row
        lbl = ttk.Label(settings, text=label)
        lbl.grid(row=row, column=0, sticky="w", pady=3)
        ent = ttk.Entry(settings, textvariable=var, width=width)
        ent.grid(row=row, column=1, sticky="w", pady=3)
        if help_key:
            add_help_button(row, help_key)
        if key:
            dynamic_widgets.setdefault(key, []).extend([(lbl, "label"), (ent, "entry")])
        row += 1
        return lbl, ent

    ttk.Label(settings, text="negative even means").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Combobox(
        settings, textvariable=neg_var, values=["over", "under"],
        width=13, state="readonly"
    ).grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "negative_even")
    row += 1

    ttk.Label(settings, text="layout").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Combobox(
        settings, textvariable=layout_var,
        values=["tutte", "shaped-tutte", "holed-tutte", "sphere-stereo", "planar", "spring", "kamada"],
        width=15, state="readonly"
    ).grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "layout")
    row += 1

    _tutte_shape_lbl = ttk.Label(settings, text="tutte shape")
    _tutte_shape_lbl.grid(row=row, column=0, sticky="w", pady=3)
    _tutte_shape_combo = ttk.Combobox(
        settings, textvariable=tutte_shape_var,
        values=list(TUTTE_SHAPE_CHOICES), width=18, state="readonly"
    )
    _tutte_shape_combo.grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "tutte_shape")
    dynamic_widgets.setdefault("tutte_shape", []).extend(
        [(_tutte_shape_lbl, "label"), (_tutte_shape_combo, "combo")]
    )
    row += 1
    _tutte_auto_aspect_chk = ttk.Checkbutton(
        settings, text="tutte auto aspect (from PCA)", variable=tutte_auto_aspect_var
    )
    _tutte_auto_aspect_chk.grid(row=row, column=1, sticky="w", pady=2)
    # Live readout of the computed auto-aspect value (updated after each preview).
    auto_aspect_value_lbl = ttk.Label(settings, text="", foreground="gray40")
    auto_aspect_value_lbl.grid(row=row, column=2, sticky="w", pady=2)
    add_help_button(row, "tutte_auto_aspect")
    dynamic_widgets.setdefault("tutte_auto_aspect", []).append((_tutte_auto_aspect_chk, "check"))
    row += 1
    add_entry("tutte aspect", tutte_aspect_var, help_key="tutte_aspect", key="tutte_aspect")
    add_entry("tutte corner radius", tutte_corner_var,
              help_key="tutte_corner_radius", key="tutte_corner")
    add_entry("tutte decompress", tutte_decompress_var,
              help_key="tutte_decompress", key="tutte_decompress")
    add_entry("tutte COM expand", tutte_com_expand_var,
              help_key="tutte_com_expand", key="tutte_com_expand")
    _tutte_auto_orient_chk = ttk.Checkbutton(
        settings, text="tutte auto orient (PCA to x)", variable=tutte_auto_orient_var
    )
    _tutte_auto_orient_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "tutte_auto_orient")
    dynamic_widgets.setdefault("tutte_auto_orient", []).append((_tutte_auto_orient_chk, "check"))
    row += 1
    add_entry("tutte orient deg", tutte_orient_var, help_key="tutte_orient", key="tutte_orient")
    _tutte_outline_chk = ttk.Checkbutton(
        settings, text="show shape outline", variable=show_tutte_outline_var
    )
    _tutte_outline_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "show_tutte_outline")
    dynamic_widgets.setdefault("show_tutte_outline", []).append((_tutte_outline_chk, "check"))
    row += 1
    _tutte_pca_chk = ttk.Checkbutton(
        settings, text="show PCA / curved axis", variable=show_tutte_pca_var
    )
    _tutte_pca_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "show_tutte_pca")
    dynamic_widgets.setdefault("show_tutte_pca", []).append((_tutte_pca_chk, "check"))
    row += 1
    add_entry("hole ratio", hole_ratio_var, help_key="hole_ratio", key="hole_ratio")
    _hole_swap_chk = ttk.Checkbutton(
        settings, text="swap inner/outer face", variable=hole_swap_var
    )
    _hole_swap_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "hole_swap")
    dynamic_widgets.setdefault("hole_swap", []).append((_hole_swap_chk, "check"))
    row += 1
    _invert_ring_chk = ttk.Checkbutton(
        settings, text="invert ring (inside-out)", variable=invert_ring_var
    )
    _invert_ring_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "invert_ring")
    dynamic_widgets.setdefault("invert_ring", []).append((_invert_ring_chk, "check"))
    row += 1
    add_entry("ring tilt deg", ring_tilt_var, help_key="ring_tilt", key="ring_tilt")
    add_entry("min separation", min_sep_var, help_key="min_sep", key="min_sep")
    add_entry("ring equalize", ring_equalize_var,
              help_key="ring_equalize", key="ring_equalize")
    add_entry("relax passes", relax_passes_var, help_key="relax_passes")
    add_entry("relax strength", relax_strength_var, help_key="relax_strength")

    ttk.Label(settings, text="y direction").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Combobox(
        settings, textvariable=ydir_var,
        values=["top-to-bottom", "bottom-to-top"],
        width=18, state="readonly"
    ).grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "y_direction")
    row += 1

    add_entry("rotate degrees", rotate_var, help_key="rotate")
    add_entry("figure size", figsize_var, help_key="figsize")
    add_entry("DT label font", font_var, help_key="font_size")
    add_entry("crossing ID font", cid_font_var, help_key="crossing_id_font_size")
    add_entry("line width", lw_var, help_key="line_width")
    add_entry("gap fraction", gap_var, help_key="gap_frac")

    ttk.Separator(settings, orient="horizontal").grid(
        row=row, column=0, columnspan=4, sticky="ew", pady=(8, 6)
    )
    row += 1

    ttk.Label(settings, text="title").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Entry(settings, textvariable=title_var).grid(
        row=row, column=1, columnspan=2, sticky="ew", pady=3
    )
    add_help_button(row, "title")
    row += 1

    ttk.Checkbutton(settings, text="show crossing IDs", variable=show_ids_var).grid(
        row=row, column=1, columnspan=2, sticky="w", pady=2
    )
    add_help_button(row, "show_crossing_ids")
    row += 1
    ttk.Checkbutton(
        settings,
        text="color crossing IDs by over-strand",
        variable=color_ids_var,
    ).grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "color_crossing_ids_by_overstrand")
    row += 1
    ttk.Checkbutton(settings, text="hide DT labels", variable=hide_labels_var).grid(
        row=row, column=1, columnspan=2, sticky="w", pady=2
    )
    add_help_button(row, "hide_labels")
    row += 1
    ttk.Checkbutton(settings, text="no arrows", variable=no_arrows_var).grid(
        row=row, column=1, columnspan=2, sticky="w", pady=2
    )
    add_help_button(row, "no_arrows")
    row += 1

    ttk.Separator(settings, orient="horizontal").grid(
        row=row, column=0, columnspan=4, sticky="ew", pady=(8, 6)
    )
    row += 1

    ttk.Label(settings, text="Crossing order").grid(row=row, column=0, sticky="nw", pady=3)
    order_text = tk.Text(settings, height=4, width=54, wrap="word")
    order_text.grid(row=row, column=1, columnspan=2, sticky="ew", pady=3)
    if initial_args.crossing_order:
        order_text.insert("1.0", initial_args.crossing_order)
    add_help_button(row, "crossing_order")
    row += 1

    ttk.Label(settings, text="Explicit map").grid(row=row, column=0, sticky="nw", pady=3)
    map_text = tk.Text(settings, height=3, width=54, wrap="word")
    map_text.grid(row=row, column=1, columnspan=2, sticky="ew", pady=3)
    if initial_args.crossing_map:
        map_text.insert("1.0", initial_args.crossing_map)
    add_help_button(row, "crossing_map")
    row += 1

    # ------------------------------------------------------------------
    # Everything below goes on the second tab (3D XYZ / sphere parameters).
    # Switching 'settings' redirects the field helpers to that tab; reset the row.
    # ------------------------------------------------------------------
    settings = tab_3d
    row = 0
    ttk.Label(settings, text="sphere layout").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Combobox(
        settings, textvariable=sphere_layout_var,
        values=["spherical-kamada", "shaped-kamada", "stereographic", "stereo-safe"],
        width=18, state="readonly"
    ).grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "sphere_layout")
    row += 1

    _surface_shape_lbl = ttk.Label(settings, text="surface shape")
    _surface_shape_lbl.grid(row=row, column=0, sticky="w", pady=3)
    _surface_shape_combo = ttk.Combobox(
        settings, textvariable=surface_shape_var,
        values=list(SURFACE_SHAPES), width=18, state="readonly"
    )
    _surface_shape_combo.grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "surface_shape")
    dynamic_widgets.setdefault("surface_shape", []).extend(
        [(_surface_shape_lbl, "label"), (_surface_shape_combo, "combo")]
    )
    row += 1
    _surface_auto_aspect_chk = ttk.Checkbutton(
        settings, text="surface auto aspect (from PCA)", variable=surface_auto_aspect_var
    )
    _surface_auto_aspect_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "surface_auto_aspect")
    dynamic_widgets.setdefault("surface_auto_aspect", []).append((_surface_auto_aspect_chk, "check"))
    row += 1
    add_entry("surface aspect", surface_aspect_var,
              help_key="surface_aspect", key="surface_aspect")
    add_entry("surface tube ratio", surface_tube_var,
              help_key="surface_tube", key="surface_tube")
    _surface_auto_orient_chk = ttk.Checkbutton(
        settings, text="surface auto orient (PCA axis)", variable=surface_auto_orient_var
    )
    _surface_auto_orient_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "surface_auto_orient")
    dynamic_widgets.setdefault("surface_auto_orient", []).append((_surface_auto_orient_chk, "check"))
    row += 1
    add_entry("surface orient deg", surface_orient_var,
              help_key="surface_orient", key="surface_orient")
    add_entry("surface tilt deg", surface_tilt_var,
              help_key="surface_tilt", key="surface_tilt")

    _direct_conn_chk = ttk.Checkbutton(
        settings,
        text="direct connecting",
        variable=direct_connecting_var,
    )
    _direct_conn_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "direct_connecting")
    dynamic_widgets.setdefault("direct_connecting", []).append((_direct_conn_chk, "check"))
    row += 1
    add_entry("sphere radius", sphere_radius_var, help_key="sphere_radius")
    add_entry("sphere extent", sphere_extent_var, help_key="sphere_extent", key="sphere_extent")
    add_entry("crossing offset", sphere_offset_var, help_key="crossing_offset")
    add_entry("crossing angle deg", sphere_angle_var,
              help_key="sphere_crossing_angle", key="sphere_crossing_angle")
    add_entry("bump fraction", sphere_bump_var, help_key="sphere_bump_frac", key="sphere_bump")
    add_entry("XYZ point spacing", xyz_spacing_var, help_key="xyz_spacing")
    ttk.Checkbutton(
        settings,
        text="final smooth 3D curve",
        variable=xyz_final_smooth_var,
    ).grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "xyz_final_smooth")
    row += 1
    add_entry("smooth window", xyz_smooth_window_var, help_key="xyz_smooth_window")
    add_entry("smooth passes", xyz_smooth_passes_var, help_key="xyz_smooth_passes")
    add_entry("clearance (0=auto)", xyz_clearance_var,
              help_key="xyz_clearance", key="xyz_clearance")
    _xyz_repair_chk = ttk.Checkbutton(
        settings,
        text="repair 3D strand clearance",
        variable=xyz_repair_var,
    )
    _xyz_repair_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "xyz_repair")
    dynamic_widgets.setdefault("xyz_repair", []).append((_xyz_repair_chk, "check"))
    row += 1
    add_entry("XYZ decimals", xyz_decimals_var, help_key="xyz_decimals")
    ttk.Checkbutton(
        settings,
        text="repeat first point to close each component",
        variable=xyz_close_var,
    ).grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "xyz_close_components")
    row += 1

    # ------------------------------------------------------------------
    # V5.0 third tab: illustrative 2D projections of the 3D curve.
    # None of these variables triggers the 2D preview; the projection is only
    # recomputed by the 'Redraw 3D projection' button (the 3D build is slow).
    # ------------------------------------------------------------------
    settings = tab_view3d
    row = 0
    ttk.Label(
        settings,
        text=("Illustrative 2D projection of the 3D curve\n"
              "(uses the '3D XYZ' tab parameters; nearer strands\n"
              "are drawn over farther ones with white halos)"),
        foreground="gray25", justify="left",
    ).grid(row=row, column=0, columnspan=3, sticky="w", pady=(2, 6))
    row += 1
    add_entry("view elevation deg", proj_elev_var, help_key="proj_elev")
    add_entry("view azimuth deg", proj_azim_var, help_key="proj_azim")
    add_entry("view roll deg", proj_roll_var, help_key="proj_roll")
    add_entry("projection line width", proj_lw_var, help_key="proj_line_width")
    add_entry("depth fade (0-1)", proj_depth_fade_var, help_key="proj_depth_fade")
    add_entry("perspective (0=ortho)", proj_perspective_var, help_key="proj_perspective")
    _proj_skel_chk = ttk.Checkbutton(
        settings, text="show mapped skeleton", variable=proj_skeleton_var
    )
    _proj_skel_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "proj_skeleton")
    row += 1
    _proj_skel_ids_chk = ttk.Checkbutton(
        settings, text="label skeleton with crossing IDs",
        variable=proj_skeleton_ids_var,
    )
    _proj_skel_ids_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "proj_skeleton_ids")
    row += 1
    _proj_skel_only_chk = ttk.Checkbutton(
        settings, text="skeleton only (framework view)",
        variable=proj_skeleton_only_var,
    )
    _proj_skel_only_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "proj_skeleton_only")
    row += 1
    _proj_chords_chk = ttk.Checkbutton(
        settings, text="framework chords (crossing links)",
        variable=proj_chords_var,
    )
    _proj_chords_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "proj_chords")
    row += 1
    _proj_dots_chk = ttk.Checkbutton(
        settings, text="crossing anchor dots",
        variable=proj_crossing_dots_var,
    )
    _proj_dots_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "proj_crossing_dots")
    row += 1
    _proj_fixed_grid_chk = ttk.Checkbutton(
        settings, text="fixed grid while rotating",
        variable=proj_fixed_grid_var,
    )
    _proj_fixed_grid_chk.grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "proj_fixed_grid")
    row += 1
    add_entry("grid density (redraw)", proj_grid_density_var, help_key="proj_grid_density")
    add_entry("grid color", proj_grid_color_var, help_key="proj_grid_color")
    add_entry("grid line width", proj_grid_lw_var, help_key="proj_grid_lw")
    add_entry("views to save", proj_views_var, width=22, help_key="proj_views")
    ttk.Checkbutton(
        settings,
        text="'Save XYZ' also saves projections",
        variable=proj_save_with_xyz_var,
    ).grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "proj_save_with_xyz")
    row += 1
    ttk.Separator(settings, orient="horizontal").grid(
        row=row, column=0, columnspan=4, sticky="ew", pady=(8, 6)
    )
    row += 1
    _proj_redraw_btn = ttk.Button(
        settings, text="Redraw 3D projection",
        command=lambda: redraw_projection(),
    )
    _proj_redraw_btn.grid(row=row, column=0, columnspan=2, sticky="ew", pady=3)
    add_help_button(row, "redraw_projection")
    row += 1
    _proj_auto_btn = ttk.Button(
        settings, text="Auto view (avoid edge-on rings)",
        command=lambda: auto_view_gui(),
    )
    _proj_auto_btn.grid(row=row, column=0, columnspan=2, sticky="ew", pady=3)
    add_help_button(row, "proj_auto_view")
    row += 1
    ttk.Button(
        settings, text="Save projection(s)",
        command=lambda: save_projections_gui(),
    ).grid(row=row, column=0, columnspan=2, sticky="ew", pady=3)
    row += 1
    ttk.Button(
        settings, text="View XYZ (interactive 3D)",
        command=lambda: view_xyz(),
    ).grid(row=row, column=0, columnspan=2, sticky="ew", pady=3)
    row += 1
    ttk.Label(
        settings,
        text=("Quick-view buttons in the projection window\n"
              "re-project the cached curve instantly; press\n"
              "'Redraw 3D projection' after changing 3D XYZ\n"
              "parameters or the DT code."),
        foreground="gray40", justify="left",
    ).grid(row=row, column=0, columnspan=3, sticky="w", pady=(4, 2))
    row += 1

    latest = {
        "args": None,
        "state": None,
        "label_coords": None,
        "crossing_xy": None,
        "base_xlim": None,
        "base_ylim": None,
        "content_bounds": None,
    }
    preview_after = {"id": None}
    preview_zoom = {"value": 1.0}

    def _preview_aspect():
        """Width/height ratio of the live preview canvas."""
        try:
            w = max(float(canvas_widget.winfo_width()), 1.0)
            h = max(float(canvas_widget.winfo_height()), 1.0)
            return w / h
        except Exception:
            return 1.0

    def _apply_preview_zoom():
        """Preview-only zoom that reduces whitespace before it ever crops content."""
        try:
            zoom = float(preview_zoom.get("value", 1.0))
        except Exception:
            zoom = 1.0
        zoom = max(0.05, min(50.0, zoom))
        preview_zoom["value"] = zoom

        # V3.10: frame the preview through the shared routine so the saved image
        # (which calls the same routine) reproduces this view exactly.
        bounds = latest.get("content_bounds")
        if apply_content_framing(ax, bounds, aspect=_preview_aspect(), zoom=zoom):
            return

        base_xlim = latest.get("base_xlim") or tuple(ax.get_xlim())
        base_ylim = latest.get("base_ylim") or tuple(ax.get_ylim())
        try:
            x0, x1 = float(base_xlim[0]), float(base_xlim[1])
            y0, y1 = float(base_ylim[0]), float(base_ylim[1])
        except Exception:
            return
        cx = 0.5 * (x0 + x1)
        cy = 0.5 * (y0 + y1)
        hx = 0.5 * abs(x1 - x0) / zoom
        hy = 0.5 * abs(y1 - y0) / zoom
        ax.set_xlim(cx - hx, cx + hx)
        ax.set_ylim(cy - hy, cy + hy)

    def _rotate_preview_by(delta_degrees):
        try:
            current = float(rotate_var.get())
        except Exception:
            current = 0.0
        rotate_var.set("%.6g" % (current + float(delta_degrees)))
        schedule_preview()

    def _zoom_preview(multiplier):
        try:
            current = float(preview_zoom.get("value", 1.0))
        except Exception:
            current = 1.0
        preview_zoom["value"] = max(0.05, min(50.0, current * float(multiplier)))
        _apply_preview_zoom()
        canvas.draw_idle()

    def _reset_preview_zoom():
        preview_zoom["value"] = 1.0
        schedule_preview()

    def _float_value(var, name):
        try:
            return float(var.get())
        except Exception:
            raise ValueError("%s must be a number." % name)

    def _int_value(var, name):
        try:
            return int(var.get())
        except Exception:
            raise ValueError("%s must be an integer." % name)

    def _set_widget_enabled(widget, kind, enabled):
        try:
            if kind == "label":
                widget.configure(foreground=("" if enabled else "#9a9a9a"))
            elif kind == "combo":
                widget.configure(state=("readonly" if enabled else "disabled"))
            elif kind == "note":
                widget.configure(foreground=("gray40" if enabled else "#c8c8c8"))
            else:  # entry / check
                widget.configure(state=("normal" if enabled else "disabled"))
        except Exception:
            pass

    def _set_dynamic(key, enabled):
        for widget, kind in dynamic_widgets.get(key, ()):
            _set_widget_enabled(widget, kind, bool(enabled))

    def _apply_dynamic_states(*_a):
        """Grey out parameter fields that do not apply to the current choices."""
        layout = layout_var.get()
        is_shaped_tutte = (layout == "shaped-tutte")
        is_holed_tutte = (layout == "holed-tutte")
        # Shape / aspect / orient controls apply to both shaped and holed tutte.
        shaped_family = is_shaped_tutte or is_holed_tutte
        tutte_family = layout in ("tutte", "shaped-tutte")
        t_shape = tutte_shape_var.get()
        t_auto_aspect = bool(tutte_auto_aspect_var.get())
        t_auto_orient = bool(tutte_auto_orient_var.get())
        _set_dynamic("tutte_shape", shaped_family)
        _set_dynamic("tutte_auto_aspect", shaped_family)
        _set_dynamic("tutte_auto_orient", shaped_family)
        _set_dynamic(
            "tutte_aspect",
            shaped_family and (not t_auto_aspect)
            and t_shape in ("ellipse", "rectangle", "rounded-rectangle"),
        )
        _set_dynamic("tutte_corner", shaped_family and t_shape == "rounded-rectangle")
        # 'tutte orient deg' tilts the shape's aspect axis from the PCA axis, so it
        # applies to any non-circular shape regardless of the auto-orient framing.
        _set_dynamic("tutte_orient", shaped_family and t_shape != "circle")
        _set_dynamic("show_tutte_outline", shaped_family)
        _set_dynamic("show_tutte_pca", shaped_family)
        # holed-tutte only.
        _set_dynamic("hole_ratio", is_holed_tutte)
        _set_dynamic("hole_swap", is_holed_tutte)
        _set_dynamic("invert_ring", is_holed_tutte)
        _set_dynamic("ring_tilt", is_holed_tutte)
        _set_dynamic("ring_equalize", is_holed_tutte or layout == "sphere-stereo")
        _set_dynamic("tutte_decompress", tutte_family)
        _set_dynamic("tutte_com_expand", tutte_family)

        sl = sphere_layout_var.get()
        is_shaped_k = (sl == "shaped-kamada")
        kamada_family = sl in ("spherical-kamada", "shaped-kamada")
        s_shape = surface_shape_var.get()
        s_auto_aspect = bool(surface_auto_aspect_var.get())
        s_auto_orient = bool(surface_auto_orient_var.get())
        _set_dynamic("surface_shape", is_shaped_k)
        _set_dynamic("surface_auto_aspect", is_shaped_k)
        _set_dynamic("surface_auto_orient", is_shaped_k)
        _set_dynamic(
            "surface_aspect",
            is_shaped_k and (not s_auto_aspect) and s_shape in ("ellipsoid", "cylinder"),
        )
        _set_dynamic("surface_tube", is_shaped_k and (not s_auto_aspect) and s_shape == "torus")
        _set_dynamic("surface_orient", is_shaped_k and (not s_auto_orient))
        _set_dynamic("surface_tilt", is_shaped_k and (not s_auto_orient))
        # Sphere knobs that only apply to the kamada family / stereographic.
        _set_dynamic("direct_connecting", kamada_family)
        _set_dynamic("sphere_crossing_angle", kamada_family)
        _set_dynamic("sphere_bump", True)
        _set_dynamic("sphere_extent", sl in ("stereographic", "stereo-safe"))
        _set_dynamic("xyz_repair", sl != "stereo-safe")
        _set_dynamic(
            "xyz_clearance",
            sl != "stereo-safe" and bool(xyz_repair_var.get()),
        )

    def collect_args():
        dt = dt_text.get("1.0", "end").strip()
        if not dt:
            raise ValueError("Please enter a signed DT code.")
        return argparse.Namespace(
            dt=dt,
            output=output_var.get().strip() or "link_diagram.svg",
            table=None,
            negative_even=neg_var.get(),
            crossing_order=order_text.get("1.0", "end").strip() or None,
            crossing_map=map_text.get("1.0", "end").strip() or None,
            show_crossing_ids=bool(show_ids_var.get()),
            color_crossing_ids_by_overstrand=bool(color_ids_var.get()),
            y_direction=ydir_var.get(),
            rotate=_float_value(rotate_var, "rotate degrees"),
            dpi=int(getattr(initial_args, "dpi", 200) or 200),
            layout=layout_var.get(),
            tutte_shape=tutte_shape_var.get(),
            tutte_aspect=_float_value(tutte_aspect_var, "tutte aspect"),
            tutte_corner_radius=_float_value(tutte_corner_var, "tutte corner radius"),
            tutte_decompress=_float_value(tutte_decompress_var, "tutte decompress"),
            tutte_com_expand=_float_value(tutte_com_expand_var, "tutte COM expand"),
            tutte_auto_aspect=bool(tutte_auto_aspect_var.get()),
            tutte_auto_orient=bool(tutte_auto_orient_var.get()),
            tutte_orient=_float_value(tutte_orient_var, "tutte orient deg"),
            show_tutte_outline=bool(show_tutte_outline_var.get()),
            show_tutte_pca=bool(show_tutte_pca_var.get()),
            hole_ratio=_float_value(hole_ratio_var, "hole ratio"),
            hole_swap=bool(hole_swap_var.get()),
            invert_ring=bool(invert_ring_var.get()),
            ring_tilt=_float_value(ring_tilt_var, "ring tilt deg"),
            min_sep=_float_value(min_sep_var, "min separation"),
            ring_equalize=_float_value(ring_equalize_var, "ring equalize"),
            relax_passes=_int_value(relax_passes_var, "relax passes"),
            relax_strength=_float_value(relax_strength_var, "relax strength"),
            title=title_var.get().strip() or None,
            font_size=_float_value(font_var, "DT label font"),
            crossing_id_font_size=_float_value(cid_font_var, "crossing ID font"),
            line_width=_float_value(lw_var, "line width"),
            gap_frac=_float_value(gap_var, "gap fraction"),
            figsize=_float_value(figsize_var, "figure size"),
            hide_labels=bool(hide_labels_var.get()),
            no_arrows=bool(no_arrows_var.get()),
            xyz_output=xyz_var.get().strip() or "link_sphere.xyz",
            write_xyz=True,
            no_image=False,
            sphere_layout=sphere_layout_var.get(),
            direct_connecting=bool(direct_connecting_var.get()),
            sphere_radius=_float_value(sphere_radius_var, "sphere radius"),
            sphere_extent=_float_value(sphere_extent_var, "sphere extent"),
            sphere_offset=_float_value(sphere_offset_var, "crossing offset"),
            sphere_bump_frac=_float_value(sphere_bump_var, "bump fraction"),
            sphere_crossing_angle=_float_value(sphere_angle_var, "crossing angle"),
            xyz_spacing=_float_value(xyz_spacing_var, "XYZ point spacing"),
            xyz_final_smooth=bool(xyz_final_smooth_var.get()),
            xyz_smooth_window=_int_value(xyz_smooth_window_var, "smooth window"),
            xyz_smooth_passes=_int_value(xyz_smooth_passes_var, "smooth passes"),
            xyz_clearance=_float_value(xyz_clearance_var, "XYZ clearance"),
            no_xyz_repair=not bool(xyz_repair_var.get()),
            xyz_decimals=_int_value(xyz_decimals_var, "XYZ decimals"),
            xyz_close_components=bool(xyz_close_var.get()),
            surface_shape=surface_shape_var.get(),
            surface_auto_orient=bool(surface_auto_orient_var.get()),
            surface_auto_aspect=bool(surface_auto_aspect_var.get()),
            surface_aspect=_float_value(surface_aspect_var, "surface aspect"),
            surface_tube=_float_value(surface_tube_var, "surface tube ratio"),
            surface_orient=_float_value(surface_orient_var, "surface orient deg"),
            surface_tilt=_float_value(surface_tilt_var, "surface tilt deg"),
            proj_elev=_float_value(proj_elev_var, "view elevation deg"),
            proj_azim=_float_value(proj_azim_var, "view azimuth deg"),
            proj_roll=_float_value(proj_roll_var, "view roll deg"),
            proj_line_width=_float_value(proj_lw_var, "projection line width"),
            proj_views=proj_views_var.get().strip() or "current",
            proj_skeleton=bool(proj_skeleton_var.get()),
            proj_skeleton_ids=bool(proj_skeleton_ids_var.get()),
            proj_save_with_xyz=bool(proj_save_with_xyz_var.get()),
            proj_depth_fade=_float_value(proj_depth_fade_var, "depth fade"),
            proj_perspective=_float_value(proj_perspective_var, "perspective"),
            proj_skeleton_only=bool(proj_skeleton_only_var.get()),
            proj_grid_density=_int_value(proj_grid_density_var, "grid density"),
            proj_grid_color=proj_grid_color_var.get().strip() or "0.82",
            proj_grid_lw=_float_value(proj_grid_lw_var, "grid line width"),
            proj_hide_chords=not bool(proj_chords_var.get()),
            proj_crossing_dots=bool(proj_crossing_dots_var.get()),
            proj_fixed_grid=bool(proj_fixed_grid_var.get()),
        )

    def show_preview_error(message):
        ax.clear()
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "Preview error:\n%s" % message,
            ha="center",
            va="center",
            transform=ax.transAxes,
            wrap=True,
        )
        canvas.draw_idle()

    def update_preview():
        preview_after["id"] = None
        try:
            ns = collect_args()
            buf = io.StringIO()
            state = compute_diagram_objects(ns, status_stream=buf)
            label_coords, crossing_xy = render_state_on_axis(ax, state, ns)
            latest["args"] = ns
            latest["state"] = state
            latest["label_coords"] = label_coords
            latest["crossing_xy"] = crossing_xy
            latest["base_xlim"] = tuple(ax.get_xlim())
            latest["base_ylim"] = tuple(ax.get_ylim())
            latest["content_bounds"] = _axis_content_bounds(ax)
            _apply_preview_zoom()
            canvas.draw_idle()
            # Report the computed auto-aspect value next to the checkbox.
            av = state.get("aspect_value") if isinstance(state, dict) else None
            if bool(tutte_auto_aspect_var.get()) and av is not None \
                    and layout_var.get() in ("shaped-tutte", "holed-tutte"):
                auto_aspect_value_lbl.configure(text="= %.2f" % float(av))
            else:
                auto_aspect_value_lbl.configure(text="")
            set_log("[preview]\n" + buf.getvalue())
        except Exception as exc:
            latest["args"] = None
            latest["state"] = None
            latest["label_coords"] = None
            latest["crossing_xy"] = None
            latest["base_xlim"] = None
            latest["base_ylim"] = None
            latest["content_bounds"] = None
            msg = str(exc)
            show_preview_error(msg)
            set_log("[preview error] %s\n" % msg)

    def schedule_preview(_event=None):
        if preview_after["id"] is not None:
            try:
                root.after_cancel(preview_after["id"])
            except Exception:
                pass
        preview_after["id"] = root.after(450, update_preview)

    def _ask_output_path():
        # Always prompt for the directory/filename (like "Save as session").
        path = filedialog.asksaveasfilename(
            title="Save image",
            defaultextension=".svg",
            initialfile=os.path.basename(output_var.get().strip() or "link_diagram.svg"),
            filetypes=[
                ("SVG", "*.svg"), ("PDF", "*.pdf"), ("PNG", "*.png"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            raise ValueError("No output image path was selected.")
        output_var.set(path)
        return path

    def _ask_table_path():
        base = output_var.get().strip() or "link_diagram.svg"
        path = filedialog.asksaveasfilename(
            title="Save table (CSV)",
            defaultextension=".csv",
            initialfile=os.path.splitext(os.path.basename(base))[0] + ".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            raise ValueError("No table path was selected.")
        return path

    def _ask_xyz_path():
        path = filedialog.asksaveasfilename(
            title="Save XYZ",
            defaultextension=".xyz",
            initialfile=os.path.basename(xyz_var.get().strip() or "link_sphere.xyz"),
            filetypes=[("XYZ coordinate file", "*.xyz"), ("Text", "*.txt"),
                       ("All files", "*.*")],
        )
        if not path:
            raise ValueError("No XYZ coordinate path was selected.")
        xyz_var.set(path)
        return path

    def _render_table_coords(state, ns):
        tmp_fig, tmp_ax = plt.subplots(figsize=(float(ns.figsize), float(ns.figsize)))
        try:
            label_coords, crossing_xy = render_state_on_axis(tmp_ax, state, ns)
        finally:
            plt.close(tmp_fig)
        return label_coords, crossing_xy

    def _current_match_view():
        """Framing that reproduces the live preview when saving an image."""
        try:
            return {"aspect": _preview_aspect(),
                    "zoom": float(preview_zoom.get("value", 1.0))}
        except Exception:
            return {"aspect": 1.0, "zoom": 1.0}

    def save_image():
        try:
            ns = collect_args()
            ns.output = _ask_output_path()
            buf = io.StringIO()
            state = compute_diagram_objects(ns, status_stream=buf)
            draw(
                state["model"],
                state["P"],
                state["centers"],
                ns.output,
                dpi=ns.dpi,
                show_crossing_ids=ns.show_crossing_ids,
                color_crossing_ids_by_overstrand=ns.color_crossing_ids_by_overstrand,
                crossing_ids=state["crossing_ids"],
                title=ns.title,
                label_fontsize=ns.font_size,
                crossing_id_fontsize=ns.crossing_id_font_size,
                line_width=ns.line_width,
                gap_frac=ns.gap_frac,
                figsize=ns.figsize,
                show_labels=not ns.hide_labels,
                arrows=not ns.no_arrows,
                match_view=_current_match_view(),
                dt_code=ns.dt,
                layout=ns.layout,
                tutte_guides=state.get("tutte_guides"),
                show_tutte_outline=bool(getattr(ns, "show_tutte_outline", False)),
                show_tutte_pca=bool(getattr(ns, "show_tutte_pca", False)),
            )
            set_log(buf.getvalue() + "[ok] wrote %s\n" % ns.output)
            messagebox.showinfo("Saved", "Wrote image:\n%s" % ns.output)
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    def save_table():
        try:
            ns = collect_args()
            ns.table = _ask_table_path()
            buf = io.StringIO()
            state = compute_diagram_objects(ns, status_stream=buf)
            label_coords, crossing_xy = _render_table_coords(state, ns)
            write_table(
                state["model"],
                label_coords,
                crossing_xy,
                ns.table,
                crossing_ids=state["crossing_ids"],
            )
            set_log(buf.getvalue() + "[ok] wrote %s\n" % ns.table)
            messagebox.showinfo("Saved", "Wrote CSV table:\n%s" % ns.table)
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    def save_xyz():
        try:
            ns = collect_args()
            ns.xyz_output = _ask_xyz_path()
            buf = io.StringIO()
            state = compute_diagram_objects(ns, status_stream=buf)
            n_points = write_spherical_xyz(
                state["model"],
                state["P"],
                state["centers"],
                ns.xyz_output,
                G=state["graph"],
                sphere_radius=ns.sphere_radius,
                sphere_extent=ns.sphere_extent,
                sphere_offset=ns.sphere_offset,
                sphere_bump_frac=ns.sphere_bump_frac,
                xyz_spacing=ns.xyz_spacing,
                close_components=ns.xyz_close_components,
                decimals=ns.xyz_decimals,
                sphere_layout=ns.sphere_layout,
                sphere_crossing_angle=ns.sphere_crossing_angle,
                direct_connecting=ns.direct_connecting,
                xyz_final_smooth=ns.xyz_final_smooth,
                xyz_smooth_window=ns.xyz_smooth_window,
                xyz_smooth_passes=ns.xyz_smooth_passes,
                surface_shape=ns.surface_shape,
                surface_auto_orient=ns.surface_auto_orient,
                surface_auto_aspect=ns.surface_auto_aspect,
                surface_aspect=ns.surface_aspect,
                surface_tube=ns.surface_tube,
                surface_orient=ns.surface_orient,
                surface_tilt=ns.surface_tilt,
                xyz_repair=not getattr(ns, 'no_xyz_repair', False),
                xyz_clearance=(getattr(ns, 'xyz_clearance', 0.0) or None),
                audit_dt=getattr(ns, 'dt', None),
                status_stream=buf,
            )
            extra = ""
            if bool(proj_save_with_xyz_var.get()):
                try:
                    written = _save_projection_files(ns, ns.xyz_output)
                    extra = "".join("[ok] wrote %s\n" % p for p in written)
                except Exception as pexc:
                    extra = "[warn] projections failed: %s\n" % pexc
            set_log(buf.getvalue() + "[ok] wrote %s (%d xyz rows)\n" % (ns.xyz_output, n_points) + extra)
            messagebox.showinfo("Saved", "Wrote XYZ coordinates:\n%s%s" % (
                ns.xyz_output,
                ("\n\n+ %d projection file(s)" % extra.count("[ok]")) if extra.startswith("[ok]") else "",
            ))
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    def save_all():
        try:
            ns = collect_args()
            ns.output = _ask_output_path()
            ns.table = _ask_table_path()
            ns.xyz_output = _ask_xyz_path()
            buf = io.StringIO()
            state = compute_diagram_objects(ns, status_stream=buf)
            label_coords, crossing_xy = draw(
                state["model"],
                state["P"],
                state["centers"],
                ns.output,
                dpi=ns.dpi,
                show_crossing_ids=ns.show_crossing_ids,
                color_crossing_ids_by_overstrand=ns.color_crossing_ids_by_overstrand,
                crossing_ids=state["crossing_ids"],
                title=ns.title,
                label_fontsize=ns.font_size,
                crossing_id_fontsize=ns.crossing_id_font_size,
                line_width=ns.line_width,
                gap_frac=ns.gap_frac,
                figsize=ns.figsize,
                show_labels=not ns.hide_labels,
                arrows=not ns.no_arrows,
                match_view=_current_match_view(),
                dt_code=ns.dt,
                layout=ns.layout,
                tutte_guides=state.get("tutte_guides"),
                show_tutte_outline=bool(getattr(ns, "show_tutte_outline", False)),
                show_tutte_pca=bool(getattr(ns, "show_tutte_pca", False)),
            )
            write_table(
                state["model"],
                label_coords,
                crossing_xy,
                ns.table,
                crossing_ids=state["crossing_ids"],
            )
            n_points = write_spherical_xyz(
                state["model"],
                state["P"],
                state["centers"],
                ns.xyz_output,
                G=state["graph"],
                sphere_radius=ns.sphere_radius,
                sphere_extent=ns.sphere_extent,
                sphere_offset=ns.sphere_offset,
                sphere_bump_frac=ns.sphere_bump_frac,
                xyz_spacing=ns.xyz_spacing,
                close_components=ns.xyz_close_components,
                decimals=ns.xyz_decimals,
                sphere_layout=ns.sphere_layout,
                sphere_crossing_angle=ns.sphere_crossing_angle,
                direct_connecting=ns.direct_connecting,
                xyz_final_smooth=ns.xyz_final_smooth,
                xyz_smooth_window=ns.xyz_smooth_window,
                xyz_smooth_passes=ns.xyz_smooth_passes,
                surface_shape=ns.surface_shape,
                surface_auto_orient=ns.surface_auto_orient,
                surface_auto_aspect=ns.surface_auto_aspect,
                surface_aspect=ns.surface_aspect,
                surface_tube=ns.surface_tube,
                surface_orient=ns.surface_orient,
                surface_tilt=ns.surface_tilt,
                xyz_repair=not getattr(ns, 'no_xyz_repair', False),
                xyz_clearance=(getattr(ns, 'xyz_clearance', 0.0) or None),
                audit_dt=getattr(ns, 'dt', None),
                status_stream=buf,
            )
            set_log(
                buf.getvalue()
                + "[ok] wrote %s\n[ok] wrote %s\n[ok] wrote %s (%d xyz rows)\n"
                % (ns.output, ns.table, ns.xyz_output, n_points)
            )
            messagebox.showinfo(
                "Saved",
                "Wrote image:\n%s\n\nWrote CSV table:\n%s\n\nWrote XYZ coordinates:\n%s"
                % (ns.output, ns.table, ns.xyz_output),
            )
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    # ------------------------------------------------------------------
    # V5.0 projection window: persistent Toplevel showing the illustrative 2D
    # projection.  The heavy 3D construction runs only in redraw_projection();
    # the quick-view buttons re-project the cached curve instantly.
    # ------------------------------------------------------------------
    proj_cache = {
        "components": None,
        "skeleton": None,
        "grid_view": None,
        "audit_messages": [],
    }
    proj_win = {
        "win": None,
        "fig": None,
        "ax": None,
        "canvas": None,
        "label": None,
        "audit_label": None,
    }

    def _build_xyz_with_skeleton(ns, state, status_messages=None):
        skel = {}
        comps3 = build_spherical_xyz_components(
            state["model"], state["P"], state["centers"], G=state["graph"],
            sphere_radius=ns.sphere_radius,
            sphere_extent=ns.sphere_extent,
            sphere_offset=ns.sphere_offset,
            sphere_bump_frac=ns.sphere_bump_frac,
            xyz_spacing=ns.xyz_spacing,
            sphere_layout=ns.sphere_layout,
            sphere_crossing_angle=ns.sphere_crossing_angle,
            direct_connecting=ns.direct_connecting,
            xyz_final_smooth=ns.xyz_final_smooth,
            xyz_smooth_window=ns.xyz_smooth_window,
            xyz_smooth_passes=ns.xyz_smooth_passes,
            surface_shape=ns.surface_shape,
            surface_auto_orient=ns.surface_auto_orient,
            surface_auto_aspect=ns.surface_auto_aspect,
            surface_aspect=ns.surface_aspect,
            surface_tube=ns.surface_tube,
            surface_orient=ns.surface_orient,
            surface_tilt=ns.surface_tilt,
            skeleton_out=skel,
            skeleton_grid=getattr(ns, "proj_grid_density", 24),
            xyz_repair=not getattr(ns, "no_xyz_repair", False),
            xyz_clearance=(getattr(ns, "xyz_clearance", 0.0) or None),
            audit_dt=getattr(ns, "dt", None),
            status_messages=status_messages,
        )
        return comps3, skel

    def _set_projection_audit(messages):
        proj_cache["audit_messages"] = list(messages or ())
        label = proj_win.get("audit_label")
        if label is None:
            return
        audit_text, audit_level = _xyz_audit_display(messages)
        label.configure(
            text=audit_text,
            foreground=_audit_foreground(audit_level),
        )

    def _ensure_proj_window():
        if proj_win["win"] is not None:
            try:
                if proj_win["win"].winfo_exists():
                    return
            except Exception:
                pass
        win = tk.Toplevel(root)
        win.title("3D projection")
        try:
            apply_tk_window_icon(win)
        except Exception:
            pass
        win.geometry("720x760")
        frame = ttk.Frame(win, padding=4)
        frame.pack(fill="both", expand=True)
        frame.rowconfigure(1, weight=1)
        frame.columnconfigure(0, weight=1)
        bar = ttk.Frame(frame)
        bar.grid(row=0, column=0, sticky="ew", pady=(0, 4))
        for j in range(8):
            bar.columnconfigure(j, weight=1)

        def _bump(var, delta):
            try:
                var.set("%.6g" % (float(var.get()) + delta))
            except Exception:
                pass
            _render_proj_from_cache()

        def _preset(name):
            e, a, r = PROJECTION_VIEW_PRESETS[name]
            proj_elev_var.set("%.6g" % e)
            proj_azim_var.set("%.6g" % a)
            proj_roll_var.set("%.6g" % r)
            _render_proj_from_cache()

        ttk.Button(bar, text="Az -15", command=lambda: _bump(proj_azim_var, -15.0)).grid(row=0, column=0, sticky="ew", padx=1)
        ttk.Button(bar, text="Az +15", command=lambda: _bump(proj_azim_var, 15.0)).grid(row=0, column=1, sticky="ew", padx=1)
        ttk.Button(bar, text="El -15", command=lambda: _bump(proj_elev_var, -15.0)).grid(row=0, column=2, sticky="ew", padx=1)
        ttk.Button(bar, text="El +15", command=lambda: _bump(proj_elev_var, 15.0)).grid(row=0, column=3, sticky="ew", padx=1)
        ttk.Button(bar, text="Top", command=lambda: _preset("top")).grid(row=0, column=4, sticky="ew", padx=1)
        ttk.Button(bar, text="Front", command=lambda: _preset("front")).grid(row=0, column=5, sticky="ew", padx=1)
        ttk.Button(bar, text="Side", command=lambda: _preset("side")).grid(row=0, column=6, sticky="ew", padx=1)
        ttk.Button(bar, text="Auto view", command=lambda: auto_view_gui()).grid(row=0, column=7, sticky="ew", padx=1)
        ttk.Button(bar, text="Redraw 3D projection", command=lambda: redraw_projection()).grid(row=0, column=8, sticky="ew", padx=1)
        bar.columnconfigure(8, weight=1)

        pfig = Figure(figsize=(6.8, 6.8), dpi=100)
        pax = pfig.add_subplot(111)
        pax.axis("off")
        pcanvas = FigureCanvasTkAgg(pfig, master=frame)
        pcanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        audit_text, audit_level = _xyz_audit_display(
            proj_cache.get("audit_messages"))
        audit_info = tk.Label(
            frame,
            text=audit_text,
            anchor="w",
            justify="left",
            wraplength=690,
            padx=6,
            pady=4,
            relief="groove",
            foreground=_audit_foreground(audit_level),
        )
        audit_info.grid(row=2, column=0, sticky="ew", pady=(4, 2))
        info = ttk.Label(frame, text="(no projection yet)", foreground="gray40")
        info.grid(row=3, column=0, sticky="w")
        proj_win.update({"win": win, "fig": pfig, "ax": pax, "canvas": pcanvas,
                         "label": info, "audit_label": audit_info, "zoom": 1.0})

        # V5.1 Chimera-style mouse navigation on the projection canvas:
        # left-drag rotates (horizontal = azimuth, vertical = elevation) with a
        # live fast redraw, releasing does the full halo render; the scroll
        # wheel zooms about the view centre.
        drag = {"on": False, "mode": "orbit", "x": 0.0, "y": 0.0,
                "e": 0.0, "a": 0.0, "r": 0.0, "last": 0.0}

        def _canvas_center():
            try:
                w = pcanvas.get_tk_widget()
                return 0.5 * float(w.winfo_width()), 0.5 * float(w.winfo_height())
            except Exception:
                return 340.0, 340.0

        def _on_press(event):
            if event.button not in (1, 3):
                return
            drag["on"] = True
            # Right-drag (or Shift+left-drag) rolls about the view axis,
            # pivoting on the canvas centre; plain left-drag orbits.
            shift = "shift" in str(getattr(event, "key", "") or "").lower()
            drag["mode"] = "roll" if (event.button == 3 or shift) else "orbit"
            drag["x"], drag["y"] = event.x, event.y
            try:
                drag["e"] = float(proj_elev_var.get())
                drag["a"] = float(proj_azim_var.get())
                drag["r"] = float(proj_roll_var.get())
            except Exception:
                drag["e"], drag["a"], drag["r"] = 25.0, -60.0, 0.0

        def _on_motion(event):
            if not drag["on"]:
                return
            now = time.time()
            if now - drag["last"] < 0.05:      # ~20 fps throttle
                return
            drag["last"] = now
            if drag["mode"] == "roll":
                cx, cy = _canvas_center()
                a0 = math.atan2(float(drag["y"]) - cy, float(drag["x"]) - cx)
                a1 = math.atan2(float(event.y) - cy, float(event.x) - cx)
                # Tk pixel y grows downward, so the sign gives the natural
                # "grab and turn" feel.
                new_r = drag["r"] + math.degrees(a1 - a0)
                new_r = ((new_r + 180.0) % 360.0) - 180.0
                proj_roll_var.set("%.6g" % new_r)
                _render_proj_from_cache(fast=True)
                return
            dx = float(event.x - drag["x"])
            dy = float(event.y - drag["y"])
            new_a = drag["a"] - 0.4 * dx
            new_e = min(89.9, max(-89.9, drag["e"] + 0.4 * dy))
            new_a = ((new_a + 180.0) % 360.0) - 180.0
            proj_azim_var.set("%.6g" % new_a)
            proj_elev_var.set("%.6g" % new_e)
            _render_proj_from_cache(fast=True)

        def _on_release(event):
            if not drag["on"]:
                return
            drag["on"] = False
            _render_proj_from_cache(fast=False)

        def _on_scroll(event):
            f = 1.2 if getattr(event, "button", "") == "up" or getattr(event, "step", 0) > 0 else (1.0 / 1.2)
            proj_win["zoom"] = min(40.0, max(0.2, float(proj_win.get("zoom", 1.0)) * f))
            _render_proj_from_cache(fast=False)

        pcanvas.mpl_connect("button_press_event", _on_press)
        pcanvas.mpl_connect("motion_notify_event", _on_motion)
        pcanvas.mpl_connect("button_release_event", _on_release)
        pcanvas.mpl_connect("scroll_event", _on_scroll)

    def _render_proj_from_cache(fast=False):
        if proj_cache["components"] is None:
            return
        _ensure_proj_window()
        try:
            e = float(proj_elev_var.get()); a = float(proj_azim_var.get())
            r = float(proj_roll_var.get()); lw = float(proj_lw_var.get())
        except Exception:
            e, a, r, lw = 25.0, -60.0, 0.0, 2.2
        try:
            df = float(proj_depth_fade_var.get())
        except Exception:
            df = 0.55
        try:
            pp = float(proj_perspective_var.get())
        except Exception:
            pp = 0.0
        try:
            glw = float(proj_grid_lw_var.get())
        except Exception:
            glw = 0.4
        grid_basis = None
        if bool(proj_fixed_grid_var.get()):
            grid_view = proj_cache.get("grid_view") or (e, a, r)
            try:
                grid_basis = _projection_basis(*grid_view)
            except Exception:
                grid_basis = _projection_basis(e, a, r)
        render_projection_on_axis(
            proj_win["ax"], proj_cache["components"],
            elev=e, azim=a, roll=r, line_width=lw,
            skeleton=proj_cache["skeleton"],
            show_skeleton=bool(proj_skeleton_var.get()),
            skeleton_ids=bool(proj_skeleton_ids_var.get()),
            depth_fade=df, perspective=pp,
            skeleton_only=bool(proj_skeleton_only_var.get()),
            grid_color=proj_grid_color_var.get().strip() or "0.82",
            grid_lw=glw,
            show_chords=bool(proj_chords_var.get()),
            show_crossing_dots=bool(proj_crossing_dots_var.get()),
            grid_basis=grid_basis,
            fast=fast,
        )
        # Re-apply the wheel zoom about the view centre.
        try:
            zf = float(proj_win.get("zoom", 1.0))
            if zf != 1.0:
                x0, x1 = proj_win["ax"].get_xlim()
                y0, y1 = proj_win["ax"].get_ylim()
                cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
                hx, hy = 0.5 * (x1 - x0) / zf, 0.5 * (y1 - y0) / zf
                proj_win["ax"].set_xlim(cx - hx, cx + hx)
                proj_win["ax"].set_ylim(cy - hy, cy + hy)
        except Exception:
            pass
        proj_win["label"].configure(
            text="elev=%.4g  azim=%.4g  roll=%.4g   left-drag = orbit, "
                 "right/shift-drag = roll, wheel = zoom%s (cached curve; "
                 "'Redraw 3D projection' after 3D XYZ / DT changes)"
                 % (e, a, r, "; grid locked" if bool(proj_fixed_grid_var.get()) else "")
        )
        try:
            (proj_win["canvas"].draw if fast else proj_win["canvas"].draw_idle)()
        except Exception:
            pass

    def auto_view_gui():
        try:
            if proj_cache["components"] is None:
                redraw_projection()
                if proj_cache["components"] is None:
                    return
            e, a = auto_projection_view(proj_cache["components"])
            proj_elev_var.set("%.6g" % e)
            proj_azim_var.set("%.6g" % a)
            _render_proj_from_cache()
            set_log("[proj] auto view: elev=%.1f azim=%.1f\n" % (e, a))
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    def redraw_projection():
        try:
            ns = collect_args()
            buf = io.StringIO()
            xyz_messages = []
            state = compute_diagram_objects(ns, status_stream=buf)
            comps3, skel = _build_xyz_with_skeleton(ns, state, xyz_messages)
            proj_cache["components"] = comps3
            proj_cache["skeleton"] = skel
            proj_cache["grid_view"] = (ns.proj_elev, ns.proj_azim, ns.proj_roll)
            n_points = sum(len(arr) for arr in comps3)
            _ensure_proj_window()
            _set_projection_audit(xyz_messages)
            _render_proj_from_cache()
            set_log(buf.getvalue() + "".join(msg + "\n" for msg in xyz_messages)
                    + "[proj] rebuilt 3D curve (%d points, %d skeleton anchors).\n"
                    % (n_points, len(skel.get("points", []))))
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    def _save_projection_files(ns, base_path):
        if proj_cache["components"] is None:
            buf = io.StringIO()
            xyz_messages = []
            state = compute_diagram_objects(ns, status_stream=buf)
            comps3, skel = _build_xyz_with_skeleton(ns, state, xyz_messages)
            proj_cache["components"] = comps3
            proj_cache["skeleton"] = skel
            proj_cache["grid_view"] = (ns.proj_elev, ns.proj_azim, ns.proj_roll)
            proj_cache["audit_messages"] = xyz_messages
        views = _parse_projection_views(ns.proj_views)
        return save_projection_images(
            proj_cache["components"], base_path, views,
            elev=ns.proj_elev, azim=ns.proj_azim, roll=ns.proj_roll,
            line_width=ns.proj_line_width,
            skeleton=proj_cache["skeleton"],
            show_skeleton=bool(ns.proj_skeleton),
            skeleton_ids=bool(ns.proj_skeleton_ids),
            figsize=ns.figsize, dpi=ns.dpi,
            depth_fade=ns.proj_depth_fade,
            perspective=ns.proj_perspective,
            skeleton_only=bool(ns.proj_skeleton_only),
            grid_color=ns.proj_grid_color,
            grid_lw=ns.proj_grid_lw,
            show_chords=not bool(ns.proj_hide_chords),
            show_crossing_dots=bool(ns.proj_crossing_dots),
            fixed_grid=bool(ns.proj_fixed_grid),
            fixed_grid_view=proj_cache.get("grid_view")
            or (ns.proj_elev, ns.proj_azim, ns.proj_roll),
        )

    def save_projections_gui():
        try:
            ns = collect_args()
            base = filedialog.asksaveasfilename(
                title="Save projections (base name; _proj_<view>.svg/.png is appended)",
                defaultextension=".xyz",
                initialfile=os.path.basename(xyz_var.get().strip() or "link_sphere.xyz"),
                filetypes=[("Base name (extension ignored)", "*.*")],
            )
            if not base:
                return
            written = _save_projection_files(ns, base)
            set_log("".join("[ok] wrote %s\n" % p for p in written))
            messagebox.showinfo("Saved", "Wrote %d projection file(s):\n%s"
                                % (len(written), "\n".join(written)))
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    def view_xyz():
        try:
            ns = collect_args()
            buf = io.StringIO()
            xyz_messages = []
            state = compute_diagram_objects(ns, status_stream=buf)
            xyz_components = build_spherical_xyz_components(
                state["model"],
                state["P"],
                state["centers"],
                G=state["graph"],
                sphere_radius=ns.sphere_radius,
                sphere_extent=ns.sphere_extent,
                sphere_offset=ns.sphere_offset,
                sphere_bump_frac=ns.sphere_bump_frac,
                xyz_spacing=ns.xyz_spacing,
                sphere_layout=ns.sphere_layout,
                sphere_crossing_angle=ns.sphere_crossing_angle,
                direct_connecting=ns.direct_connecting,
                xyz_final_smooth=ns.xyz_final_smooth,
                xyz_smooth_window=ns.xyz_smooth_window,
                xyz_smooth_passes=ns.xyz_smooth_passes,
                surface_shape=ns.surface_shape,
                surface_auto_orient=ns.surface_auto_orient,
                surface_auto_aspect=ns.surface_auto_aspect,
                surface_aspect=ns.surface_aspect,
                surface_tube=ns.surface_tube,
                surface_orient=ns.surface_orient,
                surface_tilt=ns.surface_tilt,
                xyz_repair=not getattr(ns, "no_xyz_repair", False),
                xyz_clearance=(getattr(ns, "xyz_clearance", 0.0) or None),
                audit_dt=getattr(ns, "dt", None),
                status_messages=xyz_messages,
            )
            n_points = sum(len(arr) for arr in xyz_components)
            _open_xyz_viewer_window(
                root,
                xyz_components,
                closed=True,
                title="Sphere XYZ: %s, R=%g, offset=%g" % (
                    ns.sphere_layout, ns.sphere_radius, ns.sphere_offset
                ),
                audit_messages=xyz_messages,
            )
            set_log(
                buf.getvalue()
                + "".join(msg + "\n" for msg in xyz_messages)
                + "[view] generated %d XYZ points in memory.\n" % n_points
                + "[info] The viewer closes components visually; saved-file closure is controlled by the close-component checkbox.\n"
            )
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    # ------------------------------------------------------------------
    # Save / load a full session (every setting on the panel) as JSON.
    # ------------------------------------------------------------------
    session_string_vars = {
        "output": output_var, "xyz_output": xyz_var, "negative_even": neg_var,
        "layout": layout_var, "y_direction": ydir_var, "rotate": rotate_var,
        "font_size": font_var, "crossing_id_font_size": cid_font_var,
        "line_width": lw_var, "gap_frac": gap_var, "figsize": figsize_var,
        "title": title_var,
        "sphere_layout": sphere_layout_var, "sphere_radius": sphere_radius_var,
        "sphere_extent": sphere_extent_var, "sphere_offset": sphere_offset_var,
        "sphere_bump_frac": sphere_bump_var, "sphere_crossing_angle": sphere_angle_var,
        "xyz_spacing": xyz_spacing_var, "xyz_smooth_window": xyz_smooth_window_var,
        "xyz_smooth_passes": xyz_smooth_passes_var,
        "xyz_clearance": xyz_clearance_var, "xyz_decimals": xyz_decimals_var,
        "tutte_shape": tutte_shape_var, "tutte_aspect": tutte_aspect_var,
        "tutte_corner_radius": tutte_corner_var, "tutte_decompress": tutte_decompress_var,
        "tutte_com_expand": tutte_com_expand_var, "tutte_orient": tutte_orient_var,
        "hole_ratio": hole_ratio_var, "ring_tilt": ring_tilt_var,
        "min_sep": min_sep_var, "ring_equalize": ring_equalize_var,
        "relax_passes": relax_passes_var, "relax_strength": relax_strength_var,
        "proj_elev": proj_elev_var, "proj_azim": proj_azim_var,
        "proj_roll": proj_roll_var, "proj_line_width": proj_lw_var,
        "proj_views": proj_views_var, "proj_depth_fade": proj_depth_fade_var,
        "proj_perspective": proj_perspective_var,
        "proj_grid_density": proj_grid_density_var,
        "proj_grid_color": proj_grid_color_var,
        "proj_grid_lw": proj_grid_lw_var,
        "surface_shape": surface_shape_var, "surface_aspect": surface_aspect_var,
        "surface_tube": surface_tube_var, "surface_orient": surface_orient_var,
        "surface_tilt": surface_tilt_var,
    }
    session_bool_vars = {
        "direct_connecting": direct_connecting_var,
        "xyz_final_smooth": xyz_final_smooth_var,
        "xyz_repair": xyz_repair_var,
        "xyz_close_components": xyz_close_var,
        "show_crossing_ids": show_ids_var,
        "color_crossing_ids_by_overstrand": color_ids_var,
        "hide_labels": hide_labels_var, "no_arrows": no_arrows_var,
        "tutte_auto_aspect": tutte_auto_aspect_var,
        "tutte_auto_orient": tutte_auto_orient_var,
        "show_tutte_outline": show_tutte_outline_var,
        "show_tutte_pca": show_tutte_pca_var,
        "hole_swap": hole_swap_var,
        "invert_ring": invert_ring_var,
        "surface_auto_orient": surface_auto_orient_var,
        "surface_auto_aspect": surface_auto_aspect_var,
        "proj_skeleton": proj_skeleton_var,
        "proj_skeleton_ids": proj_skeleton_ids_var,
        "proj_save_with_xyz": proj_save_with_xyz_var,
        "proj_skeleton_only": proj_skeleton_only_var,
        "proj_chords": proj_chords_var,
        "proj_crossing_dots": proj_crossing_dots_var,
        "proj_fixed_grid": proj_fixed_grid_var,
    }
    session_text_widgets = {
        "dt": dt_text, "crossing_order": order_text, "crossing_map": map_text,
    }

    def save_session():
        try:
            path = filedialog.asksaveasfilename(
                title="Save session",
                defaultextension=".json",
                filetypes=[("JSON session", "*.json"), ("All files", "*.*")],
                initialfile="dt_session.json",
            )
            if not path:
                return
            data = {
                "script": "draw_dt_original_labels",
                "version": SCRIPT_VERSION,
                "strings": {k: v.get() for k, v in session_string_vars.items()},
                "bools": {k: bool(v.get()) for k, v in session_bool_vars.items()},
                "texts": {
                    k: w.get("1.0", "end").rstrip("\n")
                    for k, w in session_text_widgets.items()
                },
                "preview_zoom": float(preview_zoom.get("value", 1.0)),
            }
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2)
            set_log("[ok] saved session to %s\n" % path)
            messagebox.showinfo("Saved", "Wrote session:\n%s" % path)
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    def load_session():
        try:
            path = filedialog.askopenfilename(
                title="Load session",
                filetypes=[("JSON session", "*.json"), ("All files", "*.*")],
            )
            if not path:
                return
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            for k, val in (data.get("strings") or {}).items():
                if k in session_string_vars:
                    session_string_vars[k].set("" if val is None else str(val))
            for k, val in (data.get("bools") or {}).items():
                if k in session_bool_vars:
                    session_bool_vars[k].set(bool(val))
            for k, val in (data.get("texts") or {}).items():
                if k in session_text_widgets:
                    w = session_text_widgets[k]
                    w.delete("1.0", "end")
                    if val:
                        w.insert("1.0", str(val))
            try:
                preview_zoom["value"] = float(data.get("preview_zoom", 1.0))
            except Exception:
                preview_zoom["value"] = 1.0
            _apply_dynamic_states()
            update_preview()
            set_log("[ok] loaded session from %s\n" % path)
            messagebox.showinfo("Loaded", "Loaded session:\n%s" % path)
        except Exception as exc:
            msg = str(exc)
            set_log("[error] %s\n" % msg)
            messagebox.showerror("Error", msg)

    ttk.Button(save_bar, text="Save image", command=save_image).grid(
        row=0, column=0, sticky="ew", padx=2
    )
    ttk.Button(save_bar, text="Save table", command=save_table).grid(
        row=0, column=1, sticky="ew", padx=2
    )
    ttk.Button(save_bar, text="Save XYZ", command=save_xyz).grid(
        row=0, column=2, sticky="ew", padx=2
    )
    # V5.0: 'View XYZ' lives on the '3D view' tab now (with the projection
    # controls), so the save row stays compact.
    ttk.Button(save_bar, text="Save as session", command=save_session).grid(
        row=0, column=3, sticky="ew", padx=2
    )
    ttk.Button(save_bar, text="Load session", command=load_session).grid(
        row=0, column=4, sticky="ew", padx=2
    )
    def redraw_2d():
        # Clear the cached 3D-torus layout (and reroll its seed) so the holed-tutte
        # diagram is recomputed from scratch even without a parameter change.
        clear_holed_cache(reroll=True)
        update_preview()

    ttk.Button(save_bar, text="Redraw 2D", command=redraw_2d).grid(
        row=0, column=5, sticky="ew", padx=2
    )

    ttk.Button(preview_control_bar, text="Zoom +", command=lambda: _zoom_preview(1.25)).grid(
        row=0, column=0, sticky="ew", padx=2
    )
    ttk.Button(preview_control_bar, text="Zoom -", command=lambda: _zoom_preview(1.0 / 1.25)).grid(
        row=0, column=1, sticky="ew", padx=2
    )
    ttk.Button(preview_control_bar, text="Reset zoom", command=_reset_preview_zoom).grid(
        row=0, column=2, sticky="ew", padx=2
    )
    ttk.Button(preview_control_bar, text="Rotate -15", command=lambda: _rotate_preview_by(-15.0)).grid(
        row=0, column=3, sticky="ew", padx=2
    )
    ttk.Button(preview_control_bar, text="Rotate +15", command=lambda: _rotate_preview_by(15.0)).grid(
        row=0, column=4, sticky="ew", padx=2
    )
    tk.Button(
        preview_control_bar,
        text="?",
        width=2,
        bg="#cfeeff",
        activebackground="#aee3ff",
        command=lambda: show_arg_help("preview_zoom"),
    ).grid(row=0, column=5, sticky="ew", padx=2)
    # Quit lives on this second button row (no extra row just for Quit).
    ttk.Button(preview_control_bar, text="Quit", command=root.destroy).grid(
        row=0, column=6, sticky="ew", padx=2
    )

    # Only variables that actually change the 2D preview trigger an automatic
    # redraw.  Sphere XYZ / surface (3D-only) settings are deliberately excluded
    # so editing them stays smooth and never re-renders the 2D diagram; the 3D
    # output picks them up on Save/View XYZ, and "Refresh 2D" forces a redraw.
    watched_vars = [
        neg_var,
        layout_var,
        ydir_var,
        rotate_var,
        font_var,
        cid_font_var,
        lw_var,
        gap_var,
        figsize_var,
        title_var,
        show_ids_var,
        color_ids_var,
        hide_labels_var,
        no_arrows_var,
        tutte_shape_var,
        tutte_aspect_var,
        tutte_corner_var,
        tutte_decompress_var,
        tutte_com_expand_var,
        tutte_auto_aspect_var,
        tutte_auto_orient_var,
        tutte_orient_var,
        show_tutte_outline_var,
        show_tutte_pca_var,
        hole_ratio_var,
        hole_swap_var,
        invert_ring_var,
        ring_tilt_var,
        min_sep_var,
        ring_equalize_var,
        relax_passes_var,
        relax_strength_var,
    ]
    for var in watched_vars:
        var.trace_add("write", lambda *_args: schedule_preview())

    # Selector variables also drive the dynamic greying of non-relevant fields.
    for var in (
        layout_var,
        tutte_shape_var,
        tutte_auto_aspect_var,
        tutte_auto_orient_var,
        sphere_layout_var,
        surface_shape_var,
        surface_auto_orient_var,
        surface_auto_aspect_var,
        xyz_repair_var,
    ):
        var.trace_add("write", lambda *_a: _apply_dynamic_states())

    for text_widget in (dt_text, order_text, map_text):
        text_widget.bind("<KeyRelease>", schedule_preview)
        text_widget.bind("<FocusOut>", schedule_preview)

    def _collapse_dt_lines(_event=None):
        # Keep the DT code on a single line: strip any line breaks (e.g. from a
        # multi-line paste) and collapse whitespace runs.
        def _do():
            s = dt_text.get("1.0", "end-1c")
            if "\n" in s or "\r" in s:
                collapsed = " ".join(s.split())
                dt_text.delete("1.0", "end")
                dt_text.insert("1.0", collapsed)
        dt_text.after_idle(_do)

    dt_text.bind("<KeyRelease>", _collapse_dt_lines, add="+")
    dt_text.bind("<<Paste>>", _collapse_dt_lines, add="+")
    dt_text.bind("<FocusOut>", _collapse_dt_lines, add="+")
    _collapse_dt_lines()  # tidy any multi-line initial/loaded value

    def _proj_style_changed(*_a):
        # Render-side styling only: re-render the cached curve (no 3D rebuild).
        if proj_cache["components"] is not None:
            _render_proj_from_cache()
    for _v in (proj_chords_var, proj_crossing_dots_var, proj_fixed_grid_var,
               proj_skeleton_var, proj_skeleton_ids_var,
               proj_skeleton_only_var):
        _v.trace_add("write", _proj_style_changed)

    _apply_dynamic_states()
    root.after(200, update_preview)
    root.mainloop()
    return 0

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.gui or len(argv) == 0:
        return run_gui(args)

    if not args.dt:
        parser.error("--dt is required unless using --gui or no arguments for GUI mode.")

    return run_pipeline(args, status_stream=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
