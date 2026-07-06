#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_dt_original_labelsV3_12.py
===============================

Draw a smooth planar oriented link diagram from a signed Dowker-Thistlethwaite
(DT) code while preserving the original traversal labels supplied by the user.

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
  python draw_dt_original_labelsV3_12.py \
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
DT_LABEL_BOX_PAD = 0.22
CROSSING_ID_BOX_PAD = 0.28
matplotlib.rcParams["font.family"] = DIAGRAM_FONT_FAMILY
matplotlib.rcParams["font.sans-serif"] = [DIAGRAM_FONT_FAMILY]
matplotlib.rcParams["font.monospace"] = [DIAGRAM_FONT_FAMILY]
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.path import Path
from matplotlib.patches import PathPatch


EXAMPLE_DT = (
    "DT: [(-8,-12,16),(-24,-22,-28,-26),"
    "(-10,-14,-2),(-20,-6,-18,-4)]"
)


GUI_HELP_TEXT = {
    "dt": "Signed Dowker-Thistlethwaite code. Multi-component links use one tuple/list per component. Example:\n\nDT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]",
    "output": "Path for the 2D diagram image. The extension controls the format: .svg, .pdf, .png, etc.",
    "table": "Optional CSV table listing each crossing, original odd/even DT labels, components, over/under label, and 2D coordinates.",
    "xyz": "Path for the spherical x y z coordinate file. The file is plain three-column coordinates. Blank lines separate link components.",
    "negative_even": "DT sign convention. Default: a negative even DT label means the even-labeled visit is over, matching the common Sage/KnotTheory convention. Choose 'under' for the opposite convention.",
    "layout": "2D preview/image layout. tutte is usually smooth and clean; planar is safest; spring and kamada are aesthetic force-directed alternatives that are audited for false crossings.",
    "y_direction": "Controls the final 2D coordinate convention. top-to-bottom is the default drawing orientation; bottom-to-top flips it.",
    "rotate": "Rotates the final 2D scheme by this many degrees. This affects the 2D drawing and the stereographic sphere layout, but not the native spherical-kamada XYZ layout.",
    "dpi": "Raster resolution for PNG output. SVG/PDF are vector formats and are mostly unaffected.",
    "figsize": "Square Matplotlib figure size in inches for the saved 2D image.",
    "font_size": "Font size for the original DT traversal labels such as 1, -8, 16. Default: 7.",
    "crossing_id_font_size": "Font size for displayed crossing IDs such as c1, c7, c14. Default: 6.",
    "line_width": "2D strand line width in Matplotlib points. Default: 2.0.",
    "gap_frac": "Under-strand gap size in the 2D image. This is a ratio of the overall 2D diagram span, not an absolute coordinate distance. Default: 0.025.",
    "sphere_layout": "XYZ sphere layout. spherical-kamada distributes the graph directly over the sphere and is best for symmetric spherical models. stereographic maps the current 2D drawing onto a sphere.",
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


def tutte_layout_connected(G, emb):
    """Barycentric/Tutte embedding for one connected planar graph."""
    faces = planar_faces(emb)
    if not faces:
        return {n: np.array([0.0, 0.0]) for n in G.nodes()}

    outer = max(faces, key=len)
    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)

    fixed = {}
    m = len(outer)
    for j, n in enumerate(outer):
        th = 2.0 * np.pi * j / max(m, 1)
        fixed[n] = np.array([np.cos(th), np.sin(th)])

    A = np.zeros((N, N))
    bx = np.zeros(N)
    by = np.zeros(N)
    for n in nodes:
        i = idx[n]
        if n in fixed:
            A[i, i] = 1.0
            bx[i], by[i] = fixed[n]
        else:
            nbrs = list(G[n])
            A[i, i] = float(len(nbrs))
            for w in nbrs:
                A[i, idx[w]] -= 1.0

    X = np.linalg.solve(A, bx)
    Y = np.linalg.solve(A, by)
    return {n: np.array([X[idx[n]], Y[idx[n]]]) for n in nodes}


def compute_positions_connected(G, layout):
    ok, emb = nx.check_planarity(G)
    if not ok:
        raise RuntimeError(
            "The crossing graph is not planar; the DT code may be non-realizable."
        )

    if layout == "tutte":
        try:
            return tutte_layout_connected(G, emb)
        except np.linalg.LinAlgError:
            # Some connected planar graphs are singular for this simple
            # barycentric solve. NetworkX planar_layout is a safe fallback.
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
    raise ValueError("Unknown layout %r" % layout)


def compute_positions(G, layout):
    """Compute graph coordinates; disconnected graph pieces are packed side-by-side."""
    if G.number_of_nodes() == 0:
        return {}

    components = [list(nodes) for nodes in nx.connected_components(G)]
    if len(components) == 1:
        return compute_positions_connected(G, layout)

    packed = {}
    x_offset = 0.0
    gap = 0.75
    for nodes in components:
        H = G.subgraph(nodes).copy()
        P = compute_positions_connected(H, layout)
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


def _add_arrows(ax, dense, centers, color, scale, n_arrows=2):
    if len(dense) < 4 or not centers:
        return
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
            arrowprops=dict(arrowstyle="-|>", color=color, lw=0, mutation_scale=18),
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
            center_data = np.asarray(entry.get("center", artist.get_position()), float)
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
            cpt = dense[starts[cidx]]
            _plot_local_curve_piece(
                ax,
                dense,
                cpt,
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
            cpt = dense[starts[cidx]]
            _plot_local_curve_piece(
                ax,
                dense,
                cpt,
                radius=gap * 1.25,
                color=info["color"],
                lw=lw,
                zorder=4,
            )

    if arrows:
        for info in curve_infos:
            _add_arrows(ax, info["dense"], crossing_xy, info["color"], span * sc)

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
        for k in sorted(centers):
            cxy = crossing_xy[k]
            txt_color = "0.25"
            edge_color = "0.60"
            if color_crossing_ids_by_overstrand:
                over_pos = over_position_for_crossing(model, k)
                over_comp = model["comp_of"][over_pos]
                txt_color = color_of(over_comp)
                edge_color = txt_color
            artist = ax.text(
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
                bbox=dict(boxstyle="circle,pad=%s" % crossing_id_box_pad,
                          fc="white", ec=edge_color, alpha=0.78),
            )
            crossing_id_entries.append({"artist": artist, "center": cxy, "crossing_index": k})

    if show_labels and show_crossing_ids:
        _move_dt_labels_away_from_crossing_ids(ax, dt_label_entries, crossing_id_entries)
        for entry in dt_label_entries:
            label_coords[entry["p"]] = np.asarray(entry["artist"].get_position(), float)

    return label_coords, crossing_xy

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
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=11, fontfamily=DIAGRAM_FONT_FAMILY)
    _maximize_axis_in_figure(ax, has_title=bool(title))
    _tighten_axis_to_content(ax, pad_frac=0.08)
    _disable_figure_clipping(fig)
    ensure_parent_dir(out_path)
    if match_view:
        # Reproduce the preview framing exactly, then save the whole figure so
        # the file matches what the live preview shows.
        apply_content_framing(ax, _axis_content_bounds(ax), aspect=aspect, zoom=zoom)
        fig.savefig(out_path, dpi=dpi, pad_inches=0.0)
    else:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pad_inches=0.03)
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
):
    """Sphere-native method: distribute the diagram directly on S^2."""
    if G is None:
        G = build_gadget_graph(model)
    raw_dirs = _kamada_3d_unit_directions(G)
    P3, centers3 = _compact_crossing_gadgets_on_sphere(
        model, raw_dirs, crossing_angle=crossing_angle
    )

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
):
    """
    Return one Nx3 array per component for the spherical diagram.

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
    if layout == "stereographic":
        xyz_components = _build_stereo_xyz_components(
            model,
            P,
            centers,
            sphere_radius=sphere_radius,
            sphere_extent=sphere_extent,
            sphere_offset=sphere_offset,
            sphere_bump_frac=sphere_bump_frac,
            xyz_spacing=xyz_spacing,
        )
    elif layout == "spherical-kamada":
        xyz_components = _build_spherical_kamada_xyz_components(
            model,
            G,
            sphere_radius=sphere_radius,
            sphere_offset=sphere_offset,
            sphere_bump_frac=sphere_bump_frac,
            xyz_spacing=xyz_spacing,
            crossing_angle=sphere_crossing_angle,
            direct_connecting=direct_connecting,
        )
    else:
        raise ValueError("Unknown --sphere-layout %r." % sphere_layout)

    return _apply_final_xyz_smoothing(
        xyz_components,
        sphere_radius=sphere_radius,
        xyz_spacing=xyz_spacing,
        enabled=bool(xyz_final_smooth),
        smooth_window=xyz_smooth_window,
        smooth_passes=xyz_smooth_passes,
    )


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
    )
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


def _open_xyz_viewer_window(parent, xyz_components, closed=True, title="Sphere XYZ preview"):
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
    P = compute_positions(G, args.layout)
    P = transform_positions(P, args.y_direction, args.rotate)
    centers = crossing_centers(model, P)

    false = audit_false_crossings(model, P, centers)
    used_layout = args.layout
    if false > 0 and args.layout != "planar":
        out.write(
            "[warn] layout '%s' produced %d false crossing(s); "
            "falling back to 'planar'.\n" % (args.layout, false)
        )
        P = compute_positions(G, "planar")
        P = transform_positions(P, args.y_direction, args.rotate)
        centers = crossing_centers(model, P)
        false = audit_false_crossings(model, P, centers)
        used_layout = "planar"

    if false > 0:
        out.write(
            "[warn] %d false crossing(s) remain; inspect the drawing carefully.\n" % false
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
    }


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
        )
        out.write("[ok] wrote %s (%d xyz rows)\n" % (args.xyz_output, n_points))

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
        choices=["tutte", "planar", "spring", "kamada"],
        default="tutte",
        help=(
            "Layout engine. 'tutte' is default. If false crossings are detected, "
            "the script falls back to 'planar'."
        ),
    )
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
        choices=["spherical-kamada", "stereographic"],
        default="spherical-kamada",
        help=(
            "XYZ sphere layout. 'spherical-kamada' spreads the graph directly "
            "on S^2 and is the V3.9 default; 'stereographic' uses the 2D diagram "
            "and inverse stereographic projection."
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
    root.title("draw_dt_original_labelsV3_12")
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
    # Left side: live preview + save buttons + status log.
    # ------------------------------------------------------------------
    left = ttk.Frame(main)
    left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))
    left.columnconfigure(0, weight=1)
    left.rowconfigure(1, weight=1)

    ttk.Label(left, text="Live preview", font=("TkDefaultFont", 11, "bold")).grid(
        row=0, column=0, sticky="w"
    )

    preview_frame = ttk.Frame(left, relief="sunken", padding=2)
    preview_frame.grid(row=1, column=0, sticky="nsew", pady=(4, 6))
    preview_frame.columnconfigure(0, weight=1)
    preview_frame.rowconfigure(0, weight=1)

    fig = Figure(figsize=(6.4, 6.4), dpi=100)
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
    for j in range(6):
        preview_control_bar.columnconfigure(j, weight=1)

    log_text = tk.Text(left, height=8, state="disabled", wrap="word")
    log_text.grid(row=4, column=0, sticky="nsew")

    def set_log(text):
        log_text.configure(state="normal")
        log_text.delete("1.0", "end")
        log_text.insert("1.0", text)
        log_text.configure(state="disabled")

    # ------------------------------------------------------------------
    # Right side: scrollable parameter panel.
    # ------------------------------------------------------------------
    right_outer = ttk.Frame(main)
    right_outer.grid(row=0, column=1, sticky="nsew")
    right_outer.columnconfigure(0, weight=1)
    right_outer.rowconfigure(0, weight=1)

    settings_canvas = tk.Canvas(right_outer, borderwidth=0, highlightthickness=0)
    settings_scroll = ttk.Scrollbar(
        right_outer, orient="vertical", command=settings_canvas.yview
    )
    settings_canvas.configure(yscrollcommand=settings_scroll.set)
    settings_canvas.grid(row=0, column=0, sticky="nsew")
    settings_scroll.grid(row=0, column=1, sticky="ns")

    settings = ttk.Frame(settings_canvas, padding=(0, 0, 6, 0))
    settings_window = settings_canvas.create_window((0, 0), window=settings, anchor="nw")

    def _settings_configure(_event=None):
        settings_canvas.configure(scrollregion=settings_canvas.bbox("all"))

    def _canvas_configure(event):
        settings_canvas.itemconfigure(settings_window, width=event.width)

    settings.bind("<Configure>", _settings_configure)
    settings_canvas.bind("<Configure>", _canvas_configure)
    settings.columnconfigure(1, weight=1)
    settings.columnconfigure(2, weight=0)
    settings.columnconfigure(3, weight=0)

    ttk.Label(
        settings, text="Parameters", font=("TkDefaultFont", 11, "bold")
    ).grid(row=0, column=0, columnspan=4, sticky="w", pady=(0, 6))

    def show_arg_help(key):
        title = key.replace("_", " ").strip().title()
        message = GUI_HELP_TEXT.get(key, "No help text is available for this parameter.")
        messagebox.showinfo(title, message)

    def add_help_button(grid_row, key):
        btn = tk.Button(
            settings,
            text="?",
            width=2,
            bg="#cfeeff",
            activebackground="#aee3ff",
            relief="raised",
            command=lambda k=key: show_arg_help(k),
        )
        btn.grid(row=grid_row, column=3, sticky="w", padx=(5, 0), pady=2)
        return btn

    row = 1
    ttk.Label(settings, text="Signed DT code").grid(row=row, column=0, sticky="nw", pady=3)
    dt_text = tk.Text(settings, height=6, width=54, wrap="word")
    dt_text.grid(row=row, column=1, columnspan=2, sticky="ew", pady=3)
    dt_text.insert("1.0", initial_args.dt if initial_args.dt else EXAMPLE_DT)
    add_help_button(row, "dt")
    row += 1

    output_var = tk.StringVar(value=initial_args.output or "link_diagram.svg")
    table_var = tk.StringVar(value=initial_args.table or "")
    xyz_var = tk.StringVar(value=initial_args.xyz_output or "link_sphere.xyz")

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

    def browse_table():
        path = filedialog.asksaveasfilename(
            title="Save CSV table as",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
        )
        if path:
            table_var.set(path)

    def browse_xyz():
        path = filedialog.asksaveasfilename(
            title="Save spherical XYZ coordinates as",
            defaultextension=".xyz",
            filetypes=[("XYZ coordinate file", "*.xyz"), ("Text", "*.txt"), ("All files", "*.*")],
        )
        if path:
            xyz_var.set(path)

    ttk.Label(settings, text="Output image").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Entry(settings, textvariable=output_var).grid(row=row, column=1, sticky="ew", pady=3)
    ttk.Button(settings, text="Browse", command=browse_output).grid(row=row, column=2, sticky="ew", pady=3)
    add_help_button(row, "output")
    row += 1

    ttk.Label(settings, text="CSV table").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Entry(settings, textvariable=table_var).grid(row=row, column=1, sticky="ew", pady=3)
    ttk.Button(settings, text="Browse", command=browse_table).grid(row=row, column=2, sticky="ew", pady=3)
    add_help_button(row, "table")
    row += 1

    ttk.Label(settings, text="Sphere XYZ").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Entry(settings, textvariable=xyz_var).grid(row=row, column=1, sticky="ew", pady=3)
    ttk.Button(settings, text="Browse", command=browse_xyz).grid(row=row, column=2, sticky="ew", pady=3)
    add_help_button(row, "xyz")
    row += 1

    neg_var = tk.StringVar(value=initial_args.negative_even)
    layout_var = tk.StringVar(value=initial_args.layout)
    ydir_var = tk.StringVar(value=initial_args.y_direction)
    rotate_var = tk.StringVar(value=str(initial_args.rotate))
    dpi_var = tk.StringVar(value=str(initial_args.dpi))
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
    xyz_decimals_var = tk.StringVar(value=str(initial_args.xyz_decimals))
    xyz_close_var = tk.BooleanVar(value=bool(initial_args.xyz_close_components))

    show_ids_var = tk.BooleanVar(value=bool(initial_args.show_crossing_ids))
    color_ids_var = tk.BooleanVar(
        value=bool(initial_args.color_crossing_ids_by_overstrand)
    )
    hide_labels_var = tk.BooleanVar(value=bool(initial_args.hide_labels))
    no_arrows_var = tk.BooleanVar(value=bool(initial_args.no_arrows))

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
        values=["tutte", "planar", "spring", "kamada"],
        width=13, state="readonly"
    ).grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "layout")
    row += 1

    ttk.Label(settings, text="y direction").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Combobox(
        settings, textvariable=ydir_var,
        values=["top-to-bottom", "bottom-to-top"],
        width=18, state="readonly"
    ).grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "y_direction")
    row += 1

    def add_entry(label, var, width=12, help_key=None):
        nonlocal row
        ttk.Label(settings, text=label).grid(row=row, column=0, sticky="w", pady=3)
        ttk.Entry(settings, textvariable=var, width=width).grid(
            row=row, column=1, sticky="w", pady=3
        )
        if help_key:
            add_help_button(row, help_key)
        row += 1

    add_entry("rotate degrees", rotate_var, help_key="rotate")
    add_entry("DPI", dpi_var, help_key="dpi")
    add_entry("figure size", figsize_var, help_key="figsize")
    add_entry("DT label font", font_var, help_key="font_size")
    add_entry("crossing ID font", cid_font_var, help_key="crossing_id_font_size")
    add_entry("line width", lw_var, help_key="line_width")
    add_entry("gap fraction", gap_var, help_key="gap_frac")

    ttk.Separator(settings, orient="horizontal").grid(
        row=row, column=0, columnspan=4, sticky="ew", pady=(8, 6)
    )
    row += 1
    ttk.Label(settings, text="Sphere XYZ", font=("TkDefaultFont", 10, "bold")).grid(
        row=row, column=0, columnspan=4, sticky="w", pady=(0, 4)
    )
    row += 1
    ttk.Label(settings, text="sphere layout").grid(row=row, column=0, sticky="w", pady=3)
    ttk.Combobox(
        settings, textvariable=sphere_layout_var,
        values=["spherical-kamada", "stereographic"],
        width=18, state="readonly"
    ).grid(row=row, column=1, sticky="w", pady=3)
    add_help_button(row, "sphere_layout")
    row += 1
    ttk.Checkbutton(
        settings,
        text="direct connecting",
        variable=direct_connecting_var,
    ).grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "direct_connecting")
    row += 1
    add_entry("sphere radius", sphere_radius_var, help_key="sphere_radius")
    add_entry("sphere extent", sphere_extent_var, help_key="sphere_extent")
    add_entry("crossing offset", sphere_offset_var, help_key="crossing_offset")
    add_entry("crossing angle deg", sphere_angle_var, help_key="sphere_crossing_angle")
    add_entry("bump fraction", sphere_bump_var, help_key="sphere_bump_frac")
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
    add_entry("XYZ decimals", xyz_decimals_var, help_key="xyz_decimals")
    ttk.Checkbutton(
        settings,
        text="repeat first point to close each component",
        variable=xyz_close_var,
    ).grid(row=row, column=1, columnspan=2, sticky="w", pady=2)
    add_help_button(row, "xyz_close_components")
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
    ttk.Label(
        settings,
        text="Optional: crossing IDs ordered by odd labels 1,3,5,...",
        foreground="gray40",
    ).grid(row=row, column=1, columnspan=2, sticky="w")
    row += 1

    ttk.Label(settings, text="Explicit map").grid(row=row, column=0, sticky="nw", pady=3)
    map_text = tk.Text(settings, height=3, width=54, wrap="word")
    map_text.grid(row=row, column=1, columnspan=2, sticky="ew", pady=3)
    if initial_args.crossing_map:
        map_text.insert("1.0", initial_args.crossing_map)
    add_help_button(row, "crossing_map")
    row += 1
    ttk.Label(
        settings,
        text="Alternative: c1=1, c7=3, c14=5, ...",
        foreground="gray40",
    ).grid(row=row, column=1, columnspan=2, sticky="w")
    row += 1

    ttk.Button(settings, text="Refresh preview", command=lambda: update_preview()).grid(
        row=row, column=1, sticky="ew", pady=(10, 4)
    )

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

    def collect_args():
        dt = dt_text.get("1.0", "end").strip()
        if not dt:
            raise ValueError("Please enter a signed DT code.")
        return argparse.Namespace(
            dt=dt,
            output=output_var.get().strip() or "link_diagram.svg",
            table=table_var.get().strip() or None,
            negative_even=neg_var.get(),
            crossing_order=order_text.get("1.0", "end").strip() or None,
            crossing_map=map_text.get("1.0", "end").strip() or None,
            show_crossing_ids=bool(show_ids_var.get()),
            color_crossing_ids_by_overstrand=bool(color_ids_var.get()),
            y_direction=ydir_var.get(),
            rotate=_float_value(rotate_var, "rotate degrees"),
            dpi=_int_value(dpi_var, "DPI"),
            layout=layout_var.get(),
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
            xyz_decimals=_int_value(xyz_decimals_var, "XYZ decimals"),
            xyz_close_components=bool(xyz_close_var.get()),
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
        path = output_var.get().strip()
        if not path:
            browse_output()
            path = output_var.get().strip()
        if not path:
            raise ValueError("No output image path was selected.")
        return path

    def _ask_table_path():
        path = table_var.get().strip()
        if not path:
            browse_table()
            path = table_var.get().strip()
        if not path:
            raise ValueError("No CSV table path was selected.")
        return path

    def _ask_xyz_path():
        path = xyz_var.get().strip()
        if not path:
            browse_xyz()
            path = xyz_var.get().strip()
        if not path:
            raise ValueError("No XYZ coordinate path was selected.")
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
            )
            set_log(buf.getvalue() + "[ok] wrote %s (%d xyz rows)\n" % (ns.xyz_output, n_points))
            messagebox.showinfo("Saved", "Wrote XYZ coordinates:\n%s" % ns.xyz_output)
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

    def view_xyz():
        try:
            ns = collect_args()
            buf = io.StringIO()
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
            )
            n_points = sum(len(arr) for arr in xyz_components)
            _open_xyz_viewer_window(
                root,
                xyz_components,
                closed=True,
                title="Sphere XYZ: %s, R=%g, offset=%g" % (
                    ns.sphere_layout, ns.sphere_radius, ns.sphere_offset
                ),
            )
            set_log(
                buf.getvalue()
                + "[view] generated %d XYZ points in memory.\n" % n_points
                + "[info] The viewer closes components visually; saved-file closure is controlled by the close-component checkbox.\n"
            )
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
    ttk.Button(save_bar, text="View XYZ", command=view_xyz).grid(
        row=0, column=3, sticky="ew", padx=2
    )
    ttk.Button(save_bar, text="Save all", command=save_all).grid(
        row=0, column=4, sticky="ew", padx=2
    )
    ttk.Button(save_bar, text="Refresh 2D", command=update_preview).grid(
        row=0, column=5, sticky="ew", padx=2
    )
    ttk.Button(save_bar, text="Quit", command=root.destroy).grid(
        row=0, column=6, sticky="ew", padx=2
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

    watched_vars = [
        output_var,
        table_var,
        xyz_var,
        neg_var,
        layout_var,
        ydir_var,
        rotate_var,
        dpi_var,
        font_var,
        cid_font_var,
        lw_var,
        gap_var,
        figsize_var,
        title_var,
        sphere_layout_var,
        direct_connecting_var,
        sphere_radius_var,
        sphere_extent_var,
        sphere_offset_var,
        sphere_bump_var,
        sphere_angle_var,
        xyz_spacing_var,
        xyz_final_smooth_var,
        xyz_smooth_window_var,
        xyz_smooth_passes_var,
        xyz_decimals_var,
        xyz_close_var,
        show_ids_var,
        color_ids_var,
        hide_labels_var,
        no_arrows_var,
    ]
    for var in watched_vars:
        var.trace_add("write", lambda *_args: schedule_preview())

    for text_widget in (dt_text, order_text, map_text):
        text_widget.bind("<KeyRelease>", schedule_preview)
        text_widget.bind("<FocusOut>", schedule_preview)

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
