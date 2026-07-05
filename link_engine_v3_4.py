#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
link_engine_v3_4.py  (V3.4)
======================

Link-diagram engine for strand-passage experiments.

What is new in V3.4
-------------------
* Imports ``draw_dt_original_labelsV3_11.py``, whose standalone GUI can use the
  optional project icon asset without requiring it at runtime.

Design
------
* The continuity engine performs strand passages (crossing changes) and safe
  Reidemeister I/II reductions while preserving original component identity by
  construction -- component index i always carries original colour
  ``component_colors[i]``.
* When SnapPy/Spherogram is available, the engine asks SnapPy for a global
  simplified diagram after each passage, then recovers component colours by
  enumerating all component permutations and orientation reversals against the
  linking matrix, reinforced by per-component Jones signatures.

SnapPy is thus a diagram/invariant oracle, never an unquestioned carrier of
component labels.  A unique invariant match necessarily agrees with the physical
strand-continuity labelling, so the two paths never disagree.

What is new in V3.2
-------------------
* Drawing backend is now ``draw_dt_original_labelsV3_11.py``.
* ``render`` follows the drawing helper's own default 2-D layout pipeline
  (default layout, top-to-bottom orientation, and the false-crossing audit that
  falls back to a planar layout) so 2-D links look exactly as the standalone
  helper draws them.
* Every diagram carries a ``_dt_labels_valid`` flag.  Any diagram freshly built
  from a DT code (the original, the direct after-passage DT, and the
  SnapPy-simplified DT) draws its own original DT traversal labels, so labels
  stay visible even after a SnapPy simplification.  The flag is cleared the
  moment a diagram is passage-changed or Reidemeister-reduced.

What is new in V3.3
-------------------
* Optional backtrack-assisted SnapPy simplification.  ``simplify('global')`` is
  a greedy heuristic that can get stuck at a non-minimal diagram (e.g. reporting
  12 crossings for a link whose minimum is 10).  ``backtrack_simplify`` can now
  repeatedly *complicate then re-simplify* to escape such plateaus, keeping the
  fewest-crossing diagram found.  It is controlled by ``backtrack_rounds`` and
  ``backtrack_steps`` and defaults to OFF (rounds=0), so behaviour is unchanged
  unless the caller enables it.

Drawing backend
---------------
V3.4 imports ``draw_dt_original_labelsV3_11.py`` for the 2-D DT
parser/model/layout/render helpers.
"""

from __future__ import annotations

import copy
import itertools
import os
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import networkx as nx

# Keep pyplot import from selecting a fragile GUI backend under Sage.
os.environ.setdefault("MPLBACKEND", "Agg")
import draw_dt_original_labelsV3_11 as D

VERSION = "3.4"
DRAWING_MODULE_NAME = getattr(D, "__name__", "draw_dt_original_labelsV3_11")

# Defaults for the backtrack-assisted SnapPy simplification (ON by default in
# the GUI and --nongui as of V3.3).
DEFAULT_BACKTRACK_ROUNDS = 200
DEFAULT_BACKTRACK_STEPS = 30

# 2-D drawing defaults: mirror the standalone helper so V3.2 draws links exactly
# as ``draw_dt_original_labelsV3_11.py`` does by default.  These fall back to the
# helper's own module-level constants when present.
DEFAULT_LAYOUT = getattr(D, "DEFAULT_LAYOUT", "tutte")
DEFAULT_Y_DIRECTION = getattr(D, "DEFAULT_Y_DIRECTION", "top-to-bottom")
DEFAULT_ROTATE = float(getattr(D, "DEFAULT_ROTATE", 0.0))


def dt_to_string(dt_code: Sequence[Sequence[int]]) -> str:
    """Return SnapPy-style ``DT: [...]`` text for a multi-component DT code."""
    pieces = []
    for comp in dt_code:
        pieces.append(repr(tuple(int(x) for x in comp)))
    return "DT: [" + ", ".join(pieces) + "]"


def parse_dt_any(value: Any) -> List[Tuple[int, ...]]:
    """Parse DT text or normalize an already materialized DT object.

    ``draw_dt_original_labelsV3_11.parse_dt`` is the parser for user text.  The
    SnapPy global branch can also receive list/tuple-like objects from
    ``DT_code()``, so this wrapper normalizes those objects without importing any
    other drawing module.
    """
    if isinstance(value, str):
        return [tuple(int(x) for x in comp) for comp in D.parse_dt(value)]
    if value is None:
        return []
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, tuple):
        value = list(value)
    if not isinstance(value, list):
        return [tuple(int(x) for x in comp) for comp in D.parse_dt(str(value))]
    if all(isinstance(x, (int, np.integer)) for x in value):
        return [tuple(int(x) for x in value)]
    out: List[Tuple[int, ...]] = []
    for comp in value:
        if hasattr(comp, "tolist"):
            comp = comp.tolist()
        if isinstance(comp, tuple):
            comp = list(comp)
        if not isinstance(comp, list):
            return [tuple(int(x) for x in c) for c in D.parse_dt(str(value))]
        out.append(tuple(int(x) for x in comp))
    return out


# --------------------------------------------------------------------------- #
#  Small data records for SnapPy/component matching
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MatchCandidate:
    """One possible matching of SnapPy-simplified components to source ones.

    ``perm[simplified_component] = source_component``.  ``signs`` records a
    possible orientation reversal pattern in the simplified component order.
    """

    perm: Tuple[int, ...]
    signs: Tuple[int, ...]


@dataclass
class ComponentMatchResult:
    status: str
    message: str
    mapping: Optional[Dict[int, int]] = None
    candidates: List[MatchCandidate] = field(default_factory=list)
    unique_permutations: List[Tuple[int, ...]] = field(default_factory=list)
    used_signatures: bool = False

    @property
    def is_unique(self) -> bool:
        return self.status == "unique" and self.mapping is not None


# --------------------------------------------------------------------------- #
#  Gauss-system diagram
# --------------------------------------------------------------------------- #
class Diagram(object):
    """
    comps      : list of components; each component is a list of visit-ids (vids)
                 in oriented cyclic order.
    visit_x    : vid -> crossing id
    visit_over : vid -> bool (True if this strand is the OVER strand here)

    Each crossing id occurs in exactly two visits.  Component identity is carried
    by ``component_colors``: component index i in this *current* diagram should
    be drawn with original component colour ``component_colors[i]``.
    """

    def __init__(self, comps, visit_x, visit_over, next_vid=0, component_colors=None):
        self.comps = comps
        self.visit_x = visit_x
        self.visit_over = visit_over
        self._next_vid = next_vid
        self._touched = False
        if component_colors is None:
            component_colors = list(range(len(comps)))
        if len(component_colors) != len(comps):
            raise ValueError("component_colors must have one entry per component")
        self.component_colors = list(component_colors)
        self.last_snap_result: Optional[Dict[str, Any]] = None
        self.display_source = "continuity"
        # True while component index still maps to a known ORIGINAL component.
        # Once a SnapPy step cannot disambiguate colours, this becomes False and
        # stays False for all descendants.
        self.colours_tracked = True
        # True when this drawing is a pristine picture of ``_dt_model`` (i.e. it
        # was just built from a DT code and not yet passage-changed / reduced).
        # While True, ``render`` draws the original DT traversal labels.  A
        # crossing change or a Reidemeister reduction clears it.
        self._dt_labels_valid = False

    # ---- construction from the drawer's DT model ----
    @classmethod
    def from_dt_model(cls, model, component_colors=None):
        comps = []
        visit_x, visit_over = {}, {}
        vid = 0
        for cp in model["comp_positions"]:
            comp = []
            for p in cp:
                visit_x[vid] = model["pos_cross"][p]
                visit_over[vid] = bool(model["over_at"][p])
                comp.append(vid)
                vid += 1
            comps.append(comp)
        return cls(comps, visit_x, visit_over, next_vid=vid,
                   component_colors=component_colors)

    @classmethod
    def from_dt_code(cls, comps_dt, negative_even="over", component_colors=None):
        comps_dt = parse_dt_any(comps_dt)
        model = D.build_model(comps_dt, negative_even=negative_even)
        d = cls.from_dt_model(model, component_colors=component_colors)
        d._dt_model = model      # keep original labels for display
        d._dt_labels_valid = True  # pristine DT drawing: show its own labels
        return d

    def copy(self):
        d = Diagram(copy.deepcopy(self.comps), dict(self.visit_x),
                    dict(self.visit_over), self._next_vid,
                    component_colors=list(self.component_colors))
        d._touched = self._touched
        d.display_source = self.display_source
        d.colours_tracked = self.colours_tracked
        d.last_snap_result = self.last_snap_result
        d._dt_labels_valid = self._dt_labels_valid
        if hasattr(self, "_dt_model"):
            d._dt_model = self._dt_model
        return d

    # ---- basic queries ----
    def crossings(self):
        return sorted(set(self.visit_x.values()))

    def visits_of(self, x):
        return [v for v in self.visit_x if self.visit_x[v] == x]

    def num_components(self):
        return len(self.comps)

    def num_visits(self):
        return sum(len(c) for c in self.comps)

    # ---- strand passage (crossing change) ----
    def crossing_change(self, x):
        vs = self.visits_of(x)
        if not vs:
            raise ValueError("Crossing c%d is not present in this diagram" % (x + 1))
        for v in vs:
            self.visit_over[v] = not self.visit_over[v]
        self._touched = True
        self._dt_labels_valid = False
        self.display_source = "continuity"
        self.last_snap_result = None

    # ---- Reidemeister I : adjacent equal visits ----
    def _reduce_R1_once(self):
        for ci, comp in enumerate(self.comps):
            n = len(comp)
            if n < 2:
                continue
            for i in range(n):
                a, b = comp[i], comp[(i + 1) % n]
                if self.visit_x[a] == self.visit_x[b] and a != b:
                    keep = [v for v in comp if v not in (a, b)]
                    self.comps[ci] = keep
                    for v in (a, b):
                        self.visit_x.pop(v, None)
                        self.visit_over.pop(v, None)
                    self._touched = True
                    self._dt_labels_valid = False
                    return True
        return False

    # ---- Reidemeister II : reducible bigon (verified against planar faces) ----
    def _reduce_R2_once(self):
        for block in self._blocks():
            info = self._block_model(block)
            if info is None:
                continue
            model, _pos2vid = info
            G = D.build_gadget_graph(model)
            try:
                ok, emb = nx.check_planarity(G)
            except Exception:  # noqa: BLE001
                ok = False
            if not ok:
                continue
            for face in D.planar_faces(emb):
                face_t = [tuple(nd) for nd in face]
                segs = [nd for nd in set(face_t) if nd[0] == "seg"]
                xs = set(nd[0] for nd in face_t if nd[0] != "seg")
                if len(segs) != 2 or len(xs) != 2:
                    continue
                p1 = segs[0][1]
                pair1 = self._arc_endpoints(model, p1)
                if pair1 is None:
                    continue
                vx, vy = pair1
                # Reducible R2 poke: same over/under sense at both crossings
                # along one boundary arc of the bigon.
                if model["over_at"][vx] == model["over_at"][vy]:
                    x, y = tuple(xs)
                    xid = model["crossings"][x]["id"]
                    yid = model["crossings"][y]["id"]
                    self._remove_crossings([xid, yid])
                    self._touched = True
                    self._dt_labels_valid = False
                    return True
        return False

    @staticmethod
    def _arc_endpoints(model, p):
        q = model["nextpos"].get(p)
        if q is None:
            return None
        cp, cq = model["pos_cross"].get(p), model["pos_cross"].get(q)
        if cp is None or cq is None or cp == cq:
            return None
        return (p, q)

    def _remove_crossings(self, xids):
        xids = set(xids)
        drop = set(v for v in self.visit_x if self.visit_x[v] in xids)
        for ci, comp in enumerate(self.comps):
            self.comps[ci] = [v for v in comp if v not in drop]
        for v in drop:
            self.visit_x.pop(v, None)
            self.visit_over.pop(v, None)

    def simplify(self, max_moves=10000):
        moves = 0
        while moves < max_moves:
            if self._reduce_R1_once():
                moves += 1
                continue
            if self._reduce_R2_once():
                moves += 1
                continue
            break
        return moves

    # ---- connected blocks (for layout of split links) ----
    def _blocks(self):
        parent = list(range(len(self.comps)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        comp_of_vid = {}
        for ci, comp in enumerate(self.comps):
            for v in comp:
                comp_of_vid[v] = ci
        for x in self.crossings():
            vs = self.visits_of(x)
            if len(vs) == 2 and vs[0] in comp_of_vid and vs[1] in comp_of_vid:
                union(comp_of_vid[vs[0]], comp_of_vid[vs[1]])
        groups = {}
        for ci in range(len(self.comps)):
            groups.setdefault(find(ci), []).append(ci)
        return list(groups.values())

    # ---- build a drawer-model for one block ----
    def _block_model(self, comp_indices):
        comps = [self.comps[ci] for ci in comp_indices]
        xset = sorted(set(self.visit_x[v] for comp in comps for v in comp))
        if not xset:
            return None  # free loop(s) only -> handled by render()
        xindex = {x: k for k, x in enumerate(xset)}

        pos = 1
        comp_positions = []
        pos_cross, pos_role, over_at = {}, {}, {}
        pos2vid = {}
        role_seen = {}
        for comp in comps:
            cp = []
            for v in comp:
                x = self.visit_x[v]
                pos_cross[pos] = xindex[x]
                over_at[pos] = self.visit_over[v]
                if xindex[x] not in role_seen:
                    role_seen[xindex[x]] = pos
                    pos_role[pos] = "o"
                else:
                    pos_role[pos] = "e"
                pos2vid[pos] = v
                cp.append(pos)
                pos += 1
            comp_positions.append(cp)

        nextpos = {}
        prevpos = {}
        comp_of = {}
        for li, cp in enumerate(comp_positions):
            L = len(cp)
            for i, p in enumerate(cp):
                nextpos[p] = cp[(i + 1) % L]
                prevpos[p] = cp[i - 1]
                comp_of[p] = li

        # Each crossing carries its "odd"/"even" visit positions so the drawing
        # helper can colour crossing IDs by the over-strand's component (the
        # helper's default look).  "odd" = first-seen visit (role 'o'), "even" =
        # the second visit (role 'e').
        crossings = [{"id": x, "odd": None, "even": None} for x in xset]
        for p in range(1, pos):
            k = pos_cross[p]
            if pos_role[p] == "o":
                crossings[k]["odd"] = p
            else:
                crossings[k]["even"] = p
        # Guard against any crossing that somehow saw only one visit in this
        # block (should not happen for a valid diagram): fall back sensibly.
        for c in crossings:
            if c["odd"] is None:
                c["odd"] = c["even"]
            if c["even"] is None:
                c["even"] = c["odd"]
        global_colours = [self.component_colors[ci] for ci in comp_indices]
        model = {
            "comp_positions": comp_positions,
            "comp_of": comp_of,
            "nextpos": nextpos,
            "prevpos": prevpos,
            "label_of": (lambda p: ""),
            "pos_cross": pos_cross,
            "pos_role": pos_role,
            "over_at": over_at,
            "crossings": crossings,
            "twon": pos - 1,
            "comp_color_index": global_colours,   # local comp -> original colour
        }
        return model, pos2vid

    # ---- signed DT code (for SnapPy / export); may be None if not DT-valid ----
    def to_signed_dt(self):
        """
        Attempt to emit a signed DT code (SnapPy format) from the current
        diagram by re-traversing components.  Returns ``(dt_string, ok)``.
        """
        num = {}
        pos = 1
        for comp in self.comps:
            for v in comp:
                num[v] = pos
                pos += 1
        twon = pos - 1
        if twon == 0:
            return dt_to_string([tuple() for _ in self.comps]), True

        even_sign = {}
        ok = True
        for x in self.crossings():
            vs = self.visits_of(x)
            labs = sorted(num[v] for v in vs)
            if len(labs) != 2:
                ok = False
                break
            lo, hi = labs
            if lo % 2 == hi % 2:
                ok = False       # not a legal DT code for this traversal
                break
            odd = lo if lo % 2 == 1 else hi
            even = hi if hi % 2 == 0 else lo
            veven = [v for v in vs if num[v] == even][0]
            over = self.visit_over[veven]
            even_sign[odd] = -even if over else even  # neg even = even over
        if not ok:
            return None, False

        comps_out = []
        for comp in self.comps:
            tup = []
            for v in comp:
                n = num[v]
                if n % 2 == 1:
                    tup.append(even_sign[n])
            comps_out.append(tuple(tup))
        return dt_to_string(comps_out), True


# --------------------------------------------------------------------------- #
#  Rendering a (possibly split) diagram with stable colours
# --------------------------------------------------------------------------- #
def color_for(global_comp_index, palette=None):
    palette = palette if palette is not None else D.DEFAULT_PALETTE
    return palette[global_comp_index % len(palette)]


def _layout_like_helper(model, G):
    """Lay out one block exactly as the drawing helper's default 2-D pipeline.

    This mirrors ``draw_dt_original_labelsV3_11.prepare_diagram``: compute the
    helper's default layout, apply its default top-to-bottom orientation and
    rotation, then run the false-crossing audit and fall back to a planar layout
    when the default one produced apparent (false) crossings.  Using the helper's
    own defaults keeps V3.2 diagrams identical to the standalone helper's output.
    """
    P = D.compute_positions(G, DEFAULT_LAYOUT)
    try:
        P = D.transform_positions(P, DEFAULT_Y_DIRECTION, DEFAULT_ROTATE)
    except Exception:  # noqa: BLE001  (older helpers may lack transform_positions)
        pass
    centers = D.crossing_centers(model, P)
    try:
        false = D.audit_false_crossings(model, P, centers)
    except Exception:  # noqa: BLE001
        false = 0
    if false > 0 and DEFAULT_LAYOUT != "planar":
        P = D.compute_positions(G, "planar")
        try:
            P = D.transform_positions(P, DEFAULT_Y_DIRECTION, DEFAULT_ROTATE)
        except Exception:  # noqa: BLE001
            pass
    return P


def render(diagram, ax, show_crossing_ids=True, show_dt_labels=False,
           palette=None, title=None):
    """
    Draw ``diagram`` onto axis ``ax``.  Colours are keyed to original component
    identities via ``diagram.component_colors``.

    Returns ``crossing_id -> (x, y)`` for click hit-testing.
    """
    ax.clear()
    palette = palette if palette is not None else D.DEFAULT_PALETTE
    blocks = diagram._blocks()

    laid = []
    free_loops = []
    for blk in blocks:
        xset = set(diagram.visit_x[v] for ci in blk for v in diagram.comps[ci])
        if xset:
            laid.append(blk)
        else:
            free_loops.extend(blk)

    tile = 2.4
    nslots = max(1, len(laid) + len(free_loops))
    ncols = max(1, int(np.ceil(np.sqrt(nslots))))
    click_targets = {}
    slot = 0

    dt_labels_model = getattr(diagram, "_dt_model", None) if show_dt_labels else None

    for blk in laid:
        info = diagram._block_model(blk)
        if info is None:
            continue
        model, pos2vid = info
        if dt_labels_model is not None and _same_as_original(diagram):
            model = _attach_dt_labels(model, pos2vid, diagram)
            do_labels = True
        else:
            do_labels = False
        G = D.build_gadget_graph(model)
        P = _layout_like_helper(model, G)
        centers = D.crossing_centers(model, P)
        r, c = divmod(slot, ncols)
        origin = (c * tile, -r * tile)
        color_of = (lambda li, model=model:
                    color_for(model["comp_color_index"][li], palette))
        _, cxy = D.render_diagram(
            ax, model, P, centers, color_of=color_of,
            show_labels=do_labels, show_crossing_ids=show_crossing_ids,
            color_crossing_ids_by_overstrand=True,
            arrows=True, origin=origin, scale_to=tile * 0.72)
        for k, xy in cxy.items():
            click_targets[model["crossings"][k]["id"]] = xy
        slot += 1

    for ci in free_loops:
        r, c = divmod(slot, ncols)
        cx, cy = c * tile, -r * tile
        th = np.linspace(0, 2 * np.pi, 200)
        rad = tile * 0.32
        original_colour = diagram.component_colors[ci]
        col = color_for(original_colour, palette)
        ax.plot(cx + rad * np.cos(th), cy + rad * np.sin(th), color=col, lw=3.0)
        ax.annotate("", xy=(cx + rad, cy + 0.02 * tile),
                    xytext=(cx + rad, cy - 0.02 * tile),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=0,
                                    mutation_scale=18))
        ax.text(cx, cy, "component %d" % (original_colour + 1),
                ha="center", va="center", fontsize=8, color="0.35")
        slot += 1

    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=11)
    ax.margins(0.05)
    return click_targets


def _same_as_original(diagram):
    """True when the drawing is a pristine picture of its own ``_dt_model``.

    V3.2: this is driven by ``_dt_labels_valid`` rather than ``_touched``.  Any
    diagram built from a DT code (original, direct after-passage, or the
    SnapPy-simplified DT) is pristine and shows its own DT traversal labels, so
    labels remain visible even after a SnapPy simplification.  The flag is
    cleared once the diagram is passage-changed or Reidemeister-reduced.
    """
    m = getattr(diagram, "_dt_model", None)
    if m is None:
        return False
    if not getattr(diagram, "_dt_labels_valid", False):
        return False
    orig_pairs = sum(len(cp) for cp in m["comp_positions"])
    cur = sum(len(c) for c in diagram.comps)
    return orig_pairs == cur


def _attach_dt_labels(model, pos2vid, diagram):
    om = diagram._dt_model
    label_of = {}
    for p, v in pos2vid.items():
        label_of[p] = om["label_of"](v + 1)
    m2 = dict(model)
    m2["label_of"] = lambda p: label_of.get(p, "")
    return m2


# --------------------------------------------------------------------------- #
#  SnapPy helpers and component-identity matching
# --------------------------------------------------------------------------- #
def _as_fraction(x) -> Fraction:
    if isinstance(x, Fraction):
        return x
    try:
        return Fraction(int(x), 1)
    except Exception:  # noqa: BLE001
        pass
    try:
        return Fraction(str(x))
    except Exception:  # noqa: BLE001
        return Fraction(float(x)).limit_denominator()


def _matrix_to_python(M) -> Optional[List[List[Fraction]]]:
    if M is None:
        return None
    try:
        nr = int(M.nrows())
        nc = int(M.ncols())
        return [[_as_fraction(M[i, j]) for j in range(nc)] for i in range(nr)]
    except Exception:  # noqa: BLE001
        pass
    try:
        rows = list(M)
        return [[_as_fraction(x) for x in row] for row in rows]
    except Exception:  # noqa: BLE001
        return None


def matrix_for_display(M: Optional[Sequence[Sequence[Fraction]]]) -> str:
    if M is None:
        return "n/a"

    def f(x):
        x = _as_fraction(x)
        return str(x.numerator) if x.denominator == 1 else "%s/%s" % (x.numerator, x.denominator)

    return "[" + "; ".join("[" + ", ".join(f(x) for x in row) + "]" for row in M) + "]"


def _snappy_linking_matrix(link) -> Optional[List[List[Fraction]]]:
    try:
        return _matrix_to_python(link.linking_matrix())
    except Exception:  # noqa: BLE001
        return None


def _snappy_component_signatures(link) -> Optional[List[Optional[Tuple[str, str, int]]]]:
    """Return inexpensive per-component signatures, if SnapPy supports them."""
    try:
        n = len(link.link_components)
    except Exception:  # noqa: BLE001
        return None
    sigs: List[Optional[Tuple[str, str, int]]] = []
    for i in range(n):
        try:
            sub = link.sublink([i])
            try:
                sub.simplify("global")
            except Exception:  # noqa: BLE001
                pass
            try:
                jones = str(sub.jones_polynomial())
            except Exception:  # noqa: BLE001
                jones = "jones:n/a"
            try:
                crossings = len(sub.crossings)
            except Exception:  # noqa: BLE001
                crossings = -1
            sigs.append(("component", jones, crossings))
        except Exception:  # noqa: BLE001
            sigs.append(None)
    return sigs


def _signature_filter_available(source_sigs, simplified_sigs, c: int) -> bool:
    if source_sigs is None or simplified_sigs is None:
        return False
    if len(source_sigs) != c or len(simplified_sigs) != c:
        return False
    return all(s is not None for s in source_sigs) and all(s is not None for s in simplified_sigs)


def match_component_order(
    source_matrix: Optional[Sequence[Sequence[Fraction]]],
    simplified_matrix: Optional[Sequence[Sequence[Fraction]]],
    source_signatures: Optional[Sequence[Any]] = None,
    simplified_signatures: Optional[Sequence[Any]] = None,
) -> ComponentMatchResult:
    """Enumerate component permutations and orientation reversals.

    Returns a mapping ``simplified component -> source component`` only when the
    permutation is unique.  Orientation signs are intentionally not used for
    colouring; they only allow linking matrices to match when SnapPy reverses a
    component orientation during simplification.
    """
    if source_matrix is None or simplified_matrix is None:
        return ComponentMatchResult(
            status="unavailable",
            message="linking matrix unavailable; cannot recover SnapPy component labels")
    c = len(source_matrix)
    if c != len(simplified_matrix):
        return ComponentMatchResult(
            status="component-count-mismatch",
            message="component count changed during SnapPy simplification")
    if c == 0:
        return ComponentMatchResult(
            status="component-count-mismatch",
            message="SnapPy returned no labelled components")
    if any(len(row) != c for row in source_matrix) or any(len(row) != c for row in simplified_matrix):
        return ComponentMatchResult(
            status="bad-matrix",
            message="linking matrices are not square with the same size")

    use_sigs = _signature_filter_available(source_signatures, simplified_signatures, c)
    matches: List[MatchCandidate] = []
    for perm in itertools.permutations(range(c)):
        if use_sigs:
            sig_ok = True
            for simp_i, src_i in enumerate(perm):
                if simplified_signatures[simp_i] != source_signatures[src_i]:
                    sig_ok = False
                    break
            if not sig_ok:
                continue
        for signs in itertools.product((1, -1), repeat=c):
            ok = True
            for i in range(c):
                for j in range(i + 1, c):
                    lhs = _as_fraction(signs[i] * signs[j]) * _as_fraction(simplified_matrix[i][j])
                    rhs = _as_fraction(source_matrix[perm[i]][perm[j]])
                    if lhs != rhs:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                matches.append(MatchCandidate(tuple(perm), tuple(signs)))

    unique_perms = sorted(set(m.perm for m in matches))
    if len(unique_perms) == 1:
        perm = unique_perms[0]
        return ComponentMatchResult(
            status="unique",
            message="unique component matching by linking matrix%s" %
                    (" + component signatures" if use_sigs else ""),
            mapping={i: perm[i] for i in range(c)},
            candidates=matches,
            unique_permutations=unique_perms,
            used_signatures=use_sigs,
        )
    if len(unique_perms) == 0:
        return ComponentMatchResult(
            status="no-match",
            message="no permutation/orientation assignment matches the linking matrices",
            candidates=matches,
            unique_permutations=unique_perms,
            used_signatures=use_sigs,
        )
    return ComponentMatchResult(
        status="ambiguous",
        message="%d component labelings have identical invariant data" % len(unique_perms),
        candidates=matches,
        unique_permutations=unique_perms,
        used_signatures=use_sigs,
    )


def _try_backtrack(link, steps):
    """Randomly complicate ``link`` by ~``steps`` moves.  Returns True if a
    ``backtrack`` method was found and called (its signature varies across
    SnapPy/Spherogram versions), False if none is available."""
    for call in (lambda: link.backtrack(num_steps=steps),
                 lambda: link.backtrack(steps),
                 lambda: link.backtrack()):
        try:
            call()
            return True
        except TypeError:
            continue
        except Exception:
            continue
    return False


def backtrack_simplify(snappy, link, mode="global", rounds=0, steps=20,
                       target=None):
    """Simplify a SnapPy link, optionally escaping local minima with backtrack.

    ``simplify('global')`` is greedy and can stop at a non-minimal diagram.  When
    ``rounds > 0`` this repeatedly complicates the diagram (``backtrack``) and
    re-simplifies, remembering the fewest-crossing diagram seen, and returns a
    link in that best state.  With ``rounds <= 0`` it is a single plain simplify
    (unchanged behaviour).  If ``target`` is given, it stops as soon as the
    crossing count reaches ``target`` (used by the per-step reconciliation to
    drive one structure down to another's crossing count efficiently).  Never
    raises: on any problem it falls back to the plainly-simplified link.
    """
    try:
        link.simplify(mode)
    except Exception:  # noqa: BLE001
        pass
    if not rounds or int(rounds) <= 0:
        return link

    best_n = len(link.crossings)
    best_dt = None
    try:
        best_dt = dt_to_string(parse_dt_any(link.DT_code()))
    except Exception:  # noqa: BLE001
        best_dt = None

    if target is not None and best_n <= int(target):
        return link

    for _ in range(int(rounds)):
        if not _try_backtrack(link, steps):
            break  # no backtrack in this build: nothing more we can do
        try:
            link.simplify(mode)
        except Exception:  # noqa: BLE001
            pass
        n = len(link.crossings)
        if n < best_n:
            best_n = n
            try:
                best_dt = dt_to_string(parse_dt_any(link.DT_code()))
            except Exception:  # noqa: BLE001
                pass
        if target is not None and best_n <= int(target):
            break

    # If the current state is not the best seen, rebuild from the best DT.
    if best_dt is not None and len(link.crossings) > best_n:
        try:
            rebuilt = snappy.Link(best_dt)
            rebuilt.simplify(mode)
            if len(rebuilt.crossings) <= best_n:
                return rebuilt
        except Exception:  # noqa: BLE001
            pass
    return link


def snappy_global_simplification(
    source_diagram: Diagram,
    dt_for_snappy: Optional[str] = None,
    negative_even: str = "over",
    mode: str = "global",
    backtrack_rounds: int = 0,
    backtrack_steps: int = 20,
) -> Dict[str, Any]:
    """Ask SnapPy for a global simplification, then try to colour it.

    ``source_diagram`` supplies the physical component colours.  ``dt_for_snappy``
    must encode the same link in the same component order as ``source_diagram``.
    The returned dictionary always contains a ``match`` key when SnapPy ran.
    If ``diagram`` is not ``None``, it is a SnapPy-simplified Diagram whose
    component colours have been matched back to the source/original components.
    """
    result: Dict[str, Any] = {
        "available": False,
        "diagram": None,
        "source_dt": dt_for_snappy,
        "mode": mode,
    }
    if dt_for_snappy is None:
        dt_for_snappy, ok = source_diagram.to_signed_dt()
        result["source_dt"] = dt_for_snappy
        if not ok:
            result["error"] = "current diagram could not be exported as a legal DT code"
            return result

    try:
        import snappy  # type: ignore  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        result["error"] = "SnapPy not available: %s" % exc
        return result

    result["available"] = True
    try:
        L0 = snappy.Link(dt_for_snappy)
        source_components = len(L0.link_components)
        result["source_components"] = source_components
        result["source_crossings"] = len(L0.crossings)
        source_matrix = _snappy_linking_matrix(L0)
        source_sigs = _snappy_component_signatures(L0)
        result["source_linking_matrix"] = source_matrix
        result["source_component_signatures"] = source_sigs
    except Exception as exc:  # noqa: BLE001
        result["error"] = "SnapPy could not read source DT: %s" % exc
        return result

    try:
        Ls = snappy.Link(dt_for_snappy)
        Ls = backtrack_simplify(snappy, Ls, mode=mode,
                                rounds=backtrack_rounds, steps=backtrack_steps)
        result["backtrack_rounds"] = int(backtrack_rounds or 0)
        result["backtrack_steps"] = int(backtrack_steps or 0)
        result["simplified_components"] = len(Ls.link_components)
        result["simplified_crossings"] = len(Ls.crossings)
        simplified_matrix = _snappy_linking_matrix(Ls)
        simplified_sigs = _snappy_component_signatures(Ls)
        result["simplified_linking_matrix"] = simplified_matrix
        result["simplified_component_signatures"] = simplified_sigs
        try:
            snap_jones = str(Ls.jones_polynomial())
        except Exception:  # noqa: BLE001
            snap_jones = "n/a"
        result["jones"] = snap_jones
        try:
            raw_dt = Ls.DT_code()
            simplified_comps = parse_dt_any(raw_dt)
            simplified_dt = dt_to_string(simplified_comps)
            result["simplified_dt"] = simplified_dt
            result["simplified_dt_comps"] = simplified_comps
            try:
                result["jones"] = str(snappy.Link(simplified_dt).jones_polynomial())
            except Exception:  # noqa: BLE001
                result["jones"] = snap_jones
        except Exception as exc:  # noqa: BLE001
            result["simplified_dt_error"] = str(exc)
            simplified_comps = None
    except Exception as exc:  # noqa: BLE001
        result["error"] = "SnapPy simplification failed: %s" % exc
        return result

    match = match_component_order(source_matrix, simplified_matrix, source_sigs, simplified_sigs)
    if source_components != result.get("simplified_components"):
        match = ComponentMatchResult(
            status="component-count-mismatch",
            message="SnapPy has %s component(s) after simplification, source has %s" %
                    (result.get("simplified_components"), source_components),
        )
    result["match"] = match

    # V3.1: ALWAYS display the SnapPy simplify('global') structure when SnapPy
    # produced a usable DT code.  If the component match is unique we recover the
    # original colours; otherwise we still show SnapPy's structure but with a
    # default colouring and flag colours_tracked = False so the GUI can note it.
    if simplified_comps is not None and len(simplified_comps) > 0:
        component_colours, tracked = _snappy_colours(
            source_diagram, simplified_comps, match)
        try:
            diag = Diagram.from_dt_code(simplified_comps,
                                        negative_even=negative_even,
                                        component_colors=component_colours)
            diag._touched = True
            diag.colours_tracked = bool(tracked) and getattr(
                source_diagram, "colours_tracked", True)
            diag.display_source = ("snappy-global"
                                   if diag.colours_tracked
                                   else "snappy-global (colours not tracked)")
            diag.last_snap_result = result
            result["diagram"] = diag
            result["colours_tracked"] = diag.colours_tracked
        except Exception as exc:  # noqa: BLE001
            result["diagram_error"] = "Could not build simplified diagram: %s" % exc
            result["colours_tracked"] = False
    else:
        result["colours_tracked"] = False
    return result


def _snappy_colours(source_diagram: Diagram, simplified_comps, match):
    """
    Decide the colour (original-component index) for each SnapPy-simplified
    component.  Returns ``(component_colours, tracked)``.

    * unique invariant match  -> map each simplified component back to the source
      component's colour; tracked = True.
    * otherwise               -> default identity colouring; tracked = False.
    """
    c = len(simplified_comps)
    source_colours = list(source_diagram.component_colors)
    if (match is not None and match.is_unique and match.mapping is not None
            and len(source_colours) == c
            and all(i in match.mapping for i in range(c))):
        return [source_colours[match.mapping[i]] for i in range(c)], True
    # Colours cannot be tracked: fall back to a stable default colouring so the
    # SnapPy structure is still drawable and clickable.
    return list(range(c)), False


def snappy_properties(diagram: Diagram, dt_for_snappy: Optional[str] = None,
                      negative_even: str = "over") -> Dict[str, Any]:
    """Compatibility wrapper returning invariants/properties for GUI reports."""
    props = {
        "components": diagram.num_components(),
        "crossings_drawn": len(diagram.crossings()),
    }
    if dt_for_snappy is None:
        dt_for_snappy, ok = diagram.to_signed_dt()
    else:
        ok = True
    props["dt_code"] = dt_for_snappy if ok else "(no legal DT code for this traversal)"
    if not ok:
        props["snappy"] = "skipped (no DT code)"
        return props
    snap = snappy_global_simplification(diagram, dt_for_snappy,
                                        negative_even=negative_even)
    props.update(snap)
    if not snap.get("available"):
        props["snappy"] = "not available"
    return props
