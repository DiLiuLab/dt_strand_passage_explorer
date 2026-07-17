#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
canonical_dt_V2_0.py  --  Canonical DT code (UP TO MIRROR) and 3-D symmetry of a
knot/link diagram.

Given a signed Dowker-Thistlethwaite (DT) code, this computes:

  * the CANONICAL DT code -- the single, reproducible representative of the
    diagram, taken as the lexicographically smallest legal DT code over all
    equivalent re-encodings (loop order, start point, direction) *and the
    whole-code sign inversion*.  Including the sign inversion means the diagram
    and its mirror image collapse to ONE canonical form: this build works
    strictly UP TO MIRROR.  (Rationale: a DT code carries over/under but not
    crossing handedness, so it cannot reliably separate a diagram from its
    mirror anyway -- see the chirality note.  V2_0 therefore stops pretending to
    and merges mirror pairs on purpose, which also makes de-duplication
    self-consistent.)
  * the SYMMETRY of the diagram, determined from the 3-D RENDERING (a symmetric
    Kamada-Kawai layout of the crossings on a sphere), not from the sign-blind
    DT re-encoding count.  Each combinatorial symmetry element is fitted to the
    3-D crossing positions by an orthogonal map; det = +1 marks a rotation, -1 a
    reflection/improper operation, and the fit residual is a built-in
    reliability check.  The result is reported as a rotation order + reflection
    count + a point-group label, all up to mirror (the natural setting, since
    the 3-D shadow is itself sign-blind).

Chirality note: a DT code carries over/under but not crossing handedness (the
writhe sign is what distinguishes a diagram from its mirror), and the over/under
sign is stated relative to the even/odd labelling, so a basepoint shift can flip
it.  Hence the DT code cannot in general separate a diagram from its mirror.
V2_0 embraces this and works up to mirror.  For genuine chirality use a signed
representation (PD code / crossing signs) or SnapPy on the oriented complement.

3-D symmetry needs NumPy and the drawing engine (draw_dt_original_labels*.py) in
the same folder; if they are missing it falls back to the sign-blind
re-encoding count and says so.

Usage
-----
    python3 canonical_dt_V2_0.py --dt "DT: [(4,6,2)]"
    python3 canonical_dt_V2_0.py --dt "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]" --ops
    python3 canonical_dt_V2_0.py                     # no arguments -> graphical interface
    python3 canonical_dt_V2_0.py --gui               # graphical interface

Notes
-----
The DT sign convention is the usual "negative even label => the even-numbered
visit is the over-strand".  The canonical search is exact; its size is
(number of loops)! x product over loops of (2 x loop length) x 2 (the mirror
flip), so it is fast for the small diagrams typical of DNA topology.
"""

import argparse
import ast
import glob
import math
import os
import re
import sys
from itertools import permutations, product

VERSION = "2.0"


# --------------------------------------------------------------------------- #
#  Parsing and the diagram model  (self-contained)
# --------------------------------------------------------------------------- #
def parse_dt(text):
    """Parse 'DT: [(-8,-12,16),(...)]' (or a bare list) into a list of tuples of
    nonzero even integers, one tuple per link component."""
    text = str(text)
    for ch in ("−", "–", "—", "‐", "‑", "－", "―"):
        text = text.replace(ch, "-")           # normalise unicode dashes to '-'
    for ch in (" ", " ", " ", " ", " ", "　"):
        text = text.replace(ch, " ")           # normalise unicode spaces
    m = re.search(r"\[.*\]", text.strip(), re.DOTALL)
    if not m:
        raise ValueError("Could not find a '[...]' list in the DT input.")
    raw = ast.literal_eval(m.group(0))
    if not isinstance(raw, (list, tuple)) or len(raw) == 0:
        raise ValueError("Empty or invalid DT code.")
    if all(isinstance(x, int) for x in raw):   # single-component shorthand [4,6,2]
        raw = [tuple(raw)]
    comps = []
    for ci, comp in enumerate(raw, 1):
        if not isinstance(comp, (list, tuple)):
            raise ValueError("Component %d is not a list/tuple." % ci)
        clean = []
        for x in comp:
            if not isinstance(x, int) or x == 0 or abs(x) % 2 != 0:
                raise ValueError("DT entries must be nonzero even integers; got %r." % (x,))
            clean.append(int(x))
        comps.append(tuple(clean))
    return comps


def build_tours(comps):
    """Return (tours, n): tours[i] is the cyclic list of (crossing_index, is_over)
    that component i passes through, and n is the number of crossings.

    Positions 1..2n are the traversal steps (= DT labels); each component occupies
    a consecutive block.  The k-th crossing (by smallest odd label) has an odd
    label and an even label; a negative even DT entry means the even-numbered
    visit is the over-strand."""
    pos = 1
    odd_partner = {}
    comp_positions = []
    for comp in comps:
        cp = []
        for signed_even in comp:
            cp.append(pos)              # odd-label traversal step
            cp.append(pos + 1)          # even-label traversal step
            odd_partner[pos] = signed_even
            pos += 2
        comp_positions.append(cp)
    twon = pos - 1
    n = twon // 2
    evens = sorted(abs(v) for v in odd_partner.values())
    if evens != list(range(2, twon + 1, 2)):
        raise ValueError("Invalid DT code: the absolute even labels must be exactly "
                         "2, 4, ..., %d with no repeats." % twon)
    pos_cross, over_at = {}, {}
    for k, odd_pos in enumerate(sorted(odd_partner)):
        signed = odd_partner[odd_pos]
        even_pos = abs(signed)
        pos_cross[odd_pos] = k
        pos_cross[even_pos] = k
        even_over = signed < 0          # convention: negative even => even visit is over
        over_at[even_pos] = even_over
        over_at[odd_pos] = not even_over
    tours = [[(pos_cross[p], over_at[p]) for p in cp] for cp in comp_positions]
    return tours, n


# --------------------------------------------------------------------------- #
#  Canonicalization
# --------------------------------------------------------------------------- #
def _variants(tour, flip):
    """All (start, direction) re-readings of one component's cyclic tour, tagged
    with (start_index, reversed?, sequence).  ``flip`` swaps over/under (mirror)."""
    base = [(c, (not o) if flip else o) for c, o in tour]
    L = len(base)
    out = []
    for reversed_flag, seq in ((0, base), (1, base[::-1])):
        for s in range(L):
            out.append((s, reversed_flag, seq[s:] + seq[:s]))
    return out


def _walk_to_dt(walk, n, bounds):
    """Re-derive a signed DT code (tuple of tuples) from a full traversal ``walk``
    (list of (crossing, is_over)); return None if that traversal is not a legal DT
    (some crossing would receive two same-parity labels)."""
    slots = [[] for _ in range(n)]
    for i, (cr, ov) in enumerate(walk):
        slots[cr].append((i + 1, ov))
    signed = {}
    for lst in slots:
        if len(lst) != 2:
            return None
        (p1, o1), (p2, o2) = lst
        if (p1 & 1) == (p2 & 1):
            return None
        if p1 & 1:
            oddp, evenp, ev = p1, p2, o2
        else:
            oddp, evenp, ev = p2, p1, o1
        signed[oddp] = (-evenp if ev else evenp)
    tup = []
    for lo, hi in bounds:
        odds = sorted(p for p in signed if (p & 1) and lo < p <= hi)
        tup.append(tuple(signed[p] for p in odds))
    return tuple(tup)


def canonicalize(comps, allow_flip=False, record_ops=False):
    """Return (canonical_tuple, symmetry_count, ops).

    canonical_tuple : the lexicographically smallest legal DT code (tuple of tuples).
    symmetry_count  : how many re-encodings reproduce it.
    ops             : if record_ops, the list of re-encodings (flip, loop_order,
                      per-loop (start, reversed?), crossing-sequence) achieving it.
    With allow_flip=True the mirror image is also allowed (used for amphichirality).
    """
    tours, n = build_tours(comps)
    C = len(tours)
    Ls = [len(t) for t in tours]
    best = None
    count = 0
    ops = []
    for flip in ((False, True) if allow_flip else (False,)):
        var = [_variants(tours[ci], flip) for ci in range(C)]
        for perm in permutations(range(C)):
            bounds, off = [], 0
            for ci in perm:
                bounds.append((off, off + Ls[ci]))
                off += Ls[ci]
            for choice in product(*[var[ci] for ci in perm]):
                walk = []
                for _, _, seq in choice:
                    walk.extend(seq)
                tup = _walk_to_dt(walk, n, bounds)
                if tup is None:
                    continue
                if best is None or tup < best:
                    best, count = tup, 1
                    if record_ops:
                        ops = [(flip, perm, tuple((choice[i][0], choice[i][1]) for i in range(C)),
                                [c for c, _ in walk])]
                elif tup == best:
                    count += 1
                    if record_ops:
                        ops.append((flip, perm, tuple((choice[i][0], choice[i][1]) for i in range(C)),
                                    [c for c, _ in walk]))
    return best, count, ops


def fmt_dt(tup):
    return "DT: [" + ", ".join("(" + ", ".join(str(x) for x in comp) + ")" for comp in tup) + "]"


def _perm_order(sigma):
    """Order of a permutation given as a dict crossing->crossing (lcm of cycle lengths)."""
    seen, order = set(), 1
    for start in sigma:
        if start in seen:
            continue
        length, x = 0, start
        while True:
            seen.add(x)
            x = sigma[x]
            length += 1
            if x == start:
                break
        order = order * length // math.gcd(order, length)
    return order


# --------------------------------------------------------------------------- #
#  3-D-rendering symmetry  (up to mirror; needs NumPy + the drawing engine)
# --------------------------------------------------------------------------- #
_DRAW_CACHE = {}


def _load_draw_module():
    """Load the highest-versioned draw_dt_original_labels*.py sitting next to us."""
    if "mod" in _DRAW_CACHE:
        return _DRAW_CACHE["mod"]
    import importlib.util
    base = os.path.dirname(os.path.abspath(__file__))
    matches = sorted(glob.glob(os.path.join(base, "draw_dt_original_labels*.py")))
    if not matches:
        raise ImportError("draw_dt_original_labels*.py not found next to canonical_dt_V2_0.py")
    path = matches[-1]
    name = os.path.splitext(os.path.basename(path))[0]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    _DRAW_CACHE["mod"] = mod
    return mod


def _crossing_centers_3d(dt_string):
    """One unit 3-D vector per crossing, from the symmetric Kamada-Kawai sphere
    layout of the diagram graph (the same layout the scorer/figure use)."""
    import numpy as np
    DDOL = _load_draw_module()
    model = DDOL.build_model(DDOL.parse_dt(dt_string))
    G = DDOL.build_gadget_graph(model)
    dirs = DDOL._kamada_3d_unit_directions(G)
    n = len(model["crossings"])
    centers = np.zeros((n, 3), float)
    counts = np.zeros(n, float)
    for node, d in dirs.items():
        cid = node[0] if (isinstance(node, tuple) and isinstance(node[0], int)) else None
        if cid is not None and 0 <= cid < n:
            centers[cid] += np.asarray(d, float)
            counts[cid] += 1.0
    counts[counts == 0] = 1.0
    centers = centers / counts[:, None]
    norms = np.linalg.norm(centers, axis=1, keepdims=True)
    norms[norms < 1e-12] = 1.0
    return centers / norms, n


def _point_group_name(rot_orders, n_refl, n_inv=0, n_S=0):
    """Schoenflies-style label from rotation orders, #mirror planes, #inversion, #S_n.
    Distinguishes Cs (a mirror) from Ci (an inversion centre) -- the fix that makes the
    Balanced clasp read Ci, not Cs."""
    n_rot = len(rot_orders) + 1          # + identity
    n = max(rot_orders) if rot_orders else 1
    order = n_rot + n_refl + n_inv + n_S
    if order <= 1:
        return "C1 (no symmetry)"
    if n_refl == 0 and n_inv == 0 and n_S == 0:
        return "C%d (pure rotation)" % n if n_rot <= 2 else "D%d/C%d-type (order %d, pure rotation)" % (n, n, order)
    if n_rot == 1 and n_inv == 1 and n_refl == 0 and n_S == 0 and order == 2:
        return "Ci (inversion centre)"
    if n_rot == 1 and n_refl == 1 and n_inv == 0 and n_S == 0 and order == 2:
        return "Cs (a single mirror plane)"
    if n == 2 and n_rot == 2 and n_refl == 2 and n_inv == 0 and order == 4:
        return "C2v (one C2 axis + 2 mirror planes)"
    if n == 2 and n_rot == 2 and n_refl == 1 and n_inv == 1 and order == 4:
        return "C2h (C2 axis + inversion centre + mirror plane)"
    if n_rot == 1 and n_refl == 0 and n_inv == 0 and n_S >= 1:
        return "S%d (rotoreflection only)" % (2 * n if n_S else 2)
    bits = []
    if n_rot: bits.append("%d rot(max C%d)" % (n_rot, n))
    if n_refl: bits.append("%d mirror" % n_refl)
    if n_inv: bits.append("inversion")
    if n_S: bits.append("%d S_n" % n_S)
    return "order-%d group (%s)" % (order, ", ".join(bits))


# --------------------------------------------------------------------------- #
#  Rotation-system + eigenvalue symmetry engine (flip=False, loop-gated).
#  A point-group operation of the on-sphere embedding always preserves inside/
#  outside, so it is a flip=False re-encoding.  We classify proper vs improper by
#  the (layout-independent) rotation system, sub-type the impropers (sigma / i /
#  S_n) by the eigenvalues of the orthogonal fit Q, and accept an element only if
#  it also maps the full 3-D LOOPS onto themselves (chamfer gate) -- this rejects
#  the coplanar-crossing false positives that fooled the old det-only method.
# --------------------------------------------------------------------------- #
def _positions_engine(comps):
    pos = 1; odd_partner = {}; comp_positions = []; comp_of = {}
    for ci, comp in enumerate(comps):
        cp = []
        for se in comp:
            odd_partner[pos] = se; cp += [pos, pos + 1]
            comp_of[pos] = ci; comp_of[pos + 1] = ci; pos += 2
        comp_positions.append(cp)
    twon = pos - 1; n = twon // 2; pos_cross = {}; over_at = {}; role = {}
    for k, op in enumerate(sorted(odd_partner)):
        s = odd_partner[op]; ep = abs(s)
        pos_cross[op] = k; pos_cross[ep] = k; role[op] = 'o'; role[ep] = 'e'
        eo = s < 0; over_at[ep] = eo; over_at[op] = not eo
    return comp_positions, pos_cross, over_at, role, comp_of, n


def _walk_to_dt_engine(walk, n, bounds):
    slots = [[] for _ in range(n)]
    for i, (cr, ov) in enumerate(walk):
        slots[cr].append((i + 1, ov))
    signed = {}
    for lst in slots:
        if len(lst) != 2:
            return None
        (p1, o1), (p2, o2) = lst
        if (p1 & 1) == (p2 & 1):
            return None
        if p1 & 1:
            oddp, evenp, ev = p1, p2, o2
        else:
            oddp, evenp, ev = p2, p1, o1
        signed[oddp] = (-evenp if ev else evenp)
    tup = []
    for lo, hi in bounds:
        odds = sorted(p for p in signed if (p & 1) and lo < p <= hi)
        tup.append(tuple(signed[p] for p in odds))
    return tuple(tup)


def _flipfalse_syms(comps, canonical):
    """All flip=False re-encodings reproducing canonical; each returned as
    (sigma_pos, rev_step)."""
    from itertools import permutations, product
    comp_positions, pos_cross, over_at, role, comp_of, n = _positions_engine(comps)
    C = len(comp_positions); Ls = [len(c) for c in comp_positions]
    base = [[(p, pos_cross[p], over_at[p]) for p in cp] for cp in comp_positions]
    def variants(ci):
        b = base[ci]; L = len(b); out = []
        for rev, seq in ((False, b), (True, b[::-1])):
            for s in range(L):
                out.append((rev, seq[s:] + seq[:s]))
        return out
    var = [variants(ci) for ci in range(C)]
    syms = []
    for perm in permutations(range(C)):
        bounds = []; off = 0
        for ci in perm:
            bounds.append((off, off + Ls[ci])); off += Ls[ci]
        for choice in product(*[var[ci] for ci in perm]):
            steps = []
            for rev, seq in choice:
                for (op, cr, ov) in seq:
                    steps.append((op, cr, ov, rev))
            if _walk_to_dt_engine([(cr, ov) for (_, cr, ov, _) in steps], n, bounds) != canonical:
                continue
            sigma_pos = {j + 1: steps[j][0] for j in range(2 * n)}
            rev_step = {j: steps[j][3] for j in range(2 * n)}
            syms.append((sigma_pos, rev_step))
    return syms, n, pos_cross, role, comp_positions


def _corner_ref(canon_dt):
    """ccw cyclic order of the 4 corner-ends at each crossing, from a Tutte layout."""
    import math
    DDOL = _load_draw_module()
    model = DDOL.build_model(DDOL.parse_dt(canon_dt))
    G = DDOL.build_gadget_graph(model)
    pos = DDOL.compute_positions_connected(G, "tutte")
    ncr = len(model["crossings"]); ref = {}
    for k in range(ncr):
        cor = {r: pos[(k, r)] for r in ("in_o", "out_o", "in_e", "out_e")}
        cen = sum(cor.values()) / 4.0
        named = {('o', 'in'): cor['in_o'], ('o', 'out'): cor['out_o'],
                 ('e', 'in'): cor['in_e'], ('e', 'out'): cor['out_e']}
        ang = {kk: math.atan2((v - cen)[1], (v - cen)[0]) for kk, v in named.items()}
        ref[k] = [kk for kk, _ in sorted(ang.items(), key=lambda kv: kv[1])]
    return ref, ncr


def _weave_loops(canon_dt, centers):
    """Full 3-D curve woven through the symmetric crossing centres (over -> outward,
    under -> inward), used as the reliability gate."""
    import math
    import numpy as np
    DDOL = _load_draw_module()
    model = DDOL.build_model(DDOL.parse_dt(canon_dt))
    pcx = model["pos_cross"]; ov = model["over_at"]
    def slerp(a, b, t):
        a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
        d = max(-1.0, min(1.0, float(a @ b))); th = math.acos(d)
        return a if th < 1e-6 else (math.sin((1 - t) * th) * a + math.sin(t * th) * b) / math.sin(th)
    loops = []
    for cp in model["comp_positions"]:
        ks = [pcx[p] for p in cp]; oo = [ov[p] for p in cp]; L = len(ks); pts = []
        for i in range(L):
            k0, k1 = ks[i], ks[(i + 1) % L]
            r0 = 1.16 if oo[i] else 0.84; r1 = 1.16 if oo[(i + 1) % L] else 0.84
            for t in np.linspace(0, 1, 20, endpoint=False):
                pts.append(slerp(centers[k0], centers[k1], t) *
                           ((r0 * (1 - t) + r1 * t) * (1 - 0.10 * math.sin(math.pi * t))))
        loops.append(np.array(pts))
    return np.vstack(loops)


def _chamfer(A, B):
    import numpy as np
    d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)
    return float(np.sqrt(d2.min(1).mean()))


def _eigvec(Q, target):
    import numpy as np
    w, V = np.linalg.eig(Q)
    for i in range(3):
        if abs(w[i].real - target) < 1e-2 and abs(w[i].imag) < 1e-3:
            v = np.real(V[:, i]); return v / (np.linalg.norm(v) + 1e-12)
    return None


def symmetry_3d(canon_comps, resid_tol=0.12, loop_tol=0.06, centers=None, loops=None):
    """True point group of the diagram's 3-D rendering.  `centers`/`loops` may be
    supplied by an external drawer (e.g. the scorer) so the returned axis/plane
    vectors come out in that drawer's coordinate frame.

    For every flip=False re-encoding (all point-group ops preserve inside/outside,
    so flip=False suffices): classify proper vs improper by the layout-independent
    ROTATION SYSTEM; sub-type the impropers by the eigenvalues of the orthogonal
    fit Q -- eig{1,1,-1}=mirror (sigma), {-1,-1,-1}=inversion (i), complex=S_n; and
    KEEP the element only if it also maps the full 3-D loops onto themselves
    (chamfer gate), which rejects the coplanar-crossing false positives that made
    the old det-only method mislabel the Balanced clasp 'Cs' instead of 'Ci'.

    Returns rotation order(s), mirror/inversion/S_n counts, a Schoenflies point
    group, the worst fit residual, and `elements`: a list of (kind, order, vec)
    with kind in {'axis','mirror','inversion','improper-axis'} and vec a unit
    rotation axis / plane normal (None for the inversion centre)."""
    import numpy as np, math
    canon = tuple(tuple(c) for c in canon_comps)
    canon_dt = fmt_dt(canon)
    if centers is None:
        centers, ncr = _crossing_centers_3d(canon_dt)
    else:
        centers = np.asarray(centers, float); ncr = len(centers)
    if loops is None:
        loops = _weave_loops(canon_dt, centers)
    ref, _ = _corner_ref(canon_dt)
    syms, n, pos_cross, role, comp_positions = _flipfalse_syms(list(canon), canon)
    posinfo = {}
    for cp in comp_positions:
        for p in cp:
            posinfo[p] = (pos_cross[p], role[p])
    cen0 = centers.mean(0); Xc = centers - cen0
    Lc = loops - loops.mean(0)
    seen = set()
    rot_orders, n_refl, n_inv, n_S, maxres = [], 0, 0, 0, 0.0
    elements = []
    for (sigma_pos, rev_step) in syms:
        cross_map = {}; M = {}
        for P in range(1, 2 * n + 1):
            Qp = sigma_pos[P]; cQ, rQ = posinfo[Qp]; cP, rP = posinfo[P]; rev = rev_step[P - 1]
            cross_map[cQ] = cP
            if not rev:
                M[(cQ, rQ, 'out')] = (cP, rP, 'out'); M[(cQ, rQ, 'in')] = (cP, rP, 'in')
            else:
                M[(cQ, rQ, 'out')] = (cP, rP, 'in'); M[(cQ, rQ, 'in')] = (cP, rP, 'out')
        p = tuple(cross_map[k] for k in range(ncr))
        if p in seen:
            continue
        seen.add(p)
        # rotation-system orientation (proper vs improper), coplanarity-immune
        votes = []
        for k in range(ncr):
            cyc = ref[k]; kk = cross_map[k]; mapped = []; good = True
            for (r, side) in cyc:
                key = (k, r, side)
                if key not in M:
                    good = False; break
                cP, rP, sp = M[key]; mapped.append((rP, sp))
            if not good:
                continue
            tgt = ref[kk]; i0 = tgt.index(mapped[0])
            votes.append(1 if tgt[(i0 + 1) % 4] == mapped[1]
                         else (-1 if tgt[(i0 - 1) % 4] == mapped[1] else 0))
        if not votes:
            continue
        proper = sum(votes) > 0
        # orthogonal fit for vectors / eigen sub-type + loop gate
        Y = np.array([Xc[p[c]] for c in range(ncr)])
        U, S, Vt = np.linalg.svd(Y.T @ Xc); Q = U @ Vt
        resid = float(np.sqrt(np.mean(np.sum((Xc @ Q.T - Y) ** 2, axis=1))))
        maxres = max(maxres, resid)
        if _chamfer(Lc @ Q.T, Lc) > loop_tol:      # must map the full loops, not just crossings
            continue
        order = _perm_order({i: p[i] for i in range(ncr)})
        if proper:
            if p == tuple(range(ncr)):
                continue                            # identity E
            rot_orders.append(order)
            elements.append(("axis", order, _eigvec(Q, +1.0)))
        else:
            evc = np.linalg.eigvals(Q)
            ev = np.sort(np.round(np.real(evc), 1))
            cplx = bool(np.any(np.abs(np.imag(evc)) > 0.1))
            if np.allclose(ev, [-1.0, -1.0, -1.0]):
                n_inv += 1; elements.append(("inversion", 2, None))
            elif list(ev) == [-1.0, 1.0, 1.0]:
                n_refl += 1; elements.append(("mirror", 2, _eigvec(Q, -1.0)))
            elif cplx:
                n_S += 1; elements.append(("improper-axis", order, _eigvec(Q, -1.0)))
            else:
                n_refl += 1; elements.append(("mirror", 2, _eigvec(Q, -1.0)))
    # de-duplicate near-parallel elements of the same kind
    uniq = []
    for kind, order, vec in elements:
        if vec is not None and any(k == kind and v is not None and
                                   abs(float(np.dot(vec, v))) > 0.98 for k, o, v in uniq):
            continue
        uniq.append((kind, order, vec))
    return {
        "order": 1 + len(rot_orders) + n_refl + n_inv + n_S,   # +1 for identity
        "rotations": len(rot_orders),
        "reflections": n_refl,
        "inversions": n_inv,
        "improper_axes": n_S,
        "rotation_orders": sorted(rot_orders),
        "point_group": _point_group_name(rot_orders, n_refl, n_inv, n_S),
        "max_residual": round(maxres, 4),
        "elements": uniq,
    }


def symmetry_elements(dt_string):
    """Public helper for external drawers (e.g. the scorer's SVG): returns the
    up-to-mirror canonical DT's typed symmetry elements and point group."""
    comps = parse_dt(dt_string)
    canonical, _, _ = canonicalize(comps, allow_flip=True)
    res = symmetry_3d([tuple(c) for c in canonical])
    return {"canonical": fmt_dt(canonical), "point_group": res["point_group"],
            "elements": res["elements"], "max_residual": res["max_residual"]}


def analyze(dt_string):
    """Full analysis of one DT code; returns a results dict."""
    comps = parse_dt(dt_string)
    n_loops = len(comps)
    in_tours, _ = build_tours(comps)
    Ls = [len(t) for t in in_tours]              # visits per loop = crossings it threads
    est = math.factorial(n_loops) * 2            # x2 for the mirror flip (up to mirror)
    for L in Ls:
        est *= 2 * L

    # UP-TO-MIRROR canonical: lex-min over re-encodings AND the whole-code sign flip,
    # so a diagram and its mirror collapse to one canonical form.
    canonical, _, _ = canonicalize(comps, allow_flip=True)
    canon_comps = [tuple(c) for c in canonical]
    tours, _ = build_tours(canon_comps)
    W_id = [c for cp in tours for c, _ in cp]
    n_cross = len(W_id) // 2

    # sign-blind re-encodings that reproduce the canonical (kept for --ops and fallback)
    _, count, ops = canonicalize(canon_comps, record_ops=True)
    op_descr = [(perm, startdirs) for (flip, perm, startdirs, walk) in ops]

    # Order of each re-encoding viewed as a crossing permutation.  This is the
    # exact combinatorial symmetry -- no geometry and no tolerance, so unlike
    # sym3d it is always available; score_diagram* names the group from it.
    orders = []
    for (_flip, _perm, _startdirs, walk) in ops:
        sigma, ok = {}, True
        for i, c in enumerate(W_id):
            if c in sigma and sigma[c] != walk[i]:
                ok = False
            sigma[c] = walk[i]
        orders.append(_perm_order(sigma) if ok and len(sigma) == n_cross else None)
    element_orders = sorted(o for o in orders if o)

    # SYMMETRY from the 3-D rendering (falls back to the sign-blind count if the
    # drawing engine / NumPy are unavailable).
    sym3d, sym_note = None, ""
    try:
        sym3d = symmetry_3d(canon_comps)
    except Exception as exc:  # noqa: BLE001
        sym_note = ("3-D symmetry unavailable (%s) -- reporting the sign-blind "
                    "re-encoding count instead." % exc)

    return {
        "input": dt_string,
        "n_crossings": n_cross,
        "n_components": n_loops,
        "input_strands": Ls,
        "canonical": fmt_dt(canonical),
        "canonical_strands": [len(t) for t in tours],
        "sym3d": sym3d,
        "sym_note": sym_note,
        "reencoding_count": count,
        "symmetry_order": count,          # same number, under the V1 key name
        "element_orders": element_orders,
        "re_encodings_searched": est,
        "operations": op_descr,
    }


# --------------------------------------------------------------------------- #
#  Reporting
# --------------------------------------------------------------------------- #
def report_text(res, list_ops=False):
    lines = []
    lines.append("Input DT code       : %s" % res["input"])
    lines.append("Crossings / loops   : %d crossings, %d components" %
                 (res["n_crossings"], res["n_components"]))
    lines.append("")
    lines.append("CANONICAL DT CODE   : %s   (UP TO MIRROR)" % res["canonical"])
    lines.append("Canonical strands   : %s   (crossings threaded by each loop, canonical order)" %
                 res["canonical_strands"])
    lines.append("(A diagram and its mirror share this canonical form -- V2_0 works up to mirror.)")
    lines.append("")
    s3 = res.get("sym3d")
    if s3 is not None:
        lines.append("SYMMETRY (from the 3-D rendering + rotation system, up to mirror):")
        lines.append("  order            : %d   (incl. identity)" % s3["order"])
        lines.append("  rotation axes    : %d   (orders %s)" %
                     (s3["rotations"], s3["rotation_orders"]))
        lines.append("  mirror planes    : %d" % s3["reflections"])
        lines.append("  inversion centre : %d" % s3.get("inversions", 0))
        lines.append("  rotoreflections  : %d" % s3.get("improper_axes", 0))
        lines.append("  point group      : %s" % s3["point_group"])
        lines.append("  fit residual     : %.4f   (0 = the 3-D layout realises every "
                     "element exactly; each element also passes the full-loop gate)" % s3["max_residual"])
        lines.append("  re-encoding check: %d of %s re-encodings reproduce the canonical"
                     % (res["reencoding_count"],
                        format(res["re_encodings_searched"] // 2, ",")))
    else:
        lines.append("SYMMETRY (sign-blind re-encoding count -- 3-D rendering unavailable):")
        if res.get("sym_note"):
            lines.append("  %s" % res["sym_note"])
        lines.append("  re-encoding count: %d" % res.get("reencoding_count", 0))
    lines.append("")
    lines.append("(searched %s equivalent encodings incl. the mirror flip; kept the smallest)"
                 % format(res["re_encodings_searched"], ","))
    if list_ops and res["operations"]:
        lines.append("")
        lines.append("Symmetry operations (loop order; per-loop start & direction):")
        for i, (perm, sd) in enumerate(res["operations"], 1):
            tag = " (identity)" if (perm == tuple(range(len(perm))) and
                                    all(s0 == 0 and d0 == 0 for s0, d0 in sd)) else ""
            parts = ", ".join("loop%d[start %d, %s]" % (perm[k], sd[k][0],
                              "reversed" if sd[k][1] else "forward")
                              for k in range(len(perm)))
            lines.append("  op%d%s: %s" % (i, tag, parts))
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
#  Graphical interface
# --------------------------------------------------------------------------- #
DEFAULT_DT = "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]"


def launch_gui(defaults=None):
    try:
        import tkinter as tk
        from tkinter import scrolledtext
        root = tk.Tk()
    except Exception as exc:  # no Tk / no display
        print("Graphical interface unavailable (%s); running on the command line instead.\n" % exc)
        dt = getattr(defaults, "dt", None) or DEFAULT_DT
        print(report_text(analyze(dt), list_ops=getattr(defaults, "ops", False)))
        return

    import threading
    import queue as _queue

    root.title("Canonical DT code & symmetry  v%s" % VERSION)
    frm = tk.Frame(root, padx=10, pady=8)
    frm.pack(fill="x")
    frm.columnconfigure(1, weight=1)

    tk.Label(frm, text="DT code", anchor="w").grid(row=0, column=0, sticky="w", pady=3)
    dt_var = tk.StringVar(value=(getattr(defaults, "dt", None) or DEFAULT_DT))
    tk.Entry(frm, textvariable=dt_var, width=70).grid(row=0, column=1, sticky="we", padx=6)

    ops_var = tk.BooleanVar(value=bool(getattr(defaults, "ops", False)))
    tk.Checkbutton(frm, text="list the symmetry operations",
                   variable=ops_var).grid(row=1, column=1, sticky="w", pady=(2, 0))

    tk.Label(frm, anchor="w", justify="left", fg="#444444",
             text="Enter a signed DT code and press Compute. V2_0 gives the canonical DT code UP TO\n"
                  "MIRROR (a diagram and its mirror share it), and the symmetry read from the 3-D\n"
                  "rendering: rotation order, mirror count and point group (fit residual 0 = exact).")\
        .grid(row=2, column=0, columnspan=2, sticky="w", pady=(6, 0))

    btns = tk.Frame(root, padx=10, pady=4)
    btns.pack(fill="x")
    run_btn = tk.Button(btns, text="Compute")
    run_btn.pack(side="left")
    tk.Button(btns, text="Quit", command=root.destroy).pack(side="right")

    out = scrolledtext.ScrolledText(root, width=96, height=22, font=("Menlo", 10))
    out.pack(fill="both", expand=True, padx=10, pady=(4, 10))

    q = _queue.Queue()

    def _poll():
        # All widget updates happen here, on the main thread, driven by the queue.
        try:
            while True:
                kind, payload = q.get_nowait()
                if kind == "clear":
                    out.delete("1.0", "end")
                elif kind == "text":
                    out.insert("end", payload)
                    out.see("end")
                elif kind == "done":
                    run_btn.config(state="normal")     # re-enable for the next run
        except _queue.Empty:
            pass
        root.after(120, _poll)

    def _run():
        dt = dt_var.get().strip() or DEFAULT_DT
        want_ops = ops_var.get()
        run_btn.config(state="disabled")
        q.put(("clear", None))
        q.put(("text", "Computing (this can take a few seconds)...\n\n"))

        def _worker():
            try:
                txt = report_text(analyze(dt), list_ops=want_ops)
            except Exception as exc:  # noqa: BLE001
                txt = "Error: %s" % exc
            q.put(("clear", None))
            q.put(("text", txt + "\n"))
            q.put(("done", None))                      # signal the main thread to re-enable

        threading.Thread(target=_worker, daemon=True).start()

    run_btn.config(command=_run)
    root.after(120, _poll)
    root.mainloop()


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    raw = list(sys.argv[1:]) if argv is None else list(argv)
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dt", default=None, help="the DT code, e.g. \"DT: [(4,6,2)]\"")
    ap.add_argument("--ops", action="store_true", help="also list the symmetry operations")
    ap.add_argument("--gui", action="store_true",
                    help="launch the graphical interface (also the default when no "
                         "arguments are given)")
    args = ap.parse_args(raw)

    if not raw or args.gui:
        launch_gui(args)
        return
    dt = args.dt or DEFAULT_DT
    print(report_text(analyze(dt), list_ops=args.ops))


if __name__ == "__main__":
    main()
