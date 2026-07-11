"""Audit a Sphere-XYZ export: does the 3D polyline realize the same link
as the DT code it was generated from?

Parses the .xyz (blank-line-separated components), closes each loop,
projects along a generic direction, extracts crossings with over/under
from interpolated depth, rebuilds the link in spherogram, and compares
hyperbolic invariants with the reference DT link.

Usage: sage -python audit_xyz.py file.xyz "DT: [...]"
"""
import argparse
import sys
import warnings

warnings.filterwarnings('ignore')
import numpy as np

spherogram = None
snappy = None


def _load_topology_modules():
    """Load optional topology dependencies only when an audit is requested."""
    global spherogram, snappy
    if spherogram is None or snappy is None:
        import spherogram as _spherogram
        import snappy as _snappy
        spherogram, snappy = _spherogram, _snappy


def load_xyz(path):
    comps, cur = [], []
    for line in open(path):
        line = line.strip()
        if not line:
            if len(cur) > 2:
                comps.append(np.array(cur))
            cur = []
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                cur.append([float(x) for x in parts[-3:]])
            except ValueError:
                pass
    if len(cur) > 2:
        comps.append(np.array(cur))
    return comps


def link_from_curve(comps, rng):
    """Build a spherogram Link from closed 3D polylines via a random
    projection.  Returns None if the projection is degenerate."""
    # random orthonormal frame
    z = rng.normal(size=3); z /= np.linalg.norm(z)
    x = np.cross(z, rng.normal(size=3)); x /= np.linalg.norm(x)
    y = np.cross(z, x)
    segs = []          # (comp, idx, P2a, P2b, deptha, depthb)
    for ci, pts in enumerate(comps):
        P = np.vstack([pts, pts[:1]])          # close the loop
        u = P @ x; v = P @ y; d = P @ z
        for i in range(len(pts)):
            segs.append((ci, i, np.array([u[i], v[i]]),
                         np.array([u[i + 1], v[i + 1]]), d[i], d[i + 1]))
    crossings = []     # (segA, segB, tA, tB, over_is_A)
    # spatial-grid prefilter: only test segment pairs whose 2D midpoints
    # fall in neighboring cells (cell size ~ max segment extent)
    mids = np.array([(s[2] + s[3]) * 0.5 for s in segs])
    ext = max(1e-9, 2.0 * max(np.linalg.norm(s[3] - s[2]) for s in segs))
    cells = {}
    keys = np.floor(mids / ext).astype(int)
    for idx, k in enumerate(map(tuple, keys)):
        cells.setdefault(k, []).append(idx)
    cand = set()
    for (kx, ky), lst in cells.items():
        neigh = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                neigh.extend(cells.get((kx + dx, ky + dy), ()))
        for a in lst:
            for b in neigh:
                if b > a:
                    cand.add((a, b))
    for a, b in cand:
        ca, ia, A0, A1, da0, da1 = segs[a]
        if True:
            cb, ib, B0, B1, db0, db1 = segs[b]
            if ca == cb and (ia == ib or abs(ia - ib) == 1 or
                             {ia, ib} == {0, len(comps[ca]) - 1}):
                continue
            r = A1 - A0; s = B1 - B0
            den = r[0] * s[1] - r[1] * s[0]
            if abs(den) < 1e-12:
                continue
            q = B0 - A0
            t = (q[0] * s[1] - q[1] * s[0]) / den
            uu = (q[0] * r[1] - q[1] * r[0]) / den
            if not (1e-9 < t < 1 - 1e-9 and 1e-9 < uu < 1 - 1e-9):
                continue
            dA = da0 + t * (da1 - da0)
            dB = db0 + uu * (db1 - db0)
            if abs(dA - dB) < 1e-9:
                return None
            crossings.append((a, b, t, uu, dA > dB))
    # order passes along each component
    from collections import defaultdict
    passes = defaultdict(list)
    for k, (a, b, t, uu, overA) in enumerate(crossings):
        ca, ia = segs[a][0], segs[a][1]
        cb, ib = segs[b][0], segs[b][1]
        passes[ca].append((ia + t, k, True, overA))
        passes[cb].append((ib + uu, k, False, not overA))
    for c in passes:
        passes[c].sort()
    Cs = {k: spherogram.Crossing(str(k)) for k in range(len(crossings))}
    info = {}
    for c in passes:
        for j, (pos, k, isA, over) in enumerate(passes[c]):
            info.setdefault(k, {})['over' if over else 'under'] = (c, j)
    # strand directions at the crossing (2D)
    def dirat(k, want_over):
        a, b, t, uu, overA = crossings[k]
        seg = segs[a] if (overA == want_over) else segs[b]
        d2 = seg[3] - seg[2]
        return d2 / np.linalg.norm(d2)
    slots = {}
    for k, rec in info.items():
        (ui, uj) = rec['under']; (oi, oj) = rec['over']
        du = dirat(k, False); do = dirat(k, True)
        cross = du[0] * do[1] - du[1] * do[0]
        s_in = 1 if cross > 0 else 3
        slots[(ui, uj)] = (Cs[k], 0, 2)
        slots[(oi, oj)] = (Cs[k], s_in, (s_in + 2) % 4)
    for c in passes:
        m = len(passes[c])
        if m == 0:
            return 'split'
        for j in range(m):
            c1, _, e1o = slots[(c, j)]
            c2, e2i, _ = slots[(c, (j + 1) % m)]
            c1[e1o] = c2[e2i]
    return spherogram.Link(list(Cs.values()))



def audit_components_against_dt(xyz_components, dt_string, attempts=4,
                                simplify_rounds=8, seed=11):
    """Library API for draw_dt_original_labels: audit a list of Nx3 closed
    polylines against the DT code they were generated from.

    Returns a dict: {'status': 'ok'|'FAIL'|'inconclusive'|'unavailable',
    'detail': str}.  Requires spherogram + snappy ('unavailable' if not)."""
    try:
        _load_topology_modules()
    except Exception as e:
        return {'status': 'unavailable',
                'detail': 'spherogram/snappy not importable (%s)' % e}
    import numpy as _np
    comps = [_np.asarray(c, float) for c in xyz_components]
    try:
        ref = spherogram.Link(dt_string)
        Mref = ref.exterior()
        vref = float(Mref.volume())
    except Exception as e:
        return {'status': 'inconclusive',
                'detail': 'reference DT link failed: %s' % str(e)[:80]}
    rng = np.random.default_rng(seed)
    last = 'no valid projection'
    for _ in range(attempts):
        L = link_from_curve(comps, rng)
        if L is None or L == 'split':
            last = 'degenerate projection'
            continue
        try:
            for _r in range(simplify_rounds):
                L.simplify('global'); L.simplify('pickup')
            if len(L.crossings) == 0:
                return {'status': 'FAIL',
                        'detail': 'curve reduces to an UNLINK; the 3D '
                                  'polyline does not realize the DT link'}
            M = L.exterior()
            for _r in range(5):
                st = M.solution_type()
                if 'positive' in st:
                    break
                M.randomize()
            v = float(M.volume())
            try:
                iso = M.is_isometric_to(Mref)
            except Exception:
                iso = None
            if iso is True:
                return {'status': 'ok',
                        'detail': 'curve link isometric to DT link '
                                  '(vol %.4f)' % v}
            if iso is False:
                return {'status': 'FAIL',
                        'detail': 'curve link vol %.4f != DT link vol %.4f '
                                  '(NOT isometric): a strand passed through '
                                  'another' % (v, vref)}
            last = ('isometry undecided (curve vol %.4f, DT vol %.4f)'
                    % (v, vref))
        except Exception as e:
            last = 'invariant error: %s' % str(e)[:70]
    return {'status': 'inconclusive', 'detail': last}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct a link from a blank-line-separated XYZ curve and "
            "compare it with the signed DT link using SnapPy isometry checks."
        )
    )
    parser.add_argument("xyz_file", help="XYZ file to audit")
    parser.add_argument("dt", help='Reference signed DT code, e.g. "DT: [(4,6,2)]"')
    parser.add_argument("--attempts", type=int, default=4,
                        help="random projection attempts (default: 4)")
    parser.add_argument("--simplify-rounds", type=int, default=8,
                        help="Spherogram simplify rounds per attempt (default: 8)")
    parser.add_argument("--seed", type=int, default=11,
                        help="random projection seed (default: 11)")
    args = parser.parse_args(argv)
    try:
        _load_topology_modules()
    except Exception as exc:
        parser.error("spherogram and snappy are required: %s" % exc)

    path, dt = args.xyz_file, args.dt
    comps = load_xyz(path)
    print(f"{len(comps)} components, {sum(len(c) for c in comps)} points")
    result = audit_components_against_dt(
        comps,
        dt,
        attempts=max(1, args.attempts),
        simplify_rounds=max(0, args.simplify_rounds),
        seed=args.seed,
    )
    prefix = {
        "ok": "[ok]",
        "FAIL": "[FAIL]",
        "inconclusive": "[WARN]",
        "unavailable": "[info]",
    }.get(result["status"], "[WARN]")
    print("%s xyz topology audit: %s" % (prefix, result["detail"]))
    return 1 if result["status"] == "FAIL" else 0


if __name__ == '__main__':
    sys.exit(main())
