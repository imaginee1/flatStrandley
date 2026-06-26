#!/usr/bin/env python3
"""
Flat "Strand"ley - a sketchpad for protein secondary-structure topology diagrams.

Inspired by Charlie Bond's TopDraw (Tcl/Tk), rebuilt in Python/Tkinter.

Conventions
-----------
* N-terminus = head = beginning of an arrow (flat end) = FULL oval of a cylinder.
* C-terminus = tail = end of an arrow (the arrowhead)  = HALF oval of a cylinder.
* Every end decoration (arrowhead / oval cap) is CENTERED on its grid square.
* At a junction the tail is drawn on top of the head, EXCEPT a loop tail meeting
  a strand head: there the strand head is painted on top.

Geometry + SVG are pure python (no Tk); only the GUI imports tkinter.
"""

import colorsys
import json
import math

# ======================================================================
#  GEOMETRY + RENDER PRIMITIVES (pure python)
# ======================================================================
UNIT = 40.0
OL = 3  # colour-run overlap (samples) to kill sub-pixel seams

DEFAULT_CFG = {
    "sheet_width": 24.0,
    "helix_width": 26.0,
    "loop_width": 7.0,
    "sheet_outline_w": 2.0,
    "helix_outline_w": 2.0,
    "outline": "#1d1d1f",
    "corner_iters": 3,
    "sample_step": 2.5,
    "head_len": 17.0,    # arrowhead length (shorter)
    "head_half": 20.0,   # arrowhead half-width (squished = wide)
    "cap_depth": 0.17,   # oval depth along the tube (smaller = more squished)
    "loops_black": False,  # render all loops black regardless of colour mode
}


def _lerp(a, b, t):
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


def _dist(a, b):
    return math.hypot(b[0] - a[0], b[1] - a[1])


def _dir(a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    L = math.hypot(dx, dy) or 1.0
    return (dx / L, dy / L)


def _chaikin(pts, iters):
    pts = [p for i, p in enumerate(pts) if i == 0 or p != pts[i - 1]]
    for _ in range(iters):
        if len(pts) < 3:
            break
        new = [pts[0]]
        for i in range(len(pts) - 1):
            p, q = pts[i], pts[i + 1]
            new.append((0.75 * p[0] + 0.25 * q[0], 0.75 * p[1] + 0.25 * q[1]))
            new.append((0.25 * p[0] + 0.75 * q[0], 0.25 * p[1] + 0.75 * q[1]))
        new.append(pts[-1])
        pts = new
    return pts


def _resample(pts, step):
    pts = [p for i, p in enumerate(pts) if i == 0 or p != pts[i - 1]]
    if len(pts) < 2:
        return pts[:], [0.0] * len(pts), 0.0
    seglens = [_dist(a, b) for a, b in zip(pts, pts[1:])]
    total = sum(seglens)
    if total == 0:
        return pts[:], [0.0] * len(pts), 0.0
    n = max(1, int(total / step))
    samples, fracs = [], []
    for k in range(n + 1):
        d = min(total, k * step)
        rem, idx = d, 0
        while idx < len(seglens) and rem > seglens[idx]:
            rem -= seglens[idx]
            idx += 1
        if idx >= len(seglens):
            samples.append(pts[-1])
        else:
            t = rem / seglens[idx] if seglens[idx] else 0.0
            samples.append(_lerp(pts[idx], pts[idx + 1], t))
        fracs.append(d / total)
    return samples, fracs, total


def _normals(samples):
    n = len(samples)
    out = []
    for j in range(n):
        a = samples[max(0, j - 1)]
        b = samples[min(n - 1, j + 1)]
        dx, dy = b[0] - a[0], b[1] - a[1]
        L = math.hypot(dx, dy) or 1.0
        out.append((-dy / L, dx / L))
    return out


def _centerline(world_pts, cfg):
    base = _chaikin(world_pts, cfg["corner_iters"])
    samples, fracs, total = _resample(base, cfg["sample_step"])
    return samples, fracs, total, _normals(samples)


def _split_runs(keep):
    runs, start = [], None
    for i, k in enumerate(keep):
        if k and start is None:
            start = i
        elif not k and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(keep)))
    return runs


# --- colour helpers -------------------------------------------------
def hex_blend(c1, c2, t):
    def h2(c):
        c = c.lstrip("#")
        if len(c) == 3:
            c = "".join(ch * 2 for ch in c)
        return tuple(int(c[i:i + 2], 16) for i in (0, 2, 4))
    r1, g1, b1 = h2(c1)
    r2, g2, b2 = h2(c2)
    return "#%02x%02x%02x" % (int(r1 + (r2 - r1) * t),
                              int(g1 + (g2 - g1) * t),
                              int(b1 + (b2 - b1) * t))


def lighten(c, t=0.5):
    return hex_blend(c, "#ffffff", t)


def darken(c, t=0.35):
    return hex_blend(c, "#000000", t)


def _offsets(samples, nrm, idxs, w):
    h = w / 2.0
    left = [(samples[j][0] + nrm[j][0] * h, samples[j][1] + nrm[j][1] * h) for j in idxs]
    right = [(samples[j][0] - nrm[j][0] * h, samples[j][1] - nrm[j][1] * h) for j in idxs]
    return left, right


def _color_runs(idxs, color_of):
    runs, a = [], 0
    for k in range(1, len(idxs)):
        if color_of(idxs[k]) != color_of(idxs[k - 1]):
            runs.append((a, k))
            a = k
    runs.append((a, len(idxs)))
    return runs


def _fill_runs(idxs, color_of, left, right, prims):
    """One filled polygon per colour run, overlapping the next run so there are
    no seams (later runs paint over the overlap)."""
    n = len(idxs)
    for (a, b) in _color_runs(idxs, color_of):
        b2 = min(b + OL, n)
        poly = left[a:b2] + right[a:b2][::-1]
        if len(poly) >= 3:
            prims.append({"t": "poly", "closed": True, "pts": poly,
                          "fill": color_of(idxs[a]), "stroke": None, "sw": 0})


def _stroke(pts, color, w, cap=None):
    p = {"t": "poly", "closed": False, "pts": pts, "fill": None, "stroke": color, "sw": w}
    if cap:
        p["cap"] = cap
    return p


# --- end caps -------------------------------------------------------
def _full_ellipse(prims, E, n, u, rmaj, depth, color, oc, ow):
    """Full oval (cylinder mouth) centred on E. Major axis = width (n),
    minor axis = depth along the tube (u)."""
    steps = 26
    pts = [(E[0] + n[0] * rmaj * math.cos(2 * math.pi * s / steps) + u[0] * depth * math.sin(2 * math.pi * s / steps),
            E[1] + n[1] * rmaj * math.cos(2 * math.pi * s / steps) + u[1] * depth * math.sin(2 * math.pi * s / steps))
           for s in range(steps)]
    prims.append({"t": "poly", "closed": True, "pts": pts,
                  "fill": lighten(color, 0.14), "stroke": None, "sw": 0})
    prims.append({"t": "poly", "closed": True, "pts": pts, "fill": None,
                  "stroke": oc, "sw": ow})


def _half_ellipse(prims, E, n, u, rmaj, depth, color, oc, ow):
    """Half oval (rounded back of a tube) centred on E, bulging along +u.
    Filled; only the curved edge gets a black stroke (the flat side is hidden)."""
    steps = 18
    arc = [(E[0] + n[0] * rmaj * math.cos(math.pi * s / steps) + u[0] * depth * math.sin(math.pi * s / steps),
            E[1] + n[1] * rmaj * math.cos(math.pi * s / steps) + u[1] * depth * math.sin(math.pi * s / steps))
           for s in range(steps + 1)]
    prims.append({"t": "poly", "closed": True, "pts": arc,
                  "fill": color, "stroke": None, "sw": 0})
    prims.append({"t": "poly", "closed": False, "pts": arc, "fill": None,
                  "stroke": oc, "sw": ow, "cap": "round"})


# --- element decomposition -----------------------------------------
def _elements(path, types, dirs):
    """Split a strand into elements: maximal runs of one type (helices and
    sheets also split on direction). Each element keeps the grid edge range
    [e0, e1) and the shared boundary vertices, so neighbours touch exactly."""
    ne = len(path) - 1
    els, i = [], 0
    while i < ne:
        t, d, j = types[i], dirs[i], i + 1
        if t == "loop":
            while j < ne and types[j] == "loop":
                j += 1
        else:
            while j < ne and types[j] == t and dirs[j] == d:
                j += 1
        els.append({"type": t, "e0": i, "e1": j})
        i = j
    return els


def _element_geometry(path, e0, e1, cfg):
    sub = [(c * UNIT, r * UNIT) for c, r in path[e0:e1 + 1]]
    return _centerline(sub, cfg)  # endpoints stay on the grid vertices


def _clip_order(n, direction, clip):
    order = list(range(n))
    if direction < 0:
        order = order[::-1]
    if clip is not None and n > 1:
        m = n - 1
        order = [order[i] for i in range(len(order)) if clip[0] <= i / m <= clip[1]]
    return order


# --- per-element emitters ------------------------------------------
def _emit_loop(prims, samples, color_of, cfg, clip):
    idxs = _clip_order(len(samples), 1, clip)
    if len(idxs) < 2:
        return
    n = len(idxs)
    for (a, b) in _color_runs(idxs, color_of):
        b2 = min(b + OL, n)
        sub = [samples[idxs[i]] for i in range(a, b2)]
        if len(sub) >= 2:
            prims.append(_stroke(sub, color_of(idxs[a]), cfg["loop_width"], "round"))


def _emit_sheet(prims, samples, nrm, color_of, cfg, direction, patch, clip,
                decorate_end=None):
    order = _clip_order(len(samples), direction, clip)
    if len(order) < 2:
        return
    ow, oc = cfg["sheet_outline_w"], cfg["outline"]
    w = cfg["sheet_width"]
    left, right = _offsets(samples, nrm, order, w)
    draw_arrow = ((not patch) or decorate_end == "C") and len(order) > 3
    draw_flatN = (not patch) or decorate_end == "N"

    if not draw_arrow:
        _fill_runs(order, color_of, left, right, prims)
        prims.append(_stroke(left, oc, ow, "round"))
        prims.append(_stroke(right, oc, ow, "round"))
        if draw_flatN:
            prims.append(_stroke([left[0], right[0]], oc, ow, "round"))
        if not patch:
            prims.append(_stroke([left[-1], right[-1]], oc, ow, "round"))
        return

    hl, hh = cfg["head_len"], cfg["head_half"]
    half = hl / 2.0
    end = order[-1]
    E = samples[end]
    u = _dir(samples[order[-2]], E)
    nb = (-u[1], u[0])
    # walk back half the head length so the arrowhead is centred on the C grid
    # vertex: midpoint(base, tip) == E (the grid vertex), tip extends past it.
    acc, bi = 0.0, 0
    for k in range(len(order) - 1, 0, -1):
        acc += _dist(samples[order[k]], samples[order[k - 1]])
        if acc >= half:
            bi = k
            break
    body = order[:bi + 1]
    bl, br = left[:bi + 1], right[:bi + 1]
    _fill_runs(body, color_of, bl, br, prims)
    # side outlines BEFORE the arrowhead so the arrowhead always sits on top
    prims.append(_stroke(bl, oc, ow, "round"))
    prims.append(_stroke(br, oc, ow, "round"))
    if draw_flatN:
        prims.append(_stroke([left[0], right[0]], oc, ow, "round"))   # N flat end (head)

    base = samples[order[bi]]
    barbL = (base[0] + nb[0] * hh, base[1] + nb[1] * hh)
    barbR = (base[0] - nb[0] * hh, base[1] - nb[1] * hh)
    tip = (E[0] + u[0] * half, E[1] + u[1] * half)
    col = color_of(end)
    prims.append({"t": "poly", "closed": True,
                  "pts": [bl[-1], barbL, tip, barbR, br[-1]],
                  "fill": col, "stroke": None, "sw": 0})
    prims.append(_stroke([bl[-1], barbL, tip, barbR, br[-1]], oc, ow, "round"))


def _emit_helix(prims, samples, nrm, color_of, cfg, direction, patch, clip,
                decorate_end=None):
    order = _clip_order(len(samples), direction, clip)
    if len(order) < 2:
        return
    ow, oc = cfg["helix_outline_w"], cfg["outline"]
    w = cfg["helix_width"]
    left, right = _offsets(samples, nrm, order, w)
    _fill_runs(order, color_of, left, right, prims)
    prims.append(_stroke(left, oc, ow, "round"))
    prims.append(_stroke(right, oc, ow, "round"))
    if patch and decorate_end is None:
        return
    depth, r = w * cfg["cap_depth"], w / 2
    if (not patch) or decorate_end == "N":
        h = order[0]                              # head (N) -> full oval, centred on vertex
        uh = _dir(samples[h], samples[order[1]])
        _full_ellipse(prims, samples[h], nrm[h], uh, r, depth, color_of(h), oc, ow)
    if (not patch) or decorate_end == "C":
        t = order[-1]                             # tail (C) -> half oval, centred on vertex
        ut = _dir(samples[order[-2]], samples[t])
        _half_ellipse(prims, samples[t], nrm[t], ut, r, depth, color_of(t), oc, ow)


def strand_prims(path, types, colors, dirs, cfg, clip=None, patch=False,
                 ss_colors=None, color_mode="custom", only_edge=None,
                 decorate_end=None):
    """Render primitives (world px) for one strand.

    Each element is built from its own sub-path so its ends land exactly on the
    grid vertices (decorations are centred there). Elements are emitted C->N so
    N-ward elements paint over their C-ward neighbour (tail-on-top rule).
    clip/only_edge/decorate_end are used to patch a local window over a crossing.
    """
    if len(path) < 2:
        return []
    els = _elements(path, types, dirs)
    out = []  # prim-lists aligned with els (N->C order)
    for el in els:
        e0, e1, t = el["e0"], el["e1"], el["type"]
        if only_edge is not None and not (e0 <= only_edge < e1):
            out.append([])
            continue
        samples, fracs, total, nrm = _element_geometry(path, e0, e1, cfg)
        if total == 0:
            out.append([])
            continue
        ne_el = e1 - e0

        def color_of(j, _e0=e0, _ne=ne_el, _fr=fracs):
            le = _e0 + min(_ne - 1, max(0, int(_fr[j] * _ne)))
            if cfg.get("loops_black") and types[le] == "loop":
                return "#000000"
            if color_mode == "ss" and ss_colors:
                return ss_colors.get(types[le], colors[le])
            return colors[le]

        d = dirs[e0]
        prims = []
        if t == "loop":
            _emit_loop(prims, samples, color_of, cfg, clip)
        elif t == "helix":
            _emit_helix(prims, samples, nrm, color_of, cfg, d, patch, clip,
                        decorate_end)
        else:
            _emit_sheet(prims, samples, nrm, color_of, cfg, d, patch, clip,
                        decorate_end)
        out.append(prims)

    # Draw order: default is C->N so N-ward elements paint over their C-ward
    # neighbour (tail-on-top). Exception: a loop immediately N-ward of a sheet
    # is drawn *under* that sheet, so the strand head paints over the loop tail.
    order = list(range(len(els) - 1, -1, -1))    # [k, k-1, ..., 0]
    pos = {idx: p for p, idx in enumerate(order)}
    for i in range(len(els) - 1):
        if els[i]["type"] == "loop" and els[i + 1]["type"] == "sheet":
            pi, pj = pos[i], pos[i + 1]
            order[pi], order[pj] = order[pj], order[pi]
            pos[i], pos[i + 1] = pj, pi

    res = []
    for idx in order:
        res.extend(out[idx])
    return res


# ----------------------------------------------------------------------
#  CROSSINGS  (overrides store the chosen top segment as (sid, vidx))
# ----------------------------------------------------------------------
def _orient_at(path, i):
    a, b, c = path[i - 1], path[i], path[i + 1]
    din = (b[0] - a[0], b[1] - a[1])
    dout = (c[0] - b[0], c[1] - b[1])
    if din[1] == 0 and dout[1] == 0 and din[0] == dout[0]:
        return "H"
    if din[0] == 0 and dout[0] == 0 and din[1] == dout[1]:
        return "V"
    return None


def _arrow_tips(path, types, dirs):
    """Arrowhead vertices of sheet elements, with the orientation of the segment
    leading into the tip (so an arrowhead landing on another chain registers as
    a crossing that can be toggled)."""
    tips = []
    for el in _elements(path, types, dirs):
        if el["type"] != "sheet":
            continue
        e0, e1, d = el["e0"], el["e1"], dirs[el["e0"]]
        if d > 0:
            tipv, a, b = e1, path[e1 - 1], path[e1]
        else:
            tipv, a, b = e0, path[e0 + 1], path[e0]
        if a[1] == b[1] and a[0] != b[0]:
            tips.append((tipv, "H"))
        elif a[0] == b[0] and a[1] != b[1]:
            tips.append((tipv, "V"))
    return tips


def compute_crossings(strands):
    table = {}
    for s in strands:
        n = len(s.path)
        for i in range(1, n - 1):
            o = _orient_at(s.path, i)
            if o:
                table.setdefault(s.path[i], []).append((s.id, i, o))
        for tipv, o in _arrow_tips(s.path, s.types, s.dirs):
            if tipv in (0, n - 1):  # a terminal arrowhead lying on another chain
                table.setdefault(s.path[tipv], []).append((s.id, tipv, o))
    out = {}
    for p, lst in table.items():
        if any(o == "H" for _, _, o in lst) and any(o == "V" for _, _, o in lst):
            out[p] = lst
    return out


def _default_top(entries):
    """Default top segment: highest strand id, and (for self-crossings) the
    more N-ward vertex of that strand."""
    msid = max(sid for sid, _, _ in entries)
    vid = min(vi for sid, vi, _ in entries if sid == msid)
    return (msid, vid)


def crossing_tops(crossings, overrides):
    out = {}
    for p, lst in crossings.items():
        entries = [(sid, vi) for sid, vi, _ in lst]
        top = overrides.get(p)
        if top not in entries:
            top = _default_top(lst)
        out[p] = top
    return out


def _chain_half_extent(cs, vidx, cfg):
    """Half-width (px) chain cs occupies across its own direction at vidx,
    accounting for a wide arrowhead landing on that vertex."""
    if any(tv == vidx for tv, _ in _arrow_tips(cs.path, cs.types, cs.dirs)):
        return cfg["head_half"]
    i = min(max(vidx - 1, 0), len(cs.types) - 1)
    t = cs.types[i]
    w = (cfg["helix_width"] if t == "helix"
         else cfg["sheet_width"] if t == "sheet" else cfg["loop_width"])
    return w / 2.0


def crossing_patches(strands, crossings, overrides, cfg, ss_colors, color_mode):
    smap = {s.id: s for s in strands}
    tops = crossing_tops(crossings, overrides)
    # every crossing vertex each strand touches (top or not) -> neighbour clamp
    cvids = {}
    for q, lst in crossings.items():
        for es, ev, _o in lst:
            cvids.setdefault(es, set()).add(ev)
    prims = []
    for p in crossings:
        sid, vidx = tops[p]
        s = smap.get(sid)
        if not s:
            continue
        els = _elements(s.path, s.types, s.dirs)
        target = None
        for el in els:                     # prefer a straight interior pass
            if el["e0"] < vidx < el["e1"]:
                target = el
                break
        if target is None:                 # otherwise an end element (arrow tip / cap)
            for el in els:
                if el["e0"] == vidx or el["e1"] == vidx:
                    target = el
                    break
        if target is None:
            continue
        e0, e1 = target["e0"], target["e1"]
        ne_el = e1 - e0
        edge = min(max(vidx - 1, e0), e1 - 1)
        floc = (vidx - e0) / ne_el

        # window must be long enough (along this chain) to fully hide whatever
        # lies beneath at p -- including a wide arrowhead -- so the whole
        # arrowhead clears the crossing instead of poking out.
        others = [(es, ev) for es, ev, _o in crossings[p]
                  if not (es == sid and ev == vidx)]
        if others:
            need_px = max(_chain_half_extent(smap[es], ev, cfg) for es, ev in others)
        else:
            need_px = _chain_half_extent(s, vidx, cfg)
        need_f = (need_px + cfg["sheet_outline_w"] + 2.0) / (ne_el * UNIT)

        # but never bleed past the midpoint toward an adjacent crossing
        inside = [v for v in cvids.get(sid, ()) if e0 <= v <= e1 and v != vidx]
        left = [v for v in inside if v < vidx]
        right = [v for v in inside if v > vidx]
        room_l = (vidx - max(left)) / 2.0 if left else (vidx - e0)
        room_r = (min(right) - vidx) / 2.0 if right else (e1 - vidx)
        wl = min(need_f, room_l / ne_el)
        wr = min(need_f, room_r / ne_el)
        clip = (max(0.0, floc - wl), min(1.0, floc + wr))

        # if the crossing is at this element's decorated terminus, draw it on top
        d = s.dirs[e0]
        dec = None
        if target["type"] == "sheet":
            if vidx == (e1 if d > 0 else e0):
                dec = "C"
        elif target["type"] == "helix":
            if vidx == e0:
                dec = "N" if d > 0 else "C"
            elif vidx == e1:
                dec = "C" if d > 0 else "N"
        prims += strand_prims(s.path, s.types, s.colors, s.dirs, cfg,
                              clip=clip, patch=True, only_edge=edge,
                              decorate_end=dec,
                              ss_colors=ss_colors, color_mode=color_mode)
    return prims


# ----------------------------------------------------------------------
#  SVG EXPORT  (transparent background by default)
# ----------------------------------------------------------------------
def _fmt(pts):
    return " ".join("%.2f,%.2f" % (x, y) for x, y in pts)


def _prim_to_svg(p):
    tag = "polygon" if p.get("closed") else "polyline"
    attrs = ['points="%s"' % _fmt(p["pts"]),
             'fill="%s"' % (p["fill"] if p["fill"] else "none")]
    if p.get("stroke"):
        attrs += ['stroke="%s"' % p["stroke"], 'stroke-width="%.2f"' % p["sw"],
                  'stroke-linejoin="round"']
        if p.get("cap") == "round":
            attrs.append('stroke-linecap="round"')
    else:
        attrs.append('stroke="none"')
    return "  <%s %s/>" % (tag, " ".join(attrs))


def project_to_svg(strands, crossings, overrides, cfg, ss_colors, color_mode,
                   bg=None, margin=30.0):
    xs, ys = [], []
    pad = max(cfg["sheet_width"], cfg["helix_width"]) + cfg["head_half"]
    for s in strands:
        for c, r in s.path:
            xs.append(c * UNIT)
            ys.append(r * UNIT)
    if not xs:
        xs, ys = [0], [0]
    minx, maxx = min(xs) - pad - margin, max(xs) + pad + margin
    miny, maxy = min(ys) - pad - margin, max(ys) + pad + margin
    w, h = maxx - minx, maxy - miny
    parts = []
    for s in strands:
        for p in strand_prims(s.path, s.types, s.colors, s.dirs, cfg,
                              ss_colors=ss_colors, color_mode=color_mode):
            parts.append(_prim_to_svg(p))
    for p in crossing_patches(strands, crossings, overrides, cfg, ss_colors, color_mode):
        parts.append(_prim_to_svg(p))
    bgrect = ('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"/>'
              % (minx, miny, w, h, bg)) if bg else ""
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<svg xmlns="http://www.w3.org/2000/svg" '
            'viewBox="%.2f %.2f %.2f %.2f" width="%.0f" height="%.0f">\n%s\n%s\n</svg>\n'
            % (minx, miny, w, h, w, h, bgrect, "\n".join(parts)))


# ======================================================================
#  GUI LAYER (tkinter)
# ======================================================================
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    HAS_TK = True
except Exception:
    HAS_TK = False

PALETTE = [
    "#1d1d1f", "#5b5b60", "#9a9aa0", "#d0d0d6",
    "#d62728", "#e3791f", "#f2c40f", "#2ca02c",
    "#1a7a3a", "#17becf", "#1f77b4", "#2a3f9e",
    "#6f42c1", "#c44ec4", "#e0529b", "#8c564b",
    "#c7a17a", "#b5651d", "#0d7c66", "#7f3f98",
    "#003f5c", "#bc5090", "#ff6361", "#ffa600",
]

SS_DEFAULTS = {"helix": "#d62728", "sheet": "#1f77b4", "loop": "#5b5b60"}


class Strand:
    _next = 1

    def __init__(self, path, stype="loop", color=None, sid=None, name=None):
        if sid is not None:
            self.id = sid
            Strand._next = max(Strand._next, sid + 1)
        else:
            self.id = Strand._next
            Strand._next += 1
        self.path = [tuple(p) for p in path]
        ne = max(0, len(self.path) - 1)
        col = color or SS_DEFAULTS.get(stype, "#5b5b60")
        self.types = [stype] * ne
        self.colors = [col] * ne
        self.dirs = [1] * ne
        self.name = name or "Chain %d" % self.id
        self.group = None

    @property
    def head(self):
        return self.path[0]

    @property
    def tail(self):
        return self.path[-1]

    def to_dict(self):
        return {"id": self.id, "path": [list(p) for p in self.path],
                "types": self.types, "colors": self.colors,
                "dirs": self.dirs, "name": self.name, "group": self.group}

    @staticmethod
    def from_dict(d):
        s = Strand(d["path"], sid=d["id"], name=d.get("name"))
        s.types = d["types"]
        s.colors = d["colors"]
        s.dirs = d.get("dirs", [1] * len(d["types"]))
        s.group = d.get("group")
        return s


def _substrand(s, v0, v1, sid=None):
    """A new Strand from vertices v0..v1 (edges v0..v1-1) of s, keeping style."""
    a = Strand(s.path[v0:v1 + 1], sid=sid, name=s.name)
    a.types = s.types[v0:v1]
    a.colors = s.colors[v0:v1]
    a.dirs = s.dirs[v0:v1]
    a.group = s.group
    return a


def split_at(s, k):
    """Split s at vertex k into two strands that still meet at vertex k."""
    ne = len(s.path) - 1
    if k <= 0 or k >= ne:
        return None
    return _substrand(s, 0, k), _substrand(s, k, ne)


def delete_vertex(s, k):
    """Remove the tile (vertex k); returns the 0-2 surviving strands."""
    ne = len(s.path) - 1
    out = []
    if k == 0:
        out.append(_substrand(s, 1, ne))
    elif k == ne:
        out.append(_substrand(s, 0, ne - 1))
    else:
        if k >= 2:
            out.append(_substrand(s, 0, k - 1))
        if k <= ne - 2:
            out.append(_substrand(s, k + 1, ne))
    return [x for x in out if len(x.path) >= 2]


def delete_edges(s, a, b):
    """Remove edges [a, b); returns the 0-2 surviving strands."""
    ne = len(s.path) - 1
    out = []
    if a >= 1:
        out.append(_substrand(s, 0, a))
    if b <= ne - 1:
        out.append(_substrand(s, b, ne))
    return [x for x in out if len(x.path) >= 2]


def next_chain_name(used):
    """Smallest 'Chain k' not present in the set `used` (so renamed/deleted
    numbers get reused)."""
    k = 1
    while ("Chain %d" % k) in used:
        k += 1
    return "Chain %d" % k


GROUP_CONFLICT = object()  # returned by _try_merge when the two groups differ


def _try_merge(a, b):
    """Combine only C-terminus -> N-terminus (preserves chain direction).
    Returns the merged Strand, None if the ends don't touch, or GROUP_CONFLICT
    if both chains belong to different (non-empty) groups."""
    if a.tail == b.head:
        path, types = a.path + b.path[1:], a.types + b.types
        colors, dirs, name = a.colors + b.colors, a.dirs + b.dirs, a.name
    elif b.tail == a.head:
        path, types = b.path + a.path[1:], b.types + a.types
        colors, dirs, name = b.colors + a.colors, b.dirs + a.dirs, a.name
    else:
        return None
    if a.group and b.group and a.group != b.group:
        return GROUP_CONFLICT
    s = Strand(path, name=name)
    s.types, s.colors, s.dirs = types, colors, dirs
    s.group = a.group or b.group        # keep the one group that exists
    return s


def _pt_seg_dist(p, a, b):
    px, py = p
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return math.hypot(px - ax, py - ay)
    t = max(0, min(1, ((px - ax) * dx + (py - ay) * dy) / L2))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


if HAS_TK:

    class PanZoomCanvas(tk.Canvas):
        def __init__(self, master, on_render=None, **kw):
            super().__init__(master, highlightthickness=0, **kw)
            self.scale_ = 1.0
            self.ox = 60.0
            self.oy = 60.0
            self.on_render = on_render
            self._pan = None
            self.bind("<ButtonPress-2>", self._pan_start)
            self.bind("<B2-Motion>", self._pan_move)
            self.bind("<MouseWheel>", self._wheel)
            self.bind("<Button-4>", lambda e: self._wheel(e, 120))
            self.bind("<Button-5>", lambda e: self._wheel(e, -120))
            self.bind("<Configure>", lambda e: self.render())

        def w2s(self, x, y):
            return (x * self.scale_ + self.ox, y * self.scale_ + self.oy)

        def s2w(self, x, y):
            return ((x - self.ox) / self.scale_, (y - self.oy) / self.scale_)

        def _pan_start(self, e):
            self._pan = (e.x, e.y, self.ox, self.oy)

        def _pan_move(self, e):
            if self._pan:
                x0, y0, ox0, oy0 = self._pan
                self.ox, self.oy = ox0 + (e.x - x0), oy0 + (e.y - y0)
                self.render()

        def _wheel(self, e, delta=None):
            d = delta if delta is not None else e.delta
            factor = 1.1 if d > 0 else 1 / 1.1
            wx, wy = self.s2w(e.x, e.y)
            self.scale_ = max(0.15, min(6.0, self.scale_ * factor))
            sx, sy = self.w2s(wx, wy)
            self.ox += e.x - sx
            self.oy += e.y - sy
            self.render()

        def reset_view(self, ox=60, oy=60, scale=1.0):
            self.ox, self.oy, self.scale_ = ox, oy, scale
            self.render()

        def render(self):
            if self.on_render:
                self.on_render(self)

        def draw_prims(self, prims):
            s = self.scale_
            for p in prims:
                pts = [self.w2s(x, y) for x, y in p["pts"]]
                flat = [v for xy in pts for v in xy]
                if len(flat) < 4:
                    continue
                if p["t"] == "poly" and p.get("closed"):
                    self.create_polygon(*flat, fill=(p["fill"] or ""),
                                        outline=(p["stroke"] or ""),
                                        width=max(0.4, p["sw"] * s))
                else:
                    kw = dict(width=max(0.4, p["sw"] * s),
                              fill=(p["stroke"] or p["fill"] or "#000"),
                              joinstyle="round")
                    if p.get("cap") == "round":
                        kw["capstyle"] = "round"
                    self.create_line(*flat, **kw)

    class ColorPicker(tk.Toplevel):
        """Friendly HSV picker: hue strip + saturation/value square + hex entry
        + live preview + eyedropper (screen pick via PIL when available)."""

        SQ = 168     # square size
        HUE_W = 22   # hue strip width

        def __init__(self, master, initial, on_ok):
            super().__init__(master)
            self.title("Choose colour")
            self.transient(master)
            self.resizable(False, False)
            self.on_ok = on_ok
            r, g, b = self._hex_rgb(initial)
            self.h, self.s, self.v = colorsys.rgb_to_hsv(r, g, b)

            body = ttk.Frame(self)
            body.pack(padx=10, pady=10)
            self.sq = tk.Canvas(body, width=self.SQ, height=self.SQ,
                                highlightthickness=1, highlightbackground="#999",
                                cursor="crosshair")
            self.sq.grid(row=0, column=0)
            self.hue = tk.Canvas(body, width=self.HUE_W, height=self.SQ,
                                 highlightthickness=1, highlightbackground="#999",
                                 cursor="crosshair")
            self.hue.grid(row=0, column=1, padx=(8, 0))
            self.sq.bind("<Button-1>", self._sq_pick)
            self.sq.bind("<B1-Motion>", self._sq_pick)
            self.hue.bind("<Button-1>", self._hue_pick)
            self.hue.bind("<B1-Motion>", self._hue_pick)

            ctl = ttk.Frame(self)
            ctl.pack(fill="x", padx=10)
            self.preview = tk.Canvas(ctl, width=40, height=24, highlightthickness=1,
                                     highlightbackground="#999")
            self.preview.pack(side="left")
            ttk.Label(ctl, text="Hex").pack(side="left", padx=(8, 2))
            self.hexvar = tk.StringVar()
            e = ttk.Entry(ctl, textvariable=self.hexvar, width=9)
            e.pack(side="left")
            e.bind("<Return>", self._hex_enter)
            ttk.Button(ctl, text="⛏", width=3, command=self._eyedrop).pack(side="left", padx=6)

            btns = ttk.Frame(self)
            btns.pack(fill="x", padx=10, pady=10)
            ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right")
            ttk.Button(btns, text="OK", command=self._ok).pack(side="right", padx=6)

            self._draw_hue()
            self._draw_square()
            self._refresh()
            self.update_idletasks()
            try:
                self.grab_set()
            except Exception:
                pass

        @staticmethod
        def _hex_rgb(h):
            h = h.lstrip("#")
            if len(h) == 3:
                h = "".join(c * 2 for c in h)
            try:
                return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255,
                        int(h[4:6], 16) / 255)
            except Exception:
                return (1.0, 0.0, 0.0)

        def _cur_hex(self):
            r, g, b = colorsys.hsv_to_rgb(self.h, self.s, self.v)
            return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))

        def _draw_hue(self):
            n = 42
            step = self.SQ / n
            for i in range(n):
                r, g, b = colorsys.hsv_to_rgb(i / n, 1, 1)
                col = "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
                self.hue.create_rectangle(0, i * step, self.HUE_W, (i + 1) * step + 1,
                                          fill=col, outline=col)

        def _draw_square(self):
            self.sq.delete("sv")
            n = 28
            step = self.SQ / n
            for i in range(n):            # saturation across x
                for j in range(n):        # value down y (top = bright)
                    s = i / (n - 1)
                    v = 1 - j / (n - 1)
                    r, g, b = colorsys.hsv_to_rgb(self.h, s, v)
                    col = "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
                    self.sq.create_rectangle(i * step, j * step, (i + 1) * step + 1,
                                             (j + 1) * step + 1, fill=col, outline=col,
                                             tags="sv")

        def _refresh(self):
            self.hexvar.set(self._cur_hex())
            self.preview.configure(bg=self._cur_hex())
            self.sq.delete("cursor")
            cx = self.s * self.SQ
            cy = (1 - self.v) * self.SQ
            self.sq.create_oval(cx - 5, cy - 5, cx + 5, cy + 5, outline="#fff",
                                width=2, tags="cursor")
            self.sq.create_oval(cx - 5, cy - 5, cx + 5, cy + 5, outline="#000",
                                tags="cursor")
            self.hue.delete("hcur")
            hy = self.h * self.SQ
            self.hue.create_rectangle(0, hy - 2, self.HUE_W, hy + 2, outline="#000",
                                      width=2, tags="hcur")

        def _sq_pick(self, e):
            self.s = min(1, max(0, e.x / self.SQ))
            self.v = min(1, max(0, 1 - e.y / self.SQ))
            self._refresh()

        def _hue_pick(self, e):
            self.h = min(1, max(0, e.y / self.SQ))
            self._draw_square()
            self._refresh()

        def _hex_enter(self, _e):
            r, g, b = self._hex_rgb(self.hexvar.get())
            self.h, self.s, self.v = colorsys.rgb_to_hsv(r, g, b)
            self._draw_square()
            self._refresh()

        def _eyedrop(self):
            try:
                from PIL import ImageGrab
            except Exception:
                messagebox.showinfo("Eyedropper",
                                    "Screen eyedropper needs the Pillow package.\n"
                                    "Install it with:  pip install pillow\n\n"
                                    "Meanwhile use the ⛏ Pick button on the panel to "
                                    "sample a colour from a chain on the canvas.")
                return
            self.withdraw()
            self.master.update()
            top = tk.Toplevel(self)
            top.attributes("-fullscreen", True)
            top.attributes("-alpha", 0.01)
            top.configure(cursor="crosshair")

            def grab(ev):
                x, y = ev.x_root, ev.y_root
                top.destroy()
                try:
                    px = ImageGrab.grab(bbox=(x, y, x + 1, y + 1)).getpixel((0, 0))
                    hexc = "#%02x%02x%02x" % (px[0], px[1], px[2])
                    self.h, self.s, self.v = colorsys.rgb_to_hsv(
                        px[0] / 255, px[1] / 255, px[2] / 255)
                    self._draw_square()
                    self._refresh()
                except Exception:
                    pass
                self.deiconify()
                try:
                    self.grab_set()
                except Exception:
                    pass
            top.bind("<Button-1>", grab)
            top.bind("<Escape>", lambda e: (top.destroy(), self.deiconify()))

        def _ok(self):
            col = self._cur_hex()
            self.destroy()
            if self.on_ok:
                self.on_ok(col)


    class App(tk.Tk):
        DRAW_TOOLS = ["helix", "sheet", "loop"]
        EDIT_TOOLS = ["select", "selgroup", "move", "split", "crossing", "delete"]
        TOOLS = DRAW_TOOLS + EDIT_TOOLS
        TOOL_LABELS = {
            "select": "Select Chain", "selgroup": "Select Group",
            "helix": "Helix (cylinder)", "sheet": "Strand (arrow)",
            "loop": "Loop (line)", "move": "Move",
            "split": "Split chain", "crossing": "Toggle crossing",
            "delete": "Delete tile",
        }

        def __init__(self):
            super().__init__()
            self.title('Flat "Strand"ley  —  protein topology sketchpad')
            self.geometry("1320x860")
            self.configure(bg="#f4f4f6")

            self.strands = []
            self.cfg = dict(DEFAULT_CFG)
            self.ss_colors = dict(SS_DEFAULTS)
            self.color_mode = tk.StringVar(value="custom")
            self.loops_black = tk.BooleanVar(value=False)
            self.cols = 30
            self.rows = 22
            self.overrides = {}
            self.crossings = {}
            self.tool = tk.StringVar(value="helix")
            self.cur_color = "#d62728"
            self.show_grid = tk.BooleanVar(value=True)
            self.sel = []                 # selected chain ids (multi-select)
            self.draw_path = None
            self.move_ref = None
            self.lin_sel = None
            self._lin_drag = None
            self.undo_stack = []          # last 10 snapshots
            self.cur_file = None
            self.dirty = False

            self._build_style()
            self._build_ui()
            self.recompute()
            self.canvas.reset_view(70, 70, 1.0)
            self.bind_all("<Control-z>", lambda e: self.undo())
            self.bind_all("<Control-Z>", lambda e: self.undo())
            self.bind_all("<Control-s>", lambda e: self.save_project())
            self.protocol("WM_DELETE_WINDOW", self.on_close)

        def _build_style(self):
            st = ttk.Style(self)
            try:
                st.theme_use("clam")
            except Exception:
                pass
            st.configure(".", background="#f4f4f6", font=("Helvetica", 10))
            st.configure("TButton", padding=4)
            st.configure("Panel.TFrame", background="#fafafb")
            st.configure("Toolbar.TFrame", background="#ffffff")
            st.configure("Head.TLabel", font=("Helvetica", 10, "bold"))

        def _build_ui(self):
            top = ttk.Frame(self, style="Toolbar.TFrame")
            top.pack(side="top", fill="x")
            self._topbar(top)

            self.hpane = ttk.PanedWindow(self, orient="horizontal")
            self.hpane.pack(fill="both", expand=True)

            self.tools_panel = ttk.Frame(self.hpane, style="Panel.TFrame", width=148)
            self._tools_panel(self.tools_panel)
            self.hpane.add(self.tools_panel, weight=0)

            self.vpane = ttk.PanedWindow(self.hpane, orient="vertical")
            self.hpane.add(self.vpane, weight=1)

            cwrap = ttk.Frame(self.vpane)
            self.canvas = PanZoomCanvas(cwrap, on_render=lambda c: self.render_main(),
                                        bg="#ffffff")
            self.canvas.pack(fill="both", expand=True)
            self.hover_lbl = tk.Label(self.canvas, text="", bg="#ffffff", fg="#333",
                                      font=("Helvetica", 10, "bold"), anchor="w",
                                      justify="left", bd=0)
            self.hover_lbl.place(x=8, y=8)
            self.vpane.add(cwrap, weight=3)

            self.lin_wrap = ttk.Frame(self.vpane)
            head = ttk.Frame(self.lin_wrap, style="Panel.TFrame")
            head.pack(fill="x")
            ttk.Label(head, text="  Linear view — drag a stretch to recolour / convert / "
                      "reverse;  click a group's box to select it, its name to rename",
                      style="Head.TLabel").pack(side="left", pady=2)
            self.lin = PanZoomCanvas(self.lin_wrap,
                                     on_render=lambda c: self.render_linear(),
                                     bg="#ffffff", height=180)
            self.lin.pack(fill="both", expand=True)
            self.lin.bind("<ButtonPress-1>", self._lin_press)
            self.lin.bind("<B1-Motion>", self._lin_motion)
            self.lin.bind("<ButtonRelease-1>", self._lin_release)
            self.vpane.add(self.lin_wrap, weight=1)

            self.props_panel = ttk.Frame(self.hpane, style="Panel.TFrame", width=240)
            self._props_panel(self.props_panel)
            self.hpane.add(self.props_panel, weight=0)

            self.canvas.bind("<ButtonPress-1>", self._press)
            self.canvas.bind("<B1-Motion>", self._motion)
            self.canvas.bind("<ButtonRelease-1>", self._release)
            self.canvas.bind("<Motion>", self._hover)

            self.status = ttk.Label(self, text="", anchor="w")
            self.status.pack(side="bottom", fill="x")
            self._set_status()

        def _topbar(self, bar):
            def b(t, c):
                ttk.Button(bar, text=t, command=c).pack(side="left", padx=2, pady=4)
            b("New", self.new_project)
            b("Open", self.open_project)
            b("Close", self.close_file)
            b("Save", self.save_project)
            b("Save As", self.save_project_as)
            b("Export SVG", self.export_svg)
            ttk.Separator(bar, orient="vertical").pack(side="left", fill="y", padx=6, pady=4)
            ttk.Button(bar, text="↶ Undo", command=self.undo).pack(side="left", padx=2, pady=4)
            ttk.Separator(bar, orient="vertical").pack(side="left", fill="y", padx=6, pady=4)
            ttk.Label(bar, text="Grid").pack(side="left")
            self.col_var = tk.IntVar(value=self.cols)
            self.row_var = tk.IntVar(value=self.rows)
            tk.Spinbox(bar, from_=4, to=300, width=4, textvariable=self.col_var,
                       command=self.resize_grid).pack(side="left", padx=2)
            ttk.Label(bar, text="x").pack(side="left")
            tk.Spinbox(bar, from_=4, to=300, width=4, textvariable=self.row_var,
                       command=self.resize_grid).pack(side="left", padx=2)
            ttk.Button(bar, text="Apply", command=self.resize_grid).pack(side="left", padx=2)
            ttk.Separator(bar, orient="vertical").pack(side="left", fill="y", padx=6, pady=4)
            ttk.Button(bar, text="Fit", command=self.fit_view).pack(side="left", padx=2)
            ttk.Button(bar, text="–", width=2, command=lambda: self._zoom(1 / 1.2)).pack(side="left")
            ttk.Button(bar, text="+", width=2, command=lambda: self._zoom(1.2)).pack(side="left", padx=2)
            ttk.Separator(bar, orient="vertical").pack(side="left", fill="y", padx=6, pady=4)
            self.t_tools = tk.BooleanVar(value=True)
            self.t_props = tk.BooleanVar(value=True)
            self.t_lin = tk.BooleanVar(value=True)
            ttk.Checkbutton(bar, text="Tools", variable=self.t_tools,
                            command=self._toggle_panels).pack(side="left", padx=2)
            ttk.Checkbutton(bar, text="Linear", variable=self.t_lin,
                            command=self._toggle_panels).pack(side="left", padx=2)
            ttk.Checkbutton(bar, text="Properties", variable=self.t_props,
                            command=self._toggle_panels).pack(side="left", padx=2)
            ttk.Checkbutton(bar, text="Grid", variable=self.show_grid,
                            command=self.render_main).pack(side="left", padx=8)

        def _tools_panel(self, p):
            ttk.Label(p, text="Drawing tools", style="Head.TLabel").pack(
                anchor="w", padx=8, pady=(8, 2))
            for t in self.DRAW_TOOLS:
                ttk.Radiobutton(p, text=self.TOOL_LABELS[t], value=t, variable=self.tool,
                                command=self._tool_changed).pack(anchor="w", padx=10)
            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6, padx=6)
            ttk.Label(p, text="Editing tools", style="Head.TLabel").pack(anchor="w", padx=8)
            for t in self.EDIT_TOOLS:
                ttk.Radiobutton(p, text=self.TOOL_LABELS[t], value=t, variable=self.tool,
                                command=self._tool_changed).pack(anchor="w", padx=10)
            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=8, padx=6)
            ttk.Button(p, text="Join selected ends",
                       command=self.join_ends).pack(fill="x", padx=8, pady=2)
            ttk.Button(p, text="Reverse whole chain",
                       command=self.reverse_selected).pack(fill="x", padx=8, pady=2)
            ttk.Button(p, text="Rename chain…",
                       command=self.rename_chain).pack(fill="x", padx=8, pady=2)
            ttk.Button(p, text="Delete selected",
                       command=self.delete_selected).pack(fill="x", padx=8, pady=2)
            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=8, padx=6)
            ttk.Label(p, text="Groups", style="Head.TLabel").pack(anchor="w", padx=8)
            ttk.Button(p, text="Create Group…",
                       command=self.create_group).pack(fill="x", padx=8, pady=2)
            ttk.Button(p, text="Add to group…",
                       command=self.add_to_group).pack(fill="x", padx=8, pady=2)
            ttk.Button(p, text="Split Group",
                       command=self.split_group).pack(fill="x", padx=8, pady=2)
            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=8, padx=6)
            ttk.Label(p, text="Shift-click: multi-select\nMiddle-drag: pan\n"
                      "Wheel: zoom\nCtrl-Z: undo   Ctrl-S: save",
                      foreground="#777").pack(anchor="w", padx=8)

        def _props_panel(self, p):
            ttk.Label(p, text="Colour", style="Head.TLabel").pack(anchor="w", padx=8, pady=(8, 2))
            sw = ttk.Frame(p)
            sw.pack(padx=8, anchor="w")
            for i, col in enumerate(PALETTE):
                c = tk.Canvas(sw, width=20, height=20, bg=col, highlightthickness=1,
                              highlightbackground="#ccc", cursor="hand2")
                c.grid(row=i // 6, column=i % 6, padx=1, pady=1)
                c.bind("<Button-1>", lambda e, col=col: self._pick_color(col))
            row = ttk.Frame(p)
            row.pack(fill="x", padx=8, pady=4)
            self.cur_swatch = tk.Canvas(row, width=46, height=18, bg=self.cur_color,
                                        highlightthickness=1, highlightbackground="#999")
            self.cur_swatch.pack(side="left")
            ttk.Button(row, text="Custom…", width=8,
                       command=self._choose_color).pack(side="left", padx=4)
            ttk.Button(row, text="⛏ Pick", width=6,
                       command=self._start_eyedropper).pack(side="left")
            ttk.Button(p, text="Apply colour to selection",
                       command=self.apply_stretch_color).pack(fill="x", padx=8, pady=2)
            ttk.Button(p, text="Reverse selected section",
                       command=self.reverse_section).pack(fill="x", padx=8, pady=2)

            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=8, padx=6)
            ttk.Label(p, text="Colour mode", style="Head.TLabel").pack(anchor="w", padx=8)
            ttk.Radiobutton(p, text="Custom colours", value="custom",
                            variable=self.color_mode, command=self.recompute).pack(anchor="w", padx=10)
            ttk.Radiobutton(p, text="By secondary structure", value="ss",
                            variable=self.color_mode, command=self.recompute).pack(anchor="w", padx=10)
            self.ss_btns = {}
            for t, lbl in [("helix", "Helix"), ("sheet", "Strand"), ("loop", "Loop")]:
                r = ttk.Frame(p)
                r.pack(fill="x", padx=12, pady=1)
                ttk.Label(r, text=lbl, width=7).pack(side="left")
                cv = tk.Canvas(r, width=28, height=16, bg=self.ss_colors[t],
                               highlightthickness=1, highlightbackground="#999", cursor="hand2")
                cv.pack(side="left")
                cv.bind("<Button-1>", lambda e, t=t: self._choose_ss(t))
                self.ss_btns[t] = cv
            ttk.Checkbutton(p, text="Show all loops black",
                            variable=self.loops_black,
                            command=self._toggle_loops_black).pack(anchor="w", padx=10, pady=(4, 0))

            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=8, padx=6)
            ttk.Label(p, text="Widths", style="Head.TLabel").pack(anchor="w", padx=8)
            self._slider(p, "Arrow", "sheet_width", 8, 44)
            self._slider(p, "Helix", "helix_width", 8, 44)
            self._slider(p, "Loop", "loop_width", 2, 18)
            self._slider(p, "Arrow edge", "sheet_outline_w", 0.5, 6)
            self._slider(p, "Helix edge", "helix_outline_w", 0.5, 6)

            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=8, padx=6)
            ttk.Label(p, text="Selected stretch", style="Head.TLabel").pack(anchor="w", padx=8)
            self.sel_label = ttk.Label(p, text="(none)", foreground="#555")
            self.sel_label.pack(anchor="w", padx=10)
            ttk.Label(p, text="Convert selection to:", foreground="#555").pack(
                anchor="w", padx=10, pady=(4, 0))
            tr = ttk.Frame(p)
            tr.pack(fill="x", padx=8, pady=2)
            for t, lbl in [("helix", "Helix"), ("sheet", "Strand"), ("loop", "Loop")]:
                ttk.Button(tr, text=lbl, width=6,
                           command=lambda t=t: self.set_stretch_type(t)).pack(side="left", padx=1)

        def _toggle_loops_black(self):
            self.cfg["loops_black"] = bool(self.loops_black.get())
            self.recompute()

        def _slider(self, parent, label, key, lo, hi):
            row = ttk.Frame(parent)
            row.pack(fill="x", padx=8, pady=1)
            ttk.Label(row, text=label, width=10).pack(side="left")
            val = tk.DoubleVar(value=self.cfg[key])
            lbl = ttk.Label(row, text="%.1f" % self.cfg[key], width=4)

            def on(v, k=key, vv=val, ll=lbl):
                self.cfg[k] = float(v)
                ll.configure(text="%.1f" % float(v))
                self.render_main()
            ttk.Scale(row, from_=lo, to=hi, variable=val, command=on).pack(
                side="left", fill="x", expand=True, padx=4)
            lbl.pack(side="left")

        def _toggle_panels(self):
            self._set_pane(self.hpane, self.tools_panel, self.t_tools.get(), 0)
            self._set_pane(self.hpane, self.props_panel, self.t_props.get(), None)
            self._set_pane(self.vpane, self.lin_wrap, self.t_lin.get(), 1)

        def _set_pane(self, pane, child, show, index):
            present = str(child) in pane.panes()
            if show and not present:
                if index is None or index >= len(pane.panes()):
                    pane.add(child)
                else:
                    pane.insert(index, child)
            elif not show and present:
                pane.forget(child)

        def _tool_changed(self):
            if self.tool.get() not in ("select", "selgroup", "move"):
                self.sel = []
                self.lin_sel = None
                self._update_sel_label()
            self.render_main()
            self.render_linear()
            self._set_status()

        def _pick_color(self, col):
            self.cur_color = col
            self.cur_swatch.configure(bg=col)
            self.apply_color()           # recolour current selection immediately

        def _choose_color(self):
            ColorPicker(self, self.cur_color, on_ok=self._pick_color)

        def _start_eyedropper(self):
            self._eyedrop = True
            self.canvas.configure(cursor="crosshair")
            self._set_status("Eyedropper: click a chain (or anywhere) to pick a colour.")

        def _do_eyedrop(self, sx, sy):
            col = None
            try:                                   # true screen pixel, if Pillow present
                from PIL import ImageGrab
                x = self.canvas.winfo_rootx() + sx
                y = self.canvas.winfo_rooty() + sy
                px = ImageGrab.grab(bbox=(x, y, x + 1, y + 1)).getpixel((0, 0))
                col = "#%02x%02x%02x" % (px[0], px[1], px[2])
            except Exception:
                col = None
            if col is None:                        # fall back to the chain under cursor
                hit = self._top_vertex_at(sx, sy)
                if hit:
                    s, k = hit
                    col = s.colors[min(max(k - 1, 0), len(s.colors) - 1)]
            self._eyedrop = False
            self.canvas.configure(cursor="")
            if col:
                self.cur_color = col
                self.cur_swatch.configure(bg=col)
                self._set_status("Picked %s" % col)
            else:
                self._set_status("Eyedropper: nothing under the cursor.")

        def _choose_ss(self, t):
            def done(c):
                self.push_undo()
                self.ss_colors[t] = c
                self.ss_btns[t].configure(bg=c)
                self.recompute()
            ColorPicker(self, self.ss_colors[t], on_ok=done)

        def _set_status(self, extra=""):
            self.status.configure(
                text="  Tool: %s   |   Chains: %d   |   Selected: %d   |   Crossings: %d   %s"
                % (self.tool.get(), len(self.strands), len(self.sel),
                   len(self.crossings), extra))

        def _update_sel_label(self):
            if not self.sel:
                self.sel_label.configure(text="(none)")
            elif len(self.sel) == 1:
                s = self._get(self.sel[0])
                self.sel_label.configure(text=(s.name if s else "(none)"))
            else:
                self.sel_label.configure(text="%d chains" % len(self.sel))

        def _get(self, sid):
            return next((s for s in self.strands if s.id == sid), None)

        def _sel_strands(self):
            return [s for s in self.strands if s.id in self.sel]

        def _dedupe_names(self):
            """Ensure chain names are unique; reused numbers fill the lowest gap."""
            seen = set()
            for s in self.strands:
                if s.name in seen:
                    s.name = next_chain_name(seen | {x.name for x in self.strands})
                seen.add(s.name)

        def _center_window(self, win):
            self.update_idletasks()
            try:
                px, py = self.winfo_rootx(), self.winfo_rooty()
                pw, ph = self.winfo_width(), self.winfo_height()
                ww, wh = win.winfo_width(), win.winfo_height()
                win.geometry("+%d+%d" % (px + (pw - ww) // 2, py + (ph - wh) // 2))
            except Exception:
                pass

        # ---- undo -------------------------------------------------
        def _snapshot(self):
            return {"strands": [s.to_dict() for s in self.strands],
                    "overrides": {"%d,%d" % p: list(t)
                                  for p, t in self.overrides.items()},
                    "sel": list(self.sel)}

        def push_undo(self):
            self.undo_stack.append(self._snapshot())
            if len(self.undo_stack) > 10:
                self.undo_stack.pop(0)
            self.dirty = True

        def undo(self):
            if not self.undo_stack:
                self._set_status("nothing to undo.")
                return
            snap = self.undo_stack.pop()
            self.strands = [Strand.from_dict(d) for d in snap["strands"]]
            self.overrides = {}
            for key, val in snap["overrides"].items():
                cx, cy = key.split(",")
                self.overrides[(int(cx), int(cy))] = tuple(val)
            self.sel = [sid for sid in snap["sel"] if self._get(sid)]
            self.lin_sel = None
            self._update_sel_label()
            self.recompute()

        # ---- colour / groups -------------------------------------
        def apply_color(self):
            """Recolour the current selection (linear stretch, else whole
            selected chains) and update immediately."""
            changed = False
            if self.lin_sel:
                sid, a, b = self.lin_sel
                s = self._get(sid)
                if s:
                    self.push_undo()
                    for k in range(a, b):
                        s.colors[k] = self.cur_color
                    changed = True
            elif self.sel:
                self.push_undo()
                for s in self._sel_strands():
                    s.colors = [self.cur_color] * len(s.colors)
                changed = True
            if changed:
                if self.color_mode.get() == "ss":
                    self._set_status("colour stored — switch to 'Custom colours' to see it.")
                self.recompute()

        def create_group(self):
            if not self.sel:
                messagebox.showinfo("Create Group", "Select one or more chains first.")
                return
            from tkinter import simpledialog
            name = simpledialog.askstring("Create Group", "Group name:", parent=self)
            if not name:
                return
            self.push_undo()
            for s in self._sel_strands():
                s.group = name
            self.recompute()

        def split_group(self):
            names = {s.group for s in self._sel_strands() if s.group}
            if not names:
                messagebox.showinfo("Split Group", "Select a chain that is in a group.")
                return
            self.push_undo()
            for s in self.strands:
                if s.group in names:
                    s.group = None
            self.recompute()

        def rename_group(self, old):
            from tkinter import simpledialog
            name = simpledialog.askstring("Rename Group", "New name:",
                                          initialvalue=old, parent=self)
            if not name:
                return
            self.push_undo()
            for s in self.strands:
                if s.group == old:
                    s.group = name
            self.recompute()

        def rename_chain(self):
            if len(self.sel) != 1:
                messagebox.showinfo("Rename chain", "Select exactly one chain to rename.")
                return
            from tkinter import simpledialog
            s = self._get(self.sel[0])
            name = simpledialog.askstring("Rename chain", "Chain name:",
                                          initialvalue=s.name, parent=self)
            if not name:
                return
            self.push_undo()
            s.name = name
            self.recompute()       # frees the old 'Chain k' for reuse

        def select_group_named(self, gname):
            self.sel = [s.id for s in self.strands if s.group == gname]
            self.lin_sel = None
            self._update_sel_label()
            self.render_main()
            self.render_linear()
            self._set_status()

        def _hover(self, e):
            sid = self._hit_strand(e.x, e.y)
            if sid is None:
                self.hover_lbl.configure(text="")
                return
            s = self._get(sid)
            grp = ("Group: %s   " % s.group) if s and s.group else ""
            self.hover_lbl.configure(text="%s%s" % (grp, s.name if s else ""))

        def _zoom(self, f):
            cx = self.canvas.winfo_width() / 2
            cy = self.canvas.winfo_height() / 2
            wx, wy = self.canvas.s2w(cx, cy)
            self.canvas.scale_ = max(0.15, min(6.0, self.canvas.scale_ * f))
            sx, sy = self.canvas.w2s(wx, wy)
            self.canvas.ox += cx - sx
            self.canvas.oy += cy - sy
            self.render_main()

        def fit_view(self):
            if not self.strands:
                self.canvas.reset_view(70, 70, 1.0)
                return
            xs = [c * UNIT for s in self.strands for c, r in s.path]
            ys = [r * UNIT for s in self.strands for c, r in s.path]
            pad = 60
            minx, maxx = min(xs) - pad, max(xs) + pad
            miny, maxy = min(ys) - pad, max(ys) + pad
            cw = max(50, self.canvas.winfo_width())
            ch = max(50, self.canvas.winfo_height())
            sc = max(0.15, min(6.0, min(cw / (maxx - minx), ch / (maxy - miny))))
            self.canvas.scale_ = sc
            self.canvas.ox = (cw - (maxx - minx) * sc) / 2 - minx * sc
            self.canvas.oy = (ch - (maxy - miny) * sc) / 2 - miny * sc
            self.render_main()

        def grid_at(self, sx, sy):
            wx, wy = self.canvas.s2w(sx, sy)
            c = max(0, min(self.cols, round(wx / UNIT)))
            r = max(0, min(self.rows, round(wy / UNIT)))
            return (int(c), int(r))

        def recompute(self):
            self.crossings = compute_crossings(self.strands)
            self.overrides = {p: t for p, t in self.overrides.items() if p in self.crossings}
            self.render_main()
            self.render_linear()
            self._set_status()

        def render_main(self):
            c = self.canvas
            c.delete("all")
            if self.show_grid.get():
                self._draw_grid(c)
            for s in self.strands:
                c.draw_prims(strand_prims(s.path, s.types, s.colors, s.dirs, self.cfg,
                                          ss_colors=self.ss_colors,
                                          color_mode=self.color_mode.get()))
            c.draw_prims(crossing_patches(self.strands, self.crossings, self.overrides,
                                          self.cfg, self.ss_colors, self.color_mode.get()))
            for s in self.strands:
                if s.id in self.sel:
                    self._draw_selection(c, s)
            for p in self.crossings:
                x, y = c.w2s(p[0] * UNIT, p[1] * UNIT)
                c.create_oval(x - 5, y - 5, x + 5, y + 5, outline="#ff8800", width=2)
            if self.draw_path and len(self.draw_path) >= 1:
                pts = [c.w2s(px * UNIT, py * UNIT) for px, py in self.draw_path]
                if len(pts) >= 2:
                    c.create_line(*[v for xy in pts for v in xy], fill="#888",
                                  width=2, dash=(4, 3))
                for px, py in pts:
                    c.create_oval(px - 3, py - 3, px + 3, py + 3, fill="#888", outline="")

        def _draw_grid(self, c):
            g = "#e7e7ec"
            for i in range(self.cols + 1):
                x0, y0 = c.w2s(i * UNIT, 0)
                x1, y1 = c.w2s(i * UNIT, self.rows * UNIT)
                c.create_line(x0, y0, x1, y1, fill=g)
            for j in range(self.rows + 1):
                x0, y0 = c.w2s(0, j * UNIT)
                x1, y1 = c.w2s(self.cols * UNIT, j * UNIT)
                c.create_line(x0, y0, x1, y1, fill=g)

        def _draw_selection(self, c, s):
            pts = [c.w2s(px * UNIT, py * UNIT) for px, py in s.path]
            if len(pts) >= 2:
                c.create_line(*[v for xy in pts for v in xy], fill="#0a84ff",
                              width=2, dash=(2, 2))
            hx, hy = pts[0]
            c.create_oval(hx - 6, hy - 6, hx + 6, hy + 6, outline="#0a84ff", width=2)
            c.create_text(hx, hy - 12, text="N", fill="#0a84ff", font=("Helvetica", 9, "bold"))

        def _press(self, e):
            if getattr(self, "_eyedrop", False):
                self._do_eyedrop(e.x, e.y)
                return
            t = self.tool.get()
            g = self.grid_at(e.x, e.y)
            if t in ("helix", "sheet", "loop"):
                self.draw_path = [g]
            elif t == "select":
                sid = self._hit_strand(e.x, e.y)
                shift = bool(e.state & 0x0001)
                if sid is None:
                    if not shift:
                        self.sel = []
                elif shift:
                    if sid in self.sel:
                        self.sel.remove(sid)
                    else:
                        self.sel.append(sid)
                else:
                    self.sel = [sid]
                self.lin_sel = None
                self._update_sel_label()
                self.render_main()
                self.render_linear()
                self._set_status()
            elif t == "selgroup":
                sid = self._hit_strand(e.x, e.y)
                if sid is not None:
                    s = self._get(sid)
                    if s and s.group:
                        self.sel = [x.id for x in self.strands if x.group == s.group]
                    else:
                        self.sel = [sid]
                    self.lin_sel = None
                    self._update_sel_label()
                    self.render_main()
                    self.render_linear()
                    self._set_status()
            elif t == "move":
                sid = self._hit_strand(e.x, e.y)
                if sid is not None:
                    # if the clicked chain is part of the current (group) selection,
                    # move every selected chain together; otherwise move just it
                    if sid in self.sel and len(self.sel) > 1:
                        moving = list(self.sel)
                    else:
                        self.sel = [sid]
                        moving = [sid]
                    self._update_sel_label()
                    self.push_undo()
                    self.move_ref = (g, [(m, [tuple(p) for p in self._get(m).path])
                                         for m in moving])
                self.render_main()
            elif t == "crossing":
                self._toggle_crossing(e.x, e.y)
            elif t == "split":
                self._split_at_click(e.x, e.y)
            elif t == "delete":
                self._delete_tile(e.x, e.y)

        def _motion(self, e):
            t = self.tool.get()
            g = self.grid_at(e.x, e.y)
            if t in ("helix", "sheet", "loop") and self.draw_path:
                self._extend_path(g)
                self.render_main()
            elif t == "move" and self.move_ref:
                start, items = self.move_ref
                dc, dr = g[0] - start[0], g[1] - start[1]
                if e.state & 0x0001:
                    if abs(dc) >= abs(dr):
                        dr = 0
                    else:
                        dc = 0
                for sid, orig in items:
                    st = self._get(sid)
                    if st:
                        st.path = [(c + dc, r + dr) for c, r in orig]
                self.render_main()

        def _release(self, e):
            t = self.tool.get()
            if t in ("helix", "sheet", "loop") and self.draw_path:
                if len(self.draw_path) >= 2:
                    self.push_undo()
                    s = Strand(self.draw_path, t, self.cur_color)
                    s.name = next_chain_name({x.name for x in self.strands})
                    self.strands.append(s)
                    self.sel = [s.id]
                    self._update_sel_label()
                self.draw_path = None
                self.recompute()
            elif t == "move" and self.move_ref:
                self.move_ref = None
                self.recompute()

        def _top_vertex_at(self, sx, sy):
            """The chain id + vertex index nearest the click, choosing the chain
            rendered on top when several overlap a cell."""
            g = self.grid_at(sx, sy)
            here = [(s, i) for s in self.strands
                    for i, v in enumerate(s.path) if v == g]
            if not here:
                return None
            if g in self.crossings:
                topsid, _ = crossing_tops(self.crossings, self.overrides)[g]
                cands = [(s, i) for (s, i) in here if s.id == topsid] or here
            else:
                cands = here
            s, i = max(cands, key=lambda si: self.strands.index(si[0]))
            return s, i

        def _split_at_click(self, sx, sy):
            hit = self._top_vertex_at(sx, sy)
            if not hit:
                return
            s, k = hit
            parts = split_at(s, k)
            if not parts:
                return
            self.push_undo()
            idx = self.strands.index(s)
            a, b = parts
            self.strands[idx:idx + 1] = [a, b]
            self._dedupe_names()
            self.sel = [a.id, b.id]
            self._update_sel_label()
            self.recompute()

        def _delete_tile(self, sx, sy):
            hit = self._top_vertex_at(sx, sy)
            if not hit:
                return
            s, k = hit
            self.push_undo()
            parts = delete_vertex(s, k)
            idx = self.strands.index(s)
            self.strands[idx:idx + 1] = parts
            self._dedupe_names()
            self.sel = [p.id for p in parts]
            self._update_sel_label()
            self.recompute()

        def _extend_path(self, g):
            path = self.draw_path
            guard = 0
            while path[-1] != g and guard < 800:
                guard += 1
                cx, cy = path[-1]
                dx, dy = g[0] - cx, g[1] - cy
                if abs(dx) >= abs(dy):
                    step = (cx + (1 if dx > 0 else -1), cy)
                else:
                    step = (cx, cy + (1 if dy > 0 else -1))
                if len(path) >= 2 and step == path[-2]:
                    path.pop()
                else:
                    path.append(step)

        def _hit_strand(self, sx, sy, thresh=12):
            wx, wy = self.canvas.s2w(sx, sy)
            best, bestd = None, thresh / self.canvas.scale_
            for s in self.strands:
                for a, b in zip(s.path, s.path[1:]):
                    d = _pt_seg_dist((wx, wy), (a[0] * UNIT, a[1] * UNIT),
                                     (b[0] * UNIT, b[1] * UNIT))
                    if d < bestd:
                        bestd, best = d, s.id
            return best

        def _toggle_crossing(self, sx, sy):
            wx, wy = self.canvas.s2w(sx, sy)
            best, bestd = None, 18 / self.canvas.scale_
            for p in self.crossings:
                d = math.hypot(wx - p[0] * UNIT, wy - p[1] * UNIT)
                if d < bestd:
                    bestd, best = d, p
            if best is None:
                return
            lst = self.crossings[best]
            entries = [(sid, vi) for sid, vi, _ in lst]
            cur = self.overrides.get(best, _default_top(lst))
            i = entries.index(cur) if cur in entries else 0
            self.push_undo()
            self.overrides[best] = entries[(i + 1) % len(entries)]
            self.render_main()

        def reverse_selected(self):
            ss = self._sel_strands()
            if not ss:
                return
            self.push_undo()
            for s in ss:
                s.path.reverse()
                s.types.reverse()
                s.colors.reverse()
                s.dirs.reverse()
            self.recompute()

        def reverse_section(self):
            if not self.lin_sel:
                messagebox.showinfo("No selection",
                                    "Drag across a stretch in the linear panel first.")
                return
            sid, a, b = self.lin_sel
            s = self._get(sid)
            if s:
                self.push_undo()
                for k in range(a, b):
                    s.dirs[k] *= -1
                self.recompute()

        def delete_selected(self):
            # linear sub-selection: split the chain and delete only that stretch
            if self.lin_sel:
                sid, a, b = self.lin_sel
                s = self._get(sid)
                if s:
                    self.push_undo()
                    parts = delete_edges(s, a, b)
                    idx = self.strands.index(s)
                    self.strands[idx:idx + 1] = parts
                    self._dedupe_names()
                    self.sel = [p.id for p in parts]
                    self.lin_sel = None
                    self._update_sel_label()
                    self.recompute()
                return
            # otherwise delete whole selected chains
            if self.sel:
                self.push_undo()
                self.strands = [s for s in self.strands if s.id not in self.sel]
                self.sel = []
                self._update_sel_label()
                self.recompute()

        def set_stretch_type(self, t):
            if self.lin_sel:
                sid, a, b = self.lin_sel
                s = self._get(sid)
                if s:
                    self.push_undo()
                    for k in range(a, b):
                        s.types[k] = t
                    self.recompute()
                    return
            ss = self._sel_strands()
            if ss:
                self.push_undo()
                for s in ss:
                    s.types = [t] * len(s.types)
                self.recompute()

        def join_ends(self):
            if len(self.sel) < 2:
                messagebox.showinfo("Join selected ends",
                                    "Select two or more chains (Shift-click) to join.")
                return
            self.push_undo()
            pool = [s for s in self.strands if s.id in self.sel]
            rest = [s for s in self.strands if s.id not in self.sel]
            conflict = False
            merged = True
            while merged:
                merged = False
                for i in range(len(pool)):
                    for j in range(len(pool)):
                        if i == j:
                            continue
                        new = _try_merge(pool[i], pool[j])
                        if new is GROUP_CONFLICT:
                            conflict = True
                            continue
                        if new:
                            new.id = pool[i].id
                            pool = [s for k, s in enumerate(pool) if k not in (i, j)]
                            pool.append(new)
                            merged = True
                            break
                    if merged:
                        break
            self.strands = rest + pool
            self.sel = [s.id for s in pool]
            self._update_sel_label()
            self.recompute()
            if conflict:
                messagebox.showinfo("", "Chains are part\nof Different Groups")

        def add_to_group(self):
            """Add the currently selected chains to an existing group (popup list)."""
            if not self.sel:
                messagebox.showinfo("Add to group", "Select one or more chains first.")
                return
            groups = sorted({s.group for s in self.strands if s.group},
                            key=lambda x: x.lower())
            if not groups:
                messagebox.showinfo("Add to group",
                                    "No groups yet — use 'Create Group…' first.")
                return
            win = tk.Toplevel(self)
            win.title("Add to group")
            win.transient(self)
            ttk.Label(win, text="Add selected chain(s) to:",
                      style="Head.TLabel").pack(padx=12, pady=(10, 6))

            def choose(g):
                self.push_undo()
                for s in self._sel_strands():
                    s.group = g
                win.destroy()
                self.recompute()
            for g in groups:
                ttk.Button(win, text=g, command=lambda g=g: choose(g)).pack(
                    fill="x", padx=12, pady=2)
            ttk.Button(win, text="Cancel", command=win.destroy).pack(pady=(6, 10))
            win.update_idletasks()
            self._center_window(win)

        def render_linear(self):
            c = self.lin
            c.delete("all")
            rowh, cellw, labx = 34, 17, 8
            indent = 30                   # deeper indent for grouped chains
            y = 10
            self._lin_layout = []
            self._lin_headers = []        # (group_name, x0, x1, y0, y1, kind)

            # order: grouped chains first (by group name), then ungrouped
            grouped = {}
            for s in self.strands:
                grouped.setdefault(s.group, []).append(s)
            order = sorted((g for g in grouped if g is not None),
                           key=lambda x: x.lower())
            sections = [(g, grouped[g]) for g in order]
            if None in grouped:
                sections.append((None, grouped[None]))

            for gi, (gname, members) in enumerate(sections):
                if gname is not None:
                    sel_all = members and all(s.id in self.sel for s in members)
                    bx0, by0 = c.w2s(labx, y)
                    # checkbox to select the whole group
                    c.create_rectangle(bx0, by0 + 2, bx0 + 16, by0 + 18,
                                       outline="#0a6", width=2,
                                       fill="#0a6" if sel_all else "#fff")
                    if sel_all:
                        c.create_line(bx0 + 3, by0 + 10, bx0 + 7, by0 + 15,
                                      fill="#fff", width=2)
                        c.create_line(bx0 + 7, by0 + 15, bx0 + 13, by0 + 4,
                                      fill="#fff", width=2)
                    self._lin_headers.append((gname, labx, labx + 16, y, y + 20, "check"))
                    # larger group name, click to rename
                    nx, ny = c.w2s(labx + 22, y)
                    c.create_text(nx, ny + 9, anchor="w", text=gname,
                                  font=("Helvetica", 13, "bold"), fill="#0a6")
                    nlen = max(60, len(gname) * 9)
                    self._lin_headers.append((gname, labx + 22, labx + 22 + nlen,
                                              y, y + 20, "rename"))
                    y += 24
                for s in members:
                    ind = indent if gname else 0
                    ne = len(s.types)
                    x0, y0 = c.w2s(labx + ind, y)
                    namecol = "#0a84ff" if s.id in self.sel else "#333"
                    c.create_text(x0, y0 + rowh / 2, anchor="w", text=s.name,
                                  font=("Helvetica", 9, "bold"), fill=namecol)
                    bx = labx + ind + 70
                    for k in range(ne):
                        col = (self.ss_colors[s.types[k]]
                               if self.color_mode.get() == "ss" else s.colors[k])
                        if self.cfg.get("loops_black") and s.types[k] == "loop":
                            col = "#000000"
                        sx0, sy0 = c.w2s(bx + k * cellw, y + 6)
                        sx1, sy1 = c.w2s(bx + (k + 1) * cellw, y + rowh - 6)
                        sel = (self.lin_sel and self.lin_sel[0] == s.id
                               and self.lin_sel[1] <= k < self.lin_sel[2])
                        c.create_rectangle(sx0, sy0, sx1, sy1, fill=col,
                                           outline="#0a84ff" if sel else "#333",
                                           width=2 if sel else 1)
                        midy = (sy0 + sy1) / 2
                        fwd = s.dirs[k] > 0
                        if s.types[k] == "sheet":
                            c.create_line(sx0 + 2, midy, sx1 - 2, midy, fill="#fff",
                                          arrow=("last" if fwd else "first"))
                        elif s.types[k] == "helix":
                            c.create_line(sx0 + 2, midy, sx1 - 2, midy, fill="#fff", width=3)
                    self._lin_layout.append((s.id, bx, y, ne, cellw, rowh))
                    y += rowh + 6
                y += 14                   # extra gap between groups

        def _lin_hit(self, sx, sy):
            wx, wy = self.lin.s2w(sx, sy)
            for (sid, bx, y, ne, cellw, rowh) in getattr(self, "_lin_layout", []):
                if y <= wy <= y + rowh:
                    k = int((wx - bx) / cellw)
                    if 0 <= k < ne:
                        return (sid, k)
            return None

        def _lin_press(self, e):
            wx, wy = self.lin.s2w(e.x, e.y)
            for (gname, x0, x1, y0, y1, kind) in getattr(self, "_lin_headers", []):
                if x0 <= wx <= x1 and y0 <= wy <= y1:
                    if kind == "check":
                        self.select_group_named(gname)
                    else:
                        self.rename_group(gname)
                    return
            hit = self._lin_hit(e.x, e.y)
            if hit:
                sid, k = hit
                self.lin_sel = (sid, k, k + 1)
                self._lin_drag = (sid, k)
                self.sel = [sid]
                self._update_sel_label()
                self.render_linear()
                self.render_main()

        def _lin_motion(self, e):
            if not self._lin_drag:
                return
            hit = self._lin_hit(e.x, e.y)
            if hit and hit[0] == self._lin_drag[0]:
                sid, k0 = self._lin_drag
                k = hit[1]
                self.lin_sel = (sid, min(k0, k), max(k0, k) + 1)
                self.render_linear()

        def _lin_release(self, e):
            self._lin_drag = None

        def apply_stretch_color(self):
            self.apply_color()

        def resize_grid(self):
            self.cols = max(4, int(self.col_var.get()))
            self.rows = max(4, int(self.row_var.get()))
            self.render_main()
            self._set_status()

        def _maybe_save(self):
            """Offer to save unsaved work. Returns True to proceed, False to abort."""
            if not self.dirty or not self.strands:
                return True
            ans = messagebox.askyesnocancel("Save changes?",
                                            "Save changes to the current diagram?")
            if ans is None:
                return False
            if ans:
                return self.save_project()
            return True

        def new_project(self):
            if not self._maybe_save():
                return
            self.strands = []
            self.overrides = {}
            self.sel = []
            self.lin_sel = None
            self.undo_stack = []
            self.cur_file = None
            self.dirty = False
            self._update_sel_label()
            self.recompute()

        def close_file(self):
            """Close the current file (prompt to save), leaving an empty canvas."""
            self.new_project()

        def _write(self, f):
            data = {"cols": self.cols, "rows": self.rows, "ss_colors": self.ss_colors,
                    "color_mode": self.color_mode.get(),
                    "loops_black": bool(self.loops_black.get()),
                    "cfg": {k: self.cfg[k] for k in
                            ("sheet_width", "helix_width", "loop_width",
                             "sheet_outline_w", "helix_outline_w")},
                    "strands": [s.to_dict() for s in self.strands],
                    "overrides": {"%d,%d" % p: list(t) for p, t in self.overrides.items()}}
            with open(f, "w") as fh:
                json.dump(data, fh, indent=1)
            self.cur_file = f
            self.dirty = False
            self._set_status("saved %s" % f)
            return True

        def save_project(self):
            if self.cur_file:
                return self._write(self.cur_file)
            return self.save_project_as()

        def save_project_as(self):
            f = filedialog.asksaveasfilename(defaultextension=".topo",
                                             filetypes=[("Flat Strandley", "*.topo"),
                                                        ("JSON", "*.json"), ("All", "*")])
            if not f:
                return False
            return self._write(f)

        def open_project(self):
            if not self._maybe_save():
                return
            f = filedialog.askopenfilename(filetypes=[("Flat Strandley", "*.topo"),
                                                      ("JSON", "*.json"), ("All", "*")])
            if not f:
                return
            with open(f) as fh:
                data = json.load(fh)
            self.cols, self.rows = data.get("cols", 30), data.get("rows", 22)
            self.col_var.set(self.cols)
            self.row_var.set(self.rows)
            self.ss_colors = data.get("ss_colors", dict(SS_DEFAULTS))
            self.color_mode.set(data.get("color_mode", "custom"))
            self.loops_black.set(bool(data.get("loops_black", False)))
            self.cfg.update(data.get("cfg", {}))
            self.cfg["loops_black"] = bool(data.get("loops_black", False))
            Strand._next = 1
            self.strands = [Strand.from_dict(d) for d in data.get("strands", [])]
            self.overrides = {}
            for key, val in data.get("overrides", {}).items():
                cx, cy = key.split(",")
                self.overrides[(int(cx), int(cy))] = tuple(val)
            for t in self.ss_btns:
                self.ss_btns[t].configure(bg=self.ss_colors[t])
            self.sel = []
            self.lin_sel = None
            self.undo_stack = []
            self.cur_file = f
            self.dirty = False
            self._update_sel_label()
            self.recompute()
            self.fit_view()

        def on_close(self):
            if self._maybe_save():
                self.destroy()

        def export_svg(self):
            f = filedialog.asksaveasfilename(defaultextension=".svg",
                                             filetypes=[("SVG", "*.svg"), ("All", "*")])
            if not f:
                return
            svg = project_to_svg(self.strands, self.crossings, self.overrides,
                                 self.cfg, self.ss_colors, self.color_mode.get())
            with open(f, "w") as fh:
                fh.write(svg)
            self._set_status("exported SVG (transparent background).")


def main():
    if not HAS_TK:
        raise SystemExit("Tkinter is required to run the GUI.")
    App().mainloop()


if __name__ == "__main__":
    main()
