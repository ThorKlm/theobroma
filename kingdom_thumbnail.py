"""Server-side SVG generation for the /browse kingdom-thumbnail.
Strict-proportional radial chart, no minimum-visible correction.
Sub-percent kingdoms remain visually invisible but listed in the legend.
"""
import math

# Kingdom display order matching the home page chip row.
KINGDOM_ORDER = ["plant", "fungi", "bacteria", "animal", "multi", "unresolved"]

# Fallback colour set if CSS variable resolution fails (renders inline-only SVG).
KINGDOM_FALLBACK = {
    "plant":      "#2d6a4f",
    "fungi":      "#9d4edd",
    "bacteria":   "#0077b6",
    "animal":     "#bc4749",
    "multi":      "#7d8597",
    "unresolved": "#aaaaaa",
}

def _polar(cx, cy, r, angle_rad):
    return (cx + r * math.cos(angle_rad), cy + r * math.sin(angle_rad))

def _arc_path(cx, cy, r_outer, r_inner, start_rad, end_rad):
    """Donut-segment path. start_rad and end_rad in radians, measured clockwise from -pi/2 (12 o'clock)."""
    large_arc = 1 if (end_rad - start_rad) > math.pi else 0
    x1, y1 = _polar(cx, cy, r_outer, start_rad)
    x2, y2 = _polar(cx, cy, r_outer, end_rad)
    x3, y3 = _polar(cx, cy, r_inner, end_rad)
    x4, y4 = _polar(cx, cy, r_inner, start_rad)
    return (f"M {x1:.2f} {y1:.2f} "
            f"A {r_outer} {r_outer} 0 {large_arc} 1 {x2:.2f} {y2:.2f} "
            f"L {x3:.2f} {y3:.2f} "
            f"A {r_inner} {r_inner} 0 {large_arc} 0 {x4:.2f} {y4:.2f} "
            f"Z")

def kingdom_thumbnail_svg(kingdom_counts, total=None, size=200, title="Kingdoms"):
    """Build the donut SVG and a legend.
    kingdom_counts: iterable of dicts with keys 'kingdom' and 'cnt', or (kingdom, count) tuples.
    Returns dict with keys 'svg', 'legend' (a list of (kingdom, count, pct) for the legend rendering).
    """
    if not kingdom_counts: return {"svg": "", "legend": []}
    # Normalise input to list of (kingdom, count)
    rows = []
    for k in kingdom_counts:
        if isinstance(k, dict):
            rows.append((k.get("kingdom") or "unknown", int(k.get("cnt") or 0)))
        else:
            rows.append((k[0] or "unknown", int(k[1] or 0)))
    if total is None: total = sum(c for _, c in rows)
    if total == 0: return {"svg": "", "legend": []}
    # Order by KINGDOM_ORDER, then by count desc for unknown
    rows_by_kingdom = {k: c for k, c in rows}
    ordered = []
    for k in KINGDOM_ORDER:
        if k in rows_by_kingdom and rows_by_kingdom[k] > 0:
            ordered.append((k, rows_by_kingdom[k]))
    for k, c in rows:
        if k not in KINGDOM_ORDER and c > 0:
            ordered.append((k, c))
    # Render donut
    cx = cy = size / 2
    r_outer = size * 0.42
    r_inner = size * 0.24
    parts = [f'<svg viewBox="0 0 {size} {size}" width="{size}" height="{size}" xmlns="http://www.w3.org/2000/svg" class="kingdom-thumbnail" style="flex-shrink:0;">']
    parts.append(f'<title>{title}: {total:,} compounds</title>')
    angle = -math.pi / 2  # start at top (12 o'clock)
    legend = []
    for kingdom, count in ordered:
        pct = 100.0 * count / total
        sweep = 2 * math.pi * count / total
        if sweep < 1e-6:
            legend.append((kingdom, count, pct))
            continue
        end = angle + sweep
        colour = KINGDOM_FALLBACK.get(kingdom, "#999")
        # Full-circle case: SVG path with identical endpoints renders as zero-length.
        # Split into two semicircles or use circle-with-hole.
        if abs(sweep - 2 * math.pi) < 1e-6:
            parts.append(f'<circle cx="{cx}" cy="{cy}" r="{r_outer}" fill="{colour}" stroke="white" stroke-width="1"><title>{kingdom}: {count:,} ({pct:.1f}%)</title></circle>')
            parts.append(f'<circle cx="{cx}" cy="{cy}" r="{r_inner}" fill="white" />')
        else:
            path = _arc_path(cx, cy, r_outer, r_inner, angle, end)
            parts.append(f'<path d="{path}" fill="{colour}" stroke="white" stroke-width="1"><title>{kingdom}: {count:,} ({pct:.1f}%)</title></path>')
        legend.append((kingdom, count, pct))
        angle = end
    # Total in the donut hole
    parts.append(f'<text x="{cx}" y="{cy}" text-anchor="middle" dominant-baseline="central" font-size="{size*0.10:.0f}" fill="#333" font-weight="600">{total:,}</text>')
    parts.append(f'<text x="{cx}" y="{cy + size*0.10:.1f}" text-anchor="middle" font-size="{size*0.06:.0f}" fill="#666">total</text>')
    parts.append('</svg>')
    return {"svg": "".join(parts), "legend": legend}
