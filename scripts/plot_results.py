#!/usr/bin/env python3
"""Abbildung der Verwandtschaftsanalyse aus mendel.tsv und king.kin0.

Die Rollen (Eltern, Kind, nicht verwandt) werden NICHT vorgegeben, sondern aus den
Daten abgeleitet:
  - Paare: Trennung an der größten Lücke der Mendel-Fehlerraten (log-Skala)
  - Trio:  das Trio mit der kleinsten Fehlerrate liefert Kind und Eltern
"""
import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

# ---------------------------------------------------------------- Farben ----
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
REL = "#2a78d6"      # Eltern-Kind (hervorgehoben)
UNREL = "#898781"    # nicht verwandt (zurückgenommen)
REL_TINT = "#cde2fb"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Inter", "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": INK,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": INK2,
    "xtick.labelcolor": INK2,
    "pdf.fonttype": 42,
})


# ------------------------------------------------------ Zahlenformatierung ----
def de(x, nd=2):
    return f"{x:.{nd}f}".replace(".", ",")


def pct(frac):
    """Anteil -> '0,015 %' mit sinnvoller Stellenzahl (deutsches Komma)."""
    v = 100 * frac
    nd = 3 if v < 0.1 else 2 if v < 1 else 1
    return de(v, nd) + "\u00a0%"


def thousands(n):
    return f"{n:,}".replace(",", ".")


# ------------------------------------------------------------ Einlesen ----
def read_mendel(path):
    pairs, trios, section, header = [], [], None, None
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("##"):
            section = "pairs" if "Paar" in s else "trios"
            header = None
            continue
        cols = s.split("\t")
        if header is None:
            header = cols
            continue
        r = dict(zip(header, cols))
        sites, err = int(r["sites"]), int(r["mendel_errors"])
        if section == "pairs":
            pairs.append(dict(a=r["probe1"], b=r["probe2"], sites=sites, err=err, rate=err / sites))
        else:
            trios.append(dict(p1=r["elter1"], p2=r["elter2"], child=r["kind"],
                              sites=sites, err=err, rate=err / sites))
    return pairs, trios


def read_king(path):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip()]
    header = [h.lstrip("#") for h in rows[0]]
    out = []
    for r in rows[1:]:
        d = dict(zip(header, r))
        out.append(dict(a=d["IID1"], b=d["IID2"], nsnp=int(d["NSNP"]),
                        ibs0=float(d["IBS0"]), kin=float(d["KINSHIP"])))
    return out


def pair_key(a, b):
    return frozenset((a, b))


def pair_label(a, b):
    return " – ".join(sorted((a, b)))


# ----------------------------------------------- Rollen aus den Daten ----
def split_pairs(pairs, min_ratio=5.0):
    """Schwelle an der größten Lücke (log) zwischen den sortierten Fehlerraten."""
    rates = sorted(max(p["rate"], 1e-7) for p in pairs)
    gaps = [(math.log(rates[i + 1] / rates[i]), i) for i in range(len(rates) - 1)]
    g, i = max(gaps)
    if g < math.log(min_ratio):
        return None
    return math.sqrt(rates[i] * rates[i + 1])


# ----------------------------------------------------------- Gestaltung ----
def style_axis(ax, grid="x"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(length=3, width=0.8, labelsize=8.5)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.7)


def panel_title(fig, ax, letter, text, x):
    """Titel als Figure-Text: x in Figure-Koordinaten, y knapp über dem Achsenkasten."""
    y = ax.get_position(original=True).y1 + 0.016
    fig.text(x, y, f"{letter}  {text}", fontsize=10.5, fontweight="bold",
             color=INK, ha="left", va="bottom")


def ring(color):
    return dict(s=78, color=color, edgecolor=SURFACE, linewidth=1.6, zorder=3)


# -------------------------------------------------------- Panel A: Netz ----
def box_edge_point(c, target, w, h, pad=0.06):
    """Punkt, an dem die Linie von c nach target den Kasten (Breite w, Höhe h) verlässt."""
    dx, dy = target[0] - c[0], target[1] - c[1]
    if dx == 0 and dy == 0:
        return c
    tx = (w / 2 + pad) / abs(dx) if dx else float("inf")
    ty = (h / 2 + pad) / abs(dy) if dy else float("inf")
    t = min(tx, ty)
    return (c[0] + dx * t, c[1] + dy * t)


def draw_network(ax, names, pairs, trio, cut):
    W, H = 1.62, 0.64
    ax.set_xlim(-2.2, 4.1)
    ax.set_ylim(-1.95, 1.45)
    ax.set_aspect("equal", adjustable="box")
    ax.set_anchor("N")
    ax.axis("off")

    related = lambda p: cut is not None and p["rate"] < cut
    rate = {pair_key(p["a"], p["b"]): p["rate"] for p in pairs}

    # Positionen: Trio als Stammbaum-Dreieck, übrige Personen rechts
    if trio:
        pos = {trio["p1"]: (-1.3, 1.0), trio["p2"]: (1.3, 1.0), trio["child"]: (0.0, -1.0)}
        rest = [n for n in names if n not in pos]
        for i, n in enumerate(rest):
            pos[n] = (3.3, 0.0 if len(rest) == 1 else 1.0 - 2.0 * i / (len(rest) - 1))
        role = {trio["p1"]: "Elternteil", trio["p2"]: "Elternteil", trio["child"]: "Kind"}
        for n in rest:
            role[n] = "nicht verwandt"
    else:
        pos, role = {}, {n: "" for n in names}
        for i, n in enumerate(names):
            ang = 2 * math.pi * i / len(names)
            pos[n] = (1.8 * math.cos(ang), 1.2 * math.sin(ang))

    # 1) nicht verwandte Paare: dünne, zurückgenommene Linien
    for p in pairs:
        if related(p):
            continue
        a, b = pos[p["a"]], pos[p["b"]]
        pa, pb = box_edge_point(a, b, W, H), box_edge_point(b, a, W, H)
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=AXIS, lw=1.2, zorder=1,
                solid_capstyle="round")

    # 2) Eltern-Kind-Paare: dicke blaue Pfeile Richtung Kind (Richtung aus dem Trio)
    for p in pairs:
        if not related(p):
            continue
        a, b = p["a"], p["b"]
        if trio and b == trio["child"]:
            src, dst = a, b
        elif trio and a == trio["child"]:
            src, dst = b, a
        else:
            src, dst = a, b
        s, d = pos[src], pos[dst]
        ps, pd = box_edge_point(s, d, W, H), box_edge_point(d, s, W, H)
        ax.annotate("", xy=pd, xytext=ps, zorder=2,
                    arrowprops=dict(arrowstyle="-|>,head_length=0.7,head_width=0.32",
                                    color=REL, lw=2.8, shrinkA=0, shrinkB=0))
        mid = ((ps[0] + pd[0]) / 2, (ps[1] + pd[1]) / 2)
        e = next(q for q in pairs if pair_key(q["a"], q["b"]) == pair_key(a, b))
        side_left = (s[0] < d[0]) if s[0] != d[0] else True
        ax.text(mid[0] + (-0.12 if side_left else 0.12), mid[1], f"Fehler {pct(e['rate'])}",
                ha="right" if side_left else "left", va="center", fontsize=8.5,
                color=INK, zorder=4,
                bbox=dict(boxstyle="round,pad=0.18", fc=SURFACE, ec="none"))

    # 3) Knoten
    for n in names:
        x, y = pos[n]
        is_child = trio and n == trio["child"]
        is_parent = trio and n in (trio["p1"], trio["p2"])
        ec = REL if (is_child or is_parent) else UNREL
        fc = REL_TINT if is_child else SURFACE
        ax.add_patch(FancyBboxPatch((x - W / 2, y - H / 2), W, H,
                                    boxstyle="round,pad=0,rounding_size=0.14",
                                    fc=fc, ec=ec, lw=2.0 if is_child else 1.4, zorder=5))
        ax.text(x, y + 0.08, n, ha="center", va="center", fontsize=9.5,
                fontweight="bold", color=INK, zorder=6)
        ax.text(x, y - 0.17, role[n], ha="center", va="center", fontsize=7.8,
                color=INK2, zorder=6)

    # 4) Fehlerrate der nicht verwandten Paare als eine Zeile je Gruppe
    unrel_rates = [p["rate"] for p in pairs if not related(p)]
    if unrel_rates:
        lo, hi = min(unrel_rates), max(unrel_rates)
        txt = f"graue Linien: nicht verwandt, Fehler {pct(lo)} bis {pct(hi)}"
        ax.text(0.5 * (ax.get_xlim()[0] + ax.get_xlim()[1]) + 0.0, -1.78, txt,
                ha="center", va="center", fontsize=8.3, color=INK2)


# ------------------------------------------------- Panel B: Paare (log) ----
def draw_pairs(ax, pairs, cut):
    rows = sorted(pairs, key=lambda p: p["rate"])
    n = len(rows)
    ys = list(range(n))
    for y, p in zip(ys, rows):
        rel = cut is not None and p["rate"] < cut
        ax.scatter(100 * p["rate"], y, **ring(REL if rel else UNREL))
        ax.text(100 * p["rate"] * 1.35, y, pct(p["rate"]), va="center", ha="left",
                fontsize=8.5, color=INK, fontweight="bold" if rel else "normal")
    ax.set_yticks(ys)
    ax.set_yticklabels([pair_label(p["a"], p["b"]) for p in rows], fontsize=8.8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.008, 40)
    ax.xaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: de(v, 2 if v < 0.1 else 1 if v < 1 else 0) + "\u00a0%"))
    ax.set_xlabel("Mendel-Fehler (Anteil der SNPs, log-Skala)", fontsize=8.8)
    ax.set_ylim(n - 0.5, -0.6)
    style_axis(ax, grid="x")
    ax.tick_params(axis="y", length=0)

    if cut is not None:
        k = sum(1 for p in rows if p["rate"] < cut)
        if 0 < k < n:
            a, b = 100 * rows[k - 1]["rate"], 100 * rows[k]["rate"]
            ymid = k - 0.5
            ax.annotate("", xy=(b, ymid), xytext=(a, ymid),
                        arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.1,
                                        shrinkA=0, shrinkB=0))
            ax.text(math.sqrt(a * b), ymid + 0.1, f"{math.floor(b / a)}-facher Abstand",
                    ha="center", va="top", fontsize=8.5, color=INK2)


# ------------------------------------------------- Panel C: Trios ----
def draw_trios(ax, trios):
    rows = sorted(trios, key=lambda t: t["rate"])
    n = len(rows)
    best = rows[0]
    for y, t in enumerate(rows):
        is_best = y == 0
        ax.scatter(100 * t["rate"], y, **ring(REL if is_best else UNREL))
    ax.axhspan(-0.5, 0.5, color=REL_TINT, alpha=0.35, lw=0, zorder=0)
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{t['p1']} + {t['p2']} → {t['child']}" for t in rows], fontsize=8.2)
    ax.invert_yaxis()
    ax.set_ylim(n - 0.5, -0.6)
    xmax = 100 * max(t["rate"] for t in rows)
    ax.set_xlim(0, math.ceil(xmax / 5.0) * 5 + 2)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v)}\u00a0%"))
    ax.set_xlabel("Mendel-Fehler im Trio (Anteil der SNPs)", fontsize=8.8)
    style_axis(ax, grid="x")
    ax.tick_params(axis="y", length=0)

    ax.text(100 * best["rate"] + 1.0, 0, pct(best["rate"]), va="center", ha="left",
            fontsize=8.5, fontweight="bold", color=INK)
    if n > 1:
        nxt = rows[1]
        ax.text(100 * nxt["rate"] + 1.0, 1, f"nächstbeste Kombination: {pct(nxt['rate'])}",
                va="center", ha="left", fontsize=8.5, color=INK2)
        ratio = nxt["rate"] / best["rate"]
        ax.annotate("", xy=(100 * nxt["rate"], 0.5), xytext=(100 * best["rate"], 0.5),
                    arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.1, shrinkA=0, shrinkB=0))
        ax.text(50 * (best["rate"] + nxt["rate"]), 0.58, f"{math.floor(ratio)}-facher Abstand",
                ha="center", va="top", fontsize=8.5, color=INK2)


# ------------------------------------------------- Panel D: KING ----
def draw_king(ax, king, pairs, cut):
    rate = {pair_key(p["a"], p["b"]): p["rate"] for p in pairs}
    pts = []
    for k in king:
        r = rate.get(pair_key(k["a"], k["b"]))
        rel = cut is not None and r is not None and r < cut
        pts.append((k, rel))

    xmax = max(k["ibs0"] for k, _ in pts)
    ymax = max(k["kin"] for k, _ in pts)
    ax.set_xlim(-0.004, xmax * 1.28)
    ax.set_ylim(-0.03, max(0.36, ymax * 1.12))

    # theoretische Erwartungswerte als ruhige Referenzlinien
    for yv, txt in ((0.25, "theoretisch Eltern–Kind: 0,25"), (0.0, "theoretisch nicht verwandt: 0")):
        ax.axhline(yv, color=AXIS, lw=0.9, zorder=1)
        ax.text(ax.get_xlim()[1], yv + 0.008, txt, ha="right", va="bottom", fontsize=8,
                color=MUTED)

    for k, rel in pts:
        ax.scatter(k["ibs0"], k["kin"], **ring(REL if rel else UNREL))

    # direkte Beschriftung mit Leitlinie: Eltern-Kind-Punkte liegen dicht beieinander
    rel_pts = sorted([k for k, rel in pts if rel], key=lambda k: k["ibs0"])
    offsets = [(46, -22), (46, 14)]          # links liegender Punkt -> Label unten, rechter -> oben
    for k, off in zip(rel_pts, offsets):
        ax.annotate(pair_label(k["a"], k["b"]), (k["ibs0"], k["kin"]),
                    xytext=off, textcoords="offset points", fontsize=8.5, color=INK,
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9, shrinkA=2, shrinkB=5))
    unrel_pts = [k for k, rel in pts if not rel]
    if unrel_pts:
        cx = sum(k["ibs0"] for k in unrel_pts) / len(unrel_pts)
        cy = sum(k["kin"] for k in unrel_pts) / len(unrel_pts)
        ax.annotate(f"{len(unrel_pts)} nicht verwandte Paare", (cx, cy),
                    xytext=(0, 24), textcoords="offset points", fontsize=8.5, color=INK2,
                    ha="center", va="bottom",
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9, shrinkA=2, shrinkB=9))

    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: de(v, 2)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: de(v, 2)))
    ax.set_xlabel("IBS0 (Anteil entgegengesetzter Homozygoter)", fontsize=8.8)
    ax.set_ylabel("KINSHIP (KING-robust)", fontsize=8.8)
    style_axis(ax, grid=None)


# ------------------------------------------------------------- Haupt ----
def save(fig, base):
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(base.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print("geschrieben:", base.with_suffix(".png"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mendel", default="results/mendel.tsv")
    ap.add_argument("--king", default="results/king.kin0")
    ap.add_argument("--outdir", default="results/abbildungen")
    a = ap.parse_args()

    pairs, trios = read_mendel(a.mendel)
    king = read_king(a.king)
    names = sorted({p["a"] for p in pairs} | {p["b"] for p in pairs})
    cut = split_pairs(pairs)
    trio = min(trios, key=lambda t: t["rate"]) if trios else None
    out = Path(a.outdir)

    # Keine Titel im Bild: Beschriftung kommt in die Abbildungsunterschrift
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    draw_network(ax, names, pairs, trio, cut)
    ax.set_xlim(-2.2, 4.35)
    save(fig, out / "abb1_verwandtschaftsnetz")

    fig, ax = plt.subplots(figsize=(7.0, 3.4), layout="constrained")
    draw_pairs(ax, pairs, cut)
    save(fig, out / "abb2_mendel_paare")

    if trios:
        fig, ax = plt.subplots(figsize=(7.4, 4.4), layout="constrained")
        draw_trios(ax, trios)
        save(fig, out / "abb3_mendel_trios")

    fig, ax = plt.subplots(figsize=(6.6, 4.6), layout="constrained")
    draw_king(ax, king, pairs, cut)
    save(fig, out / "abb4_king")


if __name__ == "__main__":
    main()
