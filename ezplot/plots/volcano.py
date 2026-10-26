"""Volcano plot in the format of the template deck (slide 3, "Differential gene expression").

x = log2FC (oriented: positive = higher in treatment), y = -log10(adjusted P). Genes with
FDR < cutoff and |log2FC| >= threshold are coloured; the rest are grey. Significant genes
are labelled (italic) with leader lines, placed by adjustText.
"""

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

from .. import style


@dataclass
class VolcanoParams:
    title: str = "Differential gene expression"
    fdr: float = 0.05
    lfc: float = 1.0
    label_mode: str = "significant"     # significant / top / custom / none
    label_top: int = 20
    label_rank: str = "fdr"             # top N by: fdr / lfc / distance
    label_genes: list = field(default_factory=list)
    neg_color: str = style.NEG
    pos_color: str = style.POS
    ns_color: str = style.NS
    neg_label: str = "Negative log$_2$FC"
    pos_label: str = "Positive log$_2$FC"
    show_thresholds: bool = True
    xlim: float = 0.0                   # 0 = automatic, symmetric
    ymax: float = 0.0                   # 0 = automatic
    point_size: float = 70
    label_fontsize: float = 13
    legend_fontsize: float = 12.5
    width: float = 8.1
    height: float = 7.0


def classify(df, p):
    sig = (df["padj"] < p.fdr) & (df["lfc"].abs() >= p.lfc)
    return np.where(sig & (df["lfc"] > 0), "pos", np.where(sig & (df["lfc"] < 0), "neg", "ns"))


def plot_volcano(genes, p=None):
    """genes: DataFrame with columns symbol, lfc (oriented), padj."""
    p = p or VolcanoParams()
    style.setup()
    df = genes[["symbol", "lfc", "padj"]].dropna().copy()
    floor = df.loc[df["padj"] > 0, "padj"].min() if (df["padj"] > 0).any() else 1e-300
    df["y"] = -np.log10(df["padj"].clip(lower=floor / 10))
    df["cls"] = classify(df, p)
    counts = df["cls"].value_counts()
    lab = _labels(df, p)

    fig = plt.figure(figsize=(p.width, p.height))
    ax = fig.add_axes([0.09, 0.10, 0.88, 0.66])

    ns = df[df["cls"] == "ns"]
    ax.scatter(ns["lfc"], ns["y"], s=p.point_size * 0.35, color=p.ns_color, alpha=0.8,
               edgecolors="none", zorder=2, rasterized=len(ns) > 3000)
    for cls, color in (("neg", p.neg_color), ("pos", p.pos_color)):
        d = df[df["cls"] == cls]
        ax.scatter(d["lfc"], d["y"], s=p.point_size, color=color, edgecolors="white",
                   linewidths=0.8, zorder=3)

    n_lab = len(lab)
    x_pad = 1.05 + (0.002 * min(n_lab, 80) if n_lab > 15 else 0.0)
    y_pad = 1.08 + (0.0025 * min(n_lab, 80) if n_lab > 15 else 0.0)
    xl = p.xlim or np.ceil(df["lfc"].abs().max() * x_pad * 2) / 2
    ax.set_xlim(-xl, xl)
    ytop = p.ymax or df["y"].max() * y_pad
    ax.set_ylim(0, max(ytop, -np.log10(p.fdr) * 1.5))

    fdr_text = None
    if p.show_thresholds:
        kw = dict(color="#888888", ls=(0, (3, 2)), lw=1.3, zorder=1)
        ax.axhline(-np.log10(p.fdr), **kw)
        if p.lfc > 0:
            ax.axvline(-p.lfc, **kw)
            ax.axvline(p.lfc, **kw)
        fdr_text = ax.text(0.995, 0.012, f"FDR < {p.fdr:g}; |log$_2$FC| ≥ {p.lfc:g}",
                           transform=ax.transAxes, ha="right", va="bottom", fontsize=13,
                           color="#6b7078", bbox=dict(fc="white", ec="none", pad=1, alpha=0.85))

    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.4)
    ax.tick_params(labelsize=14, width=1.4, length=6)
    ax.set_xlabel("log$_2$ fold change", fontsize=18, labelpad=8)
    ax.set_ylabel("$-$log$_{10}$ (adjusted $P$ value)", fontsize=18, labelpad=8)

    fig.text(0.09, 0.975, p.title, ha="left", va="top", fontsize=22, fontweight="bold",
             color=style.TEXT)
    handles = [
        Line2D([], [], ls="", marker="o", ms=9, color=p.neg_color,
               label=f"{p.neg_label} (n={counts.get('neg', 0):,})"),
        Line2D([], [], ls="", marker="o", ms=9, color=p.ns_color,
               label=f"Not significant (n={counts.get('ns', 0):,})"),
        Line2D([], [], ls="", marker="o", ms=9, color=p.pos_color,
               label=f"{p.pos_label} (n={counts.get('pos', 0):,})"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.02, 0.84, 0.97, 0.05),
               mode="expand", ncol=3, frameon=False, fontsize=p.legend_fontsize,
               handletextpad=0.3, columnspacing=0.8, borderaxespad=0)

    if len(lab):
        texts = [
            ax.text(r["lfc"], r["y"], r["symbol"], fontsize=p.label_fontsize,
                    fontstyle="italic", color=style.TEXT, zorder=5,
                    ha="center", va="center")
            for _, r in lab.iterrows()
        ]
        _place_labels(ax, texts, lab["lfc"].to_numpy(), lab["y"].to_numpy(), df, p, fdr_text)
    return fig


def _segments_intersect(p1, p2, q1, q2):
    """Vectorized check if segments p1->p2 (K, 2) intersect q1->q2 (M, 2)."""
    p1 = p1[:, None, :]
    p2 = p2[:, None, :]
    q1 = q1[None, :, :]
    q2 = q2[None, :, :]

    def ccw(a, b, c):
        return (c[..., 1] - a[..., 1]) * (b[..., 0] - a[..., 0]) > (
            b[..., 1] - a[..., 1]
        ) * (c[..., 0] - a[..., 0])

    return (ccw(p1, q1, q2) != ccw(p2, q1, q2)) & (ccw(p1, p2, q1) != ccw(p1, p2, q2))


def _seg_intersects_box(p1, p2, bx0, bx1, by0, by1):
    """Vectorized Liang-Barsky check if segments p1->p2 (K, 2) intersect boxes (M,)."""
    dx = (p2[:, 0] - p1[:, 0])[:, None]
    dy = (p2[:, 1] - p1[:, 1])[:, None]
    x0 = p1[:, 0][:, None]
    y0 = p1[:, 1][:, None]

    seg_xmin = np.minimum(p1[:, 0], p2[:, 0])[:, None]
    seg_xmax = np.maximum(p1[:, 0], p2[:, 0])[:, None]
    seg_ymin = np.minimum(p1[:, 1], p2[:, 1])[:, None]
    seg_ymax = np.maximum(p1[:, 1], p2[:, 1])[:, None]
    aabb = (
        (seg_xmax >= bx0[None, :])
        & (seg_xmin <= bx1[None, :])
        & (seg_ymax >= by0[None, :])
        & (seg_ymin <= by1[None, :])
    )
    if not np.any(aabb):
        return aabb

    t0 = np.zeros_like(aabb, dtype=float)
    t1 = np.ones_like(aabb, dtype=float)

    for p, q in (
        (-dx, x0 - bx0[None, :]),
        (dx, bx1[None, :] - x0),
        (-dy, y0 - by0[None, :]),
        (dy, by1[None, :] - y0),
    ):
        parallel = np.abs(p) < 1e-9
        outside = parallel & (q < 0)
        t0 = np.where(outside, 2.0, t0)
        safe_p = np.where(parallel, 1.0, p)
        r = q / safe_p
        t0 = np.where(~parallel & (p < 0), np.maximum(t0, r), t0)
        t1 = np.where(~parallel & (p > 0), np.minimum(t1, r), t1)

    return aabb & (t0 <= t1)


def _place_labels(ax, texts, target_x, target_y, df, p, fdr_text=None):
    """Place gene labels without overlaps while exploiting empty regions of the volcano plot."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    transform = ax.transData
    inv_transform = transform.inverted()

    ax_ext = ax.patch.get_window_extent(renderer)
    margin = 4.0
    ax_x0 = ax_ext.xmin + margin
    ax_x1 = ax_ext.xmax - margin
    ax_y0 = ax_ext.ymin + margin
    ax_y1 = ax_ext.ymax - margin
    ax_w = max(1.0, ax_x1 - ax_x0)
    ax_h = max(1.0, ax_y1 - ax_y0)

    n_labels = len(texts)
    if n_labels == 0:
        return texts, []

    widths = np.empty(n_labels)
    heights = np.empty(n_labels)
    for i, t in enumerate(texts):
        bb = t.get_window_extent(renderer)
        widths[i] = bb.width
        heights[i] = bb.height

    targets = transform.transform(np.column_stack([target_x, target_y]))
    tx = targets[:, 0]
    ty = targets[:, 1]
    x_zero_px = transform.transform([[0.0, 0.0]])[0, 0]
    x_neg_lfc_px = transform.transform([[-abs(p.lfc), 0.0]])[0, 0]
    x_pos_lfc_px = transform.transform([[abs(p.lfc), 0.0]])[0, 0]

    # 2D summed-area tables (integral images) for fast O(1) empty-space queries
    cell = 2.0
    gw = int(np.ceil(ax_w / cell)) + 2
    gh = int(np.ceil(ax_h / cell)) + 2
    grid_sig = np.zeros((gh, gw), dtype=np.float64)
    grid_ns = np.zeros((gh, gw), dtype=np.float64)
    grid_hard = np.zeros((gh, gw), dtype=np.float64)

    dpi_scale = fig.dpi / 72.0
    r_sig_px = np.sqrt(p.point_size / np.pi) * dpi_scale + 2.0
    r_ns_px = np.sqrt(p.point_size * 0.35 / np.pi) * dpi_scale + 1.0

    def stamp_points(grid_arr, pts_xy, radius_px, weight, cap=3):
        if len(pts_xy) == 0:
            return
        px_coords = transform.transform(pts_xy)
        gx = ((px_coords[:, 0] - ax_x0) / cell).astype(int)
        gy = ((px_coords[:, 1] - ax_y0) / cell).astype(int)
        valid = (gx >= 0) & (gx < gw) & (gy >= 0) & (gy < gh)
        gx, gy = gx[valid], gy[valid]
        if len(gx) == 0:
            return
        flat = gy * gw + gx
        u_flat, counts = np.unique(flat, return_counts=True)
        ugy, ugx = np.divmod(u_flat, gw)
        r_cells = int(np.ceil(radius_px / cell))
        for dy in range(-r_cells, r_cells + 1):
            for dx in range(-r_cells, r_cells + 1):
                if dx * dx + dy * dy <= (radius_px / cell + 0.3) ** 2:
                    ny = np.clip(ugy + dy, 0, gh - 1)
                    nx = np.clip(ugx + dx, 0, gw - 1)
                    np.add.at(grid_arr, (ny, nx), weight * np.minimum(counts, cap))

    sig_df = df[df["cls"] != "ns"]
    ns_df = df[df["cls"] == "ns"]
    stamp_points(grid_sig, sig_df[["lfc", "y"]].to_numpy(), r_sig_px, 1.0, cap=4)
    stamp_points(grid_ns, ns_df[["lfc", "y"]].to_numpy(), r_ns_px, 1.0, cap=2)

    if fdr_text is not None:
        fbb = fdr_text.get_window_extent(renderer)
        fx0 = int(np.clip(np.floor((fbb.xmin - 4 - ax_x0) / cell), 0, gw - 1))
        fx1 = int(np.clip(np.ceil((fbb.xmax + 4 - ax_x0) / cell), 0, gw - 1))
        fy0 = int(np.clip(np.floor((fbb.ymin - 4 - ax_y0) / cell), 0, gh - 1))
        fy1 = int(np.clip(np.ceil((fbb.ymax + 4 - ax_y0) / cell), 0, gh - 1))
        grid_hard[fy0 : fy1 + 1, fx0 : fx1 + 1] += 1.0

    def make_sat(g):
        s = np.zeros((gh + 1, gw + 1), dtype=np.float64)
        s[1:, 1:] = np.cumsum(np.cumsum(g, axis=0), axis=1)
        return s

    sat_sig = make_sat(grid_sig)
    sat_ns = make_sat(grid_ns)
    sat_hard = make_sat(grid_hard)

    def query_sat(s, x0, x1, y0, y1):
        gx0 = np.clip(np.floor((x0 - ax_x0) / cell).astype(int), 0, gw)
        gx1 = np.clip(np.ceil((x1 - ax_x0) / cell).astype(int), 0, gw)
        gy0 = np.clip(np.floor((y0 - ax_y0) / cell).astype(int), 0, gh)
        gy1 = np.clip(np.ceil((y1 - ax_y0) / cell).astype(int), 0, gh)
        return s[gy1, gx1] - s[gy0, gx1] - s[gy1, gx0] + s[gy0, gx0]

    pad_x = 4.0 if n_labels <= 45 else 2.5
    pad_y = 2.5 if n_labels <= 45 else 1.5

    gx_lin = np.arange(ax_x0 + 14, ax_x1 - 14, 8.0)
    gy_lin = np.arange(ax_y0 + 7, ax_y1 - 7, 5.0)
    ggx, ggy = np.meshgrid(gx_lin, gy_lin)
    dense_grid_cands = np.column_stack([ggx.ravel(), ggy.ravel()])

    near_radii = np.array([2.0, 5.0, 9.0, 14.0, 20.0, 28.0, 38.0], dtype=float)
    near_angles = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    cos_a = np.cos(near_angles)
    sin_a = np.sin(near_angles)

    def hemisphere_bounds(i):
        half_w = widths[i] * 0.5
        xmin_b = ax_x0 + half_w
        xmax_b = ax_x1 - half_w
        if target_x[i] < 0:
            xmax_b = min(xmax_b, x_pos_lfc_px - half_w - 2.0)
        elif target_x[i] > 0:
            xmin_b = max(xmin_b, x_neg_lfc_px + half_w + 2.0)
        return xmin_b, xmax_b

    def get_candidates(i):
        wi, hi = widths[i], heights[i]
        half_w, half_h = wi * 0.5, hi * 0.5
        txi, tyi = tx[i], ty[i]

        rx = (half_w + r_sig_px) + near_radii[:, None]
        ry = (half_h + r_sig_px) + near_radii[:, None]
        cx_loc = (txi + rx * cos_a[None, :]).ravel()
        cy_loc = (tyi + ry * sin_a[None, :]).ravel()
        cands = np.vstack([np.column_stack([cx_loc, cy_loc]), dense_grid_cands])

        xmin_b, xmax_b = hemisphere_bounds(i)
        in_bounds = (
            (cands[:, 0] >= xmin_b)
            & (cands[:, 0] <= xmax_b)
            & (cands[:, 1] - half_h >= ax_y0)
            & (cands[:, 1] + half_h <= ax_y1)
        )
        return cands[in_bounds]

    local_density = query_sat(sat_sig, tx - 35, tx + 35, ty - 25, ty + 25)
    order = np.argsort(-(local_density * 2.0 + ty * 0.1))

    placed_cx = np.full(n_labels, np.nan)
    placed_cy = np.full(n_labels, np.nan)

    def choose_best_candidate(i, active_mask):
        wi, hi = widths[i], heights[i]
        half_w, half_h = wi * 0.5, hi * 0.5
        txi, tyi = tx[i], ty[i]
        cands = get_candidates(i)
        if len(cands) == 0:
            return (
                np.clip(txi, ax_x0 + half_w, ax_x1 - half_w),
                np.clip(tyi, ax_y0 + half_h, ax_y1 - half_h),
            )

        cx = cands[:, 0]
        cy = cands[:, 1]
        x0, x1 = cx - half_w, cx + half_w
        y0, y1 = cy - half_h, cy + half_h

        own_dx = np.maximum(0.0, np.maximum(x0 - txi, txi - x1))
        own_dy = np.maximum(0.0, np.maximum(y0 - tyi, tyi - y1))
        dist_edge = np.hypot(own_dx, own_dy)
        too_close_own = dist_edge < (r_sig_px + 1.5)
        hard_ov = query_sat(sat_hard, x0, x1, y0, y1) > 0

        act_idx = np.where(active_mask)[0]
        if len(act_idx) > 0:
            acx = placed_cx[act_idx]
            acy = placed_cy[act_idx]
            aw = widths[act_idx]
            ah = heights[act_idx]
            ov_x = np.maximum(0.0, (wi + aw[None, :]) * 0.5 + pad_x - np.abs(cx[:, None] - acx[None, :]))
            ov_y = np.maximum(0.0, (hi + ah[None, :]) * 0.5 + pad_y - np.abs(cy[:, None] - acy[None, :]))
            label_ov_area = np.sum(ov_x * ov_y, axis=1)
            has_label_ov = label_ov_area > 0.0
        else:
            label_ov_area = np.zeros(len(cands))
            has_label_ov = np.zeros(len(cands), dtype=bool)

        valid_mask = (~has_label_ov) & (~too_close_own) & (~hard_ov)
        if np.any(valid_mask):
            pool = np.where(valid_mask)[0]
        elif np.any(~has_label_ov):
            pool = np.where(~has_label_ov)[0]
        else:
            pool = np.arange(len(cands))

        pcx, pcy = cx[pool], cy[pool]
        px0, px1 = x0[pool], x1[pool]
        py0, py1 = y0[pool], y1[pool]
        p_dist = dist_edge[pool]
        p_lov = label_ov_area[pool]

        sig_ov = query_sat(sat_sig, px0 - 1.0, px1 + 1.0, py0 - 1.0, py1 + 1.0)
        ns_ov = query_sat(sat_ns, px0, px1, py0, py1)

        cost = (
            p_dist * 2.0
            + 0.008 * (p_dist ** 2)
            + np.where(sig_ov > 0, 5000.0 + sig_ov * 180.0, 0.0)
            + np.where(ns_ov > 0, 120.0 + ns_ov * 12.0, 0.0)
            + p_lov * 10000.0
        )

        if target_x[i] < 0:
            cost += np.maximum(0.0, px1 - x_zero_px) * 3.0
        elif target_x[i] > 0:
            cost += np.maximum(0.0, x_zero_px - px0) * 3.0

        if len(act_idx) > 0:
            if len(pool) > 350:
                top_k = np.argpartition(cost, 350)[:350]
                pool = pool[top_k]
                pcx, pcy = cx[pool], cy[pool]
                px0, px1, py0, py1 = x0[pool], x1[pool], y0[pool], y1[pool]
                cost = cost[top_k]

            p1 = np.broadcast_to([txi, tyi], (len(pool), 2))
            p2 = np.column_stack([pcx, pcy])
            q1 = np.column_stack([tx[act_idx], ty[act_idx]])
            q2 = np.column_stack([acx, acy])

            n_cross = np.sum(_segments_intersect(p1, p2, q1, q2), axis=1)
            abx0, abx1 = acx - aw * 0.5 - 1.0, acx + aw * 0.5 + 1.0
            aby0, aby1 = acy - ah * 0.5 - 1.0, acy + ah * 0.5 + 1.0
            n_seg_box = np.sum(_seg_intersects_box(p1, p2, abx0, abx1, aby0, aby1), axis=1)
            n_box_seg = np.sum(
                _seg_intersects_box(q1, q2, px0 - 1.0, px1 + 1.0, py0 - 1.0, py1 + 1.0), axis=0
            )
            cost = cost + n_cross * 450.0 + n_seg_box * 160.0 + n_box_seg * 240.0

        best = np.argmin(cost)
        return pcx[best], pcy[best]

    # Pass 1: Greedy placement
    active = np.zeros(n_labels, dtype=bool)
    for i in order:
        placed_cx[i], placed_cy[i] = choose_best_candidate(i, active)
        active[i] = True

    # Pass 2: Re-placement sweep (each label sees all other N-1 labels)
    for i in order[::-1]:
        active[i] = False
        placed_cx[i], placed_cy[i] = choose_best_candidate(i, active)
        active[i] = True

    # Pass 3: 2-opt uncrossing & distance optimization
    nudge_offsets = np.array(
        [
            [0, 0], [-4, 0], [4, 0], [0, -3], [0, 3],
            [-8, 0], [8, 0], [0, -6], [0, 6],
            [-14, 0], [14, 0], [-6, -4], [6, -4], [-6, 4], [6, 4],
        ],
        dtype=float,
    )

    def try_fit_box(cx_base, cy_base, i, ocx, ocy, ow, oh, max_sig_ov=0.0):
        w, h = widths[i], heights[i]
        half_w, half_h = w * 0.5, h * 0.5
        xmin_b, xmax_b = hemisphere_bounds(i)
        cxs = np.clip(cx_base + nudge_offsets[:, 0], xmin_b, xmax_b)
        cys = np.clip(cy_base + nudge_offsets[:, 1], ax_y0 + half_h, ax_y1 - half_h)

        x0, x1 = cxs - half_w, cxs + half_w
        y0, y1 = cys - half_h, cys + half_h
        own_dx = np.maximum(0.0, np.maximum(x0 - tx[i], tx[i] - x1))
        own_dy = np.maximum(0.0, np.maximum(y0 - ty[i], ty[i] - y1))
        ok_own = np.hypot(own_dx, own_dy) >= (r_sig_px + 1.5)
        ok_hard = query_sat(sat_hard, x0, x1, y0, y1) == 0
        ok_sig = query_sat(sat_sig, x0 - 1.0, x1 + 1.0, y0 - 1.0, y1 + 1.0) <= max_sig_ov

        if len(ocx) > 0:
            ov = (np.abs(cxs[:, None] - ocx[None, :]) < (w + ow[None, :]) * 0.5 + pad_x) & (
                np.abs(cys[:, None] - ocy[None, :]) < (h + oh[None, :]) * 0.5 + pad_y
            )
            ok_ov = ~np.any(ov, axis=1)
        else:
            ok_ov = np.ones(len(cxs), dtype=bool)

        valid = np.where(ok_own & ok_hard & ok_sig & ok_ov)[0]
        if len(valid) == 0:
            return None
        return cxs[valid[0]], cys[valid[0]]

    for _ in range(8):
        changed = False
        for a in range(n_labels):
            for b in range(a + 1, n_labels):
                if (target_x[a] * target_x[b]) < 0:
                    continue
                p1_old = np.array([[tx[a], ty[a]]])
                p2_old = np.array([[placed_cx[a], placed_cy[a]]])
                q1_old = np.array([[tx[b], ty[b]]])
                q2_old = np.array([[placed_cx[b], placed_cy[b]]])
                crossed_before = _segments_intersect(p1_old, p2_old, q1_old, q2_old)[0, 0]

                dist_before = np.hypot(placed_cx[a] - tx[a], placed_cy[a] - ty[a]) + np.hypot(
                    placed_cx[b] - tx[b], placed_cy[b] - ty[b]
                )
                dist_after = np.hypot(placed_cx[b] - tx[a], placed_cy[b] - ty[a]) + np.hypot(
                    placed_cx[a] - tx[b], placed_cy[a] - ty[b]
                )
                if not crossed_before and dist_after >= dist_before - 3.0:
                    continue

                others = [k for k in range(n_labels) if k != a and k != b]
                ocx = placed_cx[others] if others else np.empty(0)
                ocy = placed_cy[others] if others else np.empty(0)
                ow = widths[others] if others else np.empty(0)
                oh = heights[others] if others else np.empty(0)

                cur_sig_a = query_sat(
                    sat_sig,
                    placed_cx[b] - widths[b] * 0.5 - 1.0,
                    placed_cx[b] + widths[b] * 0.5 + 1.0,
                    placed_cy[b] - heights[b] * 0.5 - 1.0,
                    placed_cy[b] + heights[b] * 0.5 + 1.0,
                )
                cur_sig_b = query_sat(
                    sat_sig,
                    placed_cx[a] - widths[a] * 0.5 - 1.0,
                    placed_cx[a] + widths[a] * 0.5 + 1.0,
                    placed_cy[a] - heights[a] * 0.5 - 1.0,
                    placed_cy[a] + heights[a] * 0.5 + 1.0,
                )

                fit_a = try_fit_box(placed_cx[b], placed_cy[b], a, ocx, ocy, ow, oh, cur_sig_a)
                if fit_a is None:
                    continue
                cxa, cya = fit_a
                fit_b = try_fit_box(
                    placed_cx[a],
                    placed_cy[a],
                    b,
                    np.append(ocx, cxa),
                    np.append(ocy, cya),
                    np.append(ow, widths[a]),
                    np.append(oh, heights[a]),
                    cur_sig_b,
                )
                if fit_b is None:
                    continue
                cxb, cyb = fit_b

                if _segments_intersect(
                    np.array([[tx[a], ty[a]]]),
                    np.array([[cxa, cya]]),
                    np.array([[tx[b], ty[b]]]),
                    np.array([[cxb, cyb]]),
                )[0, 0]:
                    continue

                placed_cx[a], placed_cy[a] = cxa, cya
                placed_cx[b], placed_cy[b] = cxb, cyb
                changed = True
        if not changed:
            break

    # Pass 4: Inward compaction (pull labels closer to their target points along empty space)
    dists_now = np.hypot(placed_cx - tx, placed_cy - ty)
    compact_order = np.argsort(dists_now)
    t_steps = np.linspace(0.15, 0.92, 12)
    lat_shifts = np.array([-8.0, 0.0, 8.0])
    for i in compact_order:
        if dists_now[i] < 18.0:
            continue
        wi, hi = widths[i], heights[i]
        half_w, half_h = wi * 0.5, hi * 0.5
        xmin_b, xmax_b = hemisphere_bounds(i)
        cur_sig = query_sat(
            sat_sig,
            placed_cx[i] - half_w - 1.0,
            placed_cx[i] + half_w + 1.0,
            placed_cy[i] - half_h - 1.0,
            placed_cy[i] + half_h + 1.0,
        )
        cur_ns = query_sat(
            sat_ns,
            placed_cx[i] - half_w,
            placed_cx[i] + half_w,
            placed_cy[i] - half_h,
            placed_cy[i] + half_h,
        )

        # Ray from target to current center
        vx = placed_cx[i] - tx[i]
        vy = placed_cy[i] - ty[i]
        v_len = max(1.0, np.hypot(vx, vy))
        px_perp = -vy / v_len
        py_perp = vx / v_len

        base_x = tx[i] + vx * t_steps[:, None] + px_perp * lat_shifts[None, :]
        base_y = ty[i] + vy * t_steps[:, None] + py_perp * lat_shifts[None, :]
        cxs = np.clip(base_x.ravel(), xmin_b, xmax_b)
        cys = np.clip(base_y.ravel(), ax_y0 + half_h, ax_y1 - half_h)

        x0, x1 = cxs - half_w, cxs + half_w
        y0, y1 = cys - half_h, cys + half_h
        own_dx = np.maximum(0.0, np.maximum(x0 - tx[i], tx[i] - x1))
        own_dy = np.maximum(0.0, np.maximum(y0 - ty[i], ty[i] - y1))
        ok = (
            (np.hypot(own_dx, own_dy) >= (r_sig_px + 1.5))
            & (query_sat(sat_hard, x0, x1, y0, y1) == 0)
            & (query_sat(sat_sig, x0 - 1.0, x1 + 1.0, y0 - 1.0, y1 + 1.0) <= cur_sig)
            & (query_sat(sat_ns, x0, x1, y0, y1) <= cur_ns)
        )

        others = [k for k in range(n_labels) if k != i]
        if others and np.any(ok):
            ocx, ocy = placed_cx[others], placed_cy[others]
            ow, oh = widths[others], heights[others]
            ov = (np.abs(cxs[:, None] - ocx[None, :]) < (wi + ow[None, :]) * 0.5 + pad_x) & (
                np.abs(cys[:, None] - ocy[None, :]) < (hi + oh[None, :]) * 0.5 + pad_y
            )
            ok &= ~np.any(ov, axis=1)

            if np.any(ok):
                idx_ok = np.where(ok)[0]
                p1 = np.broadcast_to([tx[i], ty[i]], (len(idx_ok), 2))
                p2 = np.column_stack([cxs[idx_ok], cys[idx_ok]])
                q1 = np.column_stack([tx[others], ty[others]])
                q2 = np.column_stack([ocx, ocy])
                no_cross = np.sum(_segments_intersect(p1, p2, q1, q2), axis=1) == 0
                if np.any(no_cross):
                    pick = idx_ok[np.where(no_cross)[0][0]]
                    placed_cx[i], placed_cy[i] = cxs[pick], cys[pick]

    patches = []
    final_data = inv_transform.transform(np.column_stack([placed_cx, placed_cy]))
    for i, text in enumerate(texts):
        text_mid = final_data[i]
        target = (target_x[i], target_y[i])
        text.set_verticalalignment("center")
        text.set_horizontalalignment("center")
        text.set_position(text_mid)

        wi, hi = widths[i], heights[i]
        dx = max(0.0, abs(placed_cx[i] - tx[i]) - wi * 0.5)
        dy = max(0.0, abs(placed_cy[i] - ty[i]) - hi * 0.5)
        if np.hypot(dx, dy) >= 3.5:
            arrowpatch = FancyArrowPatch(
                posA=text_mid,
                posB=target,
                patchA=text,
                transform=transform,
                arrowstyle="-",
                color="#8a8f97",
                lw=1.0,
                shrinkA=2,
                shrinkB=4,
                zorder=4,
            )
            ax.add_patch(arrowpatch)
            patches.append(arrowpatch)

    return texts, patches


def _labels(df, p):
    if p.label_mode == "none":
        return df.iloc[:0]
    if p.label_mode == "custom":
        want = {g.strip().lower() for g in p.label_genes if g.strip()}
        return df[df["symbol"].str.lower().isin(want)]
    sig = df[df["cls"] != "ns"]
    if p.label_mode == "top":
        return rank_genes(sig, p.label_rank).head(p.label_top)
    return sig.sort_values("padj").head(60)   # "significant": all, capped for legibility


def rank_genes(d, by="fdr"):
    """Order genes for labelling: smallest FDR, largest |log2FC|, or farthest from the
    origin as drawn (both axes scaled to their range, so neither dominates)."""
    if by == "lfc":
        key = -d["lfc"].abs()
    elif by == "distance":
        x = d["lfc"].abs() / (d["lfc"].abs().max() or 1)
        y = d["y"] / (d["y"].max() or 1)
        key = -(x ** 2 + y ** 2)
    else:
        key = d["padj"]
    return d.assign(_k=key).sort_values(["_k", "padj"], kind="mergesort").drop(columns="_k")
