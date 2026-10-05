"""Preranked GSEA in plain numpy (runs natively and in the browser build).

Follows fgsea (gseaParam = 1):
  ES    weighted Kolmogorov-Smirnov running sum over genes ranked by the statistic.
  null  fgseaSimple's: `nperm` random gene sets of the same size. Random sets of every
        size are prefixes of the same random permutations, so one draw serves all sizes.
  NES   ES / mean(null ES of the same sign).
  pval  (#null at least as extreme + 1) / (#null of the same sign + 1); sets beyond the
        permutation floor are refined with fgsea's adaptive multilevel splitting
        (Korotkevich et al., bioRxiv 060012), grouped by set size and sign.
  padj  Benjamini-Hochberg over all tested sets.
Validated against fgsea 1.x on real RNA-seq data: median |delta NES| ~ 0.01 and log10 p
correlation 0.97 for sets refined by multilevel splitting.
"""

import numpy as np
import pandas as pd


def _running(pos, w, n):
    """Running-sum values just after (top) and just before (bot) each hit."""
    k = pos.shape[1]
    hw = w[pos]
    nr = hw.sum(axis=1, keepdims=True)
    nr[nr == 0] = 1.0
    top = np.cumsum(hw, axis=1) / nr - (pos - np.arange(k)) / (n - k)
    return top, top - hw / nr


def _es(pos, w, n):
    """ES for gene sets given as sorted rank positions, shape (m, k)."""
    top, bot = _running(pos, w, n)
    mx, mn = top.max(axis=1), bot.min(axis=1)
    return np.where(mx > -mn, mx, mn)


def _es_pos(pos, w, n):
    """Positive deviation only (the multilevel score)."""
    return _running(pos, w, n)[0].max(axis=1)


def _random_prefixes(rng, n, kmax, nperm, chunk=256):
    """nperm x kmax matrix whose first k columns are a uniform random k-subset."""
    out = np.empty((nperm, kmax), dtype=np.int32)
    for i in range(0, nperm, chunk):
        m = min(chunk, nperm - i)
        keys = rng.random((m, n), dtype=np.float32)
        sel = np.argpartition(keys, kmax - 1, axis=1)[:, :kmax]
        order = np.argsort(np.take_along_axis(keys, sel, axis=1), axis=1)
        out[i:i + m] = np.take_along_axis(sel, order, axis=1)
    return out


def _multilevel(w, n, k, emax, rng, sample_size=101, step_frac=0.1, floor=1e-50):
    """Adaptive multilevel splitting for the positive tail of ES among size-k sets.

    Returns a list of levels (log P(score > previous level), sorted sample scores).
    Samples of a level follow the score distribution conditioned on exceeding the
    previous level, so P(score >= e) = exp(log_p) * fraction of its samples >= e.
    """
    s = sample_size
    keys = rng.random((s, n), dtype=np.float32)
    pos = np.sort(np.argpartition(keys, k - 1, axis=1)[:, :k], axis=1)
    score = _es_pos(pos, w, n)
    rows = np.arange(s)
    nsteps = max(10, int(k * step_frac))
    levels, lp = [], 0.0
    while True:
        levels.append((lp, np.sort(score)))
        thr = np.median(score)
        keep = np.flatnonzero(score > thr)
        if thr >= emax or len(keep) == 0 or lp < np.log(floor):
            break
        lp += np.log(0.5)              # P(X > sample median) = 1/2 in expectation
        pick = rng.choice(keep, s)
        pos, score = pos[pick].copy(), score[pick].copy()
        for _ in range(nsteps):        # Metropolis swaps that stay above the level
            j = rng.integers(0, k, s)
            g = rng.integers(0, n, s)
            fresh = ~(pos == g[:, None]).any(axis=1)
            new = pos.copy()
            new[rows, j] = g
            new.sort(axis=1)
            ns = _es_pos(new, w, n)
            ok = fresh & (ns > thr)
            pos[ok], score[ok] = new[ok], ns[ok]
    return levels


def _tail_p(levels, e):
    """P(score >= e) from the first level whose samples have a median >= e."""
    for lp, sc in levels:
        if np.median(sc) >= e:
            break
    frac = (len(sc) - np.searchsorted(sc, e, side="left")) / len(sc)
    return np.exp(lp) * max(frac, 1.0 / (len(sc) + 1))


def _bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def prepare_ranking(stats):
    """Series indexed by gene symbol -> sorted (descending) Series without NA/duplicates."""
    s = pd.Series(stats, dtype=float).dropna()
    s = s[s.index.notna() & (s.index.astype(str).str.strip() != "")]
    s = s[~s.index.duplicated(keep="first")]
    return s.sort_values(ascending=False, kind="mergesort")


def prerank(stats, gene_sets, min_size=15, max_size=500, nperm=1000, seed=42,
            multilevel=True, progress=None):
    """Run preranked GSEA. Returns a DataFrame sorted by padj (ties by |NES|).

    stats      Series: gene symbol -> ranking statistic (oriented log2FC)
    gene_sets  dict: name -> list of genes
    progress   optional callable(fraction, message)
    """
    def tick(f, msg):
        if progress:
            progress(f, msg)

    r = prepare_ranking(stats)
    genes = r.index.to_numpy()
    n = len(genes)
    w = np.abs(r.to_numpy())
    where = {g: i for i, g in enumerate(genes)}

    names, positions = [], []
    for name, members in gene_sets.items():
        p = np.unique([where[g] for g in members if g in where])
        if min_size <= len(p) <= max_size:
            names.append(name)
            positions.append(p)
    cols = ["pathway", "ES", "NES", "pval", "padj", "size", "leadingEdge"]
    if not names:
        return pd.DataFrame(columns=cols)

    sizes = np.array([len(p) for p in positions])
    es = np.array([_es(p[None, :], w, n)[0] for p in positions])
    lead = [";".join(_leading_edge(p, w, n, e, genes)) for p, e in zip(positions, es)]

    rng = np.random.default_rng(seed)
    uniq = np.unique(sizes)
    n_same, n_more, mean_null = (np.zeros(len(names)) for _ in range(3))
    prefixes = _random_prefixes(rng, n, int(uniq.max()), nperm)
    for j, k in enumerate(uniq):
        tick(0.6 * (j + 1) / len(uniq), f"Permutations · set size {k}")
        null = _es(np.sort(prefixes[:, :k], axis=1), w, n)
        pos_null, neg_null = null[null >= 0], null[null <= 0]
        for i in np.flatnonzero(sizes == k):
            same = pos_null if es[i] >= 0 else neg_null
            n_same[i] = len(same)
            n_more[i] = np.sum(np.abs(same) >= abs(es[i]))
            mean_null[i] = np.abs(same).mean() if len(same) else np.nan
    del prefixes
    pval = (n_more + 1) / (n_same + 1)

    # sets at the permutation floor: multilevel splitting, one run per (size, sign)
    floor_sets = np.flatnonzero(n_more < 10) if multilevel else np.array([], int)
    groups = sorted({(sizes[i], es[i] >= 0) for i in floor_sets})
    w_rev = w[::-1].copy()
    for j, (k, positive) in enumerate(groups):
        tick(0.6 + 0.4 * (j + 1) / len(groups), f"Refining small p-values · set size {k}")
        idx = [i for i in floor_sets if sizes[i] == k and (es[i] >= 0) == positive]
        e = np.abs(es[idx])
        levels = _multilevel(w if positive else w_rev, n, int(k), e.max(), rng)
        for i, ei in zip(idx, e):
            p_tail = _tail_p(levels, ei)
            p_same = n_same[i] / nperm
            if p_same > 0:
                pval[i] = min(pval[i], p_tail / p_same)

    nes = es / mean_null
    res = pd.DataFrame({"pathway": names, "ES": es, "NES": nes, "pval": pval,
                        "padj": _bh(pval), "size": sizes, "leadingEdge": lead})
    tick(1.0, "Done")
    res["_abs"] = -res["NES"].abs()
    res = res.sort_values(["padj", "pval", "_abs"], kind="mergesort").drop(columns="_abs")
    return res.reset_index(drop=True)


def _leading_edge(p, w, n, es, genes):
    top, bot = _running(p[None, :], w, n)
    if es >= 0:
        hits = p[: np.argmax(top[0]) + 1]
    else:
        hits = p[np.argmin(bot[0]):][::-1]
    return list(genes[hits])
