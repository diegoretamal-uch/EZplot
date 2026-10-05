"""Command line: `ezplot` opens the app; `ezplot batch ...` renders everything headless."""

import argparse
import socket
import sys
from pathlib import Path


def free_port(preferred=8765):
    for port in [preferred] + list(range(preferred + 1, preferred + 50)):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
                return port
            except OSError:
                continue
    return 0


def launch(port=None, browser=True):
    import shiny

    port = port or free_port()
    print(f"\n  EZplot is running at  http://127.0.0.1:{port}\n"
          "  Keep this window open while you use it; close it to quit.\n", flush=True)
    shiny.run_app("ezplot.app:app", host="127.0.0.1", port=port, launch_browser=browser,
                  log_level="warning")


def batch(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from . import export
    from .project import Project

    pr = Project()
    for f in [a.de, a.heatmap, a.gmt]:
        if f:
            kind = pr.load(f, Path(f).name)
            print(f"  {Path(f).name}: {kind}")
    if a.control:
        pr.control = a.control
        pr.treatment = a.treatment or next(g for g in pr.group_names if g != a.control)
    if a.species:
        pr.species = a.species
    if a.flip:
        pr.up_group = next(g for g in pr.group_names if g != pr.up_group)
    pr.control_label, pr.treatment_label = a.control_label, a.treatment_label
    if pr.de is not None:
        print(f"  {pr.direction_sentence()}  (auto-check agreement {pr.direction_confidence:.0%})")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    figs = {}
    if pr.de is not None:
        figs["volcano"] = pr.volcano
    if pr.heatmap is not None:
        figs["heatmap"] = pr.heatmap_fig
    keys = [k for k in (a.collections.split(",") if a.collections else pr.collections()) if k]
    if pr.de is not None:
        for k in keys:
            if k not in pr.collections():
                sys.exit(f"Unknown collection '{k}'. Available: {', '.join(pr.collections())}")
            res = pr.run_gsea(k, progress=None)
            res.to_csv(out / f"gsea_{k}.tsv", sep="\t", index=False)
            print(f"  GSEA {k}: {len(res)} sets, {(res.padj < 0.05).sum()} at FDR < 0.05")
            figs[f"gsea_{k}"] = (lambda k=k: pr.gsea_fig(k))
    width = a.width_mm or None
    for name, make in figs.items():
        fig = make()
        for fmt in a.formats.split(","):
            (out / f"{name}.{fmt}").write_bytes(export.figure_bytes(fig, fmt, a.dpi, width))
        plt.close(fig)
        print(f"  wrote {name}.{{{a.formats}}}")


def main(argv=None):
    # Windows consoles default to cp1252, which can't print "log₂FC"
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(prog="ezplot", description=__doc__)
    sub = ap.add_subparsers(dest="cmd")
    ap.add_argument("--port", type=int, default=None)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--check", action="store_true", help="verify the install and exit")
    b = sub.add_parser("batch", help="render all plots without the app")
    b.add_argument("--de", help="iDEP DE results matrix (csv)")
    b.add_argument("--heatmap", help="iDEP DEG_Heatmap_Data.csv")
    b.add_argument("--gmt", help="custom gene sets (.gmt)")
    b.add_argument("--control", help="control group name (default: auto)")
    b.add_argument("--treatment", help="treatment group name (default: the other group)")
    b.add_argument("--control-label")
    b.add_argument("--treatment-label")
    b.add_argument("--flip", action="store_true", help="override the auto-detected sign")
    b.add_argument("--species", choices=["mouse", "human"])
    b.add_argument("--collections", help="comma list: hallmarks,go_bp,go_mf,go_cc,custom")
    b.add_argument("--formats", default="png,pdf")
    b.add_argument("--dpi", type=int, default=300)
    b.add_argument("--width-mm", type=float, default=0)
    b.add_argument("--out", default="ezplot_output")
    a = ap.parse_args(argv)
    if a.check:
        import ezplot.app  # noqa: F401
        from . import __version__
        print(f"EZplot {__version__}: installation OK")
    elif a.cmd == "batch":
        batch(a)
    else:
        launch(a.port, not a.no_browser)


if __name__ == "__main__":
    main()
