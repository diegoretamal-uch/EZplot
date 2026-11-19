"""Saving figures for journals: PNG/TIFF at a chosen DPI and width, PDF/SVG with
editable text, and a zip bundle with the parameters used."""

import io as _io
import json
import zipfile
from contextlib import contextmanager
from datetime import datetime

import matplotlib.pyplot as plt
from PIL import Image

FORMATS = {"png": "PNG", "tiff": "TIFF (LZW)", "pdf": "PDF (vector)", "svg": "SVG (vector)"}
WIDTHS_MM = {"original": None, "single column (89 mm)": 89, "1.5 column (120 mm)": 120,
             "double column (183 mm)": 183}

# matplotlib's Agg renderer refuses images of 2**16 pixels or more on either side; stay
# below that with room for bbox_inches="tight" growing past the figure size
RASTER_MAX_PX = 60000


def _raster_dpi(fig, dpi):
    """`dpi`, lowered just enough that the longest side fits in RASTER_MAX_PX."""
    longest_in = max(max(fig.get_size_inches()), 1e-6)
    return max(1, min(int(dpi), int(RASTER_MAX_PX // longest_in)))


def _render_png(fig, dpi, kw):
    """PNG bytes of `fig` and the DPI actually used (lower than asked only for huge figures)."""
    eff = _raster_dpi(fig, dpi)
    for _ in range(3):
        buf = _io.BytesIO()
        try:
            fig.savefig(buf, format="png", dpi=eff, **kw)
            return buf.getvalue(), eff
        except ValueError as e:
            if "too large" not in str(e) or eff <= 1:
                raise
            eff = max(1, eff // 2)        # tight bbox grew past the estimate: retry smaller
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=eff, **kw)
    return buf.getvalue(), eff


@contextmanager
def _trusted_image_size():
    """Lift PIL's decompression-bomb guard while handling images we rendered ourselves.

    A tall heatmap at 300 dpi (~750 genes) or almost any heatmap at 1200 dpi exceeds
    PIL's 179-megapixel limit, and Image.open() would raise DecompressionBombError.
    """
    old = Image.MAX_IMAGE_PIXELS
    Image.MAX_IMAGE_PIXELS = None
    try:
        yield
    finally:
        Image.MAX_IMAGE_PIXELS = old


def figure_bytes(fig, fmt="png", dpi=300, width_mm=None):
    """Render `fig` to bytes. For raster formats `width_mm` sets the printed width
    (text shrinks with the figure); vector formats keep the design size."""
    kw = dict(bbox_inches="tight", pad_inches=0.08, facecolor="white")
    if fmt in ("pdf", "svg"):
        buf = _io.BytesIO()
        fig.savefig(buf, format=fmt, **kw)
        return buf.getvalue()
    png, eff = _render_png(fig, dpi, kw)
    if fmt != "tiff" and not width_mm:
        return png                        # matplotlib already wrote the DPI into the PNG
    with _trusted_image_size():
        img = Image.open(_io.BytesIO(png))
        img.load()
        out_dpi = eff
        if width_mm:
            out_dpi = dpi                 # printed width is fixed, so the density is as asked
            target_px = round(width_mm / 25.4 * dpi)
            if target_px != img.width:
                h = max(1, round(img.height * target_px / img.width))
                img = img.resize((target_px, h), Image.LANCZOS)
        out = _io.BytesIO()
        if fmt == "tiff":
            img.convert("RGB").save(out, format="TIFF", compression="tiff_lzw",
                                    dpi=(out_dpi, out_dpi))
        else:
            img.save(out, format="PNG", dpi=(out_dpi, out_dpi), optimize=True)
    return out.getvalue()


def _describe(exc):
    return f"{type(exc).__name__}: {exc}"


def bundle(figures, tables, params, fmts=("png", "pdf"), dpi=300, width_mm=None, errors=None):
    """figures: name -> callable returning a Figure; tables: name -> DataFrame.

    A figure or table that fails is left out and explained in errors.txt inside the zip,
    together with any `errors` the caller collected (e.g. GSEA runs that failed), so one
    problem never turns the whole download into an empty or broken file.
    """
    errs = [str(e) for e in (errors or [])]
    buf = _io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, make in figures.items():
            open_before = set(plt.get_fignums())
            try:
                fig = make()
            except Exception as exc:  # noqa: BLE001 - report it in the zip instead
                for num in set(plt.get_fignums()) - open_before:
                    plt.close(num)        # half-built figure left behind by the failure
                errs.append(f"figures/{name}: could not draw the figure ({_describe(exc)})")
                continue
            try:
                for f in fmts:
                    try:
                        data = figure_bytes(fig, f, dpi, width_mm)
                    except Exception as exc:  # noqa: BLE001
                        errs.append(f"figures/{name}.{f}: could not save ({_describe(exc)})")
                        continue
                    z.writestr(f"figures/{name}.{f}", data)
            finally:
                plt.close(fig)
        for name, df in tables.items():
            try:
                text = df.to_csv(sep="\t", index=False)
            except Exception as exc:  # noqa: BLE001
                errs.append(f"tables/{name}.tsv: could not write the table ({_describe(exc)})")
                continue
            z.writestr(f"tables/{name}.tsv", text)
        meta = {"created": datetime.now().isoformat(timespec="seconds"),
                "dpi": dpi, "width_mm": width_mm, "formats": list(fmts), **params}
        if errs:
            meta["errors"] = errs
        z.writestr("ezplot_parameters.json", json.dumps(meta, indent=2, default=str))
        if errs:
            z.writestr("errors.txt",
                       "Some items could not be produced and are missing from this zip:\n\n"
                       + "\n".join(f"- {e}" for e in errs) + "\n")
    return buf.getvalue()
