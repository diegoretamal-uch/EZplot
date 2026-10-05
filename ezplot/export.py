"""Saving figures for journals: PNG/TIFF at a chosen DPI and width, PDF/SVG with
editable text, and a zip bundle with the parameters used."""

import io as _io
import json
import zipfile
from datetime import datetime

import matplotlib.pyplot as plt
from PIL import Image

FORMATS = {"png": "PNG", "tiff": "TIFF (LZW)", "pdf": "PDF (vector)", "svg": "SVG (vector)"}
WIDTHS_MM = {"original": None, "single column (89 mm)": 89, "1.5 column (120 mm)": 120,
             "double column (183 mm)": 183}


def figure_bytes(fig, fmt="png", dpi=300, width_mm=None):
    """Render `fig` to bytes. For raster formats `width_mm` sets the printed width
    (text shrinks with the figure); vector formats keep the design size."""
    buf = _io.BytesIO()
    kw = dict(bbox_inches="tight", pad_inches=0.08, facecolor="white")
    if fmt in ("pdf", "svg"):
        fig.savefig(buf, format=fmt, **kw)
        return buf.getvalue()
    fig.savefig(buf, format="png", dpi=dpi, **kw)
    img = Image.open(_io.BytesIO(buf.getvalue()))
    if width_mm:
        target_px = round(width_mm / 25.4 * dpi)
        if target_px != img.width:
            h = round(img.height * target_px / img.width)
            img = img.resize((target_px, h), Image.LANCZOS)
    out = _io.BytesIO()
    if fmt == "tiff":
        img.convert("RGB").save(out, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    else:
        img.save(out, format="PNG", dpi=(dpi, dpi), optimize=True)
    return out.getvalue()


def bundle(figures, tables, params, fmts=("png", "pdf"), dpi=300, width_mm=None):
    """figures: name -> callable returning a Figure; tables: name -> DataFrame."""
    buf = _io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, make in figures.items():
            fig = make()
            try:
                for f in fmts:
                    z.writestr(f"figures/{name}.{f}", figure_bytes(fig, f, dpi, width_mm))
            finally:
                plt.close(fig)
        for name, df in tables.items():
            z.writestr(f"tables/{name}.tsv", df.to_csv(sep="\t", index=False))
        meta = {"created": datetime.now().isoformat(timespec="seconds"),
                "dpi": dpi, "width_mm": width_mm, "formats": list(fmts), **params}
        z.writestr("ezplot_parameters.json", json.dumps(meta, indent=2, default=str))
    return buf.getvalue()
