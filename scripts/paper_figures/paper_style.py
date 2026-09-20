"""
Shared figure style for the 2026 High-Redshift-Black-Hole-Growth paper.

Uses SciencePlots (styles 'science' and 'bright') when it is installed, and falls back to the local copy of the
'science' style (mnras_science.mplstyle) otherwise. Figures are drawn at the size at which they are printed
(single column 3.4 in, full width 7.1 in for the openjournal two-column class), so that the font sizes below are the
printed sizes. Ordered quantities (seed mass, halo mass, redshift band, eta_acc) use a sequential, colour-blind-safe,
greyscale-safe palette; categorical quantities use the SciencePlots 'bright' cycle. Every figure is written as a
vector PDF (scatter points rasterised) and as a PNG.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
COL, FULL = 3.4, 7.1                     # inches: column width and text width of the paper's class

BLUE, RED, GREEN, YELLOW, CYAN, PURPLE, GREY = "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"


def apply():
    try:
        import scienceplots  # noqa: F401
        plt.style.use(["science", "bright"])
    except Exception:
        plt.style.use(str(HERE / "mnras_science.mplstyle"))
    plt.rcParams.update({
        "axes.grid": False, "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "legend.title_fontsize": 7,
        "lines.linewidth": 1.2, "savefig.dpi": 300, "pdf.fonttype": 42,
    })


def seq(n, lo=0.05, hi=0.85):
    """n colours from the viridis map, dark to light (ordered quantities)."""
    return [plt.cm.viridis(x) for x in np.linspace(lo, hi, n)]


def save(fig, stem):
    """Write stem.pdf and stem.png (300 dpi)."""
    stem = Path(stem)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
