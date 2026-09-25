"""Shared plotting style for the locality exploration (reference palette, light mode)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # categorical slots 1-4
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
SURF = "#fcfcfb"
# diverging blue <-> gray <-> red (reference diverging pair)
DIV = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#9b1c1c"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10, "legend.fontsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e2e1dc", "grid.linewidth": 0.5, "axes.edgecolor": "#52514e",
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "savefig.dpi": 170, "savefig.bbox": "tight", "lines.linewidth": 1.8,
    "legend.frameon": False,
})
