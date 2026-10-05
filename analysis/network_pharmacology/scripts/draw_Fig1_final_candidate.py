#!/usr/bin/env python3
"""Redraw frozen Figure 1 inputs without recalculating network or enrichment results."""

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import Ellipse
import networkx as nx
import numpy as np
import pandas as pd


MODULE = Path(__file__).resolve().parents[1]
SOURCE = MODULE / "data" / "derived"
OUT = MODULE / "reproduced_figures"
DATA = SOURCE
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

COLORS = {
    "teal": "#0B5D5E",
    "cyan": "#3297AD",
    "ochre": "#C5964A",
    "coral": "#E76551",
    "burgundy": "#A5131C",
    "ink": "#202124",
    "muted": "#69727A",
    "gray": "#BCC2C5",
    "pale": "#EEF2F3",
}


def configure_fonts():
    arial = Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
    if arial.exists():
        font_manager.fontManager.addfont(str(arial))
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.55,
        "savefig.facecolor": "white",
    })


def panel_label(ax, label):
    ax.text(-0.075, 1.045, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left", color=COLORS["ink"])


def clean_axis(ax, grid=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.55)
    ax.spines["bottom"].set_linewidth(0.55)
    ax.tick_params(width=0.55, length=3, color=COLORS["ink"])
    if grid:
        ax.grid(axis="x", color="#E6E9EB", linewidth=0.45, zorder=0)


def load_frozen_inputs():
    overlap = pd.read_csv(DATA / "Fig1A_54_overlap_genes.csv")["gene"].astype(str)
    degree = pd.read_csv(DATA / "Fig1B_all_nodes_frozen_metrics.csv")
    edges = pd.read_csv(DATA / "Fig1B_all_543_STRING_edges.csv")
    top15 = pd.read_csv(DATA / "Fig1C_top15_degree.csv")
    selected = pd.read_csv(DATA / "Fig1D_selected_enrichment.csv")

    assert overlap.nunique() == 54 and len(overlap) == 54
    assert len(degree) == 54 and degree["name"].nunique() == 54
    # Preserve the archived files exactly: overlap has VEGFA whereas the frozen
    # Cytoscape/STRING node set has COL18A1. Do not silently substitute either.
    assert set(overlap) - set(degree["name"]) == {"VEGFA"}
    assert set(degree["name"]) - set(overlap) == {"COL18A1"}
    assert len(edges) == 543
    assert len(set(edges["#node1"]) | set(edges["node2"])) == 54
    pairs = edges.apply(lambda row: "|".join(sorted((row["#node1"], row["node2"]))), axis=1)
    assert not pairs.duplicated().any()
    assert len(selected) == 5

    frozen_terms = [
        "HIF-1 signaling pathway",
        "inflammatory response",
        "Interleukin-4 and Interleukin-13 signaling",
        "response to oxidative stress",
        "vascular process in circulatory system",
    ]
    assert selected["Description"].tolist() == frozen_terms
    assert selected["Category"].notna().all()
    expected = {
        "IL6": 43, "INS": 42, "IL1B": 39, "IFNG": 39, "MMP9": 39,
        "TP53": 39, "EGFR": 36, "IGF1": 35, "CXCL8": 35, "PTGS2": 35,
        "BCL2": 34, "HIF1A": 33, "ESR1": 32, "CYCS": 28, "MMP2": 28,
    }
    assert dict(zip(top15["name"], top15["Degree"])) == expected

    return overlap, degree, edges, top15, selected


def draw_overlap(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    panel_label(ax, "A")
    ax.text(0.5, 0.96, "Prespecified target overlap", ha="center", va="top",
            fontsize=9, fontweight="bold", color=COLORS["ink"])
    left = Ellipse((0.39, 0.57), .56, .42, facecolor=COLORS["coral"],
                   edgecolor="none", alpha=.30)
    right = Ellipse((0.61, 0.57), .56, .42, facecolor=COLORS["cyan"],
                    edgecolor="none", alpha=.30)
    ax.add_patch(left); ax.add_patch(right)
    ax.text(.26, .84, "C. pilosula targets", ha="center", fontsize=8,
            fontstyle="italic")
    ax.text(.74, .84, "Mountain-sickness\nassociated genes", ha="center",
            va="center", fontsize=8, linespacing=.95)
    ax.text(.22, .57, "344", ha="center", va="center", fontsize=13, fontweight="bold")
    ax.text(.50, .57, "54", ha="center", va="center", fontsize=15,
            fontweight="bold", color=COLORS["ink"])
    ax.text(.78, .57, "297", ha="center", va="center", fontsize=13, fontweight="bold")
    ax.text(.23, .30, "HERB / TCMSP", ha="center", fontsize=7, color=COLORS["muted"])
    ax.text(.77, .30, "GeneCards / OMIM", ha="center", fontsize=7, color=COLORS["muted"])
    ax.text(.50, .19, "398 putative targets        351 associated genes",
            ha="center", fontsize=6.7, color=COLORS["muted"])


def draw_network(ax, degree, edges, top15):
    graph = nx.Graph()
    for _, row in edges.iterrows():
        graph.add_edge(row["#node1"], row["node2"],
                       weight=float(row["combined_score"]))
    assert set(graph.nodes) == set(degree["name"])
    degrees = degree.set_index("name")["Degree"].to_dict()
    nx.set_node_attributes(graph, degrees, "frozen_degree")
    pos = nx.spring_layout(graph, seed=20260923, weight="weight", k=.68,
                           iterations=900, threshold=1e-5)

    weights = np.array([graph[u][v]["weight"] for u, v in graph.edges()])
    widths = .16 + .48 * (weights - weights.min()) / (weights.max() - weights.min())
    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color="#7D878F",
                           width=widths, alpha=.10)
    hubs = set(top15["name"])
    regular = [node for node in graph if node not in hubs]
    hub_nodes = [node for node in graph if node in hubs]
    nx.draw_networkx_nodes(graph, pos, ax=ax, nodelist=regular,
                           node_size=[18 + degrees[n] * 1.45 for n in regular],
                           node_color="#A8D1D2", edgecolors="white", linewidths=.35)
    hub_colors = [COLORS["burgundy"] if degrees[n] >= 39 else COLORS["coral"]
                  for n in hub_nodes]
    nx.draw_networkx_nodes(graph, pos, ax=ax, nodelist=hub_nodes,
                           node_size=[30 + degrees[n] * 2.35 for n in hub_nodes],
                           node_color=hub_colors, edgecolors="white", linewidths=.55)

    labeled_hubs = top15.head(8)["name"].tolist()
    labeled_hubs = sorted(labeled_hubs, key=lambda node: pos[node][0])
    label_groups = (labeled_hubs[:4], labeled_hubs[4:])
    all_x = np.array([xy[0] for xy in pos.values()])
    all_y = np.array([xy[1] for xy in pos.values()])
    x_span = all_x.max() - all_x.min()
    y_span = all_y.max() - all_y.min()
    for side, group in zip(("left", "right"), label_groups):
        group = sorted(group, key=lambda node: pos[node][1])
        node_y = np.array([pos[node][1] for node in group])
        targets = np.linspace(node_y.min() - .055 * y_span,
                              node_y.max() + .055 * y_span, len(group))
        x_target = (min(pos[node][0] for node in labeled_hubs) - .18 * x_span
                    if side == "left" else
                    max(pos[node][0] for node in labeled_hubs) + .18 * x_span)
        for node, target_y in zip(group, targets):
            ax.annotate(node, pos[node], xytext=(x_target, target_y),
                        textcoords="data", ha="right" if side == "left" else "left",
                        va="center", fontsize=5.7, fontstyle="italic",
                        color=COLORS["ink"], zorder=5,
                        arrowprops=dict(arrowstyle="-", color="#9EA5AA",
                                        linewidth=.35, shrinkA=1, shrinkB=2),
                        bbox=dict(boxstyle="round,pad=.10", facecolor="white",
                                  edgecolor="none", alpha=.88))
    ax.set_title("STRING protein-protein interaction network", loc="left",
                 fontsize=9, fontweight="bold", pad=2)
    ax.text(.99, .01, "54 nodes | 543 interactions\nTop 15 highlighted; top 8 labeled",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.3,
            color=COLORS["muted"])
    ax.margins(.08)
    ax.axis("off"); panel_label(ax, "B")


def draw_degree(ax, top15):
    shown = top15.iloc[::-1].copy()
    maximum = shown["Degree"].max()
    cmap = LinearSegmentedColormap.from_list("degree", ["#A8D1D2", COLORS["coral"], COLORS["burgundy"]])
    norm = Normalize(shown["Degree"].min(), maximum)
    colors = [cmap(norm(value)) for value in shown["Degree"]]
    bars = ax.barh(np.arange(len(shown)), shown["Degree"], color=colors,
                   height=.68, edgecolor="none", zorder=2)
    ax.set_yticks(np.arange(len(shown)), shown["name"], fontstyle="italic")
    ax.set_xlim(0, 47)
    ax.set_xlabel("Degree")
    ax.set_title("Highest-connectivity targets", loc="left", fontsize=9,
                 fontweight="bold", pad=4)
    for bar, value in zip(bars, shown["Degree"]):
        ax.text(value + .65, bar.get_y() + bar.get_height() / 2, f"{int(value)}",
                va="center", ha="left", fontsize=6.5, color=COLORS["ink"])
    clean_axis(ax, grid=True); panel_label(ax, "C")


def draw_enrichment(ax, selected):
    display = selected.copy()
    display["label"] = display["Description"].replace({
        "Interleukin-4 and Interleukin-13 signaling": "Interleukin-4 and\nInterleukin-13 signaling",
        "vascular process in circulatory system": "vascular process in\ncirculatory system",
    })
    display = display.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    category_colors = {
        "KEGG Pathway": COLORS["teal"],
        "GO Biological Processes": COLORS["coral"],
        "Reactome Gene Sets": COLORS["ochre"],
    }
    ax.hlines(y, 0, display["minus_log10_q"], color="#CCD1D4", linewidth=.7, zorder=1)
    sizes = 20 + display["listed_gene_count"].to_numpy() * 3.0
    for idx, row in display.iterrows():
        ax.scatter(row["minus_log10_q"], idx, s=sizes[idx],
                   color=category_colors[row["Category"]], edgecolor="white",
                   linewidth=.55, zorder=3)
        ax.text(row["minus_log10_q"] + .45, idx,
                f"{row['minus_log10_q']:.1f}  |  {int(row['listed_gene_count'])} genes",
                va="center", fontsize=6.3, color=COLORS["muted"])
    ax.set_yticks(y, display["label"])
    ax.set_xlim(0, 28)
    ax.set_xlabel("-log10(q-value)")
    ax.set_title("Selected functional enrichment terms", loc="left", fontsize=9,
                 fontweight="bold", pad=4)
    clean_axis(ax, grid=True); panel_label(ax, "D")
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=5,
                          markerfacecolor=color, markeredgecolor="white", label=label)
               for label, color in [("KEGG", COLORS["teal"]),
                                    ("GO biological process", COLORS["coral"]),
                                    ("Reactome", COLORS["ochre"])]]
    ax.legend(handles=handles, frameon=False, loc="upper left",
              bbox_to_anchor=(0, -.16), ncol=3, handletextpad=.30,
              columnspacing=.8, borderaxespad=0)


def save(fig):
    stem = OUT / "Fig1_final_candidate"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches=None,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    if shutil.which("mogrify"):
        for extension in (".png", ".tiff"):
            import subprocess
            subprocess.run(["mogrify", "-units", "PixelsPerInch", "-density", "600",
                            str(stem.with_suffix(extension))], check=True)


def main():
    configure_fonts()
    _, degree, edges, top15, selected = load_frozen_inputs()
    fig = plt.figure(figsize=(180 / 25.4, 7.05), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=(.86, 1.14), height_ratios=(.92, 1.08),
                            left=.08, right=.985, top=.965, bottom=.105,
                            wspace=.30, hspace=.34)
    draw_overlap(fig.add_subplot(grid[0, 0]))
    draw_network(fig.add_subplot(grid[0, 1]), degree, edges, top15)
    draw_degree(fig.add_subplot(grid[1, 0]), top15)
    draw_enrichment(fig.add_subplot(grid[1, 1]), selected)
    save(fig)


if __name__ == "__main__":
    main()
