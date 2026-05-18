import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

cm = 1 / 2.54

T = np.array(
    [
        0.0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.45,
        0.48,
        0.49,
        0.5,
        0.51,
        0.52,
        0.55,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
    ]
)
data10 = np.array(
    [
        2.27157,
        1.96257,
        1.62786,
        1.25086,
        0.853857,
        0.626286,
        0.39,
        0.256,
        0.02429,
        0.437,
        0.5461,
        0.759429,
        0.977429,
        1.27986,
        1.57057,
        1.83086,
        2.06343,
    ]
)
data09 = np.array(
    [
        1.88364,
        1.63614,
        1.32607,
        1.03143,
        0.704857,
        0.566553,
        0.422158,
        0.223652,
        0.0,
        0.136273,
        0.344526,
        0.644868,
        0.771071,
        1.04371,
        1.31307,
        1.52686,
        1.71857,
    ]
)
data08 = np.array(
    [
        1.54606,
        1.2896,
        1.08834,
        0.893943,
        0.725086,
        0.584368,
        0.433868,
        0.319895,
        0.0,
        0.109263,
        0.335105,
        0.500632,
        0.604229,
        0.852171,
        1.06377,
        1.21331,
        1.4256,
    ]
)
data07 = np.array(
    [
        1.23164,
        1.09414,
        0.9565,
        0.817286,
        0.642786,
        0.569763,
        0.377737,
        0.323395,
        0.0,
        0.142342,
        0.290895,
        0.405711,
        0.550714,
        0.824929,
        0.960929,
        1.10814,
        1.188,
    ]
)
data06 = np.array(
    [
        0.96646,
        0.838057,
        0.647086,
        0.445143,
        0.244686,
        0.156857,
        0.128721,
        0.097167,
        0.0,
        0.138764,
        0.162857,
        0.238743,
        0.369029,
        0.643429,
        0.831257,
        0.959886,
        1.0884,
    ]
)

mask = (T <= 1.0) & (T >= 0.0)

T = T[mask]
data10 = data10[mask]
data09 = data09[mask]
data08 = data08[mask]
data07 = data07[mask]
data06 = data06[mask]

all_data = [data10, data09, data08, data07, data06]

plt.rcParams["font.size"] = 11
plt.rcParams["font.family"] = "Times New Roman"

styles = [
    {"line": "-", "marker": "o", "color": "black", "label": "c = 1.0"},
    {"line": "--", "marker": "s", "color": "black", "label": "c = 0.9"},
    {"line": ":", "marker": "^", "color": "black", "label": "c = 0.8"},
    {"line": "-.", "marker": "x", "color": "black", "label": "c = 0.7"},
    {
        "line": (0, (8, 2, 1, 2, 1, 2, 1, 2)),
        "marker": "D",
        "color": "black",
        "label": "c = 0.6",
    },
]

fig, ax = plt.subplots(figsize=(17 * cm, 10 * cm))

ax.set_xlim(T.min(), T.max())
ax.set_ylim(0, 2.5)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

ax.grid(True, linestyle="--", alpha=0.3, color="black")
ax.tick_params(direction="out", length=5, width=1, colors="black")

ax.set_ylabel(r"Критическая температура $T_N$")
ax.set_xlabel(r"Фрустрационный параметр $r = \frac{J_2}{J_1}$")

for data, style in zip(all_data, styles):
    ax.plot(
        T,
        data,
        linestyle=style["line"],
        color=style["color"],
        marker=style["marker"],
        markersize=4,
        linewidth=1.5,
        label=style["label"],
    )

ax.legend(loc="best", fontsize=11, frameon=True, fancybox=True)

plt.tight_layout()
plt.savefig("results.png", dpi=300, bbox_inches="tight")
