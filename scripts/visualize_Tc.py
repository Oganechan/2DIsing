import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

cm = 1 / 2.54

parser = argparse.ArgumentParser()
parser.add_argument("filepath", type=str)
parser.add_argument("--output", "-o", type=str, default="results.png")

args = parser.parse_args()

plt.rcParams["font.size"] = 11
plt.rcParams["font.family"] = "Times New Roman"

data = pd.read_csv(args.filepath)

r = data["r"]
Tc = data["Tc"]

fig, ax = plt.subplots(1, 1, figsize=(17 * cm, 9 * cm), sharex=True)

ax.scatter(r, Tc, c="black")
ax.plot(r, Tc, c="black")

ax.grid(True, linestyle="--", alpha=0.3, color="black")
ax.set_xlabel(r"Фрустрационный параметр $r =  \frac{J_2}{J_1}$")
ax.set_ylabel(r"Критическая температура $T_N$")
# ax.set_xlim(0.0, 1.0)
# ax.set_ylim(0.0, 2.5)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

plt.tight_layout()
plt.savefig(args.output, dpi=300, bbox_inches="tight")
