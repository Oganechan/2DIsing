import pandas as pd
import matplotlib.pyplot as plt
import argparse

cm = 1 / 2.54

parser = argparse.ArgumentParser()
parser.add_argument("filepath", type=str)
parser.add_argument("--output", "-o", type=str, default="results.png")

args = parser.parse_args()

plt.rcParams["font.size"] = 10
plt.rcParams["font.family"] = "serif"

file = pd.read_csv(args.filepath)

T = file["T"]
Z = file["Z"]
E = file["E"]
C = file["C"]
U = file["U"]
S = file["S"]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 7), sharex=True)

ax1.scatter(T, E, c="black", s=8, linewidths=0)
ax1.set_ylabel("Mean Energy")
ax1.grid(True, linestyle="--", alpha=0.3)
ax1.set_xscale("log", base=10)

ax2.scatter(T, C, c="black", s=8, linewidths=0)
ax2.set_ylabel("Heat capacity")
ax2.grid(True, linestyle="--", alpha=0.3)
ax2.set_xscale("log", base=10)

ax3.scatter(T, S, c="black", s=8, linewidths=0)
ax3.set_ylabel("Entropy")
ax3.grid(True, linestyle="--", alpha=0.3)
ax3.set_xscale("log", base=10)

plt.tight_layout()
plt.savefig(args.output, dpi=300, bbox_inches="tight")
