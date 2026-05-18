import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse
import glob
import os


def read_snapshot(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    temp = 0
    conc = 0
    frust = 0
    for line in lines:
        if "Temperature:" in line:
            temp = float(line.split(":")[1])
        elif "Concentration:" in line:
            conc = float(line.split(":")[1])
        elif "Frustration parameter:" in line:
            frust = float(line.split(":")[1])

    data = []
    for line in lines:
        if not line.startswith("#") and line.strip():
            parts = line.split()
            if len(parts) >= 4:
                x, y, spin, magnetic = parts
                data.append([int(x), int(y), int(spin), int(magnetic)])

    data = np.array(data)

    max_x = int(data[:, 0].max()) + 1
    max_y = int(data[:, 1].max()) + 1

    lattice = np.zeros((max_x, max_y), dtype=int)

    for x, y, spin, magnetic in data:
        if magnetic:
            lattice[x, y] = spin

    return temp, conc, frust, lattice


parser = argparse.ArgumentParser()
parser.add_argument("dirpath", type=str)
parser.add_argument("-o", "--output", type=str, default="animation.gif")
args = parser.parse_args()

all_files = sorted(glob.glob(os.path.join(args.dirpath, "snapshot*.txt")))

temps = []
concs = []
frusts = []
spins_sequence = []
for f in all_files:
    temp, conc, frust, lattice = read_snapshot(f)
    temps.append(temp)
    concs.append(conc)
    frusts.append(frust)
    spins_sequence.append(lattice)

fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(
    spins_sequence[0], cmap="binary", vmin=-1, vmax=1, interpolation="nearest"
)
ax.set_title(
    f"T = {temps[0]:.3f}, c = {concs[0]:.2f}, r = {frusts[0]:.2f}", fontsize=24
)
ax.axis("off")


def update(frame):
    im.set_array(spins_sequence[frame])
    ax.set_title(
        f"T = {temps[frame]:.3f}, c = {concs[frame]:.2f}, r = {frusts[frame]:.2f}",
        fontsize=24,
    )
    return [im]


anim = FuncAnimation(fig, update, frames=len(spins_sequence), interval=100, blit=True)

anim.save(args.output, writer="pillow", fps=10)
