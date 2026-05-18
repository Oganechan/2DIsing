import pandas as pd
import matplotlib.pyplot as plt

pos = pd.read_csv("atom_positions.csv")
plt.scatter(pos['x'], pos['y'])
plt.savefig("test_positions.png")