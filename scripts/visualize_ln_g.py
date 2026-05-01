import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'serif'

file = pd.read_csv('test_wl_ln_g.csv')

key = file['key']
E = file['E']
ln_g = file['ln_g']

fig, ax = plt.subplots(1,1,figsize=(9,7), sharex=True)

ax.scatter(E, ln_g, c='black', s=8, linewidths=0)
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('test_wl_ln_g.png')