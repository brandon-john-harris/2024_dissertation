import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

csv_input= '/Users/user/ligand_rmsd.csv'  # data with headers Category Value, Rank
# Category: the PDB code (e.g. 5ek0)
# Value: the RMSD
# Rank: the Chai-1 value-assigned rank

df = pd.read_csv(csv_input, sep=',')
rank_colors = {1: '#E69F00', 2: '#56B4E9', 3: '#009E73', 4: '#F0E442', 5: '#0072B2'}  # adjust as needed

for rank, subset in df.groupby(['Rank']):
    color = rank_colors[rank]
    jittered_x = pd.factorize(subset['Category'])[0] + np.random.uniform(-0.0, 0.0, len(subset))
    plt.scatter(x=jittered_x, y=subset['Value'], label=f'Rank {rank}', color=color, marker='.', s=50, alpha=0.5)

plt.xticks(range(len(df['Category'].unique())), labels=df['Category'].unique(), rotation=-15)
# plt.xlim(-0.5, 15)

y_min, y_max = 0, 50  # adjust as needed
y_ticks = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]  # adjust as needed
plt.ylim(y_min, y_max)
plt.axhline(y=2, color='black', linestyle='--', linewidth=1)
plt.yticks(y_ticks)
plt.ylabel('Chai-1 Ligand RMSD (Å) to Structure')

plt.legend()

# plt.subplots_adjust(left=0, right=1, top=0.9, bottom=0.2) # when you need to remove white space within the plot
plt.savefig('/Users/user/save_path/ligand_rmsd_example.png', dpi=500)  # adjust as needed
exit()
