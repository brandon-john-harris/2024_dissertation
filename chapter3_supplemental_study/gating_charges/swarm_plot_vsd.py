import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import numpy as np
csv_input = '/Users/user/Data/gating_charge.csv'
# csv file contains headers Category, Subcategory, Value, Rank
# Category: VSD numbering (e.g. VSD1)
# Subcategory: Gating charge identifier (e.g. 1)
# Value: the gating charge distance, with negative indicating below hydrophobic constriction site
# Rank: The model rank assigned by Chai-1

df = pd.read_csv(csv_input, sep=',')
# df['Subcategory'] = df['Subcategory'].str.strip().str.title()
rank_markers = {'1': '.', '2': 0, '3': 1, '4': 0, '5': 1, 'Reference': '_'}  # for gating charge plot
# subcategory_colors = {'1': 'red', '2': 'orange', '3': 'yellow', '4': 'green', '5': 'blue', '6': 'purple'}
subcategory_colors = {1: '#E69F00', 2: '#56B4E9', 3: '#009E73', 4: '#F0E442', 5: '#0072B2',
                      6: '#D55E00'}  # for gating charge plot
for (subcategory, rank), subset in df.groupby(['Subcategory', 'Rank']):
    marker = rank_markers[rank]
    color = 'black' if rank == 'Reference' else subcategory_colors[subcategory]
    jittered_x = pd.factorize(subset['Category'])[0] + np.random.uniform(-0.0, 0.0, len(subset))
    plt.scatter(x=jittered_x, y=subset['Value'], label=f'Rank {rank}', marker=marker, color=color, s=50, alpha=1)

plt.xticks(range(len(df['Category'].unique())), labels=df['Category'].unique())
plt.xlim(-0.5, 3.5)  # adjust as needed

y_min, y_max = -20, 25   # adjust as needed
y_ticks = [-20, -15, -10, -5, 0, 5, 10, 15, 20, 25]   # adjust as needed
plt.ylim(y_min, y_max)
plt.yticks(y_ticks)
plt.ylabel('Cα-Cα Distance to HCS (Å)')

# plt.subplots_adjust(left=0, right=1, top=0.9, bottom=0.2) # when you need to remove white space within the plot
# plt.show()
plt.savefig('/Users/user/Data/gating_charge_example.png', dpi=500)
exit()
