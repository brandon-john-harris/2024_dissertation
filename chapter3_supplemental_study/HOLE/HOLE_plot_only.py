import matplotlib.pyplot as plt
import pandas as pd
import os

root_path = '/Users/yarovoylab/Desktop/HOLE_analysis/folder/pore_profile'  # adjust as needed

# extract results
df_files = {'black': (f'{root_path}/ref_profile.csv', 'Reference'),
          '#E69F00': (f'{root_path}/0_profile.csv', 'Rank 1'),
          '#56B4E9': (f'{root_path}/1_profile.csv', 'Rank 2'),
          '#009E73': (f'{root_path}/2_profile.csv', 'Rank 3'),
          '#D55E00': (f'{root_path}/3_profile.csv', 'Rank 4'),
          '#CC79A7': (f'{root_path}/4_profile.csv', 'Rank 5')}

for color, (df_file, label) in df_files.items():
    # convert to dataframe
    df = pd.read_csv(df_file, sep=',')

    # create plot
    df_filtered = df[df['radius'] <=7]
    plt.plot(df_filtered['radius'], df_filtered['corrected_rxcn_coord'], color=color, label=label)

x_min, x_max = 0, 7
plt.xlim(x_min, x_max)
x_ticks = [0, 1, 2, 3, 4, 5, 6, 7]
plt.xticks(x_ticks)
plt.xlabel('Pore Radius (Å)')

y_min, y_max = -40, 30
y_ticks = [-40, -30, -20, -10, 0, 10, 20, 30]
plt.ylim(y_min, y_max)
plt.yticks(y_ticks)
plt.ylabel('Transmembrane Z-Axis (Å)')

plt.legend()
plt.axvline(x=1.4, color='black', linestyle='--', linewidth=1) # water molecule radius
plt.savefig(f'{root_path}/pore_profile_plot.png', dpi=500)

exit()
