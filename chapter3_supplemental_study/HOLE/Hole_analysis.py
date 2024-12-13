import MDAnalysis as mda
from MDAnalysis.analysis import hole2
import matplotlib.pyplot as plt
import pandas as pd
import os

# load model
root_path = '/Users/user/Desktop/HOLE_analysis/'  # adjust as needed
multipdb_file = f'{root_path}/my_models.pdb'  # your PDB file containing the reference and models
cpoint_coords = [0.09, 1.69, -1.68]  # initial center point coordinates. adjust as needed
cvect_input = [0, 0, 1]
glu_center_z = [11.52, 12.05, 12.10, 11.95, 11.76]  # chai-1 selectivity filter center. Modify for your own use
# run HOLE
hole_exec = '/Users/user/anaconda3/envs/hole2/bin/hole'  # path to hole executable
sph_process = '/Users/user/anaconda3/envs/hole2/bin/sph_process'  # path to run sph_process executable

# for multiple PDB models
u = mda.Universe(multipdb_file)

ha = hole2.HoleAnalysis(u, select='protein', executable=hole_exec, cpoint=cpoint_coords, cvect=cvect_input,
                        prefix=f'{root_path}/')
ha.run()

# extract results
# colors = {'black'} # reference
colors = {'#E69F00', '#56B4E9', '#009E73', '#D55E00', '#CC79A7'}  # models
header = ['rxn_coord', 'pore_radius', 'cen_line']
for (key, recarray), color, glu_center_z in zip(ha.results.profiles.items(), colors, glu_center_z):
    # convert to dataframe
    df = pd.DataFrame.from_records(recarray, columns=header)
    df['corrected_rxcn_coord'] = df['rxn_coord'] - glu_center_z
    # save dataframe to csv
    filename = f'{key}.csv'
    outpath = f'{root_path}/pore_profile'
    os.makedirs(outpath, exist_ok=True)
    df.to_csv(f'{outpath}/{key}_profile.csv', index=False)

    # create plot
    df_filtered = df[df['radius'] <= 7]  # adjust as needed
    plt.plot(df_filtered['radius'], df_filtered['corrected_rxcn_coord'], label=key, color=color)

x_min, x_max = 0, 7  # adjust as needed
plt.xlim(x_min, x_max)
x_ticks = [0, 1, 2, 3, 4, 5, 6, 7]  # adjust as needed
plt.xticks(x_ticks)
plt.xlabel('Pore Radius (Å)')

y_min, y_max = -50, 30  # adjust as needed
y_ticks = [-50, -40, -30, -20, -10, 0, 10, 20, 30]  # adjust as needed
plt.ylim(y_min, y_max)
plt.yticks(y_ticks)
plt.ylabel('Transmembrane Z-Axis (Å)')

plt.axvline(x=1.4, color='black', linestyle='--', linewidth=1)  # water molecule radius
plt.savefig(f'{root_path}/pore_profile_plot.png', dpi=500)

# plot results
plt.plot(profile[0].radius, profile[0].rxn_coord, linestyle="-")
plt.ylabel('Pore Coordinate')
plt.xlabel('Pore Radius (Å)')
plt.show()

exit()
