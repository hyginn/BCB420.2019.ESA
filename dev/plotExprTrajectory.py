# imports the pcreode package
import pcreode
# pandas is a package used for making the handling of large data sets easier
import pandas as pd
# numpy is very common package for handling arrays and matrices
import numpy as np
from sklearn.metrics import pairwise_distances
from igraph import *
import matplotlib
from IPython.display import display, Image
from sklearn import preprocessing

import sys

file_nm = sys.argv[1]

data_raw = pd.read_csv( file_nm)
overlay = data_raw[sys.argv[2]]
data_pca = pcreode.PCA( data_raw)
data_pca.get_pca()
pca_reduced_data = data_pca.pca_set_components( 3)

dens = pcreode.Density( pca_reduced_data)

guess = dens.radius_best_guess()
density_1 = dens.get_density( radius=guess)

noise = 8.0
target = 25.0
num_runs = 3

downed, downed_ind = pcreode.Down_Sample( pca_reduced_data, density_1, noise, target)
out_graph, out_ids = pcreode.pCreode( data=pca_reduced_data, density=density_1, noise=noise,
                                      target=target, file_path="./data/temp/", num_runs=num_runs)

graph_ranks = pcreode.pCreode_Scoring( data=pca_reduced_data, file_path="./data/temp/", num_graphs=num_runs)

gid = graph_ranks[0]

analysis = pcreode.Analysis( file_path="./data/temp/", graph_id=gid, data=pca_reduced_data, density=density_1, noise=noise)

norm_dens = preprocessing.MinMaxScaler( feature_range=(2,8))
dens = norm_dens.fit_transform( analysis._density.astype( float).reshape(-1, 1))[analysis.node_data_indices]

seed = 1

upper_range=3
node_label_size=0

w_adj = pcreode.return_weighted_adj( pca_reduced_data, "./data/temp/", gid)

np.savetxt("./data/temp/weights.txt", w_adj)
np.savetxt("./data/temp/density.txt", dens)
np.savetxt("./data/temp/nodes.txt", analysis.node_data_indices)

norm_ana = preprocessing.MinMaxScaler(feature_range=(0, upper_range))
norm_ana.fit(overlay[analysis._density > analysis._noise].values.astype(np.float).reshape(-1, 1))
old_ana = norm_ana.transform(overlay[analysis._density > analysis._noise].values.astype(np.float).reshape(-1, 1))
# bin the data points to each node so that an average of closest surrounding nodes is used for overlay
bin_dist = pairwise_distances(analysis.good_cells, analysis._data[analysis.node_data_indices])
bin_assignments = np.argmin(bin_dist, axis=1)
new_ana = overlay.values[analysis.node_data_indices]
for ii in range(analysis.num_nodes):
    new_ana[ii] = np.mean(old_ana[bin_assignments == ii])
norm_1 = np.array(new_ana, dtype=float)
cl_vals_1 = [[]] * analysis.num_nodes
# colors to use for overlay
get_cl = matplotlib.cm.get_cmap('RdYlBu_r')
for jj in range(analysis.num_nodes):
    cl_vals_1[jj] = get_cl(norm_1[jj])

analysis.graph.vs["color"] = [cl_vals_1[kk] for kk in range(analysis.num_nodes)]

np.savetxt("./data/temp/colors.txt", analysis.graph.vs["color"])
