import numpy as np

# import all parameter sets from PyDREAM and pull out the unique ones, along with their counts
pvals_all = np.load('/Users/leonardharris/PycharmProjects/NERM/param_files/necro_smallest_dreamzs5720_5chain_sampledparams_0_50000.npy')[25000:,:]
pvals_unique, counts = np.unique(pvals_all, return_counts=True, axis=0)

# import the parameter sets for the mode you're interested in
pvals_modeX = np.unique(pvals_all, axis=0)

#######################################################
# mix up the order of the array to test out the code
#######################################################
max_idx = np.argsort(counts)
pvals_modeX = np.array([pvals_modeX[i] for i in max_idx])
#######################################################

mode_counts = []
for p in pvals_modeX:
    for i,u in enumerate(pvals_unique):
        if (p == u).all():
            mode_counts.append(counts[i])
            # print(i)
            continue

print(mode_counts)
