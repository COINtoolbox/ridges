import numpy as np
from astropy.io import fits
import time

import OT_distance as ot

# To get the exact same realizations as those shown in Figure 1, run this script
# 4 times with each of the following seeds:
range_random = range(25)
np.random.seed(1480)

#range_random = range(25,50)
#np.random.seed(248973)

#range_random = range(50,75)
#np.random.seed(387)

#range_random = range(75,100)
#np.random.seed(8745)

# OT/Sinkhorn hyperparameters. A very low gamma avoids excessive blurring of
# the transport plans; and since computation time hardly matters here, take
# a very large number of Sinkhorn iterations to ensure convergence.
gamma=.001
n_iter=500

save_path='Data/TransportPlans/'

ridges = np.load('Data/ridges_mask.npy')
# avoid absolute 0s
ridges[ridges==0] = 1e-10

noisridges = np.load('Data/noisridges_mask.npy')
noisridges[noisridges==0] = 1e-10

# projected mask
binmask = np.load('revisit/lores_binmask.npy')

# compute how many pixels contain ridges
mesh_points = np.sum(ridges>1e-7)  # not really meshpoints but ok
noimesh_points = np.sum(noisridges>1e-7)  # not really meshpoints but ok
unmaskd_pix = int(np.sum(binmask))

# compute cost matrix
C = ot.EuclidCost(ridges, timeit=True)

# normalize to make them histograms
ridges /= np.sum(ridges)
noisridges /= np.sum(noisridges)

# compute transoport plan between ridges and realizations of random maps
for RANDOM_IDX in range_random:
    print(' > NOW WORKING ON RANDOM REAL NUMBER {}'.format(RANDOM_IDX))
    idx = np.random.choice(unmaskd_pix, mesh_points, False)
    
    # create random map
    randomp = np.zeros(ridges.shape) + 1e-10
    randomp[np.where(binmask)[0][idx], np.where(binmask)[1][idx]] = 1

    randomp /= np.sum(randomp)

    # compute transport plan
    start = time.time()
    T_n2r = ot.logsinkhorn(noisridges.flatten(), randomp.flatten(), gamma, C, n_iter)
    np.save(save_path+'maskT_n2r_{}_morandum.npy'.format(RANDOM_IDX), T_n2r)
    print '   > N2R done, TIME ELAPSED:\t{}'.format(time.time()-start)

