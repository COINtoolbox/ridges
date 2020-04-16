import numpy as np
from astropy.io import fits
import OT_distance as ot

data_path = "Data/TransportPlans/"

def wass_from_T(T,C):
    """ Convert transport plan into Wasserstein distance."""
    return np.inner(C.flatten(),T.flatten())

# OT/Sinkhorn hyperparameters. These should match those used when
# running the "Section3_1_ComputeTransportPlans.py" script.
gamma=.001
n_iter=500
# cost matrix
ridges = np.load('Data/ridges_mask.npy')
C = ot.EuclidCost(ridges, timeit=True)

# Load converged transport plans
T_n2n = np.load(data_path+'maskT_n2n_{}.npy'.format(n_iter))

# Compute Wasserstein distances from converged transport plan
W_n2n = wass_from_T(T_n2n,C)

# Load transport plans for different random realizations, computed from 
# Section3_1_ComputeTransportPlans.py
N2Reals = []
for RANDOM_IDX in range(100):
    T = np.load(data_path+'maskT_n2r_{}_morandum.npy'.format(RANDOM_IDX))
    N2Reals += [wass_from_T(T,C)]
N2Reals = np.array(N2Reals)

# Dump random distribution to csv and print "noisy-to-noiseless"
# Wasserstein distance - the two make up the contents of figure 1
np.savetxt('RawValuesPlots/Fig1distToRand.csv', N2Reals.reshape(1,-1), 
           delimiter = ',')
print(W_n2n)

