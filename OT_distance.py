import numpy as np
import time



#############
# UTILITIES #
#############
def EuclidCost(data, divmed=False, timeit=True, trunc=False, maxtol=745.13, truncval=745.13):
    ''' Function to compute usual Euclidean distance matrix;
    assumes both images will be flattened in a np.ndarray.flatten() fashion
   
    ARGUMENTS:
        data: The input image
        divmed: (bool) if True, normalize cost matrix by median distance value (gives gamma
                    values more of a data-independant meaning)
        timeit: (bool) if True, prints the time spent computing the cost matrix
        trunc: (bool) if True, truncates cost matrix entries above a given value
        maxtol: (float) if trunc is True, value above which distances are truncated; default
                    is maximum value that avoids naught entries in Sinkhorn kernel for gamma=1
        truncval: (float) if trunc is True, value to which truncated entries are set; default
                    is maximum value that avoids naught entries in Sinkhorn kernel for gamma=1
               
    OUTPUT:
        C: (Nr*Nc x Nr*Nc np.ndarray) Distance matrix for the Euclidean distance'''
    if timeit:
        start = time.time()
    s = np.shape(data)
    N = np.size(data)
    loc = np.where(np.zeros(s)==0)
    loc = np.array(loc).T.astype(float)
    C = np.sum((loc[:,np.newaxis,:]/s[1]-loc[np.newaxis,:,:]/s[1])**2, axis=2)
    
    C *= (s[0]-1)**2*2 / np.max(C)
   
    if timeit:
        print('cost matrix computed in '+str(time.time()-start)+'s.')
    if divmed:
        C /= np.median(C)
    if trunc:
        C[C>maxtol] = truncval
    return C
        
def safelog(M, epsilon=1e-200):
    M[M==0] = epsilon
    return np.log(M)
    
def stabker(C,alpha,beta,gamma):
    ''' Function to use the stable -c_ij + alpha_i + beta_j quantity for stabilized Sinkhorn;'''
    N1, N2 = alpha.shape[0], beta.shape[0]
    return np.exp((-C + np.outer(alpha, np.ones(N1)) + np.outer(np.ones(N2),beta))/gamma)
    
##################
# TRANSPORT PLAN #
##################
def sinkhorn(mu, nu, gamma=10., C=None, n_iter=100):
    ''' Compute transport plan between input images;
    
    ARGUMENTS:
        mu: (Nx1 np.ndarray) input measure (e.g. observed lensed ark)
        nu: (Nx1 np.ndarray) input measure (e.g. current unlensed source estimate)
        gamma: (float) entropic penalty parameter; the closer to 0, the sharper the transport plan
                    (and the closer it is to the true optimal),
                    but the higher the risk of numerical instability
        C: (None or NxN np.ndarray) cost matrix; if None, computes usual Euclidean distance
                    assuming marginals are square images having been flattened in 
                    np.ndarray.flatten() fashion
        n_iter: (int) number of Sinkhorn iterations
        
    OUTPUTS:
        T: (NxN np.ndarray) approximate optimal transport plan between mu and nu, i.e. we 
                    have approximately:
        T.dot(np.ones(N)) = mu
        np.ones(N)[:,np.newaxis].dot(T) = nu
    '''
    N = mu.size
    p = mu.flatten()
    q = nu.flatten()
    # initialize kernel and scaling vectors
    if C is None:
        K = np.exp(-EuclidCost(mu)/gamma)
    else:
        K = np.exp(-C/gamma)
    a, b = np.ones(N), np.ones(N)
    # Sinkhorn updates
    for j in range(n_iter):
        a = p / K.dot(b)
        b = q / K.T.dot(a)
    T = K*a[:,np.newaxis]*b
    return T
    
    
def logsinkhorn(mu, nu, gamma=10., C=None, n_iter=100):
    ''' Compute transport plan between input images using log-domain stabilization;
    if precision is not a priority, use the (much faster) sinkhorn function with a 
    high gamma value instead
    
    ARGUMENTS:
        mu: (Nx1 np.ndarray) input measure (e.g. observed lensed ark)
        nu: (Nx1 np.ndarray) input measure (e.g. current unlensed source estimate)
        gamma: (float) entropic penalty parameter; the closer to 0, the sharper the transport plan
                    (and the closer it is to the true optimal),
                    but the higher the risk of numerical instability
        C: (None or NxN np.ndarray) cost matrix; if None, computes usual Euclidean distance
                    assuming marginals are square images having been flattened in 
                    np.ndarray.flatten() fashion
        n_iter: (int) number of Sinkhorn iterations
        
    OUTPUTS:
        T: (NxN np.ndarray) approximate optimal transport plan between mu and nu, i.e. we 
                    have approximately:
        T.dot(np.ones(N)) = mu
        np.ones(N)[:,np.newaxis].dot(T) = nu
    '''
    N = mu.size
    p = np.copy(mu.flatten())
    q = np.copy(nu.flatten())
    # add epsilon to naught entries in input measures
    if np.any(mu==0) or np.any(nu==0):
        p[p==0] = 1e-200
        q[q==0] = 1e-200
    # compute cost if needed
    if C is None:
        C = EuclidCost(mu) 
    # initialize dual scaling vectors
    alpha, beta = np.zeros(N), np.zeros(N)
    # log sinkhorn
    for j in range(1,n_iter+1):
        stab = stabker(C,alpha,beta,gamma)
        alpha = gamma * (np.log(p) - safelog(
                        np.sum(stab,axis=1))) + alpha
        stab = stabker(C,alpha,beta,gamma)
        beta = gamma * (np.log(q) - safelog(
                        np.sum(stab, axis=0))) + beta 
    return stabker(C,alpha,beta,gamma)
   
def WassDistance(mu, nu, gamma=10., C=None, n_iter=100):
    ''' Compute approximate Wasserstein distance between mu and nu.'''
    if C is None:
        C = EuclidCost(mu) 
    return np.inner(C.flatten(),sinkhorn(mu, nu, gamma, n_iter, C).flatten())
