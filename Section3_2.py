import numpy as np
from astropy.io import fits

from Section3_Downsampling import CoordsToIm, actual_to_grid, binarize_map

data_path = 'Data/'



# Load ridges
ridgecoord = np.load(data_path+'des_real_ridges.npy')

# remove weird SCMS artefacts if present
ridges = ridgecoord[ridgecoord[:,1]>-65]

# Load curvelet denoising
denois = fits.open(data_path+'denois_nonbin_fdr_nocoarse.fits')
denois_map = denois[0].data

fake_bins = CoordsToIm(np.random.rand(487,2), denois_map.shape[::-1])
grid_ridge = np.array([actual_to_grid(ra,dec,fake_bins[2],fake_bins[1]) for ra,dec in ridgecoord])

# binarize ridges
grid_size=[600, 75]
ridges = CoordsToIm(ridgecoord, grid_size=grid_size, binary=False)[0]

# Compute mismatch per threshold on the ridges (Figure 4)
fullMis = []

thresholds = [1, 5, 10, 20, 50] # in meshpoints per pixel
for minval in thresholds:

    im = np.copy(ridges)
    im[im<minval] = 0
    im[im>=minval] = 1

    # load and apply mask (note you will need to reproject the mask if
    # you change the arbitrary pixel grid)
    mask = np.load(data_path + 'hires_binmask.npy')
    mask = (1-mask).astype(bool)
    im[mask] = 0
    denois_map[mask] = 0

    # Amplitude of the curvelet map at each mesh point
    curveletEntries = np.array([denois_map[x,y] for x,y in zip(np.where(im)[0], np.where(im)[1])])

    # Compute how much of the ridges fall in regions of the curvelet map above a certain value
    curvValGrid = np.linspace(np.max(denois_map), 0, 500) # linear in the span of amplitudes of the curvelet map
    countInside = np.array([np.sum(curveletEntries >= val) for val in curvValGrid])

    # proportion of ridges that "match" curvelet for each of these values
    normInside = countInside / np.sum(im) 
    
    # "full mismatch" are places where ridges are at hard 0 areas in curvelet map
    fullMis += [1- normInside[-2]]
    
# Convert to percentage
fullMis = np.array(fullMis) * 100

# Convert thresholds to meshpoints per sq degree
raStep = (fake_bins[2][1:] - fake_bins[2][:-1])[0]
decStep = (fake_bins[1][1:] - fake_bins[1][:-1])[0]

pixArea = raStep * decStep # square degrees

# Same, per bandwidth:
bw_fullMis = []
bandwidth_names = [None, '1_25', '1_5', '2_0']
for bw in bandwidth_names:
    # Load ridges
    if bw is None: # bandwidth determined as described in Sect 2.2.1
        ridge_fn = data_path+'des_real_ridges'
    else: # multiply that cross-validated bandwidth by different values
        ridge_fn = data_path+'DifferentBW/des_real_ridges_{}'.format(bw)
    ridge_fn += '.npy'
    ridgecoord = np.load(ridge_fn)
    # remove weird SCMS artefacts if present
    ridges = ridgecoord[ridgecoord[:,1]>-65]
    
    # convert to arbitrary 2D map
    ridges = CoordsToIm(ridgecoord, grid_size=grid_size, binary=False)[0]

    # no increasing the threshold here
    minval=1

    im = np.copy(ridges)
    im[im<minval] = 0
    im[im>=minval] = 1

    # reapply mask 
    im[mask] = 0

    # and same as before:
    curveletEntries = np.array([denois_map[x,y] for x,y in zip(np.where(im)[0], np.where(im)[1])])

    curvValGrid = np.linspace(np.max(denois_map), 0, 500)
    countInside = np.array([np.sum(curveletEntries >= val) for val in curvValGrid])

    normInside = countInside / np.sum(im)
    
    bw_fullMis += [1- normInside[-2]]
    
# Convert to percentage
bw_fullMis = np.array(bw_fullMis) * 100

# Dump values used in figure 4 to .csv:
threshValues = np.vstack((thresholds / pixArea, fullMis))
np.savetxt('RawValuesPlots/Fig4ThresholdCurve.csv', threshValues, delimiter=',')

BWValues = np.vstack(([1, 1.25, 1.5, 2.], bw_fullMis))
np.savetxt('RawValuesPlots/Fig4BWCurve.csv', BWValues, delimiter=',')

# And those used in Figure 3 to .fits
redVsGreen = np.zeros((denois_map.shape))
for x,y in zip(np.where(im)[0], np.where(im)[1]):
    curvEntry = denois_map[x,y]
    if curvEntry:
        redVsGreen[x,y] = 1
    else:
        redVsGreen[x,y] = 2
        
fits.writeto('RawValuesPlots/Fig3.fits', redVsGreen)


