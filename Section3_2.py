import numpy as np
from astropy.io import fits

data_path = 'Data/'

def CoordsToIm(coords, grid_size=[400,50], mins=[-62,-62], maxs=[100,-37],
               binary=True):
    """ Convert a set of positions to a 2D image with arbitrary number of pixels,
    containing how many positions fell within each."""
    x_bins = np.linspace(mins[0], maxs[0], grid_size[0])
    y_bins = np.linspace(mins[1], maxs[1], grid_size[1])

    ra, dec = coords[:,0], coords[:,1]

    bins = np.histogram2d(dec, ra, bins=[y_bins,x_bins])
    if binary:
        bins[0][bins[0]>0] = 1
    return bins

def actual_to_grid(ra,dec,bins_x,bins_y):
    """Convert ra,dec positions to those on the gridded 2D image. (Used to overplot
    points with actual coordinates on the gridded images, eg the curvelet denoising.)"""
    x = -.5 + (ra - np.min(bins_x)) / (np.max(bins_x)-np.min(bins_x)) * (len(bins_x))
    y = -.5 + (dec - np.min(bins_y)) / (np.max(bins_y)-np.min(bins_y)) * (len(bins_y))
    return [x,y]

def binarize_map(im, count=None, min_val=None, binary=True):
    """ Convert image map to binary, keeping either the `count` highest-valued bins,
    or those containing more than `min_val` meshpoints. One of either should be a
    numerical value; if both are, `min_val` is ignored."""
    if count is None:
        if min_val is None:
            raise ValueError('Either count or min_val should be a numerical value')
        thrsh_im = np.copy(im)
        thrsh_im[thrsh_im < min_val] = 0
        if binary:
            thrsh_im[thrsh_im>0] = 1
        return thrsh_im
    else:
        sortidx = np.argsort(im.flatten())
        # count activated pixels:
        actpix = np.sum(im>0)
        thrsh_im = np.copy(im.flatten())
        thrsh_im[sortidx[:-count]] = 0
        if binary:
            thrsh_im[thrsh_im>0] = 1
        return thrsh_im.reshape(im.shape)

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


