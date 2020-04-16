""" This script performs downsampling of both the data described in section 2.3
and the ridges outputted by SCMS, for use in the experiments of Section 3.
For Section 3.1, we apply the (projected) DES mask and ensure there is the same number 
of "mesh points" (meaning active pixels in the low-res maps) in both the noisy and noiseless sims."""
import numpy as np
from astropy.io import fits

data_path = 'Data/'
output_path = 'Data/'

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
        
def main():
    ###################################################################
    ### Downsample ridges obtained from simulations for section 3.1 ###
    ###################################################################
    # load SCMS-outputted ridges
    simridges = np.load(output_path+'sim_noiseless_nonmirror_ridges.npy')
    noisimridges = np.load(output_path+'sim_noisy_ridges.npy')
    # remove weird SCMS artefacts if present
    simridges = simridges[simridges[:,1]>-65]
    noisimridges = noisimridges[noisimridges[:,1]>-65]
    
    # convert them to 2D image for Wasserstein
    grid_size=[296,37]
    simridge_bins = CoordsToIm(simridges, grid_size=grid_size, binary=False)
    noisimridge_bins = CoordsToIm(noisimridges, grid_size=grid_size, binary=False)
    ridges, noisridges = simridge_bins[0], noisimridge_bins[0]
    
    # load projected mask
    binmask = np.load('Data/lores_binmask.npy')
    
    # apply mask
    negmask = np.abs(binmask-1).astype(bool)
    ridges[negmask] = 0
    noisridges[negmask] = 0
    
    # Count activated pixels - WARNING
    act = np.sum(ridges > 0)
    noisact = np.sum(noisridges > 0)

    print(act, noisact) # WARNING: what follows assumes act<noisact
    
    # And create final, masked maps:
    im = np.copy(ridges)
    im[im>0] = 1
    nois_im = binarize_map(noisridges, count=act, binary=True)
    
    # save downsampled maps
    np.save(data_path+'ridges_mask.npy', im)
    np.save(data_path+'noisridges_mask.npy', nois_im)
    
    #################################################################################
    ### Downsample and convert real DES data to 2D maps to run curvelet denoising ###
    #################################################################################
    grid_size=[600,75] 
    
    # Load samples from real DES mass map
    desSamps = fits.open(data_path+'des-masked-noisy.fits')
    ndat = desSamps[1].data

    # Convert coordinates
    ra = ndat['ra']
    dec = ndat['dec']
    ra[ra>180] -= 360
    coords = np.array([ra,dec]).T
    
    # Convert to 2D map
    bins = CoordsToIm(coords, binary=False, grid_size=grid_size)
    
    # save
    fits.writeto(output_path+'des_nonbin.fits', bins[0], overwrite=True)
    
    
if __name__ == "__main__":
    main()
