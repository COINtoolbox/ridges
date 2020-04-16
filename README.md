# voids
Scripts to reproduce results from Moews et al., 2020

## Contents
### Section 3.2
The script `Section3_2.py` contains the code to reproduce the experiments of Section 3.2.

The curvelet-denoised map is stored in `Data/denois_nonbin_fdr_nocoarse.fits`. It can be recomputed by running:

```mr_filter -t28 -C2 -K des_nonbin.fits denois_nonbin_fdr_nocoarse```

after installing the sparse2d module of [ISAP](http://www.cosmostat.org/software/isap). `-t` selects the type of sparse representation, with `28` being curvelets; `-C2` means the choice of threshold values, for denoising, should be done using False Discovery Rate; `-K` means the coarse scale will not be readded to the reconstructed map.
