# voids
Scripts to reproduce results from Moews et al., 2020

## Contents
### Section 3.1
The script `Section3_1.py` contains the code to compute the Wasserstein distances of Section 3.1's experiment. It relies on the outputs of `Section3_1_ComputeTransportPlans.py`, which computes the full transport plans for every random map realization. These are fairly large files, and are therefore not included in this repo.

Both of these scripts import `OT_distance.py`, which contains some ad hoc optimal transport-related, pure-python functions. If you are only interested in computing optimal transport quantities, beyond the scope of this paper, we strongly recommend you look into other, maintained ressources instead of these (for instance, [POT](https://github.com/rflamary/POT)).

### Section 3.2
The script `Section3_2.py` contains the code to reproduce the experiments of Section 3.2.

The curvelet-denoised map is stored in `Data/denois_nonbin_fdr_nocoarse.fits`. It can be recomputed by running:

```mr_filter -t28 -C2 -K des_nonbin.fits denois_nonbin_fdr_nocoarse```

after installing the sparse2d module of [ISAP](http://www.cosmostat.org/software/isap). `-t` selects the type of sparse representation, with `28` being curvelets; `-C2` means the choice of threshold values, for denoising, should be chosen using False Discovery Rate; `-K` means the coarse scale will not be readded to the reconstructed map.
