[![crp6](https://img.shields.io/badge/CRP-%236-blue)](https://cosmostatistics-initiative.org/residence-programs/crp6/)
[![arXiv](https://img.shields.io/badge/arXiv-astro--ph%2F2005.08583-red)](https://arxiv.org/pdf/2005.08583.pdf) 

# <img align="right" src="images/COIN_logo_very_small.png" width="350"> Ridges in the Dark Energy Survey for cosmic trough identification

This repository holds the scripts necessary to reproduce results from [Moews et al., 2020](https://arxiv.org/abs/2005.08583).

We kindly ask you to include the [full citation](https://ui.adsabs.harvard.edu/abs/2020arXiv200508583M/exportcitation) if you use this material in your research.

## Introduction


The ridges shown in the paper were obtained by running the `filaments` function of the modified DREDGE package, with arguments:

```ridges = filaments(data, n_process = n_proc, mesh_size = 100000, convergence=.1)```

Where `data` contains is an `(n, 2)` ndarray containing 2D positions of samples drawn as described in Section 2.3.2. In the paper, these are either a random subset of those in `Data/des-masked-noisy.fits`, or from the simulatons described in Section 2.3.3 for the experiments of Section 3.1. `n_proc` is the number of parallel processes to run.

## Contents
### Preprocessing
Experiments in both Section 3.1 and 3.2 require the conversion of input data and/or ridges to 2D maps. This preprocessing step is performed by `Section3_Downsampling.py`.

The files `Data/hires_binmask.npy` and `Data/lores_binmask.npy` are projections of the DES Y1 mask `y1a1_spt_mcal_0.2_1.3_mask.fits`, found at <http://desdr-server.ncsa.illinois.edu/despublic/y1a1_files/mass_maps/>, at the arbitrary resolutions chosen in that script.

### Section 3.1
The script `Section3_1.py` contains the code to compute the Wasserstein distances of Section 3.1's experiment. It relies on the outputs of `Section3_1_ComputeTransportPlans.py`, which computes the full transport plans for every random map realization (and takes several hours to run). These are fairly large files, and are therefore not included in this repo.

Both of these scripts import `OT_distance.py`, which contains some ad hoc optimal transport-related, pure-python functions. If you are only interested in computing optimal transport quantities, beyond the scope of this paper, we strongly recommend you look into other, maintained ressources instead of these (for instance, [POT](https://github.com/rflamary/POT)).

### Section 3.2
The script `Section3_2.py` contains the code to reproduce the experiments of Section 3.2.

The curvelet-denoised map is stored in `Data/denois_nonbin_fdr_nocoarse.fits`. It can be recomputed by running:

```mr_filter -t28 -C2 -K des_nonbin.fits denois_nonbin_fdr_nocoarse```

after installing the sparse2d module of [ISAP](http://www.cosmostat.org/software/isap). `-t` selects the type of sparse representation, with `28` being curvelets; `-C2` means the choice of threshold values, for denoising, should be chosen using False Discovery Rate; `-K` means the coarse scale will not be readded to the reconstructed map.


#### Acknowledgments

This work is one of the products of [COIN Residence Program 6 (CRP#6)](https://cosmostatistics-initiative.org/residence-programs/crp6/), held in Chamonix, France, 2019 and financially supported by the Centre national de la recherche scientifique ([CNRS](http://www.cnrs.fr/)).

The [COsmostatistics INitiative (COIN)](https://cosmostatistics-initiative.org/) receives financial support from CNRS as part of its MOMENTUM programme over the 2018-2020 period, under the project Active Learning for Large Scale Sky Surveys.

This work would not be possible without intensive consultation to online platforms and discussion forums. Although it is not possible to provide a complete list of the open source material consulted in the construction of this material, **we recognize their importance and deeply thank all those who contributes to open learning platforms**.
