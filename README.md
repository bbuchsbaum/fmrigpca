# genpca.fmri

Design-free, batch-friendly dimensionality reduction for fMRI using
[`genpca`](https://github.com/bbuchsbaum/genpca) and `neuroim2`'s `NeuroVec`/`NeuroVol`.

- Voxel-wise API (`NeuroVec`)
- Parcel-wise API (`ClusteredNeuroVec`)
- Rank selection & stability helpers
- Vignettes: overview & rank/stability

Install locally:

```r
install.packages("devtools")
devtools::install_local("genpca.fmri_0.1.3.tar.gz", build_vignettes = TRUE, upgrade = "never")
```
