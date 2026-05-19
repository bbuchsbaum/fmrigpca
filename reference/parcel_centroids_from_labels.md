# Parcel centroids from a label volume

Computes centroid coordinates (i,j,k) for labeled parcels in a
`NeuroVol` label image. Labels with value 0 are ignored.

## Usage

``` r
parcel_centroids_from_labels(parc_vol, labels = NULL)
```

## Arguments

- parc_vol:

  `NeuroVol` of parcel labels (integer; 0 treated as background).

- labels:

  Optional integer vector of labels to include.

## Value

A list with `coords` (P x 3 matrix) and `labels` (length P).

## Details

For each parcel label, voxel coordinates are averaged to obtain a
centroid. Optionally restrict to a provided vector of `labels`.

## See also

Other parcel utilities:
[`aggregate_tissue_to_parcels()`](https://bbuchsbaum.github.io/fmrigpca/reference/aggregate_tissue_to_parcels.md)

## Examples

``` r
# \donttest{
pc <- parcel_centroids_from_labels(parc_vol)
#> Error: object 'parc_vol' not found
# }
```
