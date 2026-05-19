# Aggregate GM/WM/CSF tissue probabilities to parcels

Computes mean tissue probabilities within each parcel label to obtain
parcel-level summaries suitable for parcel metrics.

## Usage

``` r
aggregate_tissue_to_parcels(gm_vol, wm_vol, csf_vol, parc_vol, labels = NULL)
```

## Arguments

- gm_vol, wm_vol, csf_vol:

  `NeuroVol` tissue maps.

- parc_vol:

  `NeuroVol` of parcel labels.

- labels:

  Optional integer vector of labels to include/order.

## Value

A list with numeric vectors `gm`, `wm`, `csf`, and `labels`.

## See also

Other parcel utilities:
[`parcel_centroids_from_labels()`](https://bbuchsbaum.github.io/fmrigpca/reference/parcel_centroids_from_labels.md)

## Examples

``` r
# \donttest{
# Example requires tissue maps and parcellation
# See vignette for complete example
# }
```
