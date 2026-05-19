# Second-difference temporal smoothing penalty

Builds a T x T penalty matrix `I + lambda_t * D2' D2` where `D2` is the
second-difference operator. Promotes smoothness in time while keeping
the metric positive definite.

## Usage

``` r
make_temporal_penalty(Tlen, lambda_t = 0.3)
```

## Arguments

- Tlen:

  Integer number of time points.

- lambda_t:

  Nonnegative weight on the roughness penalty (default 0.3).

## Value

Symmetric positive definite `Matrix` of size `Tlen x Tlen`.

## See also

Other temporal metrics:
[`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md),
[`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md),
[`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md),
[`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md),
[`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)

## Examples

``` r
# \donttest{
H <- make_temporal_penalty(30, lambda_t = 0.5)
# }
```
