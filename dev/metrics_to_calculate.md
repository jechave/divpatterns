# Metrics to calculate

## Residue-dependent

| Metric | Function | Measures |
|---|---|---|
| `dr2i` | `delta_structure_dr2i` | squared displacement |
| `dh_site` | `delta_motion_dh_site` | change in fluctuation entropy |
| `dlcov` | `delta_motion_dlcov_site` | log-covariance distance |
| `size` | `delta_motion_dlcov_site` | size part of `dlcov` |
| `shape` | `delta_motion_dlcov_site` | shape part of `dlcov` |
| `orientation` | `delta_motion_dlcov_site` | orientation part of `dlcov` |

## Mode-dependent

| Metric | Function | Measures |
|---|---|---|
| `dr2n` | `delta_structure_dr2n` | squared displacement along each mode |
| `dh_mode` | `delta_motion_dh_mode` | change in fluctuation entropy along each mode |
| `nmix_mode` | `delta_motion_nmix_mode` | number of modes each mode spreads over |

## From `msamodel`, on both axes

| Metric | Definition |
|---|---|
| `lrmsd` | `log(sqrt(dr2))` |
| `nlrmsd` | `lrmsd` minus its mean |
