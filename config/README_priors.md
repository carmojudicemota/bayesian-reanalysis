# Prior registry

This directory records every prior used in the project so that the analytic
choices are inspectable in one place. The Wave 1 and core Wave 2 priors are
preregistered; the Wave 3 families were prospectively specified as planned
extensions (fixed before the final Bayesian result was produced) rather than
part of the original preregistration. Priors are split by wave, but all Wave 2
and Wave 3 families share the same tidy schema and can be row-bound into a
single table.

## Files

| File | Wave | Schema |
|------|------|--------|
| `priors_wave1.csv` | Wave 1 | `family, prior_label, rscale` |
| `priors_wave2.csv` | Wave 2 | `prior_family, prior_label, param, value` |
| `priors_wave3.csv` | Wave 3 | `prior_family, prior_label, param, value` |

