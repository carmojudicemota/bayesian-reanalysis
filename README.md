# Bayesian Reanalysis of Computationally Reproducible Studies in Psychology

This repository contains a reproducible Bayesian reanalysis of statistically significant results from studies in psychology of teaching and education.

The project separates:

1. extraction and registration of published claims;
2. frequentist reproduction or recomputation;
3. Bayesian reanalysis;
4. frequentist–Bayesian concordance classification;
5. tables, figures, and reporting outputs.

The repository is organised so that the analyses currently marked as ready can be reproduced from the raw study data, the locked R environment, and a single run-all script.

---

## Quick start

Clone the repository:

```bash
git clone https://github.com/carmojudicemota/bayesian-reanalysis-v1.git
cd bayesian-reanalysis-v1
```

Open the project in RStudio or start R from the repository root.

Install `renv` if necessary:

```r
install.packages("renv")
```

Restore the project environment:

```r
renv::restore()
```

Run the complete pipeline:

```r
source("scripts/00_run_all.R")
```

This is the recommended and canonical way to reproduce the project.

The run-all script executes the required stages in order, including:

1. rebuilding the derived claim registry;
2. running the implemented frequentist reproductions;
3. running the Bayesian analyses for claims marked `ready`;
4. combining Wave 1, Wave 2, and Wave 3 results;
5. rebuilding dataset and concordance outputs;
6. generating the project figures.

Claims marked `pending` or `held` are not treated as completed Bayesian analyses.

---

## Repository structure

```text
.
├── config/
│   └── claim_map.csv
├── data/
│   ├── raw/
│   ├── source/
│   └── derived/
├── R/
│   ├── prepare/
│   ├── reproduce/
│   ├── analyse/
│   │   ├── wave2/
│   │   └── wave3/
│   └── visualise/
├── outputs/
│   ├── reproduced/
│   ├── tables/
│   └── figures/
├── scripts/
│   ├── 00_run_all.R
│   ├── 01_run_analysis.R
│   └── 02_run_figures.R
├── renv/
├── renv.lock
└── README.md
```

### Main directories

* `config/`: project configuration, including claim readiness and analysis mappings.
* `data/raw/`: raw study datasets used in the reproductions and Bayesian analyses.
* `data/source/`: source tables and manually curated project inputs.
* `data/derived/`: generated analysis-ready data, including the claim registry.
* `R/prepare/`: functions used to prepare and refresh project inputs.
* `R/reproduce/`: study-specific frequentist reproductions and recomputations.
* `R/analyse/`: Bayesian analyses and combined analysis outputs.
* `R/analyse/wave2/`: Wave 2 Bayesian analyses.
* `R/analyse/wave3/`: Wave 3 Bayesian analyses.
* `R/visualise/`: table and figure generation.
* `outputs/reproduced/`: reproduced frequentist results.
* `outputs/tables/`: Bayesian results, dataset summaries, and concordance tables.
* `outputs/figures/`: generated figures.
* `scripts/`: top-level project entry points.

---

## Requirements

The project uses R and `renv` to record package dependencies.

The analyses were developed using R 4.5 on macOS. Other operating systems should work, although packages containing compiled code may need to be rebuilt locally.

Git is required to clone the repository.

For the analyses currently marked as ready, no external software is required beyond R and the packages recorded in `renv.lock`.

Some future Wave 3 analyses will require CmdStan. CmdStan is not installed by `renv` and is not currently required to reproduce the claims marked as ready.

---

## Clone the repository

```bash
git clone https://github.com/carmojudicemota/bayesian-reanalysis-v1.git
cd bayesian-reanalysis-v1
```

Open the project in RStudio or start R from the repository root.

Confirm that the working directory is the repository root:

```r
getwd()
```

All project paths are relative to the repository root.

---

## Restore the R environment

Install `renv` if necessary:

```r
install.packages("renv")
```

Restore the package versions recorded in `renv.lock`:

```r
renv::restore()
```

Confirm that the project library matches the lockfile:

```r
renv::status()
```

A clean environment should report that the project is consistent.

Packages containing compiled code may need to be rebuilt when moving between R versions, computer architectures, or operating systems:

```r
renv::rebuild()
```

To rebuild a specific package:

```r
renv::rebuild("package_name")
```

For example:

```r
renv::rebuild(
  packages = c(
    "RcppParallel",
    "blavaan"
  ),
  recursive = TRUE,
  type = "source",
  prompt = FALSE
)
```

This may be necessary if a compiled package was restored from an incompatible local cache.

---

## Raw data

The raw datasets required by the implemented reproductions and Bayesian analyses are stored under:

```text
data/raw/
```

Each study has its own directory:

```text
data/raw/study_XX/
```

The analysis scripts refer to these files using paths relative to the repository root.

Important: The raw files are third-party research datasets. Their inclusion in this repository is intended to support research reproducibility. Ownership remains with the original authors or data providers. Users should consult the relevant articles, data repositories, licences, and terms of use before redistributing or reusing individual datasets outside this project. The presence of a dataset in this repository should not be interpreted as a transfer of ownership or as permission for unrestricted redistribution.

---

## Run the complete pipeline

After restoring the environment, run:

```r
source("scripts/00_run_all.R")
```

This is the primary project entry point.

The run-all script should execute the project stages in their required order:

1. refresh the derived claim registry;
2. run the implemented frequentist reproductions or recomputations;
3. run the implemented Bayesian analysis waves;
4. combine Wave 1, Wave 2, and Wave 3 Bayes-factor results;
5. build the dataset summaries;
6. classify frequentist–Bayesian concordance;
7. generate the project figures.

A successful run should regenerate all derived outputs without requiring manual edits to generated files.

---

## Claim registry

The main analysis registry is generated from the project source files.

The generated registry is:

```text
data/derived/claims.csv
```

Claim readiness is controlled by:

```text
config/claim_map.csv
```

Only claims with:

```text
status = ready
```

are included in the corresponding Bayesian analysis pipeline.

Claims marked:

```text
status = pending
```

require further implementation, validation, or methodological work.

Claims marked:

```text
status = held
```

are intentionally excluded from the completed analysis, usually because they are out of scope, not reproducible, or otherwise not eligible under the project rules.

The run-all script refreshes the derived registry before running the analysis.

To refresh it independently during development:

```r
source("R/prepare/00_refresh_claims.R")
```

---

## Analysis outputs

The principal combined Bayes-factor output is:

```text
outputs/tables/bayes_factor_results.csv
```

Wave-specific outputs include:

```text
outputs/tables/bayes_factor_results_wave1.csv
outputs/tables/bayes_factor_results_wave2.csv
outputs/tables/bayes_factor_results_wave3.csv
```

Additional outputs include:

* dataset summaries;
* claim-level concordance classifications;
* study-level summaries;
* figures derived from the combined analysis results.

Generated tables are written under:

```text
outputs/tables/
```

Generated figures are written under:

```text
outputs/figures/
```

---

## Running individual stages

The run-all script is the recommended entry point. Individual stages may still be run during development or debugging.

### Refresh the claim registry

```r
source("R/prepare/00_refresh_claims.R")
```

### Run the Bayesian analysis pipeline

```r
source("scripts/01_run_analysis.R")
```

This script:

1. runs the implemented Bayesian analysis waves;
2. combines the Wave 1, Wave 2, and Wave 3 Bayes-factor results;
3. builds the dataset summaries;
4. classifies frequentist–Bayesian concordance;
5. writes the analysis tables.

### Generate figures

```r
source("scripts/02_run_figures.R")
```

The figure pipeline expects every included claim to contain exactly one primary Bayes-factor row identified by:

```text
prior_label = primary
```

Sensitivity-analysis rows use separate prior labels.

Running an individual stage assumes that all required earlier stages and inputs are already up to date. For complete reproduction, use:

```r
source("scripts/00_run_all.R")
```

---

## Reproduce individual studies

Study-specific frequentist reproductions are stored under:

```text
R/reproduce/
```

For example:

```r
source("R/reproduce/study_44.R")

reproduce_study_44()
```

The resulting files are written under:

```text
outputs/reproduced/
```

The reproduction layer verifies an original or reconstructed frequentist result.

The Bayesian analysis layer is separate and is stored under:

```text
R/analyse/
```

A Bayesian analysis may reuse the same raw-data preparation and model structure, but it does not replace or modify the frequentist reproduction.

The separation between reproduction and Bayesian reanalysis is intentional:

* `R/reproduce/` establishes the analysis input and verifies the frequentist result;
* `R/analyse/` performs the Bayesian reanalysis;
* `config/claim_map.csv` controls which verified claims are ready for inclusion.

---

## Bayesian analysis waves

### Wave 1

Wave 1 covers statistical families for which established default Bayes-factor procedures can be applied directly.

### Wave 2

Wave 2 extends the analysis to designs including:

* Welch tests;
* factorial models;
* regression models;
* repeated-measures designs;
* mixed designs.

The Wave 2 implementation is stored under:

```text
R/analyse/wave2/
```

### Wave 3

Wave 3 covers specialised statistical families, including:

* multivariate linear models;
* MANOVA;
* rank correlations;
* structural equation models;
* latent growth models;
* other analyses requiring family-specific Bayesian implementations.

The Wave 3 implementation is stored under:

```text
R/analyse/wave3/
```

Currently integrated Wave 3 studies include:

* Study 03;
* Study 18;
* Study 44.

These studies use `BFpack` generalized fractional Bayes factors for joint coefficient restrictions in multivariate linear models.

The following specialised studies remain under development:

* Study 20: Bayesian latent-difference-score SEM;
* Study 26: latent-rank correlation method;
* Study 39: latent-rank Spearman Bayes factor;
* Study 40: Bayesian latent-basis growth-model comparison;
* Study 45: pending implementation.

Pending studies remain marked as `pending` or `held` in `config/claim_map.csv` and are not included as completed Bayesian results.



---

## CmdStan

The `cmdstanr` R package is recorded in `renv.lock`, but CmdStan itself is an external dependency.

CmdStan is not currently required to reproduce the studies marked as ready. It is expected to be potentially required for the custom latent-rank implementation used by Study 39 and maybe Study 26.

When those analyses are implemented, check the local C++ toolchain:

```r
cmdstanr::check_cmdstan_toolchain()
```

Install CmdStan:

```r
cmdstanr::install_cmdstan()
```

Set or detect the installation path:

```r
cmdstanr::set_cmdstan_path()
```

Verify the installation:

```r
cmdstanr::cmdstan_path()
cmdstanr::cmdstan_version()
```

The final project release will record the exact CmdStan version used for the latent-rank analyses. Until those analyses are marked as ready, the absence of CmdStan does not prevent reproduction of the currently completed pipeline.



---

## Updating the environment

After intentionally installing or updating a package used by the project, update the lockfile:

```r
renv::snapshot()
```

Check the resulting environment status:

```r
renv::status()
```

Commit environment changes together with the code that requires them:

```bash
git add renv.lock
git commit -m "Update project dependencies"
git push
```

Do not update packages or the lockfile merely to repair an unrelated local installation unless the project intentionally adopts the new package versions.

The `renv.lock` file records R package versions. It does not record external software installations such as CmdStan itself.

---



## Project status

The repository is reproducible for the analyses currently marked as `ready`.

The canonical reproduction command is:

```r
source("scripts/00_run_all.R")
```

Full Wave 3 reproducibility is still in progress because several specialised analyses require additional validated Bayesian implementations.

Incomplete analyses are explicitly excluded from the completed pipeline through their `pending` or `held` status.
