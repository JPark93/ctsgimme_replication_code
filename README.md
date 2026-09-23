# Continuous-time subgroup GIMME: simulation study

R code for the simulation study accompanying *Classification of Heterogeneous Time Series in Continuous Time*. The scripts generate data, fit continuous-time subgroup GIMME and discrete-time S-GIMME, and produce the paper's four simulation figures.

Archived generating matrices and saved paper results are supplied; fresh fits may vary because optimization can use random restarts.

## Files

| File | Purpose |
| --- | --- |
| `01_generate_data.R` | Generate data using the supplied drift matrices and subject seeds, or generate new matrices. |
| `02_fit_models.R` | Fit both methods and export subgroup memberships and simulation metrics. |
| `03_figures.R` | Produce the four simulation figures from saved summaries or fitted memberships. |
| `R/metrics.R` | Shared metric calculations and S-GIMME empty-subgroup handling. |
| `data/truth/*.rds` | Generating matrix sets for all 500 replications of each effect condition. |
| `data/paper/` | Paper summary CSVs and compressed subject-level memberships. |

Run commands from the repository root. Each main script begins with instructions and adjustable controls.

## Reproduce the paper figures

Install the plotting packages in R:

```r
install.packages(c("ggplot2", "dplyr", "scales"))
```

Then run:

```sh
Rscript 03_figures.R
```

This uses the included results and writes four PDFs to `figures/`:

| Output | Content |
| --- | --- |
| `figure_1_ari_means.pdf` | Mean adjusted Rand index (ARI) and 95% Monte Carlo confidence intervals. |
| `figure_5_ari_above_070.pdf` | Proportion with ARI greater than .70 and 95% Wilson intervals. |
| `figure_4_ari_difference_heatmap.pdf` | Mean paired ARI difference: continuous-time minus discrete-time. |
| `figure_7_cross_rate_stability.pdf` | Correct and changing subject-pair classifications across all three sampling intervals. |

Add `--png` for 320-dpi PNG copies or `--out-dir=PATH` to change the destination. PDF output requires R with Cairo support.

To calculate the summaries from the included memberships, install `mclust` and run:

```sh
Rscript 03_figures.R --memberships=data/paper/memberships.csv.gz --out-dir=figures/recomputed
```

This also writes summary CSVs to the output directory. The compressed memberships are read directly.

## Generate data and fit models

Install the fitting dependencies in R. Continuous-time fitting uses the current CRAN `ctgimme`; discrete-time fitting uses `gimme` 0.7.18:

```r
install.packages(c("ctgimme", "OpenMx", "expm", "mclust", "remotes"),
                 repos = "https://cloud.r-project.org")
remotes::install_version("gimme", version = "0.7.18", upgrade = "never")
```

Run the full workflow:

```sh
Rscript 01_generate_data.R --workers=4
Rscript 02_fit_models.R --cores=4
Rscript 03_figures.R --memberships=results/memberships.csv --out-dir=figures/refitted
```

Both generation and fitting accept `--effects=Weak,Mid,Large`, `--replications=1:500`, and `--output-dir=PATH`. Fitting also accepts `--data-dir=PATH` and `--deltas=0.5,1,5`. `--workers` controls parallel generation across replications; `--cores` controls continuous-time fitting workers. Model conditions run sequentially.

Full generation saves complete trajectories and requires about **540 GB** of disk space and **1 GB RAM per worker**. Add `--retain-only=true` to the generation command to save only the 460 rows per subject needed for the study, reducing storage to about **2.5 GB** while preserving all fitted observations. This option still generates complete trajectories before retaining rows, so it does not reduce generation time. Complete trajectories are needed to extend the retained series length or sampling schedule.

Generation skips existing dataset/truth pairs unless `--overwrite=true` is supplied. Use a new output directory when changing settings. Fitting resumes completed jobs when inputs, settings and software versions match; failures provide a log location. The combined `results/memberships.csv` contains completed jobs in that output directory. Data directories contain `Simulations_<effect>/Dataset <rep>.RDS` and matching `Trues <rep>.RDS` files.

Use the supplied matrices to replicate the study design. `--new-truth=true --truth-seed=12345` generates new matrices. Figures from newly fitted memberships require all three effects and sampling intervals, 60 subjects per fit, and at least two matched replications per effect; the paper uses 500.

## Study settings

- Subgroup effects: **0.30, 0.60 and 0.90**; sampling intervals: **0.5, 1 and 5**.
- **500 replications per effect**, reusing trajectories at all intervals: **1,500 datasets and 4,500 condition evaluations**.
- **60 subjects**, three equally sized subgroups, six observed variables, and 200 retained observations per subject at each interval.
- Generation step .01; 100,000 observations after a 1,000-observation burn-in; identity diffusion matrix; measurement-error variance `1e-5`; process-noise covariance rounded to five decimals.
- Continuous-time fitting: group/subgroup inclusion .55; group/subgroup alpha .05; individual alpha .01; Benjamini-Hochberg correction; PAM with at most six subgroups; `ME.var = diag(1e-5, 6)` and `PE.var = diag(1, 6)`.
- Discrete-time fitting: `gimme(data, out, subgroup = TRUE)`.

Each fit's `metrics.rds` stores Type I error, power, absolute bias, relative bias, RMSE, and both methods' ARIs rounded to two decimals. Nonzero-truth coefficient metrics include diagonal drift coefficients; retained group/subgroup paths count as present even when individual tests are nonsignificant. Figure ARIs use full-precision memberships. Cross-rate summaries match subject pairs within replication across all three sampling intervals.
