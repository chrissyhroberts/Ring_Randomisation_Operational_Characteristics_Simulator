# User guide

## 1. What the model does

The model simulates ring-vaccination trials under many possible epidemiological conditions. For each condition, it repeatedly generates rings, assigns rings to immediate or delayed vaccination, applies the assumed vaccination delay, and estimates the operational VE contrast.

The main uncertainty dimensions are:

- **True VE:** the underlying biological effect of vaccination.
- **SAR:** the average secondary attack rate among direct contacts.
- **ICC:** the degree to which infection outcomes cluster within the same ring.
- **CoC risk multiplier:** the risk among contacts-of-contacts relative to direct contacts.

The operational dimensions are:

- number of rings;
- immediate:delayed allocation ratio; and
- delayed-arm vaccination day.

## 2. Software requirements

- R 4.2 or later is recommended.
- Required packages: `future`, `future.apply`, `dplyr`, `ggplot2`, and `tidyr`.

Install the packages with:

```bash
Rscript scripts/install_dependencies.R
```

## 3. Input data

The model needs a CSV containing one row per observed ring and these exact columns:

| Column | Meaning |
|---|---|
| `r_case_index_id` | Unique ring or index-case identifier |
| `direct_contact` | Number of eligible direct contacts in that ring |
| `contact_of_contact` | Number of eligible contacts-of-contacts in that ring |

The empirical rows are sampled to reproduce realistic variation in ring size. The example file in `data/example_ring_sizes.csv` is entirely synthetic.

Input rules:

- `direct_contact` must be a positive integer;
- `contact_of_contact` must be zero or a positive integer;
- missing CoC counts are converted to zero;
- rows with missing or non-positive direct-contact counts are discarded.

## 4. First run: debug mode

Run commands from the repository root so the relative input path resolves correctly:

```bash
Rscript scripts/ring_vax_modeller_v2.5.2_comprehensive.R
```

The packaged configuration uses:

```r
DEBUG <- TRUE
PARALLEL <- FALSE
N_WORKERS <- 2L
DEBUG_N_SIM <- 50L
DEBUG_N_DESIGNS <- 18L
RUN_ANALYSIS <- FALSE
```

This is a software check, not a scientific analysis. With few replicates and an incomplete design grid, estimates will be noisy and comparative grid summaries are not valid.

The public debug configuration uses sequential execution so it also works in restricted containers and continuous-integration environments. Set `PARALLEL <- TRUE` for a substantive local run after the debug check succeeds.

Confirm that:

1. all three stages complete;
2. the output CSV is created;
3. each selected design has the expected number of valid replicates; and
4. no path or dependency errors appear.

## 5. Configuring a real run

Edit only the `USER CONFIGURATION` block at the beginning of the model script.

### Run controls

```r
DEBUG <- FALSE
PARALLEL <- TRUE
N_WORKERS <- 10L
RESUME <- TRUE
SEED <- 20260821L
```

- Set `N_WORKERS` below the number of logical CPU cores if the computer must remain responsive.
- `RESUME <- TRUE` reuses completed cache and design chunks for the same run.
- The seed is used to make the parallel simulation reproducible.

### Simulation size

```r
N_SIM <- 2000L
REPLICATE_CHUNK_SIZE <- 500L
```

Increasing `N_SIM` reduces Monte Carlo noise but increases runtime and disk use. Smaller chunks create more files and scheduling overhead; larger chunks take longer to redo if interrupted.

### Run identity

Use a fixed `RUN_ID` to resume the same run. Use a new `RUN_ID` after changing the design grid or scientific assumptions.

```r
RUN_ID <- "comprehensive001"
```

The run manifest prevents an existing run identifier from silently being reused with a different grid.

### Scenario and design grids

The packaged comprehensive grid is:

```r
TRUE_VE_VALUES <- seq(0.1, 0.9, by = 0.1)
SAR_VALUES <- c(0.01, 0.02, 0.04, 0.06, 0.10)
ICC_VALUES <- c(0.05, 0.10, 0.15, 0.20, 0.30)
RING_NUMBERS <- c(500L, 750L, 900L, 1000L, 1100L, 1200L, 1500L)
ALLOCATION_RATIOS <- c(1, 2, 3)
DELAYS <- c(7L, 14L, 21L)
COC_RISK_MULTIPLIERS <- c(0.05, 0.10, 0.25, 0.50, 1.00)
```

Each additional value multiplies the number of designs. Before a large run, calculate:

```text
number of designs = VE × SAR × ICC × CoC × rings × allocations × delays
```

## 6. What the three stages mean

### Stage 1: epidemic cache

The script generates reusable epidemic scenarios for each ring-number, SAR and ICC combination. Caching avoids regenerating the same underlying outbreak for every allocation and delay choice.

### Stage 2: design simulation

Each trial design is applied to the cached epidemics. This is usually the longest stage because it contains the largest number of independent design-by-replicate chunks.

### Stage 3: aggregation and testing

The chunks for each design are combined, operating characteristics are calculated, and Pareto comparisons are made.

All three stages use a flat worker pool when parallel processing is enabled. The model does not start nested worker pools.

## 7. Monitoring progress

From another terminal, count completed design chunks with:

```bash
find runs/ring_vax_icc_RUN_ID/design_chunks -name '*.rds' | wc -l
```

Count completed design checkpoints during or after Stage 3 with:

```bash
find runs/ring_vax_icc_RUN_ID/checkpoints -name '*.rds' | wc -l
```

Replace `RUN_ID` with the actual identifier. The stage-start messages report the expected total.

## 8. Output structure

```text
runs/ring_vax_icc_RUN_ID/
├── cache/          reusable epidemic chunks
├── checkpoints/    aggregated design checkpoints
├── design_chunks/  simulated design-replicate chunks
├── logs/           stage timing summaries
├── outputs/        final operating-characteristics and Pareto tables
├── analysis/       optional built-in analysis
├── design_grid.csv
├── model_script.R
└── run_manifest.rds
```

Keep `run_manifest.rds`, the copied model script and the final output together for reproducibility.

## 9. Separate analysis scripts

### Trade-off analysis v2.6.1

This analysis asks how much robustness is obtained at each ring budget and quantifies diminishing returns. It creates:

- a design league table;
- recommendations under ring caps;
- recommendations for adequacy and robustness-retention targets;
- marginal robustness per additional 100 rings;
- prevention choices under explicit adequacy floors; and
- resource–robustness figures.

Run it with:

```bash
Rscript scripts/ring_vax_analysis_v2.6.1_tradeoffs.R INPUT.csv OUTPUT_DIRECTORY
```

### Visual decision analysis v2.6.2

This analysis compares the shortlisted designs side by side, including:

- complete-grid and focus-region SAR–ICC heatmaps;
- CoC-risk sensitivity;
- a summary of shortlisted designs;
- an all-design trade-off plot; and
- an operational design atlas.

Run it with:

```bash
Rscript scripts/ring_vax_analysis_v2.6.2_visual_decisions.R INPUT.csv OUTPUT_DIRECTORY
```

Edit `SHORTLIST_DESIGNS` near the top of the script if the operational shortlist changes.

## 10. Understanding grid adequacy

Think of every combination of true VE, SAR and ICC as a different possible world.

A design passes in one world only when all of these are true:

1. SD of the estimated VE is no more than 0.15;
2. absolute VE bias is no more than 0.05; and
3. at least 95% of simulation replicates return a valid VE estimate.

If 146 of 225 worlds pass, grid adequacy is:

```text
146 / 225 = 65%
```

This does **not** mean the design has a 65% probability of success. The grid points are deliberately given equal weight so that designs can be compared over a transparent sensitivity range. A probability statement would require a justified probability distribution over the unknown true VE, SAR and ICC.

The declared focus region repeats the same calculation over a narrower, explicitly labelled subset. It is a sensitivity analysis, not a hidden prior distribution.

## 11. CoC interpretation

The CoC multiplier is the baseline infection risk among contacts-of-contacts relative to direct contacts. For example, if direct-contact SAR is 2%:

- multiplier 0.05 implies CoC risk 0.1%;
- multiplier 0.10 implies CoC risk 0.2%; and
- multiplier 0.25 implies CoC risk 0.5%.

CoC vaccination can increase cases prevented, but it also increases vaccine use. It does not create additional direct-contact endpoint information for the primary direct-contact VE contrast.

## 12. Troubleshooting

### The input file cannot be found

Run from the repository root or change `EMPIRICAL_RING_FILE` to an absolute path.

### A resumed run reports a different grid

The existing run identifier belongs to another configuration. Restore the original configuration or select a new `RUN_ID`.

### Parallel execution fails

Confirm that `future` and `future.apply` are installed. Reduce `N_WORKERS` if memory pressure is high.

### Analysis figures are empty or incomplete

Confirm that the input is the final operating-characteristics CSV from a complete factorial production run. A debug subset is not suitable for grid-level comparative analysis.

### Cloud-synchronised folders are slow

Large runs create many small checkpoint files. A local, non-synchronised working directory may reduce filesystem overhead. Copy final outputs into the project repository after completion.

## 13. Reproducibility checklist

- Record the script version and Git commit.
- Retain the model script copied into the run directory.
- Retain `run_manifest.rds` and `design_grid.csv`.
- Record R and package versions using `sessionInfo()`.
- Use a fixed random seed.
- Give every changed grid a new run identifier.
- Report complete-grid and focus-region results separately.
- State explicitly that equally weighted grid cases are not probabilities.
