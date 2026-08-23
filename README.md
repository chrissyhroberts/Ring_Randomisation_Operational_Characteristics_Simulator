# Ring-vaccination trial design model

This repository contains an individual-level simulation model for comparing ring-vaccination trial designs when the secondary attack rate (SAR), intracluster correlation coefficient (ICC), vaccine effectiveness (VE), and contact-of-contact (CoC) risk are uncertain.

The model compares operational choices such as:

- number of rings;
- immediate:delayed allocation ratio;
- timing of delayed vaccination; and
- whether contacts-of-contacts are included.

It reports both inferential performance (precision, bias, valid estimates and Pareto status) and public-health outcomes (cases prevented, vaccine doses and cases prevented per 1,000 doses).

## Repository contents

```text
ring-vax-modeller/
├── README.md
├── USER_GUIDE.md
├── data/
│   ├── example_ring_sizes.csv
│   └── DATA_DICTIONARY.md
├── docs/
│   └── METHODS_IN_PLAIN_LANGUAGE.md
├── example_outputs/
│   └── README.md
└── scripts/
    ├── install_dependencies.R
    ├── ring_vax_modeller_v2.5.2_comprehensive.R
    ├── ring_vax_analysis_v2.6.1_tradeoffs.R
    └── ring_vax_analysis_v2.6.2_visual_decisions.R
```

## Quick start

Run these commands from the repository root.

```bash
Rscript scripts/install_dependencies.R
Rscript scripts/ring_vax_modeller_v2.5.2_comprehensive.R
```

The public copy is deliberately configured for a small debug run using the synthetic dataset. It is intended to verify installation, file paths, parallel execution and output creation before a large run is started.

Results will be written below:

```text
runs/ring_vax_icc_example_debug/
```

## Running the full comprehensive grid

Open `scripts/ring_vax_modeller_v2.5.2_comprehensive.R` and edit only the clearly marked `USER CONFIGURATION` block.

At minimum:

```r
DEBUG <- FALSE
N_WORKERS <- 10L              # adjust for the machine
EMPIRICAL_RING_FILE <- "data/your_ring_sizes.csv"
RUN_ID <- "comprehensive001"  # use a new identifier for a new grid
```

The supplied comprehensive scenario ranges remain visible at the top of the script. A full run can be large; check the calculated design and chunk counts before committing substantial compute time.

## Analysing a completed run

The analysis scripts accept the input results file and output directory as command-line arguments:

```bash
Rscript scripts/ring_vax_analysis_v2.6.1_tradeoffs.R \
  runs/ring_vax_icc_comprehensive001/outputs/ring_vax_icc_operating_characteristics_comprehensive001.csv \
  runs/ring_vax_icc_comprehensive001/analysis_tradeoffs

Rscript scripts/ring_vax_analysis_v2.6.2_visual_decisions.R \
  runs/ring_vax_icc_comprehensive001/outputs/ring_vax_icc_operating_characteristics_comprehensive001.csv \
  runs/ring_vax_icc_comprehensive001/analysis_visual_decisions
```

## Grid adequacy in one sentence

Grid adequacy is the percentage of equally weighted simulated VE–SAR–ICC settings in which a design returns a sufficiently precise, sufficiently unbiased and usually valid VE estimate. It is a sensitivity score for comparing designs, not the estimated probability that a design will succeed in the real world.

See [USER_GUIDE.md](USER_GUIDE.md) for complete instructions and [docs/METHODS_IN_PLAIN_LANGUAGE.md](docs/METHODS_IN_PLAIN_LANGUAGE.md) for a non-technical explanation of the model.

## Important cautions

- The example dataset is synthetic and contains no participant data.
- The grid points are sensitivity cases, not a probability distribution over future outbreaks.
- Model results depend on the specified assumptions and should not be treated as a substitute for statistical, epidemiological or operational review.
- Validate the code and assumptions for the intended setting before using results to determine a trial protocol.

## Licence and citation

No software licence has been assigned in this package. Add an appropriate licence and project citation before publishing the repository.

