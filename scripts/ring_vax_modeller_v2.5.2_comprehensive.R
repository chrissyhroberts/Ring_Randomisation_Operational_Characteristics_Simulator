# ============================================================
# Ring vaccination operating-characteristics model v2.5.2 (comprehensive grid)
# ============================================================
#
# Primary estimand:
#   VE_immediate_vs_delayed (an operational trial estimand, not biological VE)
#
# Scientific redesign:
#   * ICC is the primary clustering assumption.
#   * Ring-level baseline risks follow a beta distribution parameterised by
#     mean SAR and ICC. Conditional Bernoulli outcomes therefore have the
#     requested exchangeable within-ring correlation (beta-binomial model).
#   * Latent SAR variation and observed ring-level SAR CV are outputs.
#
# Execution redesign:
#   * Compact, scenario-specific caches.
#   * Three flat parallel stages: cache chunks, design-simulation chunks,
#     then design aggregation/testing. No nested worker oversubscription.
#   * Run-specific folders, checkpoints, deterministic seeds and timing logs.
# ============================================================


# ============================================================
# USER CONFIGURATION -- EDIT THIS BLOCK ONLY
# ============================================================

# Run mode
DEBUG <- TRUE
PARALLEL <- FALSE
N_WORKERS <- 2L
RESUME <- TRUE
SEED <- 20260821L

# Simulation size
N_SIM <- 2000L
REPLICATE_CHUNK_SIZE <- 500L

# Smaller settings used when DEBUG = TRUE
DEBUG_N_SIM <- 50L
DEBUG_N_DESIGNS <- 18L

# Input and output
EMPIRICAL_RING_FILE <- "data/example_ring_sizes.csv"
OUTPUT_ROOT <- "runs"

# Use a fixed name, such as "production_001", to resume the same run.
# The timestamp below creates a new run every time the script starts.
RUN_ID <- "example_debug"

# Scientific scenario grid
TRUE_VE_VALUES <- seq(0.1, 0.9, by = 0.1)
SAR_VALUES <- c(0.01, 0.02, 0.04, 0.06, 0.10)
ICC_VALUES <- c(0.05, 0.10, 0.15, 0.20, 0.30)
RING_NUMBERS <- c(
  500L,
  750L,
  900L,
  1000L,
  1100L,
  1200L,
  1500L
)

# Trial-design grid
ALLOCATION_RATIOS <- c(1, 2, 3)
DELAYS <- c(7L,14L,21L)

# Contact-of-contact assumptions
INCLUDE_COC <- TRUE
# Relative baseline risk versus a direct contact. Values span the plausible
# central range and retain 1.00 as an equal-risk upper-bound sensitivity.
COC_RISK_MULTIPLIERS <- c(0.05, 0.10, 0.25, 0.50, 1.00)

# Automated analysis
RUN_ANALYSIS <- FALSE
# The comprehensive grid creates hundreds of figure panels. Keep table outputs
# on, but skip unreadable overview figures; use a filtered run for final plots.
RUN_FIGURES <- FALSE
DECISION_PRECISION_THRESHOLD <- 0.15
PRECISION_THRESHOLDS <- c(0.10, 0.15, 0.20)
FIGURE_DPI <- 300L

# ============================================================
# END USER CONFIGURATION
# ============================================================


# Convert the visible settings above into the internal configuration object.
default_config <- function() {
  list(
    debug_mode = DEBUG,
    parallel = PARALLEL,
    n_workers = N_WORKERS,
    resume = RESUME,
    seed = SEED,

    n_sim = N_SIM,
    replicate_chunk_size = REPLICATE_CHUNK_SIZE,
    debug_n_sim = DEBUG_N_SIM,
    debug_n_designs = DEBUG_N_DESIGNS,

    empirical_ring_file = EMPIRICAL_RING_FILE,
    output_root = OUTPUT_ROOT,
    run_id = RUN_ID,

    true_ve_values = TRUE_VE_VALUES,
    sar_values = SAR_VALUES,
    icc_values = ICC_VALUES,
    ring_numbers = RING_NUMBERS,
    allocation_ratios = ALLOCATION_RATIOS,
    delays = DELAYS,

    include_coc = INCLUDE_COC,
    coc_risk_multipliers = COC_RISK_MULTIPLIERS,
    run_analysis = RUN_ANALYSIS,
    run_figures = RUN_FIGURES,
    decision_precision_threshold = DECISION_PRECISION_THRESHOLD,
    precision_thresholds = PRECISION_THRESHOLDS,
    figure_dpi = FIGURE_DPI
  )
}


validate_config <- function(config) {
  stopifnot(
    config$n_workers >= 1L,
    config$n_sim >= 1L,
    config$replicate_chunk_size >= 1L,
    all(config$true_ve_values >= 0 & config$true_ve_values <= 1),
    all(config$sar_values > 0 & config$sar_values < 1),
    all(config$icc_values > 0 & config$icc_values < 1),
    all(config$ring_numbers >= 2L),
    all(config$allocation_ratios > 0),
    all(config$delays >= 0),
    all(is.finite(config$coc_risk_multipliers)),
    all(config$coc_risk_multipliers >= 0),
    length(config$run_figures) == 1L,
    is.logical(config$run_figures),
    length(config$decision_precision_threshold) == 1L,
    is.finite(config$decision_precision_threshold),
    config$decision_precision_threshold > 0,
    all(is.finite(config$precision_thresholds)),
    all(config$precision_thresholds > 0)
  )
  invisible(config)
}


# -------------------------
# INPUT AND RUN MANAGEMENT
# -------------------------

load_empirical_rings <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Empirical ring file not found: ", path, "\n",
      "Edit EMPIRICAL_RING_FILE in the USER CONFIGURATION block."
    )
  }

  raw <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("r_case_index_id", "direct_contact", "contact_of_contact")
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop("Ring file is missing columns: ", paste(missing, collapse = ", "))
  }

  rings <- data.frame(
    ring_id = raw$r_case_index_id,
    direct_contacts = as.integer(raw$direct_contact),
    coc_contacts = as.integer(raw$contact_of_contact)
  )

  rings$coc_contacts[is.na(rings$coc_contacts)] <- 0L
  rings <- rings[
    !is.na(rings$direct_contacts) &
      rings$direct_contacts > 0L &
      rings$coc_contacts >= 0L,
    ,
    drop = FALSE
  ]

  if (!nrow(rings)) stop("No usable rings remained after input validation.")
  rings
}


script_path <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) return(NA_character_)
  normalizePath(sub("^--file=", "", arg[1]), mustWork = FALSE)
}


initialise_run <- function(config, design_grid) {
  run_dir <- file.path(config$output_root, paste0("ring_vax_icc_", config$run_id))
  paths <- list(
    run_dir = run_dir,
    cache_dir = file.path(run_dir, "cache"),
    checkpoint_dir = file.path(run_dir, "checkpoints"),
    design_chunk_dir = file.path(run_dir, "design_chunks"),
    output_dir = file.path(run_dir, "outputs"),
    analysis_dir = file.path(run_dir, "analysis"),
    log_dir = file.path(run_dir, "logs")
  )
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

  model_file <- script_path()
  if (!is.na(model_file) && file.exists(model_file)) {
    file.copy(model_file, file.path(run_dir, "model_script.R"), overwrite = TRUE)
  }

  manifest_file <- file.path(run_dir, "run_manifest.rds")
  manifest <- list(config = config, design_grid = design_grid)
  if (file.exists(manifest_file) && config$resume) {
    old <- readRDS(manifest_file)
    if (!identical(old$design_grid, design_grid)) {
      stop("The requested design grid differs from the existing run manifest.")
    }
  } else {
    saveRDS(manifest, manifest_file)
    write.csv(design_grid, file.path(run_dir, "design_grid.csv"), row.names = FALSE)
  }
  paths
}


# -------------------------
# ICC / BETA-BINOMIAL MODEL
# -------------------------

# If p_j ~ Beta(alpha, beta) and Y_ij | p_j ~ Bernoulli(p_j), then
# Corr(Y_ij, Y_i'j) = 1 / (alpha + beta + 1) for i != i'.
beta_parameters_from_icc <- function(mean_sar, icc) {
  if (mean_sar <= 0 || mean_sar >= 1) stop("mean_sar must be in (0, 1).")
  if (icc <= 0 || icc >= 1) stop("icc must be in (0, 1).")
  concentration <- 1 / icc - 1
  c(
    alpha = mean_sar * concentration,
    beta = (1 - mean_sar) * concentration
  )
}


draw_ring_sar_icc <- function(n, mean_sar, icc) {
  pars <- beta_parameters_from_icc(mean_sar, icc)
  rbeta(n, pars[["alpha"]], pars[["beta"]])
}


implied_latent_sar_cv <- function(mean_sar, icc) {
  sqrt(icc * mean_sar * (1 - mean_sar)) / mean_sar
}


cv_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || mean(x) == 0) return(NA_real_)
  stats::sd(x) / mean(x)
}


# Method-of-moments diagnostic. It removes the average binomial sampling
# component from between-ring variance before scaling by Bernoulli variance.
estimate_binary_icc <- function(cases, sizes) {
  ok <- is.finite(cases) & is.finite(sizes) & sizes > 0
  cases <- cases[ok]
  sizes <- sizes[ok]
  if (length(cases) < 2L) return(NA_real_)

  p_hat <- sum(cases) / sum(sizes)
  if (p_hat <= 0 || p_hat >= 1) return(NA_real_)
  ring_rates <- cases / sizes
  observed_between <- stats::var(ring_rates)
  sampling_component <- mean(p_hat * (1 - p_hat) / sizes)
  (observed_between - sampling_component) / (p_hat * (1 - p_hat))
}


design_effect <- function(mean_cluster_size, icc, cluster_size_cv = 0) {
  # The unequal-size extension reduces to 1 + (m - 1)ICC when CV_m = 0.
  1 + (((1 + cluster_size_cv^2) * mean_cluster_size) - 1) * icc
}


# -------------------------
# BIOLOGY AND RANDOMISATION
# -------------------------

protection_curve <- function(days) {
  approx(
    x = c(0, 7, 14, 21, 28),
    y = c(0, 0.25, 0.60, 0.80, 1),
    xout = pmax(days, 0),
    rule = 2
  )$y
}


mean_protection_over_exposure_window <- function(vaccination_day) {
  mean(protection_curve(0:28 - vaccination_day))
}


make_allocation_matrix <- function(n_sim, n_rings, ratio) {
  n_immediate <- as.integer(round(n_rings * ratio / (ratio + 1)))
  ans <- matrix(FALSE, nrow = n_sim, ncol = n_rings)
  for (i in seq_len(n_sim)) {
    ans[i, sample.int(n_rings, n_immediate)] <- TRUE
  }
  ans
}


# -------------------------
# COMPACT EPIDEMIC CACHE
# -------------------------

scenario_key <- function(n_rings, n_sim, sar, icc) {
  paste0(
    "r", n_rings,
    "_s", format(sar, scientific = FALSE, trim = TRUE),
    "_i", format(icc, scientific = FALSE, trim = TRUE),
    "_n", n_sim
  )
}


scenario_seed <- function(base_seed, n_rings, sar, icc) {
  as.integer((base_seed + n_rings * 1009 + round(sar * 1e5) * 101 +
    round(icc * 1e5) * 103) %% .Machine$integer.max)
}


cache_path <- function(paths, n_rings, n_sim, sar, icc) {
  file.path(paths$cache_dir, paste0("cache_", scenario_key(n_rings, n_sim, sar, icc), ".rds"))
}


generate_epidemic_cache <- function(
  rings, n_rings, n_sim, sar, icc, seed
) {
  set.seed(seed)
  n_source_rings <- nrow(rings)
  ring_row_index <- matrix(
    sample.int(n_source_rings, n_sim * n_rings, replace = TRUE),
    nrow = n_sim,
    ncol = n_rings
  )
  latent_sar <- matrix(
    draw_ring_sar_icc(n_sim * n_rings, sar, icc),
    nrow = n_sim,
    ncol = n_rings
  )

  direct_sizes <- matrix(
    rings$direct_contacts[ring_row_index],
    nrow = n_sim,
    ncol = n_rings
  )
  baseline_cases <- matrix(
    rbinom(length(latent_sar), size = as.vector(direct_sizes), prob = as.vector(latent_sar)),
    nrow = n_sim,
    ncol = n_rings
  )

  observed_cv <- vapply(
    seq_len(n_sim),
    function(i) cv_safe(baseline_cases[i, ] / direct_sizes[i, ]),
    numeric(1)
  )
  empirical_icc <- vapply(
    seq_len(n_sim),
    function(i) estimate_binary_icc(baseline_cases[i, ], direct_sizes[i, ]),
    numeric(1)
  )
  latent_cv <- apply(latent_sar, 1L, cv_safe)

  structure(
    list(
      ring_row_index = ring_row_index,
      latent_sar = latent_sar,
      diagnostics = data.frame(
        latent_sar_cv = latent_cv,
        observed_sar_cv = observed_cv,
        empirical_icc = empirical_icc
      ),
      metadata = list(
        n_rings = n_rings,
        n_sim = n_sim,
        sar = sar,
        icc = icc,
        implied_latent_sar_cv = implied_latent_sar_cv(sar, icc),
        seed = seed
      )
    ),
    class = "ring_vax_epidemic_cache_v2"
  )
}


make_replicate_chunks <- function(n_sim, chunk_size) {
  starts <- seq.int(1L, n_sim, by = chunk_size)
  data.frame(
    chunk_id = seq_along(starts),
    sim_start = starts,
    sim_end = pmin(starts + chunk_size - 1L, n_sim),
    n_chunk = pmin(chunk_size, n_sim - starts + 1L)
  )
}


cache_chunk_path <- function(paths, n_rings, total_n_sim, sar, icc, chunk_id) {
  file.path(
    paths$cache_dir,
    paste0(
      "cache_", scenario_key(n_rings, total_n_sim, sar, icc),
      "_chunk_", sprintf("%04d", chunk_id), ".rds"
    )
  )
}


build_one_cache_chunk <- function(
  job, cache_grid, chunk_grid, rings, total_n_sim, config, paths
) {
  scenario <- cache_grid[job$scenario_id, , drop = FALSE]
  chunk <- chunk_grid[job$chunk_id, , drop = FALSE]
  n_rings <- as.integer(scenario$rings)
  sar <- as.numeric(scenario$SAR)
  icc <- as.numeric(scenario$ICC)
  file <- cache_chunk_path(
    paths, n_rings, total_n_sim, sar, icc, chunk$chunk_id
  )
  started <- proc.time()[["elapsed"]]

  if (file.exists(file)) {
    return(data.frame(
      scenario_id = job$scenario_id,
      chunk_id = chunk$chunk_id,
      status = "hit",
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      cache_mb = file.info(file)$size / 1024^2
    ))
  }

  seed <- as.integer((
    scenario_seed(config$seed, n_rings, sar, icc) +
      as.integer(chunk$chunk_id) * 104729
  ) %% .Machine$integer.max)
  cache <- generate_epidemic_cache(
    rings = rings,
    n_rings = n_rings,
    n_sim = as.integer(chunk$n_chunk),
    sar = sar,
    icc = icc,
    seed = seed
  )
  cache$metadata$total_n_sim <- total_n_sim
  cache$metadata$chunk_id <- as.integer(chunk$chunk_id)
  cache$metadata$sim_start <- as.integer(chunk$sim_start)
  cache$metadata$sim_end <- as.integer(chunk$sim_end)

  temp <- paste0(file, ".tmp_", Sys.getpid())
  saveRDS(cache, temp, compress = FALSE)
  if (!file.rename(temp, file)) stop("Could not publish cache chunk: ", file)

  data.frame(
    scenario_id = job$scenario_id,
    chunk_id = chunk$chunk_id,
    status = "generated",
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    cache_mb = file.info(file)$size / 1024^2
  )
}


build_one_cache <- function(row, rings, n_sim, config, paths) {
  n_rings <- as.integer(row[["rings"]])
  sar <- as.numeric(row[["SAR"]])
  icc <- as.numeric(row[["ICC"]])
  file <- cache_path(paths, n_rings, n_sim, sar, icc)
  started <- proc.time()[["elapsed"]]

  if (file.exists(file)) {
    return(data.frame(
      key = scenario_key(n_rings, n_sim, sar, icc),
      status = "hit",
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      cache_mb = file.info(file)$size / 1024^2
    ))
  }

  cache <- generate_epidemic_cache(
    rings = rings,
    n_rings = n_rings,
    n_sim = n_sim,
    sar = sar,
    icc = icc,
    seed = scenario_seed(config$seed, n_rings, sar, icc)
  )
  temp <- paste0(file, ".tmp_", Sys.getpid())
  saveRDS(cache, temp, compress = FALSE)
  if (!file.rename(temp, file)) stop("Could not atomically publish cache: ", file)

  data.frame(
    key = scenario_key(n_rings, n_sim, sar, icc),
    status = "generated",
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    cache_mb = file.info(file)$size / 1024^2
  )
}


# -------------------------
# TRIAL EVALUATION
# -------------------------

row_sums_where <- function(x, where) {
  rowSums(x * where)
}


run_design <- function(design, rings, cache, config, design_seed) {
  set.seed(design_seed)
  n_sim <- cache$metadata$n_sim
  n_rings <- cache$metadata$n_rings

  direct_sizes <- matrix(
    rings$direct_contacts[cache$ring_row_index],
    nrow = n_sim,
    ncol = n_rings
  )
  coc_sizes <- matrix(
    rings$coc_contacts[cache$ring_row_index],
    nrow = n_sim,
    ncol = n_rings
  )
  if (!isTRUE(config$include_coc)) coc_sizes[] <- 0L
  immediate <- make_allocation_matrix(
    n_sim,
    n_rings,
    as.numeric(design$allocation_ratio)
  )
  delayed <- !immediate

  p_immediate <- mean_protection_over_exposure_window(0)
  p_delayed <- mean_protection_over_exposure_window(as.numeric(design$delay))
  arm_protection <- ifelse(immediate, p_immediate, p_delayed)

  direct_risk <- cache$latent_sar *
    (1 - as.numeric(design$true_VE) * arm_protection)
  direct_cases <- matrix(
    rbinom(length(direct_risk), as.vector(direct_sizes), as.vector(direct_risk)),
    nrow = n_sim,
    ncol = n_rings
  )

  coc_cases <- matrix(0L, nrow = n_sim, ncol = n_rings)
  if (isTRUE(config$include_coc)) {
    coc_risk <- pmin(
      1,
      cache$latent_sar * as.numeric(design$coc_risk_multiplier)
    ) *
      (1 - as.numeric(design$true_VE) * arm_protection)
    coc_cases <- matrix(
      rbinom(length(coc_risk), as.vector(coc_sizes), as.vector(coc_risk)),
      nrow = n_sim,
      ncol = n_rings
    )
  }

  vaccine_cases <- row_sums_where(direct_cases, immediate)
  vaccine_n <- row_sums_where(direct_sizes, immediate)
  delayed_cases <- row_sums_where(direct_cases, delayed)
  delayed_n <- row_sums_where(direct_sizes, delayed)
  risk_vaccine <- vaccine_cases / vaccine_n
  risk_delayed <- delayed_cases / delayed_n
  ve_estimate <- 1 - risk_vaccine / risk_delayed
  ve_estimate[!is.finite(ve_estimate)] <- NA_real_

  expected_no_vaccine_direct <- rowSums(cache$latent_sar * direct_sizes)
  expected_no_vaccine_coc <- rowSums(
    pmin(
      1,
      cache$latent_sar * as.numeric(design$coc_risk_multiplier)
    ) * coc_sizes
  )
  direct_cases_prevented <- expected_no_vaccine_direct - rowSums(direct_cases)
  coc_cases_prevented <- expected_no_vaccine_coc - rowSums(coc_cases)
  total_cases_prevented <- direct_cases_prevented + coc_cases_prevented
  vaccine_doses <- row_sums_where(direct_sizes + coc_sizes, immediate)

  mean_m <- mean(direct_sizes)
  cv_m <- stats::sd(as.vector(direct_sizes)) / mean_m
  de_equal <- design_effect(mean_m, as.numeric(design$ICC), 0)
  de_unequal <- design_effect(mean_m, as.numeric(design$ICC), cv_m)
  total_direct_n <- rowSums(direct_sizes)

  true_operational_ve <- 1 -
    (1 - as.numeric(design$true_VE) * p_immediate) /
    (1 - as.numeric(design$true_VE) * p_delayed)
  q <- stats::quantile(ve_estimate, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  data.frame(
    design_id = as.integer(design$design_id),
    true_biological_VE = as.numeric(design$true_VE),
    true_operational_VE = true_operational_ve,
    SAR = as.numeric(design$SAR),
    ICC = as.numeric(design$ICC),
    coc_risk_multiplier = as.numeric(design$coc_risk_multiplier),
    rings = as.integer(design$rings),
    allocation_ratio = as.numeric(design$allocation_ratio),
    delay = as.integer(design$delay),

    implied_latent_SAR_CV = cache$metadata$implied_latent_sar_cv,
    simulated_latent_SAR_CV = mean(cache$diagnostics$latent_sar_cv, na.rm = TRUE),
    observed_SAR_CV = mean(cache$diagnostics$observed_sar_cv, na.rm = TRUE),
    empirical_ICC = mean(cache$diagnostics$empirical_icc, na.rm = TRUE),

    mean_cluster_size = mean_m,
    cluster_size_CV = cv_m,
    design_effect_equal_m = de_equal,
    design_effect_unequal_m = de_unequal,
    effective_sample_size = mean(total_direct_n / de_unequal),

    VE_immediate_vs_delayed = mean(ve_estimate, na.rm = TRUE),
    sd_VE = stats::sd(ve_estimate, na.rm = TRUE),
    VE_bias = mean(ve_estimate - true_operational_ve, na.rm = TRUE),
    VE_RMSE = sqrt(mean((ve_estimate - true_operational_ve)^2, na.rm = TRUE)),
    VE_empirical_lower = q[1],
    VE_empirical_upper = q[2],
    VE_empirical_width = q[2] - q[1],
    valid_VE_replicates = sum(is.finite(ve_estimate)),

    direct_cases_prevented = mean(direct_cases_prevented),
    coc_cases_prevented = mean(coc_cases_prevented),
    total_cases_prevented = mean(total_cases_prevented),
    vaccine_doses = mean(vaccine_doses),
    cases_per_1000_doses = mean(total_cases_prevented) /
      mean(vaccine_doses) * 1000
  )
}


simulate_design_chunk <- function(design, rings, cache, config, seed) {
  set.seed(seed)
  n_sim <- cache$metadata$n_sim
  n_rings <- cache$metadata$n_rings

  direct_sizes <- matrix(
    rings$direct_contacts[cache$ring_row_index],
    nrow = n_sim,
    ncol = n_rings
  )
  coc_sizes <- matrix(
    rings$coc_contacts[cache$ring_row_index],
    nrow = n_sim,
    ncol = n_rings
  )
  if (!isTRUE(config$include_coc)) coc_sizes[] <- 0L
  immediate <- make_allocation_matrix(
    n_sim, n_rings, as.numeric(design$allocation_ratio)
  )
  delayed <- !immediate

  p_immediate <- mean_protection_over_exposure_window(0)
  p_delayed <- mean_protection_over_exposure_window(as.numeric(design$delay))
  arm_protection <- ifelse(immediate, p_immediate, p_delayed)
  direct_risk <- cache$latent_sar *
    (1 - as.numeric(design$true_VE) * arm_protection)
  direct_cases <- matrix(
    rbinom(length(direct_risk), as.vector(direct_sizes), as.vector(direct_risk)),
    nrow = n_sim,
    ncol = n_rings
  )

  coc_cases <- matrix(0L, nrow = n_sim, ncol = n_rings)
  if (isTRUE(config$include_coc)) {
    coc_risk <- pmin(
      1,
      cache$latent_sar * as.numeric(design$coc_risk_multiplier)
    ) *
      (1 - as.numeric(design$true_VE) * arm_protection)
    coc_cases <- matrix(
      rbinom(length(coc_risk), as.vector(coc_sizes), as.vector(coc_risk)),
      nrow = n_sim,
      ncol = n_rings
    )
  }

  vaccine_cases <- row_sums_where(direct_cases, immediate)
  vaccine_n <- row_sums_where(direct_sizes, immediate)
  delayed_cases <- row_sums_where(direct_cases, delayed)
  delayed_n <- row_sums_where(direct_sizes, delayed)
  ve_estimate <- 1 - (vaccine_cases / vaccine_n) / (delayed_cases / delayed_n)
  ve_estimate[!is.finite(ve_estimate)] <- NA_real_

  expected_no_vaccine_direct <- rowSums(cache$latent_sar * direct_sizes)
  expected_no_vaccine_coc <- rowSums(
    pmin(
      1,
      cache$latent_sar * as.numeric(design$coc_risk_multiplier)
    ) * coc_sizes
  )
  direct_prevented <- expected_no_vaccine_direct - rowSums(direct_cases)
  coc_prevented <- expected_no_vaccine_coc - rowSums(coc_cases)

  list(
    metrics = data.frame(
      VE = ve_estimate,
      direct_cases_prevented = direct_prevented,
      coc_cases_prevented = coc_prevented,
      total_cases_prevented = direct_prevented + coc_prevented,
      vaccine_doses = row_sums_where(direct_sizes + coc_sizes, immediate),
      total_direct_n = rowSums(direct_sizes),
      latent_sar_cv = cache$diagnostics$latent_sar_cv,
      observed_sar_cv = cache$diagnostics$observed_sar_cv,
      empirical_icc = cache$diagnostics$empirical_icc
    ),
    cluster_sufficient = c(
      n = length(direct_sizes),
      sum = sum(direct_sizes),
      sum_squares = sum(as.numeric(direct_sizes)^2)
    ),
    metadata = list(
      implied_latent_sar_cv = cache$metadata$implied_latent_sar_cv,
      chunk_id = cache$metadata$chunk_id,
      n_sim = n_sim
    )
  )
}


design_chunk_file <- function(paths, design_id, chunk_id) {
  file.path(
    paths$design_chunk_dir,
    paste0(
      "design_", design_id,
      "_chunk_", sprintf("%04d", chunk_id), ".rds"
    )
  )
}


run_one_design_chunk <- function(
  job, design_grid, chunk_grid, rings, total_n_sim, config, paths
) {
  design <- design_grid[job$design_index, , drop = FALSE]
  chunk <- chunk_grid[job$chunk_id, , drop = FALSE]
  output_file <- design_chunk_file(paths, design$design_id, chunk$chunk_id)
  started <- proc.time()[["elapsed"]]

  if (config$resume && file.exists(output_file)) {
    return(data.frame(
      design_id = design$design_id,
      chunk_id = chunk$chunk_id,
      status = "hit",
      elapsed_seconds = proc.time()[["elapsed"]] - started
    ))
  }

  cache_file <- cache_chunk_path(
    paths,
    design$rings,
    total_n_sim,
    design$SAR,
    design$ICC,
    chunk$chunk_id
  )
  if (!file.exists(cache_file)) stop("Required cache chunk is missing: ", cache_file)
  cache <- readRDS(cache_file)
  seed <- as.integer((
    config$seed + as.integer(design$design_id) * 10007 +
      as.integer(chunk$chunk_id) * 13007
  ) %% .Machine$integer.max)
  result <- simulate_design_chunk(design, rings, cache, config, seed)
  result$elapsed_seconds <- proc.time()[["elapsed"]] - started

  temp <- paste0(output_file, ".tmp_", Sys.getpid())
  saveRDS(result, temp, compress = FALSE)
  if (!file.rename(temp, output_file)) {
    stop("Could not publish design chunk: ", output_file)
  }

  data.frame(
    design_id = design$design_id,
    chunk_id = chunk$chunk_id,
    status = "generated",
    elapsed_seconds = result$elapsed_seconds
  )
}


aggregate_one_design <- function(
  i, design_grid, chunk_grid, config, paths
) {
  design <- design_grid[i, , drop = FALSE]
  checkpoint <- checkpoint_file(paths, design$design_id)
  if (config$resume && file.exists(checkpoint)) return(readRDS(checkpoint))
  started <- proc.time()[["elapsed"]]

  files <- vapply(
    chunk_grid$chunk_id,
    function(chunk_id) design_chunk_file(paths, design$design_id, chunk_id),
    character(1)
  )
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Design ", design$design_id, " is missing ", length(missing), " chunks.")
  }
  chunks <- lapply(files, readRDS)
  metrics <- do.call(rbind, lapply(chunks, `[[`, "metrics"))
  sufficient <- Reduce(`+`, lapply(chunks, `[[`, "cluster_sufficient"))

  mean_m <- sufficient[["sum"]] / sufficient[["n"]]
  cluster_variance <- (
    sufficient[["sum_squares"]] - sufficient[["sum"]]^2 / sufficient[["n"]]
  ) / (sufficient[["n"]] - 1)
  cv_m <- sqrt(max(0, cluster_variance)) / mean_m
  de_equal <- design_effect(mean_m, as.numeric(design$ICC), 0)
  de_unequal <- design_effect(mean_m, as.numeric(design$ICC), cv_m)

  p_immediate <- mean_protection_over_exposure_window(0)
  p_delayed <- mean_protection_over_exposure_window(as.numeric(design$delay))
  true_operational_ve <- 1 -
    (1 - as.numeric(design$true_VE) * p_immediate) /
    (1 - as.numeric(design$true_VE) * p_delayed)
  q <- stats::quantile(metrics$VE, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  result <- data.frame(
    design_id = as.integer(design$design_id),
    true_biological_VE = as.numeric(design$true_VE),
    true_operational_VE = true_operational_ve,
    SAR = as.numeric(design$SAR),
    ICC = as.numeric(design$ICC),
    coc_risk_multiplier = as.numeric(design$coc_risk_multiplier),
    rings = as.integer(design$rings),
    allocation_ratio = as.numeric(design$allocation_ratio),
    delay = as.integer(design$delay),

    implied_latent_SAR_CV = chunks[[1]]$metadata$implied_latent_sar_cv,
    simulated_latent_SAR_CV = mean(metrics$latent_sar_cv, na.rm = TRUE),
    observed_SAR_CV = mean(metrics$observed_sar_cv, na.rm = TRUE),
    empirical_ICC = mean(metrics$empirical_icc, na.rm = TRUE),

    mean_cluster_size = mean_m,
    cluster_size_CV = cv_m,
    design_effect_equal_m = de_equal,
    design_effect_unequal_m = de_unequal,
    effective_sample_size = mean(metrics$total_direct_n / de_unequal),

    VE_immediate_vs_delayed = mean(metrics$VE, na.rm = TRUE),
    sd_VE = stats::sd(metrics$VE, na.rm = TRUE),
    VE_bias = mean(metrics$VE - true_operational_ve, na.rm = TRUE),
    VE_RMSE = sqrt(mean((metrics$VE - true_operational_ve)^2, na.rm = TRUE)),
    VE_empirical_lower = q[1],
    VE_empirical_upper = q[2],
    VE_empirical_width = q[2] - q[1],
    valid_VE_replicates = sum(is.finite(metrics$VE)),

    direct_cases_prevented = mean(metrics$direct_cases_prevented),
    coc_cases_prevented = mean(metrics$coc_cases_prevented),
    total_cases_prevented = mean(metrics$total_cases_prevented),
    vaccine_doses = mean(metrics$vaccine_doses),
    cases_per_1000_doses = mean(metrics$total_cases_prevented) /
      mean(metrics$vaccine_doses) * 1000,
    simulation_elapsed_seconds = sum(vapply(chunks, `[[`, numeric(1), "elapsed_seconds")),
    aggregation_elapsed_seconds = proc.time()[["elapsed"]] - started
  )

  temp <- paste0(checkpoint, ".tmp_", Sys.getpid())
  saveRDS(result, temp)
  if (!file.rename(temp, checkpoint)) stop("Could not publish checkpoint: ", checkpoint)
  result
}


checkpoint_file <- function(paths, design_id) {
  file.path(paths$checkpoint_dir, paste0("design_", design_id, ".rds"))
}


run_one_design <- function(i, design_grid, rings, n_sim, config, paths) {
  design <- design_grid[i, , drop = FALSE]
  checkpoint <- checkpoint_file(paths, design$design_id)
  if (config$resume && file.exists(checkpoint)) return(readRDS(checkpoint))

  started <- proc.time()[["elapsed"]]
  cache_file <- cache_path(
    paths, design$rings, n_sim, design$SAR, design$ICC
  )
  if (!file.exists(cache_file)) stop("Required cache is missing: ", cache_file)
  cache <- readRDS(cache_file)

  result <- run_design(
    design = design,
    rings = rings,
    cache = cache,
    config = config,
    design_seed = as.integer((config$seed + design$design_id * 10007) %% .Machine$integer.max)
  )
  result$elapsed_seconds <- proc.time()[["elapsed"]] - started

  temp <- paste0(checkpoint, ".tmp_", Sys.getpid())
  saveRDS(result, temp)
  if (!file.rename(temp, checkpoint)) stop("Could not publish checkpoint: ", checkpoint)
  result
}


# -------------------------
# DESIGN GRID AND EXECUTION
# -------------------------

make_design_grid <- function(config) {
  coc_values <- if (isTRUE(config$include_coc)) {
    config$coc_risk_multipliers
  } else {
    0
  }
  grid <- expand.grid(
    true_VE = config$true_ve_values,
    SAR = config$sar_values,
    ICC = config$icc_values,
    coc_risk_multiplier = coc_values,
    rings = config$ring_numbers,
    allocation_ratio = config$allocation_ratios,
    delay = config$delays,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$design_id <- seq_len(nrow(grid))

  if (config$debug_mode) {
    if (isTRUE(config$run_analysis) && config$debug_n_designs < nrow(grid)) {
      stop(
        "Automated comparative analysis requires a complete factorial grid. ",
        "The configured grid has ", nrow(grid), " designs, but ",
        "DEBUG_N_DESIGNS is ", config$debug_n_designs, ". ",
        "Set DEBUG_N_DESIGNS <- ", nrow(grid),
        " or narrow the scenario vectors."
      )
    }
    # Keep broad ICC/ring coverage instead of taking only the first grid rows.
    ordering <- order(
      grid$ICC,
      grid$coc_risk_multiplier,
      grid$rings,
      grid$SAR,
      grid$true_VE
    )
    evenly_spaced <- unique(round(seq(1, length(ordering), length.out = min(
      config$debug_n_designs, length(ordering)
    ))))
    grid <- grid[ordering[evenly_spaced], , drop = FALSE]
  }
  rownames(grid) <- NULL
  grid
}


parallel_lapply <- function(x, fun, config) {
  if (!config$parallel) return(lapply(x, fun))
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Parallel execution requires the future and future.apply packages.")
  }
  future::plan(future::multisession, workers = config$n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  future.apply::future_lapply(x, fun, future.seed = TRUE)
}


pareto_flag <- function(df) {
  keep <- rep(TRUE, nrow(df))
  for (i in seq_len(nrow(df))) {
    dominated <- any(
      df$total_cases_prevented >= df$total_cases_prevented[i] &
        df$sd_VE <= df$sd_VE[i] &
        df$vaccine_doses <= df$vaccine_doses[i] &
        (
          df$total_cases_prevented > df$total_cases_prevented[i] |
            df$sd_VE < df$sd_VE[i] |
            df$vaccine_doses < df$vaccine_doses[i]
        ),
      na.rm = TRUE
    )
    if (dominated) keep[i] <- FALSE
  }
  keep
}


add_pareto_flags <- function(output) {
  output$pareto_optimal <- FALSE
  groups <- split(
    seq_len(nrow(output)),
    interaction(
      output$true_biological_VE,
      output$SAR,
      output$ICC,
      output$coc_risk_multiplier,
      output$delay,
      drop = TRUE
    )
  )
  for (g in groups) output$pareto_optimal[g] <- pareto_flag(output[g, , drop = FALSE])
  output
}


add_pareto_flags_parallel <- function(output, config) {
  output$pareto_optimal <- FALSE
  groups <- split(
    seq_len(nrow(output)),
    interaction(
      output$true_biological_VE,
      output$SAR,
      output$ICC,
      output$coc_risk_multiplier,
      output$delay,
      drop = TRUE
    )
  )
  flags <- parallel_lapply(
    groups,
    function(g) list(rows = g, keep = pareto_flag(output[g, , drop = FALSE])),
    config
  )
  for (x in flags) output$pareto_optimal[x$rows] <- x$keep
  output
}


# -------------------------
# AUTOMATED POST-SIMULATION ANALYSIS
# -------------------------

check_analysis_dependencies <- function(config) {
  if (!isTRUE(config$run_analysis)) return(invisible(TRUE))
  required <- c("dplyr", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "RUN_ANALYSIS is TRUE, but these packages are missing: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}


pareto_precision_prevention <- function(df) {
  keep <- rep(TRUE, nrow(df))
  for (i in seq_len(nrow(df))) {
    dominated <- any(
      df$total_cases_prevented >= df$total_cases_prevented[i] &
        df$sd_VE <= df$sd_VE[i] &
        (
          df$total_cases_prevented > df$total_cases_prevented[i] |
            df$sd_VE < df$sd_VE[i]
        ),
      na.rm = TRUE
    )
    if (dominated) keep[i] <- FALSE
  }
  keep
}


select_pragmatic_recommendations <- function(results, threshold) {
  scenario_keys <- c(
    "true_biological_VE", "SAR", "ICC", "coc_risk_multiplier", "delay"
  )
  groups <- split(
    seq_len(nrow(results)),
    do.call(interaction, c(results[scenario_keys], list(drop = TRUE)))
  )

  selected <- lapply(groups, function(rows) {
    scenario <- results[rows, , drop = FALSE]
    qualifying <- scenario[
      is.finite(scenario$sd_VE) & scenario$sd_VE <= threshold,
      ,
      drop = FALSE
    ]

    if (nrow(qualifying)) {
      minimum_rings <- min(qualifying$rings)
      candidates <- qualifying[
        qualifying$rings == minimum_rings,
        ,
        drop = FALSE
      ]
      candidates <- candidates[order(
        -candidates$total_cases_prevented,
        candidates$sd_VE,
        candidates$vaccine_doses
      ), , drop = FALSE]
      recommendation <- candidates[1L, , drop = FALSE]
      recommendation$priority <- "Pragmatic"
      recommendation$recommendation_status <- "Meets precision threshold"
      recommendation$minimum_qualifying_rings <- minimum_rings
      recommendation$meets_precision_threshold <- TRUE
    } else {
      candidates <- scenario[is.finite(scenario$sd_VE), , drop = FALSE]
      candidates <- candidates[order(
        candidates$sd_VE,
        -candidates$total_cases_prevented,
        candidates$vaccine_doses
      ), , drop = FALSE]
      recommendation <- candidates[1L, , drop = FALSE]
      recommendation$priority <- "Fallback"
      recommendation$recommendation_status <-
        "No design met threshold; best available precision"
      recommendation$minimum_qualifying_rings <- NA_integer_
      recommendation$meets_precision_threshold <- FALSE
    }

    recommendation$decision_precision_threshold <- threshold
    recommendation
  })

  dplyr::bind_rows(selected)
}


save_analysis_plot <- function(plot, filename, config, width = 11, height = 7) {
  if (!isTRUE(config$run_figures)) return(invisible(FALSE))
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = config$figure_dpi
  )
}


run_post_analysis <- function(results, analysis_dir, config) {
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  scenario_keys <- c(
    "true_biological_VE", "SAR", "ICC", "coc_risk_multiplier", "delay"
  )

  # Two-objective Pareto analysis from the legacy analysis: maximise cases
  # prevented and minimise SD(VE), evaluated within epidemiological scenarios.
  results$pareto_precision_prevention <- FALSE
  groups <- split(
    seq_len(nrow(results)),
    do.call(interaction, c(results[scenario_keys], list(drop = TRUE)))
  )
  for (g in groups) {
    results$pareto_precision_prevention[g] <-
      pareto_precision_prevention(results[g, , drop = FALSE])
  }
  pareto_only <- results[results$pareto_precision_prevention, , drop = FALSE]
  if (!nrow(pareto_only)) stop("Post-analysis found no Pareto-optimal designs.")

  pareto_summary <- pareto_only |>
    dplyr::group_by(
      true_biological_VE, SAR, ICC, coc_risk_multiplier, delay
    ) |>
    dplyr::summarise(
      n_pareto_designs = dplyr::n(),
      allocation_ratios = paste(sort(unique(allocation_ratio)), collapse = ", "),
      rings_range = paste(min(rings), max(rings), sep = "-"),
      cases_prevented_range = paste(
        round(min(total_cases_prevented)),
        round(max(total_cases_prevented)),
        sep = "-"
      ),
      sd_VE_range = paste(round(min(sd_VE), 3), round(max(sd_VE), 3), sep = "-"),
      .groups = "drop"
    )
  write.csv(
    pareto_only,
    file.path(analysis_dir, "pareto_frontier_designs.csv"),
    row.names = FALSE
  )
  write.csv(
    pareto_summary,
    file.path(analysis_dir, "pareto_frontier_summary.csv"),
    row.names = FALSE
  )

  p_pareto <- ggplot2::ggplot(
    results,
    ggplot2::aes(
      x = -sd_VE,
      y = total_cases_prevented,
      colour = factor(allocation_ratio),
      size = rings
    )
  ) +
    ggplot2::geom_point(alpha = 0.45) +
    ggplot2::geom_point(
      data = pareto_only,
      shape = 21,
      fill = "white",
      colour = "black",
      size = 3.5
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(true_biological_VE, ICC, coc_risk_multiplier),
      cols = ggplot2::vars(SAR),
      labeller = ggplot2::label_both
    ) +
    ggplot2::labs(
      x = "Inferential precision (higher = lower SD(VE))",
      y = "Total cases prevented",
      colour = "Immediate:delayed allocation",
      size = "Rings",
      title = "Pareto frontier: prevention versus inferential precision"
    ) +
    ggplot2::theme_bw()
  save_analysis_plot(
    p_pareto,
    file.path(analysis_dir, "figure_pareto_frontier.png"),
    config,
    height = 15
  )

  # Allocation comparisons are matched on VE, SAR, ICC, CoC risk, rings and
  # delay. CoC assumptions must never be merged into the same comparison.
  allocation_keys <- c(
    "true_biological_VE", "SAR", "ICC", "coc_risk_multiplier",
    "rings", "delay"
  )
  baseline <- results[results$allocation_ratio == 1, c(
    allocation_keys, "total_cases_prevented", "sd_VE", "vaccine_doses"
  ), drop = FALSE]
  names(baseline)[match(
    c("total_cases_prevented", "sd_VE", "vaccine_doses"), names(baseline)
  )] <- c("baseline_cases", "baseline_sd_VE", "baseline_doses")
  allocation_gain <- merge(
    results[results$allocation_ratio > 1, , drop = FALSE],
    baseline,
    by = allocation_keys,
    all.x = TRUE,
    sort = FALSE
  )
  allocation_gain$additional_cases <-
    allocation_gain$total_cases_prevented - allocation_gain$baseline_cases
  allocation_gain$precision_cost <-
    allocation_gain$sd_VE - allocation_gain$baseline_sd_VE
  allocation_gain$additional_doses <-
    allocation_gain$vaccine_doses - allocation_gain$baseline_doses
  allocation_gain <- allocation_gain[complete.cases(
    allocation_gain[c("baseline_cases", "additional_cases", "precision_cost")]
  ), , drop = FALSE]
  write.csv(
    allocation_gain,
    file.path(analysis_dir, "allocation_tradeoff.csv"),
    row.names = FALSE
  )

  if (nrow(allocation_gain)) {
    p_allocation <- ggplot2::ggplot(
      allocation_gain,
      ggplot2::aes(
        x = precision_cost,
        y = additional_cases,
        colour = factor(allocation_ratio),
        size = rings
      )
    ) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(true_biological_VE, ICC, coc_risk_multiplier),
        cols = ggplot2::vars(SAR),
        labeller = ggplot2::label_both
      ) +
      ggplot2::labs(
        x = "Precision cost: increase in SD(VE)",
        y = "Additional cases prevented versus 1:1",
        colour = "Allocation ratio",
        size = "Rings",
        title = "Marginal benefit and precision cost of unequal allocation"
      ) +
      ggplot2::theme_bw()
    save_analysis_plot(
      p_allocation,
      file.path(analysis_dir, "figure_allocation_tradeoff.png"),
      config,
      height = 15
    )
  }

  # Decision recommendations within each VE/SAR/ICC/CoC-risk/delay scenario.
  # pragmatic recommendation first finds the smallest ring count meeting the
  # configured precision threshold, then maximises cases prevented among the
  # qualifying allocations at that ring count. This prevents the recommendation
  # from drifting to the largest simulated trial merely because it is most
  # precise. If no design qualifies, report the best-precision fallback clearly.
  ranked <- results |>
    dplyr::group_by(
      true_biological_VE, SAR, ICC, coc_risk_multiplier, delay
    ) |>
    dplyr::mutate(
      best_precision = min(sd_VE, na.rm = TRUE),
      best_cases = max(total_cases_prevented, na.rm = TRUE)
    ) |>
    dplyr::ungroup()
  pragmatic_priority <- select_pragmatic_recommendations(
    results,
    config$decision_precision_threshold
  )
  prevention_priority <- ranked |>
    dplyr::filter(total_cases_prevented >= best_cases * 0.95) |>
    dplyr::group_by(
      true_biological_VE, SAR, ICC, coc_risk_multiplier, delay
    ) |>
    dplyr::slice_min(sd_VE, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      priority = "Prevention",
      recommendation_status = "Within 5% of maximum cases prevented",
      minimum_qualifying_rings = NA_integer_,
      meets_precision_threshold = sd_VE <= config$decision_precision_threshold,
      decision_precision_threshold = config$decision_precision_threshold
    )
  recommendations <- dplyr::bind_rows(pragmatic_priority, prevention_priority) |>
    dplyr::select(
      priority, recommendation_status, decision_precision_threshold,
      meets_precision_threshold, minimum_qualifying_rings,
      true_biological_VE, SAR, ICC, coc_risk_multiplier,
      rings, allocation_ratio, delay,
      total_cases_prevented, sd_VE, VE_empirical_width,
      effective_sample_size, vaccine_doses
    )
  write.csv(
    recommendations,
    file.path(analysis_dir, "allocation_recommendation_table.csv"),
    row.names = FALSE
  )

  # Minimum rings meeting precision thresholds, including explicit NA rows.
  rings_required <- do.call(rbind, lapply(config$precision_thresholds, function(threshold) {
    x <- results |>
      dplyr::group_by(
        true_biological_VE, SAR, ICC, coc_risk_multiplier,
        allocation_ratio, delay
      ) |>
      dplyr::summarise(
        minimum_rings = if (any(sd_VE <= threshold, na.rm = TRUE)) {
          min(rings[sd_VE <= threshold], na.rm = TRUE)
        } else {
          NA_integer_
        },
        achieved_sd_VE = if (any(sd_VE <= threshold, na.rm = TRUE)) {
          first_qualifying_ring <- min(rings[sd_VE <= threshold], na.rm = TRUE)
          min(sd_VE[rings == first_qualifying_ring], na.rm = TRUE)
        } else {
          NA_real_
        },
        .groups = "drop"
      )
    x$precision_threshold <- threshold
    x
  }))
  write.csv(
    rings_required,
    file.path(analysis_dir, "rings_required_for_precision.csv"),
    row.names = FALSE
  )

  p_precision <- ggplot2::ggplot(
    results,
    ggplot2::aes(
      x = rings,
      y = sd_VE,
      colour = factor(allocation_ratio),
      group = allocation_ratio
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(
      yintercept = config$precision_thresholds,
      linetype = "dashed",
      colour = "grey55"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(true_biological_VE, ICC, coc_risk_multiplier),
      cols = ggplot2::vars(SAR),
      labeller = ggplot2::label_both
    ) +
    ggplot2::labs(
      x = "Number of rings",
      y = "SD of operational VE estimate",
      colour = "Allocation ratio",
      title = "VE precision by number of rings and ICC"
    ) +
    ggplot2::theme_bw()
  save_analysis_plot(
    p_precision,
    file.path(analysis_dir, "figure_rings_vs_precision.png"),
    config,
    height = 15
  )

  ring_efficiency <- results |>
    dplyr::arrange(
      true_biological_VE, SAR, ICC, coc_risk_multiplier,
      delay, allocation_ratio, rings
    ) |>
    dplyr::group_by(
      true_biological_VE, SAR, ICC, coc_risk_multiplier,
      delay, allocation_ratio
    ) |>
    dplyr::mutate(
      additional_cases_per_ring =
        (total_cases_prevented - dplyr::lag(total_cases_prevented)) /
        (rings - dplyr::lag(rings))
    ) |>
    dplyr::ungroup()
  write.csv(
    ring_efficiency,
    file.path(analysis_dir, "marginal_cases_per_additional_ring.csv"),
    row.names = FALSE
  )

  # ICC-specific replacements for the legacy SAR-CV sensitivity figures.
  if (length(unique(results$ICC)) > 1L) {
    p_icc <- ggplot2::ggplot(
      results,
      ggplot2::aes(
        x = ICC,
        y = sd_VE,
        colour = factor(allocation_ratio),
        linetype = factor(rings),
        group = interaction(allocation_ratio, rings)
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(
        rows = ggplot2::vars(true_biological_VE, coc_risk_multiplier),
        cols = ggplot2::vars(SAR),
        labeller = ggplot2::label_both
      ) +
      ggplot2::labs(
        x = "Assumed ICC",
        y = "SD of operational VE estimate",
        colour = "Allocation ratio",
        linetype = "Rings",
        title = "Precision loss under increasing clustering"
      ) +
      ggplot2::theme_bw()
    save_analysis_plot(
      p_icc,
      file.path(analysis_dir, "figure_ICC_precision.png"),
      config,
      height = 10
    )

    sar_diagnostics <- unique(results[c(
      "SAR", "ICC", "rings", "empirical_ICC", "observed_SAR_CV",
      "implied_latent_SAR_CV", "simulated_latent_SAR_CV"
    )])
    write.csv(
      sar_diagnostics,
      file.path(analysis_dir, "ICC_SAR_diagnostics.csv"),
      row.names = FALSE
    )
    p_sar <- ggplot2::ggplot(
      sar_diagnostics,
      ggplot2::aes(
        x = ICC,
        y = observed_SAR_CV,
        colour = factor(rings),
        group = rings
      )
    ) +
      ggplot2::geom_hline(yintercept = 2, linetype = "dashed", colour = "grey40") +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~SAR, labeller = ggplot2::label_both) +
      ggplot2::labs(
        x = "Assumed ICC",
        y = "Simulated observed SAR-CV",
        colour = "Rings",
        title = "Observed SAR variation generated by ICC assumptions",
        subtitle = "Dashed line marks the approximate empirical CV of 200%"
      ) +
      ggplot2::theme_bw()
    save_analysis_plot(
      p_sar,
      file.path(analysis_dir, "figure_ICC_observed_SAR_CV.png"),
      config
    )
  }

  if (length(unique(results$delay)) > 1L) {
    p_delay <- ggplot2::ggplot(
      results,
      ggplot2::aes(
        x = delay,
        y = total_cases_prevented,
        colour = factor(allocation_ratio),
        group = interaction(allocation_ratio, rings, coc_risk_multiplier)
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(
        rows = ggplot2::vars(true_biological_VE, ICC, coc_risk_multiplier),
        cols = ggplot2::vars(SAR),
        labeller = ggplot2::label_both
      ) +
      ggplot2::labs(
        x = "Delayed-arm vaccination day",
        y = "Total cases prevented",
        colour = "Allocation ratio",
        title = "Impact of vaccination delay"
      ) +
      ggplot2::theme_bw()
    save_analysis_plot(
      p_delay,
      file.path(analysis_dir, "figure_delay_effect.png"),
      config,
      height = 15
    )
  }

  event_columns <- intersect(
    c(
      "true_biological_VE", "SAR", "ICC", "rings", "allocation_ratio", "delay",
      "coc_risk_multiplier",
      "direct_cases_prevented", "coc_cases_prevented", "total_cases_prevented",
      "vaccine_doses", "cases_per_1000_doses"
    ),
    names(results)
  )
  write.csv(
    results[event_columns],
    file.path(analysis_dir, "available_event_outputs.csv"),
    row.names = FALSE
  )

  rapid_columns <- intersect(
    c(
      "true_biological_VE", "true_operational_VE", "SAR", "ICC",
      "coc_risk_multiplier",
      "rings", "allocation_ratio", "delay", "VE_immediate_vs_delayed",
      "sd_VE", "VE_empirical_width", "effective_sample_size",
      "total_cases_prevented", "vaccine_doses", "cases_per_1000_doses"
    ),
    names(results)
  )
  rapid_table <- results[order(
    results$true_biological_VE, results$SAR, results$ICC,
    results$coc_risk_multiplier, results$delay,
    results$rings, results$allocation_ratio
  ), rapid_columns, drop = FALSE]
  write.csv(
    rapid_table,
    file.path(analysis_dir, "rapid_decision_table.csv"),
    row.names = FALSE
  )

  message("Automated analysis completed: ", analysis_dir)
  invisible(list(
    pareto = pareto_only,
    recommendations = recommendations,
    rings_required = rings_required,
    rapid_decision_table = rapid_table
  ))
}


main <- function(config = default_config()) {
  validate_config(config)
  check_analysis_dependencies(config)
  n_sim <- if (config$debug_mode) config$debug_n_sim else config$n_sim
  design_grid <- make_design_grid(config)
  paths <- initialise_run(config, design_grid)
  rings <- load_empirical_rings(config$empirical_ring_file)

  chunk_size <- min(config$replicate_chunk_size, n_sim)
  chunk_grid <- make_replicate_chunks(n_sim, chunk_size)

  cache_grid <- unique(design_grid[c("rings", "SAR", "ICC")])
  rownames(cache_grid) <- NULL

  # Stage 1: every scenario-replicate chunk is an independent cache job.
  cache_jobs <- expand.grid(
    scenario_id = seq_len(nrow(cache_grid)),
    chunk_id = chunk_grid$chunk_id,
    KEEP.OUT.ATTRS = FALSE
  )
  message(
    "Stage 1/3: preparing ", nrow(cache_jobs),
    " epidemic cache chunks across ",
    if (config$parallel) config$n_workers else 1L, " worker(s)..."
  )
  stage_started <- proc.time()[["elapsed"]]
  cache_timings <- parallel_lapply(
    seq_len(nrow(cache_jobs)),
    function(i) build_one_cache_chunk(
      cache_jobs[i, ], cache_grid, chunk_grid, rings, n_sim, config, paths
    ),
    config
  )
  cache_timings <- do.call(rbind, cache_timings)
  cache_timings$stage_wall_seconds <- proc.time()[["elapsed"]] - stage_started
  write.csv(cache_timings, file.path(paths$log_dir, "cache_timings.csv"), row.names = FALSE)

  # Stage 2: flatten design and replicate chunk dimensions into one queue.
  # This parallelises both designs and replicates without nested worker pools.
  design_jobs <- expand.grid(
    design_index = seq_len(nrow(design_grid)),
    chunk_id = chunk_grid$chunk_id,
    KEEP.OUT.ATTRS = FALSE
  )
  message(
    "Stage 2/3: simulating ", nrow(design_jobs),
    " design chunks across the worker pool..."
  )
  stage_started <- proc.time()[["elapsed"]]
  design_timings <- parallel_lapply(
    seq_len(nrow(design_jobs)),
    function(i) run_one_design_chunk(
      design_jobs[i, ], design_grid, chunk_grid, rings, n_sim, config, paths
    ),
    config
  )
  design_timings <- do.call(rbind, design_timings)
  design_timings$stage_wall_seconds <- proc.time()[["elapsed"]] - stage_started
  write.csv(design_timings, file.path(paths$log_dir, "design_chunk_timings.csv"), row.names = FALSE)

  # Stage 3: each design is independently aggregated and tested. Pareto
  # comparisons are subsequently parallelised within epidemiological groups.
  message(
    "Stage 3/3: aggregating and testing ", nrow(design_grid),
    " completed designs across the worker pool..."
  )
  stage_started <- proc.time()[["elapsed"]]
  results <- parallel_lapply(
    seq_len(nrow(design_grid)),
    function(i) aggregate_one_design(i, design_grid, chunk_grid, config, paths),
    config
  )
  output <- do.call(rbind, results)
  output <- add_pareto_flags_parallel(output, config)
  aggregation_wall_seconds <- proc.time()[["elapsed"]] - stage_started
  write.csv(
    data.frame(
      designs = nrow(design_grid),
      replicate_chunks = nrow(chunk_grid),
      workers = if (config$parallel) config$n_workers else 1L,
      stage_wall_seconds = aggregation_wall_seconds
    ),
    file.path(paths$log_dir, "aggregation_timing.csv"),
    row.names = FALSE
  )

  output_file <- file.path(
    paths$output_dir,
    paste0("ring_vax_icc_operating_characteristics_", config$run_id, ".csv")
  )
  pareto_file <- file.path(
    paths$output_dir,
    paste0("ring_vax_icc_pareto_frontier_", config$run_id, ".csv")
  )
  write.csv(output, output_file, row.names = FALSE)
  write.csv(output[output$pareto_optimal, , drop = FALSE], pareto_file, row.names = FALSE)
  saveRDS(output, sub("\\.csv$", ".rds", output_file))

  if (isTRUE(config$run_analysis)) {
    message("Running automated ICC design analysis...")
    run_post_analysis(output, paths$analysis_dir, config)
  }

  message("Completed. Results: ", output_file)
  invisible(output)
}


if (sys.nframe() == 0L) main()
