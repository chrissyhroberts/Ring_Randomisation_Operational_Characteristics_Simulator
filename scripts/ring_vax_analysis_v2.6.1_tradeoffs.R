# ============================================================
# Ring-vaccination comprehensive trade-off analysis v2.6.1
# ============================================================
#
# Analysis only: this script never runs or resumes simulations.
# It reads the completed comprehensive operating-characteristics table and
# reports how much robustness is obtained for different ring budgets.
#
# Primary principle:
#   1,500 rings may maximise precision and total prevention simply because it
#   is the largest trial. It is therefore reported as maximum robustness, not
#   as an automatic recommendation. Resource-constrained and efficiency-based
#   alternatives are reported alongside it.
# ============================================================


# ============================================================
# USER CONFIGURATION -- EDIT THIS BLOCK ONLY
# ============================================================

INPUT_RESULTS_FILE <- paste0(
  "runs/ring_vax_icc_comprehensive003/outputs/",
  "ring_vax_icc_operating_characteristics_comprehensive003.csv"
)

OUTPUT_DIR <- paste0(
  "runs/ring_vax_icc_comprehensive003/",
  "analysis_v2.6.1_tradeoffs"
)

N_SIM <- 2000L
FIGURE_DPI <- 300L
RUN_FIGURES <- TRUE

# Inferential adequacy criteria
PRECISION_THRESHOLD <- 0.15
MAX_ABS_BIAS <- 0.05
MIN_VALID_FRACTION <- 0.95

# The primary analysis always uses the complete grid. The focus region is a
# separately labelled sensitivity analysis, not a probability distribution.
FOCUS_TRUE_VE_VALUES <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
FOCUS_SAR_VALUES <- c(0.02, 0.04, 0.06)
FOCUS_ICC_VALUES <- c(0.05, 0.10, 0.15)
FOCUS_COC_RISK_MULTIPLIERS <- c(0.05, 0.10, 0.25)

# Operational questions to answer
RING_CAPS <- c(500L, 750L, 900L, 1000L, 1100L, 1200L, 1500L)
ABSOLUTE_COVERAGE_TARGETS_PERCENT <- c(50, 55, 60, 65, 70, 75, 80, 85)
ROBUSTNESS_RETENTION_TARGETS <- c(0.80, 0.90, 0.95)
PREVENTION_COVERAGE_FLOORS_PERCENT <- c(45, 50, 55, 60, 65)

# ============================================================
# END USER CONFIGURATION
# ============================================================


args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) INPUT_RESULTS_FILE <- args[[1L]]
if (length(args) >= 2L) OUTPUT_DIR <- args[[2L]]


check_dependencies <- function() {
  required <- c("dplyr", "ggplot2", "tidyr")
  missing <- required[!vapply(
    required, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing)) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "))
  }
}


quantile_or_na <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probability, names = FALSE, na.rm = TRUE))
}


safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}


save_plot <- function(plot, filename, width = 11, height = 7) {
  if (!isTRUE(RUN_FIGURES)) return(invisible(FALSE))
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = FIGURE_DPI,
    bg = "white"
  )
  invisible(TRUE)
}


validate_results <- function(results) {
  required <- c(
    "true_biological_VE", "SAR", "ICC", "coc_risk_multiplier",
    "rings", "allocation_ratio", "delay", "sd_VE", "VE_bias",
    "VE_RMSE", "valid_VE_replicates", "total_cases_prevented",
    "vaccine_doses", "cases_per_1000_doses"
  )
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    stop("Input results are missing columns: ", paste(missing, collapse = ", "))
  }
  if (!nrow(results)) stop("Input results contain no rows.")
  invisible(results)
}


select_best_within_cap <- function(league, cap) {
  candidates <- league[league$rings <= cap, , drop = FALSE]
  if (!nrow(candidates)) return(NULL)
  candidates <- candidates[order(
    -candidates$complete_grid_coverage_percent,
    -candidates$focus_region_coverage_percent,
    candidates$rings,
    candidates$p90_sd_VE,
    candidates$mean_relative_prevention_regret,
    candidates$median_vaccine_doses
  ), , drop = FALSE]
  selected <- candidates[1L, , drop = FALSE]
  selected$maximum_feasible_rings <- cap
  selected
}


select_for_coverage_target <- function(league, scope, target) {
  coverage_column <- if (scope == "Complete grid") {
    "complete_grid_coverage_percent"
  } else {
    "focus_region_coverage_percent"
  }
  qualifying <- league[
    is.finite(league[[coverage_column]]) & league[[coverage_column]] >= target,
    ,
    drop = FALSE
  ]
  if (!nrow(qualifying)) {
    return(data.frame(
      coverage_scope = scope,
      coverage_target_percent = target,
      target_met = FALSE,
      rings = NA_integer_,
      allocation_ratio = NA_real_,
      delay = NA_integer_,
      achieved_coverage_percent = NA_real_,
      complete_grid_coverage_percent = NA_real_,
      focus_region_coverage_percent = NA_real_,
      median_cases_prevented = NA_real_,
      median_vaccine_doses = NA_real_,
      mean_prevention_performance_percent = NA_real_
    ))
  }
  qualifying <- qualifying[order(
    qualifying$rings,
    -qualifying[[coverage_column]],
    qualifying$mean_relative_prevention_regret,
    qualifying$median_vaccine_doses,
    qualifying$p90_sd_VE
  ), , drop = FALSE]
  selected <- qualifying[1L, , drop = FALSE]
  data.frame(
    coverage_scope = scope,
    coverage_target_percent = target,
    target_met = TRUE,
    rings = selected$rings,
    allocation_ratio = selected$allocation_ratio,
    delay = selected$delay,
    achieved_coverage_percent = selected[[coverage_column]],
    complete_grid_coverage_percent =
      selected$complete_grid_coverage_percent,
    focus_region_coverage_percent = selected$focus_region_coverage_percent,
    median_cases_prevented = selected$median_cases_prevented,
    median_vaccine_doses = selected$median_vaccine_doses,
    mean_prevention_performance_percent =
      selected$mean_prevention_performance_percent
  )
}


select_for_retention_target <- function(league, scope, retention) {
  coverage_column <- if (scope == "Complete grid") {
    "complete_grid_coverage_percent"
  } else {
    "focus_region_coverage_percent"
  }
  maximum_coverage <- max(league[[coverage_column]], na.rm = TRUE)
  required_coverage <- retention * maximum_coverage
  qualifying <- league[
    is.finite(league[[coverage_column]]) &
      league[[coverage_column]] >= required_coverage,
    ,
    drop = FALSE
  ]
  qualifying <- qualifying[order(
    qualifying$rings,
    -qualifying[[coverage_column]],
    qualifying$mean_relative_prevention_regret,
    qualifying$median_vaccine_doses,
    qualifying$p90_sd_VE
  ), , drop = FALSE]
  selected <- qualifying[1L, , drop = FALSE]
  data.frame(
    coverage_scope = scope,
    robustness_retention_target = retention,
    maximum_observed_coverage_percent = maximum_coverage,
    required_coverage_percent = required_coverage,
    rings = selected$rings,
    allocation_ratio = selected$allocation_ratio,
    delay = selected$delay,
    achieved_coverage_percent = selected[[coverage_column]],
    achieved_retention = selected[[coverage_column]] / maximum_coverage,
    complete_grid_coverage_percent =
      selected$complete_grid_coverage_percent,
    focus_region_coverage_percent = selected$focus_region_coverage_percent,
    median_cases_prevented = selected$median_cases_prevented,
    median_vaccine_doses = selected$median_vaccine_doses,
    mean_prevention_performance_percent =
      selected$mean_prevention_performance_percent
  )
}


select_prevention_with_floor <- function(league, floor_percent) {
  qualifying <- league[
    is.finite(league$complete_grid_coverage_percent) &
      league$complete_grid_coverage_percent >= floor_percent,
    ,
    drop = FALSE
  ]
  if (!nrow(qualifying)) {
    return(data.frame(
      minimum_complete_grid_coverage_percent = floor_percent,
      rings = NA_integer_,
      allocation_ratio = NA_real_,
      delay = NA_integer_,
      complete_grid_coverage_percent = NA_real_,
      focus_region_coverage_percent = NA_real_,
      median_sd_VE = NA_real_,
      p90_sd_VE = NA_real_,
      worst_sd_VE = NA_real_,
      median_cases_prevented = NA_real_,
      median_cases_per_1000_doses = NA_real_,
      median_vaccine_doses = NA_real_,
      mean_prevention_performance_percent = NA_real_,
      mean_relative_prevention_regret = NA_real_
    ))
  }
  qualifying <- qualifying[order(
    qualifying$mean_relative_prevention_regret,
    -qualifying$complete_grid_coverage_percent,
    qualifying$median_vaccine_doses,
    qualifying$rings,
    qualifying$p90_sd_VE
  ), , drop = FALSE]
  selected <- qualifying[1L, , drop = FALSE]
  selected$minimum_complete_grid_coverage_percent <- floor_percent
  selected
}


main <- function() {
  check_dependencies()
  if (!file.exists(INPUT_RESULTS_FILE)) {
    stop(
      "Completed comprehensive results file not found: ",
      INPUT_RESULTS_FILE,
      "\nRun this script from the ring_vax_modeller_2 project folder or edit ",
      "INPUT_RESULTS_FILE."
    )
  }
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  results <- read.csv(INPUT_RESULTS_FILE, stringsAsFactors = FALSE)
  validate_results(results)

  operational_keys <- c("rings", "allocation_ratio", "delay")
  precision_truth_keys <- c("true_biological_VE", "SAR", "ICC")

  # CoC risk does not enter the direct-contact VE estimator. Collapse its
  # repeated Monte Carlo rows before counting precision ground truths.
  precision_by_truth <- results |>
    dplyr::group_by(
      rings, allocation_ratio, delay,
      true_biological_VE, SAR, ICC
    ) |>
    dplyr::summarise(
      sd_VE = mean(sd_VE, na.rm = TRUE),
      VE_bias = mean(VE_bias, na.rm = TRUE),
      VE_RMSE = mean(VE_RMSE, na.rm = TRUE),
      valid_fraction = mean(valid_VE_replicates / N_SIM, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      precision_pass = is.finite(sd_VE) & sd_VE <= PRECISION_THRESHOLD,
      bias_pass = is.finite(VE_bias) & abs(VE_bias) <= MAX_ABS_BIAS,
      validity_pass = is.finite(valid_fraction) &
        valid_fraction >= MIN_VALID_FRACTION,
      inferentially_adequate = precision_pass & bias_pass & validity_pass,
      in_focus_region =
        true_biological_VE %in% FOCUS_TRUE_VE_VALUES &
        SAR %in% FOCUS_SAR_VALUES &
        ICC %in% FOCUS_ICC_VALUES
    )

  precision_summary <- precision_by_truth |>
    dplyr::group_by(rings, allocation_ratio, delay) |>
    dplyr::summarise(
      n_complete_grid_ground_truths = dplyr::n(),
      n_complete_grid_adequate = sum(inferentially_adequate, na.rm = TRUE),
      complete_grid_coverage_percent =
        100 * mean(inferentially_adequate, na.rm = TRUE),
      n_focus_region_ground_truths = sum(in_focus_region),
      n_focus_region_adequate = sum(
        inferentially_adequate & in_focus_region, na.rm = TRUE
      ),
      focus_region_coverage_percent = if (any(in_focus_region)) {
        100 * mean(inferentially_adequate[in_focus_region], na.rm = TRUE)
      } else {
        NA_real_
      },
      precision_only_coverage_percent =
        100 * mean(precision_pass, na.rm = TRUE),
      bias_coverage_percent = 100 * mean(bias_pass, na.rm = TRUE),
      validity_coverage_percent = 100 * mean(validity_pass, na.rm = TRUE),
      median_sd_VE = stats::median(sd_VE, na.rm = TRUE),
      p90_sd_VE = quantile_or_na(sd_VE, 0.90),
      worst_sd_VE = safe_max(sd_VE),
      median_abs_bias = stats::median(abs(VE_bias), na.rm = TRUE),
      worst_abs_bias = safe_max(abs(VE_bias)),
      .groups = "drop"
    )

  prevention_scored <- results |>
    dplyr::group_by(
      true_biological_VE, SAR, ICC, coc_risk_multiplier
    ) |>
    dplyr::mutate(
      scenario_best_cases_prevented = max(
        total_cases_prevented, na.rm = TRUE
      ),
      relative_prevention_regret = dplyr::if_else(
        scenario_best_cases_prevented > 0,
        (scenario_best_cases_prevented - total_cases_prevented) /
          scenario_best_cases_prevented,
        0
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      in_focus_region =
        true_biological_VE %in% FOCUS_TRUE_VE_VALUES &
        SAR %in% FOCUS_SAR_VALUES &
        ICC %in% FOCUS_ICC_VALUES &
        coc_risk_multiplier %in% FOCUS_COC_RISK_MULTIPLIERS
    )

  prevention_summary <- prevention_scored |>
    dplyr::group_by(rings, allocation_ratio, delay) |>
    dplyr::summarise(
      n_prevention_ground_truths = dplyr::n(),
      median_cases_prevented = stats::median(
        total_cases_prevented, na.rm = TRUE
      ),
      p10_cases_prevented = quantile_or_na(total_cases_prevented, 0.10),
      p90_cases_prevented = quantile_or_na(total_cases_prevented, 0.90),
      median_cases_per_1000_doses = stats::median(
        cases_per_1000_doses, na.rm = TRUE
      ),
      median_vaccine_doses = stats::median(vaccine_doses, na.rm = TRUE),
      mean_relative_prevention_regret = mean(
        relative_prevention_regret, na.rm = TRUE
      ),
      p90_relative_prevention_regret = quantile_or_na(
        relative_prevention_regret, 0.90
      ),
      focus_median_cases_prevented = if (any(in_focus_region)) {
        stats::median(total_cases_prevented[in_focus_region], na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      mean_prevention_performance_percent =
        100 * (1 - mean_relative_prevention_regret)
    )

  league <- precision_summary |>
    dplyr::left_join(prevention_summary, by = operational_keys) |>
    dplyr::arrange(
      dplyr::desc(complete_grid_coverage_percent),
      p90_sd_VE,
      rings,
      median_vaccine_doses
    )

  budget_recommendations <- dplyr::bind_rows(lapply(
    RING_CAPS,
    function(cap) select_best_within_cap(league, cap)
  )) |>
    dplyr::select(
      maximum_feasible_rings,
      rings, allocation_ratio, delay,
      complete_grid_coverage_percent, focus_region_coverage_percent,
      median_sd_VE, p90_sd_VE, worst_sd_VE,
      median_cases_prevented, median_cases_per_1000_doses,
      median_vaccine_doses, mean_prevention_performance_percent,
      mean_relative_prevention_regret
    )

  coverage_target_recommendations <- dplyr::bind_rows(lapply(
    c("Complete grid", "Focus region"),
    function(scope) dplyr::bind_rows(lapply(
      ABSOLUTE_COVERAGE_TARGETS_PERCENT,
      function(target) select_for_coverage_target(league, scope, target)
    ))
  ))

  retention_recommendations <- dplyr::bind_rows(lapply(
    c("Complete grid", "Focus region"),
    function(scope) dplyr::bind_rows(lapply(
      ROBUSTNESS_RETENTION_TARGETS,
      function(target) select_for_retention_target(league, scope, target)
    ))
  ))

  prevention_constrained <- dplyr::bind_rows(lapply(
    PREVENTION_COVERAGE_FLOORS_PERCENT,
    function(floor_percent) select_prevention_with_floor(
      league, floor_percent
    )
  )) |>
    dplyr::select(
      minimum_complete_grid_coverage_percent,
      rings, allocation_ratio, delay,
      complete_grid_coverage_percent, focus_region_coverage_percent,
      median_sd_VE, p90_sd_VE, worst_sd_VE,
      median_cases_prevented, median_cases_per_1000_doses,
      median_vaccine_doses, mean_prevention_performance_percent,
      mean_relative_prevention_regret
    )

  best_by_ring <- league |>
    dplyr::group_by(rings) |>
    dplyr::arrange(
      dplyr::desc(complete_grid_coverage_percent),
      dplyr::desc(focus_region_coverage_percent),
      p90_sd_VE,
      mean_relative_prevention_regret,
      .by_group = TRUE
    ) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::arrange(rings) |>
    dplyr::mutate(
      additional_rings = rings - dplyr::lag(rings),
      complete_coverage_gain_points =
        complete_grid_coverage_percent -
          dplyr::lag(complete_grid_coverage_percent),
      focus_coverage_gain_points =
        focus_region_coverage_percent -
          dplyr::lag(focus_region_coverage_percent),
      complete_coverage_gain_per_100_rings =
        100 * complete_coverage_gain_points / additional_rings,
      focus_coverage_gain_per_100_rings =
        100 * focus_coverage_gain_points / additional_rings,
      additional_median_cases_prevented =
        median_cases_prevented - dplyr::lag(median_cases_prevented),
      additional_median_doses =
        median_vaccine_doses - dplyr::lag(median_vaccine_doses),
      marginal_cases_per_additional_1000_doses =
        1000 * additional_median_cases_prevented / additional_median_doses
    )

  criteria <- data.frame(
    item = c(
      "Inferential adequacy",
      "Complete-grid coverage",
      "Focus-region coverage",
      "Budget recommendation",
      "Coverage-target recommendation",
      "Robustness retention",
      "Prevention-constrained recommendation",
      "Ring-budget frontier",
      "Scenario weighting"
    ),
    definition = c(
      paste0(
        "SD(VE) <= ", PRECISION_THRESHOLD,
        "; abs(VE bias) <= ", MAX_ABS_BIAS,
        "; valid fraction >= ", MIN_VALID_FRACTION
      ),
      "Percentage of all VE/SAR/ICC grid points meeting inferential adequacy",
      "Percentage meeting adequacy within the explicitly declared focus values",
      "Most robust design that does not exceed the stated ring cap",
      "Smallest-ring design meeting the stated absolute coverage target",
      "Smallest-ring design retaining the stated fraction of maximum observed coverage",
      "Lowest mean prevention regret subject to the stated complete-grid coverage floor",
      "Most inferentially robust design available at each configured maximum ring count",
      "Grid points are equally weighted sensitivity cases, not estimated probabilities"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    league,
    file.path(OUTPUT_DIR, "tradeoff_design_league_table.csv"),
    row.names = FALSE
  )
  write.csv(
    budget_recommendations,
    file.path(OUTPUT_DIR, "recommendations_by_ring_cap.csv"),
    row.names = FALSE
  )
  write.csv(
    coverage_target_recommendations,
    file.path(OUTPUT_DIR, "recommendations_by_coverage_target.csv"),
    row.names = FALSE
  )
  write.csv(
    retention_recommendations,
    file.path(OUTPUT_DIR, "recommendations_by_robustness_retention.csv"),
    row.names = FALSE
  )
  write.csv(
    prevention_constrained,
    file.path(OUTPUT_DIR, "recommendations_by_prevention_constraint.csv"),
    row.names = FALSE
  )
  write.csv(
    best_by_ring,
    file.path(OUTPUT_DIR, "marginal_robustness_by_ring.csv"),
    row.names = FALSE
  )
  write.csv(
    budget_recommendations,
    file.path(OUTPUT_DIR, "ring_budget_frontier.csv"),
    row.names = FALSE
  )
  write.csv(
    criteria,
    file.path(OUTPUT_DIR, "tradeoff_analysis_definitions.csv"),
    row.names = FALSE
  )

  if (isTRUE(RUN_FIGURES)) {
    coverage_plot_data <- dplyr::bind_rows(
      best_by_ring |>
        dplyr::transmute(
          rings,
          coverage_scope = "Complete uncertainty grid",
          coverage_percent = complete_grid_coverage_percent
        ),
      best_by_ring |>
        dplyr::transmute(
          rings,
          coverage_scope = "Declared focus region",
          coverage_percent = focus_region_coverage_percent
        )
    )
    p_coverage <- ggplot2::ggplot(
      coverage_plot_data,
      ggplot2::aes(x = rings, y = coverage_percent)
    ) +
      ggplot2::geom_line(linewidth = 0.8, colour = "#2878B5") +
      ggplot2::geom_point(size = 2.5, colour = "#2878B5") +
      ggplot2::geom_text(
        ggplot2::aes(label = paste0(round(coverage_percent, 1), "%")),
        nudge_y = 2,
        size = 3.2
      ) +
      ggplot2::facet_wrap(~coverage_scope, ncol = 1) +
      ggplot2::scale_x_continuous(breaks = best_by_ring$rings) +
      ggplot2::coord_cartesian(ylim = c(0, 100)) +
      ggplot2::labs(
        x = "Maximum number of rings",
        y = "Ground truths meeting all inferential criteria (%)",
        title = "Robustness achieved at each ring budget",
        subtitle = "Each point is the most robust design available within that ring count"
      ) +
      ggplot2::theme_bw(base_size = 11)
    save_plot(
      p_coverage,
      file.path(OUTPUT_DIR, "figure_tradeoff_coverage_by_ring_budget.png"),
      height = 8
    )

    marginal_plot_data <- best_by_ring |>
      dplyr::filter(is.finite(complete_coverage_gain_per_100_rings)) |>
      dplyr::select(
        rings,
        `Complete uncertainty grid` = complete_coverage_gain_per_100_rings,
        `Declared focus region` = focus_coverage_gain_per_100_rings
      ) |>
      tidyr::pivot_longer(
        cols = -rings,
        names_to = "coverage_scope",
        values_to = "coverage_gain_per_100_rings"
      )
    p_marginal <- ggplot2::ggplot(
      marginal_plot_data,
      ggplot2::aes(
        x = factor(rings),
        y = coverage_gain_per_100_rings,
        fill = coverage_scope
      )
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_text(
        ggplot2::aes(
          label = round(coverage_gain_per_100_rings, 1)
        ),
        position = ggplot2::position_dodge(width = 0.9),
        vjust = -0.4,
        size = 3.2
      ) +
      ggplot2::labs(
        x = "New ring budget",
        y = "Coverage gain per additional 100 rings (percentage points)",
        fill = "Coverage scope",
        title = "Diminishing robustness returns from additional rings",
        subtitle = "Bars compare each ring count with the next smaller simulated option"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    save_plot(
      p_marginal,
      file.path(OUTPUT_DIR, "figure_tradeoff_marginal_robustness.png")
    )

    p_frontier <- ggplot2::ggplot(
      league,
      ggplot2::aes(
        x = rings,
        y = complete_grid_coverage_percent,
        colour = mean_prevention_performance_percent,
        shape = factor(delay),
        size = median_vaccine_doses
      )
    ) +
      ggplot2::geom_point(alpha = 0.45) +
      ggplot2::geom_point(
        data = budget_recommendations,
        shape = 21,
        fill = "white",
        colour = "black",
        stroke = 1,
        alpha = 1
      ) +
      ggplot2::geom_point(
        data = prevention_constrained,
        shape = 23,
        fill = "white",
        colour = "black",
        stroke = 1,
        alpha = 1
      ) +
      ggplot2::scale_colour_viridis_c(
        option = "C",
        name = "Mean prevention\nperformance (%)"
      ) +
      ggplot2::scale_x_continuous(breaks = sort(unique(league$rings))) +
      ggplot2::labs(
        x = "Number of rings",
        y = "Complete-grid inferential coverage (%)",
        shape = "Delay (days)",
        size = "Median vaccine doses",
        title = "Resource–robustness–prevention decision landscape",
        subtitle = paste0(
          "Circles: best design at each ring cap; diamonds: prevention choice ",
          "at each coverage floor"
        )
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    save_plot(
      p_frontier,
      file.path(OUTPUT_DIR, "figure_tradeoff_resource_frontier.png")
    )

    p_prevention <- ggplot2::ggplot(
      league,
      ggplot2::aes(
        x = complete_grid_coverage_percent,
        y = mean_prevention_performance_percent,
        colour = factor(allocation_ratio),
        shape = factor(delay),
        size = rings
      )
    ) +
      ggplot2::geom_point(alpha = 0.55) +
      ggplot2::geom_point(
        data = prevention_constrained,
        shape = 21,
        fill = "white",
        colour = "black",
        stroke = 1,
        alpha = 1
      ) +
      ggplot2::labs(
        x = "Complete-grid inferential coverage (%)",
        y = "Mean prevention performance versus scenario optimum (%)",
        colour = "Immediate:delayed allocation",
        shape = "Delay (days)",
        size = "Rings",
        title = "Prevention choices under explicit robustness constraints",
        subtitle = "Outlined designs are selected at the configured coverage floors"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    save_plot(
      p_prevention,
      file.path(OUTPUT_DIR, "figure_tradeoff_prevention_constraints.png")
    )
  }

  message("Trade-off analysis completed: ", normalizePath(
    OUTPUT_DIR, mustWork = FALSE
  ))
  invisible(list(
    league = league,
    budget_recommendations = budget_recommendations,
    coverage_target_recommendations = coverage_target_recommendations,
    retention_recommendations = retention_recommendations,
    prevention_constrained = prevention_constrained,
    best_by_ring = best_by_ring,
    ring_budget_frontier = budget_recommendations
  ))
}


if (sys.nframe() == 0L) main()
