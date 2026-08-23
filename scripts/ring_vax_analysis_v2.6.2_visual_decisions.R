# ============================================================
# Ring-vaccination visual decision analysis v2.6.2
# ============================================================
#
# Analysis only: this script never runs or resumes simulations.
# It reads the completed comprehensive operating-characteristics table and
# creates colleague-facing comparisons of:
#   1. shortlisted designs across SAR and ICC;
#   2. CoC-risk sensitivity for those shortlisted designs; and
#   3. all operational designs in a common decision landscape.
#
# Grid points are equally weighted sensitivity cases. Coverage percentages
# are not estimated probabilities that a design will succeed in practice.
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
  "analysis_v2.6.2_visual_decisions"
)

N_SIM <- 2000L
FIGURE_DPI <- 300L
RUN_FIGURES <- TRUE

# Inferential adequacy criteria
PRECISION_THRESHOLD <- 0.15
MAX_ABS_BIAS <- 0.05
MIN_VALID_FRACTION <- 0.95

# Declared focus region. This is a sensitivity subset, not a probability model.
FOCUS_TRUE_VE_VALUES <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
FOCUS_SAR_VALUES <- c(0.02, 0.04, 0.06)
FOCUS_ICC_VALUES <- c(0.05, 0.10, 0.15)
FOCUS_COC_RISK_MULTIPLIERS <- c(0.05, 0.10, 0.25)

# Designs to compare side by side. Add or remove rows as required.
SHORTLIST_DESIGNS <- data.frame(
  rings = c(1000L, 1200L, 1200L, 1500L, 1500L),
  allocation_ratio = c(1, 1, 1, 1, 1),
  delay = c(21L, 21L, 14L, 21L, 7L),
  decision_label = c(
    "Lean inference",
    "Preferred robust",
    "Balanced",
    "Maximum robustness",
    "Prevention-forward"
  ),
  stringsAsFactors = FALSE
)

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


save_plot <- function(plot, filename, width = 12, height = 8) {
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
    "valid_VE_replicates", "total_cases_prevented", "vaccine_doses",
    "cases_per_1000_doses"
  )
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    stop("Input results are missing columns: ", paste(missing, collapse = ", "))
  }
  if (!nrow(results)) stop("Input results contain no rows.")
  invisible(results)
}


design_key <- function(rings, allocation_ratio, delay) {
  paste(rings, allocation_ratio, delay, sep = "|")
}


format_design_label <- function(label, rings, allocation_ratio, delay) {
  paste0(
    label, "\n",
    rings, " rings; ", allocation_ratio, ":1; ", delay, " days"
  )
}


make_sar_icc_summary <- function(data, scope_label, focus_only = FALSE) {
  selected <- data
  if (isTRUE(focus_only)) {
    selected <- selected |>
      dplyr::filter(
        true_biological_VE %in% FOCUS_TRUE_VE_VALUES,
        SAR %in% FOCUS_SAR_VALUES,
        ICC %in% FOCUS_ICC_VALUES
      )
  }
  selected |>
    dplyr::group_by(
      design_order, design_label, rings, allocation_ratio, delay, SAR, ICC
    ) |>
    dplyr::summarise(
      n_true_VE_values = dplyr::n(),
      n_adequate = sum(inferentially_adequate, na.rm = TRUE),
      adequate_VE_percent = 100 * mean(
        inferentially_adequate, na.rm = TRUE
      ),
      median_sd_VE = stats::median(sd_VE, na.rm = TRUE),
      maximum_sd_VE = max(sd_VE, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(coverage_scope = scope_label)
}


plot_sar_icc_heatmap <- function(data, title, subtitle) {
  plot_data <- data |>
    dplyr::mutate(
      SAR_label = factor(
        format(SAR, trim = TRUE, scientific = FALSE),
        levels = format(sort(unique(SAR)), trim = TRUE, scientific = FALSE)
      ),
      ICC_label = factor(
        format(ICC, trim = TRUE, scientific = FALSE),
        levels = rev(format(
          sort(unique(ICC)), trim = TRUE, scientific = FALSE
        ))
      ),
      cell_label = paste0(round(adequate_VE_percent), "%"),
      text_colour = ifelse(adequate_VE_percent < 45, "white", "black"),
      design_label = factor(
        design_label,
        levels = unique(design_label[order(design_order)])
      )
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = SAR_label, y = ICC_label, fill = adequate_VE_percent)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = cell_label, colour = text_colour),
      size = 3
    ) +
    ggplot2::facet_wrap(~design_label, ncol = 3) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      name = "True VE values\nadequate (%)"
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::labs(
      x = "Secondary attack rate",
      y = "Intracluster correlation coefficient",
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 9),
      legend.position = "right"
    )
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

  SHORTLIST_DESIGNS$design_order <- seq_len(nrow(SHORTLIST_DESIGNS))
  SHORTLIST_DESIGNS$design_label <- format_design_label(
    SHORTLIST_DESIGNS$decision_label,
    SHORTLIST_DESIGNS$rings,
    SHORTLIST_DESIGNS$allocation_ratio,
    SHORTLIST_DESIGNS$delay
  )

  available_keys <- unique(design_key(
    results$rings, results$allocation_ratio, results$delay
  ))
  requested_keys <- design_key(
    SHORTLIST_DESIGNS$rings,
    SHORTLIST_DESIGNS$allocation_ratio,
    SHORTLIST_DESIGNS$delay
  )
  missing_shortlist <- SHORTLIST_DESIGNS[!requested_keys %in% available_keys, ]
  if (nrow(missing_shortlist)) {
    stop(
      "These shortlisted designs were not found in the results: ",
      paste(missing_shortlist$design_label, collapse = "; ")
    )
  }

  # CoC risk does not enter the direct-contact VE estimator. Collapse repeated
  # CoC rows before evaluating VE/SAR/ICC inferential ground truths.
  precision_by_truth <- results |>
    dplyr::group_by(
      rings, allocation_ratio, delay,
      true_biological_VE, SAR, ICC
    ) |>
    dplyr::summarise(
      sd_VE = mean(sd_VE, na.rm = TRUE),
      VE_bias = mean(VE_bias, na.rm = TRUE),
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
      complete_grid_coverage_percent = 100 * mean(
        inferentially_adequate, na.rm = TRUE
      ),
      focus_region_coverage_percent = if (any(in_focus_region)) {
        100 * mean(inferentially_adequate[in_focus_region], na.rm = TRUE)
      } else {
        NA_real_
      },
      median_sd_VE = stats::median(sd_VE, na.rm = TRUE),
      p90_sd_VE = quantile_or_na(sd_VE, 0.90),
      .groups = "drop"
    )

  # Prevention is allowed to vary with CoC risk. Scenario-relative performance
  # prevents high-transmission grid points from dominating simply by scale.
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
    dplyr::ungroup()

  prevention_summary <- prevention_scored |>
    dplyr::group_by(rings, allocation_ratio, delay) |>
    dplyr::summarise(
      median_cases_prevented = stats::median(
        total_cases_prevented, na.rm = TRUE
      ),
      median_vaccine_doses = stats::median(vaccine_doses, na.rm = TRUE),
      median_cases_per_1000_doses = stats::median(
        cases_per_1000_doses, na.rm = TRUE
      ),
      mean_prevention_performance_percent = 100 * (
        1 - mean(relative_prevention_regret, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  all_models <- precision_summary |>
    dplyr::left_join(
      prevention_summary,
      by = c("rings", "allocation_ratio", "delay")
    ) |>
    dplyr::mutate(
      design_key = design_key(rings, allocation_ratio, delay),
      shortlisted = design_key %in% requested_keys,
      ring_panel = factor(
        paste0(rings, " rings"),
        levels = paste0(sort(unique(rings)), " rings")
      )
    ) |>
    dplyr::arrange(rings, allocation_ratio, delay)

  shortlist_precision <- precision_by_truth |>
    dplyr::inner_join(
      SHORTLIST_DESIGNS |>
        dplyr::select(
          rings, allocation_ratio, delay,
          design_order, decision_label, design_label
        ),
      by = c("rings", "allocation_ratio", "delay")
    )

  sar_icc_complete <- make_sar_icc_summary(
    shortlist_precision,
    scope_label = "Complete uncertainty grid",
    focus_only = FALSE
  )
  sar_icc_focus <- make_sar_icc_summary(
    shortlist_precision,
    scope_label = "Declared focus region",
    focus_only = TRUE
  )
  sar_icc_summary <- dplyr::bind_rows(sar_icc_complete, sar_icc_focus) |>
    dplyr::arrange(coverage_scope, design_order, ICC, SAR)

  shortlist_summary <- all_models |>
    dplyr::filter(shortlisted) |>
    dplyr::left_join(
      SHORTLIST_DESIGNS |>
        dplyr::select(
          rings, allocation_ratio, delay,
          design_order, decision_label, design_label
        ),
      by = c("rings", "allocation_ratio", "delay")
    ) |>
    dplyr::arrange(design_order)

  coc_sensitivity <- prevention_scored |>
    dplyr::inner_join(
      SHORTLIST_DESIGNS |>
        dplyr::select(
          rings, allocation_ratio, delay,
          design_order, decision_label, design_label
        ),
      by = c("rings", "allocation_ratio", "delay")
    ) |>
    dplyr::group_by(
      design_order, design_label, rings, allocation_ratio, delay,
      coc_risk_multiplier
    ) |>
    dplyr::summarise(
      median_cases_prevented = stats::median(
        total_cases_prevented, na.rm = TRUE
      ),
      p10_cases_prevented = quantile_or_na(total_cases_prevented, 0.10),
      p90_cases_prevented = quantile_or_na(total_cases_prevented, 0.90),
      median_cases_per_1000_doses = stats::median(
        cases_per_1000_doses, na.rm = TRUE
      ),
      median_vaccine_doses = stats::median(vaccine_doses, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(design_order, coc_risk_multiplier)

  definitions <- data.frame(
    item = c(
      "Inferential adequacy",
      "SAR-ICC heatmap cell",
      "Complete heatmap",
      "Focus heatmap",
      "All-model overview",
      "All-model atlas cell",
      "CoC sensitivity",
      "Scenario weighting"
    ),
    definition = c(
      paste0(
        "SD(VE) <= ", PRECISION_THRESHOLD,
        "; abs(VE bias) <= ", MAX_ABS_BIAS,
        "; valid fraction >= ", MIN_VALID_FRACTION
      ),
      "Percentage of simulated true VE values meeting all inferential criteria at the stated SAR and ICC",
      "Uses every simulated VE, SAR and ICC value",
      "Uses only the explicitly declared focus VE, SAR and ICC values",
      "Every point is one rings/allocation/delay design; black outlines identify the shortlist",
      "Top percentage is complete-grid adequacy; bottom percentage is prevention performance",
      "Cases prevented are summarised across VE, SAR and ICC at each CoC-risk multiplier",
      "Grid points are equally weighted sensitivity cases, not estimated probabilities"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    shortlist_summary,
    file.path(OUTPUT_DIR, "shortlisted_design_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    sar_icc_summary,
    file.path(OUTPUT_DIR, "shortlist_SAR_ICC_adequacy.csv"),
    row.names = FALSE
  )
  write.csv(
    coc_sensitivity,
    file.path(OUTPUT_DIR, "shortlist_CoC_sensitivity.csv"),
    row.names = FALSE
  )
  write.csv(
    all_models,
    file.path(OUTPUT_DIR, "all_models_overview.csv"),
    row.names = FALSE
  )
  write.csv(
    definitions,
    file.path(OUTPUT_DIR, "visual_analysis_definitions.csv"),
    row.names = FALSE
  )

  if (isTRUE(RUN_FIGURES)) {
    p_complete <- plot_sar_icc_heatmap(
      sar_icc_complete,
      title = "Where shortlisted designs remain inferentially adequate",
      subtitle = paste0(
        "Each cell is the percentage of all simulated true VE values passing ",
        "the precision, bias and validity criteria"
      )
    )
    save_plot(
      p_complete,
      file.path(
        OUTPUT_DIR,
        "figure_shortlist_SAR_ICC_heatmap_complete.png"
      ),
      width = 13,
      height = 8
    )

    p_focus <- plot_sar_icc_heatmap(
      sar_icc_focus,
      title = "SAR-ICC adequacy within the declared focus region",
      subtitle = paste0(
        "Cells average only over the declared focus VE values; grid points ",
        "are sensitivity cases, not probabilities"
      )
    )
    save_plot(
      p_focus,
      file.path(
        OUTPUT_DIR,
        "figure_shortlist_SAR_ICC_heatmap_focus.png"
      ),
      width = 13,
      height = 8
    )

    point_shapes <- c("7" = 21, "14" = 24, "21" = 22)
    p_all_models <- ggplot2::ggplot(
      all_models,
      ggplot2::aes(
        x = complete_grid_coverage_percent,
        y = mean_prevention_performance_percent,
        fill = factor(allocation_ratio),
        shape = factor(delay)
      )
    ) +
      ggplot2::geom_point(
        size = 3.2,
        alpha = 0.72,
        colour = "grey35",
        stroke = 0.35
      ) +
      ggplot2::geom_point(
        data = all_models |>
          dplyr::filter(shortlisted),
        size = 5.2,
        alpha = 1,
        colour = "black",
        stroke = 1.1,
        show.legend = FALSE
      ) +
      ggplot2::facet_wrap(~ring_panel, ncol = 3) +
      ggplot2::scale_shape_manual(values = point_shapes) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
      ggplot2::labs(
        x = "Complete-grid inferential adequacy (%)",
        y = "Prevention performance versus scenario optimum (%)",
        fill = "Immediate:delayed allocation",
        shape = "Delay (days)",
        title = "All operational designs at a glance",
        subtitle = paste0(
          "Each point is one design; black outlines identify shortlisted ",
          "designs"
        )
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    save_plot(
      p_all_models,
      file.path(OUTPUT_DIR, "figure_all_models_tradeoff.png"),
      width = 12,
      height = 10
    )

    atlas_data <- all_models |>
      dplyr::mutate(
        delay_factor = factor(delay, levels = sort(unique(delay))),
        allocation_factor = factor(
          allocation_ratio,
          levels = rev(sort(unique(allocation_ratio))),
          labels = paste0(rev(sort(unique(allocation_ratio))), ":1")
        ),
        atlas_label = paste0(
          round(complete_grid_coverage_percent), "%\n",
          round(mean_prevention_performance_percent), "%"
        ),
        text_colour = ifelse(
          complete_grid_coverage_percent < 45, "white", "black"
        )
      )
    p_atlas <- ggplot2::ggplot(
      atlas_data,
      ggplot2::aes(
        x = delay_factor,
        y = allocation_factor,
        fill = complete_grid_coverage_percent
      )
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
      ggplot2::geom_tile(
        data = atlas_data |>
          dplyr::filter(shortlisted),
        fill = NA,
        colour = "black",
        linewidth = 1.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = atlas_label, colour = text_colour),
        size = 3,
        lineheight = 0.9
      ) +
      ggplot2::facet_wrap(~ring_panel, ncol = 3) +
      ggplot2::scale_fill_viridis_c(
        option = "C",
        limits = c(0, 100),
        name = "Complete-grid\nadequacy (%)"
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::labs(
        x = "Delay (days)",
        y = "Immediate:delayed allocation",
        title = "Operational design atlas",
        subtitle = paste0(
          "Cell labels: complete-grid adequacy (top) and prevention ",
          "performance (bottom); black borders mark the shortlist"
        )
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "right"
      )
    save_plot(
      p_atlas,
      file.path(OUTPUT_DIR, "figure_all_models_design_atlas.png"),
      width = 12,
      height = 10
    )

    coc_plot_data <- coc_sensitivity |>
      dplyr::mutate(
        design_label = factor(
          design_label,
          levels = unique(design_label[order(design_order)])
        )
      )
    p_coc <- ggplot2::ggplot(
      coc_plot_data,
      ggplot2::aes(
        x = coc_risk_multiplier,
        y = median_cases_prevented,
        group = design_label
      )
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = p10_cases_prevented,
          ymax = p90_cases_prevented
        ),
        fill = "#76C7C0",
        alpha = 0.22,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.8, colour = "#178F8A") +
      ggplot2::geom_point(size = 2.2, colour = "#178F8A") +
      ggplot2::facet_wrap(~design_label, ncol = 3) +
      ggplot2::scale_x_continuous(
        breaks = sort(unique(coc_plot_data$coc_risk_multiplier))
      ) +
      ggplot2::labs(
        x = "CoC risk relative to a direct contact",
        y = "Cases prevented",
        title = "How CoC-risk uncertainty changes prevention",
        subtitle = paste0(
          "Lines show medians; bands span the 10th to 90th percentiles ",
          "across VE, SAR and ICC"
        )
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    save_plot(
      p_coc,
      file.path(OUTPUT_DIR, "figure_shortlist_CoC_sensitivity.png"),
      width = 13,
      height = 8
    )
  }

  message("Visual decision analysis completed: ", normalizePath(
    OUTPUT_DIR, mustWork = FALSE
  ))
  invisible(list(
    shortlisted_designs = shortlist_summary,
    sar_icc_adequacy = sar_icc_summary,
    coc_sensitivity = coc_sensitivity,
    all_models = all_models
  ))
}


if (sys.nframe() == 0L) main()
