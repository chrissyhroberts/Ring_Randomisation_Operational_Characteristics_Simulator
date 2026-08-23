# ============================================================
# Ring vaccination design analysis
# Pareto extraction + figures + decision tables
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# Load results
# ------------------------------------------------------------

results <- read_csv("v1,18.1_final_output.csv")


# ------------------------------------------------------------
# Pareto frontier
#
# Objectives:
#   maximise total_cases_prevented
#   minimise sd_VE
# ------------------------------------------------------------

pareto_flag <- function(df){

  dominated <- rep(FALSE, nrow(df))

  for(i in seq_len(nrow(df))){

    others <- setdiff(seq_len(nrow(df)), i)

    dominated[i] <- any(

      df$total_cases_prevented[others] >=
        df$total_cases_prevented[i] &

      df$sd_VE[others] <=
        df$sd_VE[i] &

      (
        df$total_cases_prevented[others] >
          df$total_cases_prevented[i] |

        df$sd_VE[others] <
          df$sd_VE[i]
      )

    )

  }

  df$pareto_optimal <- !dominated

  df

}


# ------------------------------------------------------------
# Pareto results
# ------------------------------------------------------------

pareto_results <-

  results %>%
  group_by(
    true_biological_VE,
    SAR
  ) %>%
  group_modify(
    ~pareto_flag(.x)
  ) %>%
  ungroup()


pareto_only <-
  pareto_results %>%
  filter(pareto_optimal)


stopifnot(nrow(pareto_only) > 0)


# ------------------------------------------------------------
# Pareto summary
# ------------------------------------------------------------

pareto_summary <-

  pareto_only %>%
  group_by(
    true_biological_VE,
    SAR
  ) %>%
  summarise(

    n_pareto_designs = n(),

    allocation_ratios =
      paste(sort(unique(allocation_ratio)), collapse=", "),

    rings_range =
      paste(min(rings), max(rings), sep="-"),

    cases_prevented_range =
      paste(
        round(min(total_cases_prevented)),
        round(max(total_cases_prevented)),
        sep="-"
      ),

    sd_VE_range =
      paste(
        round(min(sd_VE),3),
        round(max(sd_VE),3),
        sep="-"
      ),

    .groups="drop"

  )

write_csv(pareto_summary,
          "pareto_frontier_summary.csv")

write_csv(pareto_only,
          "pareto_frontier_designs.csv")


# ------------------------------------------------------------
# Figure 1: Pareto frontier
# ------------------------------------------------------------

p1 <-

ggplot(
  pareto_results,
  aes(
    x=-sd_VE,
    y=total_cases_prevented,
    colour=factor(allocation_ratio),
    size=rings
  )
) +

geom_point(alpha=0.4) +

geom_point(
  data=pareto_only,
  shape=21,
  fill="white",
  colour="black",
  size=4
) +

facet_grid(
  true_biological_VE ~ SAR,
  labeller=label_both
) +

labs(
  x="Inferential precision (higher = lower SD(VE))",
  y="Total cases prevented",
  colour="Allocation ratio",
  size="Rings",
  title="Pareto frontier: ethical benefit versus inferential precision"
) +

theme_bw()

ggsave(
  "figure_pareto_frontier.png",
  p1,
  width=12,
  height=8,
  dpi=300
)



# ------------------------------------------------------------
# Figure 2: allocation trade-off
# ------------------------------------------------------------

allocation_gain <-

results %>%
group_by(
  true_biological_VE,
  SAR,
  rings
) %>%
mutate(

  baseline_cases =
    total_cases_prevented[allocation_ratio==1][1],

  baseline_sd =
    sd_VE[allocation_ratio==1][1]

) %>%
filter(allocation_ratio>1) %>%
mutate(

  additional_cases =
    total_cases_prevented-baseline_cases,

  precision_cost =
    sd_VE-baseline_sd

)


p2 <-

ggplot(
  allocation_gain,
  aes(
    x=precision_cost,
    y=additional_cases,
    colour=factor(allocation_ratio),
    size=rings
  )
) +

geom_hline(yintercept=0) +
geom_vline(xintercept=0) +
geom_point(alpha=0.6) +

facet_grid(
  true_biological_VE ~ SAR,
  labeller=label_both
) +

labs(
  x="Precision cost (increase in SD(VE))",
  y="Additional cases prevented versus 1:1",
  colour="Allocation ratio",
  size="Rings",
  title="Marginal ethical benefit of unequal allocation"
) +

theme_bw()

ggsave(
  "figure_allocation_tradeoff.png",
  p2,
  width=12,
  height=8,
  dpi=300
)



# ------------------------------------------------------------
# Figure 3: allocation benefit
# ------------------------------------------------------------

p3 <-

results %>%

ggplot(
  aes(
    factor(allocation_ratio),
    total_cases_prevented,
    fill=factor(allocation_ratio)
  )
) +

geom_boxplot() +

facet_grid(
  true_biological_VE ~ SAR,
  labeller=label_both
) +

labs(
  x="Allocation ratio",
  y="Total cases prevented",
  title="Public health benefit of unequal allocation"
) +

theme_bw()

ggsave(
  "figure_allocation_cases.png",
  p3,
  width=12,
  height=8,
  dpi=300
)



# ------------------------------------------------------------
# Figure 4: delay sensitivity
# ------------------------------------------------------------

p4 <-

results %>%

ggplot(
  aes(
    delay,
    total_cases_prevented,
    colour=factor(allocation_ratio),
    group=interaction(allocation_ratio,rings)
  )
) +

geom_line(alpha=0.5) +

facet_grid(
  true_biological_VE ~ SAR,
  labeller=label_both
) +

labs(
  x="Vaccination delay",
  y="Cases prevented",
  colour="Allocation ratio",
  title="Impact of vaccination delay"
) +

theme_bw()

ggsave(
  "figure_delay_effect.png",
  p4,
  width=12,
  height=8,
  dpi=300
)



# ------------------------------------------------------------
# Decision table
# ------------------------------------------------------------

decision_table <-

results %>%

group_by(
  true_biological_VE,
  SAR
) %>%

mutate(
  best_precision=min(sd_VE),
  best_cases=max(total_cases_prevented)
) %>%

ungroup()


inference_priority <-

decision_table %>%
filter(sd_VE <= best_precision*1.05) %>%
group_by(true_biological_VE,SAR) %>%
slice_max(total_cases_prevented,n=1) %>%
ungroup() %>%
mutate(priority="Inference")


prevention_priority <-

decision_table %>%
filter(total_cases_prevented >= best_cases*0.95) %>%
group_by(true_biological_VE,SAR) %>%
slice_min(sd_VE,n=1) %>%
ungroup() %>%
mutate(priority="Prevention")


recommendations <-

bind_rows(
  inference_priority,
  prevention_priority
) %>%

select(
  priority,
  true_biological_VE,
  SAR,
  rings,
  allocation_ratio,
  delay,
  total_cases_prevented,
  sd_VE
)

write_csv(
  recommendations,
  "allocation_recommendation_table.csv"
)


print(pareto_summary)
print(recommendations)


# ============================================================
# VE precision / rings required analysis
# ============================================================

# Define precision thresholds
precision_thresholds <- c(
  0.02,
  0.05,
  0.10
)

rings_required <-

  purrr::map_dfr(
    precision_thresholds,
    function(threshold){

      results %>%
        filter(delay == min(delay)) %>%
        group_by(
          true_biological_VE,
          SAR,
          allocation_ratio
        ) %>%
        filter(sd_VE <= threshold) %>%
        summarise(
          precision_threshold = threshold,
          minimum_rings = min(rings),
          achieved_sd_VE = min(sd_VE),
          .groups = "drop"
        )

    }
  )

write_csv(
  rings_required,
  "rings_required_for_precision.csv"
)


# Plot SD(VE) as a function of number of rings

p_precision <-

results %>%
filter(delay == min(delay)) %>%
ggplot(
  aes(
    x=rings,
    y=sd_VE,
    colour=factor(allocation_ratio)
  )
) +

geom_line() +
geom_point() +

geom_hline(
  yintercept=c(0.02,0.05,0.10),
  linetype="dashed"
) +

facet_grid(
  true_biological_VE ~ SAR,
  labeller=label_both
) +

labs(
  x="Number of rings",
  y="SD of VE estimate",
  colour="Allocation ratio",
  title="Precision of VE estimation by number of rings"
) +

theme_bw()


ggsave(
  "figure_rings_vs_precision.png",
  p_precision,
  width=12,
  height=8,
  dpi=300
)


# ============================================================
# Event burden diagnostics
# ============================================================

event_summary_columns <- intersect(
  c(
    "cases",
    "total_cases",
    "cases_observed",
    "direct_cases",
    "coc_cases",
    "direct_cases_prevented",
    "coc_cases_prevented"
  ),
  names(results)
)

write_csv(
  results %>%
    select(all_of(event_summary_columns)),
  "available_event_outputs.csv"
)


# ============================================================
# Marginal benefit per additional ring
# ============================================================

ring_efficiency <-

results %>%
filter(delay == min(delay)) %>%
arrange(
  true_biological_VE,
  SAR,
  allocation_ratio,
  rings
) %>%
group_by(
  true_biological_VE,
  SAR,
  allocation_ratio
) %>%
mutate(
  additional_cases_per_ring =
    (total_cases_prevented - lag(total_cases_prevented)) /
    (rings - lag(rings))
) %>%
ungroup()


write_csv(
  ring_efficiency,
  "marginal_cases_per_additional_ring.csv"
)


p_efficiency <-

ring_efficiency %>%
filter(!is.na(additional_cases_per_ring)) %>%
ggplot(
  aes(
    x=rings,
    y=additional_cases_per_ring,
    colour=factor(allocation_ratio)
  )
) +

geom_line() +
geom_point() +

facet_grid(
  true_biological_VE ~ SAR,
  labeller=label_both
) +

labs(
  x="Number of rings",
  y="Additional cases prevented per additional ring",
  colour="Allocation ratio",
  title="Diminishing public health returns with increasing ring numbers"
) +

theme_bw()


ggsave(
  "figure_marginal_ring_benefit.png",
  p_efficiency,
  width=12,
  height=8,
  dpi=300
)


# ============================================================
# End
# ============================================================


# ============================================================
# Optional sensitivity analysis: latent between-ring SAR heterogeneity
# ============================================================

# This section is intentionally optional.
#
# The core model assumes SAR is identical across rings.
# If future simulation outputs include:
#
#   SAR_CV
#
# representing latent between-ring heterogeneity, this section
# analyses how many rings are required for VE precision under
# increasing heterogeneity.
#
# Observed ring-level attack-rate CV should not be used directly:
# sparse event counts inflate empirical CV.

if("SAR_CV" %in% names(results)){

  write_csv(
    results %>%
      group_by(
        true_biological_VE,
        SAR,
        SAR_CV,
        allocation_ratio
      ) %>%
      summarise(
        minimum_rings_sd05 =
          ifelse(any(sd_VE <= 0.05),
                 min(rings[sd_VE <= 0.05]),
                 NA),
        minimum_rings_sd02 =
          ifelse(any(sd_VE <= 0.02),
                 min(rings[sd_VE <= 0.02]),
                 NA),
        best_sd_VE = min(sd_VE),
        .groups="drop"
      ),
    "rings_required_by_SAR_heterogeneity.csv"
  )


  p_sar_cv <-

    results %>%
    filter(delay == min(delay)) %>%
    group_by(
      true_biological_VE,
      SAR,
      SAR_CV,
      allocation_ratio
    ) %>%
    summarise(
      minimum_rings =
        ifelse(any(sd_VE <= 0.05),
               min(rings[sd_VE <= 0.05]),
               NA),
      .groups="drop"
    ) %>%

    ggplot(
      aes(
        x=SAR_CV,
        y=minimum_rings,
        colour=factor(allocation_ratio)
      )
    ) +

    geom_point(size=3) +
    geom_line(aes(group=allocation_ratio)) +

    facet_grid(
      true_biological_VE ~ SAR,
      labeller=label_both
    ) +

    labs(
      x="Latent between-ring SAR heterogeneity (CV)",
      y="Minimum rings required for SD(VE) <= 0.05",
      colour="Allocation ratio",
      title="Impact of latent SAR heterogeneity on required ring numbers"
    ) +

    theme_bw()


  ggsave(
    "figure_SAR_heterogeneity_rings_required.png",
    p_sar_cv,
    width=12,
    height=8,
    dpi=300
  )


  p_sar_cv_cost <-

    results %>%
    filter(delay == min(delay)) %>%
    group_by(
      true_biological_VE,
      SAR,
      SAR_CV
    ) %>%
    summarise(
      best_precision=min(sd_VE),
      .groups="drop"
    ) %>%

    ggplot(
      aes(
        x=SAR_CV,
        y=best_precision
      )
    ) +

    geom_point() +
    geom_line() +

    facet_grid(
      true_biological_VE ~ SAR,
      labeller=label_both
    ) +

    labs(
      x="Latent between-ring SAR heterogeneity (CV)",
      y="Best achievable SD(VE)",
      title="Precision loss under increasing SAR heterogeneity"
    ) +

    theme_bw()


  ggsave(
    "figure_SAR_heterogeneity_precision.png",
    p_sar_cv_cost,
    width=12,
    height=8,
    dpi=300
  )


} else {

  message(
    "No SAR_CV column found. Skipping latent SAR heterogeneity analysis."
  )

}


# ============================================================
# End
# ============================================================
