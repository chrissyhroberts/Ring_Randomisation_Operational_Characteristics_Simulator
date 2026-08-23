# ============================================================
# Ring vaccination simulator v1.19.1
#
# Full operating characteristics framework
#
# Primary estimand:
#   VE_immediate_vs_delayed
#
# This is NOT biological VE.
# TRUE biological VE is a scenario assumption.
#
# Evaluates:
#   - direct contact inference
#   - CoC operational benefit
#   - vaccine requirements
#   - ethics/inference trade-offs
# ============================================================


# -------------------------
# CONFIGURATION
# -------------------------

DEBUG_MODE <- FALSE

# v1.18.1 execution options
PARALLEL <- FALSE
N_WORKERS <- 6

RESUME <- TRUE

# -------------------------
# RUN-SPECIFIC OUTPUT FOLDER
# -------------------------

RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")

RUN_DIR <- file.path(
  "runs",
  paste0("ring_vax_run_", RUN_ID)
)

file.copy(
  "ring_vax_modeller_v1.20.0_run_managed.R",
  file.path(RUN_DIR, "model_script.R")
)

CACHE_DIR <- file.path(RUN_DIR, "cache")
CHECKPOINT_DIR <- file.path(RUN_DIR, "checkpoints")
OUTPUT_RUN_DIR <- file.path(RUN_DIR, "outputs")

dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_RUN_DIR, recursive = TRUE, showWarnings = FALSE)

DEBUG_N_DESIGNS <- 5
N_SIM <- 2000
DEBUG_N_SIM <- 50

EMPIRICAL_RING_FILE <- "participants_per_index_R_210621.0311.csv"

TRUE_VE_VALUES <- c(0.3,0.5,0.7)
SAR_VALUES <- c(0.02,0.05,0.10)
SAR_CV_VALUES <- c(0,0.25,0.5)

RING_NUMBERS <- c(400,600,1000,1500)
ALLOCATION_RATIOS <- c(1,1.5,2,3)
DELAYS <- c(9,14,21,28)

INCLUDE_COC <- TRUE
COC_RISK_MULTIPLIER <- 1

set.seed(20260821)


# -------------------------
# PRE-WARM EPIDEMIC CACHE
# -------------------------
# Generate shared epidemic structures before parallel execution.
# This prevents multiple workers simultaneously generating the same cache.

warm_keys <- unique(
  design_grid[, c("rings","SAR","SAR_CV")]
)

invisible(
  apply(
    warm_keys,
    1,
    function(z){

      get_epidemic_cache(
        rings,
        as.numeric(z[["rings"]]),
        n_run,
        as.numeric(z[["SAR"]]),
        as.numeric(z[["SAR_CV"]])
      )

    }
  )
)

if(PARALLEL){

  if(!requireNamespace("future", quietly=TRUE) ||
     !requireNamespace("future.apply", quietly=TRUE) ||
     !requireNamespace("progressr", quietly=TRUE)){
    stop("Install future, future.apply and progressr")
  }

  library(future)
  library(future.apply)
  library(progressr)

  plan(multisession, workers=N_WORKERS)
}


# -------------------------
# LOAD RINGS
# -------------------------

load_empirical_rings <- function(path){

  raw <- read.csv(path)

  rings <- data.frame(
    ring_id=raw$r_case_index_id,
    direct_contacts=raw$direct_contact,
    coc_contacts=raw$contact_of_contact
  )

  subset(
    rings,
    !is.na(direct_contacts) &
      direct_contacts > 0
  )
}


sample_rings <- function(rings,n){
  rings[sample(seq_len(nrow(rings)), n, replace=TRUE),]
}




# -------------------------
# SAR HETEROGENEITY
# -------------------------

# Draw ring-specific SAR values while preserving the requested
# population mean SAR. CV=0 reproduces the original model.
draw_ring_SAR <- function(n, mean_SAR, CV){

  if(CV == 0){
    return(rep(mean_SAR,n))
  }

  variance <- (CV * mean_SAR)^2

  max_variance <- mean_SAR * (1-mean_SAR)

  if(variance >= max_variance){
    variance <- 0.99 * max_variance
  }

  common <- mean_SAR*(1-mean_SAR)/variance - 1

  rbeta(
    n,
    mean_SAR * common,
    (1-mean_SAR) * common
  )
}

# -------------------------
# RANDOMISATION
# -------------------------

randomise_rings <- function(n,ratio){

  arms <- c(
    rep("vaccine", round(ratio*2)),
    rep("delayed",2)
  )

  sample(rep(arms,length.out=n))
}


# -------------------------
# BIOLOGY
# -------------------------

protection_curve <- function(days){

  approx(
    x=c(0,7,14,21,28),
    y=c(0,0.25,0.60,0.80,1),
    xout=pmax(days,0),
    rule=2
  )$y
}


simulate_contacts <- function(n,vaccination_day,true_VE,SAR){

  if(n==0) return(0)

  exposure_day <- sample(0:28,n,replace=TRUE)

  protection <- protection_curve(
    exposure_day-vaccination_day
  )

  risk <- SAR*(1-true_VE*protection)

  # Preserve heterogeneous risks while avoiding one Bernoulli draw
  # per contact by grouping identical probabilities.
  probs <- table(round(risk, 10))

  sum(
    mapply(
      function(p, m) rbinom(1, m, p),
      as.numeric(names(probs)),
      as.numeric(probs)
    )
  )
}


simulate_population <- function(ring,
                                arm,
                                delay,
                                true_VE,
                                SAR){

  vacc_day <- ifelse(
    arm=="vaccine",
    0,
    delay
  )

  direct_cases <- sum(
    simulate_contacts(
      ring$direct_contacts,
      vacc_day,
      true_VE,
      SAR
    )
  )

  coc_cases <- 0

  if(INCLUDE_COC){

    coc_cases <- sum(
      simulate_contacts(
        ring$coc_contacts,
        vacc_day,
        true_VE,
        SAR*COC_RISK_MULTIPLIER
      )
    )
  }


  data.frame(
    arm=arm,
    direct_contacts=ring$direct_contacts,
    coc_contacts=ring$coc_contacts,
    direct_cases=direct_cases,
    coc_cases=coc_cases
  )
}


# -------------------------
# ANALYSIS
# -------------------------

analyse_trial <- function(dat){

  vaccine <- subset(dat,arm=="vaccine")
  delayed <- subset(dat,arm=="delayed")

  risk_vaccine <-
    sum(vaccine$direct_cases) /
    sum(vaccine$direct_contacts)

  risk_delayed <-
    sum(delayed$direct_cases) /
    sum(delayed$direct_contacts)


  VE_immediate_vs_delayed <-
    1-risk_vaccine/risk_delayed


  expected_direct <-
    risk_delayed *
    sum(dat$direct_contacts)

  expected_coc <-
    (sum(delayed$coc_cases)/
       sum(delayed$coc_contacts)) *
    sum(dat$coc_contacts)


  total_observed <-
    sum(dat$direct_cases)+
    sum(dat$coc_cases)


  data.frame(
    VE_immediate_vs_delayed=
      VE_immediate_vs_delayed,

    direct_cases_prevented=
      expected_direct-sum(dat$direct_cases),

    coc_cases_prevented=
      expected_coc-sum(dat$coc_cases),

    total_cases_prevented=
      expected_direct+
      expected_coc-
      total_observed,

    vaccine_doses=
      sum(vaccine$direct_contacts)+
      sum(vaccine$coc_contacts)
  )
}


# -------------------------
# DESIGN RUNNER
# -------------------------

# Cache reusable trial structure.
# Ring composition, randomisation and latent SAR are independent of
# vaccine efficacy and delay, so do not regenerate them for every scenario.


# ============================================================
# v1.19.2 epidemic cache engine
#
# Expensive transmission structures are generated once.
# Intervention scenarios (VE, delay, allocation) are evaluated
# against cached epidemic structures.
# ============================================================

EPIDEMIC_CACHE_FILE <- file.path(
  CACHE_DIR,
  "epidemic_cache_v1.19.8.rds"
)
VALIDATE_CACHE <- FALSE

# Cache diagnostics
cache_status <- function(){

  mem_keys <- ls(envir = epidemic_cache_store)

  disk_files <- list.files(
    path = CACHE_DIR,
    pattern = "^cache_.*\\.rds$",
    full.names = TRUE
  )

  list(
    memory_entries = length(mem_keys),
    memory_keys = mem_keys,
    disk_entries = length(disk_files),
    disk_files = disk_files
  )

}


generate_epidemic_cache <- function(
  rings,
  n_rings,
  n_sim,
  SAR,
  SAR_CV
){

  lapply(seq_len(n_sim), function(i){

    sampled <- sample_rings(
      rings,
      n_rings
    )

    sampled$ring_SAR <- draw_ring_SAR(
      n_rings,
      SAR,
      SAR_CV
    )

    sampled

  })

}


get_epidemic_cache <- function(
  rings,
  n_rings,
  n_sim,
  SAR,
  SAR_CV,
  cache_env = epidemic_cache_store
){

  key <- paste(
    n_rings,
    n_sim,
    SAR,
    SAR_CV,
    sep="_"
  )

  if(exists(key, envir=cache_env)){

    message("CACHE HIT (memory): ", key)

    return(
      get(key, envir=cache_env)
    )

  }

  disk_cache <- file.path(
    CACHE_DIR,
    paste0("cache_", key, ".rds")
  )

  if(file.exists(disk_cache)){

    message("CACHE HIT (disk): ", disk_cache)

    cache <- readRDS(disk_cache)

    assign(
      key,
      cache,
      envir = cache_env
    )

    return(cache)

  }

  message("CACHE MISS -> GENERATING: ", key)

  # Generate once
  cache <- generate_epidemic_cache(
    rings,
    n_rings,
    n_sim,
    SAR,
    SAR_CV
  )

  # Store in worker-local cache
  assign(
    key,
    cache,
    envir = cache_env
  )

  # Also persist to disk so parallel workers / later sessions can reuse it
  cache_file <- file.path(
    CACHE_DIR,
    paste0("cache_", key, ".rds")
  )

  try(
    saveRDS(
      cache,
      file = cache_file
    ),
    silent = TRUE
  )

  message("CACHE WRITE: ", cache_file)

  cache

}


epidemic_cache_store <- new.env(
  parent=emptyenv()
)


run_design <- function(
  rings,
  n_rings,
  n_sim,
  true_VE,
  SAR,
  SAR_CV,
  allocation_ratio,
  delay,
  scenario_cache=NULL
){

  if(is.null(scenario_cache)){

    scenario_cache <- get_epidemic_cache(
      rings,
      n_rings,
      n_sim,
      SAR,
      SAR_CV
    )

  }


  sims <- lapply(
    scenario_cache,
    function(sampled_epidemic){

      # Allocation changes only at intervention stage
      allocated <- sampled_epidemic

      allocated$arm <- randomise_rings(
        nrow(allocated),
        allocation_ratio
      )


      dat <- do.call(
        rbind,
        lapply(
          seq_len(nrow(allocated)),
          function(j){

            simulate_population(
              allocated[j,],
              allocated$arm[j],
              delay,
              true_VE,
              allocated$ring_SAR[j]
            )

          }
        )
      )

      analyse_trial(dat)

    }
  )


  x <- do.call(
    rbind,
    sims
  )


  data.frame(

    true_biological_VE=true_VE,
    SAR=SAR,
    SAR_CV=SAR_CV,
    rings=n_rings,
    allocation_ratio=allocation_ratio,
    delay=delay,

    VE_immediate_vs_delayed=
      mean(x$VE_immediate_vs_delayed,na.rm=TRUE),

    sd_VE=
      sd(x$VE_immediate_vs_delayed,na.rm=TRUE),

    direct_cases_prevented=
      mean(x$direct_cases_prevented),

    coc_cases_prevented=
      mean(x$coc_cases_prevented),

    total_cases_prevented=
      mean(x$total_cases_prevented),

    vaccine_doses=
      mean(x$vaccine_doses),

    cases_per_1000_doses=
      mean(x$total_cases_prevented) /
      mean(x$vaccine_doses) *
      1000,

    precision=
      1/(1+sd(x$VE_immediate_vs_delayed))

  )

}



# -------------------------
# EXECUTION
# -------------------------

rings <- load_empirical_rings(
  EMPIRICAL_RING_FILE
)


design_grid <- expand.grid(
  true_VE=TRUE_VE_VALUES,
  SAR=SAR_VALUES,
  SAR_CV=SAR_CV_VALUES,
  rings=RING_NUMBERS,
  allocation_ratio=ALLOCATION_RATIOS,
  delay=DELAYS
)

design_grid$design_id <- seq_len(nrow(design_grid))


if(DEBUG_MODE){

  # Stratified debug subset:
  # ensure debug runs include all SAR_CV scenarios rather than only
  # the first rows produced by expand.grid()

  n_per_group <- ceiling(
    DEBUG_N_DESIGNS / length(unique(design_grid$SAR_CV))
  )

  design_grid <- design_grid %>%
    dplyr::group_by(
      SAR_CV
    ) %>%
    dplyr::slice_head(
      n = n_per_group
    ) %>%
    dplyr::ungroup()

  n_run <- DEBUG_N_SIM

}else{

  n_run <- N_SIM

}


checkpoint_file <- function(id){

  file.path(
    CHECKPOINT_DIR,
    paste0("design_",id,".rds")
  )

}


completed <- function(){

  files <- list.files(
    CHECKPOINT_DIR,
    pattern="design_.*\\.rds"
  )

  if(length(files)==0)
    return(integer())

  as.numeric(
    gsub(
      "design_|\\.rds",
      "",
      files
    )
  )
}


completed_ids <- completed()


# ------------------------------------------------------------
# v1.19.3 cache architecture
#
# Single cache layer only. Epidemic structures are cached before
# intervention evaluation. Allocation is deliberately applied later
# so the same epidemic realisations are compared across allocation
# strategies.
# ------------------------------------------------------------

run_one_design <- function(i){

  if(RESUME && i %in% completed_ids){

    return(readRDS(checkpoint_file(i)))

  }


  d <- design_grid[i,]


  cache <- get_epidemic_cache(
    rings,
    d$rings,
    n_run,
    d$SAR,
    d$SAR_CV
  )

  result <- run_design(
    rings,
    d$rings,
    n_run,
    d$true_VE,
    d$SAR,
    d$SAR_CV,
    d$allocation_ratio,
    d$delay,
    scenario_cache=cache
  )


  result$design_id <- i


  saveRDS(
    result,
    checkpoint_file(i)
  )


  result
}


if(PARALLEL){

  handlers(global=TRUE)
  handlers("progress")

  results <- with_progress({

    p <- progressor(
      along=seq_len(nrow(design_grid))
    )

    future_lapply(
      seq_len(nrow(design_grid)),
      function(i){

        x <- run_one_design(i)
        p()
        x

      },
      future.seed=TRUE
    )

  })

}else{

  results <- lapply(
    seq_len(nrow(design_grid)),
    run_one_design
  )

}


final_output <- do.call(rbind,results)


# -------------------------
# PARETO FRONTIER
# -------------------------

pareto_flag <- function(df){

  keep <- rep(TRUE,nrow(df))

  for(i in seq_len(nrow(df))){

    dominated <- any(
      df$total_cases_prevented >= df$total_cases_prevented[i] &
      df$precision >= df$precision[i] &
      df$vaccine_doses <= df$vaccine_doses[i] &
      (
        df$total_cases_prevented > df$total_cases_prevented[i] |
        df$precision > df$precision[i] |
        df$vaccine_doses < df$vaccine_doses[i]
      )
    )

    if(dominated)
      keep[i] <- FALSE

  }

  keep
}


# Pareto frontier is calculated within epidemiological scenarios.
# Do not compare different SAR/VE/delay scenarios against each other.
final_output$pareto_optimal <- FALSE

scenario_groups <- split(
  seq_len(nrow(final_output)),
  interaction(
    final_output$true_biological_VE,
    final_output$SAR,
    final_output$SAR_CV,
    final_output$delay,
    drop = TRUE
  )
)

for(g in scenario_groups){

  final_output$pareto_optimal[g] <-
    pareto_flag(final_output[g, ])

}


write.csv(
  final_output,
  file.path(
    OUTPUT_RUN_DIR,
    paste0("ring_vax_operating_characteristics_", RUN_ID, ".csv")
  ),
  row.names=FALSE
)


write.csv(
  subset(final_output,pareto_optimal),
  file.path(
    OUTPUT_RUN_DIR,
    paste0("ring_vax_pareto_frontier_", RUN_ID, ".csv")
  ),
  row.names=FALSE
)


saveRDS(
  list(
    RUN_ID = RUN_ID,
    timestamp = Sys.time(),
    parameters = list(
      TRUE_VE_VALUES = TRUE_VE_VALUES,
      SAR_VALUES = SAR_VALUES,
      SAR_CV_VALUES = SAR_CV_VALUES,
      RING_NUMBERS = RING_NUMBERS,
      N_SIM = N_SIM,
      DEBUG_MODE = DEBUG_MODE
    )
  ),
  file.path(RUN_DIR, "config.rds")
)

print(final_output)



old_results


