required_packages <- c(
  "future",
  "future.apply",
  "dplyr",
  "ggplot2",
  "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}

