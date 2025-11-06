# wm-gwas_source.R
# 
# Copyright (c) Andrea Estandia, 2025, except where indicated
# Date Created: 2025-10-06


# --------------------------------------------------------------------------
# REQUIRES
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(showtext)

})


text_size = 11
# --------------------------------------------------------------------------
# PATHS
# --------------------------------------------------------------------------

data_path <- file.path(getwd(), "data")
reports_path <- file.path(getwd(), "reports")
figures_path <- file.path(getwd(), "reports", "plots")

if (!dir.exists(data_path)) {
  dir.create(data_path, recursive = TRUE)
}

if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}

if (!dir.exists(reports_path)) {
  dir.create(reports_path, recursive = TRUE)
}

font_add_google("Roboto Condensed", "roboto_condensed")
showtext_auto() 


'%!in%' <- function(x,y)!('%in%'(x,y))
CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 