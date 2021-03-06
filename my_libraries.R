
pacman::p_load_current_gh("tidyverse/tidyverse",
                          "tidymodels/tidymodels",
                          "sfirke/janitor",
                          "ddsjoberg/gtsummary",
                          "mayoverse/arsenal",
                          "rstudio/gt",
                          "jeroen/curl"
                          ) |>
  suppressMessages()

pacman::p_load(rcompanion) |> suppressMessages()

pacman::p_load("readxl") |> suppressMessages()

# File with useful functions for descriptives:
# source(file = "descriptives.R")
# Disabled when working remotely
