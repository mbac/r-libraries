# Descriptive statistics according to MB

# Summary function for scalar variables
suppressMessages(library(tidyverse))
suppressMessages(library(rcompanion))
suppressMessages(library(DescTools))
pacman::p_load_gh("tidymodels/infer")

# This function uses the infer package to produce median/BCa confidence
# intervals. Order of arguments is mandated by gtsummary package…
myBootCI <- function(data, value,
                     group = NULL,
                     statistic,
                     ci_type = "bias-corrected",
                     output = "tibble",
                     digits = 0,
                     ...) {
  value <- ensym(value)

  boot <- data |>
    specify(response = !!value) |>
    generate(reps = 10000, type = "bootstrap") |>
    calculate(stat = statistic) |>
    suppressWarnings()

  point <- data |> observe(
    response = !!value,
    stat = statistic)

  ci <- get_confidence_interval(boot,
                                type = ci_type,
                                point_estimate = point
  ) |> suppressWarnings() |> suppressMessages()

  # Compose dataframe and round
  result_df <- add_column(ci, point, .before = 1) |>
    map_dfc(~ myround(.x, digits = digits))

  if (output == "tibble")
    result <- result_df

  else if (output == "string")
    result <- glue::glue("{result_df$stat} ({result_df$lower_ci} — {result_df$upper_ci})")

  return(result)
}

# This function is a wrapper to the `infer` package; here we compute diff in
# means/medians and their bootstrapped CIs. Requires grouping variable


# order A string vector of specifying the order in which the levels of the
# explanatory variable should be ordered for subtraction (or division for
# ratio-based statistics), where order = c("first", "second") means ("first" -
# "second"), or the analogue for ratios. Needed for inference on difference in
# means, medians, proportions, ratios, t, and z statistics.

mydiffCI <- function(data,
                     variable,
                     by,
                     group = NULL,
                     type = NULL,
                     conf.level = 0.95,
                     statistic = "median",
                     triage = NULL,
                     ci_type = "bias-corrected",
                     output = "tibble",
                     digits = 0,
                     na.rm = TRUE,
                     ...) {

  # If a statistic is not specified, look for global variable
  if (is.null(statistic)) {
    # If a global was set…
    if (exists(default_stat) & is.character(default_stat))
      statistic <- default_stat
    else {
      warning("Attention! argument statistic was not specified, and
                    there's no global default_stat.

                    ARBITRARILY DEFAULTING TO DIFF OF MEDIANS")
      statistic <- "median"
    }
  }

  # We need an order for variables in case of differences
  if (is.null(triage)) {
    tmp <- data[[by]]
    triage <- unique(tmp)

    warning("ATTENTION! No triage was given for diff operations. Default is
               first_appearing - second_appearing")
  }

  # Defuse variables per tidyverse style
  variable <- ensym(variable)
  # data <- ensym(data)
  by <- ensym(by)

  # print(print(as.list(class(match.call()))))

  # change the function arg into one compatible with infer syntax
  stat <- glue::glue("diff in {statistic}s")
  # Statistic name which may be reported in the printed table
  test_name = glue::glue("Bootstrapped {statistic} (95% confidence intervals)")
  # Recover dot-dot-dot args. This is a list.
  varargs <- rlang::exprs(...)


  # Calculate point estimate:
  point <- data |>
    #specify(expr(!!value ~ expr(!!group))) |>
    specify(expr(!!variable ~ !!by)) |>
    calculate(stat = stat, order = triage, na.rm = TRUE) |>
    suppressWarnings()

  boot <- data |>
    specify(expr(!!variable ~ !!by)) |>
    generate(reps = 1000, type = "bootstrap") |>
    calculate(stat = stat, order = triage, na.rm = TRUE) |>
    suppressWarnings()


  result_df <- get_confidence_interval(boot,
                                       type = ci_type,
                                       point_estimate = point
  ) |> suppressWarnings() |> suppressMessages()

  # Compose dataframe and round
  result_df <- result_df |>
    mutate(# Include statistic name (sort of)
      method = test_name,
      # Include point estimate in confidence intervals
      estimate = point$stat,
      .before = 1)


  if (output == "tibble") {
    result <- result_df |>
      rename(
        conf.low = lower_ci,
        conf.high = upper_ci)

    return(result)
  }

  else if (output == "string")
    result <- glue::glue("{result_df$stat} ({result_df$lower_ci} — {result_df$upper_ci})")

}


# Takes a tibble with interesting data, extracts numeric variables and calculates
# bootstrap medians and CI; if there is a group, set var name (string)
medianizer <- function (data,
                        group_var = NULL) {
  result <- tibble()

  # Only select numeric, non-factor data
  data <- select(
    data,
    where(is.numeric) &
      !where(is.instant) &
      !where(is.timespan) &
      !where(is.POSIXct) &
      !where(is.factor) &
      # Double escaped regexp matching:
      # regexp: id(._-) or (._-)id, to eliminate obvious id variables
      !matches("(^id(\\.|\\-|\\_)?)|((\\.|\\-|\\_)+id$)")
  )

  # To better catch errors, create an empty tibble
  # with same columns as groupwiseMedian()
  dummy <- tibble(
    .id = as.character(NA),
    n = as.numeric(NA),
    Median = as.numeric(NA),
    Boot.median = as.numeric(NA),
    Conf.level = as.numeric(NA),
    Bca.lower = as.numeric(NA),
    Bca.upper = as.numeric(NA)
  )
  # Error-safe function
  gMedian <- safely(groupwiseMedian, dummy)

  for (i in seq_along(data)) {
    curr_var <- names(data[i])
    loop_result <- gMedian(
      data = data,
      var = curr_var,
      boot = T,
      exact = T,
      group = group_var,
      digits = 1,
      R = 10000
    )
    # Because we're using a "safely" function with a dummy result, if
    # there's an error, we get a $result with NAs (dummy above) and a $error
    # with the text. For now, let's just focus on $result.
    result <- loop_result$result |>
      # Now, only select useful variables. We cannot use the specific name
      # because one of the values just might not be there
      select(
        contains(
          c("median",
            "bca",
            "boot",
            "exact",
            "basic"
          )
        )|
          matches(
            "^n$"
          )

      ) |>
      # This adds the "variable name" column
      mutate(var = curr_var) |>
      # This adds current variable to the result
      bind_rows(result)
  }

  result <- result |>
    # Variable name first
    select(var, n, everything()) |>
    # Reverse row order so they match the order of variables provided in data
    arrange(-row_number())
  return(result)


}

############################
#
# Convert Logit to Probability
#
############################
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# Printable values ========================

# from broman package
myround <- function (x, digits = 1, p=FALSE, printable = FALSE)
{
  # If it's a p-value, just do this:
  if (p) {
    if (round(x, 3) > 0.001) {
      x <- myround(x, 3, p = F, printable = F)
      if (printable)
        x <- paste("=", as.character(x))
    }
    else {
      x <- 0.001
      if (printable)
        x <- paste("≤", as.character(x))
    }
    return(x)
  }

  else if (digits <= 1)
    return(round(x, digits))

  else if (length(digits) > 1) {
    digits <- digits[1]
    warning("Using only digits[1]")
  }

  number <- format(x, nsmall = digits) |>
    as.double()

  number
}

mybindesc <- function (data,
                       fact,
                       label = NULL,
                       sort = FALSE,
                       digits = 0
) {

  if (is.null(label) | !is.character(label))
    return("Sorry, you have to specify a factor level for now (as string)")


  counts <- data |> dplyr::group_by({{fact}}) |>
    dplyr::summarise(count = n()) |>
    dplyr::mutate(pct = count/sum(count) * 100,
                  pct = myround(pct, digits = digits))

  # Search and count desired value occurrences
  if (!is.null(label)) {
    x <- counts |> filter({{fact}} == label)
    return(glue::glue('
                          {x$count} ({x$pct}%)')
    )
  }

}

