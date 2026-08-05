# core_design.R -- turn declared roles/contrasts/factorial/covariates into
# the derived factors and numeric variables every test module models
# against. This is where "study strings live in JSON, not in code" is
# actually enforced: no module below this file ever compares a group value
# to a literal string.

flex_derive_time <- function(df, session_col = "session") {
  levels_sorted <- sort(unique(df[[session_col]]))
  df$time_f <- factor(df[[session_col]], levels = levels_sorted)
  df$time_numeric <- as.numeric(df$time_f) - 1
  df
}

# Maps group values to an "intervention"/"control" factor per
# design.contrasts. Rows whose group value is in neither list are dropped
# (with a warning) rather than silently coerced -- mirrors
# 01_primary_lmm.R's unknown_groups warning.
flex_derive_contrast <- function(df, cfg, group_col = "group") {
  contrasts <- cfg$design$contrasts
  if (is.null(contrasts) || is.null(contrasts$intervention) || is.null(contrasts$control)) {
    stop("design.contrasts.intervention/control not set in config -- required for this test")
  }
  intervention_groups <- unlist(contrasts$intervention)
  control_group <- unlist(contrasts$control)

  known <- c(intervention_groups, control_group)
  unknown <- setdiff(unique(df[[group_col]]), known)
  if (length(unknown)) {
    warning(sprintf("participants have group values not in design.contrasts: %s (these subjects will be dropped)",
                     paste(unknown, collapse = ", ")))
  }
  df <- df[df[[group_col]] %in% known, , drop = FALSE]
  df$intervention <- factor(ifelse(df[[group_col]] %in% intervention_groups, "intervention", "control"),
                             levels = c("control", "intervention"))
  df$style <- factor(df[[group_col]], levels = c(control_group, intervention_groups))
  df
}

# Splits the intervention arms along a named factorial axis (e.g. "social",
# "duration") per design.factorial.<axis>.<level>: [group values]. Rows
# whose group doesn't appear under any level of this axis are dropped
# (factorial designs are only defined within the arms that carry the axis
# -- e.g. control has no "social context", matching 13_factorial_social_duration.R).
flex_derive_factorial <- function(df, cfg, axis_name, group_col = "group") {
  factorial_cfg <- cfg$design$factorial
  if (is.null(factorial_cfg) || is.null(factorial_cfg[[axis_name]])) {
    stop(sprintf("design.factorial.%s not set in config -- required for this test", axis_name))
  }
  axis <- factorial_cfg[[axis_name]]
  level_names <- names(axis)
  lookup <- setNames(rep(level_names, times = vapply(axis, length, integer(1))), unlist(axis))

  df[[axis_name]] <- unname(lookup[df[[group_col]]])
  n_dropped <- sum(is.na(df[[axis_name]]))
  if (n_dropped) {
    df <- df[!is.na(df[[axis_name]]), , drop = FALSE]
  }
  df[[axis_name]] <- factor(df[[axis_name]], levels = level_names)
  df
}

# Applies every design.covariates entry, creating the transformed column
# under `as` (or the original column name if `as` is omitted).
#   transform "z"       -> scale(as.numeric(col))
#   transform "rank_z"  -> scale(rank(as.numeric(col))) -- see
#                          01_primary_lmm.R's age_z rationale: caps leverage
#                          from a sparse tail without dropping subjects.
#   transform "none"    -> as.numeric(col), unchanged
#   type "factor"       -> as.factor(col), transform ignored
flex_apply_covariates <- function(df, cfg) {
  covariates <- cfg$design$covariates
  if (is.null(covariates)) return(df)
  for (spec in covariates) {
    col <- spec$col
    target <- spec$as %||% col
    if (!col %in% names(df)) stop(sprintf("design.covariates references column '%s' not present in participants", col))
    if (identical(spec$type, "factor")) {
      df[[target]] <- factor(df[[col]])
      next
    }
    transform <- spec$transform %||% "none"
    numeric_col <- suppressWarnings(as.numeric(df[[col]]))
    df[[target]] <- switch(transform,
      "z" = as.numeric(scale(numeric_col)),
      "rank_z" = as.numeric(scale(rank(numeric_col))),
      "none" = numeric_col,
      stop(sprintf("unknown covariate transform '%s' for column '%s'", transform, col))
    )
  }
  df
}

# profile$transform: "log" or "identity" -- the log(volume) vs raw-scale
# thickness choice documented in 01_primary_lmm.R / 29_hippo_thickness_lmm.R.
flex_apply_profile_transform <- function(values, profile) {
  switch(profile$transform,
    "log" = log(values),
    "identity" = values,
    stop(sprintf("unknown profile transform '%s'", profile$transform))
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
