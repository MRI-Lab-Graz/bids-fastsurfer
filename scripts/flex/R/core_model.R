# core_model.R -- formula assembly from a measure's profile + the study's
# declared covariates, and the random-slope -> random-intercept singularity
# fallback used by every LMM-based test (verbatim from 01_primary_lmm.R).

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

flex_covariate_terms <- function(cfg) {
  covariates <- cfg$design$covariates
  if (is.null(covariates)) return(character(0))
  vapply(covariates, function(spec) spec$as %||% spec$col, character(1))
}

# Assembles the fixed-effects RHS (everything except the response and
# random-effects term) from a profile + the covariate list: the crossed
# factor term (default "intervention * time_f"), hemisphere (only when the
# profile pools hemispheres), every declared covariate, and etiv_z (only
# when the profile's etiv_covariate is TRUE) -- this is the single place
# that replaces the seven near-duplicate LMM scripts' hand-written formulas.
flex_build_rhs <- function(profile, covariate_terms, factor_term = "intervention * time_f") {
  terms <- factor_term
  if (identical(profile$hemisphere, "pooled")) terms <- c(terms, "hemisphere")
  if (length(covariate_terms)) terms <- c(terms, covariate_terms)
  if (isTRUE(profile$etiv_covariate)) terms <- c(terms, "etiv_z")
  paste(terms, collapse = " + ")
}

# Fits response ~ rhs + (1 + time_numeric | subject_id), falling back to
# (1 | subject_id) if that model fails to fit or is singular -- verbatim
# logic from 01_primary_lmm.R's fit_one_roi()/fit_re().
flex_fit_lmm_with_fallback <- function(data, rhs, response = "y") {
  fit_re <- function(re_formula) {
    f <- stats::as.formula(paste(response, "~", rhs, "+", re_formula))
    lmerTest::lmer(f, data = data, REML = TRUE)
  }
  model <- tryCatch(fit_re("(1 + time_numeric | subject_id)"), error = function(e) NULL)
  re_structure <- "random_slope"
  if (is.null(model) || lme4::isSingular(model, tol = 1e-4)) {
    model <- tryCatch(fit_re("(1 | subject_id)"), error = function(e) NULL)
    re_structure <- "random_intercept_only"
  }
  list(model = model, re_structure = re_structure)
}

flex_emmeans_contrast <- function(model, spec, method = "revpairwise") {
  em <- tryCatch(emmeans::emmeans(model, spec), error = function(e) NULL)
  if (is.null(em)) return(NULL)
  as.data.frame(emmeans::contrast(em, method = method))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
