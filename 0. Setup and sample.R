# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 0.1 Packages -----------------------------------------------------------

pkgs <- c("haven", "dplyr", "survey")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(haven)
library(dplyr)
library(survey)


# ---- 0.2 Paths (edit data_dir if the .dta file is stored elsewhere) ----------

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")
dir.create(out_dir, showWarnings = FALSE)


# ---- 0.3 Read the prepared Stata file ---------------------------------------

dat <- read_dta(file.path(data_dir, "RLS Public Use File_2_Clean.dta"))
stopifnot(nrow(dat) == 36908)


# ---- 0.4 Convert categorical covariates to factors --------------------------

fac_vars <- c("rel_curr9", "rel_raised9", "partyid3", "marital4",
              "partner_rel", "sp_relimp", "income4", "educrec",
              "ideo", "agegrp", "usgen", "woman")

clean_factor <- function(x) {
  f <- as_factor(x, levels = "labels")
  levels(f) <- sub("^-?\\d+\\.\\s*", "", levels(f))
  droplevels(f)
}

dat <- dat %>% mutate(across(all_of(fac_vars), clean_factor))

# Remaining Stata value labels are removed so that numeric variables are
# plain numbers in all later stages.
dat <- zap_labels(dat)

# Check: print the categories of each factor
print(lapply(dat[fac_vars], levels))


# ---- 0.5 LPA indicators and analytic-sample flag ----------------------------

indices <- c("religbehav", "religbelif", "spirbehav", "spirbelif")

# Expected ranges: 9-46, 4-16, 5-25, 5-18
print(sapply(dat[indices], range, na.rm = TRUE))

# Analytic sample: complete on all four indices (equivalent to complete on all
# 23 index items, because each index requires all of its items)
dat$lpa_sample <- complete.cases(dat[, indices])
dat$lpa_f <- factor(dat$lpa_sample, levels = c(TRUE, FALSE),
                    labels = c("Included", "Excluded"))

stopifnot(sum(dat$lpa_sample) == 31470)
stopifnot(sum(!is.na(dat$half)) == 31470)   # agrees with the Stata split-half flag

table(dat$lpa_f)


# ---- 0.6 Missingness table --------------------------------------------------

items <- c("attndperrlsr", "attndonrlsr", "synmemr", "prayfreqr", "gracefreqr",
           "prac_ar", "prac_br", "prac_cr", "prac_dr",
           "relimpr", "relperr", "gdcertainr", "afterlife",
           "spiract_ar", "spiract_br", "spiract_cr", "spiract_dr",
           "spiract_fr",
           "spirperr", "soulr", "spirworld2", "exp_ar", "exp_dr")

covariates <- c("agegrp", "woman", "educrec", "income4", "ideo", "partyid3",
                "usgen", "marital4", "rel_curr9", "rel_raised9",
                "child_relig", "ch_attend", "ch_relschool", "ch_religed",
                "partner_rel", "sp_relimp", "wellbeing",
                "lifedir", "lifespir")

miss_vars <- c(items, indices, covariates)

miss_tab <- data.frame(
  variable    = miss_vars,
  block       = c(rep("Index item", length(items)),
                  rep("LPA index", length(indices)),
                  rep("Covariate / outcome", length(covariates))),
  n_missing   = sapply(dat[miss_vars], function(v) sum(is.na(v))),
  row.names   = NULL
)
miss_tab$pct_missing <- round(100 * miss_tab$n_missing / nrow(dat), 2)

write.csv(miss_tab, file.path(out_dir, "stage0_missingness.csv"), row.names = FALSE)
print(miss_tab)


# ---- 0.7 Survey design ------------------------------------------------------
# Full-sample weight plus 500 Rao-Wu rescaling bootstrap replicate weights,
# as documented in the RLS methodology report.

stopifnot(length(grep("^repwt_", names(dat))) == 500)

des <- svrepdesign(data             = dat,
                   weights          = ~weight,
                   repweights       = "repwt_[0-9]+",
                   type             = "bootstrap",
                   combined.weights = TRUE,
                   mse              = TRUE)


# ---- 0.8 Comparison of included and excluded respondents --------------------
# Categorical covariates: weighted column percentages, unweighted Ns,
# Rao-Scott F test, and Cramer's V computed from the weighted table rescaled
# to the unweighted N (an approximate, descriptive effect size).
# Continuous covariates: weighted means, design-based t test, and a
# standardised mean difference (difference / weighted SD of the full sample).

cat_vars  <- c("agegrp", "woman", "educrec", "income4", "ideo", "partyid3",
               "usgen", "marital4", "rel_curr9", "rel_raised9", "partner_rel")
cont_vars <- c("child_relig", "ch_attend", "wellbeing")

compare_cat <- function(v) {
  f    <- as.formula(paste0("~", v, " + lpa_f"))
  tabw <- svytable(f, des)
  pct  <- round(100 * prop.table(tabw, margin = 2), 1)
  n    <- table(dat[[v]], dat$lpa_f)[rownames(pct), , drop = FALSE]
  test <- svychisq(f, des, statistic = "F")
  
  tab_scaled <- tabw / sum(tabw) * sum(n)
  chi2 <- suppressWarnings(chisq.test(tab_scaled, correct = FALSE)$statistic)
  v_cramer <- sqrt(as.numeric(chi2) / (sum(n) * (min(dim(tabw)) - 1)))
  
  k <- nrow(pct)
  data.frame(
    variable     = v,
    category     = rownames(pct),
    n_included   = as.vector(n[, "Included"]),
    pct_included = as.vector(pct[, "Included"]),
    n_excluded   = as.vector(n[, "Excluded"]),
    pct_excluded = as.vector(pct[, "Excluded"]),
    diff_pp      = as.vector(pct[, "Included"] - pct[, "Excluded"]),
    p_value      = c(signif(test$p.value, 3), rep(NA, k - 1)),
    effect_size  = c(round(v_cramer, 3), rep(NA, k - 1)),
    effect_type  = c("Cramer's V", rep(NA, k - 1))
  )
}

compare_cont <- function(v) {
  m    <- svyby(as.formula(paste0("~", v)), ~lpa_f, des, svymean, na.rm = TRUE)
  est  <- coef(m)
  sdw  <- sqrt(as.numeric(coef(svyvar(as.formula(paste0("~", v)), des, na.rm = TRUE))))
  test <- svyttest(as.formula(paste0(v, " ~ lpa_f")), des)
  n    <- tapply(!is.na(dat[[v]]), dat$lpa_f, sum)
  
  data.frame(
    variable     = v,
    category     = "Mean",
    n_included   = n[["Included"]],
    pct_included = round(est[["Included"]], 2),
    n_excluded   = n[["Excluded"]],
    pct_excluded = round(est[["Excluded"]], 2),
    diff_pp      = round(est[["Included"]] - est[["Excluded"]], 2),
    p_value      = signif(test$p.value, 3),
    effect_size  = round((est[["Included"]] - est[["Excluded"]]) / sdw, 3),
    effect_type  = "Standardised mean difference"
  )
}

comparison <- bind_rows(lapply(cat_vars, compare_cat),
                        lapply(cont_vars, compare_cont))

# For continuous variables, the "pct" columns contain weighted means and
# diff_pp contains the difference in means.
write.csv(comparison, file.path(out_dir, "stage0_sample_comparison.csv"),
          row.names = FALSE, na = "")
print(comparison)


# ---- 0.9 Save for later stages ----------------------------------------------

saveRDS(dat, file.path(out_dir, "rls_stage0.rds"))

sessionInfo()