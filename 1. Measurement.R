# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 1.0 Packages, paths, data ----------------------------------------------

pkgs <- c("lavaan", "semTools", "psych", "dplyr")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(lavaan)
library(semTools)
library(psych)
library(dplyr)

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")

dat <- readRDS(file.path(out_dir, "rls_stage0.rds"))
lpa <- dat %>% filter(lpa_sample)          # analytic sample, n = 31,470
stopifnot(nrow(lpa) == 31470)

h1 <- lpa %>% filter(half == 1)            # n = 15,805 (EFA)
h2 <- lpa %>% filter(half == 2)            # n = 15,665 (CFA)

rb_items <- c("attndperrlsr", "attndonrlsr", "synmemr", "prayfreqr", "gracefreqr",
              "prac_ar", "prac_br", "prac_cr", "prac_dr")
rl_items <- c("relimpr", "relperr", "gdcertainr", "afterlife")
sb_items <- c("spiract_ar", "spiract_br", "spiract_cr", "spiract_dr",
              "spiract_fr")
sl_items <- c("spirperr", "soulr", "spirworld2", "exp_ar", "exp_dr")
items    <- c(rb_items, rl_items, sb_items, sl_items)


# ---- 1.1 Item distributions (analytic sample) -------------------------------

item_dist <- do.call(rbind, lapply(items, function(i) {
  tb <- table(lpa[[i]])
  data.frame(item = i, category = names(tb), n = as.vector(tb),
             pct = round(100 * as.vector(prop.table(tb)), 1))
}))
write.csv(item_dist, file.path(out_dir, "stage1_item_distributions.csv"), row.names = FALSE)


# ---- 1.2 Number of factors (half 1) -----------------------------------------
# Polychoric correlation matrix, its eigenvalues, and parallel analysis.
# Parallel analysis compares the observed eigenvalues with those from random
# data of the same size; factors are retained while observed > random.

pc_h1 <- lavCor(h1[items], ordered = items)

eig <- data.frame(factor = seq_along(items), eigenvalue = round(eigen(pc_h1)$values, 3))
write.csv(eig, file.path(out_dir, "stage1_eigenvalues.csv"), row.names = FALSE)
print(head(eig, 8))

set.seed(20261016)
pa <- fa.parallel(pc_h1, n.obs = nrow(h1), fa = "fa", fm = "minres",
                  n.iter = 50, plot = TRUE)
print(pa$nfact)


# ---- 1.3 Exploratory factor analysis (half 1) -------------------------------
# 3 to 6 factors, oblique (oblimin) rotation, WLSMV.

efa_syntax <- function(k) {
  paste0(paste0('efa("efa")*f', 1:k, collapse = " + "),
         " =~ ", paste(items, collapse = " + "))
}

fit_measures <- c("chisq.scaled", "df.scaled", "pvalue.scaled",
                  "cfi.scaled", "tli.scaled", "rmsea.scaled",
                  "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled",
                  "cfi.robust", "tli.robust", "rmsea.robust", "srmr")

fm_row <- function(fit, label) {
  data.frame(model = label, as.list(round(unclass(fitMeasures(fit, fit_measures)), 3)))
}

efa_fits <- lapply(3:6, function(k) {
  cfa(efa_syntax(k), data = h1, ordered = items,
      estimator = "WLSMV", rotation = "oblimin")
})
names(efa_fits) <- paste0("efa_", 3:6)

efa_fit_table <- do.call(rbind, lapply(names(efa_fits), function(n) fm_row(efa_fits[[n]], n)))
write.csv(efa_fit_table, file.path(out_dir, "stage1_efa_fit.csv"), row.names = FALSE)
print(efa_fit_table)

# Standardised loadings and factor correlations for each solution
for (n in names(efa_fits)) {
  std <- lavInspect(efa_fits[[n]], "std")
  load <- round(std$lambda, 3)
  write.csv(load, file.path(out_dir, paste0("stage1_", n, "_loadings.csv")))
  write.csv(round(std$psi, 3), file.path(out_dir, paste0("stage1_", n, "_factor_cor.csv")))
  cat("\n====", n, "loadings (|loading| >= .30 shown) ====\n")
  load_print <- load
  load_print[abs(load_print) < .30] <- NA
  print(load_print, na.print = "")
}


# ---- 1.4 Confirmatory factor analysis (half 2) ------------------------------

cfa_model <- paste0(
  "RelBehav   =~ ", paste(rb_items, collapse = " + "), "\n",
  "RelBelief  =~ ", paste(rl_items, collapse = " + "), "\n",
  "SpirBehav  =~ ", paste(sb_items, collapse = " + "), "\n",
  "SpirBelief =~ ", paste(sl_items, collapse = " + ")
)

fit_cfa <- cfa(cfa_model, data = h2, ordered = items, estimator = "WLSMV")

cfa_fit_table <- fm_row(fit_cfa, "cfa_4factor")
write.csv(cfa_fit_table, file.path(out_dir, "stage1_cfa_fit.csv"), row.names = FALSE)
print(cfa_fit_table)

cfa_std <- standardizedSolution(fit_cfa) %>%
  filter(op == "=~" | (op == "~~" & lhs != rhs & lhs %in% c("RelBehav", "RelBelief", "SpirBehav", "SpirBelief")))
write.csv(cfa_std, file.path(out_dir, "stage1_cfa_std_solution.csv"), row.names = FALSE)
print(cfa_std[, c("lhs", "op", "rhs", "est.std", "se", "ci.lower", "ci.upper")])

# Largest modification indices (for diagnosis only; the model is not modified)
mi <- modindices(fit_cfa, sort. = TRUE, maximum.number = 20)
write.csv(mi, file.path(out_dir, "stage1_cfa_modindices.csv"), row.names = FALSE)
print(mi[, c("lhs", "op", "rhs", "mi", "sepc.all")])


# ---- 1.4b Alternative structures suggested by the half-1 EFA (half 2) -------
# Decision rules, fixed before fitting:
#  (a) best robust CFI, TLI, RMSEA and SRMR;
#  (b) every factor correlation has an upper 95% CI limit below .85;
#  (c) every item loads at least .50 on its own factor.

pi_items <- c(rb_items, "relimpr", "relperr")                 # practice and identity
id_items <- c("relimpr", "relperr")                           # religious identity
bel7     <- c("gdcertainr", "afterlife", "soulr", "spirworld2", "spirperr",
              "exp_ar", "exp_dr")                             # belief, with experience
bel5     <- c("gdcertainr", "afterlife", "soulr", "spirworld2", "spirperr")
ex_items <- c("exp_ar", "exp_dr")                             # spiritual experience

mk <- function(...) {
  f <- list(...)
  paste(paste(names(f), "=~", sapply(f, paste, collapse = " + ")), collapse = "\n")
}

alt_models <- list(
  m4_original   = cfa_model,
  m3            = mk(PracID = pi_items, Belief = bel7, SpirPrac = sb_items),
  m4_identity   = mk(RelPrac = rb_items, RelID = id_items, Belief = bel7, SpirPrac = sb_items),
  m4_experience = mk(PracID = pi_items, Belief = bel5, SpirPrac = sb_items, SpirExp = ex_items)
)

alt_fits <- lapply(alt_models, function(m)
  cfa(m, data = h2, ordered = items, estimator = "WLSMV"))

alt_fit_table <- do.call(rbind, lapply(names(alt_fits), function(n) fm_row(alt_fits[[n]], n)))
write.csv(alt_fit_table, file.path(out_dir, "stage1_alt_cfa_fit.csv"), row.names = FALSE)
print(alt_fit_table)

for (n in names(alt_fits)) {
  s  <- standardizedSolution(alt_fits[[n]])
  lv <- lavNames(alt_fits[[n]], "lv")
  cat("\n====", n, ": factor correlations ====\n")
  print(s[s$op == "~~" & s$lhs != s$rhs & s$lhs %in% lv,
          c("lhs", "rhs", "est.std", "ci.lower", "ci.upper")], row.names = FALSE)
  cat("\n====", n, ": loadings below .50 ====\n")
  print(s[s$op == "=~" & s$est.std < .50, c("lhs", "rhs", "est.std")], row.names = FALSE)
  cat("\n====", n, ": 10 largest modification indices ====\n")
  print(modindices(alt_fits[[n]], sort. = TRUE, maximum.number = 10)[, c("lhs", "op", "rhs", "mi", "sepc.all")],
        row.names = FALSE)
}


# ---- 1.5 Reliability (half 2) -----------------------------------------------
# Omega: reliability of each raw sum score, taking into account that the items
# have a small number of answer categories (Green & Yang, 2009).
# Cronbach's alpha: on the raw items, as in the Stata file, same half sample.

omega_ord <- compRelSEM(fit_cfa, ord.scale = TRUE)

alpha_raw <- sapply(list(RelBehav = rb_items, RelBelief = rl_items,
                         SpirBehav = sb_items, SpirBelief = sl_items),
                    function(v) psych::alpha(h2[v], warnings = FALSE)$total$raw_alpha)

reliability <- data.frame(factor         = names(omega_ord),
                          omega          = round(unlist(omega_ord), 3),
                          cronbach_alpha = round(alpha_raw[names(omega_ord)], 3),
                          row.names = NULL)
write.csv(reliability, file.path(out_dir, "stage1_reliability.csv"), row.names = FALSE)
print(reliability)



# ---- 1.5b Reliability in the full analytic sample (n = 31,470) --------------
# Same measures as 1.5 (ordinal omega and Cronbach's alpha), calculated on the
# whole analytic sample, which is the sample the sum scores are used on.

fit_cfa_all <- cfa(cfa_model, data = lpa, ordered = items, estimator = "WLSMV")
omega_all   <- compRelSEM(fit_cfa_all, ord.scale = TRUE)

alpha_all <- sapply(list(RelBehav = rb_items, RelBelief = rl_items,
                         SpirBehav = sb_items, SpirBelief = sl_items),
                    function(v) psych::alpha(lpa[v], warnings = FALSE)$total$raw_alpha)

rel_all <- data.frame(factor         = names(omega_all),
                      n              = nrow(lpa),
                      omega          = round(unlist(omega_all), 3),
                      cronbach_alpha = round(alpha_all[names(omega_all)], 3),
                      row.names = NULL)
write.csv(rel_all, file.path(out_dir, "stage1_reliability_full_sample.csv"), row.names = FALSE)
print(rel_all)


# ---- 1.5c Reliability of the wellbeing index (analytic sample) --------------
# Cronbach's alpha of the three wellbeing items (general happiness, self-rated
# health, quality of family life), for respondents in the analytic sample with
# all three items.

wb_items <- c("wb_happy", "wb_health", "wb_family")
wb_ok    <- complete.cases(lpa[, wb_items])
wb_alpha <- psych::alpha(lpa[wb_ok, wb_items], warnings = FALSE)$total$raw_alpha
cat("Wellbeing: Cronbach's alpha =", round(wb_alpha, 3), " n =", sum(wb_ok), "\n")


# ---- 1.6 Measurement invariance across nine religious traditions -------------
# Full analytic sample, groups = rel_curr9 (respondents with missing tradition
# are excluded from this step only).
#
# Step A: check that every response category of every item occurs in every
# group. An ordinal multi-group model cannot be estimated if a category is
# empty in any group. If empty cells are found, the script stops here and
# reports them, so that a decision can be made explicitly.

inv <- lpa %>% filter(!is.na(rel_curr9)) %>% mutate(tradition = droplevels(rel_curr9))
print(table(inv$tradition))

cell_counts <- do.call(rbind, lapply(items, function(i) {
  tb <- table(factor(inv[[i]]), inv$tradition)
  data.frame(item = i, category = rep(rownames(tb), ncol(tb)),
             group = rep(colnames(tb), each = nrow(tb)),
             n = as.vector(tb))
}))
write.csv(cell_counts, file.path(out_dir, "stage1_invariance_cell_counts.csv"), row.names = FALSE)

empty_cells <- cell_counts %>% filter(n == 0)
sparse_cells <- cell_counts %>% filter(n > 0 & n < 10)
cat("\nEmpty cells:\n");  print(empty_cells)
cat("\nSparse cells (1 to 9 cases):\n"); print(sparse_cells)

if (nrow(empty_cells) > 0) {
  stop("Empty response categories in at least one tradition (listed above). ",
       "Stop here and decide how to handle them before running the invariance models.")
}

# Step B: invariance sequence for ordinal indicators (Wu & Estabrook, 2016):
#   configural -> equal thresholds -> + equal loadings -> + equal intercepts
# Changes in fit are judged by delta CFI <= .010 and delta RMSEA <= .015
# (Chen, 2007). With N = 31,000, chi-square difference tests will be
# significant for trivial differences and are reported for completeness only.

fit_mi <- function(eq) {
  measEq.syntax(configural.model = cfa_model, data = inv, ordered = items,
                parameterization = "theta", ID.fac = "std.lv",
                ID.cat = "Wu.Estabrook.2016", group = "tradition",
                group.equal = eq, estimator = "WLSMV", return.fit = TRUE)
}

fit_config <- fit_mi("")
fit_thresh <- fit_mi("thresholds")
fit_metric <- fit_mi(c("thresholds", "loadings"))
fit_scalar <- fit_mi(c("thresholds", "loadings", "intercepts"))

mi_list <- list(configural = fit_config, thresholds = fit_thresh,
                metric = fit_metric, scalar = fit_scalar)

mi_table <- do.call(rbind, lapply(names(mi_list), function(n) fm_row(mi_list[[n]], n)))
mi_table$delta_cfi   <- c(NA, diff(mi_table$cfi.scaled))
mi_table$delta_rmsea <- c(NA, diff(mi_table$rmsea.scaled))
write.csv(mi_table, file.path(out_dir, "stage1_invariance_fit.csv"), row.names = FALSE)
print(mi_table)

print(lavTestLRT(fit_config, fit_thresh, fit_metric, fit_scalar))


# ---- 1.7 Save ---------------------------------------------------------------

saveRDS(list(efa = efa_fits, cfa = fit_cfa, invariance = mi_list),
        file.path(out_dir, "stage1_fits.rds"))

sessionInfo()