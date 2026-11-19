# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 9.0 Packages, paths, data ----------------------------------------------

library(dplyr)
library(survey)

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")

dat  <- readRDS(file.path(out_dir, "rls_stage0.rds"))
post <- readRDS(file.path(out_dir, "stage4_posteriors.rds"))

pp <- as.data.frame(post$A5)
names(pp)[names(pp) == "modal"] <- "profile"
pp$respid  <- as.numeric(post$respid)
dat$respid <- as.numeric(dat$respid)
dat <- left_join(dat, pp, by = "respid")

# Spouse has a different religion (as in Stage 7)
dat$sp_diffrel <- factor(as.numeric(dat$partner_rel == "Spouse/partner different religion"),
                         levels = 0:1, labels = c("No", "Yes"))

lp    <- dat %>% filter(lpa_sample)
k     <- 5
pcols <- paste0("p", 1:k)
reps  <- grep("^repwt_[0-9]+$", names(lp), value = TRUE)
P     <- as.matrix(lp[, pcols])
s     <- lp$profile
w     <- lp$weight

cat_vars <- c("rel_curr9", "rel_raised9", "agegrp", "woman", "educrec", "income4",
              "ideo", "partyid3", "usgen", "sp_relimp", "sp_diffrel")
num_vars <- c("child_relig", "ch_attend", "ch_relschool", "ch_religed")

# Categorical covariates stored as numbers are turned into factors
for (v in cat_vars) if (!is.factor(lp[[v]])) lp[[v]] <- factor(lp[[v]])
lp$prof <- factor(s, levels = 1:k)

des <- svrepdesign(data = lp, weights = ~weight, repweights = "repwt_[0-9]+",
                   type = "bootstrap", combined.weights = TRUE, mse = TRUE)

cramer_v <- function(tab) {
  p <- tab / sum(tab)
  e <- outer(rowSums(p), colSums(p))
  sqrt(sum((p - e)^2 / e) / (min(dim(p)) - 1))
}


# ---- 9.1 Bivariate associations with profile membership ---------------------

# (a) Categorical covariates
col_pct <- list()
cat_res <- do.call(rbind, lapply(cat_vars, function(v) {
  x  <- lp[[v]]
  ok <- !is.na(x)
  xo <- droplevels(x[ok])
  tab_prop  <- rowsum(P[ok, ] * w[ok], xo)                       # categories x profiles
  tab_modal <- xtabs(w[ok] ~ xo + factor(s[ok], levels = 1:k))
  col_pct[[v]] <<- data.frame(variable = v, category = rownames(tab_prop),
                              round(100 * sweep(tab_prop, 2, colSums(tab_prop), "/"), 1),
                              row.names = NULL)
  rs <- svychisq(as.formula(paste0("~", v, " + prof")), des, statistic = "F")
  data.frame(variable = v, n = sum(ok), categories = nlevels(xo),
             F = round(unname(rs$statistic), 1),
             df1 = round(unname(rs$parameter[1]), 2),
             df2 = round(unname(rs$parameter[2]), 1),
             p = unname(rs$p.value),
             V_proportional = round(cramer_v(tab_prop), 3),
             V_modal = round(cramer_v(tab_modal), 3),
             row.names = NULL)
}))

# (b) Numeric covariates (upbringing scores)
num_res <- do.call(rbind, lapply(num_vars, function(v) {
  y  <- lp[[v]]
  ok <- !is.na(y)
  wp <- P[ok, ] * w[ok]
  m_prof <- colSums(wp * y[ok]) / colSums(wp)
  gm  <- weighted.mean(y[ok], w[ok])
  tss <- sum(w[ok] * (y[ok] - gm)^2)
  bss <- sum(colSums(wp) * (m_prof - gm)^2)
  fit <- svyglm(as.formula(paste(v, "~ prof")), des[ok, ])
  tt  <- regTermTest(fit, ~prof)
  data.frame(variable = v, n = sum(ok),
             F = round(as.numeric(tt$Ftest), 1), df1 = tt$df, df2 = tt$ddf,
             p = as.numeric(tt$p),
             eta_proportional = round(sqrt(bss / tss), 3),
             t(setNames(round(m_prof, 2), paste0("mean_p", 1:k))),
             row.names = NULL)
}))

# Bonferroni correction over all tests in 9.1
m_tests <- nrow(cat_res) + nrow(num_res)
cat_res$p_bonferroni <- pmin(1, cat_res$p * m_tests)
num_res$p_bonferroni <- pmin(1, num_res$p * m_tests)
cat_res[, c("p", "p_bonferroni")] <- signif(cat_res[, c("p", "p_bonferroni")], 3)
num_res[, c("p", "p_bonferroni")] <- signif(num_res[, c("p", "p_bonferroni")], 3)

cat_res <- cat_res[order(-cat_res$V_proportional), ]
num_res <- num_res[order(-num_res$eta_proportional), ]

write.csv(cat_res, file.path(out_dir, "stage9_bivariate_categorical.csv"), row.names = FALSE)
write.csv(num_res, file.path(out_dir, "stage9_bivariate_numeric.csv"), row.names = FALSE)
cat("\n==== Categorical covariates by profile, ranked by Cramer's V ====\n")
print(cat_res, row.names = FALSE)
cat("\n==== Upbringing scores by profile (weighted means), ranked by eta ====\n")
print(num_res, row.names = FALSE)

# (c) Distribution of each categorical covariate within each profile
col_pct <- do.call(rbind, col_pct)
names(col_pct)[3:(k + 2)] <- paste0("profile_", 1:k)
write.csv(col_pct, file.path(out_dir, "stage9_covariates_within_profiles.csv"), row.names = FALSE)
cat("\n==== Distribution of each covariate within each profile (column %) ====\n")
print(col_pct, row.names = FALSE)


# ---- 9.2 Incremental validity: wellbeing ------------------------------------
# Three weighted least-squares models of wellbeing, all with the same
# demographic controls (age group, gender, education, income) and the same
# respondents:
#   (a) tradition: current religious tradition
#   (b) profile:   profile membership (BCH-corrected, as in Stage 8)
#   (c) both
# R-squared of each, and the gain from adding profile to tradition (c - a)
# and tradition to profile (c - b), with 95% confidence intervals from the
# 500 replicate weights.

covs <- c("agegrp", "woman", "educrec", "income4")
cc   <- complete.cases(lp[, c("wellbeing", "rel_curr9", covs)])

D <- t(sapply(1:k, function(c) sapply(1:k, function(j) sum(w * P[, c] * (s == j)) / sum(w * P[, c]))))
H <- solve(D)[s, , drop = FALSE]

ex <- lp[rep(which(cc), each = k), c("wellbeing", "rel_curr9", covs)]
ex$prof <- factor(rep(1:k, times = sum(cc)), levels = c(2, 1, 3, 4, 5))
ex   <- droplevels(ex)
hvec <- as.vector(t(H[cc, , drop = FALSE]))
y    <- ex$wellbeing
Xa   <- model.matrix(reformulate(c("rel_curr9", covs)), data = ex)
Xb   <- model.matrix(reformulate(c("prof", covs)), data = ex)
Xc   <- model.matrix(reformulate(c("prof", "rel_curr9", covs)), data = ex)
yr   <- lp$wellbeing[cc]

r2 <- function(X, wt_resp) {
  wv  <- rep(wt_resp[cc], each = k) * hvec
  b   <- solve(crossprod(X, X * wv), crossprod(X, wv * y))
  ssr <- sum(wv * (y - X %*% b)^2)
  wr  <- wt_resp[cc]
  tss <- sum(wr * (yr - weighted.mean(yr, wr))^2)
  1 - ssr / tss
}
r2_all <- function(wt) {
  a <- r2(Xa, wt); b <- r2(Xb, wt); c <- r2(Xc, wt)
  c(tradition = a, profile = b, both = c,
    profile_beyond_tradition = c - a, tradition_beyond_profile = c - b)
}

r2_est <- r2_all(w)
r2_rep <- sapply(reps, function(r) r2_all(lp[[r]]))
r2_se  <- sqrt(rowSums((r2_rep - r2_est)^2) / (length(reps) - 1))

r2_tab <- data.frame(model = names(r2_est),
                     r_squared = round(r2_est, 4),
                     ci_low = round(r2_est - 1.96 * r2_se, 4),
                     ci_high = round(r2_est + 1.96 * r2_se, 4),
                     row.names = NULL)
write.csv(r2_tab, file.path(out_dir, "stage9_incremental_validity_wellbeing.csv"), row.names = FALSE)
cat("\n==== Wellbeing: R-squared with tradition, profile, and both (demographics in all) ====\n")
cat("Respondents:", sum(cc), "\n")
print(r2_tab, row.names = FALSE)


# ---- 9.3 Session information ------------------------------------------------

sessionInfo()