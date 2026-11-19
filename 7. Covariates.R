# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 7.0 Packages, paths, data ----------------------------------------------

library(dplyr)

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")

dat  <- readRDS(file.path(out_dir, "rls_stage0.rds"))
post <- readRDS(file.path(out_dir, "stage4_posteriors.rds"))

pp <- as.data.frame(post$A5)
names(pp)[names(pp) == "modal"] <- "profile"
pp$respid <- as.numeric(post$respid)
dat$respid <- as.numeric(dat$respid)
dat <- left_join(dat, pp, by = "respid")

k      <- 5
pcols  <- paste0("p", 1:k)
ref    <- 2                                   # reference profile
reps   <- grep("^repwt_[0-9]+$", names(dat), value = TRUE)

# Spouse has a different religion (0 = same religion or no spouse/partner)
dat$sp_diffrel <- as.numeric(dat$partner_rel == "Spouse/partner different religion")

blocks <- list(
  M1 = c("agegrp", "woman", "educrec", "income4", "ideo", "partyid3", "usgen"),
  M2 = c("rel_raised9", "child_relig", "ch_attend", "ch_relschool", "ch_religed"),
  M3 = c("sp_relimp", "sp_diffrel")
)
vars <- list(M1 = blocks$M1,
             M2 = c(blocks$M1, blocks$M2),
             M3 = c(blocks$M1, blocks$M2, blocks$M3))


# ---- 7.1 Analytic sample for the covariate models ---------------------------

cv <- dat %>% filter(lpa_sample)
cv <- cv[complete.cases(cv[, vars$M3]), ]
cv <- droplevels(cv)
cat("Respondents in the covariate models:", nrow(cv),
    "of", sum(dat$lpa_sample), "in the profile analysis\n")

# Sparse cells: estimated number of respondents (sum of profile
# probabilities) in each category x profile. Results for categories with
# fewer than 20 estimated members in a profile are unreliable and should be
# reported with that caution.
fac_cov <- vars$M3[sapply(cv[vars$M3], is.factor)]
sparse <- do.call(rbind, lapply(fac_cov, function(v) {
  tab <- rowsum(as.matrix(cv[, paste0("p", 1:5)]), cv[[v]])
  data.frame(variable = v, category = rownames(tab), round(tab, 1), row.names = NULL)
}))
write.csv(sparse, file.path(out_dir, "stage7_category_by_profile.csv"), row.names = FALSE)
cat("\nCategories with fewer than 20 estimated members in at least one profile:\n")
print(sparse[apply(sparse[, paste0("p", 1:5)], 1, min) < 20, ], row.names = FALSE)


# ---- 7.2 Classification error matrix ----------------------------------------
# D[c, s] = weighted share of true profile c that is assigned to profile s.

P  <- as.matrix(cv[, pcols])
s  <- cv$profile
w  <- cv$weight
D  <- t(sapply(1:k, function(c) sapply(1:k, function(j) sum(w * P[, c] * (s == j)) / sum(w * P[, c]))))
dimnames(D) <- list(paste0("true_", 1:k), paste0("assigned_", 1:k))
write.csv(round(D, 4), file.path(out_dir, "stage7_classification_matrix.csv"))
cat("\nShare of each true profile assigned to each profile:\n"); print(round(D, 3))


# ---- 7.3 Weighted multinomial logistic regression, ML three-step ----------
# Likelihood for respondent i, assigned to profile s_i:
#   L_i = sum over profiles c of P(c | covariates_i) * D[c, s_i]
# i.e. the model for true profile membership, combined with the known
# probability of being assigned to s_i when the true profile is c. Weighted by
# the survey weight (rescaled to average 1). Fitted by Newton-Raphson with
# step halving.
# A very small ridge penalty (lambda = 1, equivalent to a weak prior with
# standard deviation 1 on each log relative risk ratio, intercepts excluded)
# keeps estimates finite when a small category has no members of a profile.
# With 27,000 respondents it changes well-estimated coefficients only in the
# third decimal; it is reported in the Methods section.
lambda <- 1

fit_ml3 <- function(X, wt, start = NULL, tol = 1e-7, maxit = 200) {
  nonref <- setdiff(1:k, ref)
  p  <- ncol(X); q <- length(nonref)
  wt <- wt / mean(wt)
  pen <- rep(as.numeric(colnames(X) != "(Intercept)"), q)   # penalised terms
  b  <- if (is.null(start)) rep(0, p * q) else start
  Ds <- t(D)[s, , drop = FALSE]               # n x k: D[c, s_i] for each c
  ev <- function(b) {
    eta <- matrix(0, nrow(X), k)
    eta[, nonref] <- X %*% matrix(b, ncol = q)
    m   <- apply(eta, 1, max)
    Pr  <- exp(eta - m); Pr <- Pr / rowSums(Pr)
    L   <- rowSums(Pr * Ds)
    A   <- Pr * Ds / L                        # probability of true profile c
    list(ll = sum(wt * log(L)) - 0.5 * lambda * sum(pen * b^2), P = Pr, A = A)
  }
  cur <- ev(b); conv <- 1
  for (it in seq_len(maxit)) {
    P <- cur$P; A <- cur$A
    G <- as.vector(crossprod(X, wt * (A[, nonref] - P[, nonref]))) - lambda * pen * b
    # Curvature: exact second derivative, or the multinomial information
    # matrix if the exact one is not usable at the current point
    build_H <- function(exact) {
      H <- matrix(0, p * q, p * q)
      for (aa in 1:q) for (cc in aa:q) {
        ja <- nonref[aa]; jc <- nonref[cc]
        v <- P[, ja] * ((aa == cc) - P[, jc])
        if (exact) v <- v - A[, ja] * ((aa == cc) - A[, jc])
        blk <- crossprod(X, X * (wt * v))
        H[(aa - 1) * p + 1:p, (cc - 1) * p + 1:p] <- blk
        H[(cc - 1) * p + 1:p, (aa - 1) * p + 1:p] <- t(blk)
      }
      H + diag(lambda * pen, p * q)
    }
    step <- tryCatch(solve(build_H(TRUE), G), error = function(e) NULL)
    if (is.null(step) || sum(step * G) <= 0) step <- solve(build_H(FALSE), G)
    t_ <- 1
    repeat {
      new <- ev(b + t_ * step)
      if (new$ll >= cur$ll - 1e-10 || t_ < 1e-4) break
      t_ <- t_ / 2
    }
    b <- b + t_ * step; cur <- new
    if (max(abs(t_ * step)) < tol) { conv <- 0; break }
  }
  list(b = b, ll = cur$ll + 0.5 * lambda * sum(pen * b^2), conv = conv)
}

run_model <- function(mname) {
  f <- reformulate(vars[[mname]])
  X <- model.matrix(f, data = cv)
  t0  <- Sys.time()
  full <- fit_ml3(X, w)
  if (full$conv != 0) warning(mname, ": full-sample estimation did not converge")
  
  # Intercept-only model, for a descriptive pseudo R-squared
  null <- fit_ml3(X[, 1, drop = FALSE], w)
  
  # Re-estimate with each replicate weight (classification matrix kept fixed)
  rep_b <- sapply(reps, function(r) fit_ml3(X, cv[[r]], start = full$b)$b)
  # Bootstrap replicate variance, as in survey::svrepdesign(type = "bootstrap",
  # mse = TRUE): squared deviations from the full-sample estimate / (R - 1)
  se <- sqrt(rowSums((rep_b - full$b)^2) / (length(reps) - 1))
  
  nonref <- setdiff(1:k, ref)
  out <- data.frame(term    = rep(colnames(X), times = k - 1),
                    profile = rep(nonref, each = ncol(X)),
                    b = full$b, se = se)
  out$rrr      <- exp(out$b)
  out$ci_low   <- exp(out$b - 1.96 * out$se)
  out$ci_high  <- exp(out$b + 1.96 * out$se)
  out$p        <- 2 * pnorm(-abs(out$b / out$se))
  out <- out[out$term != "(Intercept)", ]
  out[, c("rrr", "ci_low", "ci_high")] <- round(out[, c("rrr", "ci_low", "ci_high")], 2)
  out$p <- signif(out$p, 3)
  
  message(mname, ": done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
  list(table = out,
       fit = data.frame(model = mname, n = nrow(cv), parameters = length(full$b),
                        loglik_weighted = round(full$ll, 1),
                        pseudo_r2 = round(1 - full$ll / null$ll, 3)))
}

res <- lapply(names(vars), run_model)
names(res) <- names(vars)


# ---- 7.4 Results tables -----------------------------------------------------
# One table per model. Rows: covariate category; one block of rows per
# profile (compared with profile 2). Wide versions show the relative risk
# ratio and its confidence interval for each profile side by side.

for (mname in names(res)) {
  tab <- res[[mname]]$table
  write.csv(tab, file.path(out_dir, paste0("stage7_rrr_", mname, ".csv")), row.names = FALSE)
  tab$cell <- sprintf("%.2f (%.2f-%.2f)", tab$rrr, tab$ci_low, tab$ci_high)
  wide <- reshape(tab[, c("term", "profile", "cell")], idvar = "term",
                  timevar = "profile", direction = "wide")
  names(wide) <- sub("^cell\\.", "profile_", names(wide))
  cat("\n==== ", mname, ": relative risk ratios (95% CI), reference = profile ", ref, " ====\n", sep = "")
  print(wide, row.names = FALSE)
}

fit_tab <- do.call(rbind, lapply(res, `[[`, "fit"))
write.csv(fit_tab, file.path(out_dir, "stage7_model_fit.csv"), row.names = FALSE)
print(fit_tab, row.names = FALSE)


# ---- 7.5 Collinearity check and average marginal effects (model M3) ---------
# (a) Collinearity among the covariates of the full model (M3), with special
#     attention to the upbringing block (raised-in tradition, childhood
#     religiosity, childhood attendance, religious schooling, religious
#     education). Shown as (i) weighted correlations among the four
#     upbringing scores and (ii) generalised variance inflation factors
#     (GVIF; Fox & Monette, 1992) for every covariate. GVIF^(1/(2*df)) above
#     about 2.2 (equivalent to a VIF of 5) signals a collinearity problem.
# (b) Average predicted probabilities: for each category of each covariate,
#     the predicted percentage of respondents in each profile if everyone in
#     the sample had that category, keeping their other characteristics as
#     they are (average marginal effects in percentage points). For 0/1
#     indicators, everyone at 0 versus everyone at 1; for other numeric
#     covariates (scores), the effect of a one-unit increase. Weighted; 95% confidence intervals from
#     the 500 replicate weights (model M3 re-estimated 500 times; this takes
#     about as long as M3 did in section 7.3). Set ame_se to FALSE to skip the
#     confidence intervals and get the estimates in a few seconds.

ame_se <- TRUE

# (a)(i) Weighted correlations among the upbringing scores
up_num <- c("child_relig", "ch_attend", "ch_relschool", "ch_religed")
up_cor <- round(cov.wt(as.matrix(cv[, up_num]), wt = w, cor = TRUE)$cor, 3)
write.csv(up_cor, file.path(out_dir, "stage7_upbringing_correlations.csv"))
cat("\n==== Weighted correlations among the upbringing scores ====\n")
print(up_cor)

# (a)(ii) Generalised variance inflation factors for all M3 covariates.
# GVIF depends only on the covariates, so any outcome can be used in the
# auxiliary linear model; religious behaviour is used here.
gv <- car::vif(lm(reformulate(vars$M3, response = "religbehav"), data = cv, weights = w))
gv <- as.data.frame(gv)
gv$covariate <- rownames(gv)
names(gv) <- c("GVIF", "df", "GVIF_adjusted", "covariate")
gv <- gv[, c("covariate", "df", "GVIF", "GVIF_adjusted")]
gv[, c("GVIF", "GVIF_adjusted")] <- round(gv[, c("GVIF", "GVIF_adjusted")], 2)
write.csv(gv, file.path(out_dir, "stage7_gvif.csv"), row.names = FALSE)
cat("\n==== Generalised variance inflation factors (model M3) ====\n")
cat("GVIF_adjusted = GVIF^(1/(2*df)); values above about 2.2 signal a problem\n")
print(gv, row.names = FALSE)

# (b) Average predicted probabilities
f3   <- reformulate(vars$M3)
X3   <- model.matrix(f3, data = cv)
fit3 <- fit_ml3(X3, w)
rep3 <- if (ame_se) sapply(reps, function(r) fit_ml3(X3, cv[[r]], start = fit3$b)$b) else NULL

nonref <- setdiff(1:k, ref)
pred_P <- function(b, X) {
  eta <- matrix(0, nrow(X), k)
  eta[, nonref] <- X %*% matrix(b, ncol = k - 1)
  m <- apply(eta, 1, max)
  Pr <- exp(eta - m); Pr / rowSums(Pr)
}
avg_P <- function(b, X, wt) colSums(pred_P(b, X) * wt) / sum(wt)   # length k

# Counterfactual design matrices: each category of each factor covariate,
# and each numeric covariate at its observed value and +1
scen <- list()
for (v in vars$M3) {
  if (is.factor(cv[[v]])) {
    for (lv in levels(cv[[v]])) {
      nd <- cv; nd[[v]] <- factor(lv, levels = levels(cv[[v]]))
      scen[[length(scen) + 1]] <- list(variable = v, level = lv,
                                       is_ref = lv == levels(cv[[v]])[1],
                                       X = model.matrix(f3, data = nd))
    }
  } else if (all(cv[[v]] %in% c(0, 1))) {
    # 0/1 indicator: everyone at 0 (reference) versus everyone at 1
    for (val in 0:1) {
      nd <- cv; nd[[v]] <- val
      scen[[length(scen) + 1]] <- list(variable = v, level = as.character(val),
                                       is_ref = val == 0, X = model.matrix(f3, data = nd))
    }
  } else {
    nd0 <- cv; nd1 <- cv; nd1[[v]] <- nd1[[v]] + 1
    scen[[length(scen) + 1]] <- list(variable = v, level = "observed", is_ref = TRUE,
                                     X = model.matrix(f3, data = nd0))
    scen[[length(scen) + 1]] <- list(variable = v, level = "+1", is_ref = FALSE,
                                     X = model.matrix(f3, data = nd1))
  }
}

ame <- do.call(rbind, lapply(unique(sapply(scen, `[[`, "variable")), function(v) {
  sv  <- scen[sapply(scen, `[[`, "variable") == v]
  ref_s <- sv[[which(sapply(sv, `[[`, "is_ref"))[1]]]
  do.call(rbind, lapply(sv, function(sc) {
    est  <- avg_P(fit3$b, sc$X, w)
    dref <- est - avg_P(fit3$b, ref_s$X, w)
    if (ame_se) {
      r_est  <- sapply(seq_along(reps), function(i) avg_P(rep3[, i], sc$X, cv[[reps[i]]]))
      r_dref <- r_est - sapply(seq_along(reps), function(i) avg_P(rep3[, i], ref_s$X, cv[[reps[i]]]))
      se_est <- sqrt(rowSums((r_est - est)^2) / (length(reps) - 1))
      se_d   <- sqrt(rowSums((r_dref - dref)^2) / (length(reps) - 1))
    } else { se_est <- se_d <- rep(NA, k) }
    data.frame(variable = v, level = sc$level, profile = 1:k,
               pred_pct = round(100 * est, 1),
               pred_ci_low = round(100 * (est - 1.96 * se_est), 1),
               pred_ci_high = round(100 * (est + 1.96 * se_est), 1),
               diff_pp = round(100 * dref, 1),
               diff_ci_low = round(100 * (dref - 1.96 * se_d), 1),
               diff_ci_high = round(100 * (dref + 1.96 * se_d), 1),
               row.names = NULL)
  }))
}))
write.csv(ame, file.path(out_dir, "stage7_average_marginal_effects_M3.csv"), row.names = FALSE)

# Readable layout: predicted % in each profile, one row per category.
# For factor covariates the first category is the reference; diff_pp in the
# saved file gives each category's difference from it, with its CI.
ame_wide <- reshape(ame[, c("variable", "level", "profile", "pred_pct")],
                    idvar = c("variable", "level"), timevar = "profile", direction = "wide")
names(ame_wide) <- sub("^pred_pct\\.", "profile_", names(ame_wide))
cat("\n==== Average predicted percentage in each profile (model M3) ====\n")
print(ame_wide, row.names = FALSE)


# ---- 7.6 Session information ------------------------------------------------

sessionInfo()