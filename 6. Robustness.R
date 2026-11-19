# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 6.0 Packages, paths, data, helper functions ----------------------------

library(OpenMx)
library(tidySEM)
library(dplyr)

mxOption(NULL, "Number of Threads", max(1, parallel::detectCores() - 1))
mxOption(NULL, "Calculate Hessian", "No")
mxOption(NULL, "Standard Errors", "No")

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")
rob_dir  <- file.path(out_dir, "stage6_models")
dir.create(rob_dir, showWarnings = FALSE)

dat <- readRDS(file.path(out_dir, "rls_stage0.rds"))
lpa <- dat %>% filter(lpa_sample)
stopifnot(nrow(lpa) == 31470)

indices <- c("religbehav", "religbelif", "spirbehav", "spirbelif")
sds     <- sapply(lpa[indices], sd)

# Main solution (Stage 4): profile probabilities and assigned profiles
post <- readRDS(file.path(out_dir, "stage4_posteriors.rds"))
stopifnot(identical(post$respid, lpa$respid))
main_pp    <- post$A5[, 1:5]
main_modal <- post$A5[, "modal"]

n_starts <- 5

# Fits a type A model from n_starts starting points and keeps the best one
fit_best <- function(data, k, seed) {
  fits <- lapply(seq_len(n_starts), function(i) {
    set.seed(seed + i)
    tryCatch({
      m <- mx_profiles(data = data, classes = k, variances = "equal",
                       covariances = "zero", run = FALSE)
      mxTryHard(m, extraTries = 5, silent = TRUE, verbose = FALSE,
                bestInitsOutput = FALSE)
    }, error = function(e) NULL)
  })
  fits <- fits[!sapply(fits, is.null)]
  if (length(fits) == 0) stop("no starting point produced a solution")
  m2ll <- sapply(fits, function(m) m$output$Minus2LogLikelihood)
  best <- fits[[which.min(m2ll)]]
  attr(best, "starts_best") <- sum(abs(m2ll - min(m2ll)) < 0.01)
  best
}

# Loads a saved model, or fits and saves it
fit_or_load <- function(name, data, k, seed = 20261023) {
  f <- file.path(rob_dir, paste0(name, ".rds"))
  if (file.exists(f)) {
    message(name, ": already saved, loaded")
    return(readRDS(f))
  }
  t0 <- Sys.time()
  m  <- fit_best(data, k, seed)
  saveRDS(m, f)
  message(name, ": done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
          " minutes (", attr(m, "starts_best"), " of ", n_starts,
          " starting points reached the best solution)")
  m
}

# Profile probabilities (respondents x profiles)
post_of <- function(m) {
  pp <- class_prob(m, type = "individual")$individual
  pp[, 1:(ncol(pp) - 1), drop = FALSE]
}

# Probability-weighted averages of the four original indices, per profile
w_means <- function(pp, data) {
  t(sapply(seq_len(ncol(pp)), function(j) sapply(indices, function(v) {
    ok <- !is.na(data[[v]])
    sum(pp[ok, j] * data[[v]][ok]) / sum(pp[ok, j])
  })))
}

main_means <- w_means(main_pp, lpa)

# All orderings of 1..n
perms <- function(n) {
  if (n == 1) return(matrix(1))
  p <- perms(n - 1)
  do.call(rbind, lapply(1:n, function(i) cbind(i, ifelse(p >= i, p + 1, p))))
}

# Pairs the profiles of a check with the main profiles (closest averages).
# Returns, for each main profile, the number of the matching check profile.
match_to_main <- function(cand_means) {
  kc <- nrow(cand_means); km <- nrow(main_means)
  z  <- function(x) sweep(x, 2, sds, "/")
  if (kc >= km) {
    best <- NULL; best_d <- Inf
    for (sel in combn(kc, km, simplify = FALSE)) {
      P <- perms(km)
      for (r in seq_len(nrow(P))) {
        idx <- sel[P[r, ]]
        d <- sum((z(cand_means[idx, , drop = FALSE]) - z(main_means))^2)
        if (d < best_d) { best_d <- d; best <- idx }
      }
    }
  } else {
    best <- rep(NA, km); best_d <- Inf
    P <- perms(km)
    for (r in seq_len(nrow(P))) {
      idx <- P[r, 1:kc]
      d <- sum((z(cand_means) - z(main_means[idx, , drop = FALSE]))^2)
      if (d < best_d) { best_d <- d; best <- rep(NA, km); best[idx] <- 1:kc }
    }
  }
  best
}

# Compares one check with the main solution
compare <- function(label, m, data, same_people = NULL) {
  pp  <- post_of(m)
  cm  <- w_means(pp, data)
  mt  <- match_to_main(cm)
  ft  <- as.data.frame(table_fit(m))
  ok  <- !is.na(mt)
  diff_sd <- abs(sweep(cm[mt[ok], , drop = FALSE] - main_means[ok, , drop = FALSE], 2, sds, "/"))
  agree <- NA
  if (!is.null(same_people)) {
    modal_c <- apply(pp[same_people, , drop = FALSE], 1, which.max)
    relab   <- match(modal_c, mt)          # check profile -> main profile number
    agree   <- round(mean(!is.na(relab) & relab == main_modal), 3)
  }
  summary_row <- data.frame(check = label, n = nrow(data), profiles = ncol(pp),
                            BIC = round(ft$BIC, 1), entropy = round(ft$Entropy, 3),
                            prob_min = round(ft$prob_min, 3), n_min = round(ft$n_min, 3),
                            starts_best = attr(m, "starts_best"),
                            max_diff_sd = round(max(diff_sd), 2), agreement = agree)
  means_long <- do.call(rbind, lapply(which(ok), function(j) {
    data.frame(check = label, main_profile = j,
               share = round(mean(apply(pp, 1, which.max) == mt[j]), 3),
               t(round(cm[mt[j], ], 2)), row.names = NULL)
  }))
  list(summary = summary_row, means = means_long)
}

results <- list()
results[["main"]] <- list(
  summary = data.frame(check = "Main solution", n = nrow(lpa), profiles = 5,
                       BIC = NA, entropy = NA, prob_min = NA, n_min = NA,
                       starts_best = NA, max_diff_sd = 0, agreement = 1),
  means = data.frame(check = "Main solution", main_profile = 1:5,
                     share = round(as.vector(table(factor(main_modal, levels = 1:5))) / nrow(lpa), 3),
                     round(main_means, 2), row.names = NULL))


# ---- 6.1 Split-half replication ---------------------------------------------

for (h in 1:2) {
  dh <- as.data.frame(lpa[lpa$half == h, indices])
  m  <- fit_or_load(paste0("half", h, "_k5"), dh, 5)
  results[[paste0("half", h)]] <- compare(paste0("Half ", h, " only"), m, dh)
}


# ---- 6.2 Version 2 indices (three indices) ----------------------------------
# Religious practice: 9 religious behaviour items + importance of religion +
#   seeing oneself as religious
# Belief: belief in God, afterlife, seeing oneself as spiritual, soul,
#   spirit world, and the two experience items
# Spiritual practice: the five spiritual behaviour items
# Raw sums, not rescaled. Each respondent in the analytic sample has all 23
# items, so the sample is the same 31,470 respondents.

v2 <- data.frame(
  relprac  = rowSums(lpa[, c("attndperrlsr", "attndonrlsr", "synmemr", "prayfreqr", "gracefreqr",
                             "prac_ar", "prac_br", "prac_cr", "prac_dr", "relimpr", "relperr")]),
  belief   = rowSums(lpa[, c("gdcertainr", "afterlife", "spirperr", "soulr", "spirworld2",
                             "exp_ar", "exp_dr")]),
  spirprac = rowSums(lpa[, c("spiract_ar", "spiract_br", "spiract_cr", "spiract_dr", "spiract_fr")])
)
stopifnot(!anyNA(v2))

v2_fit <- list()
for (k in 4:6) {
  m <- fit_or_load(paste0("v2_k", k), v2, k)
  v2_fit[[paste0("k", k)]] <- m
  results[[paste0("v2_k", k)]] <- compare(paste0("Version 2 indices, ", k, " profiles"), m,
                                          as.data.frame(lpa[, indices]),
                                          same_people = seq_len(nrow(lpa)))
}


# ---- 6.3 Missing data: all respondents with at least one index --------------

has_any <- rowSums(!is.na(dat[, indices])) > 0
df_fiml <- as.data.frame(dat[has_any, indices])
in_lpa  <- which(dat$lpa_sample[has_any])        # rows of df_fiml in the analytic sample
stopifnot(all(as.numeric(dat$respid[has_any][in_lpa]) == as.numeric(lpa$respid)))

# Starting values: the main 5-profile solution (Stage 3). The usual automatic
# starting values take many hours when some scores are missing; starting from
# the main solution avoids that. mxTryHard then also tries perturbed starts.
f_fiml <- file.path(rob_dir, "fiml_k5.rds")
if (file.exists(f_fiml)) {
  m <- readRDS(f_fiml)
  message("fiml_k5: already saved, loaded")
} else {
  t0   <- Sys.time()
  main <- readRDS(file.path(out_dir, "stage3_models", "typeA_k5.rds"))
  m    <- mx_profiles(data = df_fiml, classes = 5, variances = "equal",
                      covariances = "zero", run = FALSE, run_mixture_starts = FALSE)
  p    <- omxGetParameters(main)
  m    <- omxSetParameters(m, labels = names(p), values = p)
  set.seed(20261023)
  m    <- mxTryHard(m, extraTries = 5, silent = TRUE, verbose = FALSE,
                    bestInitsOutput = FALSE)
  attr(m, "starts_best") <- NA
  saveRDS(m, f_fiml)
  message("fiml_k5: done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
}
results[["fiml"]] <- compare("All available data (FIML), 5 profiles", m, df_fiml,
                             same_people = in_lpa)



# ---- 6.4 Adjacent solutions: 4 and 6 profiles -------------------------------
# The 4- and 6-profile models from Stage 3 (type A), compared with the main
# 5-profile solution. Shows which main profiles merge (4 profiles) and which
# profile splits or is added (6 profiles). The full table lists every profile
# of each adjacent model with the main profile it matches ("none" = no match).

adj <- list(k4 = readRDS(file.path(out_dir, "stage3_models", "typeA_k4.rds")),
            k6 = readRDS(file.path(out_dir, "stage3_models", "typeA_k6.rds")))

adj_full <- list()
for (nm in names(adj)) {
  if (is.null(attr(adj[[nm]], "starts_best"))) attr(adj[[nm]], "starts_best") <- NA
  kk <- sub("k", "", nm)
  results[[paste0("adjacent_", nm)]] <- compare(paste0("Adjacent solution, ", kk, " profiles"),
                                                adj[[nm]], as.data.frame(lpa[, indices]),
                                                same_people = seq_len(nrow(lpa)))
  pp_a <- post_of(adj[[nm]])
  cm_a <- w_means(pp_a, lpa)
  mt_a <- match_to_main(cm_a)
  adj_full[[nm]] <- data.frame(model = paste(kk, "profiles"),
                               profile = seq_len(nrow(cm_a)),
                               matches_main_profile = sapply(seq_len(nrow(cm_a)), function(j) {
                                 hit <- which(mt_a == j); if (length(hit)) hit else "none" }),
                               share = round(colMeans(pp_a), 3),
                               round(cm_a, 2), row.names = NULL)
}
adj_full <- do.call(rbind, adj_full); rownames(adj_full) <- NULL
write.csv(adj_full, file.path(out_dir, "stage6_adjacent_profiles.csv"), row.names = FALSE)
cat("\n==== Adjacent solutions: every profile, with the main profile it matches ====\n")
print(adj_full, row.names = FALSE)


# ---- 6.5 Weighted estimation ------------------------------------------------
# The main model was estimated without survey weights. Here the same model
# (5 profiles, equal spreads, indices unrelated within profiles) is
# re-estimated with the survey weight (rescaled to average 1), by weighted
# maximum likelihood (pseudo-likelihood), using the EM algorithm, starting
# from the main solution. (OpenMx does not apply case weights to mixture
# models, so the estimation is written out here; it is the standard EM
# algorithm for this model, with each respondent counted in proportion to his
# or her weight.)

em_weighted <- function(X, w, pp, maxit = 5000, tol = 1e-8) {
  X <- as.matrix(X); w <- w / mean(w); n <- nrow(X); kk <- ncol(pp)
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    # M-step: shares, profile means, common variances
    wp   <- pp * w
    pis  <- colSums(wp) / sum(w)
    mu   <- t(wp) %*% X / colSums(wp)                       # kk x indicators
    s2   <- colSums(Reduce(`+`, lapply(1:kk, function(j)
      wp[, j] * sweep(X, 2, mu[j, ])^2))) / sum(w)
    # E-step: profile probabilities
    logd <- sapply(1:kk, function(j) log(pis[j]) +
                     rowSums(dnorm(X, matrix(mu[j, ], n, ncol(X), byrow = TRUE),
                                   matrix(sqrt(s2), n, ncol(X), byrow = TRUE), log = TRUE)))
    mx   <- apply(logd, 1, max)
    lse  <- mx + log(rowSums(exp(logd - mx)))
    pp   <- exp(logd - lse)
    ll   <- sum(w * lse)
    if (abs(ll - ll_old) < tol * abs(ll)) break
    ll_old <- ll
  }
  list(pp = pp, mu = mu, pis = pis, s2 = s2, ll = ll, iterations = it)
}

ew <- em_weighted(lpa[, indices], lpa$weight, main_pp)
cat("\nWeighted EM converged in", ew$iterations, "iterations\n")

cm_w <- w_means(ew$pp, lpa)
mt_w <- match_to_main(cm_w)
modal_w <- apply(ew$pp, 1, which.max)
avepp_w <- sapply(1:5, function(j) mean(ew$pp[modal_w == j, j]))
ent_w   <- 1 + sum(ew$pp * log(ew$pp + 1e-12)) / (nrow(ew$pp) * log(5))
diff_w  <- abs(sweep(cm_w[mt_w, , drop = FALSE] - main_means, 2, sds, "/"))

results[["weighted"]] <- list(
  summary = data.frame(check = "Survey-weighted estimation, 5 profiles", n = nrow(lpa),
                       profiles = 5, BIC = NA, entropy = round(ent_w, 3),
                       prob_min = round(min(avepp_w), 3),
                       n_min = round(min(table(factor(modal_w, levels = 1:5))) / nrow(lpa), 3),
                       starts_best = NA, max_diff_sd = round(max(diff_w), 2),
                       agreement = round(mean(match(modal_w, mt_w) == main_modal), 3)),
  means = do.call(rbind, lapply(1:5, function(j) data.frame(
    check = "Survey-weighted estimation, 5 profiles", main_profile = j,
    share = round(ew$pis[mt_w[j]], 3),          # weighted population share
    t(round(cm_w[mt_w[j], ], 2)), row.names = NULL))))



# ---- 6.6 Summary ------------------------------------------------------------

summary_tab <- do.call(rbind, lapply(results, `[[`, "summary"))
means_tab   <- do.call(rbind, lapply(results, `[[`, "means"))
rownames(summary_tab) <- NULL; rownames(means_tab) <- NULL

write.csv(summary_tab, file.path(out_dir, "stage6_summary.csv"), row.names = FALSE)
write.csv(means_tab,   file.path(out_dir, "stage6_profile_means.csv"), row.names = FALSE)

print(summary_tab, row.names = FALSE)
print(means_tab, row.names = FALSE)




# ---- 6.7 Session information ------------------------------------------------

sessionInfo()
