# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 4.0 Packages, paths, data ----------------------------------------------

library(OpenMx)
library(tidySEM)
library(dplyr)

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")
fit_dir  <- file.path(out_dir, "stage3_models")

dat <- readRDS(file.path(out_dir, "rls_stage0.rds"))
lpa <- dat %>% filter(lpa_sample)
stopifnot(nrow(lpa) == 31470)

indices <- c("religbehav", "religbelif", "spirbehav", "spirbelif")
labels  <- c(religbehav = "Religious behaviour (9-46)",
             religbelif = "Religious belief (4-16)",
             spirbehav  = "Spiritual behaviour (5-25)",
             spirbelif  = "Spiritual belief (5-18)")
ranges  <- list(religbehav = c(9, 46), religbelif = c(4, 16),
                spirbehav  = c(5, 25), spirbelif  = c(5, 18))
df <- as.data.frame(lpa[, indices])

models <- list(A5 = readRDS(file.path(fit_dir, "typeA_k5.rds")),
               C3 = readRDS(file.path(fit_dir, "typeC_k3.rds")))

# Profile probabilities for each respondent, with profiles renumbered by
# average religious behaviour (highest = profile 1)
get_post <- function(m) {
  pp <- class_prob(m, type = "individual")$individual
  k  <- ncol(pp) - 1
  pp <- pp[, 1:k, drop = FALSE]
  wmean <- colSums(pp * df$religbehav) / colSums(pp)
  ord   <- order(wmean, decreasing = TRUE)
  pp    <- pp[, ord, drop = FALSE]
  colnames(pp) <- paste0("p", 1:k)
  pp
}
posts <- lapply(models, get_post)
modal <- lapply(posts, function(pp) apply(pp, 1, which.max))


# ---- 4.1 Profile sizes ------------------------------------------------------
# share_post: estimated share of the sample in each profile (sum of
# probabilities); n_modal and share_modal: respondents assigned to the profile
# they most probably belong to. Unweighted; weighted shares follow in Stage 7.

sizes <- do.call(rbind, lapply(names(models), function(nm) {
  pp <- posts[[nm]]; k <- ncol(pp)
  data.frame(model       = nm,
             profile     = 1:k,
             share_post  = round(colSums(pp) / nrow(pp), 3),
             n_modal     = as.vector(table(factor(modal[[nm]], levels = 1:k))),
             share_modal = round(as.vector(table(factor(modal[[nm]], levels = 1:k))) / nrow(pp), 3))
}))
write.csv(sizes, file.path(out_dir, "stage4_sizes.csv"), row.names = FALSE)
print(sizes, row.names = FALSE)


# ---- 4.2 Classification quality ---------------------------------------------
# Rows: the profile people were assigned to. Columns: their average
# probability of belonging to each profile. The diagonal shows how clearly
# people in each profile are assigned (.80 or higher is the Stage 3 rule).

for (nm in names(models)) {
  pp <- posts[[nm]]; k <- ncol(pp)
  avg <- t(sapply(1:k, function(j) colMeans(pp[modal[[nm]] == j, , drop = FALSE])))
  dimnames(avg) <- list(paste0("assigned_", 1:k), paste0("prob_", 1:k))
  avg <- round(avg, 3)
  write.csv(avg, file.path(out_dir, paste0("stage4_classification_", nm, ".csv")))
  cat("\n==== Model", nm, ": average assignment probabilities ====\n")
  print(avg)
}


# ---- 4.3 Profile averages and spreads ---------------------------------------
# Each respondent counts towards each profile in proportion to his or her
# probability of belonging to it.

desc <- do.call(rbind, lapply(names(models), function(nm) {
  pp <- posts[[nm]]; k <- ncol(pp)
  do.call(rbind, lapply(1:k, function(j) {
    w <- pp[, j]
    data.frame(model   = nm,
               profile = j,
               index   = indices,
               mean    = round(sapply(indices, function(v) sum(w * df[[v]]) / sum(w)), 2),
               sd      = round(sapply(indices, function(v) {
                 m <- sum(w * df[[v]]) / sum(w)
                 sqrt(sum(w * (df[[v]] - m)^2) / sum(w)) }), 2),
               row.names = NULL)
  }))
}))
write.csv(desc, file.path(out_dir, "stage4_profile_means.csv"), row.names = FALSE)

# Wide layout for reading: one row per profile, one column per index mean
means_wide <- reshape(desc[, c("model", "profile", "index", "mean")],
                      idvar = c("model", "profile"), timevar = "index",
                      direction = "wide")
names(means_wide) <- sub("^mean\\.", "", names(means_wide))
means_wide <- merge(means_wide, sizes[, c("model", "profile", "share_post")],
                    by = c("model", "profile"))
print(means_wide, row.names = FALSE)


# ---- 4.4 Share at the lowest and highest possible score ---------------------
# Among respondents assigned to each profile. Shows whether a profile is made
# up mainly of people at the floor or ceiling of an index.

ends <- do.call(rbind, lapply(names(models), function(nm) {
  k <- ncol(posts[[nm]])
  do.call(rbind, lapply(1:k, function(j) {
    sub <- df[modal[[nm]] == j, , drop = FALSE]
    data.frame(model   = nm,
               profile = j,
               index   = indices,
               pct_at_min = round(100 * sapply(indices, function(v) mean(sub[[v]] == ranges[[v]][1])), 1),
               pct_at_max = round(100 * sapply(indices, function(v) mean(sub[[v]] == ranges[[v]][2])), 1),
               row.names = NULL)
  }))
}))
write.csv(ends, file.path(out_dir, "stage4_floor_ceiling.csv"), row.names = FALSE)
print(ends, row.names = FALSE)


# ---- 4.5 Profile plot -------------------------------------------------------
# One panel per index, on its own original scale (no rescaling).
# Bars show profile averages; the vertical axis covers the full possible range.

for (nm in names(models)) {
  k <- ncol(posts[[nm]])
  png(file.path(out_dir, paste0("stage4_profiles_", nm, ".png")),
      width = 2200, height = 1600, res = 200)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (v in indices) {
    d <- desc[desc$model == nm & desc$index == v, ]
    bp <- barplot(d$mean, names.arg = paste("P", d$profile), ylim = ranges[[v]],
                  xpd = FALSE, main = labels[[v]], ylab = "Average score",
                  col = grey(seq(0.3, 0.85, length.out = k)), border = NA)
    abline(h = ranges[[v]][1])
  }
  dev.off()
}


# ---- 4.6 Save profile probabilities for later stages ------------------------

post_out <- list(respid = lpa$respid,
                 A5 = cbind(posts$A5, modal = modal$A5),
                 C3 = cbind(posts$C3, modal = modal$C3))
saveRDS(post_out, file.path(out_dir, "stage4_posteriors.rds"))


# ---- 4.7 Additional classification and description checks (model A5) --------
# (a) Odds of correct classification (Nagin, 2005): for each profile, the odds
#     of correct assignment divided by the odds expected from the profile's
#     size alone. Values above 5 indicate good classification.
# (b) Local independence: correlations among the four indices within each
#     profile, with each respondent weighted by his or her probability of
#     belonging to the profile. The selected model assumes these are zero;
#     values below .30 in absolute size are usually treated as small.
# (c) Profile means with 95% confidence intervals: weighted by the survey
#     weight, with standard errors from the 500 replicate weights (profile
#     probabilities held fixed). Also shown on a standardised metric: the
#     difference from the weighted mean of the analytic sample, in weighted
#     standard deviations.

pp5 <- posts$A5
k5  <- ncol(pp5)

# (a) Odds of correct classification
avepp <- sapply(1:k5, function(j) mean(pp5[modal$A5 == j, j]))
pi_j  <- colSums(pp5) / nrow(pp5)
occ   <- data.frame(profile = 1:k5,
                    avg_prob = round(avepp, 3),
                    share = round(pi_j, 3),
                    occ = round((avepp / (1 - avepp)) / (pi_j / (1 - pi_j)), 1))
write.csv(occ, file.path(out_dir, "stage4_occ_A5.csv"), row.names = FALSE)
cat("\n==== Odds of correct classification (model A5) ====\n")
print(occ, row.names = FALSE)

# (b) Within-profile correlations
li <- do.call(rbind, lapply(1:k5, function(j) {
  r  <- cov.wt(df, wt = pp5[, j], cor = TRUE)$cor
  ut <- which(upper.tri(r), arr.ind = TRUE)
  data.frame(profile = j,
             pair = paste(indices[ut[, 1]], indices[ut[, 2]], sep = " - "),
             r = round(r[ut], 3))
}))
write.csv(li, file.path(out_dir, "stage4_within_profile_cor_A5.csv"), row.names = FALSE)
li_wide <- reshape(li, idvar = "pair", timevar = "profile", direction = "wide")
names(li_wide) <- sub("^r\\.", "profile_", names(li_wide))
cat("\n==== Correlations among the indices within each profile (model A5) ====\n")
print(li_wide, row.names = FALSE)

# (c) Weighted profile means with confidence intervals
reps <- grep("^repwt_[0-9]+$", names(lpa), value = TRUE)
wmeans <- function(wt) {
  sapply(indices, function(v) sapply(1:k5, function(j)
    sum(wt * pp5[, j] * df[[v]]) / sum(wt * pp5[, j])))
}
m_full <- wmeans(lpa$weight)                      # profiles x indices
m_reps <- lapply(reps, function(r) wmeans(lpa[[r]]))
m_se   <- sqrt(Reduce(`+`, lapply(m_reps, function(m) (m - m_full)^2)) / (length(reps) - 1))

grand_m  <- sapply(indices, function(v) weighted.mean(df[[v]], lpa$weight))
grand_sd <- sapply(indices, function(v) sqrt(sum(lpa$weight * (df[[v]] - grand_m[v])^2) / sum(lpa$weight)))

ci_tab <- do.call(rbind, lapply(1:k5, function(j) data.frame(
  profile = j, index = indices,
  mean    = round(m_full[j, ], 2),
  ci_low  = round(m_full[j, ] - 1.96 * m_se[j, ], 2),
  ci_high = round(m_full[j, ] + 1.96 * m_se[j, ], 2),
  standardised = round((m_full[j, ] - grand_m) / grand_sd, 2),
  row.names = NULL)))
write.csv(ci_tab, file.path(out_dir, "stage4_profile_means_ci_A5.csv"), row.names = FALSE)
cat("\n==== Weighted profile means with 95% CI, and standardised means (model A5) ====\n")
print(ci_tab, row.names = FALSE)



# ---- 4.8 Session information -------------------------------------------------

sessionInfo()