# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 8.0 Packages, paths, data ----------------------------------------------

library(dplyr)

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")

dat  <- readRDS(file.path(out_dir, "rls_stage0.rds"))
post <- readRDS(file.path(out_dir, "stage4_posteriors.rds"))

pp <- as.data.frame(post$A5)
names(pp)[names(pp) == "modal"] <- "profile"
pp$respid  <- as.numeric(post$respid)
dat$respid <- as.numeric(dat$respid)
dat <- left_join(dat, pp, by = "respid")

lp   <- dat %>% filter(lpa_sample)
k    <- 5
reps <- grep("^repwt_[0-9]+$", names(lp), value = TRUE)
P    <- as.matrix(lp[, paste0("p", 1:k)])
s    <- lp$profile
w    <- lp$weight

# Classification error matrix: D[c, s] = weighted share of true profile c
# assigned to profile s; BCH weights: solve(D)[s_i, c]
D <- t(sapply(1:k, function(c) sapply(1:k, function(j) sum(w * P[, c] * (s == j)) / sum(w * P[, c]))))
H <- solve(D)[s, , drop = FALSE]

# Replicate standard error: squared deviations from the full-sample estimate
# divided by (R - 1), as in svrepdesign(type = "bootstrap", mse = TRUE)
rep_se <- function(full, rep_mat) sqrt(rowSums((rep_mat - full)^2) / (ncol(rep_mat) - 1))

# BCH-corrected mean of y in each profile, for weight vector wt
bch_means <- function(y, wt) {
  ok <- !is.na(y)
  sapply(1:k, function(c) sum(wt[ok] * H[ok, c] * y[ok]) / sum(wt[ok] * H[ok, c]))
}


# ---- 8.1 Lifetime change in religiousness and spirituality ------------------
# Percentage of each profile giving each answer (each profile's percentages
# add up to 100 for each question).

change_labels <- list(
  lifedir  = c("1" = "Become more religious", "2" = "Become less religious",
               "3" = "Some of both", "4" = "Stayed about the same"),
  lifespir = c("1" = "Become more spiritual", "2" = "Become less spiritual",
               "3" = "Some of both", "4" = "Stayed about the same")
)

change_tab <- do.call(rbind, lapply(names(change_labels), function(v) {
  do.call(rbind, lapply(names(change_labels[[v]]), function(code) {
    y   <- ifelse(is.na(lp[[v]]), NA, as.numeric(lp[[v]] == as.numeric(code)))
    est <- bch_means(y, w)
    rm  <- sapply(reps, function(r) bch_means(y, lp[[r]]))
    se  <- rep_se(est, rm)
    data.frame(question = v, answer = change_labels[[v]][[code]], profile = 1:k,
               pct = round(100 * est, 1), se = round(100 * se, 1),
               ci_low = round(100 * (est - 1.96 * se), 1),
               ci_high = round(100 * (est + 1.96 * se), 1))
  }))
}))
write.csv(change_tab, file.path(out_dir, "stage8_lifetime_change.csv"), row.names = FALSE)

change_wide <- reshape(change_tab[, c("question", "answer", "profile", "pct")],
                       idvar = c("question", "answer"), timevar = "profile",
                       direction = "wide")
names(change_wide) <- sub("^pct\\.", "profile_", names(change_wide))
cat("\n==== Lifetime change: percentage of each profile giving each answer ====\n")
print(change_wide, row.names = FALSE)


# ---- 8.2 Wellbeing by profile -----------------------------------------------

wb     <- lp$wellbeing
wb_est <- bch_means(wb, w)
wb_rep <- sapply(reps, function(r) bch_means(wb, lp[[r]]))
wb_se  <- rep_se(wb_est, wb_rep)

wb_tab <- data.frame(profile = 1:k, n_with_wellbeing = sapply(1:k, function(c) sum(s == c & !is.na(wb))),
                     mean = round(wb_est, 2), se = round(wb_se, 3),
                     ci_low = round(wb_est - 1.96 * wb_se, 2),
                     ci_high = round(wb_est + 1.96 * wb_se, 2))
write.csv(wb_tab, file.path(out_dir, "stage8_wellbeing_means.csv"), row.names = FALSE)
cat("\n==== Wellbeing (0-10) by profile ====\n")
print(wb_tab, row.names = FALSE)

# All pairwise differences, with Holm adjustment for 10 comparisons
pairs <- t(combn(k, 2))
pw <- data.frame(profile_a = pairs[, 1], profile_b = pairs[, 2])
pw$difference <- wb_est[pw$profile_a] - wb_est[pw$profile_b]
pw$se <- sapply(seq_len(nrow(pw)), function(i) {
  d_full <- pw$difference[i]
  d_rep  <- wb_rep[pw$profile_a[i], ] - wb_rep[pw$profile_b[i], ]
  sqrt(sum((d_rep - d_full)^2) / (length(reps) - 1))
})
pw$p      <- 2 * pnorm(-abs(pw$difference / pw$se))
pw$p_holm <- p.adjust(pw$p, method = "holm")
# Standardised difference: difference / weighted SD of wellbeing
wsd <- sqrt(sum(w * (wb - weighted.mean(wb, w, na.rm = TRUE))^2, na.rm = TRUE) /
              sum(w[!is.na(wb)]))
pw$d <- pw$difference / wsd
pw[, c("difference", "se", "d")] <- round(pw[, c("difference", "se", "d")], 3)
pw[, c("p", "p_holm")] <- signif(pw[, c("p", "p_holm")], 3)
write.csv(pw, file.path(out_dir, "stage8_wellbeing_pairwise.csv"), row.names = FALSE)
cat("\n==== Wellbeing: pairwise differences between profiles ====\n")
print(pw, row.names = FALSE)


# ---- 8.3 Wellbeing adjusted for religious tradition and demographics --------
# Weighted least-squares regression on data expanded to one row per
# respondent per profile, with weight = survey weight x BCH weight.
# Profile 2 is the reference. Covariates: current religious tradition, age
# group, gender, education, income. A joint test asks whether the profiles
# differ in wellbeing once tradition and demographics are held constant.

covs <- c("rel_curr9", "agegrp", "woman", "educrec", "income4")
cc   <- complete.cases(lp[, c("wellbeing", covs)])
ex   <- lp[rep(which(cc), each = k), c("wellbeing", covs)]
ex$prof <- factor(rep(1:k, times = sum(cc)), levels = c(2, 1, 3, 4, 5))
Hc   <- H[cc, , drop = FALSE]
hvec <- as.vector(t(Hc))                       # BCH weight for each expanded row
X    <- model.matrix(reformulate(c("prof", covs)), data = droplevels(ex))
y    <- ex$wellbeing

wls <- function(wt_resp) {
  wv <- rep(wt_resp[cc], each = k) * hvec
  solve(crossprod(X, X * wv), crossprod(X, wv * y))[, 1]
}
b_full <- wls(w)
b_rep  <- sapply(reps, function(r) wls(lp[[r]]))
V      <- (b_rep - b_full) %*% t(b_rep - b_full) / (length(reps) - 1)
se     <- sqrt(diag(V))

adj <- data.frame(term = names(b_full), estimate = round(b_full, 3), se = round(se, 3),
                  ci_low = round(b_full - 1.96 * se, 3), ci_high = round(b_full + 1.96 * se, 3),
                  p = signif(2 * pnorm(-abs(b_full / se)), 3), row.names = NULL)
write.csv(adj, file.path(out_dir, "stage8_wellbeing_adjusted.csv"), row.names = FALSE)
cat("\n==== Wellbeing adjusted for tradition and demographics (reference = profile 2) ====\n")
cat("Respondents:", sum(cc), "\n")
print(adj[grepl("^prof", adj$term), ], row.names = FALSE)

# Joint test that all four profile differences are zero
pt   <- grep("^prof", names(b_full))
wald <- as.numeric(t(b_full[pt]) %*% solve(V[pt, pt]) %*% b_full[pt])
cat("\nJoint test of profile differences: chi-square =", round(wald, 1),
    "on", length(pt), "df, p =", signif(pchisq(wald, length(pt), lower.tail = FALSE), 3), "\n")

# Same joint test for religious tradition, holding profile constant
tt    <- grep("^rel_curr9", names(b_full))
wald2 <- as.numeric(t(b_full[tt]) %*% solve(V[tt, tt]) %*% b_full[tt])
cat("Joint test of tradition differences: chi-square =", round(wald2, 1),
    "on", length(tt), "df, p =", signif(pchisq(wald2, length(tt), lower.tail = FALSE), 3), "\n")


# ---- 8.4 Session information ------------------------------------------------

sessionInfo()