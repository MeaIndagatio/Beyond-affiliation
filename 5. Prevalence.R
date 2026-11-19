# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 5.0 Packages, paths, data ----------------------------------------------

library(dplyr)
library(survey)

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")

dat  <- readRDS(file.path(out_dir, "rls_stage0.rds"))
post <- readRDS(file.path(out_dir, "stage4_posteriors.rds"))

# Attach each respondent's profile probabilities (p1 to p5) and assigned
# profile; respondents outside the analytic sample get missing values.
pp <- as.data.frame(post$A5)
pp$respid <- post$respid
names(pp)[names(pp) == "modal"] <- "profile"
dat <- left_join(dat, pp, by = "respid")
stopifnot(sum(!is.na(dat$p1)) == 31470)

pcols <- paste0("p", 1:5)

des <- svrepdesign(data             = dat,
                   weights          = ~weight,
                   repweights       = "repwt_[0-9]+",
                   type             = "bootstrap",
                   combined.weights = TRUE,
                   mse              = TRUE)
des_lpa <- subset(des, lpa_sample)


# ---- 5.1 Weighted share of each profile -------------------------------------

sh <- svymean(reformulate(pcols), des_lpa)
ci <- confint(sh)
shares <- data.frame(profile   = 1:5,
                     pct       = round(100 * as.numeric(coef(sh)), 1),
                     ci_low    = round(100 * ci[, 1], 1),
                     ci_high   = round(100 * ci[, 2], 1),
                     pct_unweighted = round(100 * colMeans(dat[dat$lpa_sample, pcols]), 1),
                     row.names = NULL)
write.csv(shares, file.path(out_dir, "stage5_weighted_shares.csv"), row.names = FALSE)
print(shares, row.names = FALSE)


# ---- 5.2 Profile make-up of each religious tradition ------------------------
# Each row is a tradition; the five percentages add up to 100.

by_trad <- as.data.frame(svyby(reformulate(pcols), ~rel_curr9, des_lpa, svymean, na.rm = TRUE))
k   <- length(pcols)
trad_names <- as.character(by_trad[[1]])
est <- round(100 * as.matrix(by_trad[, 2:(k + 1)]), 1)
se  <- round(100 * as.matrix(by_trad[, (k + 2):(2 * k + 1)]), 1)
dimnames(est) <- list(trad_names, pcols)
dimnames(se)  <- list(trad_names, paste0("se_", pcols))
n_trad <- table(dat$rel_curr9[dat$lpa_sample])

within_trad <- data.frame(tradition = trad_names,
                          n = as.vector(n_trad[trad_names]),
                          est, se, row.names = NULL)
write.csv(within_trad, file.path(out_dir, "stage5_profiles_within_traditions.csv"), row.names = FALSE)
print(within_trad[, c("tradition", "n", pcols)], row.names = FALSE)


# ---- 5.3 Religious-tradition make-up of each profile ------------------------
# Each column is a profile; the nine percentages in a column add up to 100.
# For profile j, each respondent's weight (and each replicate weight) is
# multiplied by his or her probability of belonging to profile j.

lpa_dat <- dat[dat$lpa_sample & !is.na(dat$rel_curr9), ]
rep_cols <- grep("^repwt_[0-9]+$", names(lpa_dat), value = TRUE)

trad_in_profile <- lapply(1:k, function(j) {
  d <- lpa_dat
  d$w_j <- d$weight * d[[pcols[j]]]
  d[rep_cols] <- d[rep_cols] * d[[pcols[j]]]
  dj <- svrepdesign(data = d, weights = ~w_j, repweights = "repwt_[0-9]+",
                    type = "bootstrap", combined.weights = TRUE, mse = TRUE)
  m <- svymean(~rel_curr9, dj)
  data.frame(tradition = sub("^rel_curr9", "", names(coef(m))),
             profile   = j,
             pct       = round(100 * as.numeric(coef(m)), 1),
             se        = round(100 * as.numeric(SE(m)), 1),
             row.names = NULL)
})
trad_in_profile <- do.call(rbind, trad_in_profile)
write.csv(trad_in_profile, file.path(out_dir, "stage5_traditions_within_profiles.csv"), row.names = FALSE)

trad_wide <- reshape(trad_in_profile[, c("tradition", "profile", "pct")],
                     idvar = "tradition", timevar = "profile", direction = "wide")
names(trad_wide) <- sub("^pct\\.", "profile_", names(trad_wide))
print(trad_wide, row.names = FALSE)


# ---- 5.4 Profiles by religious tradition raised in --------------------------
# Same two tables as 5.2 and 5.3, for the tradition respondents were raised
# in (rel_raised9) instead of their current tradition.
# (a) Profile make-up of each raised-in tradition: each row adds up to 100.
# (b) Raised-in tradition make-up of each profile: each column adds up to 100.

# (a)
by_raised <- as.data.frame(svyby(reformulate(pcols), ~rel_raised9, des_lpa, svymean, na.rm = TRUE))
raised_names <- as.character(by_raised[[1]])
est_r <- round(100 * as.matrix(by_raised[, 2:(k + 1)]), 1)
se_r  <- round(100 * as.matrix(by_raised[, (k + 2):(2 * k + 1)]), 1)
dimnames(est_r) <- list(raised_names, pcols)
dimnames(se_r)  <- list(raised_names, paste0("se_", pcols))
n_raised <- table(dat$rel_raised9[dat$lpa_sample])

within_raised <- data.frame(raised_in = raised_names,
                            n = as.vector(n_raised[raised_names]),
                            est_r, se_r, row.names = NULL)
write.csv(within_raised, file.path(out_dir, "stage5_profiles_within_raised.csv"), row.names = FALSE)
cat("\n==== Profile make-up of each raised-in tradition (row %) ====\n")
print(within_raised[, c("raised_in", "n", pcols)], row.names = FALSE)

# (b)
lpa_dat_r <- dat[dat$lpa_sample & !is.na(dat$rel_raised9), ]

raised_in_profile <- do.call(rbind, lapply(1:k, function(j) {
  d <- lpa_dat_r
  d$w_j <- d$weight * d[[pcols[j]]]
  d[rep_cols] <- d[rep_cols] * d[[pcols[j]]]
  dj <- svrepdesign(data = d, weights = ~w_j, repweights = "repwt_[0-9]+",
                    type = "bootstrap", combined.weights = TRUE, mse = TRUE)
  m <- svymean(~rel_raised9, dj)
  data.frame(raised_in = sub("^rel_raised9", "", names(coef(m))),
             profile   = j,
             pct       = round(100 * as.numeric(coef(m)), 1),
             se        = round(100 * as.numeric(SE(m)), 1),
             row.names = NULL)
}))
write.csv(raised_in_profile, file.path(out_dir, "stage5_raised_within_profiles.csv"), row.names = FALSE)

raised_wide <- reshape(raised_in_profile[, c("raised_in", "profile", "pct")],
                       idvar = "raised_in", timevar = "profile", direction = "wide")
names(raised_wide) <- sub("^pct\\.", "profile_", names(raised_wide))
cat("\n==== Raised-in tradition make-up of each profile (column %) ====\n")
print(raised_wide, row.names = FALSE)



# ---- 5.5 Session information ------------------------------------------------

sessionInfo()