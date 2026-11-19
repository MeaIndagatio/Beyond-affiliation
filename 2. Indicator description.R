# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 2.0 Packages, paths, data ----------------------------------------------

pkgs <- c("dplyr", "survey", "psych")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr)
library(survey)
library(psych)

data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")

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


# ---- 2.1 Unweighted descriptive statistics ----------------------------------
# Kurtosis is excess kurtosis (0 = normal distribution). Stata reports
# kurtosis without subtracting 3, so Stata value = this value + 3.

desc <- do.call(rbind, lapply(indices, function(v) {
  x <- lpa[[v]]
  data.frame(index        = v,
             n            = length(x),
             possible_min = ranges[[v]][1],
             possible_max = ranges[[v]][2],
             observed_min = min(x),
             observed_max = max(x),
             mean         = round(mean(x), 2),
             sd           = round(sd(x), 2),
             median       = median(x),
             skewness     = round(skew(x), 2),
             kurtosis     = round(kurtosi(x), 2),
             pct_at_min   = round(100 * mean(x == ranges[[v]][1]), 1),
             pct_at_max   = round(100 * mean(x == ranges[[v]][2]), 1))
}))
write.csv(desc, file.path(out_dir, "stage2_descriptives.csv"), row.names = FALSE)
print(desc, row.names = FALSE)


# ---- 2.2 Weighted means and standard deviations -----------------------------
# Same survey design as Stage 0 (full-sample weight plus 500 replicate
# weights), restricted to the analytic sample.

des <- svrepdesign(data             = dat,
                   weights          = ~weight,
                   repweights       = "repwt_[0-9]+",
                   type             = "bootstrap",
                   combined.weights = TRUE,
                   mse              = TRUE)
des_lpa <- subset(des, lpa_sample)

w_desc <- do.call(rbind, lapply(indices, function(v) {
  f  <- as.formula(paste0("~", v))
  m  <- svymean(f, des_lpa)
  vr <- svyvar(f, des_lpa)
  data.frame(index         = v,
             weighted_mean = round(as.numeric(coef(m)), 2),
             se_mean       = round(as.numeric(SE(m)), 3),
             weighted_sd   = round(sqrt(as.numeric(coef(vr))), 2))
}))
write.csv(w_desc, file.path(out_dir, "stage2_weighted_descriptives.csv"), row.names = FALSE)
print(w_desc, row.names = FALSE)


# ---- 2.3 Correlations among the four indices --------------------------------

cor_unw <- round(cor(lpa[indices]), 3)
cor_w   <- round(cov.wt(lpa[indices], wt = lpa$weight, cor = TRUE)$cor, 3)

write.csv(cor_unw, file.path(out_dir, "stage2_correlations_unweighted.csv"))
write.csv(cor_w,   file.path(out_dir, "stage2_correlations_weighted.csv"))
cat("\nUnweighted correlations\n"); print(cor_unw)
cat("\nWeighted correlations\n");   print(cor_w)


# ---- 2.4 Histograms ---------------------------------------------------------
# Saved to a file rather than drawn in the Plots pane, to avoid the
# "figure margins too large" error.

png(file.path(out_dir, "stage2_histograms.png"), width = 2000, height = 1600, res = 200)
par(mfrow = c(2, 2))
for (v in indices) {
  r <- ranges[[v]]
  hist(lpa[[v]], breaks = seq(r[1] - 0.5, r[2] + 0.5, 1),
       main = labels[[v]], xlab = "Score", ylab = "Respondents",
       col = "grey75", border = "white")
}
dev.off()


# ---- 2.5 Weighted means by religious tradition ------------------------------
# Respondents with missing tradition are excluded from this table only.

by_trad <- svyby(reformulate(indices), ~rel_curr9, des_lpa, svymean, na.rm = TRUE)
by_trad <- as.data.frame(by_trad)
n_trad  <- table(lpa$rel_curr9)
by_trad$n <- as.vector(n_trad[as.character(by_trad$rel_curr9)])
by_trad[, -1] <- round(by_trad[, -1], 2)
by_trad <- by_trad[, c(1, ncol(by_trad), 2:(ncol(by_trad) - 1))]
names(by_trad) <- c("rel_curr9", "n", indices, paste0("se_", indices))

write.csv(by_trad, file.path(out_dir, "stage2_means_by_tradition.csv"), row.names = FALSE)
print(by_trad, row.names = FALSE)


# ---- 2.6 Session information ------------------------------------------------

sessionInfo()