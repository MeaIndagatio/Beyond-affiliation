# ==============================================================================
# Latent profile analysis of religious and spiritual belief and practice,
# 2023-24 Religious Landscape Study
# ==============================================================================


# ---- 3.0 Packages, paths, data ----------------------------------------------

pkgs <- c("tidySEM", "OpenMx", "dplyr")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(OpenMx)
library(tidySEM)
library(dplyr)

# Use all but one processor core
mxOption(NULL, "Number of Threads", max(1, parallel::detectCores() - 1))

mxOption(NULL, "Calculate Hessian", "No")
mxOption(NULL, "Standard Errors", "No")


data_dir <- "path/to/data"
out_dir  <- file.path(data_dir, "R output")
fit_dir  <- file.path(out_dir, "stage3_models")

dat <- readRDS(file.path(out_dir, "rls_stage0.rds"))
lpa <- dat %>% filter(lpa_sample)
stopifnot(nrow(lpa) == 31470)

indices <- c("religbehav", "religbelif", "spirbehav", "spirbelif")
df <- as.data.frame(lpa[, indices])

max_k <- 8

specs <- data.frame(
  type        = c("A", "B", "C", "D"),
  variances   = c("equal", "varying", "equal", "varying"),
  covariances = c("zero", "zero", "equal", "varying"),
  stringsAsFactors = FALSE
)

# Number of different starting points per model. A model can stop at a
# solution that is not the best one; fitting it from several starting points
# and keeping the best guards against this. The number of starting points that
# reach the best solution is recorded ("starts_best"); if only one start
# reaches it, the solution is less trustworthy.
n_starts <- 10   # first run: all models, saved in stage3_models

# Fits one model from n_starts starting points and keeps the best solution.
fit_best <- function(k, variances, covariances, seed) {
  fits <- lapply(seq_len(n_starts), function(i) {
    set.seed(seed + i)
    tryCatch({
      m <- mx_profiles(data = df, classes = k, variances = variances,
                       covariances = covariances, run = FALSE)
      mxTryHard(m, extraTries = 5, silent = TRUE, verbose = FALSE,
                bestInitsOutput = FALSE)
    }, error = function(e) NULL)
  })
  fits <- fits[!sapply(fits, is.null)]
  if (length(fits) == 0) stop("no starting point produced a solution")
  m2ll <- sapply(fits, function(m) m$output$Minus2LogLikelihood)
  best <- fits[[which.min(m2ll)]]
  attr(best, "starts_best") <- sum(abs(m2ll - min(m2ll)) < 0.01)
  attr(best, "starts_ok")   <- length(fits)
  best
}

# Fits one model type with 1 to max_k profiles; each model is saved to its own
# file. Models that already have a file are skipped.
run_type <- function(type) {
  s <- specs[specs$type == type, ]
  for (k in 1:max_k) {
    f <- file.path(fit_dir, paste0("type", type, "_k", k, ".rds"))
    if (file.exists(f)) {
      message("Type ", type, ", ", k, " profiles: already saved, skipped")
      next
    }
    t0 <- Sys.time()
    m <- tryCatch(fit_best(k, s$variances, s$covariances, seed = 20261016 + 100 * k),
                  error = function(e) e)
    mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    if (inherits(m, "error")) {
      message("Type ", type, ", ", k, " profiles: FAILED after ", mins,
              " minutes: ", conditionMessage(m))
      next
    }
    saveRDS(m, f)
    message("Type ", type, ", ", k, " profiles: done in ", mins, " minutes (",
            attr(m, "starts_best"), " of ", attr(m, "starts_ok"),
            " starting points reached the best solution)")
  }
}


# ---- 3.1 Model type A -------------------------------------------------------

# run_type("A")


# ---- 3.2 Model type B -------------------------------------------------------

# run_type("B")


# ---- 3.3 Model type C -------------------------------------------------------

# run_type("C")


# ---- 3.4 Model type D -------------------------------------------------------

# run_type("D")   # not estimated: see Methods (same problem as type B)


fit_dir <- file.path(out_dir, "stage3_models_100starts")
dir.create(fit_dir, showWarnings = FALSE)

# Check: the four models that determined the selection, re-estimated
# from 100 starting values
n_starts <- 100
to_run <- data.frame(type = c("A", "A", "A", "C"), k = c(4, 5, 6, 3))

for (i in seq_len(nrow(to_run))) {
  type <- to_run$type[i]; k <- to_run$k[i]
  s <- specs[specs$type == type, ]
  f <- file.path(fit_dir, paste0("type", type, "_k", k, ".rds"))
  if (file.exists(f)) next
  m <- fit_best(k, s$variances, s$covariances, seed = 20261016 + 100 * k)
  saveRDS(m, f)
  message(type, k, ": ", attr(m, "starts_best"), " of ",
          attr(m, "starts_ok"), " starts reached the best solution")
}







# ---- 3.5 Fit table for all models -------------------------------------------
# Can be run at any time; it uses whichever models have been saved so far.
# status: 0 or 1 = the estimation finished normally; 5 or 6 = a solution was
# found but the program could not confirm it is the best possible; other
# values = a problem. starts_best = number of the 10 starting points that
# reached the best solution.

load_type <- function(type) {
  files <- file.path(fit_dir, paste0("type", type, "_k", 1:max_k, ".rds"))
  files <- files[file.exists(files)]
  mods  <- lapply(files, readRDS)
  names(mods) <- sub("\\.rds$", "", basename(files))
  mods
}

all_mods <- lapply(specs$type, load_type)
names(all_mods) <- specs$type
message("All saved models loaded")

fit_rows <- list()
for (type in specs$type) {
  for (nm in names(all_mods[[type]])) {
    m  <- all_mods[[type]][[nm]]
    t0 <- Sys.time()
    ft <- as.data.frame(table_fit(m))
    fit_rows[[nm]] <- data.frame(type        = type,
                                 profiles    = ft$Classes,
                                 LL          = round(ft$LL, 1),
                                 params      = ft$Parameters,
                                 AIC         = round(ft$AIC, 1),
                                 BIC         = round(ft$BIC, 1),
                                 saBIC       = round(ft$saBIC, 1),
                                 entropy     = round(ft$Entropy, 3),
                                 prob_min    = round(ft$prob_min, 3),
                                 n_min       = round(ft$n_min, 3),
                                 status      = m$output$status$code,
                                 starts_best = attr(m, "starts_best"))
    message(nm, ": ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " seconds")
  }
}
fit_table <- do.call(rbind, fit_rows)
rownames(fit_table) <- NULL



write.csv(fit_table, file.path(out_dir, "stage3_fit.csv"), row.names = FALSE)
print(fit_table, row.names = FALSE)


# ---- 3.6 Plot of BIC by number of profiles ----------------------------------

png(file.path(out_dir, "stage3_bic.png"), width = 1800, height = 1400, res = 200)
plot(NULL, xlim = c(1, max_k), ylim = range(fit_table$BIC, na.rm = TRUE),
     xlab = "Number of profiles", ylab = "BIC (lower is better)",
     main = "BIC by number of profiles and model type")
for (i in seq_along(specs$type)) {
  d <- fit_table[fit_table$type == specs$type[i], ]
  lines(d$profiles, d$BIC, type = "b", pch = i, lty = i)
}
legend("topright", legend = paste("Type", specs$type), pch = 1:4, lty = 1:4, bty = "n")
dev.off()


# ---- 3.7 Session information ------------------------------------------------

sessionInfo()