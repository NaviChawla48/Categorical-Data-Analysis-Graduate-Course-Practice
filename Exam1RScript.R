# Math 536 Exam 1 - DARPA Drone Distance Error Analysis
library(ggplot2)
library(boot)
library(lmtest)
library(car)

# ── 1. Load Data ─────────────────────────────────────────────
dat <- read.csv("C:/Users/navic/Downloads/Categorical Data Analysis/exam1.csv", header = TRUE)
colnames(dat) <- c("distance", "error")

# ── 2. Exploratory Plot ──────────────────────────────────────
p_raw <- ggplot(dat, aes(x = distance, y = error)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  labs(
    title    = "Drone Position Error vs. True Distance",
    subtitle = "Raw scale — oscillating behavior visible",
    x        = "True Distance (km)",
    y        = "Position Error (m)"
  ) +
  theme_bw(base_size = 12)

p_log <- ggplot(dat, aes(x = log(distance), y = error)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 1) +
  labs(
    title    = "Drone Position Error vs. log(Distance)",
    subtitle = "Log-transformed predictor linearizes the mean trend",
    x        = "log(True Distance) (log km)",
    y        = "Position Error (m)"
  ) +
  theme_bw(base_size = 12)

# ── 3. Fit Model: error ~ log(distance) ─────────────────────
# log(x) as predictor — cannot log(y) since errors are negative
dat$log_dist <- log(dat$distance)
fit <- lm(error ~ log_dist, data = dat)
summary(fit)

# ── 4. Assumption Checks ─────────────────────────────────────

# 4a. Residuals vs Fitted
p_rvf <- ggplot(data.frame(fitted = fitted(fit), resid = resid(fit)),
                aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick") +
  geom_smooth(method = "loess", se = FALSE, color = "darkorange", linewidth = 0.8) +
  labs(title = "Residuals vs Fitted",
       x = "Fitted Values", y = "Residuals") +
  theme_bw(base_size = 12)

# 4b. Normal Q-Q
p_qq <- ggplot(data.frame(resid = resid(fit)), aes(sample = resid)) +
  stat_qq(alpha = 0.4, size = 0.8, color = "steelblue") +
  stat_qq_line(color = "firebrick", linewidth = 1) +
  labs(title = "Normal Q-Q Plot of Residuals",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_bw(base_size = 12)

# 4c. Scale-Location (sqrt |standardized residuals| vs fitted)
p_sl <- ggplot(data.frame(fitted = fitted(fit),
                           sqrt_abs_resid = sqrt(abs(rstandard(fit)))),
               aes(x = fitted, y = sqrt_abs_resid)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "darkorange", linewidth = 0.8) +
  labs(title = "Scale-Location",
       x = "Fitted Values", y = expression(sqrt("|Standardized Residuals|"))) +
  theme_bw(base_size = 12)

# 4d. Residuals vs Leverage
p_lev <- ggplot(data.frame(leverage  = hatvalues(fit),
                             std_resid = rstandard(fit)),
                aes(x = leverage, y = std_resid)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick") +
  geom_vline(xintercept = 2 * 2 / nrow(dat),   # 2p/n threshold
             linetype = "dotted", color = "gray50") +
  labs(title = "Residuals vs Leverage",
       x = "Leverage", y = "Standardized Residuals") +
  theme_bw(base_size = 12)

# 4e. Formal tests
bp_test   <- bptest(fit)          # Breusch-Pagan: homoscedasticity
sw_test <- shapiro.test(resid(fit))
# Shapiro-Wilk max n=5000; subsample if needed
dw_test   <- dwtest(fit)          # Durbin-Watson: independence

cat("\n── Breusch-Pagan Test (Homoscedasticity) ──\n"); print(bp_test)
cat("\n── Shapiro-Wilk Test (Normality) ──\n");          print(sw_test)
cat("\n── Durbin-Watson Test (Independence) ──\n");      print(dw_test)

# ── 5. Nonparametric Bootstrap CIs ──────────────────────────
# Strategy: resample (x, y) pairs (case resampling bootstrap)
# This is fully nonparametric — no distributional assumptions needed
# Handles heteroscedasticity and non-normality automatically

set.seed(536)
B <- 10000  # number of bootstrap replicates

# New points to predict at
new_x <- data.frame(log_dist = log(c(0.1, 1, 10)))
x_vals <- c(0.1, 1, 10)

# Bootstrap function: resample rows, refit, predict
boot_fn <- function(data, indices) {
  d   <- data[indices, ]
  fit_b <- lm(error ~ log_dist, data = d)
  predict(fit_b, newdata = new_x)
}

boot_out <- boot(data = dat, statistic = boot_fn, R = B)

# Extract bootstrap CIs for each prediction point (BCa method — best)
ci_01  <- boot.ci(boot_out, type = "bca", index = 1)  # 0.1 km
ci_1   <- boot.ci(boot_out, type = "bca", index = 2)  # 1 km
ci_10  <- boot.ci(boot_out, type = "bca", index = 3)  # 10 km

# Point estimates from original model
pred_pts <- predict(fit, newdata = new_x)

# ── 6. Summary Table ─────────────────────────────────────────
results <- data.frame(
  Distance_km        = x_vals,
  Predicted_Error_m  = round(pred_pts, 4),
  CI_Lower_95        = round(c(ci_01$bca[4], ci_1$bca[4], ci_10$bca[4]), 4),
  CI_Upper_95        = round(c(ci_01$bca[5], ci_1$bca[5], ci_10$bca[5]), 4)
)
print(results)

# ── 7. Prediction Plot with Bootstrap CI band ────────────────
# Build a smooth CI band across full distance range
dist_seq    <- exp(seq(log(min(dat$distance)), log(max(dat$distance)), length.out = 300))
newdat_seq  <- data.frame(log_dist = log(dist_seq))

# Bootstrap CI band (slower — use 1000 reps for the band)
set.seed(536)
B_band <- 2000
boot_band <- matrix(NA, nrow = B_band, ncol = 300)
for (b in seq_len(B_band)) {
  idx      <- sample(nrow(dat), replace = TRUE)
  fit_b    <- lm(error ~ log_dist, data = dat[idx, ])
  boot_band[b, ] <- predict(fit_b, newdata = newdat_seq)
}
ci_lo <- apply(boot_band, 2, quantile, 0.025)
ci_hi <- apply(boot_band, 2, quantile, 0.975)

band_df <- data.frame(distance = dist_seq,
                      fit      = predict(fit, newdata = newdat_seq),
                      lo       = ci_lo,
                      hi       = ci_hi)

pred_df <- data.frame(
  distance = x_vals,
  pred     = pred_pts,
  lo       = c(ci_01$bca[4], ci_1$bca[4], ci_10$bca[4]),
  hi       = c(ci_01$bca[5], ci_1$bca[5], ci_10$bca[5]),
  label    = c("0.1 km", "1 km", "10 km")
)

p_final <- ggplot() +
  geom_point(data = dat, aes(x = distance, y = error),
             alpha = 0.2, size = 0.6, color = "steelblue") +
  geom_ribbon(data = band_df, aes(x = distance, ymin = lo, ymax = hi),
              fill = "firebrick", alpha = 0.2) +
  geom_line(data = band_df, aes(x = distance, y = fit),
            color = "firebrick", linewidth = 1) +
  geom_point(data = pred_df, aes(x = distance, y = pred),
             color = "darkorange", size = 3, shape = 18) +
  geom_errorbar(data = pred_df, aes(x = distance, ymin = lo, ymax = hi),
                color = "darkorange", width = 0.2, linewidth = 0.8) +
  geom_text(data = pred_df, aes(x = distance, y = hi + 0.4, label = label),
            size = 3.5, color = "darkorange") +
  scale_x_log10() +
  labs(
    title    = "Fitted Model with 95% Bootstrap Confidence Band",
    subtitle = "OLS fit with log(distance) predictor; CIs via nonparametric bootstrap (BCa, B=10,000)",
    x        = "True Distance (km, log scale)",
    y        = "Position Error (m)"
  ) +
  theme_bw(base_size = 12)

# ── 8. Bootstrap distribution plots (transparency) ───────────
boot_df <- data.frame(
  value    = c(boot_out$t[,1], boot_out$t[,2], boot_out$t[,3]),
  Distance = rep(c("0.1 km", "1 km", "10 km"), each = B)
)

p_boot_dist <- ggplot(boot_df, aes(x = value, fill = Distance)) +
  geom_histogram(bins = 60, alpha = 0.7, color = "white") +
  facet_wrap(~Distance, scales = "free") +
  labs(
    title = "Bootstrap Distributions of Predicted Mean Error",
    x     = "Predicted Error (m)",
    y     = "Count"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# Print all plots
print(p_raw)
print(p_log)
print(p_rvf)
print(p_qq)
print(p_sl)
print(p_lev)
print(p_final)
print(p_boot_dist)
