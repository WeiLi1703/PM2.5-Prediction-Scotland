# Project 66.2: Developing predictive models for PM2.5 concentrations

# 0. Packages, Conflicts & Reproducibility

## 0.1 Core Data Manipulation & Conflict Management
library(tidyverse)    # ggplot2, dplyr, tidyr, readr, purrr, etc.
library(conflicted)   # explicit conflict resolution
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::recode)

## 0.2 Date–Time and Spatial Utilities
library(lubridate)    # date/time handling
library(sf)           # simple features for spatial data

## 0.3 Visualization Extensions
library(ggthemes)     # additional themes for ggplot2
library(ggpubr)       # publication-ready ggplot utilities
library(GGally)       # ggpairs scatterplot matrix
library(corrplot)     # correlation heatmaps
library(scales)       # number/percent formatting & rescaling helpers

## 0.4 Modeling & Diagnostics
library(car)          # vif(), durbinWatsonTest()
library(lmtest)       # bptest()
library(performance)  # check_model(), model performance metrics
library(broom)        # tidy(), augment()
library(forecast)     # Acf() for autocorrelation diagnostics
library(sandwich)     # robust / cluster-robust covariance

## 0.5 Evaluation Metrics & Helpers
library(Metrics)      # rmse(), mae(), etc.

## 0.6 Regularized & Tree Models
library(glmnet)       # LASSO / elastic net
library(ranger)       # Random Forest (fast)

## 0.7 Misc
library(e1071)        # skewness()

## 0.8 Global random seed: ensure full script reproducibility
set.seed(2025)

# 1. Read Data & Initial Cleaning
pm25 <- read_csv("Data - PM25.csv") %>%
  mutate(
    Year          = as.integer(Year),
    MonthNum      = as.integer(Month),
    Month         = factor(Month, levels = 1:12, labels = month.abb),
    Date          = make_date(Year, MonthNum, 15),
    Easting       = as.numeric(Easting),
    Northing      = as.numeric(Northing),
    PM25_modelled = as.numeric(PM25_modelled),
    Site          = factor(Site),
    UR            = factor(UR),
    Type          = factor(Type)
  )

## 1.1 Missing-value summary
pm25 %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to="variable", values_to="n_missing") %>%
  arrange(desc(n_missing)) %>% print()

## 1.2 Drop rows with missing in analysis variables
pm25_clean <- pm25 %>%
  drop_na(PM25, PM25_modelled,
          Hurs, Psl, Rainfall, SfcWind, Sun, Tas,
          NDVI, PopDen, nearest_road_distance,
          Easting, Northing,
          Year, Month, UR, Type, Site)

## 1.3 Lump rare Type levels (keep top-3, others → “Other”)
pm25_clean <- pm25_clean %>%
  mutate(Type = fct_lump_n(Type, n = 3, other_level = "Other"))

## 1.4 Z-score standardize continuous covariates
cont_vars <- c("Hurs","Psl","Rainfall","SfcWind","Sun",
               "Tas","NDVI","PopDen","nearest_road_distance",
               "PM25_modelled","Easting","Northing")

pm25_ready <- pm25_clean %>%
  mutate(across(all_of(cont_vars), ~ as.numeric(scale(.))))

## 1.5 Assess necessity of log-transform

### 1.5.1 Compute skewness before & after
orig_skew <- e1071::skewness(pm25_ready$PM25,          na.rm = TRUE)
log_skew  <- e1071::skewness(log(pm25_ready$PM25 + 1), na.rm = TRUE)
message(sprintf("Skewness — original: %.3f; log-transformed: %.3f", orig_skew, log_skew))

### 1.5.2 Plot histograms to compare distribution shape
df_skew <- pm25_ready %>%
  select(PM25) %>%
  mutate(log_PM25 = log(PM25 + 1)) %>%
  pivot_longer(c(PM25, log_PM25), names_to = "Scale", values_to = "Value")

ggplot(df_skew, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~ Scale, scales = "free_x",
             labeller = as_labeller(c(PM25 = "Original PM2.5",
                                      log_PM25 = "log(PM2.5 + 1)"))) +
  labs(title = "Distribution: Original vs Log-Transformed PM2.5",
       x = NULL, y = "Count") +
  theme_minimal()

## 1.6 Apply log-transform for downstream models
pm25_ready <- pm25_ready %>%
  mutate(log_PM25 = log(PM25 + 1))

## 1.7 Row-count audit
rows_raw   <- nrow(pm25)
rows_clean <- nrow(pm25_clean)
rows_ready <- nrow(pm25_ready)

cat(sprintf(
  "Row audit — raw: %d | after drop_na: %d | after scaling/log: %d\n",
  rows_raw, rows_clean, rows_ready
))

# 2. Exploratory Data Analysis

## 2.1 Boxplots by Month & Year
ggplot(pm25_ready, aes(Month, PM25)) +
  geom_boxplot(fill="lightblue") +
  labs(title="Monthly PM2.5 Distribution", y=expression(PM[2.5]~(µg/m^3))) +
  theme_minimal()

ggplot(pm25_ready, aes(factor(Year), PM25)) +
  geom_boxplot(fill="lightgreen") +
  labs(title="Yearly PM2.5 Distribution", x="Year", y=expression(PM[2.5]~(µg/m^3))) +
  theme_minimal()

## 2.2 Spatial mean PM2.5 per site
site_mean <- pm25_ready %>%
  group_by(Site, Easting, Northing) %>%
  summarise(mean_PM25 = mean(PM25), .groups="drop")

ggplot(site_mean, aes(Easting, Northing, color=mean_PM25)) +
  geom_point(size=3) +
  scale_color_viridis_c(name=expression(Mean~PM[2.5])) +
  labs(title="Spatial Mean PM2.5") +
  theme_void() + theme(legend.position="right")

## 2.3 Scatterplot matrix
vars_all  <- c(cont_vars, "PM25")
pm_small  <- pm25_ready %>% select(all_of(vars_all))
GGally::ggpairs(
  pm_small,
  diag  = list(continuous = wrap("densityDiag")),
  upper = list(continuous = wrap("points", alpha=0.3, size=0.5)),
  lower = list(continuous = wrap("smooth", method="loess", se=FALSE))
) + theme_minimal()

## 2.4 Correlation heatmap
corr_mat <- cor(pm_small, use="pairwise.complete.obs")
corrplot(corr_mat, method = "color", type   = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "white", number.cex = 0.7)

corr_tbl <- corr_mat %>%
  round(2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  tibble::as_tibble()

print(corr_tbl)

corr_long <- corr_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation")

corr_with_PM25 <- corr_long %>%
  dplyr::filter(Var2 == "PM25") %>%
  dplyr::arrange(dplyr::desc(abs(Correlation)))

print(corr_with_PM25)

## 2.5 Numerical summaries by Month & Year
pm25_ready %>%
  group_by(Month) %>%
  summarise(mean_PM25=mean(PM25), sd_PM25=sd(PM25),
            median_PM25=median(PM25), IQR_PM25=IQR(PM25), .groups="drop") %>% print()
pm25_ready %>%
  group_by(Year) %>%
  summarise(mean_PM25=mean(PM25), sd_PM25=sd(PM25),
            median_PM25=median(PM25), IQR_PM25=IQR(PM25), .groups="drop") %>% print()

# 3. Feature Screening & Selection

## 3.1 Correlation with response
corr_with_y <- corr_mat[,"PM25"] %>%
  enframe(name="var", value="corr") %>%
  arrange(desc(abs(corr)))
print(corr_with_y)

## 3.2 Variance Inflation Factors on full set
vif_full <- lm(PM25 ~ Hurs + Psl + Rainfall + SfcWind + Sun +
                 Tas + NDVI + PopDen + nearest_road_distance +
                 PM25_modelled + Easting + Northing + Year +
                 factor(Month) + UR + Type,
               data = pm25_ready)
print(vif(vif_full))

## 3.3 Identify and drop highly collinear pairs
high_corr_pairs <- which(abs(corr_mat) > 0.8 & abs(corr_mat) < 1, arr.ind=TRUE)
print(high_corr_pairs)
drops <- intersect(c("Sun","Northing"), colnames(pm_small))
pm25_reduced <- pm25_ready %>% select(-all_of(drops))

### 3.3.1 Collinearity drop audit — rows/cols change and dropped variables
rows_reduced <- nrow(pm25_reduced)
cols_before  <- colnames(pm25_ready)
cols_after   <- colnames(pm25_reduced)
dropped_cols <- setdiff(cols_before, cols_after)

num_pairs <- if (length(high_corr_pairs)) nrow(as.data.frame(high_corr_pairs)) else 0L

cat(sprintf(
  "Collinearity audit — rows: %d (unchanged); cols: %d -> %d; dropped: %s\n",
  rows_reduced,
  length(cols_before), length(cols_after),
  ifelse(length(dropped_cols) == 0, "none", paste(dropped_cols, collapse = ", "))
))
cat(sprintf("High-correlation pairs detected (|r| > 0.80): %d\n", num_pairs))

## 3.4 Verify no remaining collinearity
vars2    <- c(setdiff(cont_vars, drops), "PM25")
corr_mat2<- pm25_reduced %>% select(all_of(vars2)) %>%
  cor(use="pairwise.complete.obs")
print(which(abs(corr_mat2) > 0.8 & abs(corr_mat2) < 1, arr.ind=TRUE))
vif_reduced <- lm(PM25 ~ Hurs + Psl + Rainfall + SfcWind +
                    Tas + NDVI + PopDen + nearest_road_distance +
                    PM25_modelled + Easting + Year +
                    factor(Month) + UR + Type,
                  data = pm25_reduced)
print(vif(vif_reduced))

## 3.5 Automatic feature selection via stepwise AIC (log-scale)
null_mod  <- lm(log_PM25 ~ 1, data = pm25_reduced)
full_mod  <- lm(log_PM25 ~ Hurs + Psl + Rainfall + SfcWind +
                  Tas + NDVI + PopDen + nearest_road_distance +
                  PM25_modelled + Easting + Year +
                  factor(Month) + UR + Type,
                data = pm25_reduced)
step_mod  <- step(null_mod, scope=list(lower=null_mod, upper=full_mod),
                  direction="both", trace=FALSE)
best_formula <- formula(step_mod)
print(best_formula)
print(AIC(step_mod))

## 3.6 Global factor levels for harmonization
month_lv <- levels(pm25_reduced$Month)
type_lv  <- levels(pm25_reduced$Type)
ur_lv    <- levels(pm25_reduced$UR)

# 4. Site-based Train/Test split (80% sites train, 20% test)
set.seed(2025)

all_sites  <- unique(pm25_reduced$Site)
n_sites    <- length(all_sites)
n_test     <- max(1L, floor(0.2 * n_sites))
test_sites <- sample(all_sites, n_test)
train_data <- pm25_reduced %>% filter(!Site %in% test_sites)
test_data  <- pm25_reduced %>% filter( Site %in% test_sites)

cat(sprintf(
  "[Split] %d sites total | %d train sites, %d test sites\n", 
  n_sites, n_sites - length(test_sites), length(test_sites)
))
cat(sprintf(
  "[Split] %d rows in train, %d rows in test\n", 
  nrow(train_data), nrow(test_data)
))
stopifnot(length(intersect(train_data$Site, test_data$Site)) == 0)

use_aligned_intervals <- TRUE
alpha_pi <- 0.05

conformal_q <- function(y_true, y_hat, alpha = 0.05) {
  err <- abs(y_true - y_hat)
  err <- err[is.finite(err)]
  n   <- length(err)
  if (n == 0) return(NA_real_)
  k <- ceiling((n + 1) * (1 - alpha))
  k <- max(1L, min(k, n))
  sort(err, na.last = NA)[k]
}

apply_q <- function(y_hat, q) {
  if (!is.finite(q)) {
    tibble::tibble(pred = y_hat, lwr = NA_real_, upr = NA_real_)
  } else {
    tibble::tibble(
      pred = y_hat,
      lwr  = pmax(y_hat - q, 0),
      upr  = y_hat + q
    )
  }
}

# 5. Linear Model (log-scale)
set.seed(2025)
lm_model <- lm(best_formula, data = train_data)
summary(lm_model)

## 5.1 Duan smearing factor (from log-scale residuals on training LM)
resid_log_lm <- residuals(lm_model)
smear_factor <- mean(exp(resid_log_lm), na.rm = TRUE)
message(sprintf("LM Duan smearing factor = %.4f", smear_factor))

## 5.2 Training In-sample metrics (LM)
fit_log_tr <- predict(lm_model)
S_train    <- mean(exp(resid(lm_model)), na.rm = TRUE)
pred_tr    <- S_train * exp(fit_log_tr) - 1

if (use_aligned_intervals) {
  # Split-conformal: same interval logic for all models
  q_lm_train <- conformal_q(train_data$PM25, pred_tr, alpha = alpha_pi)
  pi_tr      <- apply_q(pred_tr, q_lm_train)
  lm_train_insample <- tibble(
    RMSE     = Metrics::rmse(train_data$PM25, pi_tr$pred),
    MAE      = Metrics::mae( train_data$PM25, pi_tr$pred),
    R2       = cor(train_data$PM25, pi_tr$pred, use = "complete.obs")^2,
    Coverage = mean(train_data$PM25 >= pi_tr$lwr & train_data$PM25 <= pi_tr$upr, na.rm = TRUE),
    Width    = mean(pi_tr$upr - pi_tr$lwr, na.rm = TRUE)
  )
} else {
  # Fallback: classical parametric PI on log scale, then smear back
  pi_train <- suppressWarnings(tryCatch(
    predict(lm_model, interval = "prediction", level = 1 - alpha_pi),
    error = function(e) NULL
  ))
  if (!is.null(pi_train)) {
    fit_log <- pi_train[, "fit"]; lwr_log <- pi_train[, "lwr"]; upr_log <- pi_train[, "upr"]
    pred_bt <- S_train * exp(fit_log) - 1
    lwr_bt  <- S_train * exp(lwr_log) - 1
    upr_bt  <- S_train * exp(upr_log) - 1
    lm_train_insample <- tibble(
      RMSE     = Metrics::rmse(train_data$PM25, pred_bt),
      MAE      = Metrics::mae( train_data$PM25, pred_bt),
      R2       = cor(train_data$PM25, pred_bt, use = "complete.obs")^2,
      Coverage = mean(train_data$PM25 >= lwr_bt & train_data$PM25 <= upr_bt, na.rm = TRUE),
      Width    = mean(upr_bt - lwr_bt, na.rm = TRUE)
    )
  } else {
    # If interval prediction failed, compute point metrics only
    lm_train_insample <- tibble(
      RMSE = Metrics::rmse(train_data$PM25, pred_tr),
      MAE  = Metrics::mae( train_data$PM25, pred_tr),
      R2   = cor(train_data$PM25, pred_tr, use = "complete.obs")^2,
      Coverage = NA_real_, Width = NA_real_
    )
  }
}

print(lm_train_insample)

# 6. Leave-One-Site-Out Cross-Validation (LM)

## 6.1 Harmonise factor levels in train_data
train_data <- train_data %>%
  mutate(
    Month = factor(as.character(Month), levels = month_lv),
    Type  = factor(as.character(Type),  levels = type_lv),
    UR    = factor(as.character(UR),    levels = ur_lv)
  )
train_sites_lm <- unique(train_data$Site)

set.seed(2025)

cv_results <- purrr::map_df(train_sites_lm, function(site_heldout) {
  tryCatch({
    tr <- dplyr::filter(train_data, Site != site_heldout)
    te <- dplyr::filter(train_data, Site == site_heldout)
    
    if (nrow(te) < 2) {
      message("Skipping site ", site_heldout, " (n_te = ", nrow(te), ")")
      return(tibble(
        Site     = site_heldout,
        RMSE     = NA_real_, MAE = NA_real_, R2 = NA_real_,
        Coverage = NA_real_, Width = NA_real_
      ))
    }
    
    # Fit fold model with stepwise-selected formula
    fit <- lm(best_formula, data = tr)
    
    # Ensure factor levels for predict
    vars_in_fit <- all.vars(terms(fit))
    xl <- list()
    if ("Month" %in% vars_in_fit) xl$Month <- month_lv
    if ("Type"  %in% vars_in_fit) xl$Type  <- type_lv
    if ("UR"    %in% vars_in_fit) xl$UR    <- ur_lv
    if (length(xl)) fit$xlevels <- modifyList(fit$xlevels, xl)
    
    # Fold-specific Duan smearing factor (based on log-scale residuals)
    S_fold <- mean(exp(residuals(fit)), na.rm = TRUE)
    
    # Point predictions on train/test (log -> original + smearing)
    fit_log_tr <- predict(fit, newdata = tr)
    pred_tr    <- (exp(fit_log_tr) - 1) * S_fold
    
    fit_log_te <- predict(fit, newdata = te)
    pred_te    <- (exp(fit_log_te) - 1) * S_fold
    
    # Intervals
    if (use_aligned_intervals) {
      # Conformal: symmetric residual quantile on the fold's training set
      q_fold <- conformal_q(tr$PM25, pred_tr, alpha_pi)
      pi_te  <- apply_q(pred_te, q_fold)
    } else {
      # Parametric PI on log-scale, then smearing back — with robust guards
      if (df.residual(fit) <= 0) {
        warning(sprintf("[LM][LOSO CV][Site=%s] df.residual <= 0; cannot form parametric PI. Using point prediction only.",
                        site_heldout))
        pi_te <- tibble(pred = pred_te, lwr = NA_real_, upr = NA_real_)
      } else {
        err_msg <- NULL
        pi_int <- suppressWarnings(
          tryCatch(
            predict(fit, newdata = te, interval = "prediction", level = 1 - alpha_pi),
            error = function(e) { err_msg <<- conditionMessage(e); NULL }
          )
        )
        if (is.null(pi_int)) {
          warning(sprintf("[LM][LOSO CV][Site=%s] Interval prediction failed — %s. Falling back to point prediction only.",
                          site_heldout, err_msg))
          pi_te <- tibble(pred = pred_te, lwr = NA_real_, upr = NA_real_)
        } else {
          lwr_te <- (exp(pi_int[, "lwr"]) - 1) * S_fold
          upr_te <- (exp(pi_int[, "upr"]) - 1) * S_fold
          pi_te  <- tibble(pred = pred_te, lwr = lwr_te, upr = upr_te)
        }
      }
    }
    
    obs  <- te$PM25
    good <- is.finite(pi_te$pred) & is.finite(obs)
    
    tibble(
      Site     = site_heldout,
      RMSE     = sqrt(mean((pi_te$pred[good] - obs[good])^2, na.rm = TRUE)),
      MAE      = mean(abs(pi_te$pred[good] - obs[good]), na.rm = TRUE),
      R2       = if (sum(good) >= 2) cor(pi_te$pred[good], obs[good])^2 else NA_real_,
      Coverage = mean(obs[good] >= pi_te$lwr[good] & obs[good] <= pi_te$upr[good], na.rm = TRUE),
      Width    = mean(pi_te$upr[good] - pi_te$lwr[good], na.rm = TRUE)
    )
  }, error = function(e) {
    # Outer catch: never let a single fold crash the whole CV
    warning(sprintf("[LM][LOSO CV][Site=%s] Fold failed — %s",
                    site_heldout, conditionMessage(e)))
    tibble(
      Site     = site_heldout,
      RMSE     = NA_real_, MAE = NA_real_, R2 = NA_real_,
      Coverage = NA_real_, Width = NA_real_
    )
  })
})

site_n <- train_data %>% count(Site, name = "n_te")

cv_diag_lm <- cv_results %>%
  dplyr::select(Site, Coverage, Width) %>%
  dplyr::mutate(PI_failed = is.na(Coverage)) %>%
  dplyr::left_join(site_n, by = "Site") %>%
  dplyr::arrange(dplyr::desc(PI_failed), n_te)

print(cv_diag_lm, n = Inf)
cat(sprintf("[LM][LOSO] PI failed in %d / %d folds\n",
            sum(cv_diag_lm$PI_failed, na.rm = TRUE),
            nrow(cv_diag_lm)))

## 6.2 Summarise CV metrics (LM)
cv_summary <- cv_results %>%
  summarise(across(c(RMSE, MAE, R2, Coverage, Width), ~ mean(.x, na.rm = TRUE)))
print(cv_summary)

# 7. Evaluate on hold-out test data

## 7.1 Harmonise factor levels in test set
test_std <- test_data %>%
  mutate(
    Month = factor(as.character(Month), levels = month_lv),
    Type  = factor(as.character(Type),  levels = type_lv),
    UR    = factor(as.character(UR),    levels = ur_lv)
  )

## 7.2 Hold-out metrics (LM)
fit_log_test <- predict(lm_model, newdata = test_std)
S_train      <- mean(exp(resid(lm_model)), na.rm = TRUE)
pred_test    <- (exp(fit_log_test) - 1) * S_train

if (use_aligned_intervals) {
  # Use training residual quantiles (conformal)
  fit_log_tr <- predict(lm_model)
  pred_tr    <- (exp(fit_log_tr) - 1) * S_train
  q_train    <- conformal_q(train_data$PM25, pred_tr, alpha_pi)
  pi_test    <- apply_q(pred_test, q_train)
} else {
  # Fallback: parametric PI with error capture
  err_msg <- NULL
  pi_int <- suppressWarnings(
    tryCatch(
      predict(lm_model, newdata = test_std, interval = "prediction", level = 1 - alpha_pi),
      error = function(e) { err_msg <<- conditionMessage(e); NULL }
    )
  )
  if (is.null(pi_int)) {
    warning(sprintf("[LM][Holdout] Interval prediction failed — %s. Using point prediction only.", err_msg))
    pi_test <- tibble(pred = pred_test, lwr = NA_real_, upr = NA_real_)
  } else {
    lwr <- (exp(pi_int[, "lwr"]) - 1) * S_train
    upr <- (exp(pi_int[, "upr"]) - 1) * S_train
    pi_test <- tibble(pred = pred_test, lwr = lwr, upr = upr)
  }
}

holdout_metrics <- tibble(
  RMSE     = Metrics::rmse(test_std$PM25,  pi_test$pred),
  MAE      = Metrics::mae( test_std$PM25,  pi_test$pred),
  R2       = cor(test_std$PM25,  pi_test$pred, use = "complete.obs")^2,
  Coverage = mean(test_std$PM25 >= pi_test$lwr & test_std$PM25 <= pi_test$upr, na.rm = TRUE),
  Width    = mean(pi_test$upr - pi_test$lwr, na.rm = TRUE)
)
print(holdout_metrics)

# 8. Diagnostics (LM)
set.seed(2025)

## 8.1 Check linear model assumptions
check_model(lm_model)

## 8.2 Durbin–Watson test for autocorrelation
durbinWatsonTest(lm_model)

## 8.3 Breusch–Pagan test for heteroskedasticity
bptest(lm_model)

## 8.4 Variance Inflation Factors for multicollinearity
print(vif(lm_model))

## 8.5 Coefficient estimates with 95% CI
coefs <- tidy(lm_model, conf.int = TRUE)
print(coefs)

## 8.6 Robust / Cluster-robust SE for coefficient inference
.lm_mf <- model.frame(lm_model)
.idx    <- as.integer(rownames(.lm_mf))
.cluster_vec <- train_data$Site[.idx]

vc_hc3 <- sandwich::vcovHC(lm_model, type = "HC3")
robust_hc3 <- lmtest::coeftest(lm_model, vcov = vc_hc3)
print(robust_hc3)

vc_cl  <- sandwich::vcovCL(lm_model, cluster = .cluster_vec)
robust_cl <- lmtest::coeftest(lm_model, vcov = vc_cl)
print(robust_cl)

# 9. Residual diagnostics plots (LM)
resid_df <- augment(lm_model)

## 9.1 Histogram of residuals
ggplot(resid_df, aes(x = .resid)) +
  geom_histogram(bins = 30, fill = "gray", color = "white") +
  labs(title = "Histogram of Residuals", x = "Residual", y = "Count") +
  theme_minimal()

## 9.2 Residuals vs. fitted values
ggplot(resid_df, aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs. Fitted", x = "Fitted values", y = "Residuals") +
  theme_minimal()

## 9.3 Q–Q plot of residuals
ggplot(resid_df, aes(sample = .resid)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q–Q Plot of Residuals") +
  theme_minimal()

# 10. Temporal Autocorrelation in Residuals
Acf(resid(lm_model), main = "ACF of Model Residuals")

# 11. LASSO REGRESSION (log-scale)

## 11.0 Use best_formula for LASSO too
lasso_formula <- best_formula

## 11.1 Harmonise factor levels
train_data_lasso <- train_data %>%
  mutate(
    Month = factor(as.character(Month), levels = month_lv),
    Type  = factor(as.character(Type),  levels = type_lv),
    UR    = factor(as.character(UR),    levels = ur_lv)
  )

## 11.2 Design matrices for glmnet
X_train_lasso <- model.matrix(lasso_formula, data = train_data_lasso)[ , -1, drop = FALSE]
y_train_lasso <- train_data_lasso$log_PM25

## 11.3 LOSO-style folds for tuning
train_sites_lasso  <- unique(train_data_lasso$Site)
foldid_train_lasso <- purrr::map_int(train_data_lasso$Site, ~ which(train_sites_lasso == .x))

## 11.4 cv.glmnet with LOSO folds
set.seed(2025)
cv_lasso <- cv.glmnet(
  x = X_train_lasso, y = y_train_lasso,
  alpha = 1, family = "gaussian",
  foldid = foldid_train_lasso,
  standardize = FALSE,
  nlambda = 200
)
lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se
message(sprintf("LASSO: lambda.min = %.6f; lambda.1se = %.6f", lambda_min, lambda_1se))

## 11.5 Final LASSO on full training (lambda_1se)
set.seed(2025)
lasso_fit <- glmnet(
  x = X_train_lasso, y = y_train_lasso,
  alpha = 1, family = "gaussian",
  lambda = lambda_1se,
  standardize = FALSE
)

## 11.6 design alignment (reusable by LASSO/RF)
align_design <- function(Xref, Xnew) {
  ref_names <- colnames(Xref)
  miss <- setdiff(ref_names, colnames(Xnew))
  if (length(miss) > 0) {
    add <- matrix(0, nrow = nrow(Xnew), ncol = length(miss))
    colnames(add) <- miss
    Xnew <- cbind(Xnew, add)
  }
  extra <- setdiff(colnames(Xnew), ref_names)
  if (length(extra) > 0) Xnew <- Xnew[, !(colnames(Xnew) %in% extra), drop = FALSE]
  Xnew <- Xnew[, ref_names, drop = FALSE]
  Xnew
}

## 11.7 Training — In-sample metrics (LASSO)
pred_log_train_lasso <- as.numeric(predict(lasso_fit, newx = X_train_lasso, s = lambda_1se))
pred_pm_train_lasso  <- exp(pred_log_train_lasso) - 1

if (use_aligned_intervals) {
  # Conformal q from training residuals (original scale)
  q_train_lasso  <- conformal_q(train_data_lasso$PM25, pred_pm_train_lasso, alpha_pi)
  pi_train_lasso <- apply_q(pred_pm_train_lasso, q_train_lasso)  # tibble(pred, lwr, upr)
} else {
  # Fallback: symmetric band on log residuals, then back-transform
  fitted_log_train_lasso <- pred_log_train_lasso
  resid_log_train_lasso  <- y_train_lasso - fitted_log_train_lasso
  q975_train_lasso       <- as.numeric(quantile(abs(resid_log_train_lasso), probs = 0.975, na.rm = TRUE))
  lwr_pm_train_lasso     <- exp(pred_log_train_lasso - q975_train_lasso) - 1
  upr_pm_train_lasso     <- exp(pred_log_train_lasso + q975_train_lasso) - 1
  pi_train_lasso <- tibble(pred = pred_pm_train_lasso, lwr = lwr_pm_train_lasso, upr = upr_pm_train_lasso)
}

lasso_train_insample <- tibble(
  RMSE     = Metrics::rmse(train_data_lasso$PM25, pi_train_lasso$pred),
  MAE      = Metrics::mae( train_data_lasso$PM25, pi_train_lasso$pred),
  R2       = cor(train_data_lasso$PM25, pi_train_lasso$pred, use = "complete.obs")^2,
  Coverage = mean(train_data_lasso$PM25 >= pi_train_lasso$lwr & train_data_lasso$PM25 <= pi_train_lasso$upr, na.rm = TRUE),
  Width    = mean(pi_train_lasso$upr - pi_train_lasso$lwr, na.rm = TRUE)
)
print(lasso_train_insample)

## 11.8 Training — LOSO CV metrics (LASSO)
set.seed(2025)
lasso_cv_results <- purrr::map_df(train_sites_lasso, function(site_heldout) {
  tr <- dplyr::filter(train_data_lasso, Site != site_heldout)
  te <- dplyr::filter(train_data_lasso, Site == site_heldout)
  
  if (nrow(te) < 2) {
    message("Skipping site ", site_heldout, " in LASSO-CV (n_te=", nrow(te), ")")
    return(tibble(
      Site     = site_heldout,
      RMSE     = NA_real_, MAE = NA_real_, R2 = NA_real_,
      Coverage = NA_real_, Width = NA_real_
    ))
  }
  
  X_tr <- model.matrix(lasso_formula, data = tr)[ , -1, drop = FALSE]
  y_tr <- tr$log_PM25
  X_te <- model.matrix(lasso_formula, data = te)[ , -1, drop = FALSE]
  X_te <- align_design(X_tr, X_te)
  
  # Fit fold model at tuned lambda
  fit_cv <- glmnet(
    x = X_tr, y = y_tr,
    alpha = 1, family = "gaussian",
    lambda = lambda_1se,
    standardize = FALSE
  )
  
  # Test predictions (log → original)
  pred_log_te <- as.numeric(predict(fit_cv, newx = X_te, s = lambda_1se))
  pred_pm_te  <- exp(pred_log_te) - 1
  
  # Training predictions in this fold (for conformal q)
  pred_log_tr <- as.numeric(predict(fit_cv, newx = X_tr, s = lambda_1se))
  pred_pm_tr  <- exp(pred_log_tr) - 1
  
  # Intervals
  if (use_aligned_intervals) {
    q_fold <- conformal_q(tr$PM25, pred_pm_tr, alpha_pi)
    pi_te  <- apply_q(pred_pm_te, q_fold)  # tibble(pred, lwr, upr)
  } else {
    # Fallback: symmetric band on log residuals, then back-transform
    resid_log_tr <- y_tr - pred_log_tr
    q975         <- as.numeric(quantile(abs(resid_log_tr), probs = 0.975, na.rm = TRUE))
    lwr_pm <- exp(pred_log_te - q975) - 1
    upr_pm <- exp(pred_log_te + q975) - 1
    pi_te  <- tibble(pred = pred_pm_te, lwr = lwr_pm, upr = upr_pm)
  }
  
  obs  <- te$PM25
  good <- is.finite(pi_te$pred) & is.finite(obs)
  
  tibble(
    Site     = site_heldout,
    RMSE     = sqrt(mean((pi_te$pred[good] - obs[good])^2, na.rm = TRUE)),
    MAE      = mean(abs(pi_te$pred[good] - obs[good]), na.rm = TRUE),
    R2       = if (sum(good) >= 2) cor(pi_te$pred[good], obs[good])^2 else NA_real_,
    Coverage = mean(obs[good] >= pi_te$lwr[good] & obs[good] <= pi_te$upr[good], na.rm = TRUE),
    Width    = mean(pi_te$upr[good] - pi_te$lwr[good], na.rm = TRUE)
  )
})

lasso_cv_summary <- lasso_cv_results %>%
  summarise(across(c(RMSE, MAE, R2, Coverage, Width), ~ mean(.x, na.rm = TRUE)))
print(lasso_cv_summary)

## 11.9 Hold-out — Test set metrics (LASSO)
X_test_lasso  <- model.matrix(lasso_formula, data = test_std)[ , -1, drop = FALSE]
X_test_lassoA <- align_design(X_train_lasso, X_test_lasso)
pred_log_test_lasso <- as.numeric(predict(lasso_fit, newx = X_test_lassoA, s = lambda_1se))
pred_pm_test_lasso  <- exp(pred_log_test_lasso) - 1

if (use_aligned_intervals) {
  pi_test_lasso <- apply_q(pred_pm_test_lasso, q_train_lasso)
} else {
  lwr_pm_test_lasso <- exp(pred_log_test_lasso - q975_train_lasso) - 1
  upr_pm_test_lasso <- exp(pred_log_test_lasso + q975_train_lasso) - 1
  pi_test_lasso <- tibble(pred = pred_pm_test_lasso, lwr = lwr_pm_test_lasso, upr = upr_pm_test_lasso)
}

lasso_holdout <- tibble(
  RMSE     = Metrics::rmse(test_std$PM25, pi_test_lasso$pred),
  MAE      = Metrics::mae( test_std$PM25, pi_test_lasso$pred),
  R2       = cor(test_std$PM25, pi_test_lasso$pred, use = "complete.obs")^2,
  Coverage = mean(test_std$PM25 >= pi_test_lasso$lwr & test_std$PM25 <= pi_test_lasso$upr, na.rm = TRUE),
  Width    = mean(pi_test_lasso$upr - pi_test_lasso$lwr, na.rm = TRUE)
)
print(lasso_holdout)

# 12. RANDOM FOREST (log-scale)

## 12.0 Use best_formula for RF too
rf_formula <- best_formula

## 12.1 Harmonise factor levels
train_data_rf <- train_data %>%
  mutate(
    Month = factor(as.character(Month), levels = month_lv),
    Type  = factor(as.character(Type),  levels = type_lv),
    UR    = factor(as.character(UR),    levels = ur_lv)
  )

## 12.2 Design space (matrix -> data.frame for ranger)
X_train_rf <- model.matrix(rf_formula, data = train_data_rf)[ , -1, drop = FALSE]
df_train_rf <- data.frame(log_PM25 = train_data_rf$log_PM25, X_train_rf)

## 12.3 Global hyperparameter tuning via OOB (log-scale RMSE)
p_rf <- ncol(X_train_rf)
mtry_grid <- unique(pmax(1, pmin(p_rf, round(c(sqrt(p_rf), 2*sqrt(p_rf), p_rf/3, p_rf/2)))))
min_node_grid <- c(5, 10)
num_trees <- 800

grid <- expand.grid(mtry = mtry_grid, min.node.size = min_node_grid)
grid$RMSE_log <- NA_real_

set.seed(2025)
for (i in seq_len(nrow(grid))) {
  set.seed(2025)
  fit_try <- ranger(
    dependent.variable.name = "log_PM25",
    data = df_train_rf,
    num.trees = num_trees,
    mtry = grid$mtry[i],
    min.node.size = grid$min.node.size[i],
    oob.error = TRUE, seed = 2025
  )
  grid$RMSE_log[i] <- sqrt(fit_try$prediction.error)
}

best_idx <- which.min(grid$RMSE_log)
best_mtry <- grid$mtry[best_idx]
best_min_node <- grid$min.node.size[best_idx]
message(sprintf("RF tuned: mtry=%d, min.node.size=%d (OOB RMSE_log=%.4f)",
                best_mtry, best_min_node, grid$RMSE_log[best_idx]))

## 12.4 Final RF on full training
set.seed(2025)
rf_fit <- ranger(
  dependent.variable.name = "log_PM25",
  data = df_train_rf,
  num.trees = num_trees,
  mtry = best_mtry,
  min.node.size = best_min_node,
  importance = "impurity",
  oob.error = TRUE, seed = 2025
)

### 12.4.1 Final RF (quantile forest for intervals)
q_probs <- c(0.025, 0.5, 0.975)

set.seed(2025)
rf_fit_q <- ranger(
  dependent.variable.name = "log_PM25",
  data = df_train_rf,
  num.trees = num_trees,
  mtry = best_mtry,
  min.node.size = best_min_node,
  quantreg = TRUE,
  oob.error = TRUE,
  seed = 2025
)

## 12.5 Training — In-sample metrics
q50_train <- predict(
  rf_fit_q, data = df_train_rf, type = "quantiles", quantiles = 0.5
)$predictions
if (is.null(dim(q50_train))) q50_train <- matrix(q50_train, ncol = 1)
pred_pm_train_rf <- exp(q50_train[, 1]) - 1

if (isTRUE(use_aligned_intervals)) {
  if (!exists("q_train_rf") || !is.numeric(q_train_rf) || !is.finite(q_train_rf) || length(q_train_rf) != 1) {
    if (!exists("rf_mean_full")) {
      rf_mean_full <- ranger(
        dependent.variable.name = "log_PM25",
        data = df_train_rf,
        num.trees = num_trees,
        mtry = best_mtry,
        min.node.size = best_min_node,
        oob.error = TRUE,
        seed = 2025
      )
    }
    pred_pm_oob_full <- exp(rf_mean_full$predictions) - 1
    q_train_rf <- conformal_q(train_data_rf$PM25, pred_pm_oob_full, alpha_pi)
  }
  pi_train_rf <- apply_q(pred_pm_train_rf, q_train_rf)
  
} else {
  q_train_all <- predict(
    rf_fit_q, data = df_train_rf, type = "quantiles", quantiles = c(0.025, 0.975)
  )$predictions
  if (is.null(dim(q_train_all))) q_train_all <- matrix(q_train_all, ncol = 2)
  lwr_pm_train_rf <- exp(q_train_all[, 1]) - 1
  upr_pm_train_rf <- exp(q_train_all[, 2]) - 1
  pi_train_rf <- tibble::tibble(pred = pred_pm_train_rf, lwr = lwr_pm_train_rf, upr = upr_pm_train_rf)
}

rf_train_insample <- tibble::tibble(
  RMSE     = Metrics::rmse(train_data_rf$PM25,  pi_train_rf$pred),
  MAE      = Metrics::mae( train_data_rf$PM25,  pi_train_rf$pred),
  R2       = cor(train_data_rf$PM25,            pi_train_rf$pred, use = "complete.obs")^2,
  Coverage = mean(train_data_rf$PM25 >= pi_train_rf$lwr & train_data_rf$PM25 <= pi_train_rf$upr, na.rm = TRUE),
  Width    = mean(pi_train_rf$upr - pi_train_rf$lwr, na.rm = TRUE)
)
print(rf_train_insample)

## 12.6 Training — LOSO CV metrics (RF, quantile forest)
train_sites_rf <- unique(train_data_rf$Site)
set.seed(2025)
rf_cv_results <- purrr::map_df(train_sites_rf, function(site_heldout) {
  tr <- dplyr::filter(train_data_rf, Site != site_heldout)
  te <- dplyr::filter(train_data_rf, Site == site_heldout)
  
  if (nrow(te) < 2) {
    message("Skipping site ", site_heldout, " in RF-CV (n_te=", nrow(te), ")")
    return(tibble(
      Site     = site_heldout,
      RMSE     = NA_real_, MAE = NA_real_, R2 = NA_real_,
      Coverage = NA_real_, Width = NA_real_
    ))
  }
  
  X_tr <- model.matrix(rf_formula, data = tr)[, -1, drop = FALSE]
  X_te <- model.matrix(rf_formula, data = te)[, -1, drop = FALSE]
  X_te <- align_design(X_tr, X_te)
  df_tr <- data.frame(log_PM25 = tr$log_PM25, X_tr)
  df_te <- data.frame(X_te)
  
  err_fit <- NULL
  fit_cv_q <- tryCatch(
    ranger(
      dependent.variable.name = "log_PM25",
      data            = df_tr,
      num.trees       = num_trees,
      mtry            = best_mtry,
      min.node.size   = best_min_node,
      quantreg        = TRUE,
      seed            = 2025
    ),
    error = function(e) { err_fit <<- conditionMessage(e); NULL }
  )
  if (is.null(fit_cv_q)) {
    warning(sprintf("[RF][LOSO CV][Site=%s] Quantile forest fit failed — %s.", site_heldout, err_fit))
    return(tibble(
      Site     = site_heldout,
      RMSE     = NA_real_, MAE = NA_real_, R2 = NA_real_,
      Coverage = NA_real_, Width = NA_real_
    ))
  }
  
  err_fit_mean <- NULL
  fit_cv_mean <- tryCatch(
    ranger(
      dependent.variable.name = "log_PM25",
      data            = df_tr,
      num.trees       = num_trees,
      mtry            = best_mtry,
      min.node.size   = best_min_node,
      oob.error       = TRUE,
      seed            = 2025
    ),
    error = function(e) { err_fit_mean <<- conditionMessage(e); NULL }
  )
  if (is.null(fit_cv_mean)) {
    warning(sprintf("[RF][LOSO CV][Site=%s] Mean forest (OOB) fit failed — %s.", site_heldout, err_fit_mean))
    return(tibble(
      Site     = site_heldout,
      RMSE     = NA_real_, MAE = NA_real_, R2 = NA_real_,
      Coverage = NA_real_, Width = NA_real_
    ))
  }
  
  pred_pm_oob <- exp(fit_cv_mean$predictions) - 1
  q_fold_rf   <- conformal_q(tr$PM25, pred_pm_oob, alpha_pi)
  
  err_qte <- NULL
  q_te <- tryCatch(
    predict(fit_cv_q, data = df_te, type = "quantiles", quantiles = q_probs)$predictions,
    error = function(e) { err_qte <<- conditionMessage(e); NULL }
  )
  if (is.null(q_te)) {
    warning(sprintf("[RF][LOSO CV][Site=%s] Test quantile prediction failed — %s.", site_heldout, err_qte))
    return(tibble(
      Site     = site_heldout,
      RMSE     = NA_real_, MAE = NA_real_, R2 = NA_real_,
      Coverage = NA_real_, Width = NA_real_
    ))
  }
  if (is.null(dim(q_te))) q_te <- matrix(q_te, ncol = length(q_probs))
  colnames(q_te) <- c("q025", "q50", "q975")
  
  pred_pm_te <- exp(q_te[, "q50"]) - 1
  
  if (use_aligned_intervals) {
    pi_te <- apply_q(pred_pm_te, q_fold_rf)
  } else {
    lwr_pm <- exp(q_te[, "q025"]) - 1
    upr_pm <- exp(q_te[, "q975"]) - 1
    pi_te  <- tibble(pred = pred_pm_te, lwr = lwr_pm, upr = upr_pm)
  }
  
  obs  <- te$PM25
  good <- is.finite(pi_te$pred) & is.finite(obs)
  
  tibble(
    Site     = site_heldout,
    RMSE     = sqrt(mean((pi_te$pred[good] - obs[good])^2, na.rm = TRUE)),
    MAE      = mean(abs(pi_te$pred[good] - obs[good]), na.rm = TRUE),
    R2       = if (sum(good) >= 2) cor(pi_te$pred[good], obs[good])^2 else NA_real_,
    Coverage = mean(obs[good] >= pi_te$lwr[good] & obs[good] <= pi_te$upr[good], na.rm = TRUE),
    Width    = mean(pi_te$upr[good] - pi_te$lwr[good], na.rm = TRUE)
  )
})

rf_cv_summary <- rf_cv_results %>%
  dplyr::summarise(
    dplyr::across(c(RMSE, MAE, R2, Coverage, Width), ~ mean(.x, na.rm = TRUE))
  )
print(rf_cv_summary)

## 12.7 Hold-out — Test set metrics (RF, quantile forest)
rf_mean_full <- ranger(
  dependent.variable.name = "log_PM25",
  data          = df_train_rf,
  num.trees     = num_trees,
  mtry          = best_mtry,
  min.node.size = best_min_node,
  oob.error     = TRUE,
  seed          = 2025
)
pred_pm_oob_full <- exp(rf_mean_full$predictions) - 1
q_train_rf <- conformal_q(train_data_rf$PM25, pred_pm_oob_full, alpha_pi)
if (!is.finite(q_train_rf)) q_train_rf <- NA_real_

X_test_rf  <- model.matrix(rf_formula, data = test_std)[ , -1, drop = FALSE]
X_test_rfA <- align_design(X_train_rf, X_test_rf)
df_test_rf <- data.frame(X_test_rfA)

q_test_rf <- tryCatch(
  predict(rf_fit_q, data = df_test_rf, type = "quantiles", quantiles = 0.5)$predictions,
  error = function(e) { warning(sprintf("[RF][Holdout] q50 prediction failed — %s", conditionMessage(e))); NULL }
)
if (is.null(q_test_rf)) {
  pred_pm_test_rf <- rep(NA_real_, nrow(df_test_rf))
} else {
  if (is.null(dim(q_test_rf))) q_test_rf <- matrix(q_test_rf, ncol = 1)
  pred_pm_test_rf <- exp(q_test_rf[,1]) - 1
}

if (use_aligned_intervals) {
  pi_test_rf <- apply_q(pred_pm_test_rf, q_train_rf)
} else {
  pred_q_test <- tryCatch(
    predict(rf_fit_q, data = df_test_rf, type = "quantiles", quantiles = c(0.025, 0.975))$predictions,
    error = function(e) { warning(sprintf("[RF][Holdout] interval prediction failed — %s", conditionMessage(e))); NULL }
  )
  if (is.null(pred_q_test)) {
    pi_test_rf <- tibble(pred = pred_pm_test_rf, lwr = NA_real_, upr = NA_real_)
  } else {
    if (is.null(dim(pred_q_test))) pred_q_test <- matrix(pred_q_test, ncol = 2)
    lwr_pm_test_rf <- exp(pred_q_test[,1]) - 1
    upr_pm_test_rf <- exp(pred_q_test[,2]) - 1
    pi_test_rf <- tibble(pred = pred_pm_test_rf, lwr = lwr_pm_test_rf, upr = upr_pm_test_rf)
  }
}

rf_holdout <- tibble(
  RMSE     = Metrics::rmse(test_std$PM25, pi_test_rf$pred),
  MAE      = Metrics::mae( test_std$PM25, pi_test_rf$pred),
  R2       = cor(test_std$PM25,  pi_test_rf$pred, use = "complete.obs")^2,
  Coverage = mean(test_std$PM25 >= pi_test_rf$lwr & test_std$PM25 <= pi_test_rf$upr, na.rm = TRUE),
  Width    = mean(pi_test_rf$upr - pi_test_rf$lwr, na.rm = TRUE)
)
print(rf_holdout)

## 12.8 RF tuning visualization

### 12.8.1 Ensure `grid` contains mtry, min.node.size, and RMSE_log
stopifnot(all(c("mtry","min.node.size","RMSE_log") %in% names(grid)))

### 12.8.2 Mark the best hyperparameter combination
grid$IsBest <- with(grid, mtry == best_mtry & min.node.size == best_min_node)

### 12.8.3 Line plot: OOB RMSE vs mtry across different min.node.size values
p_line <- ggplot(grid, aes(x = mtry, y = RMSE_log, group = factor(min.node.size))) +
  geom_line(aes(linetype = factor(min.node.size))) +
  geom_point(size = 2) +
  geom_point(data = subset(grid, IsBest), color = "red", size = 3) +
  geom_text(
    data = subset(grid, IsBest),
    aes(label = sprintf("best: mtry=%d\nmin.node=%d\nRMSE=%.3f", mtry, min.node.size, RMSE_log)),
    vjust = -1, size = 3.2
  ) +
  labs(
    title = "RF OOB RMSE vs mtry",
    subtitle = "Lines indicate different min.node.size values",
    x = "mtry",
    y = "OOB RMSE (log scale target)",
    linetype = "min.node.size"
  ) +
  theme_minimal(base_size = 12)

print(p_line)

### 12.8.4 Heatmap: OOB RMSE over the mtry × min.node.size grid
p_heat <- ggplot(grid, aes(x = factor(mtry), y = factor(min.node.size), fill = RMSE_log)) +
  geom_tile() +
  geom_point(
    data = subset(grid, IsBest),
    aes(x = factor(mtry), y = factor(min.node.size)),
    shape = 21, stroke = 0.8, size = 4, fill = NA, color = "red"
  ) +
  geom_text(
    data = subset(grid, IsBest),
    aes(label = sprintf("best\n%.3f", RMSE_log)),
    vjust = 1.6, size = 3.1, color = "red"
  ) +
  scale_fill_viridis_c(name = "OOB RMSE") +
  labs(
    title = "RF Tuning Heatmap",
    subtitle = "Lower is better (log-scale target)",
    x = "mtry",
    y = "min.node.size"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_heat)

## 12.9 RF variable importance

### 12.9.1 Impurity-based importance (from rf_fit)
rf_impurity_tbl <- as_tibble(
  sort(rf_fit$variable.importance, decreasing = TRUE),
  rownames = "term"
) %>%
  rename(importance = value) %>%
  mutate(term = forcats::fct_reorder(term, importance))

topN <- 20
rf_impurity_top <- rf_impurity_tbl %>% slice_head(n = topN)

p_rf_impurity <- ggplot(rf_impurity_top, aes(x = term, y = importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = sprintf("Random Forest Variable Importance (Impurity-based) — Top %d", topN),
    x = NULL, y = "Importance (impurity decrease)"
  ) +
  theme_minimal(base_size = 12)

print(p_rf_impurity)

### 12.9.2 Permutation-based importance (refit with same tuned hyperparameters)
set.seed(2025)
rf_fit_perm <- ranger::ranger(
  dependent.variable.name = "log_PM25",
  data = df_train_rf,
  num.trees = num_trees,
  mtry = best_mtry,
  min.node.size = best_min_node,
  importance = "permutation",
  oob.error = TRUE, seed = 2025
)

rf_perm_tbl <- as_tibble(
  sort(rf_fit_perm$variable.importance, decreasing = TRUE),
  rownames = "term"
) %>%
  rename(importance = value) %>%
  mutate(term = forcats::fct_reorder(term, importance))

rf_perm_top <- rf_perm_tbl %>% slice_head(n = topN)

p_rf_perm <- ggplot(rf_perm_top, aes(x = term, y = importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = sprintf("Random Forest Variable Importance (Permutation-based) — Top %d", topN),
    x = NULL, y = "Importance (permutation drop)"
  ) +
  theme_minimal(base_size = 12)

print(p_rf_perm)

# 13. COMPARISON TABLE
add_label <- function(model, split, tbl) {
  tbl %>% mutate(Model = model, Split = split) %>% relocate(Model, Split)
}
metrics_comparison <- bind_rows(
  add_label("Linear Model",  "Training-InSample", lm_train_insample),
  add_label("Linear Model",  "Training-LOSO",     cv_summary),
  add_label("Linear Model",  "Holdout-Test",      holdout_metrics),
  
  add_label("LASSO",         "Training-InSample", lasso_train_insample),
  add_label("LASSO",         "Training-LOSO",     lasso_cv_summary),
  add_label("LASSO",         "Holdout-Test",      lasso_holdout),
  
  add_label("Random Forest", "Training-InSample", rf_train_insample),
  add_label("Random Forest", "Training-LOSO",     rf_cv_summary),
  add_label("Random Forest", "Holdout-Test",      rf_holdout)
)

alpha_pi <- 0.05
K_cal <- 5
set.seed(2025)

cal_sites <- sample(unique(train_data$Site))
folds <- split(cal_sites, cut(seq_along(cal_sites), breaks = K_cal, labels = FALSE))

q_from_res <- function(err, alpha = alpha_pi) {
  err <- err[is.finite(err)]
  n <- length(err)
  if (n == 0) return(NA_real_)
  k <- ceiling((n + 1) * (1 - alpha))
  k <- max(1L, min(k, n))
  sort(err, na.last = NA)[k]
}

strip_factor_wrappers <- function(fml) {
  f <- paste(deparse(fml), collapse = "")
  f <- gsub("factor\\s*\\(([^)]+)\\)", "\\1", f)
  as.formula(f)
}
best_formula_cv_lm <- strip_factor_wrappers(best_formula)

safe_predict_lm <- function(fit, newdata) {
  Terms <- terms(fit)
  
  mf <- model.frame(Terms, newdata, xlev = fit$xlevels, na.action = na.pass)
  X  <- model.matrix(Terms, mf, contrasts.arg = fit$contrasts)
  
  beta <- coef(fit)
  beta[!is.finite(beta)] <- 0
  
  cn <- colnames(X); bn <- names(beta)
  
  miss <- setdiff(cn, bn)
  if (length(miss)) beta <- c(beta, stats::setNames(rep(0, length(miss)), miss))
  
  beta <- beta[cn]
  drop(X %*% beta)
}

abs_err_lm <- c()
for (k in seq_along(folds)) {
  te_sites <- folds[[k]]
  tr <- dplyr::filter(train_data, !Site %in% te_sites)
  te <- dplyr::filter(train_data,  Site %in% te_sites)
  
  if (nrow(tr) < 5L || nrow(te) == 0L) next
  
  tr <- tr %>% mutate(
    Month = factor(as.character(Month), levels = month_lv),
    Type  = factor(as.character(Type),  levels = type_lv),
    UR    = factor(as.character(UR),    levels = ur_lv)
  )
  te <- te %>% mutate(
    Month = factor(as.character(Month), levels = month_lv),
    Type  = factor(as.character(Type),  levels = type_lv),
    UR    = factor(as.character(UR),    levels = ur_lv)
  )
  
  fit <- lm(best_formula_cv_lm, data = tr)
  
  vars_in_fit <- all.vars(terms(fit))
  xl <- list()
  if ("Month" %in% vars_in_fit) xl$Month <- month_lv
  if ("Type"  %in% vars_in_fit) xl$Type  <- type_lv
  if ("UR"    %in% vars_in_fit) xl$UR    <- ur_lv
  if (length(xl)) fit$xlevels <- modifyList(fit$xlevels, xl)
  
  S_fold <- mean(exp(residuals(fit)), na.rm = TRUE)
  
  yhat_te <- safe_predict_lm(fit, te)
  pred_te <- (exp(yhat_te) - 1) * S_fold
  
  abs_err_lm <- c(abs_err_lm, abs(te$PM25 - pred_te))
}
q_cv_lm <- q_from_res(abs_err_lm, alpha_pi)

abs_err_lasso <- c()
for (k in seq_along(folds)) {
  te_sites <- folds[[k]]
  tr <- dplyr::filter(train_data_lasso, !Site %in% te_sites)
  te <- dplyr::filter(train_data_lasso,  Site %in% te_sites)
  
  X_tr <- model.matrix(lasso_formula, data = tr)[, -1, drop = FALSE]
  y_tr <- tr$log_PM25
  X_te <- model.matrix(lasso_formula, data = te)[, -1, drop = FALSE]
  X_te <- align_design(X_tr, X_te)
  
  fit_cv <- glmnet(x = X_tr, y = y_tr, alpha = 1, family = "gaussian",
                   lambda = lambda_1se, standardize = FALSE)
  
  pred_te <- exp(as.numeric(predict(fit_cv, newx = X_te, s = lambda_1se))) - 1
  abs_err_lasso <- c(abs_err_lasso, abs(te$PM25 - pred_te))
}
q_cv_lasso <- q_from_res(abs_err_lasso, alpha_pi)

abs_err_rf <- c()
for (k in seq_along(folds)) {
  te_sites <- folds[[k]]
  tr <- dplyr::filter(train_data_rf, !Site %in% te_sites)
  te <- dplyr::filter(train_data_rf,  Site %in% te_sites)
  
  X_tr <- model.matrix(rf_formula, data = tr)[, -1, drop = FALSE]
  X_te <- model.matrix(rf_formula, data = te)[, -1, drop = FALSE]
  X_te <- align_design(X_tr, X_te)
  df_tr <- data.frame(log_PM25 = tr$log_PM25, X_tr)
  df_te <- data.frame(X_te)
  
  fit_cv_q <- ranger(
    dependent.variable.name = "log_PM25",
    data = df_tr,
    num.trees = num_trees,
    mtry = best_mtry,
    min.node.size = best_min_node,
    quantreg = TRUE,
    seed = 2025
  )
  
  q50_te <- predict(fit_cv_q, data = df_te, type = "quantiles", quantiles = 0.5)$predictions
  if (is.null(dim(q50_te))) q50_te <- matrix(q50_te, ncol = 1)
  pred_te <- exp(q50_te[, 1]) - 1
  
  abs_err_rf <- c(abs_err_rf, abs(te$PM25 - pred_te))
}
q_cv_rf <- q_from_res(abs_err_rf, alpha_pi)

pi_test_lm_cv     <- apply_q(pred_test,            q_cv_lm)
pi_test_lasso_cv  <- apply_q(pred_pm_test_lasso,   q_cv_lasso)
pi_test_rf_cv     <- apply_q(pred_pm_test_rf,      q_cv_rf)

lm_holdout_cv <- tibble(
  RMSE     = Metrics::rmse(test_std$PM25, pi_test_lm_cv$pred),
  MAE      = Metrics::mae( test_std$PM25, pi_test_lm_cv$pred),
  R2       = cor(test_std$PM25, pi_test_lm_cv$pred, use = "complete.obs")^2,
  Coverage = mean(test_std$PM25 >= pi_test_lm_cv$lwr & test_std$PM25 <= pi_test_lm_cv$upr, na.rm = TRUE),
  Width    = mean(pi_test_lm_cv$upr - pi_test_lm_cv$lwr, na.rm = TRUE)
)

lasso_holdout_cv <- tibble(
  RMSE     = Metrics::rmse(test_std$PM25, pi_test_lasso_cv$pred),
  MAE      = Metrics::mae( test_std$PM25, pi_test_lasso_cv$pred),
  R2       = cor(test_std$PM25, pi_test_lasso_cv$pred, use = "complete.obs")^2,
  Coverage = mean(test_std$PM25 >= pi_test_lasso_cv$lwr & test_std$PM25 <= pi_test_lasso_cv$upr, na.rm = TRUE),
  Width    = mean(pi_test_lasso_cv$upr - pi_test_lasso_cv$lwr, na.rm = TRUE)
)

rf_holdout_cv <- tibble(
  RMSE     = Metrics::rmse(test_std$PM25, pi_test_rf_cv$pred),
  MAE      = Metrics::mae( test_std$PM25, pi_test_rf_cv$pred),
  R2       = cor(test_std$PM25, pi_test_rf_cv$pred, use = "complete.obs")^2,
  Coverage = mean(test_std$PM25 >= pi_test_rf_cv$lwr & test_std$PM25 <= pi_test_rf_cv$upr, na.rm = TRUE),
  Width    = mean(pi_test_rf_cv$upr - pi_test_rf_cv$lwr, na.rm = TRUE)
)

metrics_comparison <- bind_rows(
  metrics_comparison,
  add_label("Linear Model",  "Holdout-Test (CV-Conformal)",  lm_holdout_cv),
  add_label("LASSO",         "Holdout-Test (CV-Conformal)",  lasso_holdout_cv),
  add_label("Random Forest", "Holdout-Test (CV-Conformal)",  rf_holdout_cv)
) %>%
  mutate(
    Split = factor(Split, levels = c("Training-InSample","Training-LOSO","Holdout-Test","Holdout-Test (CV-Conformal)")),
    Model = factor(Model, levels = c("Linear Model","LASSO","Random Forest"))
  ) %>%
  arrange(Split, Model) %>%
  select(Model, Split, RMSE, MAE, R2, Coverage, Width)

print(metrics_comparison)

# 14. VISUALIZATION — Side-by-side comparison plot
metrics_long <- metrics_comparison %>%
  tidyr::pivot_longer(
    cols = c(RMSE, MAE, R2, Coverage, Width),
    names_to = "Metric", values_to = "Value"
  ) %>%
  dplyr::mutate(
    Split  = factor(Split, levels = c("Training-InSample","Training-LOSO","Holdout-Test","Holdout-Test (CV-Conformal)")),
    Model  = factor(Model, levels = c("Linear Model","LASSO","Random Forest")),
    Metric = factor(Metric, levels = c("RMSE","MAE","R2","Coverage","Width")),
    MetricLabel = dplyr::recode(
      Metric,
      RMSE     = "RMSE \u2193",
      MAE      = "MAE \u2193",
      R2       = "R\u00B2 \u2191",
      Coverage = "Coverage \u2191",
      Width    = "Interval Width \u2193"
    ),
    LabelText = dplyr::case_when(
      Metric %in% c("Coverage","R2") ~ ifelse(is.na(Value), "NA", scales::percent(Value, accuracy = 0.1)),
      TRUE                            ~ ifelse(is.na(Value), "NA", scales::number(Value, accuracy = 0.01))
    )
  )

model_colors <- c(
  "Linear Model"  = "#4E79A7",
  "LASSO"         = "#F28E2B",
  "Random Forest" = "#59A14F"
)
theme_vis <- theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold"),
    strip.background   = element_rect(fill = "#F3F4F5", color = NA),
    axis.text.x        = element_text(face = "bold"),
    plot.title.position= "plot",
    plot.margin        = margin(8, 8, 8, 8),
    plot.caption       = element_text(size = 9, color = "#5f6368")
  )
p_comp <- ggplot(metrics_long, aes(x = Model, y = Value, fill = Model)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = LabelText), vjust = -0.35, size = 3, fontface = "bold") +
  facet_grid(MetricLabel ~ Split, scales = "free_y") +
  scale_fill_manual(values = model_colors, breaks = levels(metrics_long$Model)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    title    = "Model Performance Comparison",
    subtitle = "Training-LOSO vs Training-InSample vs Holdout-Test (with CV-Conformal)",
    caption  = "Arrows indicate preferred direction (↑ higher is better, ↓ lower is better).",
    x = NULL, y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_vis
print(p_comp)

# 15. VISUALIZATION — Heatmap summary
safe_rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) rep(0.5, length(x))
  else scales::rescale(x, to = c(0,1), from = rng)
}
metrics_heat <- metrics_long %>%
  dplyr::group_by(Metric, Split) %>%
  dplyr::mutate(
    Score01 = dplyr::case_when(
      Metric %in% c("R2","Coverage") ~ safe_rescale(Value),
      TRUE                            ~ safe_rescale(-Value)
    )
  ) %>%
  dplyr::ungroup()
p_heat <- ggplot(metrics_heat, aes(x = Model, y = MetricLabel, fill = Score01)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = LabelText), size = 3, fontface = "bold") +
  facet_grid(~ Split) +
  scale_fill_viridis_c(name = "Normalized score\n(higher is better)", limits = c(0,1), na.value = "grey90") +
  labs(
    title    = "Model Performance — Heatmap Summary",
    subtitle = "0–1 normalized within each Metric × Split",
    caption  = "Labels show original units (%, RMSE, etc.).",
    x = NULL, y = NULL
  ) +
  theme_vis + theme(panel.grid = element_blank())
print(p_heat)

# 16. IMPUTE MISSING PM2.5 WITH THE BEST MODEL (RF, quantile forest)
rows_na <- which(is.na(pm25$PM25))
n_na    <- length(rows_na)
n_total <- nrow(pm25)
message(sprintf("PM25 missing: %d / %d rows", n_na, n_total))

if (n_na == 0) {
  message("No missing PM25 found. Skipping imputation & export.")
} else {
  # Freeze the scaler using TRAINING (observed) data statistics
  stopifnot(exists("pm25_clean"), exists("cont_vars"))
  scaler_mean <- sapply(pm25_clean[cont_vars], function(x) mean(x, na.rm = TRUE))
  scaler_sd   <- sapply(pm25_clean[cont_vars], function(x) sd(x,   na.rm = TRUE))
  scaler_sd[!is.finite(scaler_sd) | scaler_sd == 0] <- 1
  
  scale_with_training <- function(df) {
    for (nm in names(scaler_mean)) {
      if (nm %in% names(df)) {
        df[[nm]] <- (as.numeric(df[[nm]]) - scaler_mean[[nm]]) / scaler_sd[[nm]]
      }
    }
    df
  }
  
  # Build the subset to impute, applying the SAME pre-processing pipeline
  stopifnot(exists("month_lv"), exists("type_lv"), exists("ur_lv"))
  
  data_to_impute <- pm25 %>%
    dplyr::slice(rows_na) %>%
    dplyr::mutate(
      # enforce the same Type grouping: keep known top-3 levels, others -> "Other"
      Type  = as.character(Type),
      Type  = ifelse(Type %in% setdiff(type_lv, "Other"), Type, "Other"),
      Type  = factor(Type,  levels = type_lv),
      
      # harmonize factor levels
      Month = factor(as.character(Month), levels = month_lv),
      UR    = factor(as.character(UR),    levels = ur_lv),
      
      # simple type hygiene
      Year  = as.integer(Year)
    ) %>%
    scale_with_training()
  
  # Build design matrix from RHS-only terms and align to RF training matrix
  stopifnot(exists("rf_formula"), exists("X_train_rf"), exists("align_design"))
  
  # Use RHS-only terms so model.matrix won't look for 'log_PM25'
  mm_rf_rhs <- delete.response(terms(rf_formula))
  
  X_imp <- model.matrix(mm_rf_rhs, data = data_to_impute, na.action = na.pass)
  
  # Drop intercept column if present (to match how X_train_rf was built with [,-1])
  if ("(Intercept)" %in% colnames(X_imp)) {
    X_imp <- X_imp[, setdiff(colnames(X_imp), "(Intercept)"), drop = FALSE]
  }
  
  X_impA <- align_design(X_train_rf, X_imp)
  df_imp <- data.frame(X_impA)
  
  # Predict with the FULL-TRAINING quantile forest
  stopifnot(exists("rf_fit_q"))
  q_pred_imp <- predict(
    rf_fit_q, data = df_imp, type = "quantiles", quantiles = c(0.025, 0.5, 0.975)
  )$predictions
  if (is.null(dim(q_pred_imp))) q_pred_imp <- matrix(q_pred_imp, ncol = 3)
  colnames(q_pred_imp) <- c("log_q025", "log_q50", "log_q975")
  
  imputed_tbl <- tibble::tibble(
    .row_id    = rows_na,
    PM25_pred  = exp(q_pred_imp[, "log_q50"])  - 1,
    PM25_lwr   = exp(q_pred_imp[, "log_q025"]) - 1,
    PM25_upr   = exp(q_pred_imp[, "log_q975"]) - 1
  )
  
  # Merge predictions back to the original data frame
  pm25_filled <- pm25 %>%
    dplyr::mutate(.row_id = dplyr::row_number()) %>%
    dplyr::left_join(imputed_tbl, by = ".row_id") %>%
    dplyr::mutate(
      PM25_filled = dplyr::if_else(is.na(PM25), PM25_pred, PM25),
      ImputedFlag = is.na(PM25)
    ) %>%
    dplyr::select(-.row_id)
  
  # Export a single CSV (UTF-8)
  out_csv <- "Data - PM25_imputed.csv"
  readr::write_csv(pm25_filled, out_csv)
  
  # Audit summary
  message(sprintf(
    "Imputed %d PM25 values using RF quantile median. File written: %s",
    n_na, out_csv
  ))
  message("Quick check (first 5 imputed rows):")
  print(
    pm25_filled %>%
      dplyr::filter(ImputedFlag) %>%
      dplyr::select(PM25_filled, PM25_pred, PM25_lwr, PM25_upr) %>%
      dplyr::slice_head(n = 5)
  )
}

p_dist <- bind_rows(
  pm25 %>% filter(!is.na(PM25)) %>% transmute(flag = "Observed", value = PM25),
  pm25_filled %>% filter(ImputedFlag) %>% transmute(flag = "Imputed", value = PM25_filled)
) %>%
  ggplot(aes(value, fill = flag)) +
  geom_histogram(alpha = .5, position = "identity", bins = 30) +
  labs(title = "PM2.5: Observed vs Imputed", x = expression(PM[2.5]), y = "Count") +
  theme_minimal()

print(p_dist)

p_ts <- pm25_filled %>%
  mutate(
    Month = as.character(Month),
    Month = ifelse(Month %in% month.abb, Month, month.abb[pmax(1, pmin(12, as.integer(Month)))]),
    Month = factor(Month, levels = month.abb)
  ) %>%
  group_by(Month, ImputedFlag) %>%
  summarise(m = mean(if_else(ImputedFlag, PM25_filled, PM25), na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(Month, m, group = ImputedFlag, linetype = ImputedFlag)) +
  geom_line() + geom_point() +
  labs(title = "Monthly mean PM2.5: observed vs imputed", y = expression(bar(PM[2.5]))) +
  theme_minimal()

print(p_ts)
