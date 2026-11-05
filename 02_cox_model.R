library(survival)
library(tidyverse)
library(gtsummary)

# ---- Data preparation (same as script 01) ----
data(cancer, package = "survival")
df <- gbsg %>%
  mutate(
    treatment = factor(hormon, levels = c(0, 1),
                       labels = c("No hormone therapy", "Hormone therapy")),
    menopause = factor(meno, levels = c(0, 1),
                       labels = c("Pre-menopausal", "Post-menopausal")),
    grade     = factor(grade),
    rfs_years = rfstime / 365.25
  )

# Make receptor values easier to read (per 100 units instead of per 1)
df$pgr100 <- df$pgr / 100
df$er100  <- df$er / 100

# Model 1: compare treatment groups only
cox1 <- coxph(Surv(rfs_years, status) ~ treatment, data = df)

# Model 2: treatment plus other patient factors
cox2 <- coxph(Surv(rfs_years, status) ~ treatment + age + menopause +
                size + grade + nodes + pgr100 + er100, data = df)

# Hazard ratios with 95% confidence intervals
exp(cbind(HR = coef(cox1), confint(cox1)))
exp(cbind(HR = coef(cox2), confint(cox2)))

# Nice table for the report
tbl_cox <- tbl_regression(cox2, exponentiate = TRUE)
tbl_cox

# Test whether each hazard ratio stays constant over time
ph_test <- cox.zph(cox2)
ph_test

# Plot for treatment only (first variable in the model)
plot(ph_test[1])

# Red line = the model's average treatment effect (log hazard ratio)
abline(h = coef(cox2)[1], col = "red")


# Stratify by grade and menopause (they broke the assumption)
# strata() = "let each group have its own baseline risk"
cox3 <- coxph(Surv(rfs_years, status) ~ treatment + age + size + nodes +
                pgr100 + er100 + strata(grade, menopause), data = df)

# Hazard ratios with 95% CI
exp(cbind(HR = coef(cox3), confint(cox3)))

# Check the assumption again
cox.zph(cox3)