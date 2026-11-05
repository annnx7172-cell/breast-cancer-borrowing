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

# Keep women whose 5-year outcome is known:
#   had an event within 5 years, OR were followed for more than 5 years
df5 <- subset(df, (status == 1 & rfs_years <= 5) | rfs_years > 5)

# Outcome: 1 = recurrence/death within 5 years, 0 = no
df5$recur5 <- ifelse(df5$status == 1 & df5$rfs_years <= 5, 1, 0)
table(df5$recur5)

# Logistic regression
logit <- glm(recur5 ~ treatment + age + menopause + size + grade +
               nodes + pgr100 + er100, data = df5, family = binomial)

# Odds ratios with 95% CI
exp(cbind(OR = coef(logit), confint(logit)))

# Multicollinearity check
library(car)
vif(logit)


# ---- Subgroup 1: menopause ----
# Treatment HR in pre-menopausal women only
cox_pre <- coxph(Surv(rfs_years, status) ~ treatment,
                 data = subset(df, menopause == "Pre-menopausal"))
exp(cbind(HR = coef(cox_pre), confint(cox_pre)))

# Treatment HR in post-menopausal women only
cox_post <- coxph(Surv(rfs_years, status) ~ treatment,
                  data = subset(df, menopause == "Post-menopausal"))
exp(cbind(HR = coef(cox_post), confint(cox_post)))

# Interaction test: does the treatment effect differ by menopause?
cox_int1 <- coxph(Surv(rfs_years, status) ~ treatment * menopause, data = df)
summary(cox_int1)

# ---- Subgroup 2: lymph nodes ----
# Make two groups: 1-3 nodes vs 4 or more
df$nodes_grp <- ifelse(df$nodes <= 3, "1-3 nodes", "4+ nodes")

cox_low <- coxph(Surv(rfs_years, status) ~ treatment,
                 data = subset(df, nodes_grp == "1-3 nodes"))
exp(cbind(HR = coef(cox_low), confint(cox_low)))

cox_high <- coxph(Surv(rfs_years, status) ~ treatment,
                  data = subset(df, nodes_grp == "4+ nodes"))
exp(cbind(HR = coef(cox_high), confint(cox_high)))

# Interaction test: does the treatment effect differ by nodes?
cox_int2 <- coxph(Surv(rfs_years, status) ~ treatment * nodes_grp, data = df)
summary(cox_int2)


# Small helper: pulls the treatment HR and 95% CI out of a model
get_hr <- function(model) {
  exp(c(coef(model)[1], confint(model)[1, ]))
}

# Collect all results into one table
forest <- rbind(
  get_hr(cox1),      # all patients
  get_hr(cox_pre),   # pre-menopausal
  get_hr(cox_post),  # post-menopausal
  get_hr(cox_low),   # 1-3 nodes
  get_hr(cox_high)   # 4+ nodes
)
forest <- as.data.frame(forest)
names(forest) <- c("hr", "lower", "upper")
forest$group <- c("All patients", "Pre-menopausal", "Post-menopausal",
                  "1-3 nodes", "4+ nodes")

# Keep the groups in this order from top to bottom
forest$group <- factor(forest$group, levels = rev(forest$group))
forest

# Draw the forest plot
ggplot(forest, aes(x = hr, y = group)) +
  geom_point(size = 3) +                                    # the HR
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +  # the 95% CI
  geom_vline(xintercept = 1, linetype = "dashed") +         # line of no effect
  scale_x_log10() +                                         # ratio scale
  labs(x = "Hazard ratio (hormone therapy vs none, 95% CI)", y = NULL,
       title = "Effect of hormone therapy across subgroups") +
  theme_minimal()

ggsave("output/forest_plot.png", width = 8, height = 5)