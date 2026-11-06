library(survival)                                   # survival models
library(dplyr)                                      # data handling

trial_full <- gbsg %>%
  transmute(treat     = hormon,                      # 1 = hormone therapy, 0 = none
            age,                                     # age in years
            menopause = ifelse(meno == 0, "Pre", "Post"),   # same coding as the old controls
            nodes,                                   # number of positive nodes
            rfs_years = rfstime / 365.25,            # days to years
            status)                                  # 1 = recurrence or death

hr_table <- function(fit) {                          # small helper: pulls out the treatment HR, 95% CI and p-value
  ci <- summary(fit)$conf.int["treat", c(1, 3, 4)]   # HR, lower and upper CI
  p  <- summary(fit)$coefficients["treat", "Pr(>|z|)"]   # p-value
  round(c(HR = ci[1], lower = ci[2], upper = ci[3], p = p), 3)
}

# 1. Benchmark: the full trial (the "true answer")
fit_full <- coxph(Surv(rfs_years, status) ~ treat + age + menopause + nodes,
                  data = trial_full)                 # Cox model adjusted for the 3 shared variables
hr_table(fit_full)

# 2. Pretend small trial: all treated women + a random third of the controls
set.seed(2026)                                       # fix the random sample so you get the same result each time
small_ctrl  <- trial_full %>% filter(treat == 0) %>% slice_sample(prop = 1/3)   # keep a random 1/3 of the 440 controls
small_trial <- bind_rows(trial_full %>% filter(treat == 1), small_ctrl)          # add them to all 246 treated women

table(small_trial$treat)                             # check: 246 treated, about 147 controls

fit_small <- coxph(Surv(rfs_years, status) ~ treat + age + menopause + nodes,
                   data = small_trial)               # same model on the small trial
hr_table(fit_small)


max_fu <- max(trial_full$rfs_years)                  # longest follow-up in GBSG (about 7.3 years)

old_borrow <- old_ctrl %>%
  transmute(treat     = 0,                            # old patients are all controls
            age, menopause, nodes,                    # the 3 shared variables
            status    = ifelse(rfs_years > max_fu, 0, status),   # an event after the GBSG window becomes censored
            rfs_years = pmin(rfs_years, max_fu))      # cut follow-up at the GBSG maximum

weights <- c(0, 0.25, 0.5, 1)                         # how much each old patient counts
results <- data.frame()                               # empty table to collect the results

for (w in weights) {                                  # repeat the analysis for each weight
  if (w == 0) {
    dat <- small_trial %>% mutate(wt = 1)             # weight 0 = no borrowing (coxph can't take zero weights, so leave the old patients out)
  } else {
    dat <- bind_rows(small_trial %>% mutate(wt = 1),  # current patients count fully
                     old_borrow  %>% mutate(wt = w))  # old controls count as a fraction w
  }
  fit <- coxph(Surv(rfs_years, status) ~ treat + age + menopause + nodes,
               data = dat, weights = wt, robust = TRUE)   # weighted Cox model with robust standard errors
  results <- rbind(results, c(weight = w, hr_table(fit)))  # add this weight's HR, CI and p to the table
}

names(results) <- c("weight", "HR", "lower", "upper", "p")   # tidy column names
results    
# compare each row with the benchmark HR 0.69



ttp_data <- bind_rows(small_ctrl %>% mutate(source = "current"),   # the 146 current controls in the small trial
                      old_borrow %>% mutate(source = "old"))       # the 745 old controls (follow-up cut at the GBSG maximum)

lr   <- survdiff(Surv(rfs_years, status) ~ source, data = ttp_data)   # log-rank test: do current and old controls differ?
p_lr <- 1 - pchisq(lr$chisq, df = 1)                  # turn the chi-square statistic into a p-value
p_lr                                                  # small p = the outcomes differ

w_ttp <- ifelse(p_lr > 0.05, 1, 0)                    # rule: borrow 100% if no evidence of a difference, otherwise borrow 0%
w_ttp# 0 means test-then-pool refuses to borrow

fit_hr <- function(dat) {                             # fit the weighted Cox model and return the treatment HR and p-value
  fit <- coxph(Surv(rfs_years, status) ~ treat + age + menopause + nodes,
               data = dat, weights = wt, robust = TRUE)
  s <- summary(fit)
  c(hr = unname(s$conf.int["treat", 1]),              # hazard ratio for treatment
    p  = unname(s$coefficients["treat", "Pr(>|z|)"])) # its p-value
}

analyse <- function(trt, cur, hist) {                 # run all 4 methods on one simulated trial
  trial <- bind_rows(trt, cur) %>% mutate(wt = 1)     # the small trial: everyone counts fully
  none  <- fit_hr(trial)                              # method 1: no borrowing
  half  <- fit_hr(bind_rows(trial, hist %>% mutate(wt = 0.5)))   # method 2: old controls count as half
  full  <- fit_hr(bind_rows(trial, hist %>% mutate(wt = 1)))     # method 3: old controls count fully
  lr    <- survdiff(Surv(rfs_years, status) ~ grp,
                    data = bind_rows(cur %>% mutate(grp = "cur"),
                                     hist %>% mutate(grp = "hist")))   # log-rank test: current vs old controls
  ttp   <- if (1 - pchisq(lr$chisq, 1) > 0.05) full else none        # method 4: test-then-pool
  data.frame(method = c("none", "half", "full", "ttp"),
             hr = c(none["hr"], half["hr"], full["hr"], ttp["hr"]),
             p  = c(none["p"],  half["p"],  full["p"],  ttp["p"]))
}

ctrl_all    <- trial_full %>% filter(treat == 0)      # all 440 GBSG controls
treated_all <- trial_full %>% filter(treat == 1)      # all 246 treated women


n_sim <- 200                                           # start with 20 to test the code quickly, then change to 500
set.seed(123)                                         # make the simulation repeatable
sim <- list()                                         # empty list to collect each repeat's results

for (i in 1:n_sim) {
  idx  <- sample(nrow(ctrl_all))                      # shuffle the 440 controls
  cur  <- ctrl_all[idx[1:146], ]                      # 146 current controls for the small trial
  fake <- ctrl_all[idx[147:292], ] %>% mutate(treat = 1)   # 146 controls pretending to be treated (true HR = 1)
  exch <- ctrl_all[idx[293:440], ]                    # 148 leftover controls = truly similar "historical" data
  
  sim[[i]] <- bind_rows(
    analyse(treated_all, cur, old_borrow) %>% mutate(setting = "benefit", hist = "conflicting"),
    analyse(treated_all, cur, exch)       %>% mutate(setting = "benefit", hist = "exchangeable"),
    analyse(fake, cur, old_borrow)        %>% mutate(setting = "null",    hist = "conflicting"),
    analyse(fake, cur, exch)              %>% mutate(setting = "null",    hist = "exchangeable")
  ) %>% mutate(sim = i)                               # label which repeat this was
}

sim_df <- bind_rows(sim)                              # stack all repeats into one table

sim_summary <- sim_df %>%
  group_by(setting, hist, method) %>%                 # one row per setting, data source and method
  summarise(median_HR   = round(median(hr), 2),       # typical HR across repeats
            claim_benefit = mean(p < 0.05 & hr < 1),  # % of trials concluding "therapy works" (power, or type I error in the null setting)
            claim_harm    = mean(p < 0.05 & hr > 1),  # % of trials concluding "therapy harms"
            .groups = "drop")
sim_summary



library(ggplot2)

library(tidyr)                                       # for pivot_longer()

sim_long <- sim_summary %>%
  pivot_longer(c(claim_benefit, claim_harm),
               names_to = "conclusion", values_to = "prop") %>%
  mutate(conclusion = dplyr::recode(conclusion,             # dplyr:: tells R to use dplyr's recode, not car's
                                    claim_benefit = "Claims benefit",
                                    claim_harm    = "Claims harm"),
         method = factor(method,
                         levels = c("none", "half", "full", "ttp"),
                         labels = c("None", "50%", "100%", "Test-then-pool")),
         hist = dplyr::recode(hist, conflicting  = "Conflicting (TAX316)",
                              exchangeable = "Similar controls"),
         setting = dplyr::recode(setting, benefit = "True benefit (HR 0.69)",
                                 null    = "No true effect (HR 1)"))

p_sim <- ggplot(sim_long, aes(x = method, y = prop, fill = conclusion)) +
  geom_col() +                                       # stacked bars: green = benefit, red = harm
  facet_grid(setting ~ hist) +                       # 2×2 grid
  geom_hline(yintercept = 0.05, linetype = "dashed") +   # 5% reference line
  scale_fill_manual(values = c("Claims benefit" = "#1b9e77",
                               "Claims harm"    = "#d95f02")) +
  scale_y_continuous(labels = scales::percent) +     # show the y-axis as percentages
  labs(x = NULL, y = "% of 200 simulated trials", fill = NULL,
       title = "What each borrowing method concludes") +
  theme_minimal() + theme(legend.position = "bottom")

p_sim
ggsave("output/simulation_results.png", p_sim, width = 9, height = 6)


sim_long <- sim_long %>%
  mutate(setting = factor(setting, levels = c("True benefit (HR 0.69)",
                                              "No true effect (HR 1)")),   # put the benefit row on top
         method  = factor(method, labels = c("None", "50%", "100%", "Test-\nthen-pool")))   # \n splits the label over two lines

p_sim <- p_sim %+% sim_long                                  # redraw the same plot with the updated data
p_sim
ggsave("output/simulation_results.png", p_sim, width = 9, height = 6)
