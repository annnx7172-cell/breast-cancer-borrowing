library(survival)
data(cancer, package = "survival")

dim(gbsg)    # should print 686 11
head(gbsg)   # first 6 women
str(gbsg)    # type of each column
?gbsg        # help page
 
library(tidyverse)

df <- gbsg %>%
  mutate(
    treatment = factor(hormon, levels = c(0, 1),
                       labels = c("No hormone therapy", "Hormone therapy")),
    menopause = factor(meno, levels = c(0, 1),
                       labels = c("Pre-menopausal", "Post-menopausal")),
    grade     = factor(grade),
    rfs_years = rfstime / 365.25     # days -> years, easier to read
  )

table(df$treatment)          # how many in each group?
table(df$status)             # how many events vs censored?
summary(df$rfs_years)        # follow-up time in years


library(gtsummary)

tbl1 <- df %>%
  select(treatment, age, menopause, size, grade, nodes, pgr, er) %>%
  tbl_summary(
    by = treatment,
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",
    label = list(age ~ "Age (years)", menopause ~ "Menopause status",
                 size ~ "Tumour size (mm)", grade ~ "Tumour grade",
                 nodes ~ "Positive lymph nodes",
                 pgr ~ "Progesterone receptor (fmol/l)",
                 er ~ "Oestrogen receptor (fmol/l)")
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_labels()

tbl1

# Save Table 1
tbl1 %>% as_gt() %>% gt::gtsave("output/table1.html")

# Kaplan-Meier curves
library(ggsurvfit)

km <- survfit2(Surv(rfs_years, status) ~ treatment, data = df)
km                                   # median recurrence-free time per group

ggsurvfit(km) +
  add_confidence_interval() +
  add_risktable() +
  add_pvalue(caption = "Log-rank {p.value}") +
  labs(x = "Years since surgery",
       y = "Recurrence-free survival probability")

ggsave("output/km_treatment.png", width = 8, height = 6)

# Log-rank test and 5-year survival
survdiff(Surv(rfs_years, status) ~ treatment, data = df)
summary(km, times = 5)

# 1. True median follow-up (reverse Kaplan-Meier: "censoring" becomes the event)
survfit(Surv(rfs_years, 1 - status) ~ 1, data = df)

# 2. ANOVA: does tumour size differ across grades?
summary(aov(size ~ grade, data = df))

# 3. Curves by number of positive nodes
df <- df %>%
  mutate(node_group = cut(nodes, breaks = c(0, 3, 9, Inf),
                          labels = c("1-3", "4-9", "10+")))

km_nodes <- survfit2(Surv(rfs_years, status) ~ node_group, data = df)
ggsurvfit(km_nodes) +
  add_risktable() +
  add_pvalue(caption = "Log-rank {p.value}") +
  labs(x = "Years since surgery", y = "Recurrence-free survival probability")
ggsave("output/km_nodes.png", width = 8, height = 6)
km_nodes