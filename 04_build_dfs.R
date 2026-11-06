library(haven)                                   # reads SAS files
library(dplyr)                                   # filter(), group_by(), summarise()
library(survival)                                # survival functions

adds <- read_sas("data/adds.sas7bdat")           # reload the disposition file (use data/tax316/... if that's your folder)
dim(adds)                                        # check it loaded: expect about 13,000 rows and 45 columns
dfs_events <- c("BREAST CANCER RELAPSE", "SECOND PRIMARY MALIGNANCY", "DEATH")

first_event <- adds %>%
  filter(DSDECOD %in% dfs_events, !is.na(DSSTDY)) %>%
  group_by(RUSUBJID) %>%
  summarise(event_day = min(DSSTDY))

nrow(first_event)

last_contact <- adds %>%
  group_by(RUSUBJID) %>%                                # one group per patient
  summarise(last_day = max(c(DSSTDY, RFENDY), na.rm = TRUE))   # latest day she was seen

dfs <- last_contact %>%
  left_join(first_event, by = "RUSUBJID") %>%           # add the event day (NA = no event)
  mutate(status    = ifelse(is.na(event_day), 0, 1),    # 1 = event, 0 = censored
         dfs_day   = ifelse(status == 1, event_day, last_day),   # event day, or last contact if censored
         rfs_years = dfs_day / 365.25)                  # days to years

nrow(first_event)                                       # patients with an event (expect about 300)
nrow(dfs)                                               # all patients (expect 746)
sum(dfs$status)                                         # number of events (should equal nrow(first_event))
summary(dfs$rfs_years)                                  # follow-up times (check there are no negative or zero values)
dfs %>% filter(rfs_years <= 0)        # show every patient with a zero or negative time


dfs <- dfs %>% filter(rfs_years > 0)          # drop the 1 patient with no follow-up (last contact before day 1)
nrow(dfs)                                     # should now be 745

old_ctrl <- addm %>%
  transmute(RUSUBJID,                          # patient ID, used to join the two tables
            age       = AGE,                   # age in years
            menopause = ifelse(PRE_POST %in% c(1, 3), "Pre", "Post"),   # same menopause mapping as before
            nodes     = NBLYNOI,               # number of positive nodes
            source    = "TAX316 (old)") %>%    # label the source
  inner_join(dfs %>% select(RUSUBJID, rfs_years, status),   # add DFS time and status; inner_join keeps only the 745 patients in both tables
             by = "RUSUBJID")

nrow(old_ctrl)                                # expect 745

km_old <- survfit(Surv(rfs_years, status) ~ 1, data = old_ctrl)   # Kaplan–Meier curve for the old controls
summary(km_old, times = c(5, 10))             # DFS at 5 and 10 years, with 95% CIs



library(ggsurvfit)                                   # for survival plots

cur_ctrl_surv <- gbsg %>%
  filter(hormon == 0) %>%                            # current controls: untreated GBSG women
  transmute(age, 
            menopause = ifelse(meno == 0, "Pre", "Post"),
            nodes,
            source    = "GBSG (current)",
            rfs_years = rfstime / 365.25,            # days to years
            status)                                  # 1 = recurrence or death, 0 = censored

all_ctrl <- bind_rows(cur_ctrl_surv,
                      old_ctrl %>% select(-RUSUBJID))   # stack current and old controls (drop the ID, which GBSG doesn't have)

km_ctrl <- survfit2(Surv(rfs_years, status) ~ source, data = all_ctrl)   # one KM curve per source

p_km <- km_ctrl %>%
  ggsurvfit() +                                      # draw the curves
  add_confidence_interval() +                        # shaded 95% CI bands
  add_risktable() +                                  # numbers at risk below the plot
  labs(x = "Years", y = "Recurrence/disease-free survival",
       title = "Old vs current control groups")

p_km
ggsave("output/km_controls.png", p_km, width = 7, height = 5)   # save the plot

survdiff(Surv(rfs_years, status) ~ source, data = all_ctrl)      # log-rank test: do the two control curves differ?
p_km7 <- p_km + coord_cartesian(xlim = c(0, 7))     # show only the first 7 years, where both groups have follow-up
p_km7
ggsave("output/km_controls.png", p_km7, width = 7, height = 6)   # overwrite the saved plot