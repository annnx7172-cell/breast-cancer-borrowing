library(survival)                               # for the gbsg data

old_ctrl <- addm %>%                            # start from the TAX316 data
  transmute(                                    # keep only the new columns we create
    age       = AGE,                            # age in years
    menopause = ifelse(PRE_POST %in% c(1, 3),   # codes 1 and 3 count as pre-menopausal...
                       "Pre", "Post"),          # ...codes 2 and 4 as post-menopausal
    nodes     = NBLYNOI,                        # number of positive lymph nodes
    meno_code = PRE_POST,                       # keep the original code for the sensitivity check
    source    = "TAX316 (old)"                  # label where these patients came from
  )

cur_ctrl <- gbsg %>%                            # start from the main study
  filter(hormon == 0) %>%                       # keep only the untreated women (the current controls)
  transmute(
    age       = age,                            # age in years
    menopause = ifelse(meno == 0, "Pre", "Post"),   # gbsg codes 0 = pre, 1 = post
    nodes     = nodes,                          # number of positive lymph nodes
    meno_code = NA,                             # no original code needed here
    source    = "GBSG (current)"                # label the source
  )

controls <- bind_rows(cur_ctrl, old_ctrl)       # stack both groups into one table

controls %>% count(source, menopause)           # check: 440 GBSG + 746 TAX316, split by menopause


library(gtsummary)                              # makes publication-style summary tables

tbl_controls <- controls %>%                    # start from the combined controls table
  select(source, age, menopause, nodes) %>%     # keep only the variables to compare
  tbl_summary(by = source) %>%                  # one column per source, with median (IQR) or n (%)
  add_p()                                       # add a p-value testing whether the two sources differ

tbl_controls                                    # show the table in the Viewer

tbl_controls %>%
  as_gt() %>%                                   # convert to a gt table so it can be saved
  gt::gtsave("output/table3_controls.html")     # save it next to your other tables


library(DescTools)                                  # has HotellingsT2Test()

controls <- controls %>%
  mutate(post = ifelse(menopause == "Post", 1, 0))  # turn menopause into a number (1 = post, 0 = pre) so it can go in the test

hot <- HotellingsT2Test(cbind(age, nodes, post) ~ source,   # compare the average of age, nodes and post together between the two sources
                        data = controls)
hot                                                 # print the T² statistic, F value, degrees of freedom and p-value


cur_vars <- controls %>%
  filter(source == "GBSG (current)") %>%           # current controls only
  select(age, nodes, post)                          # the three shared variables

old_vars <- controls %>%
  filter(source == "TAX316 (old)") %>%              # old controls only
  select(age, nodes, post)

centre <- colMeans(cur_vars)                        # the "typical" current control: mean of each variable
covmat <- cov(cur_vars)                             # how the variables spread and relate in the current controls

d2 <- mahalanobis(old_vars, center = centre, cov = covmat)   # squared distance of each old patient from the typical current control

cutoff <- qchisq(0.975, df = 3)                     # cut-off: distances above the 97.5th percentile of chi-square with 3 df count as unusual (about 9.35)

sum(d2 > cutoff)                                    # how many old patients are unusual
mean(d2 > cutoff) * 100                             # the same as a percentage


old_flagged <- old_vars[d2 > cutoff, ]       # keep only the 20 unusual old patients
summary(old_flagged)                         # what makes them unusual (probably very high node counts)

d2_cur <- mahalanobis(cur_vars, center = centre, cov = covmat)   # distances for the current controls themselves
mean(d2_cur > cutoff) * 100                  # % of current controls flagged by the same cut-off


library(ggplot2)                                    # for plotting

pca <- prcomp(controls[, c("age", "nodes", "post")],   # run PCA on the three shared variables for all 1,186 controls
              scale. = TRUE)                        # standardise first so age (in years) doesn't dominate nodes and post

summary(pca)                                        # % of total variation captured by each principal component

pca_df <- data.frame(pca$x[, 1:2],                  # each patient's scores on PC1 and PC2
                     source = controls$source)      # attach where each patient came from

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = source)) +   # PC1 on x, PC2 on y, coloured by source
  geom_point(alpha = 0.4) +                         # semi-transparent points so overlapping areas show up
  stat_ellipse(level = 0.95) +                      # a 95% ellipse around each group's cloud
  labs(title = "Old vs current controls (PCA)", colour = "Source") +
  theme_minimal()

p_pca                                               # show the plot
ggsave("output/pca_controls.png", p_pca, width = 7, height = 5)   # save it with your other outputs

round(pca$rotation, 2)    # how much each variable contributes to each PC (large absolute values = main drivers)
