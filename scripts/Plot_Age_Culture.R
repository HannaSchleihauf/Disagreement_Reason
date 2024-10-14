
############################################################################
# PLOT AGE CULTURE
############################################################################
# Model with interaction
age.culture.model <-
  glmer(resp_mat ~ condition.disagreement.code + (culture +
                                                    z.age)^2 + test.trial.2.code +
          (1 + test.trial.2.code || dyad.nr),
        data = t.data, control = contr, family = binomial(link = "logit")
  )
boot_res_age_culture <-
  boot.glmm.pred(
    model.res = age.culture.model, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("culture", "z.age")
  )
# Model without interaction
age.culture.model.2 <-
  glmer(resp_mat ~ condition.disagreement.code + (culture +
                                                    z.age) + test.trial.2.code +
          (1 + test.trial.2.code || dyad.nr),
        data = t.data, control = contr, family = binomial(link = "logit")
  )
boot_res_age_culture.2 <-
  boot.glmm.pred(
    model.res = age.culture.model.2, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("culture", "z.age")
  )

# Aggregate Data
ydata$z.age <- scale(ydata$exact.age.children.dyad)
ydata.agg.age.culture <- ydata %>%
  group_by(dyad.nr, culture, z.age, condition) %>%
  summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
  ungroup()
# breaks for x axis
breaks <- c(
  (4 - mean(ydata$exact.age.children.dyad)) /
    sd(ydata$exact.age.children.dyad),
  (5 - mean(ydata$exact.age.children.dyad)) /
    sd(ydata$exact.age.children.dyad),
  (6 - mean(ydata$exact.age.children.dyad)) /
    sd(ydata$exact.age.children.dyad),
  (7 - mean(ydata$exact.age.children.dyad)) /
    sd(ydata$exact.age.children.dyad),
  (8 - mean(ydata$exact.age.children.dyad)) /
    sd(ydata$exact.age.children.dyad),
  (9 - mean(ydata$exact.age.children.dyad)) /
    sd(ydata$exact.age.children.dyad),
  (10 - mean(ydata$exact.age.children.dyad)) /
    sd(ydata$exact.age.children.dyad)
)# Plot Age Culture
exp1_plot_age_culture.2 <-
  ggplot() +
  geom_point(
    data = ydata.agg.age.culture,
    aes(x = z.age, y = mean.resp),
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.00, h = 0.015)
  ) +
  geom_ribbon(
    data = boot_res_age_culture.2$ci.predicted,
    aes(ymin = lower.cl, ymax = upper.cl, x = z.age, group = culture, fill = culture),
    alpha = 0.1
  ) +
  geom_textline(
    data = boot_res_age_culture.2$ci.predicted,
    aes(x = z.age, y = fitted, group = culture, color = culture, label = culture),
    hjust = c(0.3)
  ) +
  scale_x_continuous(
    name = "Age",
    breaks = breaks,
    labels = c(4, 5, 6, 7, 8, 9, 10)
  ) +
  scale_y_continuous(
    name = "Proportion of reasons given",
    labels = scales::percent
  ) +
  # labs(title = "Interaction effect of age and culture") +
  scale_color_npg() +
  scale_fill_npg() +
  theme_classic() +
  my_theme +
  theme(legend.position = "none")
exp1_plot_age_culture.2
