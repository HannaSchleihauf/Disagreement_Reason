# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 1: We expect that children from all three cultural
# settings are more likely to give reasons in the disagreement relative to the
# agreement condition (i.e., we expect no cultural differences in Study 1).

# September 2023

############################################################################
# INSTALL OR LOAD PACKAGES & SOURCE FUNCTIONS
############################################################################
# Vector of required packages
required_packages <- c(
  "xlsx", "tidyverse", "emmeans", "effects", "ggdist", "dplyr",
  "lme4", "directlabels", "ggsci", "geomtextpath", "exactRankTests"
)
# Load or install packages
lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})
# Function files to source
function_files <- c(
  "./functions/diagnostic_fcns.r", "./functions/glmm_stability.r",
  "./functions/drop1_para.r", "./functions/boot_glmm.r",
  "./functions/paired_data_functions.r"
)
# Source function files
lapply(function_files, source)

############################################################################
# R SETTINGS
############################################################################
options(scipen = 999)

############################################################################
# DATA
############################################################################
xdata <-
  read.csv(
    file =
      "./data/Disagreement_Reason_Coding_Sheet_All_Three.csv",
    head = TRUE, stringsAsFactors = TRUE
  )
str(xdata)

############################################################################
# DATA WRANGELING
############################################################################

xdata <- xdata %>% mutate(
  # Make. trial a factor
  test.trial = as.factor(test.trial),
  # Calculate the exact age
  exact.age.perceptual.child =
    (as.Date(test.date, "%m/%d/%y") -
      as.Date(birth.date.perceptual.child, "%m/%d/%y")) / 365.25,
  exact.age.testimonial.child =
    (as.Date(test.date, "%m/%d/%y") -
      as.Date(birth.date.testimonial.child, "%m/%d/%y")) / 365.25,
  exact.age.perceptual.child =
    ifelse(is.na(exact.age.perceptual.child),
      age.perceptual.child, exact.age.perceptual.child),
  exact.age.testimonial.child =
    ifelse(is.na(exact.age.testimonial.child),
      age.testimonial.child,exact.age.testimonial.child),
  exact.age.children.dyad =
    (exact.age.perceptual.child + exact.age.testimonial.child) / 2,
  # Create new variable with max reason giving number = 2
  perceptual.reason.given.dyad.nr =
      ifelse(perceptual.reason.given == "yes", 1, 0),
  testimonial.reason.given.dyad.nr =
      ifelse(testimonial.reason.given == "yes", 1, 0),
  nr.reason.dyad.max.2 =
      perceptual.reason.given.dyad.nr + testimonial.reason.given.dyad.nr,
  nr.reason.not.given.dyad.max.2 = 2 - nr.reason.dyad.max.2,
  # Create new numeric variable for the reliability coding
  perceptual.reason.given.dyad.nr.reli = case_when(
      perceptual.reason.given.reli == "yes" ~ 1,
      perceptual.reason.given.reli == "no" ~ 0,
      TRUE ~ NA_real_),
  testimonial.reason.given.dyad.nr.reli = case_when(
      testimonial.reason.given.reli == "yes" ~ 1,
      testimonial.reason.given.reli == "no" ~ 0,
      TRUE ~ NA_real_),
  nr.reason.dyad.max.2.reli =
      perceptual.reason.given.dyad.nr.reli +
        testimonial.reason.given.dyad.nr.reli,
  nr.reason.not.given.dyad.max.2.reli = 2 - nr.reason.dyad.max.2.reli)

############################################################################
# DO ANY TRIALS NEED TO BE EXCLUDED?
############################################################################
tapply(xdata$duration.inaudibles, xdata$culture, mean)
tapply(xdata$duration.inaudibles, xdata$culture, sd)

############################################################################
# INTERRATRER RELIABILITY
############################################################################
# extract only the rows for which we have a reliability coding
interrater <-
  xdata %>%
  filter(!is.na(nr.reason.dyad.max.2.reli))
# percentage of reliability coding 25%
nrow(interrater) / nrow(xdata)

# number of dyads for which the reliability coding was not the same
sum(interrater$nr.reason.dyad.max.2.reli !=
  interrater$nr.reason.dyad.max.2)

# Cohen's Kappa
library(irr)
kappa_results <-
  kappa2(cbind(
    interrater$nr.reason.dyad.max.2.reli,
    interrater$nr.reason.dyad.max.2
  ), "unweighted")
print(kappa_results)

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
# Use xx.fe.re function to determine the full random slope structure
xx.fe.re <- fe.re.tab(
  fe.model = "nr.reason.dyad.max.2 ~ condition*culture*exact.age.children.dyad +
              gender.dyad + test.trial",
  re = "(1|dyad.nr)",
  other.vars = "nr.reason.not.given.dyad.max.2",
  data = xdata
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

# Create a response matrix
resp_mat <-
  cbind(t.data$nr.reason.dyad.max.2, t.data$nr.reason.not.given.dyad.max.2)
t.data$olre <-
  as.factor(1:nrow(t.data))

# Center dummy variables
t.data$condition.disagreement.code <-
  t.data$condition.disagreement - mean(t.data$condition.disagreement)
t.data$culture.USA.code <-
  t.data$culture.USA - mean(t.data$culture.USA)
t.data$culture.Kenya.code <-
  t.data$culture.Kenya - mean(t.data$culture.Kenya)
t.data$gender.dyad.male.code <-
  t.data$gender.dyad.male - mean(t.data$gender.dyad.male)
t.data$test.trial.2.code <-
  t.data$test.trial.2 - mean(t.data$test.trial.2)
t.data$z.age <- scale(t.data$exact.age.children.dyad)

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full <- glmer(resp_mat ~
  condition * culture * z.age + test.trial + gender.dyad +
  (1 + test.trial.2.code || dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)
red <- glmer(resp_mat ~
  (condition + culture + z.age)^2 + test.trial + gender.dyad +
  (1 + test.trial.2.code || dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)
red2 <- glmer(resp_mat ~
  (condition + culture)^2 + z.age + test.trial + gender.dyad +
  (1 + test.trial.2.code || dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)
red3 <- glmer(resp_mat ~
  condition + (culture + z.age)^2 + test.trial + gender.dyad +
  (1 + test.trial.2.code || dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)
main <- glmer(resp_mat ~
  (condition + culture) + z.age + test.trial + gender.dyad +
  (1 + test.trial.2.code || dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)
null <- glmer(resp_mat ~ 1 +
  (1 + test.trial.2.code || dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)

summary(full)$varcor # the random slope can not be estimated
ranef.diagn.plot(full)
xx <- round(summary(full)$coefficients, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full)
library(car)
vif(main)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full, contr = contr, use = c("dyad.nr"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))
# Model stablility only with models that did not give a warning message
m.stab.b$detailed <- m.stab.b$detailed %>%
  filter(warnings == "none")
cbind(
  apply(m.stab.b$detailed, 2, min),
  apply(m.stab.b$detailed, 2, max)
)

############################################################################
# MODEL COMPARISONS
############################################################################
round(anova(full, null, test = "Chisq"), 3)
round(drop1(full, test = "Chisq"), 3)
round(drop1(red, test = "Chisq"), 3)
round(drop1(red2, test = "Chisq"), 3)
round(drop1(red3, test = "Chisq"), 3)
round(drop1(main, test = "Chisq"), 3)

# Quick look at the results
plot(effect("condition:culture", red2), type = "response")

# Coefficients of the full model
summary(full)$coefficients

emm1 <- emmeans(red2, ~ condition * culture)
summary(emm1, type = "response")
summary(pairs(regrid(emm1))[c(7, 9, 14, 1, 10, 15)],
  type = "response", adjust = "fdr"
)

############################################################################
# NON-PARAMETRIC WILCOX TEST
############################################################################
t.data.China <- t.data %>% filter(culture == "China")
t.data.Kenya <- t.data %>% filter(culture == "Kenya")
t.data.USA <- t.data %>% filter(culture == "USA")

wilcox.test(nr.reason.dyad.max.2 ~
  condition, data = t.data.China, paired = FALSE)
wilcox.test(nr.reason.dyad.max.2 ~
  condition, data = t.data.Kenya, paired = FALSE)
wilcox.test(nr.reason.dyad.max.2 ~
  condition, data = t.data.USA, paired = FALSE)


############################################################################
# PLOTTING PREPERATION
############################################################################
# Long data set for plotting
ydata <-
  xdata %>%
  pivot_longer(
    cols = c(
      perceptual.reason.given.dyad.nr,
      testimonial.reason.given.dyad.nr
    ),
    names_to = "type.reason",
    values_to = "reason.given"
  )

############################################################################
# BOOTSTRAP FULL MODEL
############################################################################
boot_full <-
  boot.glmm.pred(
    model.res = full, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("condition", "culture", "z.age", "test.trial", "gender.dyad")
  )
round(boot_full$ci.estimates, 3)
as.data.frame(round(boot_full$ci.estimates, 3))
m.stab.plot(round(boot_full$ci.estimates, 3))
boot_full$ci.predicted

############################################################################
# SETTING PLOT THEME
############################################################################
my_theme <-
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.background = element_blank(),
    legend.title = element_blank(),
    # legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13), # , family = "Arial Narrow"),
    strip.text.x = element_text(size = 15),
    # text = element_text(family = "Arial Narrow"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    )
  )

############################################################################
# PLOT CULTURE CONDITION
############################################################################
# Boot straps for plotting without interaction
plot.model.culture.condition <-
  glmer(resp_mat ~ (condition * culture) + z.age + test.trial.2.code +
    (1 + test.trial.2.code || dyad.nr),
  data = t.data, control = contr, family = binomial(link = "logit")
  )

boot_res_condition_culture.2 <-
  boot.glmm.pred(
    model.res = plot.model.culture.condition, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("condition", "culture")
  )

# Aggregated data for plotting
ydata.agg <- ydata %>%
  group_by(dyad.nr, condition, culture) %>%
  summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
  ungroup() %>%
  na.omit()

exp1_plot_condition.culture <-
  ggplot(
    data = ydata.agg,
    aes(x = factor(condition,
      levels = c("agreement", "disagreement")
    ), y = mean.resp)) +
  geom_point(
    data = ydata.agg,
    aes(x = condition, color = condition), size = 1.5,
    alpha = .4, position = position_jitter(w = 0.25, h = 0.025)) +
  geom_violin(
    data = ydata.agg,
    aes(x = condition, y = mean.resp, fill = condition),
    position = position_nudge(x = 0.00),
    alpha = .2, color = NA) +
  geom_errorbar(
    data = boot_res_condition_culture.2$ci.predicted,
    aes(
      x = as.numeric(condition), y = fitted,
      ymin = lower.cl, ymax = upper.cl, color = condition),
    width = 0.1, linewidth = 1) +
  geom_point(
    data = boot_res_condition_culture.2$ci.predicted,
    aes(x = as.numeric(condition), y = fitted, fill = condition),
    size = 3.5) +
  geom_point(
    data = boot_res_condition_culture.2$ci.predicted,
    aes(x = as.numeric(condition), y = fitted, fill = condition),
    size = 3.5) +
  facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_discrete(
    limits = c("agreement", "disagreement"),
    name = element_blank(),
    labels = c("Agreement\ncondition", "Disagreement\ncondition")) +
  scale_y_continuous(
    name = "Proportion of reasons given",
    labels = scales::percent) +
  # labs(title = "A) Reason-giving as a function of condition and culture") +
  scale_color_lancet() +
  scale_fill_lancet() +
  theme_light() +
  my_theme +
  theme(axis.text.x = element_text(size = 11),
        legend.position = "none",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 10, unit = "pt"))
exp1_plot_condition.culture

# Checking whethe rthe estimates are similar to the means
tapply(ydata.agg$mean.resp, list(ydata.agg$culture, ydata.agg$condition), FUN = mean)

############################################################################
# PLOT  AGE
############################################################################
# Bootstrap
age.model <-
  glmer(resp_mat ~ condition.disagreement.code +
    (culture.Kenya.code + culture.USA.code) +
    z.age + test.trial.2.code +
    (1 + test.trial.2.code || dyad.nr),
  data = t.data, control = contr, family = binomial(link = "logit")
  )

boot_res_age <-
  boot.glmm.pred(
    model.res = age.model, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("z.age")
  )

age.culture.model <-
  glmer(resp_mat ~ condition.disagreement.code +
    culture +
    z.age + test.trial.2.code +
    (1 + test.trial.2.code || dyad.nr),
  data = t.data, control = contr, family = binomial(link = "logit")
  )

boot_res_age_culture <-
  boot.glmm.pred(
    model.res = age.culture.model, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("z.age")
  )

# age.condition.model <-
#   glmer(resp_mat ~ (condition + (culture.Kenya.code + culture.USA.code) +
#     z.age) + test.trial.2.code +
#     (1 + test.trial.2.code || dyad.nr),
#   data = t.data, control = contr, family = binomial(link = "logit")
#   )
# boot_res_age_condition <-
#   boot.glmm.pred(
#     model.res = age.condition.model, excl.warnings = TRUE,
#     nboots = 1000, para = FALSE, level = 0.95,
#     use = c("condition", "z.age", "culture")
#   )

# Aggregate Data
ydata$z.age <- scale(ydata$exact.age.children.dyad)
ydata.agg.age <- ydata %>%
  group_by(condition, age.group) %>% # group_by(dyad.nr, condition, z.age) %>%
  summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
  ungroup()
ydata.agg.age$z.age <- as.numeric(scale(ydata.agg.age$age.group))

# Plot Age
set.seed(1)
exp1_plot_age <-
  ggplot() +
  geom_point(
    data = ydata.agg.age,
    aes(x = z.age, y = mean.resp, color = condition), # z.age
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.00, h = 0.015)
  ) +
  # geom_ribbon(
  #   data = boot_res_age_condition$ci.predicted,
  #   aes(ymin = lower.cl, ymax = upper.cl, x = z.age, fill = condition),
  #   alpha = 0.1
  # ) +
  # geom_line(
  #   data = boot_res_age_condition$ci.predicted,
  #   aes(x = z.age, y = fitted, color = condition)
  # ) +
  # geom_textline(
  #   data = boot_res_age_culture$ci.predicted,
  #   aes(x = z.age, y = fitted, group = culture,
  #       label = culture), #color = culture
  #   hjust = c(0.3), linewidth = 0.5, color = "gray30"
  # ) +
  # geom_ribbon(
  #   data = boot_res_age_culture$ci.predicted,
  #   aes(ymin = lower.cl, ymax = upper.cl,
  #       x = z.age, group = culture), #, fill = culture
  #   alpha = 0.2, fill = "gray50"
  # ) +
  geom_ribbon(
    data = boot_res_age$ci.predicted,
    aes(ymin = lower.cl, ymax = upper.cl, x = z.age),
    alpha = 0.5, fill = "grey"
  ) +
  geom_line(
    data = boot_res_age$ci.predicted, aes(x = z.age, y = fitted),
    color = "black", linewidth = 1
  ) +
  # facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_continuous(
    name = "Age",
    breaks = c(
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
    ),
    labels = c(5, 6, 7, 8, 9, 10)
  ) +
  scale_y_continuous(
    name = "Proportion of reasons given",
    labels = scales::percent
  ) +
  scale_color_lancet() +
  scale_fill_lancet() +
  theme_light() +
  my_theme +
  theme(plot.margin = margin(t = 25, r = 10, b = 10, l = 10, unit = "pt"))
exp1_plot_age

# exp1_plot_age +
#   scale_color_manual(values = c(
#     "gray50", "gray50",
#     "#0073C2FF", "#FFA319FF", "#868686FF"
#   ))

############################################################################
# PLOT TRIAL
############################################################################
ydata$z.trial <- scale(as.numeric(ydata$test.trial))
ydata.agg.trial <- ydata %>%
  group_by(dyad.nr, condition, test.trial) %>%
  summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
  ungroup()

trial.model.culture <-
  glmer(resp_mat ~ condition.disagreement.code +
    culture +
    z.age + test.trial +
    (1 + test.trial.2.code || dyad.nr),
  data = t.data, control = contr, family = binomial(link = "logit")
  )

trial.model <-
  glmer(resp_mat ~ (condition.disagreement.code +
    culture.Kenya.code + culture.USA.code) +
    z.age + test.trial +
    (1 + test.trial.2.code || dyad.nr),
  data = t.data, control = contr, family = binomial(link = "logit")
  )

boot_res_trial_culture <-
  boot.glmm.pred(
    model.res = trial.model, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("test.trial", "culture")
  )

boot_res_trial <-
  boot.glmm.pred(
    model.res = trial.model, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("test.trial")
  )

exp1_plot_trial <-
  ggplot() +
  geom_point(
    data = ydata.agg.trial,
    aes(x = test.trial, y = mean.resp, color = condition),
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.2, h = 0.05)
  ) +
  geom_violin(
    data = ydata.agg.trial,
    aes(x = test.trial, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "lightgrey", color = NA, alpha = .2, linewidth = 0
  ) +
  geom_errorbar(
    data = boot_res_trial$ci.predicted,
    aes(
      x = as.numeric(test.trial), y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "gray40",
    width = 0.1, size = 1
  ) +
  geom_point(
    data = boot_res_trial$ci.predicted,
    aes(x = as.numeric(test.trial), y = fitted),
    color = "black", size = 3.5
  ) +
  # facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_discrete(
    limits = c("1", "2"),
    name = "Test Trial",
    labels = c("1", "2")
  ) +
  scale_y_continuous(
    name = "Proportion of reasons given",
    labels = scales::percent
  ) +
  scale_color_lancet() +
  scale_fill_lancet() +
  # labs(title = "Reason-giving as a function of test trial") +
  theme_light() +
  my_theme +
  theme(
    legend.position = "none"
  ) +
  theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10, unit = "pt"))
exp1_plot_trial

############################################################################
# COMBINE PLOTS
############################################################################
library(ggpubr)
# theme_set(theme_pubr())
figure <- ggarrange(exp1_plot_age, exp1_plot_trial,
  labels = c(
    "(a)",
    "(b)"
  ),
  # label.x = 0,    # X position 0 for left
  # label.y = 1,    # Y position 1 for top
  #hjust = -1, # hjust = 0 for left alignment
  #vjust = c(0, 0), # hjust = 0 for left alignment
  ncol = 1, nrow = 2
)
figure
ggexport(figure, filename = "figure.png", width = 713, height = 993, res = 400)

save.image("RData_files/Hypothesis_1_&_2_&_3.RData")
load("RData_files/Hypothesis_1_&_2_&_3.RData")

ftable(box.chosen ~ dominant.child.1, xdata)
ftable(condition ~ box.chosen, xdata)
levels(xdata$dominant.child.1)
