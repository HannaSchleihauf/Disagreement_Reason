# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 1: We expect that children from all three cultural
# settings are more likely to give reasons in the disagreement relative to the
# agreement condition (i.e., we expect no cultural differences in Study 1).
# In this analysis we look at the NUMBER OF REASONS GIVEN

# September 2023

############################################################################
# PACKAGES & FUNCTIONS
############################################################################
# library("readxl")
library("xlsx")
library("tidyverse")
library("emmeans")
library("effects")
library("ggdist")
library("lme4")
library("directlabels")
library("ggsci")
library("geomtextpath")

source("./functions/diagnostic_fcns.r")
source("./functions/glmm_stability.r")
source("./functions/drop1_para.r")
source("./functions/boot_glmm.r")

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
    head = TRUE, stringsAsFactors = TRUE, na.strings = c("NA", "n/a")
  )
str(xdata)

############################################################################
# DATA WRANGELING
############################################################################
# Make test trial a factor
xdata$test.trial <-
  as.factor(xdata$test.trial)
# Calculate the exact age
xdata <- xdata %>% mutate(
  exact.age.perceptual.child =
    (as.Date(test.date, "%m/%d/%y") -
       as.Date(birth.date.perceptual.child, "%m/%d/%y")) / 365.25,
  exact.age.testimonial.child =
    (as.Date(test.date, "%m/%d/%y") -
       as.Date(birth.date.testimonial.child, "%m/%d/%y")) / 365.25,
  exact.age.perceptual.child =
    ifelse(is.na(exact.age.perceptual.child),
           age.perceptual.child,
           exact.age.perceptual.child
    ),
  exact.age.testimonial.child =
    ifelse(is.na(exact.age.testimonial.child),
           age.testimonial.child,
           exact.age.testimonial.child
    ),
  exact.age.children.dyad =
    (exact.age.perceptual.child + exact.age.testimonial.child) / 2
)

# Create new variable with max reason giving number = 2
xdata <-
  xdata %>%
  mutate(
    overall.nr.reasons =
      nr.of.perceptual.reasons + nr.of.testimonial.reasons,
    overall.nr.reasons.reli =
      nr.of.perceptual.reasons.reli + nr.of.testimonial.reasons.reli)

range(xdata$overall.nr.reasons, na.rm = T)

############################################################################
# DO ANY TRIALS NEED TO BE EXCLUDED?
############################################################################
max(xdata$duration.inaudibles) #no trials need to be excluded
#how many inaudible seconds did we have on average
tapply(xdata$duration.inaudibles, xdata$culture, mean)
tapply(xdata$duration.inaudibles, xdata$culture, sd)

############################################################################
# INTERRATRER RELIABILITY
############################################################################
# extract only the rows for which we have a reliability coding
interrater <-
  xdata %>%
  filter(!is.na(overall.nr.reasons.reli))
nrow(interrater) / nrow(xdata) # percentage of reliability coding 25%

# number of dyads for which the reliability coding was not the same
sum(interrater$overall.nr.reasons.reli !=
      interrater$overall.nr.reasons)

# Possion Kappa
calculatePoissonAgreement <- function(counts1, counts2) {
  # Ensure the input vectors are of the same length
  if(length(counts1) != length(counts2)) {
    stop("The two sets of counts must have the same length.")
  }
  # Calculate the expected mean for each count as the average of the two observations
  expectedMean <- (counts1 + counts2) / 2
  # Calculate the absolute differences from the expected mean
  diff1 <- abs(counts1 - expectedMean)
  diff2 <- abs(counts2 - expectedMean)
  # Calculate the total difference
  totalDiff <- sum(diff1 + diff2)
  # Calculate a simple similarity score as inverse of total difference normalized by the number of observations
  # This is a basic and illustrative metric of agreement
  similarityScore <- 1 - (totalDiff / (2 * sum(expectedMean)))
  return(similarityScore)
}

agreementScore <- calculatePoissonAgreement(interrater$overall.nr.reasons.reli,
                                            interrater$overall.nr.reasons)


############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "overall.nr.reasons ~ condition*culture*exact.age.children.dyad +
                   gender.dyad + test.trial",
  re = "(1|dyad.nr)",
  data = xdata
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)
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
# FITTING THE MODEL AS PREREGISTEred_add
############################################################################
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000000))

full_add <- glmer(overall.nr.reasons ~
                condition * culture * z.age + test.trial + gender.dyad +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = negative.binomial(theta = 1.2)
)
red_add <- glmer(overall.nr.reasons ~
               (condition + culture + z.age)^2 + test.trial + gender.dyad +
               (1 + test.trial.2.code || dyad.nr),
             data = t.data, control = contr,
             family = negative.binomial(theta = 1.2)
)
red_add2 <- glmer(overall.nr.reasons ~
                (condition + culture)^2 + z.age + test.trial + gender.dyad +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = negative.binomial(theta = 1.2)
)
red_add3 <- glmer(overall.nr.reasons ~
                condition + (culture + z.age)^2 + test.trial + gender.dyad +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = negative.binomial(theta = 1.2)
)
main_add <- glmer(overall.nr.reasons ~
                (condition + culture) + z.age + test.trial + gender.dyad +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = negative.binomial(theta = 1.2)
)
null_add <- glmer(overall.nr.reasons ~ 1 +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = negative.binomial(theta = 1.2)
)

summary(full_add)$varcor # the random slope can not be estimated
ranef.diagn.plot(full_add)
xx = round(summary(full_add)$coefficients, 3)
xx

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full_add)

library(car)
vif(main_add)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full_add, contr = contr, use = c("dyad.nr"))
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
round(anova(full_add, null_add, test = "Chisq"), 3)
round(drop1(full_add, test = "Chisq"), 3)
round(drop1(red_add, test = "Chisq"), 3)
round(drop1(red_add2, test = "Chisq"), 3)
round(drop1(red_add3, test = "Chisq"), 3)
round(drop1(main_add, test = "Chisq"), 3)

plot(effect("culture:z.age", red_add), type = "response")

# Coefficients of the full_add model
round(summary(full_add)$coefficients, 3)

emm1 <- emmeans(red_add2, ~ condition * culture)
summary(emm1, type = "response")
summary(pairs(regrid(emm1))[c(7, 9, 14, 1, 10, 15)],
        type = "response", adjust = "BH"
)

############################################################################
# PLOTTING PREPERATION
############################################################################
# Long data set for plotting
# ydata <-
#   xdata %>%
#   pivot_longer(
#     cols = c(
#       perceptual.reason.given.dyad.nr,
#       testimonial.reason.given.dyad.nr
#     ),
#     names_to = "type.reason",
#     values_to = "reason.given"
#   )

############################################################################
# BOOTSTRAP full_add MODEL
############################################################################
boot_full_add <-
  boot.glmm.pred(
    model.res = full_add, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("condition", "culture", "z.age", "test.trial", "gender.dyad")
  )
round(boot_full_add$ci.estimates, 3)
as.data.frame(round(boot_full_add$ci.estimates, 3))
m.stab.plot(round(boot_full_add$ci.estimates, 3))
boot_full_add$ci.predicted

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
  glmer(overall.nr.reasons ~
          (condition*culture) + z.age + test.trial.2.code +
          (1 + test.trial.2.code || dyad.nr),
        data = t.data, control = contr, family = poisson
  )

boot_res_condition_culture.2 <-
  boot.glmm.pred(
    model.res = plot.model.culture.condition, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("condition", "culture")
  )

# Aggregated data for plotting
ydata.agg <- xdata %>%
  group_by(dyad.nr, condition, culture) %>%
  summarise(mean.resp = mean(overall.nr.reasons, na.rm = T)) %>%
  ungroup() %>%
  na.omit()

exp1_plot_condition.culture <-
  ggplot(
    data = ydata.agg,
    aes(x = factor(condition,
                   levels = c("agreement", "disagreement")
    ), y = mean.resp)
  ) +
  geom_point(
    data = ydata.agg,
    aes(x = condition, color = condition), size = 1.5,
    alpha = .4, position = position_jitter(w = 0.25, h = 0.025)
  ) +
  geom_violin(
    data = ydata.agg,
    aes(x = condition, y = mean.resp, fill = condition),
    position = position_nudge(x = 0.00),
    alpha = .2, color = NA
  ) +
  geom_errorbar(
    data = boot_res_condition_culture.2$ci.predicted,
    aes(
      x = as.numeric(condition), y = fitted,
      ymin = lower.cl, ymax = upper.cl, color = condition
    ),
    width = 0.1, linewidth = 1
  ) +
  geom_point(
    data = boot_res_condition_culture.2$ci.predicted,
    aes(x = as.numeric(condition), y = fitted, fill = condition),
    size = 3.5
  ) +
  geom_point(
    data = boot_res_condition_culture.2$ci.predicted,
    aes(x = as.numeric(condition), y = fitted, fill = condition),
    size = 3.5
  ) +
  facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_discrete(
    limits = c("agreement", "disagreement"),
    name = element_blank(),
    labels = c("Agreement\ncondition", "Disagreement\ncondition")
  ) +
  scale_y_continuous(
    name = "Average number of reasons given per dyad"
  ) +
  # labs(title = "A) Reason-giving as a function of condition and culture") +
  scale_color_lancet() +
  scale_fill_lancet() +
  theme_light() +
  my_theme +
  theme(axis.text.x = element_text(size = 11),
        legend.position = "none")

exp1_plot_condition.culture



############################################################################
# PLOT  AGE & CULTURE
############################################################################
# Bootstrap
age.culture.model <-
  glmer(overall.nr.reasons ~ condition.disagreement.code +
          culture*z.age + test.trial.2.code +
          (1 + test.trial.2.code || dyad.nr),
        data = t.data, control = contr, family = poisson
  )

boot_res_age_culture <-
  boot.glmm.pred(
    model.res = age.culture.model, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("culture", "z.age")
  )

# Aggregate Data
xdata$z.age <- scale(xdata$exact.age.children.dyad)
ydata.agg.age <- xdata %>%
  group_by(dyad.nr, culture, z.age) %>%
  summarise(mean.resp = mean(overall.nr.reasons, na.rm = T)) %>%
  ungroup()

# Plot Age
library("viridis")
exp1_plot_age <-
  ggplot() +
  geom_point(
    data = ydata.agg.age,
    aes(x = z.age, y = mean.resp, color = culture),
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.00, h = 0.015)
  ) +
geom_ribbon(
  data = boot_res_age_culture$ci.predicted,
  aes(ymin = lower.cl, ymax = upper.cl, x = z.age, fill = culture),
  alpha = 0.2
) +
geom_textline(
  data = boot_res_age_culture$ci.predicted,
  aes(x = z.age, y = fitted, group = culture,
      label = culture, color = culture),
  hjust = c(0.2), linewidth = 0.5,
) +
  #facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_continuous(
    name = "Age",
    breaks = c(
      (5 - mean(xdata$exact.age.children.dyad)) /
        sd(xdata$exact.age.children.dyad),
      (6 - mean(xdata$exact.age.children.dyad)) /
        sd(xdata$exact.age.children.dyad),
      (7 - mean(xdata$exact.age.children.dyad)) /
        sd(xdata$exact.age.children.dyad),
      (8 - mean(xdata$exact.age.children.dyad)) /
        sd(xdata$exact.age.children.dyad),
      (9 - mean(xdata$exact.age.children.dyad)) /
        sd(xdata$exact.age.children.dyad),
      (10 - mean(xdata$exact.age.children.dyad)) /
        sd(xdata$exact.age.children.dyad)
    ),
    labels = c(5, 6, 7, 8, 9, 10)
  ) +
  scale_y_continuous(
    name = "Average number of reasons given per dyad"
  ) +
  scale_color_jco() +
  scale_fill_jco() +
  theme_light() +
  my_theme +
  theme(legend.position = "none")
exp1_plot_age

exp1_plot_age +
  scale_color_manual(values=c("#0073C2FF", "#FFA319FF", "#868686FF")) +
  scale_fill_manual(values=c("#0073C2FF", "#FFA319FF", "#868686FF"))

############################################################################
# PLOT TRIAL
############################################################################
ydata$z.trial <- scale(as.numeric(ydata$test.trial))
ydata.agg.trial <- ydata %>%
  group_by(dyad.nr, condition, test.trial) %>%
  summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
  ungroup()

trial.model.culture <-
  glmer(overall.nr.reasons ~ condition.disagreement.code +
          culture +
          z.age + test.trial +
          (1 + test.trial.2.code || dyad.nr),
        data = t.data, control = contr, family = binomial(link = "logit")
  )

trial.model <-
  glmer(overall.nr.reasons ~ (condition.disagreement.code +
                      culture.Kenya.code + culture.USA.code) +
          z.age + test.trial +
          (1 + test.trial.2.code || dyad.nr),
        data = t.data, control = contr, family = binomial(link = "logit")
  )

boot_res_trial_culture <-
  boot.glmm.pred_add(
    model.res = trial.model, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("test.trial", "culture")
  )

boot_res_trial <-
  boot.glmm.pred_add(
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
    data = boot_res_trial$ci.pred_addicted,
    aes(
      x = as.numeric(test.trial), y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "gray40",
    width = 0.1, size = 1
  ) +
  geom_point(
    data = boot_res_trial$ci.pred_addicted,
    aes(x = as.numeric(test.trial), y = fitted),
    color = "black", size = 3.5
  ) +
  #facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_discrete(
    limits = c("1", "2"),
    name = "Test Trial",
    labels = c("Test trial 1", "Test trial 2")
  ) +
  scale_y_continuous(
    name = "Proportion of reasons given",
    labels = scales::percent
  ) +
  scale_color_lancet() +
  scale_fill_lancet() +
  # labs(title = "Reason-giving as a function of test trial") +
  theme_light() +
  my_theme
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
                    hjust = -1, # hjust = 0 for left alignment
                    vjust = 2.5, # hjust = 0 for left alignment
                    ncol = 2, nrow = 1
)
figure
ggexport(figure, filename = "figure.png", width = 713, height = 993, res = 400)



save.image("RData_files/additional_analyses.RData")
load("RData_files/additional_analyses.RData")

ftable(box.chosen ~ dominant.child.1, xdata)
ftable(condition ~ box.chosen, xdata)
levels(xdata$dominant.child.1)
