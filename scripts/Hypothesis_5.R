# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 5: DV: meta-talk per dyad (yes/no)
# H5A: We expect that children will engage in more meta-talk (reasoning about
# reasons) in the disagreement than in the agreement condition.
# H5B: In the disagreement condition, we also expect older children to be more
# likely to engage in meta-talk than younger children.

# December 2023

############################################################################
# PACKAGES & FUNCTIONS
############################################################################
library("tidyverse")
library("effects")
library("lme4")

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
    head = TRUE, stringsAsFactors = TRUE
  )
str(xdata)

# How many dyads do we have per culture
table(xdata$culture) / 2

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
       as.Date(birth.date.testimony.child, "%m/%d/%y")) / 365.25,
  exact.age.perceptual.child =
    ifelse(is.na(exact.age.perceptual.child),
           age.perceptual.child,
           exact.age.perceptual.child
    ),
  exact.age.testimonial.child =
    ifelse(is.na(exact.age.testimonial.child),
           age.testimony.child,
           exact.age.testimonial.child
    ),
  exact.age.children.dyad =
    (exact.age.perceptual.child + exact.age.testimonial.child) / 2
)

# Creating a meta-talk variable that entails both types of metatalk codings
xdata <-
  xdata %>%
  mutate(
    metatalk.dyad.combined = ifelse(
      metatalk.dyad.explicit == "yes", "yes", "no")) %>%
  mutate(metatalk.dyad.combined =
           ifelse(metatalk.dyad.2 == "yes", "yes", metatalk.dyad.combined))



############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "metatalk.dyad.combined ~ condition*culture*exact.age.children.dyad +
                   gender.dyad + test.trial",
  re = "(1|dyad.nr)",
  data = xdata
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

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

table(xdata$metatalk.dyad.combined, xdata$age.group, xdata$culture)
table(xdata$metatalk.dyad.combined, xdata$culture)
table(xdata$metatalk.dyad.combined, xdata$condition)


############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
# contr =
# glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=10000000))
contr <-
glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
# contr <-
#   glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 10000000))

full <- glmer(metatalk.dyad.combined ~ condition * culture * z.age +
                test.trial + gender.dyad +
                (1 | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

library(brms)
full.brm <- brm(as.numeric(metatalk.dyad.combined) - 1  ~ condition * culture * z.age +
                  test.trial + gender.dyad +
                  (1 + test.trial.2.code | dyad.nr),
                data = t.data,
                family = 'bernoulli',
                prior = c(
                  set_prior("normal(0, 1)", class = "sd"),
                  set_prior("normal(0, 2.5)", class = "b")
                )
)
conditional_effects(full.brm)

red <- glmer(metatalk.dyad.combined ~ (condition + culture + z.age)^2  +
               test.trial + gender.dyad +
               (1 | dyad.nr),
             data = t.data, control = contr,
             family = binomial(link = "logit")
)
red2 <- glmer(metatalk.dyad.combined ~ (condition + culture)^2 + z.age  +
                test.trial + gender.dyad +
                (1 | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)
red3 <- glmer(metatalk.dyad.combined ~ condition + (culture + z.age)^2 +
                test.trial + gender.dyad +
                (1 | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)
main <- glmer(metatalk.dyad.combined ~ (condition + culture) +  z.age  +
                test.trial + gender.dyad +
                (1 | dyad.nr) ,
              data = t.data, control = contr,
              family = binomial(link = "logit")
)
null <- glmer(metatalk.dyad.combined ~ 1 +
                (1 | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

summary(full)$varcor #the random slope can not be estimated
ranef.diagn.plot(full)
round(summary(full)$coefficients, 3)

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

plot(effect("culture:z.age", full), type = "response")
plot(effect("condition", red2), type = "response")

############################################################################
# PLOTTING PREPERATION
############################################################################
# Long data set for plotting
ydata <-
  xdata %>%
  pivot_longer(
    cols = c(perceptual.reason.given.dyad.nr,
             testimonial.reason.given.dyad.nr),
    names_to = "type.reason",
    values_to = "reason.given"
  )

############################################################################
# SETTING PLOT THEME
############################################################################
my_theme <-
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification='left',
    legend.background = element_blank(),
    #legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size=15),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13), #, family = "Arial Narrow"),
    strip.text.x = element_text(size = 15),

    #text = element_text(family = "Arial Narrow"),

    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    ))

# Boot straps for plotting
red.plot <- glmer(metatalk.dyad.combined ~ condition * culture + z.age +
                test.trial.2.code + gender.dyad.male.code +
                (1 | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

boot_red.plot <-
  boot.glmm.pred(
    model.res = red.plot, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("condition", "culture")
  )

# Boot straps for plotting without interaction
plot.model.culture.condition <-
  glmer(resp_mat ~ (condition + culture) + z.age + test.trial.2.code +
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
                   levels = c("agreement", "disagreement")), y = mean.resp)
  ) +
  geom_point(
    data = ydata.agg %>% filter(condition == "agreement"),
    aes(x = condition), color = "green4", size = 1.5,
    alpha = .4, position = position_jitter(w = 0.25, h = 0.025)
  ) +
  geom_point(
    data = ydata.agg %>% filter(condition == "disagreement"),
    aes(x = condition), color = "red3", size = 1.5,
    alpha = .4, position = position_jitter(w = 0.25, h = 0.025)
  ) +
  geom_violin(
    data = ydata.agg %>% filter(condition == "agreement"),
    aes(x = condition, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "chartreuse3", alpha = .2, color = NA
  ) +
  geom_violin(
    data = ydata.agg %>% filter(condition == "disagreement"),
    aes(x = condition, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "orangered3", alpha = .2, color = NA
  ) +
  geom_errorbar(
    data = boot_res_condition_culture$ci.predicted %>%
      filter(condition == "agreement"),
    aes(
      x = as.numeric(condition), y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "green4",
    width = 0.1, linewidth = 1
  ) +
  geom_errorbar(
    data = boot_res_condition_culture$ci.predicted %>%
      filter(condition == "disagreement"),
    aes(
      x = as.numeric(condition), y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "red3",
    width = 0.1, linewidth = 1
  ) +
  geom_point(
    data = boot_res_condition_culture$ci.predicted %>%
      filter(condition == "agreement"),
    aes(x = as.numeric(condition), y = fitted),
    color = "green4", size = 3.5
  ) +
  geom_point(
    data = boot_res_condition_culture$ci.predicted %>%
      filter(condition == "disagreement"),
    aes(x = as.numeric(condition), y = fitted),
    color = "red3", size = 3.5
  ) +
  facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_discrete(
    limits = c("agreement", "disagreement"),
    name = element_blank(),
    labels = c("Agreement\ncondition", "Disagreement\ncondition")
  ) +
  scale_y_continuous(
    name = "Proportion of reasons given",
    labels = scales::percent
  ) +
  #labs(title = "A) Reason-giving as a function of condition and culture") +
  #theme_light() +
  my_theme +
  theme(axis.text.x = element_text(size = 11))

exp1_plot_condition.culture

#Checking whethe rthe estimates are similar to the means
tapply(ydata.agg$mean.resp, list(ydata.agg$culture, ydata.agg$condition), FUN = mean)

