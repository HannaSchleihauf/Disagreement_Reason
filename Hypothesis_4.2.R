# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 5: Box chosen

# September 2023

############################################################################
# PACKAGES & FUNCTIONS
############################################################################
library("lme4")
library("xlsx")
library("tidyverse")
library("emmeans")
library("effects")
library("ggdist")
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
xdata <- read.csv(
  file = "./data/Disagreement_Reason_Coding_Sheet_All_Three.csv",
  head = TRUE, stringsAsFactors = TRUE
)
str(xdata)

############################################################################
# DATA WRANGELING
############################################################################
xdata$test.trial <-
  as.factor(xdata$test.trial)
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
xdata$z.age <- scale(xdata$exact.age.children.dyad)
xdata$box.chosen <- relevel(xdata$box.chosen, ref = "testimonial")
ftable(box.chosen ~ culture + condition, xdata)

############################################################################
# Disagreement analysis
############################################################################

xdata_disagreement <- xdata %>%
  filter(condition == "disagreement")
xdata_disagreement$box.chosen <- droplevels(xdata_disagreement$box.chosen)

# MODEL PREPERATION
xx.fe.re <- fe.re.tab(
  fe.model =
    "box.chosen ~ culture*exact.age.children.dyad + dominant.child.1 +
                   gender.dyad + test.trial",
  re = "(1|dyad.nr)",
  data = xdata_disagreement
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

# Center dummy variables
t.data$culture.USA.code <-
  t.data$culture.USA - mean(t.data$culture.USA)
t.data$culture.Kenya.code <-
  t.data$culture.Kenya - mean(t.data$culture.Kenya)
t.data$gender.dyad.male.code <-
  t.data$gender.dyad.male - mean(t.data$gender.dyad.male)
t.data$test.trial.2.code <-
  t.data$test.trial.2 - mean(t.data$test.trial.2)
t.data$z.age <- scale(t.data$exact.age.children.dyad)

# FITTING THE MODEL AS PREREGISTERED
# contr =
# glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=10000000))
# contr <-
# glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
contr <-
  glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 10000000))

full <- glmer(box.chosen ~ culture * dominant.child.1 +
                (1 + test.trial.2.code | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

red <- glmer(box.chosen ~ (culture + dominant.child.1)^2 +
                (1 + test.trial.2.code | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

main <- glmer(box.chosen ~ (culture + dominant.child.1) +
                (1 + test.trial.2.code | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

null <- glmer(box.chosen ~ 1 +
                (1 + test.trial.2.code | dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

summary(full)$varcor
ranef.diagn.plot(full)
round(summary(full)$coefficients, 3)

# CHECKING ASSUMPTIONS
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

# MODEL COMPARISONS
round(anova(full, null, test = "Chisq"), 3)
round(drop1(full, test = "Chisq"), 3)
round(drop1(main, test = "Chisq"), 3)

plot(effect("culture:dominant.child.1", full), type = "response")
plot(effect("culture", main), type = "response")

# Coefficients of the full model
summary(full)$coefficients

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

############################################################################
# Plot Culture*Age Effect
############################################################################

ydata.agg.box.chosen.dominant <- xdata_disagreement %>%
  group_by(dyad.nr, culture, dominant.child.1) %>%
  summarise(mean.resp =
              mean((as.numeric(box.chosen) - 1), na.rm = T)) %>%
  na.omit() %>%
  ungroup()

boot_res_full <-
  boot.glmm.pred(
    model.res = full, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("culture", "dominant.child.1")
  )

plot_box.chosen.dominant <-
  ggplot() +
  geom_point(
    data = ydata.agg.box.chosen.dominant,
    aes(x = dominant.child.1, y = mean.resp, color = culture),
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.15, h = 0.03)
  ) +
  geom_violin(
    data = ydata.agg.box.chosen.dominant,
    aes(x = dominant.child.1, y = mean.resp, fill = culture),
    alpha = .2, color = NA
  ) +
  geom_errorbar(
    data = boot_res_full$ci.predicted,
    aes(
      x = dominant.child.1, y = fitted,
      ymin = lower.cl, ymax = upper.cl, color = culture
    ),
    width = 0.1, linewidth = 1
  ) +
  geom_point(
    data = boot_res_full$ci.predicted,
    aes(x = dominant.child.1, y = fitted, color = culture),
    size = 3.5
  ) +
  facet_wrap(~culture, nrow = 1, switch = "x") +
  scale_x_discrete(
    limits = c("testimonial", "perceptual"),
    name = "Evidence received by the dominant child",
    labels = c("testimonial", "perceptual")
  ) +
  scale_y_continuous(
    name = "Proportion of perceptual box choices",
    labels = scales::percent
  ) +
  # labs(title = "Interaction effect of age and culture") +
  theme_classic() +
  my_theme +
  theme(legend.position = "none")
plot_box.chosen.dominant
