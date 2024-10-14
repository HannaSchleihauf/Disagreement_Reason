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
library("ggsci")
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
  head = TRUE, stringsAsFactors = TRUE, na.strings = c("NA", "")
)
str(xdata)

?read.csv

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
xdata$z.age <- scale(xdata$exact.age.children.dyad)
xdata$box.chosen <- relevel(xdata$box.chosen, ref = "testimonial")
ftable(box.chosen ~ culture + condition, xdata)

############################################################################
# INTERRATRER RELIABILITY
############################################################################

# extract only the rows for which we have a reliability coding
interrater <-
  xdata %>%
  filter(!is.na(dominant.child.1.reli))

nrow(interrater) / nrow(xdata) # percentage of reliability coding 25%

# number of dyads for which the reliability coding was not the same
sum(interrater$dominant.child.1.reli !=
      interrater$dominant.child.1)

which(interrater$dominant.child.1.reli !=
      interrater$dominant.child.1)

# Cohen's Kappa
library(irr)
kappa_results <-
  kappa2(cbind(
    interrater$dominant.child.1.reli,
    interrater$dominant.child.1
  ), "unweighted")
print(kappa_results)

############################################################################
# Disagreement analysis
############################################################################

xdata_disagreement <- xdata %>%
  filter(condition == "disagreement")
xdata_disagreement$box.chosen <-
  droplevels(xdata_disagreement$box.chosen)
xdata_disagreement <-
  xdata_disagreement %>%
  mutate(followed.evidence.perceived.by.dominant =
           ifelse(dominant.child.1 == box.chosen, 1, 0))

# MODEL PREPERATION
xx.fe.re <- fe.re.tab(
  fe.model =
    "followed.evidence.perceived.by.dominant ~ culture*exact.age.children.dyad +
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

full_hyp4.2 <- glmer(followed.evidence.perceived.by.dominant ~
                culture * z.age +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

main_hyp4.2 <- glmer(followed.evidence.perceived.by.dominant ~
               culture + z.age +
               (1 + test.trial.2.code || dyad.nr),
             data = t.data, control = contr,
             family = binomial(link = "logit")
)

null_hyp4.2 <- glmer(followed.evidence.perceived.by.dominant ~
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

summary(full_hyp4.2)$varcor
ranef.diagn.plot(full_hyp4.2)
round(summary(full_hyp4.2)$coefficients, 3)

# CHECKING ASSUMPTIONS
overdisp.test(full_hyp4.2)
library(car)
vif(main_hyp4.2)
# Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full_hyp4.2, contr = contr, use = c("dyad.nr"))
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
round(anova(full_hyp4.2, null_hyp4.2, test = "Chisq"), 3)
round(drop1(full_hyp4.2, test = "Chisq"), 3)
round(drop1(main_hyp4.2, test = "Chisq"), 3)

plot(effect("culture", main_hyp4.2), type = "response")

# Coefficients of the full_hyp4.2 model
round(summary(full_hyp4.2)$coefficients, 3)

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
    axis.text = element_text(size = 13), #family = "Arial Narrow"),
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
# Plot Culture Effect
############################################################################

ydata.agg.box.chosen.dominant <- xdata_disagreement %>%
  group_by(dyad.nr, culture) %>%
  summarise(mean.resp =
              mean(followed.evidence.perceived.by.dominant, na.rm = T)) %>%
  na.omit() %>%
  ungroup()

boot_res_full_hyp4.2 <-
  boot.glmm.pred(
    model.res = full_hyp4.2, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("culture")
  )
round(boot_res_full_hyp4.2$ci.estimates, 3)
as.data.frame(round(boot_res_full_hyp4.2$ci.estimates, 3))
m.stab.plot(round(boot_res_full_hyp4.2$ci.estimates, 3))
boot_res_full_hyp4.2$ci.predicted


plot_box.chosen.dominant <-
  ggplot() +
  geom_hline(yintercept=0.5, linetype="dashed", color = "gray70") +
  geom_point(
    data = ydata.agg.box.chosen.dominant,
    aes(x = culture, y = mean.resp, color = culture),
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.15, h = 0.03)
  ) +
  geom_violin(
    data = ydata.agg.box.chosen.dominant,
    aes(x = culture, y = mean.resp, fill = culture),
    alpha = .2, color = NA
  ) +
  geom_errorbar(
    data = boot_res_full_hyp4.2$ci.predicted,
    aes(
      x = culture, y = fitted,
      ymin = lower.cl, ymax = upper.cl, color = culture
    ),
    width = 0.1, linewidth = 1
  ) +
  geom_point(
    data = boot_res_full_hyp4.2$ci.predicted,
    aes(x = culture, y = fitted, color = culture),
    size = 3.5
  ) +

  scale_y_continuous(
    name = "Proportion of trials in which dyads followed\nthe evidence received by the dominant child",
    labels = scales::percent
  ) +
  scale_color_jco() +
  scale_fill_jco() +
  # labs(title = "Interaction effect of age and culture") +
  theme_light() +
  my_theme +
  theme(legend.position = "none")
plot_box.chosen.dominant

save.image("RData_files/Hypothesis_4.2.RData")
load("RData_files/Hypothesis_4.2.RData")
