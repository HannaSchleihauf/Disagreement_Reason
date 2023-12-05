# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 3: In the disagreement condition, we expect to find
# differences between cultures and age groups in how often dyads give the
# testimony-based or the perceptual reason until they come to a decision.
# For example, since younger children in China might find authority-based
# reasons more convincing, and the testimony-based reason is an authority-based
# reason, they might repeat the testimony-based reason several times to convince
# their partner.

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
library("glmmTMB")
library("DHARMa")
source("./functions/diagnostic_fcns.r")
source("./functions/glmm_stability.r")
source("./functions/drop1_para.r")
source("./functions/boot_glmm.r")

source("./functions/boot_glmmTMB.r")
source("./functions/glmmTMB_stability.r")

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

############################################################################
# PIVOT_LONGER
############################################################################
nr.data <-
  pivot_longer(
    data = xdata,
    cols = c(
      "nr.of.testimonial.reasons",
      "nr.of.perceptual.reasons"
    ),
    names_to = "type_of_reason",
    values_to = "number_of_reasons"
  )

hist(nr.data$number_of_reasons)

nr.data <- subset(nr.data, nr.data$condition == "disagreement")

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "number_of_reasons ~
                     type_of_reason*culture*exact.age.children.dyad +
                   gender.dyad + test.trial",
  re = "(1|dyad.nr)",
  other.vars = c("age.group"),
  data = nr.data
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

## Center dummy variables
t.data$culture.USA.code <-
  t.data$culture.USA -
  mean(t.data$culture.USA)
t.data$culture.Kenya.code <-
  t.data$culture.Kenya -
  mean(t.data$culture.Kenya)
t.data$gender.dyad.male.code <-
  t.data$gender.dyad.male -
  mean(t.data$gender.dyad.male)
t.data$test.trial.2.code <-
  t.data$test.trial.2 -
  mean(t.data$test.trial.2)
t.data$type_of_reason.nr.of.testimonial.reasons.code <-
  t.data$type_of_reason.nr.of.testimonial.reasons -
  mean(t.data$type_of_reason.nr.of.testimonial.reasons)
t.data$z.age <-
  scale(t.data$exact.age.children.dyad)

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
  glmmTMBControl(
    optCtrl = list(iter.max = 200000, eval.max = 200000),
    profile = FALSE, collect = FALSE
  )

full <- glmmTMB(number_of_reasons ~
  # condition*
  type_of_reason *
  culture *
  z.age +
  (1 | dyad.nr),
data = t.data, control = contr,
family = nbinom1(link = "log")
)

full_resi <- simulateResiduals(full)
plot(full_resi)

red <- glmmTMB(number_of_reasons ~
  # condition*
  (type_of_reason +
    culture +
    z.age)^2 +
  (1 | dyad.nr),
data = t.data, control = contr,
family = nbinom1(link = "log")
)

red2 <- glmmTMB(number_of_reasons ~
  # condition*
  type_of_reason +
  (culture +
    z.age)^2 +
  (1 | dyad.nr),
data = t.data, control = contr,
family = nbinom1(link = "log")
)

main <- glmmTMB(number_of_reasons ~
  # condition*
  (type_of_reason +
    culture +
    z.age) +
  (1 | dyad.nr),
data = t.data, control = contr,
family = nbinom1(link = "log")
)

null <- glmmTMB(number_of_reasons ~
  1 +
  (1 | dyad.nr),
data = t.data, control = contr,
family = nbinom1(link = "log")
)
summary(full)$varcor
ranef.diagn.plot(full)
round(summary(full)$coefficients$cond, 3)

############################################################################
# CHECKING ASSUMPTIONS
############################################################################
overdisp.test(full)
library(car)
vif(main)
## Checking model stability
full.stab <- glmmTMB.stab(
  model.res = full, contr = contr,
  ind.cases = F, para = F,
  data = t.data, use = "dyad.nr",
  n.cores = c("all"), save.path = NULL,
  load.lib = T, lib.loc = .libPaths()
)
full.stab$detailed$converged
round(full.stab$summary[, -1], 3)
m.stab.plot(full.stab$summary[, -1])

############################################################################
# MODEL COMPARISONS
############################################################################
round(anova(full, null, test = "Chisq"), 3)
round(drop1(full, test = "Chisq"), 3)
round(drop1(red, test = "Chisq"), 3)
round(drop1(red2, test = "Chisq"), 3)
round(drop1(main, test = "Chisq"), 3)

plot(effect("culture:z.age", red2), type = "response")


############################################################################
# PLOT CULTURE AGE
############################################################################
plot.red2 <- glmmTMB(number_of_reasons ~
  # condition*
  type_of_reason.nr.of.testimonial.reasons.code +
  (culture +
    z.age)^2 +
  (1 | dyad.nr),
data = t.data, control = contr,
family = nbinom1(link = "log")
)

boot.culture.age <- boot.glmmTMB(
  model.res = plot.red2, data = t.data, excl.non.conv = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("culture", "z.age"),
  contr = contr, circ.var.name = NULL, circ.var = NULL,
  n.cores = c("all-1"), save.path = NULL,
  load.lib = T, lib.loc = .libPaths(), set.all.effects.2.zero = F
)
round(boot.culture.age$ci.estimates$fe, 3)
boot.culture.age$ci.fitted$number_of_reasons <-
  boot.culture.age$ci.fitted$fitted


plot_culture_age <-
  ggplot(
    data = t.data,
    aes(
      x = z.age,
      y = number_of_reasons
    )
  ) +
  geom_point(
    aes(color = type_of_reason),
    data = t.data,
    size = 1.5,
    alpha = .4, position = position_jitter(w = 0.08, h = 0.25)
  ) +
  geom_ribbon(
    data = boot.culture.age$ci.fitted,
    aes(ymin = lower.cl, ymax = upper.cl, x = z.age),
    alpha = 0.5, fill = "grey"
  ) +
  geom_line(
    data = boot.culture.age$ci.fitted,
    aes(x = z.age, y = number_of_reasons),
    color = "darkgrey"
  ) +
  facet_grid(~culture) +
  scale_color_manual(
    values = c("seagreen3", "palevioletred"),
    labels = c(
      "Perceptual reason ",
      "Testimonial reason       "
    ),
    name = NULL
  ) +
  scale_x_continuous(
    name = "Age",
    breaks = c(
      (4 - mean(t.data$exact.age.children.dyad)) /
        sd(t.data$exact.age.children.dyad),
      (5 - mean(t.data$exact.age.children.dyad)) /
        sd(t.data$exact.age.children.dyad),
      (6 - mean(t.data$exact.age.children.dyad)) /
        sd(t.data$exact.age.children.dyad),
      (7 - mean(t.data$exact.age.children.dyad)) /
        sd(t.data$exact.age.children.dyad),
      (8 - mean(t.data$exact.age.children.dyad)) /
        sd(t.data$exact.age.children.dyad),
      (9 - mean(t.data$exact.age.children.dyad)) /
        sd(t.data$exact.age.children.dyad),
      (10 - mean(t.data$exact.age.children.dyad)) /
        sd(t.data$exact.age.children.dyad)
    ),
    labels = c(4, 5, 6, 7, 8, 9, 10)
  ) +
  scale_y_continuous(
    name = "Number a reason has been given",
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  ) +
  labs(title = "Number a reason has been given as a function of culture and age") +
  theme_light() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "lightgrey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    )
  )

plot_culture_age

save.image("RData_files/Hypothesis_4.RData")
