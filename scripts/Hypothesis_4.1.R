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

full <- glmer(box.chosen ~ culture * z.age +
  (1 + test.trial.2.code | dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)

main <- glmer(box.chosen ~ (culture + z.age) +
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

plot(effect("culture:z.age", full), type = "response")
plot(effect("culture", main), type = "response")

# Coefficients of the full model
summary(full)$coefficients

############################################################################
# SETTING PLOT THEME
############################################################################
my_theme <-
  theme(
    legend.position = "none",
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
# Plot Culture Effect
############################################################################
levels(xdata_disagreement$box.chosen)

ydata.agg.box.chosen <- xdata_disagreement %>%
  group_by(dyad.nr, culture) %>%
  summarise(mean.resp = mean((as.numeric(box.chosen) - 1), na.rm = T)) %>%
  ungroup()

main <- glmer(box.chosen ~ (culture + z.age) +
  (1 + test.trial.2.code | dyad.nr),
data = t.data, control = contr,
family = binomial(link = "logit")
)

boot_res_main <-
  boot.glmm.pred(
    model.res = main, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("culture")
  )

plot_box.chosen <-
  ggplot() +
  geom_point(
    data = ydata.agg.box.chosen,
    aes(x = culture, y = mean.resp, color = culture),
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.2, h = 0.05)
  ) +
  geom_violin(
    data = ydata.agg.box.chosen,
    aes(x = culture, y = mean.resp, fill = culture),
    alpha = .2, linewidth = 0
  ) +
  geom_errorbar(
    data = boot_res_main$ci.predicted,
    aes(
      x = culture, y = fitted,
      ymin = lower.cl, ymax = upper.cl,
      color = culture
    ),
    width = 0.1, linewidth = 1
  ) +
  geom_point(
    data = boot_res_main$ci.predicted,
    aes(x = culture, y = fitted, color = culture),
    size = 2.5
  ) +
  scale_y_continuous(
    name = "Proportion of perceptual box choices",
    labels = scales::percent
  ) +
  theme_classic() +
  my_theme
plot_box.chosen

############################################################################
# Plot Culture*Age Effect
############################################################################

ydata.agg.box.chosen.age <- xdata_disagreement %>%
  group_by(dyad.nr, culture, z.age) %>%
  summarise(mean.resp =
              mean((as.numeric(box.chosen) - 1), na.rm = T)) %>%
  ungroup()

boot_res_full <-
  boot.glmm.pred(
    model.res = full, excl.warnings = TRUE,
    nboots = 1000, para = FALSE, level = 0.95,
    use = c("culture", "z.age")
  )

plot_box.chosen.age <-
  ggplot() +
  geom_point(
    data = ydata.agg.box.chosen.age,
    aes(x = z.age, y = mean.resp, color = culture),
    size = 1.5, alpha = .4,
    position = position_jitter(w = 0.00, h = 0.015)
  ) +
  geom_ribbon(
    data = boot_res_full$ci.predicted,
    aes(ymin = lower.cl, ymax = upper.cl,
        x = z.age, group = culture, fill = culture),
    alpha = 0.1
  ) +
  geom_textline(
    data = boot_res_full$ci.predicted,
    aes(x = z.age, y = fitted, group = culture,
        color = culture, label = culture),
    hjust = c(0.3)
  ) +
  scale_x_continuous(
    name = "Age",
    breaks = breaks,
    labels = c(4, 5, 6, 7, 8, 9, 10)
  ) +
  scale_y_continuous(
    name = "Proportion of perceptual box choices",
    labels = scales::percent
  ) +
  # labs(title = "Interaction effect of age and culture") +
  theme_classic() +
  my_theme +
  theme(legend.position = "none")
plot_box.chosen.age


############################################################################
# AGREEMENT ANALYSIS
############################################################################

xdata_agreement <- xdata %>%
  filter(condition == "agreement")
xdata_agreement$box.chosen <- droplevels(xdata_agreement$box.chosen)
xdata$box.chosen <- relevel(xdata$box.chosen, ref = "uncued")

xx <-
  xdata_agreement %>%
  group_by(culture, dyad.nr, box.chosen) %>%
  tally()

# MODEL PREPERATION
xx.fe.re <- fe.re.tab(
  fe.model =
    "box.chosen ~ culture*exact.age.children.dyad + dominant.child.1 +
                   gender.dyad + test.trial",
  re = "(1|dyad.nr)",
  data = xdata_agreement
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

table(t.data$box.chosen, t.data$culture)

# Since children in China ALWAYS choose the cure dbox, we have a complete
# separation issue. To tackle this problem, we create 100 data sets, in
# each of which we randomly change one data point at a time in the cells in
# which we have complete separation. We then fit the model and check
# how many of these 1000 models converged, and average the results over
# all converged models.

# function to keep warnings
library(kyotil)

# determining in which cells we will be changing one data point at a time
to.change <- (1:nrow(t.data))[(t.data$culture == "China")]

# creating empty variables to store the results
n.sims = 100
all.res <- data.frame(
  n.sim = c(1:n.sims),
  to.change.1 = NA,
  to.change.2 = NA,
  # full warning
  full.warnings = NA,
  # all.full.coeffs
  all.full.intercept = NA,
  all.full.cultureKenya = NA,
  all.full.cultureUSA = NA,
  all.full.z.age = NA,
  all.full.z.age.cultureKenya = NA,
  all.full.z.age.cultureUSA = NA,
  # null warnings
  null.warnings = NA,
  # test.full.null
  test.full.null.Chisq = NA,
  test.full.null.Df = NA,
  test.full.null.Pr..Chisq = NA,
  # test.2.way.int
  test.2.way.int.Chisq = NA,
  test.2.way.int.Chi.Df = NA,
  test.2.way.int.Pr..Chisq = NA,
  test.2.way.int.n.opt.warnings = NA,
  test.2.way.int.n.fun.warnings = NA,
  # red warnings
  main.warnings = NA,
  # red coefficients
  all.main.intercept = NA,
  all.main.cultureKenya = NA,
  all.main.cultureUSA = NA,
  all.main.z.age = NA,
  # main model comparisons
  test.main.z.age.Chisq = NA,
  test.main.z.age.Chi.Df = NA,
  test.main.z.age.Pr..Chisq = NA,
  test.main.z.age.n.opt.warnings = NA,
  test.main.z.age.n.fun.warnings = NA,
  #
  test.main.culture.Chisq = NA,
  test.main.culture.Chi.Df = NA,
  test.main.culture.Pr..Chisq = NA,
  test.main.culture.n.opt.warnings = NA,
  test.main.culture.n.fun.warnings = NA,
  # assumptions
  overdispersion.test.full = NA,
  colliniarity.test.z.age = NA,
  colliniarity.test.culture = NA
)

boot.full.values <- c()
boot.plot.1.values <- c()
boot.plot.2.values <- c()
boot.full.estimates <- c()
boot.plot.1.estimates <- c()
boot.plot.2.estimates <- c()

# starting the loop
contr <-
  glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 10000000))

t.data$box.chosen <-
  relevel(t.data$box.chosen, ref = "uncued")
t.data$new.resp <-
  as.numeric(t.data$box.chosen) - 1

for (i in 16:length(to.change)) { # i=16
  t.data$new.resp <-
    as.numeric(t.data$box.chosen) - 1
  set.seed(i)
  xx <- to.change[i]
  xx <- sample(to.change, 2, replace = FALSE)
  new.1 <- xx[1]
  new.2 <- xx[2]
  t.data$new.resp[new.1] <- ifelse(t.data$new.resp[new.1] == 1, 0, 1)
  t.data$new.resp[new.2] <- ifelse(t.data$new.resp[new.2] == 1, 0, 1)
  all.res$to.change.1[i] <- new.1
  all.res$to.change.2[i] <- new.2
  # Full model
  full1 <- keepWarnings(glmer(new.resp ~ culture * z.age +
    (1 + test.trial.2.code | dyad.nr),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  ))
  if (length(full1$warnings) == 0) {
    full <- full1$value
    all.res$full.warnings[i] <- "no"
    all.res$all.full.intercept[i] <-
      summary(full)$coefficients["(Intercept)", 1]
    all.res$all.full.cultureKenya[i] <-
      summary(full)$coefficients["cultureKenya", 1]
    all.res$all.full.cultureUSA[i] <-
      summary(full)$coefficients["cultureUSA", 1]
    all.res$all.full.z.age[i] <-
      summary(full)$coefficients["z.age", 1]
    all.res$all.full.z.age.cultureKenya[i] <-
      summary(full)$coefficients["cultureKenya:z.age", 1]
    all.res$all.full.z.age.cultureUSA[i] <-
      summary(full)$coefficients["cultureUSA:z.age", 1]

    # full-red model comparisons
    tests.full <- drop1p(model.res = full, contr = contr)
    test.2.way.int <-
      as.data.frame(tests.full$drop1.res[2, c(
        "Chisq",
        "Chi.Df",
        "Pr..Chisq.",
        "n.opt.warnings",
        "n.fun.warnings"
      )])
    all.res$test.2.way.int.Chisq[i] <-
      test.2.way.int$Chisq
    all.res$test.2.way.int.Chi.Df[i] <-
      test.2.way.int$Chi.Df
    all.res$test.2.way.int.Pr..Chisq[i] <-
      test.2.way.int$Pr..Chisq.
    all.res$test.2.way.int.n.opt.warnings[i] <-
      test.2.way.int$n.opt.warnings
    all.res$test.2.way.int.n.fun.warnings[i] <-
      test.2.way.int$n.fun.warnings

    # assumptions
    all.res$overdispersion.test.full[i] <-
      overdisp.test(full)[1, c("dispersion.parameter")]

    # boot full model
    boot.full <- boot.glmm.pred(
      model.res = full, excl.warnings = T, nboots = 10,
      para = F, resol = 100, level = 0.95, use = c("culture", "z.age")
    )

    boot.full.values <-
      rbind(boot.full.values, data.frame(
        culture = boot.full$ci.predicted$culture,
        z.age = boot.full$ci.predicted$z.age,
        fitted = boot.full$ci.predicted$fitted,
        lower.cl = boot.full$ci.predicted$lower.cl,
        upper.cl = boot.full$ci.predicted$upper.cl
      ))
    boot.full.estimates <-
      rbind(boot.full.estimates, data.frame(
        term = rownames(boot.full$ci.estimates),
        orig = boot.full$ci.estimates$orig,
        X2.5. = boot.full$ci.estimates$X2.5.,
        X97.5. = boot.full$ci.estimates$X97.5.
      ))
  } else {
    all.res$full.warnings[i] <- "yes"
  }

  # Null model
  null1 <- keepWarnings(glmer(new.resp ~ 1 +
    (1 + test.trial.2.code | dyad.nr),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  ))

  if (length(full1$warnings) == 0 & length(null1$warnings) == 0) {
    null <- null1$value
    all.res$null.warnings[i] <- "no"

    # full null model comparisons
    test.full.null <-
      as.data.frame(anova(null, full,
        test = "Chisq"
      ))["full", c(
        "Chisq",
        "Df",
        "Pr(>Chisq)"
      )]
    all.res$test.full.null.Chisq[i] <-
      test.full.null$Chisq
    all.res$test.full.null.Df[i] <-
      test.full.null$Df
    all.res$test.full.null.Pr..Chisq[i] <-
      test.full.null$`Pr(>Chisq)`
  } else {
    all.res$null.warnings[i] <- "yes"
  }

  # Red model
  main1 <- keepWarnings(glmer(new.resp ~ culture + z.age +
    (1 + test.trial.2.code | dyad.nr),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  ))

  if (length(main1$warnings) == 0) {
    main <- main1$value
    all.res$main.warnings[i] <- "no"

    # check colliniarity
    all.res$colliniarity.test.culture[i] <- vif(main)[1, 3]
    all.res$colliniarity.test.z.age[i] <- vif(main)[2, 3]

    # main coefficients
    all.res$all.main.intercept[i] <-
      summary(main)$coefficients["(Intercept)", 1]
    all.res$all.main.cultureKenya[i] <-
      summary(main)$coefficients["cultureKenya", 1]
    all.res$all.main.cultureUSA[i] <-
      summary(main)$coefficients["cultureUSA", 1]
    all.res$all.main.z.age[i] <-
      summary(main)$coefficients["z.age", 1]

    # main effects
    tests.main <-
      drop1p(model.res = main, contr = contr)

    test.main.z.age <-
      as.data.frame(tests.main$drop1.res[3, c(
        "Chisq",
        "Chi.Df",
        "Pr..Chisq.",
        "n.opt.warnings",
        "n.fun.warnings"
      )])
    all.res$test.main.z.age.Chisq[i] <-
      test.main.z.age$Chisq
    all.res$test.main.z.age.Chi.Df[i] <-
      test.main.z.age$Chi.Df
    all.res$test.main.z.age.Pr..Chisq[i] <-
      test.main.z.age$Pr..Chisq.
    all.res$test.main.z.age.n.opt.warnings[i] <-
      test.main.z.age$n.opt.warnings
    all.res$test.main.z.age.n.fun.warnings[i] <-
      test.main.z.age$n.fun.warnings
    #
    test.main.culture <-
      as.data.frame(tests.main$drop1.res[2, c(
        "Chisq",
        "Chi.Df",
        "Pr..Chisq.",
        "n.opt.warnings",
        "n.fun.warnings"
      )])
    all.res$test.main.culture.Chisq[i] <-
      test.main.culture$Chisq
    all.res$test.main.culture.Chi.Df[i] <-
      test.main.culture$Chi.Df
    all.res$test.main.culture.Pr..Chisq[i] <-
      test.main.culture$Pr..Chisq.
    all.res$test.main.culture.n.opt.warnings[i] <-
      test.main.culture$n.opt.warnings
    all.res$test.main.culture.n.fun.warnings[i] <-
      test.main.culture$n.fun.warnings
  } else {
    all.res$main.warnings[i] <- "yes"
  }

  print(i)
}

## Evaluation of the results -----------------------------------------------
str(all.res)
# how many models did converge
# full
sum(all.res$full.warnings == "no", na.rm = T) # 45 of the models converged
# null
sum(all.res$null.warnings == "no", na.rm = T) # 39 of the models converged
# main
sum(all.res$main.warnings == "no", na.rm = T) # 47 of the models converged

#remove models with overdispersion
all.res <-
  all.res %>%
  filter(overdispersion.test.full < 10 &
         overdispersion.test.full > 0.02 &
          main.warnings == "no" &
          null.warnings == "no" &
          test.full.null.Pr..Chisq != 1.0000 &
          test.2.way.int.Pr..Chisq != 1.0000 &
          test.main.z.age.Pr..Chisq != 1.0000 &
          test.main.culture.Pr..Chisq != 1.0000)

# Means of full-null-comparisons
# Chisq
round(mean(all.res$test.full.null.Chisq, na.rm = T), 10)
round(range(all.res$test.full.null.Chisq, na.rm = T), 10)
# DF
range(all.res$test.full.null.Df, na.rm = T)
# p-value
round(mean(all.res$test.full.null.Pr..Chisq, na.rm = T), 10)
round(range(all.res$test.full.null.Pr..Chisq, na.rm = T), 10)
# percent of significant models
(sum(all.res$test.full.null.Pr..Chisq <= 0.051, na.rm = T) /
       sum(all.res$null.warnings == "no"))

# Assumption tests
# overdispersion parameter
round(mean(all.res$overdispersion.test.full, na.rm = T), 10)
round(range(all.res$overdispersion.test.full, na.rm = T), 10)
# colliniarity
round(mean(all.res$colliniarity.test.z.age, na.rm = T), 10)
round(range(all.res$colliniarity.test.z.age, na.rm = T), 10)
round(mean(all.res$colliniarity.test.culture, na.rm = T), 10)
round(range(all.res$colliniarity.test.culture, na.rm = T), 10)

# Means of two-way interactions
# age  * culture
sum(all.res$test.2.way.z.age.culture.n.opt.warnings == 0 &
      all.res$test.2.way.z.age.culture.n.fun.warnings == 0, na.rm = T)
# Chisq
round(mean(all.res$test.2.way.int.Chisq, na.rm = T), 10)
round(range(all.res$test.2.way.int.Chisq, na.rm = T), 10)
# DF
range(all.res$test.2.way.int.Chi.Df, na.rm = T)
# p-value
round(mean(all.res$test.2.way.int.Pr..Chisq, na.rm = T), 10)
# percent of non-significant models
(sum(all.res$test.2.way.int.Pr..Chisq <= 0.051, na.rm = T) /
       sum(all.res$main.warnings == "no"))
round(range(all.res$test.2.way.int.Pr..Chisq, na.rm = T), 10)

# main effect age group
sum(all.res$test.main.z.age.n.opt.warnings == 0 &
      all.res$test.main.z.age.n.fun.warnings == 0, na.rm = T)
# Chisq
round(mean(all.res$test.main.z.age.Chisq, na.rm = T), 10)
round(range(all.res$test.main.z.age.Chisq, na.rm = T), 10)
# DF
range(all.res$test.main.z.age.Chi.Df, na.rm = T)
# p-value
round(mean(all.res$test.main.z.age.Pr..Chisq, na.rm = T), 10)
# percent of significant models
sum(all.res$test.main.z.age.Pr..Chisq <= 0.051, na.rm = T) /
       sum(all.res$main.warnings == "no")
round(range(all.res$test.main.z.age.Pr..Chisq, na.rm = T), 10)

# main effect culture
sum(all.res$test.main.culture.n.opt.warnings == 0 &
      all.res$test.main.culture.n.fun.warnings == 0, na.rm = T)
# Chisq
round(mean(all.res$test.main.culture.Chisq, na.rm = T), 10)
round(range(all.res$test.main.culture.Chisq, na.rm = T), 10)
# DF
range(all.res$test.main.culture.Chi.Df, na.rm = T)
# p-value
round(mean(all.res$test.main.culture.Pr..Chisq, na.rm = T), 10)
# percent of significant models
(sum(all.res$test.main.culture.Pr..Chisq <= 0.051, na.rm = T) /
    sum(all.res$main.warnings == "no"))
round(range(all.res$test.main.culture.Pr..Chisq, na.rm = T), 10)

ftable(t.data$box.chosen ~ t.data$culture)

# boot.full fitted values and confidence intervals
boot.full <- mapply(FUN = tapply, X = boot.full.values[, -1], MoreArgs = list(INDEX = boot.full.values$culture, FUN = mean)) # fitted
xx <- round(boot.full, 3)
boot.full <- mapply(FUN = tapply, X = boot.full.values[, -1], MoreArgs = list(INDEX = boot.full.values$term, FUN = min)) # lower ci
yy <- round(boot.full, 3)
boot.full <- mapply(FUN = tapply, X = boot.full.values[, -1], MoreArgs = list(INDEX = boot.full.values$term, FUN = max)) # upper ci
zz <- round(boot.full, 3)
