# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 1_alternative: We expect that children from all three cultural
# settings are more likely to give reasons in the disagreement relative to the
# agreement condition (i.e., we expect no cultural differences in Study 1).

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

ydata =
subset(xdata, xdata$reason.giving.dyad == "yes" &
         xdata$testimonial.reason.given == "no" &
         xdata$perceptual.reason.given == "no")

levels(xdata$reason.giving.dyad)
xdata$reason_giving_dyad_nr <-
  as.numeric(xdata$reason.giving.dyad) - 1 #1 is yes


############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "reason_giving_dyad_nr ~ condition*culture*exact.age.children.dyad +
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

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr <-
glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
#glmerControl(optimizer = "nlminbwrap", optCtrl = list(maxfun = 10000000))
contr = glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 10000000))
contr = glmerControl(optimizer = "nloptwrap",
                     optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD",  maxit = 1e9))

full <- glmer(reason_giving_dyad_nr ~ condition * culture * z.age +
                test.trial + gender.dyad +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit"))

red <- glmer(reason_giving_dyad_nr ~ (condition + culture + z.age)^2 +
               test.trial + gender.dyad +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit"))

red2 <- glmer(reason_giving_dyad_nr ~ (condition + culture)^2 + z.age +
                test.trial + gender.dyad +
               (1 + test.trial.2.code || dyad.nr),
             data = t.data, control = contr,
             family = binomial(link = "logit"))

main <- glmer(reason_giving_dyad_nr ~ (condition + culture + z.age) +
                test.trial + gender.dyad +
               (1 + test.trial.2.code || dyad.nr),
             data = t.data, control = contr,
             family = binomial(link = "logit"))

null <- glmer(reason_giving_dyad_nr ~ 1 +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit"))

summary(full)$varcor
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
round(drop1(main, test = "Chisq"), 3)

plot(effect("condition:culture:z.age", full), type = "response")

# Coefficients of the full model
summary(full)$coefficients

emm1 <- emmeans(red2, ~ condition * culture)
summary(emm1, type = "response")
summary(pairs(regrid(emm1))[c(7, 9, 14)],
        type = "response", adjust = "fdr")
