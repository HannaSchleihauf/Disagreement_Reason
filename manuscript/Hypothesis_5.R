# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 2: When only looking at the disagreement condition,
# we expect to find that children’s tendency to give reasons increases with age.

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
options(scipen=999)

############################################################################
# DATA
############################################################################
xdata <-  read.csv(file="./data/Disagreement_Reason_Coding_Sheet_All_Three.csv",
                   head = TRUE, stringsAsFactors = TRUE)
str(xdata)

############################################################################
# DATA WRANGELING
############################################################################
xdata$test.trial <-
  as.factor(xdata$test.trial)
xdata <- xdata %>% mutate(
  exact.age.perceptual.child =
    (as.Date(test.date, "%m/%d/%y") -
       as.Date(birth.date.perceptual.child, "%m/%d/%y"))/ 365.25,
  exact.age.testimonial.child =
    (as.Date(test.date, "%m/%d/%y") -
       as.Date(birth.date.testimony.child, "%m/%d/%y"))/ 365.25,
  exact.age.perceptual.child =
    ifelse(is.na(exact.age.perceptual.child),
           age.perceptual.child,
           exact.age.perceptual.child),
  exact.age.testimonial.child =
    ifelse(is.na(exact.age.testimonial.child),
           age.testimony.child,
           exact.age.testimonial.child),
  exact.age.children.dyad =
    (exact.age.perceptual.child + exact.age.testimonial.child) / 2
)

xdata_disagreement <- xdata%>%
  filter(condition == "disagreement")
xdata_disagreement$box.chosen <- droplevels(xdata_disagreement$box.chosen)

table(xdata_disagreement$culture, xdata_disagreement$box.chosen)

############################################################################
# PREPARE DATAFRAME FOR MODEL FITTING
############################################################################
xx.fe.re <- fe.re.tab(
  fe.model =
    "box.chosen ~ culture*exact.age.children.dyad +
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

t.data$box.chosen <- relevel(t.data$box.chosen, ref = "testimonial")
levels(t.data$box.chosen)
t.data$box.chosen.nr <-
  as.numeric(t.data$box.chosen) - 1 #1 is perceptual

############################################################################
# FITTING THE MODEL AS PREREGISTERED
############################################################################
contr = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=10000000))
#contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

full <- glmer(box.chosen.nr ~ culture * z.age +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

red <- glmer(box.chosen.nr ~ culture + z.age +
                (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr,
              family = binomial(link = "logit")
)

null <- glmer(box.chosen.nr ~ 1 +
               (1 + test.trial.2.code || dyad.nr),
             data = t.data, control = contr,
             family = binomial(link = "logit")
)


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

plot(effect("culture:z.age", full), type = "response")

# Coefficients of the full model
summary(full)$coefficients

emm1 <- emmeans(full, ~ z.age * culture)
summary(emm1, type = "response")
summary(pairs(regrid(emm1)),
        type = "response", adjust = "fdr")
