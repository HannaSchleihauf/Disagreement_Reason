library("lme4")
#library("readxl")
library("xlsx")
library("tidyverse")

library("emmeans")
library("effects")
library("ggdist")

source("./functions/diagnostic_fcns.r")
source("./functions/glmm_stability.r")
source("./functions/drop1_para.r")
source("./functions/boot_glmm.r")


#us & kenya data
xdata <-  read.csv(file="./Disagreement_Reason_Coding_Sheet_All_Three.csv", 
                   head = TRUE, stringsAsFactors = TRUE)
str(xdata)

xdata_disagreement <- 
xdata %>%
filter(condition == "disagreement")

xdata_disagreement$box.chosen <- 
 droplevels(xdata_disagreement$box.chosen)

ftable(box.chosen ~ 
        culture + 
        age.group, xdata_disagreement)

xdata_disagreement$test.trial <- 
 as.factor(xdata_disagreement$test.trial)

xdata_disagreement <- xdata_disagreement %>% 
 mutate(
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

xdata_disagreement$box.chosen <- 
 relevel(xdata_disagreement$box.chosen, ref = "testimonial")

## Prepare data for model fitting------------------------------------
xx.fe.re=fe.re.tab(fe.model = 
                    "box.chosen ~ culture*exact.age.children.dyad + 
                    test.trial",
                   re = "(1|dyad.nr)", 
                   data = xdata_disagreement
)
xx.fe.re$summary
t.data=xx.fe.re$data 
str(t.data) 

t.data$box.chosen.nr <- 
 as.numeric(t.data$box.chosen) - 1

## Center dummy variables (necessary for the random effects in the model and potentially plotting)
t.data$culture.USA.code = t.data$culture.USA - mean(t.data$culture.USA)
t.data$culture.Kenya.code = t.data$culture.Kenya - mean(t.data$culture.Kenya)
t.data$test.trial.2.code = t.data$test.trial.2 - mean(t.data$test.trial.2)
t.data$z.age = scale(t.data$exact.age.children.dyad)

## Fitting the model as pre-registered
contr <- 
 # glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
 glmerControl(optimizer="nlminbwrap", optCtrl=list(maxfun=10000000))


full <- glmer(box.chosen.nr ~ culture*z.age + 
               (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr, 
              family = binomial(link = "logit"))
main <- glmer(box.chosen.nr ~ culture + z.age + 
               (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr, 
              family = binomial(link = "logit"))
null <- glmer(box.chosen.nr ~ 1 + 
               (1 + test.trial.2.code || dyad.nr),
              data = t.data, control = contr, 
              family = binomial(link = "logit"))


summary(full)$varcor
ranef.diagn.plot(full)
round(summary(full)$coefficients, 3)

## Assumptions 
overdisp.test(full)
library(car)
vif(main)
## Checking model stability
m.stab.b = 
 glmm.model.stab(model.res=full, contr=contr, use=c("dyad.nr")) 
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))

m.stab.b$detailed <- m.stab.b$detailed  %>%
 filter(warnings == "none")
cbind(apply(m.stab.b$detailed,2,min),
      apply(m.stab.b$detailed,2,max))


## Model comparison
round(anova(full, null, test="Chisq"), 3)
round(drop1(full, test="Chisq"), 3)
round(drop1(red, test="Chisq"), 3)
round(drop1(main, test="Chisq"), 3)

plot(effect("culture:z.age", full), type = "response")

