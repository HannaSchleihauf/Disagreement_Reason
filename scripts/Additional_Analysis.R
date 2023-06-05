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



# xdata <- xdata%>%
#  filter(excluded !=1)

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

xdata <- xdata %>% mutate(
 overall.nr.reasons = nr.of.perceptual.reasons + nr.of.testimonial.reasons)

hist(xdata$overall.nr.reasons)


## Prepare data for model fitting------------------------------------
xx.fe.re=fe.re.tab(fe.model = 
                    "nr.reason.dyad.max.2 ~ condition*culture*exact.age.children.dyad + 
                   gender.dyad + test.trial",
                   re = "(1|dyad.nr)", 
                   other.vars = "nr.reason.not.given.dyad.max.2", 
                   data = xdata
)
xx.fe.re$summary
t.data=xx.fe.re$data 
str(t.data) 

resp_mat <- 
 cbind(t.data$nr.reason.dyad.max.2,
       t.data$nr.reason.not.given.dyad.max.2)
t.data$olre <- 
 as.factor(1:nrow(t.data))


## Center dummy variables (necessary for the random effects in the model and potentially plotting)
t.data$condition.disagreement.code = t.data$condition.disagreement - mean(t.data$condition.disagreement)
t.data$culture.USA.code = t.data$culture.USA - mean(t.data$culture.USA)
t.data$culture.Kenya.code = t.data$culture.Kenya - mean(t.data$culture.Kenya)
t.data$gender.dyad.male.code = t.data$gender.dyad.male - mean(t.data$gender.dyad.male)
t.data$test.trial.2.code = t.data$test.trial.2 - mean(t.data$test.trial.2)
t.data$z.age = scale(t.data$exact.age.children.dyad)

## Fitting the model as pre-registered
contr <- 
 # glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
 glmerControl(optimizer="nlminbwrap", optCtrl=list(maxfun=10000000))


table(xdata$metatalk.dyad.explicit)
