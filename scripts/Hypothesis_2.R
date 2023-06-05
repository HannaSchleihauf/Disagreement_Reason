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

# create new variable with max reason giving number = 2
xdata <-
 xdata %>%
 mutate(
  perceptual.reason.given.dyad.nr =
   ifelse(perceptual.reason.given == "yes", 1, 0)
 ) %>%
 mutate(
  testimonial.reason.given.dyad.nr =
   ifelse(testimonial.reason.given == "yes", 1, 0)
 ) %>%
 mutate(
  nr.reason.dyad.max.2 =
   perceptual.reason.given.dyad.nr + testimonial.reason.given.dyad.nr
 ) %>%
 mutate(nr.reason.not.given.dyad.max.2 = 2 - nr.reason.dyad.max.2)

xdata_disagreement <- xdata%>%
  filter(condition == "disagreement")
str(xdata_disagreement)


## Prepare data for model fitting------------------------------------
xx.fe.re=fe.re.tab(fe.model = 
                    "nr.reason.dyad.max.2 ~ culture*exact.age.children.dyad + 
                   gender.dyad + test.trial",
                   re = "(1|dyad.nr)", 
                   other.vars = "nr.reason.not.given.dyad.max.2", 
                   data = xdata_disagreement
)
xx.fe.re$summary
t.data_disagreement=xx.fe.re$data 
str(t.data_disagreement) 

resp_mat <- 
 cbind(t.data_disagreement$nr.reason.dyad.max.2,
       t.data_disagreement$nr.reason.not.given.dyad.max.2)
t.data_disagreement$olre <- 
 as.factor(1:nrow(t.data_disagreement))


## Center dummy variables (necessary for the random effects in the model and potentially plotting)
t.data_disagreement$culture.USA.code = t.data_disagreement$culture.USA - mean(t.data_disagreement$culture.USA)
t.data_disagreement$culture.Kenya.code = t.data_disagreement$culture.Kenya - mean(t.data_disagreement$culture.Kenya)
t.data_disagreement$gender.dyad.male.code = t.data_disagreement$gender.dyad.male - mean(t.data_disagreement$gender.dyad.male)
t.data_disagreement$test.trial.2.code = t.data_disagreement$test.trial.2 - mean(t.data_disagreement$test.trial.2)
t.data_disagreement$z.age = scale(t.data_disagreement$exact.age.children.dyad)

## Fitting the model as pre-registered
contr <- 
 # glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
 glmerControl(optimizer="nlminbwrap", optCtrl=list(maxfun=10000000))

full <- glmer(resp_mat ~ culture*z.age*test.trial + 
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data_disagreement, control = contr, 
              family = binomial(link = "logit"))
red <- glmer(resp_mat ~ (culture+z.age+test.trial)^2 +
              (1 + test.trial.2.code || dyad.nr) + (1 | olre),
             data = t.data_disagreement, control = contr, 
             family = binomial(link = "logit"))
red2 <- glmer(resp_mat ~ culture + z.age*test.trial + 
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data_disagreement, control = contr, 
              family = binomial(link = "logit"))
main <- glmer(resp_mat ~ culture + z.age + test.trial + 
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data_disagreement, control = contr, 
              family = binomial(link = "logit"))
null <- glmer(resp_mat ~ 1 +
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data_disagreement, control = contr, 
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
round(drop1(red2, test="Chisq"), 3)
round(drop1(main, test="Chisq"), 3)

plot(effect("culture:z.age:test.trial", full), type = "response")

## Plotting Preparations
# Long data set for plotting 
ydata_disagreement <-
 xdata_disagreement %>%
 pivot_longer(
  cols = c(perceptual.reason.given.dyad.nr, testimonial.reason.given.dyad.nr),
  names_to = "type.reason",
  values_to = "reason.given"
 ) 

## Plot Culture Condition 
# Boot straps for plotting

# ydata_disagreement$z.age <- 
#  scale(ydata_disagreement$exact.age.children.dyad)
# ydata.agg.age <- ydata_disagreement %>%
#  group_by(dyad.nr, culture, trest.trial, z.age) %>%
#  summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
#  ungroup()

boot_res_H2 <-
 boot.glmm.pred(
  model.res = full, excl.warnings = TRUE,
  nboots = 1000, para = FALSE, level = 0.95, 
  use = c("condition", "culture", "z.age")
 )

t.data_disagreement$percentage = 
 t.data_disagreement$nr.reason.dyad.max.2 / 2

# Aggregated data for plotting 

 ggplot() +
 
 geom_point(data = t.data_disagreement,
            aes(x = z.age, 
                y = percentage), color = "darkorange", size = 1.5, 
            alpha = .4, position = position_jitter(w = 0.2, h = 0.015)) + 

  geom_ribbon(data = boot_res_H2$ci.predicted, 
              aes(ymin = lower.cl, ymax = upper.cl, x = z.age), 
              alpha = 0.3, fill = "grey") +
  geom_line(data=boot_res_H2$ci.predicted, 
            aes(x=z.age, y=fitted), 
            color = "darkgrey") +
  
 
 facet_grid(culture ~ test.trial) +
  
  scale_x_continuous(name="Age", 
                     breaks = c((4-mean(t.data_disagreement$exact.age.children.dyad))/
                                 sd(t.data_disagreement$exact.age.children.dyad),
                                (5-mean(t.data_disagreement$exact.age.children.dyad))/
                                 sd(t.data_disagreement$exact.age.children.dyad), 
                                (6-mean(t.data_disagreement$exact.age.children.dyad))/
                                 sd(t.data_disagreement$exact.age.children.dyad),
                                (7-mean(t.data_disagreement$exact.age.children.dyad))/
                                 sd(t.data_disagreement$exact.age.children.dyad), 
                                (8-mean(t.data_disagreement$exact.age.children.dyad))/
                                 sd(t.data_disagreement$exact.age.children.dyad), 
                                (9-mean(t.data_disagreement$exact.age.children.dyad))/
                                 sd(t.data_disagreement$exact.age.children.dyad), 
                                (10-mean(t.data_disagreement$exact.age.children.dyad))/
                                 sd(t.data_disagreement$exact.age.children.dyad)),
                     labels = c(4, 5, 6, 7, 8, 9, 10)) +
 scale_y_continuous(name = "Proportion of reasons given", 
                    labels = scales::percent) +
 
 labs(title="Interaction effect of culture, trial number, and age") +
 theme_light() + 
 theme(panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title = element_text(size = 15),
       axis.text = element_text(size = 13),
       strip.text.x = element_text(size = 13),
       plot.margin=unit(c(1,1,1,1), "cm"),
       axis.title.y.left = element_text(vjust=3),
       plot.title = element_text(color="black", size=15, face="bold"), 
       axis.title.x = element_text(margin = 
                                    margin(t = 10, r = 0, b = 0, l = 0)))


