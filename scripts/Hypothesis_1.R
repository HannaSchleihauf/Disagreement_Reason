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
xdata <-  read.csv(file="./data/Disagreement_Reason_Coding_Sheet_All_Three.csv",
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

full <- glmer(resp_mat ~ condition*culture*z.age + test.trial +
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data, control = contr,
              family = binomial(link = "logit"))
red <- glmer(resp_mat ~ (condition + culture + z.age)^2 + test.trial +
              (1 + test.trial.2.code || dyad.nr) + (1 | olre),
             data = t.data, control = contr,
             family = binomial(link = "logit"))
red2 <- glmer(resp_mat ~ (condition + culture)^2 + z.age + test.trial +
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data, control = contr,
              family = binomial(link = "logit"))
main <- glmer(resp_mat ~ (condition + culture) + z.age + test.trial +
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data, control = contr,
              family = binomial(link = "logit"))
null <- glmer(resp_mat ~ 1 +
               (1 + test.trial.2.code || dyad.nr) + (1 | olre),
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
round(drop1(red2, test="Chisq"), 3)
round(drop1(main, test="Chisq"), 3)

plot(effect("condition:culture", red2), type = "response")


# coefficients of the full model
summary(full)$coefficients

emm1 <- emmeans(red2, ~ condition * culture)
summary(emm1, type = "response") # no effect of condition
summary(pairs(regrid(emm1)), type = "response")

## Plotting Preparations
# Long data set for plotting
ydata <-
 xdata %>%
 pivot_longer(
  cols = c(perceptual.reason.given.dyad.nr, testimonial.reason.given.dyad.nr),
  names_to = "type.reason",
  values_to = "reason.given"
 )

## Plot Culture Condition
# Boot straps for plotting
interaction.plot.model <-
 glmer(resp_mat ~ (condition + culture)^2  + z.age + test.trial.2.code +
        (1 + test.trial.2.code || dyad.nr) + (1 | olre),
       data = t.data, control = contr, family = binomial(link = "logit"))

boot_res_condition_culture <-
 boot.glmm.pred(
  model.res = interaction.plot.model, excl.warnings = TRUE,
  nboots = 1000, para = FALSE, level = 0.95,
  use = c("condition", "culture")
 )

# Aggregated data for plotting
ydata.agg <- ydata %>%
 group_by(dyad.nr, condition, culture) %>%
 summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
 ungroup()

exp1_plot_condition.culture <-
 ggplot(data = ydata.agg,
        aes(x = factor(condition,
            levels=c("agreement", "disagreement")),
            y = mean.resp)) +

 geom_point(data = ydata.agg %>% filter(condition == "agreement"),
            aes(x = condition), color = "dodgerblue", size = 1.5,
            alpha = .4, position = position_jitter(w = 0.2, h = 0.015)) +
 geom_point(data = ydata.agg %>% filter(condition == "disagreement"),
            aes(x = condition), color = "darkorange", size = 1.5,
            alpha = .4, position = position_jitter(w = 0.2, h = 0.015)) +

 geom_violin(data = ydata.agg %>% filter(condition == "agreement"),
             aes(x = condition, y = mean.resp),
             position = position_nudge(x = 0.00),
             fill = "dodgerblue", alpha = .2) +
 geom_violin(data = ydata.agg %>% filter(condition == "disagreement"),
             aes(x = condition, y = mean.resp),
             position = position_nudge(x = 0.00),
             fill = "darkorange", alpha = .2) +

 geom_errorbar(data = boot_res_condition_culture$ci.predicted %>%
                filter(condition == "agreement"),
               aes(x = as.numeric(condition) + 0.25, y = fitted,
                   ymin = lower.cl, ymax = upper.cl), color = "dodgerblue",
               width = 0.1, size = 1) +
 geom_errorbar(data = boot_res_condition_culture$ci.predicted %>%
                filter(condition == "disagreement"),
               aes(x = as.numeric(condition) + 0.25, y = fitted,
                   ymin = lower.cl, ymax = upper.cl), color = "darkorange",
               width = 0.1, size = 1) +

 geom_point(data = boot_res_condition_culture$ci.predicted %>%
             filter(condition == "agreement"),
            aes(x = as.numeric(condition) + 0.25, y = fitted),
            color = "dodgerblue", size = 2.5) +
 geom_point(data = boot_res_condition_culture$ci.predicted %>%
             filter(condition == "disagreement"),
            aes(x = as.numeric(condition) + 0.25, y = fitted),
            color = "darkorange", size = 2.5) +

 facet_wrap(~culture, switch = "x") +

 scale_x_discrete(limits = c("agreement", "disagreement"),
                  name = "Condition",
                  labels = c("Agreement\ncondition", "Disagreement\ncondition")) +
 scale_y_continuous(name = "Proportion of reasons given",
                    labels = scales::percent) +

 labs(title="Interaction effect of condition and culture") +
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
exp1_plot_condition.culture


# Age Plot
#Bootstrap
age.model <-
 glmer(resp_mat ~ (condition.disagreement.code +
                  culture.USA.code)^2  +
               z.age + test.trial.2.code +
        (1 + test.trial.2.code || dyad.nr) + (1 | olre),
              data = t.data, control = contr, family = binomial(link = "logit"))
boot_res_age <-
 boot.glmm.pred(
  model.res = age.model, excl.warnings = TRUE,
  nboots = 1000, para = FALSE, level = 0.95,
  use = c("z.age")
 )

# Aggregate Data
ydata$z.age <- scale(ydata$exact.age.children.dyad)
ydata.agg.age <- ydata %>%
 group_by(dyad.nr, condition, z.age) %>%
 summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
 ungroup()

# Plot Age
exp1_plot_age <-
 ggplot() +
 geom_point(data = ydata.agg.age,
            aes(x = z.age, y = mean.resp, color = condition),
            size = 1.5, alpha = .4,
            position = position_jitter(w = 0.00, h = 0.015)) +
 scale_color_manual(values=c("dodgerblue", "darkorange"),
                    labels = c("Agreement condition",
                               "Disagreement condition       "),
                    name = NULL) +

 geom_ribbon(data = boot_res_age$ci.predicted,
             aes(ymin = lower.cl, ymax = upper.cl, x = z.age),
             alpha = 0.3, fill = "grey") +

 geom_line(data=boot_res_age$ci.predicted, aes(x=z.age, y=fitted),
           color = "darkgrey") +

 scale_x_continuous(name="Age",
                    breaks = c((5-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (6-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (7-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (8-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (9-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (10-mean(ydata$exact.age.children.dyad))/
                               sd(ydata$exact.age.children.dyad)),
                    labels = c(5, 6, 7, 8, 9, 10)) +
 scale_y_continuous(name = "Proportion of reasons given",
                    labels = scales::percent) +
 labs(title="Effect of age") +
 theme_classic() +
 theme(legend.position="top",
       legend.direction = "horizontal",
       legend.background = element_blank(),
       legend.box.background = element_rect(colour = "lightgrey"),
       axis.title = element_text(size = 15),
       axis.text = element_text(size = 13),
       plot.margin=unit(c(1,1,1,1), "cm"),
       axis.title.y.left = element_text(vjust=3),
       plot.title = element_text(color="black", size=15, face="bold"),
       axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
exp1_plot_age


# Trial
ydata$z.trial <- scale(as.numeric(ydata$test.trial))
ydata.agg.trial <- ydata %>%
 group_by(dyad.nr, condition, test.trial) %>%
 summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
 ungroup()

trial.model <-
 glmer(resp_mat ~ (condition.disagreement.code +
                    culture.USA.code)^2  +
        z.age + test.trial +
        (1 + test.trial.2.code || dyad.nr) + (1 | olre),
       data = t.data, control = contr, family = binomial(link = "logit"))

boot_res_trial <-
 boot.glmm.pred(
  model.res = trial.model, excl.warnings = TRUE,
  nboots = 1000, para = FALSE, level = 0.95,
  use = c("test.trial")
 )


exp1_plot_trial <-
 ggplot() +
 geom_point(data = ydata.agg.trial,
            aes(x = test.trial, y = mean.resp, color = condition),
            size = 1.5, alpha = .4,
            position = position_jitter(w = 0.2, h = 0.05)) +
 scale_color_manual(values=c("dodgerblue", "darkorange"),
                    labels = c("Agreement condition",
                               "Disagreement condition       "),
                    name = NULL) +

 geom_violin(data = ydata.agg.trial %>% filter(test.trial == "1"),
             aes(x = test.trial, y = mean.resp),
             position = position_nudge(x = 0.00),
             fill = "lightgrey", alpha = .2, linewidth = 0) +
 geom_violin(data = ydata.agg.trial %>% filter(test.trial == "2"),
             aes(x = test.trial, y = mean.resp),
             position = position_nudge(x = 0.00),
             fill = "lightgrey", alpha = .2, linewidth = 0) +

 geom_errorbar(data = boot_res_trial$ci.predicted %>%
                filter(test.trial == "1"),
               aes(x = as.numeric(test.trial) + 0.25, y = fitted,
                   ymin = lower.cl, ymax = upper.cl), color = "darkgrey",
               width = 0.1, size = 1) +
 geom_errorbar(data = boot_res_trial$ci.predicted %>%
                filter(test.trial == "2"),
               aes(x = as.numeric(test.trial) + 0.25, y = fitted,
                   ymin = lower.cl, ymax = upper.cl), color = "darkgrey",
               width = 0.1, size = 1) +

 geom_point(data =boot_res_trial$ci.predicted %>%
             filter(test.trial == "1"),
            aes(x = as.numeric(test.trial) + 0.25, y = fitted),
            color = "darkgrey", size = 2.5) +
 geom_point(data = boot_res_trial$ci.predicted %>%
             filter(test.trial == "2"),
            aes(x = as.numeric(test.trial) + 0.25, y = fitted),
            color = "darkgrey", size = 2.5) +

 scale_x_discrete(limits = c("1", "2"),
                  name = "Test Trial",
                  labels = c("Test trial 1", "Test trial 2")) +
 scale_y_continuous(name = "Proportion of reasons given",
                    labels = scales::percent) +

 labs(title="Effect of test trial") +
  theme_classic() +
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "lightgrey"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title.y.left = element_text(vjust=3),
        plot.title = element_text(color="black", size=15, face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
exp1_plot_trial



# Age * Culture Plot
#Bootstrap
age.culture.model <-
 glmer(resp_mat ~ condition.disagreement.code +
                    culture*z.age + test.trial.2.code +
        (1 + test.trial.2.code || dyad.nr) + (1 | olre),
       data = t.data, control = contr, family = binomial(link = "logit"))
boot_res_age_culture <-
 boot.glmm.pred(
  model.res = age.culture.model, excl.warnings = TRUE,
  nboots = 1000, para = FALSE, level = 0.95,
  use = c("z.age", "culture")
 )

# Aggregate Data
ydata$z.age <- scale(ydata$exact.age.children.dyad)
ydata.agg.age.culture <- ydata %>%
 group_by(dyad.nr, culture, z.age, condition) %>%
 summarise(mean.resp = mean(reason.given, na.rm = T)) %>%
 ungroup()

# Plot Age Culture
exp1_plot_age_culture <-
 ggplot() +
 geom_point(data = ydata.agg.age.culture,
            aes(x = z.age, y = mean.resp, color = condition),
            size = 1.5, alpha = .4,
            position = position_jitter(w = 0.00, h = 0.015)) +
 scale_color_manual(values=c("dodgerblue", "darkorange"),
                    labels = c("Agreement condition",
                               "Disagreement condition       "),
                    name = NULL) +

 geom_ribbon(data = boot_res_age_culture$ci.predicted,
             aes(ymin = lower.cl, ymax = upper.cl, x = z.age),
             alpha = 0.3, fill = "grey") +

 geom_line(data=boot_res_age_culture$ci.predicted,
           aes(x=z.age, y=fitted),
           color = "darkgrey") +

 facet_grid(~ culture) +


 scale_x_continuous(name="Age",
                    breaks = c((4-mean(t.data_disagreement$exact.age.children.dyad))/
                                sd(t.data_disagreement$exact.age.children.dyad),
                               (5-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (6-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (7-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (8-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (9-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad),
                               (10-mean(ydata$exact.age.children.dyad))/
                                sd(ydata$exact.age.children.dyad)),
                    labels = c(4, 5, 6, 7, 8, 9, 10)) +
 scale_y_continuous(name = "Proportion of reasons given",
                    labels = scales::percent) +
 labs(title="Interaction effect of age and culture") +
 theme_classic() +
 theme(legend.position="top",
       legend.direction = "horizontal",
       legend.background = element_blank(),
       legend.box.background = element_rect(colour = "lightgrey"),
       axis.title = element_text(size = 15),
       axis.text = element_text(size = 13),
       plot.margin=unit(c(1,1,1,1), "cm"),
       axis.title.y.left = element_text(vjust=3),
       plot.title = element_text(color="black", size=15, face="bold"),
       axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
exp1_plot_age_culture


save.image("First_Peek")





