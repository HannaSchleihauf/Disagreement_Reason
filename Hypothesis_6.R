# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 6A: We expect that children will engage in more meta-talk
# (reasoning about reasons) in the disagreement than in the agreement condition.


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

table(xdata$metatalk.dyad.explicit)


#We do not have enough variance to run this analysis.
