# Study Name: How Disagreement influences Reason-Giving
# Authors: Hanna Schleihauf & Antonia Langenhoff
# Analysis: Hypotheses 7


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

xdata_disagreement <- xdata%>%
  filter(condition == "disagreement")

table(xdata_disagreement$decision.strategy)
table(xdata_disagreement$decision.strategy, xdata_disagreement$age.group)
table(xdata_disagreement$decision.strategy, xdata_disagreement$culture)
