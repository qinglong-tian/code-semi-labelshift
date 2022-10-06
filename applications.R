library(tidyverse)
library(reshape2)
library(nnet)

source("application_func.R")
dat_dir <- "real_dat"
dat3 <- readr::read_csv(paste(dat_dir, "/data3.csv", sep = ""))

data <-
  dat3 %>% select(
    diasbp_min,
    glucose_min,
    resprate_max,
    sysbp_min,
    temp_min,
    temp_max,
    hematocrit_mean,
    hematocrit_min,
    platelets_mean,
    platelets_min,
    redbloodcell_mean,
    redbloodcell_min,
    spo2_min,
    urea_n_min,
    urea_n_max,
    urea_n_mean,
    insurance,
    age,
    subject_id,
    admittime,
    sofa
  ) %>% mutate(admittime = as.Date(admittime, "%m/%d/%Y")) %>% filter(age < 65,
                                                                      sofa <= 10) %>% group_by(subject_id) %>% slice(which.min(admittime)) %>% na.omit %>% mutate(R = if_else(insurance %in% c("Medicare", "Medicaid"), 0, 1))

rSofaLogistic <-
  glm(R ~ poly(sofa, 2),
      family = binomial(link = "logit"),
      data = data)

sDat <- data %>% filter(R == 1) %>% select(-R,-insurance)
tDat <- data %>% filter(R == 0) %>% select(-R,-insurance)

Compute_GLM_Pr_R1_Given_Y(sDat, tDat, rSofaLogistic) -> df
df <- df %>% mutate(SOFA = as.numeric(SOFA))
df %>% melt(
  measure.vars = c("rProb", "rProbPred"),
  variable.name = "Type",
  value.name = "Prob"
) -> df
df %>% ggplot(aes(x = SOFA, y = Prob, col = Type)) + geom_point() + geom_line()

#
# sDat %>% mutate(sofa = factor(sofa, levels = paste(0:10))) %>% mutate(sofa = relevel(sofa, ref = "0")) -> sDat1
# multinom_model <- multinom(sofa ~ ., data = sDat1)
# predict(multinom_model, newdata = sDat[1:5, ], type = "probs") -> probMat
# apply(probMat, 1, cumsum)

Pred_Sofa_Given_X_S(sDat, tDat) -> probList
probListXs <- probList$prob_list_xs
probListXt <- probList$prob_list_xt
piVal <- nrow(sDat) / (nrow(sDat) + nrow(tDat))

optim(
  0,
  Compute_S_Eff_Sum_App,
  sDat = sDat,
  piVal = piVal,
  prob_list_xs = probListXs,
  prob_list_xt = probListXt,
  method = "Brent",
  lower = -5,
  upper = 5
)
