library(tidyverse)
library(parallel)
source("estimation_functions_binary.R")

dat_dir <- "real_dat"
dat3 <-
  read.table(paste(dat_dir, "/data3.csv", sep = ""),
             sep = ",",
             header = T)
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
  ) %>% mutate(admittime = as.Date(admittime, "%m/%d/%Y")) %>%
  mutate(Y = ifelse(sofa >= 3, 1, 0), .before = diasbp_min) %>%
  filter(age < 65) %>% group_by(subject_id) %>% slice(which.min(admittime)) %>%
  na.omit() %>% mutate(R = if_else(insurance %in% c("Medicare", "Medicaid"), 0, 1))
data %>% ungroup(subject_id) %>% select(-age,-subject_id,-admittime,-insurance,-sofa) -> data

# winequality <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv", sep = ";")
# data <- winequality %>%
#   mutate(Y = ifelse(quality <= 4, 1, 0),
#          .before = fixed.acidity) %>% 
#   mutate(R = ifelse(alcohol < 10.5, 1, 0)) %>% 
#   select(-quality)

sDat <- data %>% filter(R == 1) %>% select(-R) %>% as.data.frame()
tDat <-
  data %>% filter(R == 0) %>% select(-R,-Y) %>% as.data.frame()

logitReg <- glm(Y ~ ., data = sDat, family = binomial)
rfProbs <- predict(logitReg, newdata = sDat, type = "response")
predClass <- rfProbs > 0.5
table(sDat$Y, predClass)
rfProbt <- predict(logitReg, newdata = tDat, type = "response")

piVal <- nrow(sDat) / nrow(data)
trueBetaVal <- glm(R ~ Y, data = data, family = binomial)

initBeta <- 0
optim(
  par = initBeta,
  fn = Compute_S_Eff_Sum,
  method = "BFGS",
  sData = sDat,
  piVal = piVal,
  yFittedGivenX = rfProbs,
  yFittedGivenX_tDat = rfProbt
) -> optimOut # Beta proposed
c_ps <- E_S_RHO_Binary(optimOut$par, sDat)
Compute_SE_Binary(
  betaVal = optimOut$par,
  sData = sDat,
  c_ps = c_ps,
  piVal = piVal,
  yFittedGivenX = rfProbs,
  yFittedGivenX_tDat = rfProbt
) -> se # Variance Beta proposed

Estimate_Beta_Lipton(rfProbs, rfProbt, sDat, initBeta) -> naiveBeta # Beta Naive

initTheta <- 0.2
optim(
  par = initTheta,
  fn = Estimate_Theta_Sum_Eff,
  method = "Brent",
  betaVal = optimOut$par,
  yFittedGivenX = rfProbs,
  yFittedGivenX_t = rfProbt,
  piVal = piVal,
  sData = sDat,
  lower = 0,
  upper = 0.9
) -> thetaOut
thetaOut$par # Theta proposed

tDat_full <- data %>% filter(R == 0)
mean(tDat_full$Y) # Theta oracle

optim(
  par = initTheta,
  fn = Estimate_Theta_Sum_Naive,
  method = "Brent",
  betaVal = naiveBeta,
  sData = sDat,
  piVal = piVal,
  lower = 0,
  upper = 1
) -> thetaNaiveOut
thetaNaiveOut$par # Theta naive

B <- 200
rexpVec_List <- lapply(1:B, function(x) {
  rexp(nrow(data))
})
mclapply(rexpVec_List, function(rexpVec)
{
  Estimate_Beta_Lipton_Pert(rfProbs, rfProbt, sDat, initBeta, rexpVec) -> o1
  Estimate_optim_Theta_Naive_Pert(naiveBeta, sDat, piVal, rexpVec) -> o2
  Estimate_Theta_Pert_Optim(optimOut$par, rfProbs, rfProbt, piVal, sDat, rexpVec) -> o3
  return(list(
    beta = o1,
    theta = o2,
    theta_p = o3
  ))
},
mc.cores = detectCores()) -> beta_vector

unlist(beta_vector) %>% matrix(ncol = 3, byrow = T) -> pert_SE
apply(pert_SE, 2, sd)

