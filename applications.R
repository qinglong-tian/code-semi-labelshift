library(tidyverse)
library(reshape2)
library(nnet)
library(fastGHQuad)

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
# winequality <-
#   read.csv(
#     "http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv",
#     sep = ";"
#   )
# sofa <- winequality$quality
# rd <- runif(nrow(winequality))
# winequality <- cbind(winequality[, -ncol(winequality)], sofa, rd)
# data <-
#   winequality %>% mutate(R = ifelse(sofa <= 5, rd < 0.15, rd < 0.75))
# data <- data[, -ncol(data)+1]

rSofaLogistic <-
  glm(R ~ sofa,
      family = binomial(link = "logit"),
      data = data)
data <- data %>% ungroup(subject_id)

sDat <-
  data %>% filter(R == 1) %>% select(-R, -insurance, -admittime,-subject_id, -age)
tDat <-
  data %>% filter(R == 0) %>% select(-R, -insurance, -admittime,-subject_id, -age)

# sDat <-
#   data %>% filter(R == 1) %>% select(-R)
# tDat <-
#   data %>% filter(R == 0) %>% select(-R)

Compute_GLM_Pr_R1_Given_Y(sDat, tDat, rSofaLogistic) -> df
df <- df %>% mutate(SOFA = as.numeric(SOFA))
df %>% melt(
  measure.vars = c("rProb", "rProbPred"),
  variable.name = "Type",
  value.name = "Prob"
) -> df
df %>% ggplot(aes(x = SOFA, y = Prob, col = Type)) + geom_point() + geom_line()

Fit_Sofa_Score(sDat) -> pyxs
piVal <- nrow(sDat) / (nrow(sDat) + nrow(tDat))
ghDat <- gaussHermiteData(50)

# # #
# Y <- sDat$sofa
# sDat1 <- sDat[,-ncol(sDat)]
# sDat1 <- cbind(Y, sDat1)
# tDat1 <- tDat[,-ncol(tDat)]
# 
# coef_y_x_s_hat <- coef(pyxs)
# sigma_y_x_s_hat <- sigma(pyxs)
# Mu_Y_S_hat <- mean(sDat1[, "Y"])
# Sig_Y_S_hat <- sd(sDat1[, "Y"])
# 
# GH_Materials <-
#   list(
#     y_vec = sDat1[, "Y"],
#     mu = Mu_Y_S_hat,
#     sigma = Sig_Y_S_hat,
#     xList = ghDat$x,
#     wList = ghDat$w
#   )
# #
# # Estimate_Beta_App_Wrapper <-
# #   function(betaVal,
# #            sData,
# #            tData,
# #            piVal,
# #            tDat_ext,
# #            coef_y_x_s,
# #            sigma_y_x_s,
# #            ispar,
# #            parameters,
# #            xList,
# #            wList,
# #            weights = FALSE)
# #   {
# #     beta_rho <- c(betaVal, 0)
# #     EstimateBetaFunc_CPP(
# #       beta_rho,
# #       sData,
# #       tData,
# #       piVal,
# #       tDat_ext,
# #       coef_y_x_s,
# #       sigma_y_x_s,
# #       ispar,
# #       parameters,
# #       xList,
# #       wList,
# #       weights
# #     )
# #   }
# #
# #
# optim(
#   c(0,0),
#   EstimateBetaFunc_CPP,
#   sData = as.matrix(sDat1),
#   tData = as.matrix(tDat1),
#   piVal = piVal,
#   tDat_ext = as.matrix(tDat1),
#   coef_y_x_s = coef_y_x_s_hat,
#   sigma_y_x_s = sigma_y_x_s_hat,
#   ispar = F,
#   parameters = GH_Materials,
#   xList = ghDat$x,
#   wList = ghDat$w,
#   weights = F
# ) -> fit_old
#
# EstimateBetaCovMat_CPP(
#   fit_old$par,
#   sData = as.matrix(sDat1),
#   tData = as.matrix(tDat1),
#   piVal = piVal,
#   tDat_ext = as.matrix(tDat1),
#   coef_y_x_s = coef_y_x_s_hat,
#   sigma_y_x_s = sigma_y_x_s_hat,
#   ispar = F,
#   parameters = GH_Materials,
#   xList = ghDat$x,
#   wList = ghDat$w
# )
#
# #

optim(
  -0.1,
  Compute_S_Eff_Sum_App,
  pyxs = pyxs,
  ghDat = ghDat,
  xMat_s = sDat[,-ncol(sDat)],
  xMat_t = tDat[,-ncol(tDat)],
  sDat = sDat,
  piVal = piVal,
  method = "Brent",
  lower = -5,
  upper = 5
) -> fit1
fit1
c_ps <- E_S_RHO_App(fit1$par, sDat)
optim(4, Compute_Eff_Theta_Sum_App,
      betaVal = fit1$par,
      pyxs = pyxs,
      ghDat = ghDat,
      xMat_t = tDat[,-ncol(tDat)],
      c_ps = c_ps,
      piVal = piVal,
      xMat_s = sDat[,-ncol(sDat)],
      sDat = sDat,
      method = "Brent",
      lower = 0,
      upper = 10
      )

optim(
  -0.1,
  Compute_S_Eff_Sum_Naive_App,
  pyxs = pyxs,
  ghDat = ghDat,
  xMat_s = sDat[,-ncol(sDat)],
  xMat_t = tDat[,-ncol(tDat)],
  sDat = sDat,
  piVal = piVal,
  method = "Brent",
  lower = -2,
  upper = 2
) -> fit2
fit2
c_ps <- E_S_RHO_App(fit2$par, sDat)
optim(
  4,
  Compute_Eff_Theta_Naive_Sum_App,
  betaVal = fit2$par,
  c_ps = c_ps,
  piVal = piVal,
  sDat = sDat,
  method = "Brent",
  lower = 0,
  upper = 10
)

rSofaLogistic

Compute_Beta_Var_App(
  fit1$par,
  pyxs,
  ghDat,
  as.matrix(sDat[, -ncol(sDat)]),
  as.matrix(tDat[, -ncol(tDat)]),
  E_S_RHO_App(fit1$par, sDat),
  piVal,
  sDat
) -> var1
fit1$par - 1.96 * sqrt(var1)

Compute_Beta_Var_App(
  fit2$par,
  pyxs,
  ghDat,
  as.matrix(sDat[, -ncol(sDat)]),
  as.matrix(tDat[, -ncol(tDat)]),
  E_S_RHO_App(fit2$par, sDat),
  piVal,
  sDat
) -> var2

fit2$par - 1.96 * sqrt(var2)
