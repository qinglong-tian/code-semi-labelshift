library(tidyverse)
library(reshape2)
library(fastGHQuad)
library(CondIndTests)
library(parallel)
source("application_func.R")

B <- 200
set.seed(37)
## Read in the data
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
data <-
  data %>% ungroup(subject_id) %>% select(-insurance, -admittime,-subject_id, -age)

# CondIndTest(as.matrix(data[, -c(ncol(data), ncol(data) - 1)]),
#             as.matrix(data[, ncol(data)]),
#             as.matrix(data[, ncol(data)-1]),
#             verbose = T)

sDat <-
  data %>% filter(R == 1) %>% select(-R)
tDat <-
  data %>% filter(R == 0) %>% select(-R)
wilcox.test(sDat$sofa, tDat$sofa, paired = F)

# Find true value
rSofaLogistic <-
  glm(R ~ sofa,
      family = binomial(link = "logit"),
      data = data)

# 
# Compute_GLM_Pr_R1_Given_Y(sDat, tDat, rSofaLogistic) -> df
# df <- df %>% mutate(SOFA = as.numeric(SOFA))
# df %>% melt(
#   measure.vars = c("rProb", "rProbPred"),
#   variable.name = "Type",
#   value.name = "Prob"
# ) -> df
# df %>% ggplot(aes(x = SOFA, y = Prob, col = Type)) + geom_point() + geom_line()

## Estimation
Fit_Sofa_Score(sDat) -> pyxs
piVal <- nrow(sDat) / (nrow(sDat) + nrow(tDat))
ghDat <- gaussHermiteData(12)

# Proposed Method
optim(
  -0.1,
  Compute_S_Eff_Sum_App,
  pyxs = pyxs,
  ghDat = ghDat,
  xMat_s = sDat[, -ncol(sDat)],
  xMat_t = tDat[, -ncol(tDat)],
  sDat = sDat,
  piVal = piVal,
  method = "Brent",
  lower = -5,
  upper = 5
) -> fit1

c_ps <- E_S_RHO_App(fit1$par, sDat)
Compute_Beta_Var_App(fit1$par,
                     pyxs,
                     ghDat,
                     as.matrix(sDat[, -ncol(sDat)]),
                     as.matrix(tDat[, -ncol(tDat)]),
                     c_ps,
                     piVal,
                     sDat) -> var1

optim(
  4,
  Compute_Eff_Theta_Sum_App,
  betaVal = fit1$par,
  pyxs = pyxs,
  ghDat = ghDat,
  xMat_t = tDat[, -ncol(tDat)],
  c_ps = c_ps,
  piVal = piVal,
  xMat_s = sDat[, -ncol(sDat)],
  sDat = sDat,
  method = "Brent",
  lower = 0,
  upper = 10
) -> proposed_theta_est

rexpVecList <- lapply(1:B, function(x) {
  rexp(nrow(data))
})

mclapply(rexpVecList, function(rexpVec)
{
  tryCatch({
    optim(
      -0.1,
      Compute_S_Eff_Pert_Sum_App,
      pyxs = pyxs,
      ghDat = ghDat,
      xMat_s = sDat[,-ncol(sDat)],
      xMat_t = tDat[,-ncol(tDat)],
      sDat = sDat,
      piVal = piVal,
      rexpVec = rexpVec,
      method = "Brent",
      lower = -5,
      upper = 5
    ) -> fit1
    para <- fit1$par
    c_ps <- E_S_RHO_App(para, sDat)
    optim(
      4,
      Compute_Eff_Theta_Sum_Pert_App,
      betaVal = para,
      pyxs = pyxs,
      ghDat = ghDat,
      xMat_t = tDat[,-ncol(tDat)],
      c_ps = c_ps,
      piVal = piVal,
      xMat_s = sDat[,-ncol(sDat)],
      sDat = sDat,
      rexpVec = rexpVec,
      method = "Brent",
      lower = 0,
      upper = 18
    ) -> fit_theta
    theta_b <- fit_theta$par
  },
  error = function(e)
  {
    return(NULL)
  },
  warning = function(w)
  {
    return(NULL)
  })
  
  return(theta_b)
},
mc.cores = detectCores()) -> proposed_theta_B


# Naive Method
optim(
  -0.1,
  Compute_S_Eff_Sum_Naive_App,
  pyxs = pyxs,
  ghDat = ghDat,
  xMat_s = sDat[, -ncol(sDat)],
  xMat_t = tDat[, -ncol(tDat)],
  sDat = sDat,
  piVal = piVal,
  method = "Brent",
  lower = -5,
  upper = 5
) -> fit2
c_ps <- E_S_RHO_App(fit2$par, sDat)
optim(
  4,
  Compute_Eff_Theta_Naive_Sum_App,
  betaVal = fit2$par,
  c_ps = c_ps,
  piVal = piVal,
  sDat = sDat,
  method = "Brent",
  lower = -3,
  upper = 20
) -> naive_theta_est

mclapply(rexpVecList, function(rexpVec)
{
  tryCatch({
    optim(
      -0.1,
      Compute_S_Eff_Sum_Naive_Pert_App,
      pyxs = pyxs,
      ghDat = ghDat,
      xMat_s = sDat[, -ncol(sDat)],
      xMat_t = tDat[, -ncol(tDat)],
      sDat = sDat,
      piVal = piVal,
      rexpVec = rexpVec,
      method = "Brent",
      lower = -2,
      upper = 2
    ) -> fit22
    c_ps <- E_S_RHO_App(fit22$par, sDat)
    optim(
      4,
      Compute_Eff_Theta_Naive_Sum_Pert_App,
      betaVal = fit22$par,
      c_ps = c_ps,
      piVal = piVal,
      sDat = sDat,
      rexpVec = rexpVec,
      method = "Brent",
      lower = 0,
      upper = 10
    ) -> naive_thetaa
  },
  error = function(e)
  {
    return(NULL)
  },
  warning = function(w)
  {
    return(NULL)
  })
  
  return(naive_thetaa$par)
},
mc.cores = detectCores()) -> naive_theta_B

Compute_Beta_Var_App(
  fit2$par,
  pyxs,
  ghDat,
  as.matrix(sDat[,-ncol(sDat)]),
  as.matrix(tDat[,-ncol(tDat)]),
  E_S_RHO_App(fit2$par, sDat),
  piVal,
  sDat
) -> var2

print_results(fit1, proposed_theta_est, var1, proposed_theta_B)

print_results(fit2, naive_theta_est, var2, naive_theta_B)

-rSofaLogistic$coefficients[2]
mean(tDat$sofa)
