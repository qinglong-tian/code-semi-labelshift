# Prior info and useful snippets
###########################################
library(tidyverse)
Mu_Y_T <- 1.5

read_in_data <- function(filename, dir)
{
  rt <- str_match(filename, "ratio(.*?)_n")[2] %>% as.numeric()
  n <- str_match(filename, "_n(.*?)_.RDS")[2] %>% as.numeric()
  dat <- readRDS(paste(dir, filename, sep = ""))
  
  return(list(rt = rt, dat = dat, n = n))
}

remove_outliers_2 <- function(vec)
{
  vec[!is.na(vec) & !is.finite(vec)] -> vec
  one_iqr <- IQR(vec, na.rm = T)
  lwrq <- quantile(vec, 0.25, na.rm = T) - 2 * one_iqr
  uprq <- quantile(vec, 0.75, na.rm = T) + 2 * one_iqr
  
  removed <- is.nan(vec) | is.na(vec) | (vec < lwrq) | (vec > uprq)
  out <- vec[!removed]
  
  return(list(cleaned = out, removed = removed))
}

compute_mse <- function(vec, trueval)
{
  mean((trueval-vec)^2)
}

dat_dir <- "dat2/"
filenames <- list.files(dat_dir)
###########################################
out <- NULL
for (filename in filenames)
{
  Dat <- read_in_data(filename, dat_dir)
  
  rt <- Dat$rt
  n <- Dat$n
  rtDat <- Dat$dat
  
  richInfo <- rtDat$Rich
  
  # Effcient Theta
  sapply(richInfo, function(x) {
    x$ThetaHat
  }) -> thetaHatVec
  eff_outlier <- remove_outliers_2(thetaHatVec)
  thetaHatVec <- eff_outlier$cleaned
  removed_eff <- eff_outlier$removed
  sdTheta <- sd(thetaHatVec)
  meanTheta <- mean(thetaHatVec)
  compute_mse(thetaHatVec, Mu_Y_T) -> mse_eff
  
  # Naive Theta
  sapply(richInfo, function(x) {
    x$ThetaNaive
  }) -> thetaHatNaiveVec
  naive_outlier <- remove_outliers_2(thetaHatNaiveVec)
  thetaHatNaiveVec <- naive_outlier$cleaned
  removed_naive <- naive_outlier$removed
  sdNaiveTheta <- sd(thetaHatNaiveVec)
  meanNaiveTheta <- mean(thetaHatNaiveVec)
  compute_mse(thetaHatNaiveVec, Mu_Y_T) -> mse_naive
  
  # Use formula: SE
  sapply(richInfo, function(x) {
    x$SdHat
  }) -> sdHatVec
  sdHatVec <- sdHatVec[!removed_eff]
  
  # Use formula: Coverage
  sapply(richInfo, function(x) {
    x$CP
  }) -> CPVec
  CPVec <- CPVec[!removed_eff]
  
  # Use Perturbation: SE
  sdPertVec <- rtDat$Sd
  sdPertVec <- sdPertVec[!removed_eff]
  remove_outliers_2(sdPertVec) -> pert_outliers
  removed_pert <- pert_outliers$removed
  sdPertVec <- pert_outliers$cleaned
  
  # Bias of the efficient estimator
  biasEff <- meanTheta - Mu_Y_T
  # SE of the efficient estimator
  seEff <- sdTheta
  
  # Bias of the naive estimator
  biasNaive <- meanNaiveTheta - Mu_Y_T
  # SE of the naive estimator
  seNaive <- sdNaiveTheta
  
  # Coverage of the Semi-inference
  cpFmla <- mean(CPVec)
  # Mean Sd of Semi-inference
  sdMeanSemi <- median(sdHatVec)
  
  # Coverage using perturbation
  lwbPert <- thetaHatVec[!removed_pert] - 1.96 * sdPertVec
  upbPert <- thetaHatVec[!removed_pert] + 1.96 * sdPertVec
  cpPert <- mean(lwbPert < Mu_Y_T & upbPert > Mu_Y_T)
  # Mean Sd of perturbation
  sdMeanPert <- median(sdPertVec)
  
  out <-
    rbind(
      out,
      c(
        rt,
        n,
        biasEff,
        seEff,
        biasNaive,
        seNaive,
        cpFmla,
        sdMeanSemi,
        cpPert,
        sdMeanPert,
        mse_eff,
        mse_naive
      )
    )
}
colnames(out) <-
  c(
    "ratio",
    "n",
    "BiasEff",
    "SEEff",
    "BiasNaive",
    "SENaive",
    "CPFormula",
    "SDFormula",
    "CPPert",
    "SDPert",
    "MSEEff",
    "MSENaive"
  )

out <- as.data.frame(out)
