################################################
# Simulation results are saved in dat1/
################################################
library(parallel)
library(Rcpp)
library(tidyverse)
library(fastGHQuad)
source("data_generating_functions.R")
sourceCpp("fast_estimation_functions.cpp")
################################################
# Data Generating Parameters
Mu_YX <- c(2, 1, 1, 1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1,-1] <- ar1_cor(3, 0.3)
SigMat_YX <- SigMat_YX + 0.1 * diag(4)
SigMat_YX[1, 1] <- 1.44
Mu_Y_T <- 1.5
Sig_Y_T <- 1.5
Mu_Y_S <- Mu_YX[1]
Sig_Y_S <- sqrt(SigMat_YX[1, 1])

trueBetaRho <-
  Compute_Rho_Parameters(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)
yx_dist <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
coef_y_x_s_true <- c(yx_dist$Beta0, yx_dist$Beta1)
var_y_x_s_true <- yx_dist$VarYX
sigma_y_x_s_true <- sqrt(var_y_x_s_true)

################################################
beta_rho <- trueBetaRho
B1 <- 100 # Monte-Carlo Sample Size
B2 <- 500

gh_num <- 12
ghxw <- gaussHermiteData(gh_num)
xList <- ghxw$x
wList <- ghxw$w
################################################
set.seed(888)
n_vector <- c(1000)
mnratio_vec <- c(1)

t0 <- Sys.time()
################################################
for (n in n_vector)
{
  for (mnratio in mnratio_vec)
  {
    m <- mnratio * n
    rexpVec_list <- mclapply(1:B2, function(x)
    {
      rexp(n + m)
    },
    mc.cores = detectCores())
    
    data_list_mc <- mclapply(1:B1, function(x)
    {
      Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
    },
    mc.cores = detectCores())
    
    results_output <- mclapply(data_list_mc, function(dataList)
    {
      sData <- dataList$sDat
      tData <- dataList$tDat
      piVal <- n / (n + m)
      ispar <- T
      
      fityx <- lm(Y ~ ., data = as.data.frame(sData))
      coef_y_x_s_hat <- coef(fityx)
      sigma_y_x_s_hat <- sigma(fityx)
      
      Mu_Y_S_hat <- mean(sData[, "Y"])
      Sig_Y_S_hat <- sd(sData[, "Y"])
      
      parameters <-
        list(
          y_vec = sData[, "Y"],
          mu = Mu_Y_S_hat,
          sigma = Sig_Y_S_hat,
          xList = xList,
          wList = wList
        )
      
      betaHat <-
        optim(
          beta_rho,
          EstimateBetaFunc_CPP,
          sData = sData,
          tData = tData,
          piVal = piVal,
          tDat_ext = tData,
          coef_y_x_s = coef_y_x_s_hat,
          sigma_y_x_s = sigma_y_x_s_hat,
          ispar = ispar,
          parameters = parameters,
          xList = xList,
          wList = wList
        )$par
      
      betaSd1 <-
        EstimateBetaVarFunc_CPP(
          betaHat,
          sData,
          tData,
          piVal,
          tData,
          coef_y_x_s_hat,
          sigma_y_x_s_hat,
          ispar,
          parameters,
          xList,
          wList
        )
      
      CI1 <- matrix(betaHat,
                    nrow = 2,
                    ncol = 2,
                    byrow = T)
      Sd1 <- matrix(betaSd1,
                    nrow = 2,
                    ncol = 2,
                    byrow = T)
      Sd1[1,] <- -1.96 * Sd1[1,]
      Sd1[2,] <- 1.96 * Sd1[2,]
      CI1 <- CI1 + Sd1
      
      CP1 <-
        c(
          trueBetaRho[1] >= CI1[1, 1] &
            trueBetaRho[1] <= CI1[2, 1],
          trueBetaRho[2] >= CI1[1, 2] &
            trueBetaRho[2] <= CI1[2, 2]
        )
      
      betaHatNaive <-
        optim(
          beta_rho,
          Estimate_Naive_Beta,
          sData = sData,
          tData = tData,
          piVal = piVal,
          ispar = ispar,
          parameters = parameters
        )$par
      
      # Compute the variance of the naive estimator
      DerivMat <-
        Compute_Derivative_Mat(betaHatNaive, sData, tData, ispar, parameters)
      EMM <-
        Compute_Var_m(betaHatNaive, sData, tData, ispar, parameters, piVal)
      betaSd2 <-
        sqrt(diag(solve(DerivMat) %*% EMM %*% t(solve(DerivMat)) / (n + m)))
      CI2 <- matrix(betaHatNaive,
                    nrow = 2,
                    ncol = 2,
                    byrow = T)
      Sd2 <- matrix(betaSd2,
                    nrow = 2,
                    ncol = 2,
                    byrow = T)
      
      Sd2[1,] <- -1.96 * Sd2[1,]
      Sd2[2,] <- +1.96 * Sd2[2,]
      CI2 <- CI2 + Sd2
      
      CP2 <-
        c(
          trueBetaRho[1] >= CI2[1, 1] &
            trueBetaRho[1] <= CI2[2, 1],
          trueBetaRho[2] >= CI2[1, 2] &
            trueBetaRho[2] <= CI2[2, 2]
        )
      
      return(
        list(
          TrueBeta = trueBetaRho,
          BetaHatEff = betaHat,
          BetaHatNaive = betaHatNaive,
          SdEff = betaSd1,
          SdNaive = betaSd2,
          CPEff = CP1,
          CPNaive = CP2
        )
      )
    },
    mc.cores = detectCores())
    
    results_output_pert <-
      mclapply(data_list_mc, function(dataList) {
        sData <- dataList$sDat
        tData <- dataList$tDat
        piVal <- n / (n + m)
        ispar <- T
        
        fityx <- lm(Y ~ ., data = as.data.frame(sData))
        coef_y_x_s_hat <- coef(fityx)
        sigma_y_x_s_hat <- sigma(fityx)
        
        Mu_Y_S_hat <- mean(sData[, "Y"])
        Sig_Y_S_hat <- sd(sData[, "Y"])
        
        parameters <-
          list(
            y_vec = sData[, "Y"],
            mu = Mu_Y_S_hat,
            sigma = Sig_Y_S_hat,
            xList = xList,
            wList = wList
          )
        
        betaVec <- matrix(nrow = B2, ncol = 2)
        for (j in 1:B2)
        {
          rexpVec <- rexpVec_list[[j]]
          betaHatPert <-
            optim(
              beta_rho,
              EstimateBetaFuncPert_CPP,
              sData = sData,
              tData = tData,
              piVal = piVal,
              tDat_ext = tData,
              coef_y_x_s = coef_y_x_s_hat,
              sigma_y_x_s = sigma_y_x_s_hat,
              ispar = ispar,
              parameters = parameters,
              xList = xList,
              wList = wList,
              rexpVec = rexpVec
            )
          betaVec[j,] <- betaHatPert$par
        }
        
        return(betaVec)
      },
      mc.cores = detectCores())
    
    outList <-
      list(Out1 = results_output, Out2 = results_output_pert)
    saveRDS(outList,
            file = paste("dat3/n_", n, "_ratio_", mnratio, "_.RDS", sep = ""))
  }
}

Sys.time() - t0
