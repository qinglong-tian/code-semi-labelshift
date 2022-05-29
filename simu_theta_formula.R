################################################
library(parallel)
library(Rcpp)
library(tidyverse)
library(fastGHQuad)
source("data_generating_functions.R")
sourceCpp("fast_estimation_functions.cpp")
################################################
# Data Generation #
Mu_YX <- c(2, 1, 1, 1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1, -1] <- ar1_cor(3, 0.3)
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
B1 <- 1000 # Monte-Carlo Sample Size
B2 <- 500
################################################
# Factors #
n_vec <- c(200, 400, 600, 800, 1000)
mnratio_vec <- c(0.5, 1, 1.5)


for (n in n_vec)
{
  for (mnratio in mnratio_vec)
  {
    m <- mnratio * n
    
    rexpVec_list <- mclapply(1:B2, function(x)
    {
      rexp(n + m)
    },
    mc.cores = detectCores())
    
    beta_rho <- trueBetaRho
    # B2 <- 10 # Bootstrap (Perturbation Size)
    
    gh_num <- 10
    ghxw <- gaussHermiteData(gh_num)
    xList <- ghxw$x
    wList <- ghxw$w
    ################################################
    set.seed(888)
    data_list_mc <- mclapply(1:B1, function(x) {
      Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
    }, mc.cores = detectCores())
    ################################################
    results_output <- mclapply(data_list_mc, function(dataList) {
      sData <- dataList$sDat
      tData <- dataList$tDat
      piVal <- n / (n + m)
      ispar <- T
      
      fityx <- lm(Y ~ ., data = as.data.frame(sData))
      coef_y_x_s_hat <- coef(fityx)
      sigma_y_x_s_hat <- sigma(fityx)
      
      Mu_Y_S_hat <- mean(sData[, "Y"])
      Sig_Y_S_hat <- sd(sData[, "Y"])
      GH_Materials <-
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
          parameters = GH_Materials,
          xList = xList,
          wList = wList
        )
      betaHat <- betaHat$par
      
      num_of_target <- m
      yVec <- sData[, "Y"]
      rhoValSource <-
        exp(c(cbind(yVec, yVec ^ 2) %*% matrix(betaHat, ncol = 1)))
      c_ps <- E_S_RHO_CPP(betaHat, ispar, GH_Materials)
      e_s_rho_x_ext <-
        E_S_RHO_X_CPP(betaHat,
                      1,
                      tData,
                      coef_y_x_s_hat,
                      sigma_y_x_s_hat,
                      xList,
                      wList)
      e_s_rho2_x_ext <-
        E_S_RHO_X_CPP(betaHat,
                      2,
                      tData,
                      coef_y_x_s_hat,
                      sigma_y_x_s_hat,
                      xList,
                      wList)
      e_t_tau <-
        E_T_TAU_CPP(
          e_s_rho_x = e_s_rho_x_ext,
          e_s_rho2_x = e_s_rho2_x_ext,
          c_ps = c_ps,
          piVal = piVal
        )
      tau_x_external <-
        COMPUTE_TAU_CPP(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal)
      
      xMatAll <- rbind(sData[, -1], tData)
      e_s_rho_x_all <-
        E_S_RHO_X_CPP(betaHat,
                      1,
                      xMatAll,
                      coef_y_x_s_hat,
                      sigma_y_x_s_hat,
                      xList,
                      wList)
      e_s_rho2_x_all <-
        E_S_RHO_X_CPP(betaHat,
                      2,
                      xMatAll,
                      coef_y_x_s_hat,
                      sigma_y_x_s_hat,
                      xList,
                      wList)
      tau_x_internal_all <-
        COMPUTE_TAU_CPP(e_s_rho_x_all, e_s_rho2_x_all, c_ps, piVal)
      e_s_rho2_psi_x_internal_all <-
        E_S_RHO2_PSI_X_CPP(betaHat,
                           xMatAll,
                           coef_y_x_s_hat,
                           sigma_y_x_s_hat,
                           xList,
                           wList)
      e_s_rho2_psi_x_external <-
        E_S_RHO2_PSI_X_CPP(betaHat,
                           tData,
                           coef_y_x_s_hat,
                           sigma_y_x_s_hat,
                           xList,
                           wList)
      MatInv <-
        EstimateBetaCovMat_CPP(
          betaHat,
          sData,
          tData,
          piVal,
          tData,
          coef_y_x_s_hat,
          sigma_y_x_s_hat,
          ispar,
          GH_Materials,
          xList,
          wList
        )
      e_s_rho2_psi_x_internal_source <-
        E_S_RHO2_PSI_X_CPP(betaHat,
                           sData[, -1],
                           coef_y_x_s_hat,
                           sigma_y_x_s_hat,
                           xList,
                           wList)
      e_s_rho2_x_internal_source <-
        E_S_RHO_X_CPP(betaHat,
                      2,
                      sData[, -1],
                      coef_y_x_s_hat,
                      sigma_y_x_s_hat,
                      xList,
                      wList)
      e_s_rho_x_internal_source <-
        E_S_RHO_X_CPP(betaHat,
                      1,
                      sData[, -1],
                      coef_y_x_s_hat,
                      sigma_y_x_s_hat,
                      xList,
                      wList)
      tau_for_x_source <-
        COMPUTE_TAU_CPP(
          e_s_rho_x = e_s_rho_x_internal_source,
          e_s_rho2_x = e_s_rho2_x_internal_source,
          c_ps = c_ps,
          piVal = piVal
        )
      
      SEff <-
        ComputeEfficientScore_CPP(
          betaHat,
          sData,
          tData,
          piVal,
          tData,
          coef_y_x_s_hat,
          sigma_y_x_s_hat,
          ispar,
          GH_Materials,
          xList,
          wList
        )
      
      thetaHat <-
        COMPUTE_THETA_CPP(
          num_of_target,
          piVal,
          rhoValSource,
          c_ps,
          yVec,
          betaHat,
          e_t_tau,
          tau_x_internal_all,
          tau_x_external,
          e_s_rho2_psi_x_internal_all,
          e_s_rho2_x_all,
          e_s_rho2_psi_x_external,
          e_s_rho2_x_ext
        )
      
      thetaNHat <-
        COMPUTE_THETA_NAIVE(betaHat, rhoValSource, yVec, c_ps)
      
      Phi_Theta <-
        COMPUTE_EFFICIENT_IF_FOR_THETA_CPP(
          thetaHat,
          num_of_target,
          piVal,
          rhoValSource,
          c_ps,
          yVec,
          betaHat,
          e_t_tau,
          tau_x_internal_all,
          tau_x_external,
          e_s_rho2_psi_x_internal_all,
          e_s_rho2_x_all,
          e_s_rho2_psi_x_external,
          e_s_rho2_x_ext,
          MatInv,
          ispar,
          GH_Materials,
          sData,
          coef_y_x_s_hat,
          sigma_y_x_s_hat,
          e_s_rho2_psi_x_internal_source,
          e_s_rho2_x_internal_source,
          tau_for_x_source,
          e_s_rho_x_internal_source,
          xList,
          wList,
          SEff
        )
      
      var_est <- mean(Phi_Theta ^ 2)
      sd_est <- sqrt(var_est / (n + m))
      CI_Lower <- thetaHat - 1.96 * sd_est
      CI_Upper <- thetaHat + 1.96 * sd_est
      CP <- (Mu_Y_T > CI_Lower) & (Mu_Y_T < CI_Upper)
      
      return(
        list(
          ThetaHat = thetaHat,
          SdHat = sd_est,
          Bias = thetaHat - Mu_Y_T,
          CP = CP,
          ThetaNaive = thetaNHat
        )
      )
    },
    mc.cores = detectCores())
    
    ############################
    # Perturbation
    ############################
    
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
        GH_Materials <-
          list(
            y_vec = sData[, "Y"],
            mu = Mu_Y_S_hat,
            sigma = Sig_Y_S_hat,
            xList = xList,
            wList = wList
          )
        
        thetaHatVec <- numeric(B2)
        for (i in 1:B2)
        {
          rexpVec <- rexpVec_list[[i]]
          betaHat <-
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
              parameters = GH_Materials,
              xList = xList,
              wList = wList,
              rexpVec = rexpVec
            )
          betaHat <- betaHat$par
          
          num_of_target <- m
          yVec <- sData[, "Y"]
          rhoValSource <-
            exp(c(cbind(yVec, yVec ^ 2) %*% matrix(betaHat, ncol = 1)))
          c_ps <- E_S_RHO_CPP(betaHat, ispar, GH_Materials)
          e_s_rho_x_ext <-
            E_S_RHO_X_CPP(betaHat,
                          1,
                          tData,
                          coef_y_x_s_hat,
                          sigma_y_x_s_hat,
                          xList,
                          wList)
          e_s_rho2_x_ext <-
            E_S_RHO_X_CPP(betaHat,
                          2,
                          tData,
                          coef_y_x_s_hat,
                          sigma_y_x_s_hat,
                          xList,
                          wList)
          e_t_tau <-
            E_T_TAU_CPP(
              e_s_rho_x = e_s_rho_x_ext,
              e_s_rho2_x = e_s_rho2_x_ext,
              c_ps = c_ps,
              piVal = piVal
            )
          tau_x_external <-
            COMPUTE_TAU_CPP(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal)
          
          xMatAll <- rbind(sData[, -1], tData)
          e_s_rho_x_all <-
            E_S_RHO_X_CPP(betaHat,
                          1,
                          xMatAll,
                          coef_y_x_s_hat,
                          sigma_y_x_s_hat,
                          xList,
                          wList)
          e_s_rho2_x_all <-
            E_S_RHO_X_CPP(betaHat,
                          2,
                          xMatAll,
                          coef_y_x_s_hat,
                          sigma_y_x_s_hat,
                          xList,
                          wList)
          tau_x_internal_all <-
            COMPUTE_TAU_CPP(e_s_rho_x_all, e_s_rho2_x_all, c_ps, piVal)
          e_s_rho2_psi_x_internal_all <-
            E_S_RHO2_PSI_X_CPP(betaHat,
                               xMatAll,
                               coef_y_x_s_hat,
                               sigma_y_x_s_hat,
                               xList,
                               wList)
          e_s_rho2_psi_x_external <-
            E_S_RHO2_PSI_X_CPP(betaHat,
                               tData,
                               coef_y_x_s_hat,
                               sigma_y_x_s_hat,
                               xList,
                               wList)
          MatInv <-
            EstimateBetaCovMat_CPP(
              betaHat,
              sData,
              tData,
              piVal,
              tData,
              coef_y_x_s_hat,
              sigma_y_x_s_hat,
              ispar,
              GH_Materials,
              xList,
              wList
            )
          e_s_rho2_psi_x_internal_source <-
            E_S_RHO2_PSI_X_CPP(betaHat,
                               sData[, -1],
                               coef_y_x_s_hat,
                               sigma_y_x_s_hat,
                               xList,
                               wList)
          e_s_rho2_x_internal_source <-
            E_S_RHO_X_CPP(betaHat,
                          2,
                          sData[, -1],
                          coef_y_x_s_hat,
                          sigma_y_x_s_hat,
                          xList,
                          wList)
          e_s_rho_x_internal_source <-
            E_S_RHO_X_CPP(betaHat,
                          1,
                          sData[, -1],
                          coef_y_x_s_hat,
                          sigma_y_x_s_hat,
                          xList,
                          wList)
          tau_for_x_source <-
            COMPUTE_TAU_CPP(
              e_s_rho_x = e_s_rho_x_internal_source,
              e_s_rho2_x = e_s_rho2_x_internal_source,
              c_ps = c_ps,
              piVal = piVal
            )
          
          SEff <-
            ComputeEfficientScore_CPP(
              betaHat,
              sData,
              tData,
              piVal,
              tData,
              coef_y_x_s_hat,
              sigma_y_x_s_hat,
              ispar,
              GH_Materials,
              xList,
              wList
            )
          
          thetaHat <-
            COMPUTE_THETA_Pert_CPP(
              num_of_target,
              piVal,
              rhoValSource,
              c_ps,
              yVec,
              betaHat,
              e_t_tau,
              tau_x_internal_all,
              tau_x_external,
              e_s_rho2_psi_x_internal_all,
              e_s_rho2_x_all,
              e_s_rho2_psi_x_external,
              e_s_rho2_x_ext,
              rexpVec = rexpVec
            )
          thetaHatVec[i] <- thetaHat
        }
        return(sd(thetaHatVec))
      },
      mc.cores = detectCores())
    output <-
      list(Rich = results_output, Sd = unlist(results_output_pert))
    saveRDS(output,
            file = paste("theta_output_ratio", mnratio, "_n", n, "_.RDS", sep = ""))
  }
}

################################################
# Extract Information from the Simulation Results
################################################

library(dplyr)
library(stringr)
Mu_Y_T <- 1.5

read_in_data <- function(filename, dir)
{
  rt <- str_match(filename, "ratio(.*?)_n")[2] %>% as.numeric()
  n <- str_match(filename, "_n(.*?)_.RDS")[2] %>% as.numeric()
  dat <- readRDS(paste(dir, filename, sep = ""))
  
  return(list(rt = rt, dat = dat, n = n))
}

dat_dir <- "dat2/"
filenames <- list.files(dat_dir)
out <- NULL

for (filename in filenames)
{
  Dat <- read_in_data(filename, dat_dir)
  
  rt <- Dat$rt
  n <- Dat$n
  rtDat <- Dat$dat
  
  richInfo <- rtDat$Rich
  sapply(richInfo, function(x) {
    x$ThetaHat
  }) -> thetaHatVec
  
  sapply(richInfo, function(x) {
    x$ThetaNaive
  }) -> thetaHatNaiveVec
  
  sapply(richInfo, function(x) {
    x$SdHat
  }) -> sdHatVec
  
  sapply(richInfo, function(x) {
    x$CP
  }) -> CPVec
  
  sdPertVec <- rtDat$Sd
  
  # Bias of the efficient estimator
  biasEff <- mean(thetaHatVec) - Mu_Y_T
  # SE of the efficient estimator
  seEff <- sd(thetaHatVec)
  
  # Bias of the naive estimator
  biasNaive <- mean(thetaHatNaiveVec) - Mu_Y_T
  # SE of the naive estimator
  seNaive <- sd(thetaHatNaiveVec)
  
  # Coverage of the Semi-inference
  cpFmla <- mean(CPVec)
  # Mean Sd of Semi-inference
  sdMeanSemi <- mean(sdHatVec)
  
  # Coverage using perturbation
  lwbPert <- thetaHatVec - 1.96 * sdPertVec
  upbPert <- thetaHatVec + 1.96 * sdPertVec
  cpPert <-
    sum(lwbPert < Mu_Y_T & upbPert > Mu_Y_T) / length(thetaHatVec)
  # Mean Sd of perturbation
  sdMeanPert <- mean(sdPertVec)
  
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
        sdMeanPert
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
    "SDPert"
  )

################################################
# Plot Simulation Results
################################################
library(ggplot2)
# Bias
outBias <- out[, c(1, 2, 3, 5)] %>% as.data.frame
outBias <-
  reshape2::melt(
    outBias,
    id.vars = c("ratio", "n"),
    value.name = "Bias",
    measure.vars = c("BiasEff", "BiasNaive"),
    variable.name = "Method"
  )
outBias$Method <-
  factor(
    outBias$Method,
    levels = c("BiasEff", "BiasNaive"),
    labels = c("Efficient", "Naive")
  )
outBias %>% ggplot(aes(x = n, y = Bias)) + geom_line(aes(col = Method, linetype = Method)) +
  geom_point(aes(col = Method, shape = Method)) + facet_wrap( ~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + geom_hline(aes(yintercept = 0), linetype = "dashed")

# SE
outSE <- out[, c(1, 2, 4, 6, 8, 10)] %>% as.data.frame
outSE <-
  reshape2::melt(
    outSE,
    id.vars = c("ratio", "n"),
    value.name = "SE",
    measure.vars = c("SEEff", "SENaive", "SDFormula", "SDPert"),
    variable.name = "Method"
  )
outSE$Method <-
  factor(
    outSE$Method,
    levels = c("SEEff", "SENaive", "SDFormula", "SDPert"),
    labels = c(
      "Efficient (Emp.)",
      "Naive (Emp.)",
      "Eff. Formula (Est.)",
      "Eff. Perturbation (Est.)"
    )
  )
outSE %>% ggplot(aes(x = n, y = SE)) + geom_line(aes(col = Method, linetype = Method)) +
  geom_point(aes(col = Method, shape = Method)) + facet_wrap( ~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + ylab("Empirical/Estimated Standard Deviation")

# Coverage: Formula vs. Perturbation
outCP <- out[, c(1, 2, 7, 9)] %>% as.data.frame
outCP <-
  reshape2::melt(
    outCP,
    id.vars = c("ratio", "n"),
    value.name = "CP",
    measure.vars = c("CPFormula", "CPPert"),
    variable.name = "Method"
  )
outCP$Method <-
  factor(
    outCP$Method,
    levels = c("CPFormula", "CPPert"),
    labels = c("Formula", "Perturbation")
  )
outCP %>% ggplot(aes(x = n, y = CP)) + geom_line(aes(col = Method, linetype = Method)) +
  geom_point(aes(col = Method, shape = Method)) + facet_wrap( ~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + geom_hline(aes(yintercept = 0.95), linetype = "dashed") + ylab("Coverage Probability")
