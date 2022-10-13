library(e1071)
library(caret)
# Just a glm
Fit_Binary_Y_With_X <- function(YXMat, probit = F)
{
  link_name <- ifelse(probit, "probit", "logit")
  fit <-
    glm(Y ~ .,
        data = as.data.frame(YXMat),
        family = binomial(link = link_name))
  return(fit)
}

# Functions to derive the efficient estimator
E_S_RHO_Binary <- function(betaVal, sData)
{
  yVec <- sData[, "Y"]
  prY1 <- mean(yVec)
  prY0 <- 1 - prY1
  out <- prY0 * exp(betaVal * 0) + prY1 * exp(betaVal * 1)
  
  return(out)
}

Generate_Y_Given_X <-
  function(xMats, xMatt, sData, Method = "logit")
  {
    fglm <-
      glm(Y ~ .,
          family = binomial(link = Method),
          data = as.data.frame(sData))
    probVecs <-
      predict.glm(fglm, as.data.frame(xMats), type = "response")
    probVect <-
      predict.glm(fglm, as.data.frame(xMatt), type = "response")
    
    return(list(prob_s = probVecs,
                prob_t = probVect))
  }

E_S_RHO_Given_X_Binary <- function(betaVal, yFittedGivenX, pwr)
{
  prY1 <- yFittedGivenX
  prY0 <- 1 - prY1
  out <-
    prY0 * exp(pwr * betaVal * 0) + prY1 * exp(pwr * betaVal * 1)
  
  return(out)
}

E_T_RHO_Given_X_Binary <-
  function(betaVal, yFittedGivenX)
  {
    es2 <-
      E_S_RHO_Given_X_Binary(betaVal, yFittedGivenX, 2)
    es <-
      E_S_RHO_Given_X_Binary(betaVal, yFittedGivenX, 1)
    
    return(es2 / es)
  }

Compute_Tau_Binary <-
  function(betaVal,
           yFittedGivenX,
           c_ps,
           piVal)
  {
    e_t_rho_given_x <-
      E_T_RHO_Given_X_Binary(betaVal, yFittedGivenX)
    tmp <- e_t_rho_given_x / piVal / c_ps
    
    return(tmp / (tmp + 1 / (1 - piVal)))
  }

# For computing function S
E_S_RHO_Y_Given_X_Binary <-
  function(betaVal, yFittedGivenX)
  {
    probVec <- yFittedGivenX
    pr1 <- probVec
    pr1 * 1 * exp(1 * betaVal) + (1 - pr1) * 0 * exp(0 * betaVal)
  }

E_S_RHO_Y_Binary <- function(betaVal, sData)
{
  yVec <- sData[, "Y"]
  rho <- exp(yVec * betaVal)
  mean(rho * yVec)
}

Compute_S_Binary <-
  function(betaVal,
           yFittedGivenX,
           sData,
           c_ps)
  {
    e_s_rho_y_given_x <-
      E_S_RHO_Y_Given_X_Binary(betaVal, yFittedGivenX)
    e_s_rho_given_x <-
      E_S_RHO_Given_X_Binary(betaVal, yFittedGivenX, 1)
    e_s_rho_y <- E_S_RHO_Y_Binary(betaVal, sData)
    
    e_s_rho_y_given_x / e_s_rho_given_x - e_s_rho_y / c_ps
  }

E_T_Tau_S_Binary <-
  function(betaVal,
           yFittedGivenX_tDat,
           sData,
           c_ps,
           piVal)
  {
    SVec <-
      Compute_S_Binary(betaVal, yFittedGivenX_tDat, sData, c_ps)
    rhoVec <-
      Compute_Tau_Binary(betaVal, yFittedGivenX_tDat, c_ps, piVal)
    mean(SVec * rhoVec)
  }

E_T_Tau_Binary <-
  function(betaVal, yFittedGivenX_tDat, c_ps, piVal)
  {
    mean(Compute_Tau_Binary(betaVal, yFittedGivenX_tDat, c_ps, piVal))
  }

Compute_B_Binary <-
  function(betaVal,
           yFittedGivenX,
           yFittedGivenX_tDat,
           sData,
           c_ps,
           piVal)
  {
    tau_x <-
      Compute_Tau_Binary(betaVal, yFittedGivenX, c_ps, piVal)
    s_x <- Compute_S_Binary(betaVal, yFittedGivenX, sData, c_ps)
    e_t_tau_s <-
      E_T_Tau_S_Binary(betaVal, yFittedGivenX_tDat, sData, c_ps, piVal)
    e_t_tau <-
      E_T_Tau_Binary(betaVal, yFittedGivenX_tDat, c_ps, piVal)
    return(-1 * (1 - piVal) * (1 - tau_x) * (s_x - e_t_tau_s / (e_t_tau -
                                                                  1)))
  }

Compute_S_Eff_Binary <-
  function(betaVal,
           sData,
           c_ps,
           piVal,
           yFittedGivenX,
           yFittedGivenX_tDat)
  {
    c_ps <- E_S_RHO_Binary(betaVal, sData)
    # Labeled Data
    yVec <- sData[, "Y"]
    rhoVec <- exp(betaVal * yVec)
    mult1 <- 1 / piVal * rhoVec / c_ps
    
    b1 <-
      Compute_B_Binary(betaVal,
                       yFittedGivenX,
                       yFittedGivenX_tDat,
                       sData,
                       c_ps,
                       piVal)
    
    # Unlabeled Data
    mult2 <- -1 / (1 - piVal)
    b2 <-
      Compute_B_Binary(betaVal,
                       yFittedGivenX_tDat,
                       yFittedGivenX_tDat,
                       sData,
                       c_ps,
                       piVal)
    
    S_Eff <- c(mult1 * b1, mult2 * b2)
    
    return(as.numeric(S_Eff))
  }

Compute_S_Eff_Sum <-
  function(betaVal,
           sData,
           piVal,
           yFittedGivenX,
           yFittedGivenX_tDat)
  {
    c_ps <- E_S_RHO_Binary(betaVal, sData)
    S_Eff <-
      Compute_S_Eff_Binary(betaVal,
                           sData,
                           c_ps,
                           piVal,
                           yFittedGivenX,
                           yFittedGivenX_tDat)
    return(mean(S_Eff) ^ 2)
  }

Compute_SE_Binary <-
  function(betaVal,
           sData,
           c_ps,
           piVal,
           yFittedGivenX,
           yFittedGivenX_tDat)
  {
    SEffVec <-
      Compute_S_Eff_Binary(betaVal,
                           sData,
                           c_ps,
                           piVal,
                           yFittedGivenX,
                           yFittedGivenX_tDat)
    outVar <- 1 / mean(SEffVec ^ 2) / length(SEffVec)
    sqrt(outVar)
  }

Compute_Lipton_Score <- function(betaVal, probVecS, probVecT, sData)
{
  yVec <- sData[, "Y"]
  n <- length(probVecS)
  m <- length(probVecT)
  piVal <- n / (n + m)
  c_ps <- E_S_RHO_Binary(betaVal, sData)
  
  score <- 0
  for (i in 1:n)
  {
    rhoVal <- exp(betaVal * yVec[i])
    score <- score + 1 / piVal * rhoVal / c_ps * probVecS[i]
  }
  for (i in 1:m)
  {
    score <- score - 1 / (1 - piVal) * probVecT[i]
  }
  
  return((score / (n + m)) ^ 2)
}

Compute_Lipton_Score_Pert <-
  function(betaVal,
           probVecS,
           probVecT,
           sData,
           rexpVec)
  {
    yVec <- sData[, "Y"]
    n <- length(probVecS)
    m <- length(probVecT)
    piVal <- n / (n + m)
    c_ps <- E_S_RHO_Binary(betaVal, sData)
    
    score <- 0
    for (i in 1:n)
    {
      rhoVal <- exp(betaVal * yVec[i])
      score <-
        score + 1 / piVal * rhoVal / c_ps * probVecS[i] * rexpVec[i]
    }
    for (i in 1:m)
    {
      score <- score - 1 / (1 - piVal) * probVecT[i] * rexpVec[n + i]
    }
    return((score / (n + m)) ^ 2)
  }


Estimate_Beta_Lipton <-
  function(probVecS, probVecT, sData, initBeta)
  {
    opt <-
      optim(
        initBeta,
        Compute_Lipton_Score,
        probVecS = probVecS,
        probVecT = probVecT,
        sData = sData,
        method = "Brent",
        lower = -3,
        upper = 3
      )
    return(opt$par)
  }

Estimate_Beta_Lipton_Pert <-
  function(probVecS,
           probVecT,
           sData,
           initBeta,
           rexpVec)
  {
    opt <-
      optim(
        initBeta,
        Compute_Lipton_Score_Pert,
        probVecS = probVecS,
        probVecT = probVecT,
        sData = sData,
        rexpVec = rexpVec,
        method = "Brent",
        lower = -3,
        upper = 3
      )
    return(opt$par)
  }

Compute_B_Matrix_Binary <-
  function(betaVal,
           yFittedGivenX,
           yFittedGivenX_t,
           piVal,
           sData)
  {
    c_ps <- E_S_RHO_Binary(betaVal, sData)
    Tau_X <- Compute_Tau_Binary(betaVal, yFittedGivenX, c_ps, piVal)
    Tau_X_t <-
      Compute_Tau_Binary(betaVal, yFittedGivenX_t, c_ps, piVal)
    E_t_Tau <- mean(Tau_X_t)
    return((E_t_Tau - Tau_X) / (1 - E_t_Tau))
  }

E_S_Rho2_Phi_Given_X_Binary <- function(betaVal, yFittedGivenX)
{
  pr1 <- yFittedGivenX
  out <- exp(2 * betaVal * 1) * 1 * pr1
  return(out)
}

Compute_A_Mat_Binary <-
  function(betaVal,
           yFittedGivenX,
           yFittedGivenX_t,
           piVal,
           sData)
  {
    c_ps <- E_S_RHO_Binary(betaVal, sData)
    Tau_x <- Compute_Tau_Binary(betaVal, yFittedGivenX, c_ps, piVal)
    Tau_X_t <-
      Compute_Tau_Binary(betaVal, yFittedGivenX_t, c_ps, piVal)
    E_t_Tau <- mean(Tau_X_t)
    firstTerm <- (1 - Tau_x) / (1 - E_t_Tau)
    
    E_t_Tau_Phi_X_t <-
      E_S_Rho2_Phi_Given_X_Binary(betaVal, yFittedGivenX_t) / E_S_RHO_Given_X_Binary(betaVal, yFittedGivenX_t, 2)
    secondTerm <- mean(Tau_X_t * E_t_Tau_Phi_X_t)
    
    E_t_Tau_Phi_X <-
      E_S_Rho2_Phi_Given_X_Binary(betaVal, yFittedGivenX) / E_S_RHO_Given_X_Binary(betaVal, yFittedGivenX, 2)
    thirdTerm <- -Tau_x * E_t_Tau_Phi_X
    
    A <- firstTerm * secondTerm + thirdTerm
    
    return(A)
  }

Compute_Phi_Eff_Theta <-
  function(thetaVal,
           betaVal,
           yFittedGivenX,
           yFittedGivenX_t,
           piVal,
           sData)
  {
    c_ps <- E_S_RHO_Binary(betaVal, sData)
    yVec <- sData[, "Y"]
    rhoVec <- exp(betaVal * yVec)
    Compute_A_Mat_Binary(betaVal, yFittedGivenX, yFittedGivenX_t, piVal, sData) -> A_x
    Compute_A_Mat_Binary(betaVal, yFittedGivenX_t, yFittedGivenX_t, piVal, sData) -> A_x_t
    
    Compute_B_Matrix_Binary(betaVal, yFittedGivenX, yFittedGivenX_t, piVal, sData) -> B_x
    Compute_B_Matrix_Binary(betaVal, yFittedGivenX_t, yFittedGivenX_t, piVal, sData) -> B_x_t
    
    1 / piVal * rhoVec / c_ps * (yVec - thetaVal + A_x - B_x * thetaVal) -> part1
    - 1 / (1 - piVal) * (A_x_t - B_x_t * thetaVal) -> part2
    
    return(c(part1, part2))
  }

Compute_Phi_Theta <- function(thetaVal, betaVal, sData, piVal)
{
  yVec <- sData[, "Y"]
  c_ps <- E_S_RHO_Binary(betaVal, sData)
  rhoVec <- exp(betaVal * yVec)
  rhoVec / c_ps * (yVec - thetaVal) / piVal
}

Estimate_Theta_Sum_Eff <- function(thetaVal,
                                   betaVal,
                                   yFittedGivenX,
                                   yFittedGivenX_t,
                                   piVal,
                                   sData)
{
  Compute_Phi_Eff_Theta(thetaVal,
                        betaVal,
                        yFittedGivenX,
                        yFittedGivenX_t,
                        piVal,
                        sData) -> seff
  return(mean(seff) ^ 2)
}

Estiamte_Theta_Sum_Eff_Pert <- function(thetaVal,
                                        betaVal,
                                        yFittedGivenX,
                                        yFittedGivenX_t,
                                        piVal,
                                        sData,
                                        rexpVec)
{
  Compute_Phi_Eff_Theta(thetaVal,
                        betaVal,
                        yFittedGivenX,
                        yFittedGivenX_t,
                        piVal,
                        sData) -> seff
  mean(seff * rexpVec) ^ 2
}

Estimate_Theta_Pert_Optim <- function(betaVal,
                                      yFittedGivenX,
                                      yFittedGivenX_t,
                                      piVal,
                                      sData,
                                      rexpVec)
{
  optim(
    initTheta,
    Estiamte_Theta_Sum_Eff_Pert,
    betaVal = betaVal,
    yFittedGivenX = yFittedGivenX,
    yFittedGivenX_t = yFittedGivenX_t,
    piVal = piVal,
    sData = sData,
    rexpVec = rexpVec,
    method = "Brent",
    lower = 0,
    upper = 1
  ) -> optimOut
  return(optimOut$par)
}

Estimate_Theta_Sum_Naive <- function(thetaVal, betaVal, sData, piVal)
{
  seff <- Compute_Phi_Theta(thetaVal, betaVal, sData, piVal)
  return(mean(seff) ^ 2)
}

Estimate_Theta_Sum_Naive_Pert <-
  function(thetaVal, betaVal, sData, piVal, rexpVec)
  {
    seff <- Compute_Phi_Theta(thetaVal, betaVal, sData, piVal)
    return(mean(seff * rexpVec[1:length(seff)]) ^ 2)
  }

Estimate_optim_Theta_Naive_Pert <-
  function(betaVal, sData, piVal, rexpVec)
  {
    optim(
      initTheta,
      Estimate_Theta_Sum_Naive_Pert,
      betaVal = betaVal,
      sData = sData,
      piVal = piVal,
      rexpVec = rexpVec,
      method = "Brent",
      lower = 0,
      upper = 1
    ) -> opt
    opt$par
  }
