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

Generate_Y_Given_X <- function(xMat, sData, Method)
{
  if (Method == "logit" | Method == "probit")
  {
    fglm <-
      glm(Y ~ .,
          family = binomial(link = Method),
          data = as.data.frame(sData))
    probVec <-
      predict.glm(fglm, as.data.frame(xMat), type = "response")
  }
  
  return(probVec)
}

E_S_RHO_Given_X_Binary_ <- function(betaVal, yFittedGivenX, pwr)
{
  prY1 <- yFittedGivenX
  prY0 <- 1 - prY1
  out <-
    prY0 * exp(pwr * betaVal * 0) + prY1 * exp(pwr * betaVal * 1)
  
  return(out)
}

E_S_RHO_Given_X_Binary <-
  function(betaVal, xMat, sData, pwr, Method)
  {
    yFittedGivenX <- Generate_Y_Given_X(xMat, sData, Method)
    out <- E_S_RHO_Given_X_Binary_(betaVal, yFittedGivenX, pwr)
    
    return(out)
  }

E_T_RHO_Given_X_Binary <-
  function(betaVal, xMat, sData, Method)
  {
    es2 <-
      E_S_RHO_Given_X_Binary(betaVal, xMat, sData, 2, Method)
    es <-
      E_S_RHO_Given_X_Binary(betaVal, xMat, sData, 1, Method)
    
    return(es2 / es)
  }

Compute_Tau_Binary <-
  function(betaVal,
           xMat,
           sData,
           c_ps,
           piVal,
           Method)
  {
    e_t_rho_given_x <-
      E_T_RHO_Given_X_Binary(betaVal, xMat, sData, Method)
    tmp <- e_t_rho_given_x / piVal / c_ps
    
    return(tmp / (tmp + 1 / (1 - piVal)))
  }

# For computing function S
E_S_RHO_Y_Given_X_Binary <-
  function(betaVal, xMat, sData, Method)
  {
    probVec <- Generate_Y_Given_X(xMat, sData, Method)
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
           xMat,
           sData,
           c_ps,
           Method)
  {
    e_s_rho_y_given_x <-
      E_S_RHO_Y_Given_X_Binary(betaVal, xMat, sData, Method)
    e_s_rho_given_x <-
      E_S_RHO_Given_X_Binary(betaVal, xMat, sData, 1, Method)
    e_s_rho_y <- E_S_RHO_Y_Binary(betaVal, sData)
    
    e_s_rho_y_given_x / e_s_rho_given_x - e_s_rho_y / c_ps
  }

E_T_Tau_S_Binary <-
  function(betaVal,
           xMat_t,
           sData,
           c_ps,
           piVal,
           Method)
  {
    SVec <-
      Compute_S_Binary(betaVal, xMat_t, sData, c_ps, Method)
    rhoVec <-
      Compute_Tau_Binary(betaVal, xMat_t, sData, c_ps, piVal, Method)
    mean(SVec * rhoVec)
  }

E_T_Tau_Binary <-
  function(betaVal, xMat_t, sData, c_ps, piVal, Method)
  {
    mean(Compute_Tau_Binary(betaVal, xMat_t, sData, c_ps, piVal, Method))
  }

Compute_B_Binary <-
  function(betaVal,
           xMat,
           sData,
           xMat_t,
           c_ps,
           piVal,
           Method)
  {
    tau_x <-
      Compute_Tau_Binary(betaVal, xMat, sData, c_ps, piVal, Method)
    s_x <- Compute_S_Binary(betaVal, xMat, sData, c_ps, Method)
    e_t_tau_s <-
      E_T_Tau_S_Binary(betaVal, xMat_t, sData, c_ps, piVal, Method)
    e_t_tau <-
      E_T_Tau_Binary(betaVal, xMat_t, sData, c_ps, piVal, Method)
    return(-1 * (1 - piVal) * (1 - tau_x) * (s_x - e_t_tau_s / (e_t_tau -
                                                                  1)))
  }

Compute_S_Eff_Binary <-
  function(betaVal, tData, sData, c_ps, piVal, Method)
  {
    # Labeled Data
    yVec <- sData[, "Y"]
    rhoVec <- exp(betaVal * yVec)
    mult1 <- 1 / piVal * rhoVec / c_ps
    b1 <-
      Compute_B_Binary(betaVal, sData[,-1], sData, tData, c_ps, piVal, Method)
    
    # Unlabeled Data
    mult2 <- -1 / (1 - piVal)
    b2 <-
      Compute_B_Binary(betaVal, tData, sData, tData, c_ps, piVal, Method)
    
    S_Eff <- c(mult1 * b1, mult2 * b2)
    
    return(as.numeric(S_Eff))
  }

Compute_S_Eff_Sum <-
  function(betaVal, tData, sData, c_ps, piVal, Method)
  {
    S_Eff <-
      Compute_S_Eff_Binary(betaVal, tData, sData, c_ps, piVal, Method)
    return(mean(S_Eff) ^ 2)
  }
