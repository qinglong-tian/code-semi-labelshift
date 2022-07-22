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


Generate_Y_Given_X <- function(xMats, xMatt, sData, Method = "logit")
{
  if (Method == "logit" | Method == "probit")
  {
    fglm <-
      glm(Y ~ .,
          family = binomial(link = Method),
          data = as.data.frame(sData))
    probVecs <-
      predict.glm(fglm, as.data.frame(xMats), type = "response")
    probVect <- 
      predict.glm(fglm, as.data.frame(xMatt), type = "response")
  }
  else if (Method == "nb")
  {
    nbFit <-
      train(sData[, -1],
            as.factor(sData[, "Y"]),
            'nb',
            trControl = trainControl(method = 'cv', number = 10))
    probVecs <- predict(nbFit, xMats, type = "prob")[, 2]
    probVect <- predict(nbFit, xMatt, type = "prob")[, 2]
  }
  else if (Method == "rf")
  {
    control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
    metric <- "Accuracy"
    sData <- as.data.frame(sData)
    sData[,"Y"] <- as.factor(sData[,"Y"])
    rfFit <- train(Y~.,
                   data = sData,
                   method = 'rf',
                   metric=metric, tuneLength=15, trControl=control)
    probVecs <- predict(rfFit, xMats, type = "prob")[, 2]
    probVect <- predict(rfFit, xMatt, type = "prob")[, 2]
  }
  return(list(
    prob_s = probVecs,
    prob_t = probVect
  ))
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
  function(betaVal, sData, c_ps, piVal, yFittedGivenX, yFittedGivenX_tDat)
  {
    # Labeled Data
    yVec <- sData[, "Y"]
    rhoVec <- exp(betaVal * yVec)
    mult1 <- 1 / piVal * rhoVec / c_ps
    
    b1 <-
      Compute_B_Binary(betaVal, yFittedGivenX, yFittedGivenX_tDat, sData, c_ps, piVal)
    
    # Unlabeled Data
    mult2 <- -1 / (1 - piVal)
    b2 <-
      Compute_B_Binary(betaVal, yFittedGivenX_tDat, yFittedGivenX_tDat, sData, c_ps, piVal)
    
    S_Eff <- c(mult1 * b1, mult2 * b2)
    
    return(as.numeric(S_Eff))
  }

Compute_S_Eff_Sum <-
  function(betaVal, sData, c_ps, piVal, yFittedGivenX, yFittedGivenX_tDat)
  {
    S_Eff <-
      Compute_S_Eff_Binary(betaVal, sData, c_ps, piVal, yFittedGivenX, yFittedGivenX_tDat)
    return(mean(S_Eff) ^ 2)
  }
