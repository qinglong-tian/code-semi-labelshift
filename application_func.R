Compute_Marginal_of_SOFA <- function(sData, tData)
{
  tSofa <- unique(tData$sofa)
  sSofa <- unique(sData$sofa)
  
  if (any(tSofa %in% sSofa == FALSE))
  {
    stop("The support of target data is not a subset of that of the source data.")
  }
  
  pSofaS <- table(sData$sofa) / length(sData$sofa)
  pSofat <- table(tData$sofa) / length(tData$sofa)
  sofaClassesS <- dimnames(pSofaS)[[1]]
  sofaClassesT <- dimnames(pSofat)[[1]]
  
  dfOut <- data.frame()
  for (i in 1:length(sofaClassesS))
  {
    sofaclass <- sofaClassesS[i]
    pt <- ifelse(sofaclass %in% sofaClassesT, pSofat[sofaclass], 0)
    ps <- pSofaS[sofaclass]
    ratio <- pt / ps
    dfOut[i, 1] <- sofaclass
    dfOut[i, 2] <- ratio
    dfOut[i, 3] <- ps
    dfOut[i, 4] <- pt
  }
  colnames(dfOut) <- c("SOFA", "Ratio", "Ps", "Pt")
  return(dfOut)
}

Compute_Pr_R1_Given_Y <- function(sData, tData)
{
  num_of_s <- nrow(sData)
  num_of_t <- nrow(tData)
  
  rMarginalRatio <- num_of_t / num_of_s
  df <- Compute_Marginal_of_SOFA(sData, tData)
  dff <- df %>% mutate(rProb = 1 / (Ratio * rMarginalRatio + 1))
  
  return(dff)
}

Compute_GLM_Pr_R1_Given_Y <- function(sData, tData, rSofaLogistic)
{
  dff <- Compute_Pr_R1_Given_Y(sData, tData)
  as.numeric(dff$SOFA) -> sofa_num
  predict(rSofaLogistic,
          newdata = data.frame(sofa = sofa_num),
          type = "response") -> predProb
  predProb <- unname(predProb)
  dff %>% mutate(rProbPred = predProb) -> dff
  
  return(dff)
}

Fit_Sofa_Score <- function(sDat)
{
  sDat %>% mutate(sofa = factor(sofa, levels = paste(0:10))) %>% mutate (sofa = relevel(sofa, ref = "10")) -> sDat1
  multinom_model <- multinom(sofa ~ ., data = sDat1)
  return(multinom_model)
}

E_S_RHO_App <- function(betaVal, sDat)
{
  sofaScore <- sDat$sofa
  mean(exp(betaVal * sofaScore))
}

Pred_Sofa_Given_X_S <- function(sDat, tDat)
{
  multi_nominal <- Fit_Sofa_Score(sDat)
  prob_xs <- predict(multi_nominal, newdata = sDat, type = "probs")
  prob_xt <- predict(multi_nominal, newdata = tDat, type = "probs")
  
  return(list(prob_list_xs = prob_xs,
              prob_list_xt = prob_xt))
}

E_S_RHO_Given_X_App <- function(betaVal, prob_list_x, pwr)
{
  matrix(
    exp(pwr * betaVal * (0:10)),
    ncol = length(0:10),
    nrow = nrow(prob_list_x),
    byrow = T
  ) -> rhoMat
  rowSums(rhoMat * prob_list_x)
}

E_T_RHO_Given_X_App <- function(betaVal, prob_list_x)
{
  es2 <- E_S_RHO_Given_X_App(betaVal, prob_list_x, 2)
  es <- E_S_RHO_Given_X_App(betaVal, prob_list_x, 1)
  
  es2 / es
}

Compute_Tau_App <- function(betaVal, prob_list_x, c_ps, piVal)
{
  e_t_rho_given_x <- E_T_RHO_Given_X_App(betaVal, prob_list_x)
  tmp <- e_t_rho_given_x / piVal / c_ps
  
  tmp / (tmp + 1 / (1 - piVal))
}

E_S_RHO_Y_Given_X_App <- function(betaVal, prob_list_x)
{
  matrix((0:10) * exp((0:10) * betaVal),
         ncol = length(0:10),
         nrow = nrow(prob_list_x),
         byrow = T
  ) -> rhoMat
  
  rowSums(rhoMat * prob_list_x)
}

E_S_RHO_Y_App <- function(betaVal, sDat)
{
  sofaScore <- sDat$sofa
  rhoScore <- exp(sofaScore * betaVal)
  mean(sofaScore * rhoScore)
}

Compute_S_App <- function(betaVal, prob_list_x, sDat, c_ps)
{
  e_s_rho_y_given_x <- E_S_RHO_Y_Given_X_App(betaVal, prob_list_x)
  e_s_rho_given_x <- E_S_RHO_Given_X_App(betaVal, prob_list_x, 1)
  e_s_rho_y <- E_S_RHO_Y_App(betaVal, sDat)
  
  e_s_rho_y_given_x / e_s_rho_given_x - e_s_rho_y / c_ps
}

E_T_Tau_S_App <- function(betaVal, prob_list_xt, sDat, c_ps, piVal)
{
  SVec <- Compute_S_App(betaVal, prob_list_xt, sDat, c_ps)
  rhoVec <- Compute_Tau_App(betaVal, prob_list_xt, c_ps, piVal)
  mean(SVec * rhoVec)
}

E_T_Tau_App <- function(betaVal, prob_list_xt, c_ps, piVal)
{
  mean(Compute_Tau_App(betaVal, prob_list_xt, c_ps, piVal))
}

Compute_B_App <-
  function(betaVal,
           prob_list_x,
           prob_list_xt,
           sDat,
           c_ps,
           piVal)
  {
    tau_x <- Compute_Tau_App(betaVal, prob_list_x, c_ps, piVal)
    s_x <- Compute_S_App(betaVal, prob_list_x, sDat, c_ps)
    e_t_tau_s <-
      E_T_Tau_S_App(betaVal, prob_list_xt, sDat, c_ps, piVal)
    e_t_tau <- E_T_Tau_App(betaVal, prob_list_xt, c_ps, piVal)
    (-1) * (1 - piVal) * (1 - tau_x) * (s_x - e_t_tau_s / (e_t_tau - 1))
  }

Compute_S_Eff_App <-
  function(betaVal,
           sDat,
           c_ps,
           piVal,
           prob_list_xs,
           prob_list_xt)
  {
    c_ps <- E_S_RHO_App(betaVal, sDat)
    sofaScore <- sDat$sofa
    rhoVec <- exp(betaVal * sofaScore)
    
    mult1 <- 1 / piVal * rhoVec / c_ps
    b1 <-
      Compute_B_App(betaVal, prob_list_xs, prob_list_xt, sDat, c_ps, piVal)
    
    mult2 <- -1 / (1 - piVal)
    b2 <-
      Compute_B_App(betaVal, prob_list_xt, prob_list_xt, sDat, c_ps, piVal)
    
    S_Eff <- c(mult1 * b1, mult2 * b2)
    
    return(as.numeric(S_Eff))
  }

Compute_S_Eff_Sum_App <-
  function(betaVal,
           sDat,
           piVal,
           prob_list_xs,
           prob_list_xt)
  {
    c_ps <- E_S_RHO_App(betaVal, sDat)
    S_Eff <-
      Compute_S_Eff_App(betaVal, sDat, c_ps, piVal, prob_list_xs, prob_list_xt)
    mean(S_Eff) ^ 2
  }
