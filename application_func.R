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
  lm(sofa ~ ., data = sDat) -> fit
  return(fit)
}

E_S_RHO_App <- function(betaVal, sDat)
{
  sofaScore <- sDat$sofa
  mean(exp(betaVal * sofaScore))
}

E_S_RHO_Given_X_App <- function(betaVal, pyxs, pwr, ghDat, xMat)
{
  coef_yx_s <- coef(pyxs)
  sigma_yx_s <- sigma(pyxs)
  
  xGH <- ghDat$x
  wGH <- ghDat$w
  
  muVec <- cbind(1, xMat) %*% matrix(c(coef_yx_s), ncol = 1)
  sofaMat <- matrix(muVec, nrow = nrow(xMat), ncol = length(xGH))
  addedMat <-
    matrix(
      sqrt(2) * sigma_yx_s * xGH,
      ncol = length(xGH),
      nrow = nrow(xMat),
      byrow = T
    )
  sofaMat <- sofaMat + addedMat
  wMat <-
    matrix(wGH,
           nrow = nrow(xMat),
           ncol = length(wGH),
           byrow = T)
  rowSums(exp(pwr * betaVal * sofaMat) * wMat) / sqrt(pi)
}

E_T_RHO_Given_X_App <- function(betaVal, pyxs, ghDat, xMat)
{
  es2 <- E_S_RHO_Given_X_App(betaVal, pyxs, 2, ghDat, xMat)
  es <- E_S_RHO_Given_X_App(betaVal, pyxs, 1, ghDat, xMat)
  
  es2 / es
}

Compute_Tau_App <- function(betaVal, pyxs, ghDat, xMat, c_ps, piVal)
{
  e_t_rho_given_x <- E_T_RHO_Given_X_App(betaVal, pyxs, ghDat, xMat)
  tmp <- e_t_rho_given_x / piVal / c_ps
  
  tmp / (tmp + 1 / (1 - piVal))
}

E_S_RHO_Y_Given_X_App <-
  function(betaVal, pyxs, ghDat, xMat, pwr = 1)
  {
    coef_yx_s <- coef(pyxs)
    sigma_yx_s <- sigma(pyxs)
    
    xGH <- ghDat$x
    wGH <- ghDat$w
    
    muVec <- cbind(1, xMat) %*% matrix(coef_yx_s, ncol = 1)
    sofaMat <- matrix(muVec, nrow = nrow(xMat), ncol = length(xGH))
    addedMat <-
      matrix(
        sqrt(2) * sigma_yx_s * xGH,
        ncol = length(xGH),
        nrow = nrow(xMat),
        byrow = T
      )
    sofaMat <- sofaMat + addedMat
    wMat <-
      matrix(wGH,
             nrow = nrow(xMat),
             ncol = length(wGH),
             byrow = T)
    rowSums(sofaMat * exp(pwr * betaVal * sofaMat) * wMat) / sqrt(pi)
  }

E_S_RHO_Y_App <- function(betaVal, sDat)
{
  sofaScore <- sDat$sofa
  rhoScore <- exp(sofaScore * betaVal)
  mean(sofaScore * rhoScore)
}

Compute_S_App <- function(betaVal, pyxs, ghDat, xMat, c_ps, sDat)
{
  e_s_rho_y_given_x <-
    E_S_RHO_Y_Given_X_App(betaVal, pyxs, ghDat, xMat)
  e_s_rho_given_x <-
    E_S_RHO_Given_X_App(betaVal, pyxs, 1, ghDat, xMat)
  e_s_rho_y <- E_S_RHO_Y_App(betaVal, sDat)
  
  e_s_rho_y_given_x / e_s_rho_given_x - e_s_rho_y / c_ps
}

E_T_Tau_S_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat_t,
           c_ps,
           sDat,
           piVal)
  {
    SVec <- Compute_S_App(betaVal, pyxs, ghDat, xMat_t, c_ps, sDat)
    rhoVec <-
      Compute_Tau_App(betaVal, pyxs, ghDat, xMat_t, c_ps, piVal)
    mean(SVec * rhoVec)
  }

E_T_Tau_App <- function(betaVal, pyxs, ghDat, xMat_t, c_ps, piVal)
{
  mean(Compute_Tau_App(betaVal, pyxs, ghDat, xMat_t, c_ps, piVal))
}

Compute_B_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat,
           xMat_t,
           c_ps,
           piVal,
           sDat)
  {
    tau_x <- Compute_Tau_App(betaVal, pyxs, ghDat, xMat, c_ps, piVal)
    s_x <- Compute_S_App(betaVal, pyxs, ghDat, xMat, c_ps, sDat)
    e_t_tau_s <-
      E_T_Tau_S_App(betaVal, pyxs, ghDat, xMat_t, c_ps, sDat, piVal)
    e_t_tau <-
      E_T_Tau_App(betaVal, pyxs, ghDat, xMat_t, c_ps, piVal)
    (-1) * (1 - piVal) * (1 - tau_x) * (s_x - e_t_tau_s / (e_t_tau - 1))
  }

Compute_S_Eff_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat_s,
           xMat_t,
           c_ps,
           piVal,
           sDat)
  {
    c_ps <- E_S_RHO_App(betaVal, sDat)
    sofaScore <- sDat$sofa
    rhoVec <- exp(betaVal * sofaScore)
    
    mult1 <- 1 / piVal * rhoVec / c_ps
    b1 <-
      Compute_B_App(betaVal,
                    pyxs,
                    ghDat,
                    xMat_s,
                    xMat_t,
                    c_ps,
                    piVal,
                    sDat)
    
    mult2 <- -1 / (1 - piVal)
    b2 <-
      Compute_B_App(betaVal,
                    pyxs,
                    ghDat,
                    xMat_t,
                    xMat_t,
                    c_ps,
                    piVal,
                    sDat)
    
    S_Eff <- c(mult1 * b1, mult2 * b2)
    
    return(as.numeric(S_Eff))
  }

Compute_S_Eff_Naive_App <- function(betaVal,
                                    pyxs,
                                    ghDat,
                                    xMat_s,
                                    xMat_t,
                                    c_ps,
                                    piVal,
                                    sDat)
{
  c_ps <- E_S_RHO_App(betaVal, sDat)
  sofaScore <- sDat$sofa
  rhoVec <- exp(betaVal * sofaScore)
  
  mult1 <- 1 / piVal * rhoVec / c_ps
  mult2 <- -1 / (1 - piVal)
  
  b1 <- predict(pyxs, newdata = as.data.frame(xMat_s))
  b2 <- predict(pyxs, newdata = as.data.frame(xMat_t))
  
  S_Eff <- c(mult1 * b1, mult2 * b2)
  
  return(as.numeric(S_Eff))
}

Compute_S_Eff_Sum_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat_s,
           xMat_t,
           piVal,
           sDat)
  {
    xMat_s <- as.matrix(xMat_s)
    xMat_t <- as.matrix(xMat_t)
    
    c_ps <- E_S_RHO_App(betaVal, sDat)
    S_Eff <-
      Compute_S_Eff_App(betaVal,
                        pyxs,
                        ghDat,
                        xMat_s,
                        xMat_t,
                        c_ps,
                        piVal,
                        sDat)
    mean(S_Eff) ^ 2
  }

Compute_S_Eff_Pert_Sum_App <- function(betaVal,
                                       pyxs,
                                       ghDat,
                                       xMat_s,
                                       xMat_t,
                                       piVal,
                                       sDat,
                                       rexpVec)
{
  xMat_s <- as.matrix(xMat_s)
  xMat_t <- as.matrix(xMat_t)
  
  c_ps <- E_S_RHO_App(betaVal, sDat)
  
  Compute_S_Eff_App(betaVal,
                    pyxs,
                    ghDat,
                    xMat_s,
                    xMat_t,
                    c_ps,
                    piVal,
                    sDat) -> seff
  mean(seff * rexpVec) ^ 2
}


Compute_S_Eff_Sum_Naive_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat_s,
           xMat_t,
           piVal,
           sDat)
  {
    xMat_s <- as.matrix(xMat_s)
    xMat_t <- as.matrix(xMat_t)
    c_ps <- E_S_RHO_App(betaVal, sDat)
    S_Eff <-
      Compute_S_Eff_Naive_App(betaVal,
                              pyxs,
                              ghDat,
                              xMat_s,
                              xMat_t,
                              c_ps,
                              piVal,
                              sDat)
    mean(S_Eff) ^ 2
    
  }

Compute_S_Eff_Sum_Naive_Pert_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat_s,
           xMat_t,
           piVal,
           sDat,
           rexpVec)
  {
    xMat_s <- as.matrix(xMat_s)
    xMat_t <- as.matrix(xMat_t)
    c_ps <- E_S_RHO_App(betaVal, sDat)
    S_Eff <-
      Compute_S_Eff_Naive_App(betaVal,
                              pyxs,
                              ghDat,
                              xMat_s,
                              xMat_t,
                              c_ps,
                              piVal,
                              sDat)
    mean(S_Eff * rexpVec) ^ 2
    
  }

Compute_Beta_Var_App <- function(betaVal,
                                 pyxs,
                                 ghDat,
                                 xMat_s,
                                 xMat_t,
                                 c_ps,
                                 piVal,
                                 sDat,
                                 proposed = T)
{
  if (!proposed)
  {
    seff <- Compute_S_Eff_Naive_App(betaVal,
                                    pyxs,
                                    ghDat,
                                    xMat_s,
                                    xMat_t,
                                    c_ps,
                                    piVal,
                                    sDat)
  }
  else
  {
    seff <- Compute_S_Eff_App(betaVal,
                              pyxs,
                              ghDat,
                              xMat_s,
                              xMat_t,
                              c_ps,
                              piVal,
                              sDat)
  }
  1 / mean(seff ^ 2) / length(seff)
}

Compute_BVec_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat_t,
           c_ps,
           piVal,
           xMat)
  {
    E_T_Tau_App(betaVal, pyxs, ghDat, xMat_t, c_ps, piVal) -> e_t_tau
    Compute_Tau_App(betaVal, pyxs, ghDat, xMat, c_ps, piVal) -> tau_x
    (e_t_tau - tau_x) / (1 - e_t_tau)
  }

Compute_AVec_Phi_Y_App <-
  function(betaVal,
           pyxs,
           ghDat,
           xMat_t,
           c_ps,
           piVal,
           xMat)
  {
    tau_x_t <-
      Compute_Tau_App(betaVal, pyxs, ghDat, xMat_t, c_ps, piVal)
    e_t_tau <- mean(tau_x_t)
    tau_x <-
      Compute_Tau_App(betaVal, pyxs, ghDat, xMat, c_ps, piVal)
    e_s_rho_2_y_given_x_t <-
      E_S_RHO_Y_Given_X_App(betaVal, pyxs, ghDat, xMat_t, pwr = 2)
    e_s_rho_2_given_x_t <-
      E_S_RHO_Given_X_App(betaVal, pyxs, 2, ghDat, xMat_t)
    mean(e_s_rho_2_y_given_x_t / e_s_rho_2_given_x_t * tau_x_t) -> e_t_tau_xx
    
    e_s_rho_2_y_given_x <-
      E_S_RHO_Y_Given_X_App(betaVal, pyxs, ghDat, xMat, pwr = 2)
    e_s_rho_2_given_x <-
      E_S_RHO_Given_X_App(betaVal, pyxs, 2, ghDat, xMat)
    t_tau_xx <- tau_x * e_s_rho_2_y_given_x / e_s_rho_2_given_x
    
    A_x <-
      (1 - tau_x) / (1 - e_t_tau) * e_t_tau_xx - tau_x * t_tau_xx
    
    return(A_x)
  }

Compute_Eff_Theta_App <- function(theta,
                                  betaVal,
                                  pyxs,
                                  ghDat,
                                  xMat_t,
                                  c_ps,
                                  piVal,
                                  xMat_s,
                                  sDat)
{
  xMat_s <- as.matrix(xMat_s)
  xMat_t <- as.matrix(xMat_t)
  
  sofaScore <- sDat$sofa
  AVec_s <- Compute_AVec_Phi_Y_App(betaVal,
                                   pyxs,
                                   ghDat,
                                   xMat_t,
                                   c_ps,
                                   piVal,
                                   xMat_s)
  AVec_t <- Compute_AVec_Phi_Y_App(betaVal,
                                   pyxs,
                                   ghDat,
                                   xMat_t,
                                   c_ps,
                                   piVal,
                                   xMat_t)
  BVec_s <- Compute_BVec_App(betaVal,
                             pyxs,
                             ghDat,
                             xMat_t,
                             c_ps,
                             piVal,
                             xMat_s)
  BVec_t <- Compute_BVec_App(betaVal,
                             pyxs,
                             ghDat,
                             xMat_t,
                             c_ps,
                             piVal,
                             xMat_t)
  part1 <-
    1 / piVal * exp(betaVal * sofaScore) / c_ps * (sofaScore - theta + AVec_s -
                                                     BVec_s * theta)
  part2 <- -1 / (1 - piVal) * (AVec_t - BVec_t * theta)
  c(part1, part2)
}

Compute_Eff_Theta_Sum_App <- function(theta,
                                      betaVal,
                                      pyxs,
                                      ghDat,
                                      xMat_t,
                                      c_ps,
                                      piVal,
                                      xMat_s,
                                      sDat)
{
  Compute_Eff_Theta_App(theta,
                        betaVal,
                        pyxs,
                        ghDat,
                        xMat_t,
                        c_ps,
                        piVal,
                        xMat_s,
                        sDat) -> theta_eff
  mean(theta_eff) ^ 2
}

Compute_Eff_Theta_Sum_Pert_App <- function(theta,
                                           betaVal,
                                           pyxs,
                                           ghDat,
                                           xMat_t,
                                           c_ps,
                                           piVal,
                                           xMat_s,
                                           sDat,
                                           rexpVec)
{
  Compute_Eff_Theta_App(theta,
                        betaVal,
                        pyxs,
                        ghDat,
                        xMat_t,
                        c_ps,
                        piVal,
                        xMat_s,
                        sDat) -> theta_eff
  mean(theta_eff * rexpVec) ^ 2
}

Compute_Eff_Theta_Naive_Sum_App <- function(theta,
                                            betaVal,
                                            c_ps,
                                            piVal,
                                            sDat)
{
  sofaScore <- sDat$sofa
  1 / piVal * exp(betaVal * sofaScore) / c_ps * (sofaScore - theta) -> eff
  mean(eff) ^ 2
}

Compute_Eff_Theta_Naive_Sum_Pert_App <- function(theta,
                                                 betaVal,
                                                 c_ps,
                                                 piVal,
                                                 sDat,
                                                 rexpVec)
{
  sofaScore <- sDat$sofa
  1 / piVal * exp(betaVal * sofaScore) / c_ps * (sofaScore - theta) -> eff
  mean(eff * rexpVec[1:length(eff)]) ^ 2
}

print_results <- function(fit_beta, fit_theta, var_beta, theta_B)
{
  betaHat <- fit_beta$par
  betaSE <- sqrt(var_beta)
  betaLB <- betaHat - 1.96 * betaSE
  betaUB <- betaHat + 1.96 * betaSE
  
  thetaHat <- fit_theta$par
  theta_b_vec <- unlist(theta_B)
  thetaSE <- sd(theta_b_vec)
  thetaLB <- thetaHat - 1.96 * thetaSE
  thetaUB <- thetaHat + 1.96 * thetaSE
  
  outMat <- matrix(nrow = 2, ncol = 4)
  colnames(outMat) <- c("Estimate", "SE", "Lower", "Upper")
  rownames(outMat) <- c("beta", "theta")
  outMat[1,] <- c(betaHat, betaSE, betaLB, betaUB)
  outMat[2,] <- c(thetaHat, thetaSE, thetaLB, thetaUB)
  
  return(outMat)
}
