ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

Generate_Y_X_Marginal <- function(n, Mu_YX, SigMat_YX)
{
  if (length(Mu_YX) != nrow(SigMat_YX))
  {
    stop("Wrong Dimensions!")
  }
  yxOut <- MASS::mvrnorm(n, Mu_YX, SigMat_YX)
  colnames(yxOut) <- c("Y", paste("X", 1:(length(Mu_YX)-1), sep = ""))
  return(yxOut)
}

Compute_Y_Given_X <- function(Mu_YX, SigMat_YX)
{
  Mu_Y <- Mu_YX[1]
  Mu_X <- matrix(Mu_YX[-1], ncol = 1)
  
  SigYY <- SigMat_YX[1,1]
  SigXX <- SigMat_YX[-1,-1]
  SigYX <- matrix(SigMat_YX[1,-1], nrow = 1)
  SigXY <- matrix(SigMat_YX[-1,1], ncol = 1)
  
  Beta_0 <- Mu_Y-SigYX%*%solve(SigXX)%*%Mu_X
  Beta_1 <- SigYX%*%solve(SigXX)
  Sd <- SigYY-SigYX%*%solve(SigXX)%*%SigXY
  
  return(list(Beta0=Beta_0,
              Beta1=t(Beta_1),
              VarYX = Sd))
}

Compute_X_Given_Y <- function(Mu_YX, SigMat_YX)
{
  Mu_Y <- Mu_YX[1]
  Mu_X <- matrix(Mu_YX[-1], ncol = 1)
  
  SigYY <- SigMat_YX[1,1]
  SigXX <- SigMat_YX[-1,-1]
  SigYX <- matrix(SigMat_YX[1,-1], nrow = 1)
  SigXY <- matrix(SigMat_YX[-1,1], ncol = 1)
  
  Beta_0 <- Mu_X-SigXY%*%solve(SigYY)%*%Mu_Y
  Beta_1 <- SigXY%*%solve(SigYY)
  SdMat <- SigXX-SigXY%*%solve(SigYY)%*%SigYX
  
  return(list(Beta0=t(Beta_0),
              Beta1=t(Beta_1),
              VarXY=SdMat))
}

Generate_Y_Marginal_Target <- function(m, Mu_Y_T, Sig_Y_T)
{
  rnorm(m, Mu_Y_T, Sig_Y_T)
}

Compute_Rho_Parameters <- function(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)
{
  Mu_Y_S <- Mu_YX[1]
  Var_Y_S <- SigMat_YX[1,1]
  Var_Y_T <- Sig_Y_T^2
  
  Beta_Y <- Mu_Y_T/Var_Y_T-Mu_Y_S/Var_Y_S
  Beta_Y2 <- (1/Var_Y_S-1/Var_Y_T)/2
  
  return(c(Beta_Y, Beta_Y2))
}

Generate_X_Given_Y <- function(YVec, X_Given_Y_List)
{
  n <- length(YVec)
  Beta0 <- X_Given_Y_List$Beta0
  Beta1 <- X_Given_Y_List$Beta1
  xDim <- length(Beta0)
  yMat <- matrix(YVec, ncol = xDim, nrow = n)
  
  Beta0 <- matrix(c(Beta0), ncol = xDim, nrow = n, byrow = T)
  Beta1 <- matrix(c(Beta1), ncol = xDim, nrow = n, byrow = T)
  
  Mu_X_Given_Y <- Beta0+Beta1*yMat
  SdMat <- X_Given_Y_List$VarXY
  
  outMat <- matrix(0, nrow = n, ncol = xDim)
  for (i in 1:n)
  {
    tempMu <- Mu_X_Given_Y[i,]
    outMat[i,] <- MASS::mvrnorm(1, tempMu, SdMat)
  }
  colnames(outMat) <- paste("X", 1:xDim, sep = "")
  
  return(outMat)
}

Fit_YX_Model_Source <- function(xyDat)
{
  lFit <- lm(Y~., data = as.data.frame(xyDat))
  SdYX <- sigma(lFit)
  VarYX <- SdYX^2
  Beta0 <- coef(lFit)[1]
  Beta1 <- coef(lFit)[-1]
  
  return(list(Beta0 = matrix(Beta0, ncol = 1, nrow = 1),
              Beta1 = matrix(Beta1, ncol = 1, nrow = length(Beta1)),
              VarYX = matrix(VarYX, ncol = 1, nrow = 1)))
}

Generate_Y_Given_X <- function(YX_List, Xt, xyDat)
{
  m <- nrow(Xt)
  
  Beta0 <- c(YX_List$Beta0)
  Beta1 <- YX_List$Beta1
  MuVec <- Beta0+Xt%*%Beta1
  VarYX <- YX_List$VarYX
  
  addTerm <- matrix(rnorm(m, 0, sqrt(VarYX)), ncol = 1)
  yCol <- MuVec+addTerm
  colnames(yCol) <- "Y"
  YXOut <- cbind(yCol, Xt)
  YXOut <- rbind(xyDat, YXOut)
  
  return(YXOut)
}

Generate_Dat <- function(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
{
  dat <- Generate_Y_X_Marginal(n = n, Mu_YX = Mu_YX, SigMat_YX = SigMat_YX)
  XY_List <- Compute_X_Given_Y(Mu_YX = Mu_YX, SigMat_YX = SigMat_YX)
  
  yTarget <- Generate_Y_Marginal_Target(m = m, Mu_Y_T = Mu_Y_T, Sig_Y_T = Sig_Y_T)
  xTarget <- Generate_X_Given_Y(yTarget, XY_List)
  
  return(list(sDat = dat,
              tDat = xTarget))
}

