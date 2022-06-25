# Generate Binary Y Data

Generate_X_Joint <- function(n, Mu_X, Sigma_X)
{
  xOut <- MASS::mvrnorm(n, Mu_X, Sigma_X)
  colnames(xOut) <- paste("X", 1:length(Mu_X), sep = "")
  
  return(xOut)
}

Generate_Binary_Y_Given_X <- function(xMat, alphaVec, probit = F)
{
  xMat1 <- cbind(1, xMat)
  oddVec <- xMat1 %*% matrix(alphaVec, ncol = 1)
  if (!probit)
  {
    probVec <- 1 / (1 + exp(-oddVec))
  }
  else
  {
    probVec <- pnorm(oddVec)
  }
  Y <- as.numeric(runif(nrow(xMat1)) < probVec)
  YXMat <- cbind(Y, xMat)
  
  return(YXMat)
}

Generate_R_Given_Y <- function(yVec, gammaVec)
{
  yMat <- cbind(1, yVec)
  gammaVec <- matrix(gammaVec, ncol = 1)
  oddVec <- yMat %*% gammaVec
  probVec <- 1 / (1 + exp(-oddVec))
  RVec <- as.numeric(runif(length(yVec)) < probVec)
  
  return(RVec)
}

Generate_Binary_Data <- function(n, Mu_X, Sigma_X, alphaVec, gammaVec, probit = F)
{
  xMat <- Generate_X_Joint(n, Mu_X, Sigma_X)
  YXMat <- Generate_Binary_Y_Given_X(xMat, alphaVec, probit)
  
  yVec <- YXMat[,1]
  R <- Generate_R_Given_Y(yVec, gammaVec)
  
  YXMat1 <- YXMat[R == 1, ]
  YXMat0 <- YXMat[R == 0, ]
  
  return(list(sData = YXMat1, tData = YXMat0[,-1], tData_full = YXMat0))
}
