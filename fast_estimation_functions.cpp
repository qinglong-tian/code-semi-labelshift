#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector E_S_RHO_X_CPP(NumericVector beta_rho, double rho_pwr,
                            NumericMatrix x_mat_no_intercept, NumericVector coef_y_x_s, double sigma_y_x_s,
                            NumericVector xList, NumericVector wList) {
  int i, j;
  int num_row = x_mat_no_intercept.nrow(), num_col = x_mat_no_intercept.ncol(), ghpoints = xList.length();
  NumericVector mean_y_x_s(num_row), e_s_rho_x(num_row);
  
  for (i=0; i<num_row; i++) {
    mean_y_x_s(i) = coef_y_x_s(0);
    for (j=0; j<num_col; j++) {
      mean_y_x_s(i) += coef_y_x_s(j+1)*x_mat_no_intercept(i, j);
    }//
  }
  
  double yTemp;
  for (i=0; i<num_row; i++) {
    yTemp = 0;
    for (j=0; j<ghpoints; j++) {
      yTemp = sqrt(2)*sigma_y_x_s*xList(j)+mean_y_x_s(i);
      e_s_rho_x(i) += exp(rho_pwr*yTemp*beta_rho(0)+rho_pwr*yTemp*yTemp*beta_rho(1))*wList(j)/sqrt(3.1415926);
    }
  }
  
  return e_s_rho_x;
}

// [[Rcpp::export]]
double E_S_RHO_CPP(NumericVector beta_rho, bool ispar, List parameters) {
  double e_s_rho=0;
  
  if (ispar) {
    double Mu_Y_S = parameters["mu"], Sig_Y_S = parameters["sigma"];
    NumericVector xList = parameters["xList"], wList = parameters["wList"];
    int ghpoint = xList.length();
    
    double yTemp;
    for (int i=0; i<ghpoint; i++) {
      yTemp = sqrt(2)*Sig_Y_S*xList(i)+Mu_Y_S;
      e_s_rho += exp(yTemp*beta_rho(0)+yTemp*yTemp*beta_rho(1))*wList(i);
    }
    e_s_rho /= sqrt(3.1415926);
  }
  else {
    NumericVector y_vec = parameters[0];
    int yLen = y_vec.length();
    
    for (int i=0; i<yLen; i++) {
      e_s_rho += exp(beta_rho(0)*y_vec(i)+beta_rho(1)*y_vec(i)*y_vec(i));
    }
    e_s_rho /= yLen;
  }
  
  return e_s_rho;
}

// [[Rcpp::export]]
NumericVector COMPUTE_TAU_CPP(NumericVector e_s_rho_x, NumericVector e_s_rho2_x, double c_ps, double piVal) {
  int num_of_x = e_s_rho_x.length();
  double tmp;
  NumericVector tauX(num_of_x);
  for (int i=0; i<num_of_x; i++) {
    tmp = e_s_rho2_x(i)/e_s_rho_x(i)/c_ps/piVal;
    tauX(i) = tmp/(tmp+1/(1-piVal));
  }
  
  return tauX;
}

// [[Rcpp::export]]
double E_T_TAU_CPP(NumericVector e_s_rho_x, NumericVector e_s_rho2_x, double c_ps, double piVal) {
  NumericVector tau_vec = COMPUTE_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal);
  double mean_tau = mean(tau_vec);
  
  return mean_tau;
}

// [[Rcpp::export]]
NumericVector E_T_D_LOG_RHO_DIV_D_BETA_CPP(NumericVector beta_rho, bool ispar, List parameters, double c_ps) {
  
  NumericVector outBeta(2);
  
  if (ispar) {
    double Mu_Y_S = parameters["mu"], Sig_Y_S = parameters["sigma"];
    NumericVector xList = parameters["xList"], wList = parameters["wList"];
    int ghnum = xList.length();
    
    double yTemp, rhoTmp;
    for (int i=0; i<ghnum; i++) {
      yTemp = sqrt(2)*Sig_Y_S*xList(i)+Mu_Y_S;
      rhoTmp = exp(yTemp*beta_rho(0)+yTemp*yTemp*beta_rho(1));
      outBeta(0) += wList(i)*rhoTmp*yTemp;
      outBeta(1) += wList(i)*rhoTmp*yTemp*yTemp;
    }
    outBeta(0) /= (c_ps*sqrt(3.1415926));
    outBeta(1) /= (c_ps*sqrt(3.1415926));
  }
  else
  {
    NumericVector y_s_external = parameters[0];
    int y_num = y_s_external.length();
    double yCurrent, rhoCurrent;
    for (int i=0; i<y_num; i++) {
      
      yCurrent = y_s_external(i);
      rhoCurrent = exp(yCurrent*beta_rho(0)+yCurrent*yCurrent*beta_rho(1));
      
      outBeta(0) += rhoCurrent*yCurrent;
      outBeta(1) += rhoCurrent*yCurrent*yCurrent;
      
    }
    outBeta(0) /= (y_num*c_ps);
    outBeta(1) /= (y_num*c_ps);
  }
  
  return outBeta;
}

// [[Rcpp::export]]
NumericMatrix E_T_D_LOG_RHO_DIV_D_BETA_X_CPP(NumericVector beta_rho, NumericMatrix x_mat_no_intercept,
                                             NumericVector coef_y_x_s, double sigma_y_x_s, NumericVector e_s_rho_x,
                                             NumericVector xList, NumericVector wList) {
  int num_of_x = x_mat_no_intercept.nrow(), num_of_covariate = x_mat_no_intercept.ncol(), ghnum = xList.length();
  int i, j;
  NumericMatrix out_mat(num_of_x, 2);
  NumericVector y_x_mean(num_of_x);
  
  for (i=0; i<num_of_x; i++) {
    y_x_mean(i) = coef_y_x_s(0);
    for (j=0; j<num_of_covariate; j++) {
      y_x_mean(i) += coef_y_x_s(j+1)*x_mat_no_intercept(i,j);
    }
  }
  
  double yTmp, rhoSrc, tmp1, tmp2;
  for (i=0; i<num_of_x; i++) {
    tmp1 = 0;
    tmp2 = 0;
    for (j=0; j<ghnum; j++) {
      yTmp = sqrt(2)*sigma_y_x_s*xList(j)+y_x_mean(i);
      rhoSrc = exp(yTmp*beta_rho(0)+yTmp*yTmp*beta_rho(1));
      tmp1 += wList(j)*rhoSrc*yTmp;
      tmp2 += wList(j)*rhoSrc*yTmp*yTmp;
    }
    tmp1 /= (sqrt(3.1415926)*e_s_rho_x(i));
    tmp2 /= (sqrt(3.1415926)*e_s_rho_x(i));
    
    out_mat(i,0) = tmp1;
    out_mat(i,1) = tmp2;
  }
  
  return out_mat;
}

// [[Rcpp::export]]
NumericMatrix Compute_S_CPP(NumericVector beta_rho, NumericMatrix x_mat_no_intercept,
                            NumericVector coef_y_x_s, double sigma_y_x_s, NumericVector e_s_rho_x,
                            NumericVector xList, NumericVector wList,
                            bool ispar, List parameters, double c_ps) {
  NumericMatrix e_t_d_log_rho_div_d_beta_x = E_T_D_LOG_RHO_DIV_D_BETA_X_CPP(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s,
                                                                            e_s_rho_x, xList, wList);
  NumericVector e_t_d_log_rho_div_d_beta = E_T_D_LOG_RHO_DIV_D_BETA_CPP(beta_rho, ispar, parameters, c_ps);
  int num_of_x = e_t_d_log_rho_div_d_beta_x.nrow();
  for (int i=0; i< num_of_x; i++) {
    for (int j=0; j<2; j++) {
      e_t_d_log_rho_div_d_beta_x(i,j) -= e_t_d_log_rho_div_d_beta(j);
    }
  }
  
  return e_t_d_log_rho_div_d_beta_x;
}

// [[Rcpp::export]]
NumericMatrix ComputeEfficientScore_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                                        double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                                        bool ispar, List parameters, NumericVector xList, NumericVector wList) {
  int num_of_source = sData.nrow(), num_of_target = tData.nrow(), num_of_cov = tData.ncol();
  int num_of_total = num_of_source+num_of_target;
  
  NumericVector y_internal = sData( _ , 0);
  NumericVector rho_internal(num_of_source);
  for (int i=0; i<num_of_source; i++) {
    rho_internal(i) = exp(y_internal(i)*beta_rho(0)+y_internal(i)*y_internal(i)*beta_rho(1));
  }
  double c_ps = E_S_RHO_CPP(beta_rho, ispar, parameters);
  
  NumericMatrix xMatAll(num_of_total, num_of_cov);
  for (int i=0; i<num_of_total; i++) {
    for (int j=0; j<num_of_cov; j++) {
      if (i < num_of_source)
      {
        xMatAll(i,j) = sData(i, j+1);
      }
      else
      {
        xMatAll(i,j) = tData(i-num_of_source, j);
      }
    }
  }
  
  NumericVector multiplier(num_of_total);
  for(int i=0; i<num_of_total; i++) {
    if(i<num_of_source)
    {
      multiplier(i) = rho_internal(i)/c_ps/piVal;
    }
    else
    {
      multiplier(i) = -1/(1-piVal);
    }
  }
  
  NumericVector e_s_rho2_x = E_S_RHO_X_CPP(beta_rho, 2, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector e_s_rho_x = E_S_RHO_X_CPP(beta_rho, 1, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector tau_x_internal = COMPUTE_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal);
  NumericMatrix s_x_internal = Compute_S_CPP(beta_rho, xMatAll, coef_y_x_s, sigma_y_x_s, e_s_rho_x, xList, wList,
                                             ispar, parameters, c_ps);

  NumericVector e_s_rho_x_ext = E_S_RHO_X_CPP(beta_rho, 1, tDat_ext, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector e_s_rho2_x_ext = E_S_RHO_X_CPP(beta_rho, 2, tDat_ext, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector tau_x_external = COMPUTE_TAU_CPP(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal);
  NumericMatrix s_x_external = Compute_S_CPP(beta_rho, tDat_ext, coef_y_x_s, sigma_y_x_s, e_s_rho_x_ext, xList, wList,
                                             ispar, parameters, c_ps);

  double e_t_tau = mean(tau_x_external);
  int num_of_ext = tDat_ext.nrow();
  NumericVector e_t_tau_s(2);
  for(int i=0; i<num_of_ext; i++) {
    e_t_tau_s(0) += tau_x_external(i)*s_x_external(i,0);
    e_t_tau_s(1) += tau_x_external(i)*s_x_external(i,1);
  }
  e_t_tau_s(0) /= num_of_ext;
  e_t_tau_s(1) /= num_of_ext;

  double b2ndpart1, b2ndpart2;
  NumericMatrix b1_x(num_of_total, 2);
  for(int i=0; i<num_of_total; i++) {
    b2ndpart1 = s_x_internal(i,0)-e_t_tau_s(0)/(e_t_tau-1);
    b2ndpart2 = s_x_internal(i,1)-e_t_tau_s(1)/(e_t_tau-1);

    b1_x(i,0) = -multiplier(i)*(1-piVal)*(1-tau_x_internal(i))*b2ndpart1;
    b1_x(i,1) = -multiplier(i)*(1-piVal)*(1-tau_x_internal(i))*b2ndpart2;
  }
  
  return b1_x;
}

// [[Rcpp::export]]
double EstimateBetaFunc_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                            double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                            bool ispar, List parameters, NumericVector xList, NumericVector wList, bool weights = false) {
  double out;
  NumericMatrix Seff = ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s,
                                                 ispar, parameters, xList, wList);
  int num_of_row = Seff.nrow();
  
  NumericVector rexpVec(num_of_row);
  if (!weights) {
    for (int i=0; i<num_of_row; i++) {
      rexpVec(i) = 1;
    }
  }
  else {
    for (int i=0; i<num_of_row; i++) {
      rexpVec(i) = R::rexp(1);
    }
  }
  
  double sumTmp1 = 0, sumTmp2 = 0;
  for(int i = 0; i<num_of_row; i++) {
    sumTmp1 += rexpVec(i)*Seff(i,0);
    sumTmp2 += rexpVec(i)*Seff(i,1);
  }
  sumTmp1 /= num_of_row;
  sumTmp2 /= num_of_row;
  
  out = sumTmp1*sumTmp1+sumTmp2*sumTmp2;

  return out;
}

// For the perturbation
// [[Rcpp::export]]
double EstimateBetaFuncPert_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                               double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                               bool ispar, List parameters, NumericVector xList, NumericVector wList, NumericVector rexpVec) {
  double out;
  NumericMatrix Seff = ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s,
                                                 ispar, parameters, xList, wList);
  int num_of_row = Seff.nrow();
  double sumTmp1 = 0, sumTmp2 = 0;
  
  for(int i = 0; i<num_of_row; i++) {
    sumTmp1 += rexpVec(i)*Seff(i,0);
    sumTmp2 += rexpVec(i)*Seff(i,1);
  }
  
  out = sumTmp1*sumTmp1+sumTmp2*sumTmp2;
  
  return out;
}

// [[Rcpp::export]]
NumericVector EstimateBetaVarFunc_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                                  double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                                  bool ispar, List parameters, NumericVector xList, NumericVector wList) {
  NumericMatrix Seff = ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s,
                                                 ispar, parameters, xList, wList);
  NumericMatrix tmpMat(2,2);
  NumericVector outVec(2);
  
  int num_of_row = Seff.nrow();
  for(int i=0; i<num_of_row; i++) {
    tmpMat(0,0) += Seff(i,0)*Seff(i,0);
    tmpMat(0,1) += Seff(i,0)*Seff(i,1);
    tmpMat(1,0) += Seff(i,0)*Seff(i,1);
    tmpMat(1,1) += Seff(i,1)*Seff(i,1);
  }
  
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      tmpMat(i,j) /= num_of_row;
    }
  }
  
  double det = tmpMat(0,0)*tmpMat(1,1)-tmpMat(1,0)*tmpMat(0,1);
  outVec(0) = sqrt(tmpMat(1,1)/det/num_of_row);
  outVec(1) = sqrt(tmpMat(0,0)/det/num_of_row);

  return outVec;
}

// [[Rcpp::export]]
NumericMatrix EstimateBetaCovMat_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                                     double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                                     bool ispar, List parameters, NumericVector xList, NumericVector wList) {
  NumericMatrix Seff = ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s,
                                                 ispar, parameters, xList, wList);
  NumericMatrix tmpMat(2,2);
  NumericMatrix outVec(2,2);
  
  int num_of_row = Seff.nrow();
  for(int i=0; i<num_of_row; i++) {
    tmpMat(0,0) += Seff(i,0)*Seff(i,0);
    tmpMat(0,1) += Seff(i,0)*Seff(i,1);
    tmpMat(1,0) += Seff(i,0)*Seff(i,1);
    tmpMat(1,1) += Seff(i,1)*Seff(i,1);
  }
  
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      tmpMat(i,j) /= num_of_row;
    }
  }
  
  double det = tmpMat(0,0)*tmpMat(1,1)-tmpMat(1,0)*tmpMat(0,1);
  outVec(0,0) = tmpMat(1,1)/det;
  outVec(1,1) = tmpMat(0,0)/det;
  outVec(0,1) = -tmpMat(1,0)/det;
  outVec(1,0) = -tmpMat(0,1)/det;
  
  return outVec;
}

// [[Rcpp::export]]

NumericVector COMPUTE_B_CPP(NumericVector tau_for_x, double e_t_tau, double c_ps, double piVal) {
  int num_of_x = tau_for_x.length();
  NumericVector Bx(num_of_x);
  for (int i=0; i<num_of_x; i++) {
    Bx(i) = (e_t_tau-tau_for_x(i))/(1-e_t_tau);
  }
  return Bx;
}

// [[Rcpp::export]]

NumericVector E_S_RHO2_PSI_X_CPP(NumericVector beta_rho, NumericMatrix x_mat_no_intercept, NumericVector coef_y_x_s,
                                 double sigma_y_x_s, NumericVector xList, NumericVector wList) {
  int num_of_x = x_mat_no_intercept.nrow();
  int num_of_p = x_mat_no_intercept.ncol();

  NumericVector mean_y_x_s(num_of_x);
  for(int i=0; i<num_of_x; i++) {
    mean_y_x_s(i) = coef_y_x_s(0);
    for(int j=0; j<num_of_p; j++) {
      mean_y_x_s(i) += x_mat_no_intercept(i,j)*coef_y_x_s(j+1);
    }
  }

  int num_of_gh = xList.length();
  double yTmp;
  NumericVector e_s_rho2_psi_x(num_of_x);
  for(int i=0; i<num_of_x; i++){
    for(int j=0; j<num_of_gh; j++){
      yTmp = mean_y_x_s(i)+sqrt(2)*sigma_y_x_s*xList(j);
      e_s_rho2_psi_x(i) += exp(2*yTmp*beta_rho(0)+2*yTmp*yTmp*beta_rho(1))*yTmp*wList(j);
    }
    e_s_rho2_psi_x(i) /=sqrt(3.1415926);
  }

  return e_s_rho2_psi_x;
}

// [[Rcpp::export]]

NumericVector COMPUTE_A_CPP(NumericVector beta_rho, double e_t_tau,
                            NumericVector tau_x_internal, NumericVector tau_x_external,
                            NumericVector e_s_rho2_psi_x_internal, NumericVector e_s_rho2_x_internal,
                            NumericVector e_s_rho2_psi_x_external, NumericVector e_s_rho2_x_external) {
  int num_of_x = tau_x_internal.length();

  NumericVector thirdPart(num_of_x);
  for(int i=0; i<num_of_x; i++) {
    thirdPart(i) = tau_x_internal(i)*e_s_rho2_psi_x_internal(i)/e_s_rho2_x_internal(i);
  }

  int num_of_ext_x = tau_x_external.length();
  double secondPart=0;
  for(int i=0; i<num_of_ext_x; i++) {
    secondPart += tau_x_external(i)*e_s_rho2_psi_x_external(i)/e_s_rho2_x_external(i);
  }
  secondPart /= num_of_ext_x;

  NumericVector firstPart(num_of_x);
  for(int i=0; i<num_of_x; i++) {
    firstPart(i) = (1-tau_x_internal(i))/(1-e_t_tau)*secondPart;
  }

  NumericVector AVec(num_of_x);
  for(int i=0; i<num_of_x; i++) {
    AVec(i) = firstPart(i)-thirdPart(i);
  }

  return AVec;
}

// [[Rcpp::export]]

NumericVector COMPUTE_H_CPP(NumericVector beta_rho, double theta, double c_ps, NumericMatrix sData, double e_t_tau, NumericVector tau_x_external, double piVal,
  NumericVector coef_y_x_s, double sigma_y_x_s, NumericVector e_s_rho2_psi_x_external, NumericVector e_s_rho2_x_external,
  NumericVector e_s_rho2_psi_x_internal, NumericVector e_s_rho2_x_internal, NumericVector tau_for_x) {
  NumericVector yVec = sData( _ , 0);
  int num_of_source_x = sData.nrow();
  NumericVector multiplier(num_of_source_x);
  double yTmp;
  for(int i=0; i<num_of_source_x; i++) {
    yTmp = yVec(i);
    multiplier(i) = 1/piVal/c_ps*exp(yTmp*beta_rho(0)+yTmp*yTmp*beta_rho(1));
  }

  NumericVector AVec = COMPUTE_A_CPP(beta_rho, e_t_tau, tau_for_x, tau_x_external, e_s_rho2_psi_x_internal, e_s_rho2_x_internal, e_s_rho2_psi_x_external, e_s_rho2_x_external);
  NumericVector BVec = COMPUTE_B_CPP(tau_for_x, e_t_tau, c_ps, piVal);

  NumericVector hVec(num_of_source_x);
  for (int i=0; i<num_of_source_x; i++) {
    hVec(i) = multiplier(i)*(yVec(i)-theta+AVec(i)-BVec(i)*theta);
  }

  return hVec;
}

// [[Rcpp::export]]

NumericVector E_S_RHO_Y_CPP(NumericVector beta_rho, bool ispar, List parameters, NumericMatrix sData) {
  NumericVector outVec(3);
  double yTmp;
  double rhoVal;

  if(ispar) {
    double mu_y_s = parameters["mu"], sigma_y_s = parameters["sigma"];
    NumericVector xList = parameters["xList"], wList = parameters["wList"];

    int ghnum = xList.length();
    for (int i=1; i<ghnum; i++) {
      yTmp = sqrt(2)*sigma_y_s*xList(i)+mu_y_s;
      rhoVal = exp(yTmp*beta_rho(0)+yTmp*yTmp*beta_rho(1));
      outVec(0) += wList(i)*rhoVal*yTmp;
      outVec(1) += wList(i)*rhoVal*yTmp*yTmp;
      outVec(2) += wList(i)*rhoVal*yTmp*yTmp*yTmp;
    }
    outVec(0) /= sqrt(3.1415926);
    outVec(1) /= sqrt(3.1415926);
    outVec(2) /= sqrt(3.1415926);
  }
  else {
    NumericVector yVec = sData( _ , 0);
    int num_of_source = yVec.length();
    for (int i=0; i<num_of_source; i++) {
      yTmp = yVec(i);
      rhoVal = exp(yTmp*beta_rho(0)+yTmp*yTmp*beta_rho(1));
      outVec(0) += rhoVal*yTmp;
      outVec(1) += rhoVal*yTmp*yTmp;
      outVec(2) += rhoVal*yTmp*yTmp*yTmp;
    }
    outVec(0) /= num_of_source;
    outVec(1) /= num_of_source;
    outVec(2) /= num_of_source;
  }

  return outVec;
}

// [[Rcpp::export]]

NumericVector COMPUTE_D_MULTIPLIER(NumericVector beta_rho, bool ispar, List parameters, NumericMatrix sData,
                                   double theta, double piVal, double c_ps,
                                   double e_t_tau, NumericVector tau_x_external, NumericVector coef_y_x_s,
                                   double sigma_y_x_s, NumericVector e_s_rho2_psi_x_external,
                                   NumericVector e_s_rho2_x_external,
                                   NumericVector e_s_rho2_psi_x_internal_source, NumericVector e_s_rho2_x_internal_source, NumericVector tau_for_x_source,
                                   NumericVector e_s_rho_x_source, NumericVector xList, NumericVector wList) {
  NumericVector e_s_rho_y_123 = E_S_RHO_Y_CPP(beta_rho, ispar, parameters, sData);
  NumericVector DMult(2);
  
  double e_t_rho_y_1 = e_s_rho_y_123(0)/c_ps, e_t_rho_y_2 = e_s_rho_y_123(1)/c_ps, e_t_rho_y_3 = e_s_rho_y_123(2)/c_ps;
  DMult(0) = e_t_rho_y_2-theta*e_t_rho_y_1;
  DMult(1) = e_t_rho_y_3-theta*e_t_rho_y_2;

  int num_of_xy = sData.ncol();
  NumericMatrix sData_x_only = sData( _ , Range(1, num_of_xy-1));
  NumericVector HVec = COMPUTE_H_CPP(beta_rho, theta, c_ps, sData, e_t_tau, tau_x_external, piVal, coef_y_x_s, sigma_y_x_s, e_s_rho2_psi_x_external, e_s_rho2_x_external, e_s_rho2_psi_x_internal_source, e_s_rho2_x_internal_source, tau_for_x_source);
  NumericMatrix SVec = Compute_S_CPP(beta_rho, sData_x_only, coef_y_x_s, sigma_y_x_s, e_s_rho_x_source, xList, wList, ispar, parameters, c_ps);

  int num_of_x_source = sData.nrow();
  NumericVector rhoVal(num_of_x_source);
  double yTmp;
  for (int i=0; i<num_of_x_source; i++) {
    yTmp = sData(i,0);
    rhoVal(i) = exp(yTmp*beta_rho(0)+yTmp*yTmp*beta_rho(1));
  }

  NumericVector HS_Mean(2);
  for(int i=0; i<num_of_x_source; i++) {
    HS_Mean(0) += rhoVal(i)*HVec(i)*SVec(i,0);
    HS_Mean(1) += rhoVal(i)*HVec(i)*SVec(i,1);
  }
  HS_Mean(0) /= num_of_x_source*c_ps;
  HS_Mean(1) /= num_of_x_source*c_ps;

  DMult(0) -= (1-piVal)*HS_Mean(0);
  DMult(1) -= (1-piVal)*HS_Mean(1);

  return DMult;
}

// [[Rcpp::export]]

NumericVector COMPUTE_D_CPP(NumericVector beta_rho, bool ispar, List parameters, NumericMatrix sData,
                                   double theta, double piVal, double c_ps,
                                   double e_t_tau, NumericVector tau_x_external, NumericVector coef_y_x_s,
                                   double sigma_y_x_s, NumericVector e_s_rho2_psi_x_external,
                                   NumericVector e_s_rho2_x_external,
                                   NumericVector e_s_rho2_psi_x_internal_source, NumericVector e_s_rho2_x_internal_source, NumericVector tau_for_x_source,
                                   NumericVector e_s_rho_x_source, NumericVector xList, NumericVector wList,
                                   NumericMatrix MatInv) {
  NumericVector DMult = COMPUTE_D_MULTIPLIER(beta_rho, ispar, parameters, sData,
                                             theta, piVal, c_ps,
                                             e_t_tau, tau_x_external, coef_y_x_s,
                                             sigma_y_x_s, e_s_rho2_psi_x_external,
                                             e_s_rho2_x_external, e_s_rho2_psi_x_internal_source,
                                             e_s_rho2_x_internal_source, tau_for_x_source,
                                             e_s_rho_x_source, xList, wList);
  NumericVector DVec(2);
  DVec(0) = DMult(0)*MatInv(0,0)+DMult(1)*MatInv(1,0);
  DVec(1) = DMult(0)*MatInv(1,0)+DMult(1)*MatInv(1,1);

  return DVec;
}

// [[Rcpp::export]]

NumericVector COMPUTE_EFFICIENT_IF_FOR_THETA_CPP(double theta,
                                                 int num_of_target,
                                                 double piVal, NumericVector rhoValSource,
                                                 double c_ps, NumericVector yVec,
                                                 NumericVector beta_rho, double e_t_tau,
                                                 NumericVector tau_x_internal_all,
                                                 NumericVector tau_x_external,
                                                 NumericVector e_s_rho2_psi_x_internal_all,
                                                 NumericVector e_s_rho2_x_internal_all,
                                                 NumericVector e_s_rho2_psi_x_external,
                                                 NumericVector e_s_rho2_x_external,
                                                 NumericMatrix MatInv,
                                                 bool ispar, List parameters,
                                                 NumericMatrix sData,
                                                 NumericVector coef_y_x_s, double sigma_y_x_s,
                                                 NumericVector e_s_rho2_psi_x_internal_source,
                                                 NumericVector e_s_rho2_x_internal_source,
                                                 NumericVector tau_for_x_source,
                                                 NumericVector e_s_rho_x_source,
                                                 NumericVector xList, NumericVector wList,
                                                 NumericMatrix SEff) {
  int num_of_source =  yVec.length();
  NumericVector firstPart(num_of_source+num_of_target);
  for(int i=0; i<num_of_source; i++) {
    firstPart(i) = 1/piVal*rhoValSource(i)/c_ps*(yVec(i)-theta);
  }
  
  NumericVector secondPart(num_of_source+num_of_target);
  NumericVector AVec = COMPUTE_A_CPP(beta_rho, e_t_tau,
                                     tau_x_internal_all, tau_x_external,
                                     e_s_rho2_psi_x_internal_all, e_s_rho2_x_internal_all,
                                     e_s_rho2_psi_x_external, e_s_rho2_x_external);
  NumericVector BVec = COMPUTE_B_CPP(tau_x_internal_all, e_t_tau, c_ps, piVal);
  for(int i=0; i<num_of_source; i++) {
    secondPart(i) = rhoValSource(i)/c_ps/piVal*(AVec(i)-theta*BVec(i));
  }
  for(int i=0; i<num_of_target; i++) {
    secondPart(num_of_source+i) = -1/(1-piVal)*(AVec(num_of_source+i)-theta*BVec(num_of_source+i));
  }

  NumericVector thirdPart(num_of_source+num_of_target);
  NumericVector DMat = COMPUTE_D_CPP(beta_rho, ispar, parameters, sData,
                                     theta, piVal, c_ps, e_t_tau,
                                     tau_x_external, coef_y_x_s, sigma_y_x_s,
                                     e_s_rho2_psi_x_external, e_s_rho2_x_external,
                                     e_s_rho2_psi_x_internal_source, e_s_rho2_x_internal_source, tau_for_x_source,
                                     e_s_rho_x_source, xList, wList, MatInv);
  for(int i=0; i<num_of_source+num_of_target; i++) {
    thirdPart(i) = SEff(i,0)*DMat(0)+SEff(i,1)*DMat(1);
  }

  NumericVector IF(num_of_source+num_of_target);
  for(int i=0; i<num_of_source+num_of_target; i++) {
    IF(i) = firstPart(i)+secondPart(i)+thirdPart(i);
  }

  return IF;
}

// [[Rcpp::export]]
double COMPUTE_THETA_CPP(int num_of_target,
                        double piVal, NumericVector rhoValSource,
                        double c_ps, NumericVector yVec,
                        NumericVector beta_rho, double e_t_tau,
                        NumericVector tau_x_internal_all,
                        NumericVector tau_x_external,
                        NumericVector e_s_rho2_psi_x_internal_all,
                        NumericVector e_s_rho2_x_internal_all,
                        NumericVector e_s_rho2_psi_x_external,
                        NumericVector e_s_rho2_x_external) {
  int num_of_source = yVec.length(), i;
  double numerator = 0, denominator = 0;
  for (i=0; i<num_of_source; i++) {
    numerator += rhoValSource(i)/c_ps/piVal*yVec(i);
  }
  NumericVector secondPart(num_of_source+num_of_target);
  NumericVector AVec = COMPUTE_A_CPP(beta_rho, e_t_tau,
                                     tau_x_internal_all, tau_x_external,
                                     e_s_rho2_psi_x_internal_all, e_s_rho2_x_internal_all,
                                     e_s_rho2_psi_x_external, e_s_rho2_x_external);
  for (i=0; i<num_of_source; i++) {
    numerator += 1/piVal*rhoValSource(i)/c_ps*AVec(i);
  }
  for (i=0; i<num_of_target; i++) {
    numerator -= 1/(1-piVal)*AVec(i+num_of_source);
  }

  NumericVector BVec = COMPUTE_B_CPP(tau_x_internal_all, e_t_tau, c_ps, piVal);
  for (i=0; i<num_of_source; i++) {
    denominator += 1/piVal*rhoValSource(i)/c_ps;
    denominator += 1/piVal*rhoValSource(i)/c_ps*BVec(i);
  }
  for (i=0; i<num_of_target; i++) {
    denominator -= 1/(1-piVal)*BVec(i+num_of_source);
  }

  double Theta = numerator/denominator;

  return Theta;
}

// [[Rcpp::export]]

double COMPUTE_THETA_NAIVE(NumericVector betaHat, NumericVector rhoValSource,
                           NumericVector yVec, double c_ps) {
  double out = 0;
  int num_of_s = yVec.length();
  for (int i=0; i<num_of_s; i++) {
    out += yVec(i)*rhoValSource(i);
  }
  out /= num_of_s*c_ps;
  
  return out;
}

// [[Rcpp::export]]
double COMPUTE_THETA_Pert_CPP(int num_of_target,
                              double piVal,
                              NumericVector rhoValSource,
                              double c_ps,
                              NumericVector yVec,
                              NumericVector beta_rho,
                              double e_t_tau,
                              NumericVector tau_x_internal_all,
                              NumericVector tau_x_external,
                              NumericVector e_s_rho2_psi_x_internal_all,
                              NumericVector e_s_rho2_x_internal_all,
                              NumericVector e_s_rho2_psi_x_external,
                              NumericVector e_s_rho2_x_external,
                              NumericVector rexpVec) {
  int num_of_source = yVec.length(), i;
  double numerator = 0, denominator = 0;
  for (i=0; i<num_of_source; i++) {
    numerator += rhoValSource(i)/c_ps/piVal*yVec(i)*rexpVec(i);
  }
  
  NumericVector AVec = COMPUTE_A_CPP(beta_rho, e_t_tau,
                                     tau_x_internal_all, tau_x_external,
                                     e_s_rho2_psi_x_internal_all, e_s_rho2_x_internal_all,
                                     e_s_rho2_psi_x_external, e_s_rho2_x_external);
  for (i=0; i<num_of_source; i++) {
    numerator += 1/piVal*rhoValSource(i)/c_ps*AVec(i)*rexpVec(i);
  }
  for (i=0; i<num_of_target; i++) {
    numerator -= 1/(1-piVal)*AVec(i+num_of_source)*rexpVec(i+num_of_source);
  }
  
  NumericVector BVec = COMPUTE_B_CPP(tau_x_internal_all, e_t_tau, c_ps, piVal);
  for (i=0; i<num_of_source; i++) {
    denominator += 1/piVal*rhoValSource(i)/c_ps*rexpVec(i);
    denominator += 1/piVal*rhoValSource(i)/c_ps*BVec(i)*rexpVec(i);
  }
  for (i=0; i<num_of_target; i++) {
    denominator -= 1/(1-piVal)*BVec(i+num_of_source)*rexpVec(i+num_of_source);
  }
  
  double Theta = numerator/denominator;
  
  return Theta;
}

// [[Rcpp::export]]
NumericMatrix Compute_Naive_S(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData, double c_ps, double piVal)
{
  int num_of_source = sData.nrow();
  int num_of_target = tData.nrow();
  
  int i;
  double rho_y;
  NumericMatrix outMat(num_of_source+num_of_target, 2);
  
  for(i=0;i<num_of_source;i++)
  {
    rho_y = exp(beta_rho(0)*sData(i,0)+beta_rho(1)*sData(i,0)*sData(i,0));
    outMat(i,0) = rho_y/piVal/c_ps*sData(i,1);
    outMat(i,1) = rho_y/piVal/c_ps*sData(i,2);
  }
  
  for(i=0;i<num_of_target;i++)
  {
    outMat(i+num_of_source,0) = -1/(1-piVal)*tData(i,0);
    outMat(i+num_of_source,1) = -1/(1-piVal)*tData(i,1);
  }
  
  return outMat;
}

// [[Rcpp::export]]
double Estimate_Naive_Beta(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData, double piVal, bool ispar, List parameters)
{
  double c_ps = E_S_RHO_CPP(beta_rho, ispar, parameters);
  NumericMatrix Naive_Mat = Compute_Naive_S(beta_rho, sData, tData, c_ps, piVal);
  int num_of_total = Naive_Mat.nrow();
  double total1 = 0, total2 = 0;
  for(int i=0; i<num_of_total; i++)
  {
    total1 += Naive_Mat(i,0);
    total2 += Naive_Mat(i,1);
  }
  total1 /= num_of_total;
  total2 /= num_of_total;
  
  double out = total1*total1+total2*total2;
  
  return out;
}
