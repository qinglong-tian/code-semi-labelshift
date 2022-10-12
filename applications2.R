source("estimation_functions_binary.R")

dat_dir <- "real_dat"
dat3 <-
  read.table(paste(dat_dir, "/data3.csv", sep = ""),
             sep = ",",
             header = T)
data <-
  dat3 %>% select(
    diasbp_min,
    glucose_min,
    resprate_max,
    sysbp_min,
    temp_min,
    temp_max,
    hematocrit_mean,
    hematocrit_min,
    platelets_mean,
    platelets_min,
    redbloodcell_mean,
    redbloodcell_min,
    spo2_min,
    urea_n_min,
    urea_n_max,
    urea_n_mean,
    insurance,
    age,
    subject_id,
    admittime,
    sofa
  ) %>% mutate(admittime = as.Date(admittime, "%m/%d/%Y")) %>%
  mutate(Y = ifelse(sofa >= 4, 1, 0), .before = diasbp_min) %>%
  filter(age < 65) %>% group_by(subject_id) %>% slice(which.min(admittime)) %>%
  na.omit() %>% mutate(R = if_else(insurance %in% c("Medicare", "Medicaid"), 0, 1))
data %>% ungroup(subject_id) %>% select(-age,-subject_id,-admittime,-insurance,-sofa) -> data
sDat <- data %>% filter(R == 1) %>% select(-R) %>% as.data.frame()
tDat <-
  data %>% filter(R == 0) %>% select(-R,-Y) %>% as.data.frame()

logitReg <- glm(Y ~ ., data = sDat, family = binomial)
rfProbs <- predict(logitReg, newdata = sDat, type = "response")
rfProbt <- predict(logitReg, newdata = tDat, type = "response")
piVal <- nrow(sDat) / nrow(data)

trueBetaVal <- glm(R ~ Y, data = data, family = binomial)

initBeta <- 0
optim(
  par = trueBeta,
  fn = Compute_S_Eff_Sum,
  method = "BFGS",
  sData = sDat,
  piVal = piVal,
  yFittedGivenX = rfProbs,
  yFittedGivenX_tDat = rfProbt
) -> optimOut
c_ps <- E_S_RHO_Binary(optimOut$par, sDat)
Compute_SE_Binary(
  betaVal = optimOut$par,
  sData = sDat,
  c_ps = c_ps,
  piVal = piVal,
  yFittedGivenX = rfProbs,
  yFittedGivenX_tDat = rfProbt
) -> se
Estimate_Beta_Lipton(rfProbs, rfProbt, sDat, initBeta)

B <- 100
rexpVec_List <- lapply(1:B, function(x) {
  rexp(nrow(data))
})
mclapply(rexpVec_List, function(rexpVec)
{
  Estimate_Beta_Lipton_Pert(rfProbs, rfProbt, sDat, initBeta, rexpVec)
},
mc.cores = detectCores()) -> beta_vector
sd(unlist(beta_vector))
