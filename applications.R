library(tidyverse)
library(reshape2)
source("application_func.R")
dat_dir <- "real_dat"
dat3 <- readr::read_csv(paste(dat_dir, "/data3.csv", sep = ""))

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
  ) %>% mutate(admittime = as.Date(admittime, "%m/%d/%Y")) %>% filter(age < 65,
                                                                      sofa <= 10) %>% group_by(subject_id) %>% slice(which.min(admittime)) %>% na.omit %>% mutate(R = if_else(insurance %in% c("Medicare", "Medicaid"), 0, 1))

sDat <- data %>% filter(R == 1)
tDat <- data %>% filter(R == 0)

rSofaLogistic <-
  glm(R ~ poly(sofa, 2),
      family = binomial(link = "logit"),
      data = data)
Compute_GLM_Pr_R1_Given_Y(sDat, tDat, rSofaLogistic) -> df
df <- df %>% mutate(SOFA = as.numeric(SOFA))
df %>% melt(
  measure.vars = c("rProb", "rProbPred"),
  variable.name = "Type",
  value.name = "Prob"
) -> df
df %>% ggplot(aes(x = SOFA, y = Prob, col = Type)) + geom_point() + geom_line()
