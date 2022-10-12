
icu_data <- read.csv("data3.csv", header=T)

###Jiwei: the code below generates data4.csv that was previously used for analyzing missing data

## select columns
# demo: 8, 9, 10, 11
# chart: 39, 42, 45, 48, 51, 54
# lab: 15, 18, 21, 24, 27, 30, 32, 36, 56
# metrics: 58, 59

ind_demo <- c(8, 9, 10)
ind_chart <- c(39, 42, 45, 48, 51, 54)
ind_lab <- c(15, 18, 20, 24, 27, 30, 32, 36, 56)
ind_metrics <- c(58, 59)

# only include "ind" columns
icu_v1 <- icu_data[, c(ind_demo, ind_chart, ind_lab, ind_metrics)]
# remove NAs
icu_v2 <- icu_v1[complete.cases(icu_v1[, -which(colnames(icu_v1) == "albumin_mean")]), ]
# age 25-40 and married
if_age <- icu_v2$age>=25 & icu_v2$age<=40
if_marital <- icu_v2$marital_status=="MARRIED"
icu_v3 <- icu_v2[if_age & if_marital, -which(colnames(icu_v2) == "marital_status")]

# scale "score" variables
icu_v3$albumin_mean <- scale(icu_v3$albumin_mean)
icu_v3$sapsii <- scale(icu_v3$sapsii)
icu_v3$sofa <- scale(icu_v3$sofa)

# factorize
icu_v3$gender <- as.numeric(icu_v3$gender == "F")