library(tidyverse)
source("summarize_theta_results.R")
# Bias
outBias <- out[, c(1, 2, 3, 5)] %>% as.data.frame
outBias <-
  reshape2::melt(
    outBias,
    id.vars = c("ratio", "n"),
    value.name = "Bias",
    measure.vars = c("BiasEff", "BiasNaive"),
    variable.name = "Method"
  )
outBias$Method <-
  factor(
    outBias$Method,
    levels = c("BiasEff", "BiasNaive"),
    labels = c("Efficient", "Naive")
  )
outBias %>% ggplot(aes(x = n, y = Bias)) + geom_line(aes(col = Method, linetype = Method)) +
  geom_point(aes(col = Method, shape = Method)) + facet_wrap(~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + geom_hline(aes(yintercept = 0), linetype = "dashed")

# SE
outSE <- out[, c(1, 2, 4, 6, 8, 10)] %>% as.data.frame
outSE <-
  reshape2::melt(
    outSE,
    id.vars = c("ratio", "n"),
    value.name = "SE",
    measure.vars = c("SEEff", "SENaive", "SDFormula", "SDPert"),
    variable.name = "Method"
  )
outSE$Method <-
  factor(
    outSE$Method,
    levels = c("SEEff", "SENaive", "SDFormula", "SDPert"),
    labels = c(
      "Proposed (Empirical)",
      "Naive (Empirical)",
      "Eff. Formula (Est.)",
      "Eff. Perturbation (Est.)"
    )
  )
outSE %>% ggplot(aes(x = n, y = SE)) + geom_line(aes(col = Method, linetype = Method)) +
  geom_point(aes(col = Method, shape = Method)) + facet_wrap(~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + ylab("Standard Error")

# Coverage: Formula vs. Perturbation
outCP <- out[, c(1, 2, 7, 9)] %>% as.data.frame
outCP <-
  reshape2::melt(
    outCP,
    id.vars = c("ratio", "n"),
    value.name = "CP",
    measure.vars = c("CPFormula", "CPPert"),
    variable.name = "Method"
  )
outCP$Method <-
  factor(
    outCP$Method,
    levels = c("CPFormula", "CPPert"),
    labels = c("Formula", "Perturbation")
  )
outCP %>% ggplot(aes(x = n, y = CP)) + geom_line(aes(col = Method, linetype = Method)) +
  geom_point(aes(col = Method, shape = Method)) + facet_wrap(~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + geom_hline(aes(yintercept = 0.95), linetype = "dashed") + ylab("Coverage Probability") +
  ylim(c(0.9, 1))

# MSE Ratio
out %>% mutate(MSE_Ratio = MSEEff / MSENaive) %>% ggplot(aes(x = n, y = MSE_Ratio)) +
  geom_line() + geom_point() + facet_wrap( ~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  )))+ylim(c(0,1))+geom_hline(yintercept = 1, linetype = "dashed")+ylab("MSE Ratio")
