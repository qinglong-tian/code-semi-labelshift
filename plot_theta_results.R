library(tidyverse)
source("summarize_theta_results.R")
library(ggpubr)
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
    labels = c("Proposed", "Naive")
  )
outBias %>% ggplot(aes(x = n, y = Bias)) + geom_line(aes(col = Method, linetype = Method)) +
  geom_point(aes(col = Method, shape = Method)) + facet_wrap(~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + geom_hline(aes(yintercept = 0), linetype = "dashed") -> biasfig2

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
    levels = c("SENaive", "SEEff", "SDFormula", "SDPert"),
    labels = c(
      "Naive",
      "Proposed",
      "Eff. Formula (Est.)",
      "Eff. Perturbation (Est.)"
    )
  )


outSE %>% filter(
  Method %in% c("Proposed", "Naive")
) %>% ggplot(
  aes(x = n, y = SE)
) + geom_line(
  aes(col = Method, linetype = Method)
) + geom_point(
  aes(col = Method, shape = Method)
) + facet_wrap( ~ ratio, nrow = 1, labeller = labeller(ratio = c(
    "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
  ))) + ylab("Standard Deviation")

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

# MSE
out %>% select(n, ratio, MSEEff, MSENaive) %>% melt(
  id.vars = c("ratio", "n"),
  value.name = "MSE",
  measure.vars = c("MSEEff", "MSENaive"),
  variable.name = "Method"
) %>% mutate(Method = factor(
  Method,
  levels = c("MSEEff", "MSENaive"),
  labels = c("Proposed", "Naive")
)) %>% ggplot(aes(x = n, y = MSE, col = Method)) + geom_line(aes(linetype = Method)) + geom_point(aes(shape = Method)) + facet_wrap(~ ratio, nrow = 1, labeller = labeller(ratio = c(
  "0.5" = "m/n=0.5", "1" = "m/n=1", "1.5" = "m/n=1.5"
))) -> msefig2
ggarrange(biasfig2, msefig2, ncol = 1, labels = c("(a)", "(b)"))

# Table

out %>% select(ratio, n, BiasEff, SEEff, SDPert, CPPert) %>% arrange(ratio, n) -> latexTable2
