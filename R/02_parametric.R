##########################################################################
## 02_parametric.R
## Parametric models with left truncation:
##  1) six baselines with no covariates (the 2016 exponential, Weibull, log-normal and
##     gamma fits, plus Gompertz and log-logistic), compared by AIC
##  2) log-normal and Weibull fitted separately by BMI group, to overlay on the KM curves
##  3) log-normal proportional hazards models: BMI only and all covariates
## Run from the repository root:  source("R/02_parametric.R")
##########################################################################
rm(list=ls())
source("R/00_functions.R")

# prepare the data
d <- readData()
talpha <- d$talpha
time <- d$time
event <- d$event
bmi <- d$bmi

x1 <- makeX("bmi", d)    # BMI only
x <- makeX(covars, d)    # all covariates


##########################################################################
# 1) baselines without covariates
##########################################################################
fits0 <- lapply(names(base), function(dist) fitPH(dist, time, talpha, event))
names(fits0) <- names(base)

aic0 <- data.frame(dist = names(fits0),
                   npar = sapply(fits0, function(f) f$npar),
                   logLik = sapply(fits0, function(f) f$loglik),
                   AIC = sapply(fits0, function(f) -2*f$loglik + 2*f$npar),
                   convergence = sapply(fits0, function(f) f$convergence), row.names = NULL)
aic0 <- aic0[order(aic0$AIC), ]
print(aic0, digits = 7)

# check: the exponential MLE has a closed form, exits / time at risk
cat("check exponential rate:", exp(fits0$exponential$estimate), sum(event)/sum(time - talpha), "\n")


##########################################################################
# 2) log-normal and Weibull by BMI group (no PH constraint), for the KM overlays
##########################################################################
bygroup <- list()
for (g in 0:1) {
  s <- bmi == g
  for (dist in c("lognormal", "weibull")) {
    bygroup[[paste(dist, g, sep = ".")]] <- fitPH(dist, time[s], talpha[s], event[s])
  }
}
print(do.call(rbind, lapply(names(bygroup), function(n)
  data.frame(fit = n, par = paste(names(bygroup[[n]]$estimate), collapse = ", "),
             est = paste(round(bygroup[[n]]$estimate, 4), collapse = ", "), logLik = bygroup[[n]]$loglik))))


##########################################################################
# 3) log-normal proportional hazards models
##########################################################################
# log L_i = delta_i [log h0(t_i) + xb_i] - e^{xb_i} [H0(t_i) - H0(talpha_i)]
# which is the log of the 2016 likelihood
#   [f(t) S(t)^(e^xb - 1) e^xb]^delta [S(t)^e^xb]^(1 - delta) / S(talpha)^e^xb
lnormPH.fit1 <- fitPH("lognormal", time, talpha, event, x1)
lnormPH.fit2 <- fitPH("lognormal", time, talpha, event, x)

print(round(report(lnormPH.fit2), 5))
cat(sprintf("log-likelihood: BMI only %.5f, all covariates %.5f\n", lnormPH.fit1$loglik, lnormPH.fit2$loglik))

# check against the 2016 maxLik fit (Table 1, log-normal column)
lnorm2016 <- c(a = 0.57090, b = 1.28906, age = -0.01777, bmi = 0.18400, logLik = -3381.60085)
cat("check vs 2016:\n")
print(round(rbind(y2016 = lnorm2016,
                  y2026 = c(report(lnormPH.fit2)[c("a", "b", "Bage", "Bbmi"), "est"], lnormPH.fit2$loglik)), 5))

saveRDS(list(fits0 = fits0, aic0 = aic0, bygroup = bygroup,
             lnormPH.fit1 = lnormPH.fit1, lnormPH.fit2 = lnormPH.fit2), "output/02_parametric.rds")
write.csv(aic0, "output/parametric_baselines_aic.csv", row.names = F)
