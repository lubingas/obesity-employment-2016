##########################################################################
## 04_diagnostics.R
## Residual diagnostics for the all-covariate models:
##  1) Cox-Snell residuals: overall fit of the Cox model and of the log-normal PH model
##  2) martingale residuals against age: is age linear on the log-hazard scale?
##  3) deviance residuals against the linear predictor: outliers
##  4) Schoenfeld residuals: does the BMI effect change with time (PH assumption)?
## Needs output/02_parametric.rds and output/03_cox.rds.
## Run from the repository root:  source("R/04_diagnostics.R")
##########################################################################
rm(list=ls())
source("R/00_functions.R")
library(survival)   ## only to check our answers

# prepare the data
d <- readData()
talpha <- d$talpha
time <- d$time
event <- d$event
bmi <- d$bmi
age <- d$age
x <- makeX(covars, d)

cx <- readRDS("output/03_cox.rds")
pm <- readRDS("output/02_parametric.rds")


##########################################################################
# 1) residuals from our Cox fit, and from the log-normal PH fit
##########################################################################
xb.cox <- linpred(cx$cox.fit2$estimate, x)
res.cox <- phResiduals(xb.cox, cx$H0.cox, time, talpha, event)

# log-normal PH: same residuals, with the fitted parametric H0(t)
ln <- pm$lnormPH.fit2$estimate
xb.ln <- linpred(ln[-(1:2)], x)
H0.ln <- function(t) base$lognormal$H(ln[1:2], t)
res.ln <- phResiduals(xb.ln, H0.ln, time, talpha, event)

# check: martingale and deviance residuals match residuals() from coxph() (same ties)
coxph.fit2 <- coxph(Surv(talpha, time, event) ~ age + sex + black + other + wasmarried + nevermarried +
                      publicinsur + privateinsur + lowinc + midinc + highinc + bmi,
                    data = d, ties = "breslow")
cat("check martingale vs residuals(coxph):", max(abs(res.cox$martingale - resid(coxph.fit2, "martingale"))), "\n")
cat("check deviance vs residuals(coxph):  ", max(abs(res.cox$deviance - resid(coxph.fit2, "deviance"))), "\n")


##########################################################################
# 2) Cox-Snell plots: Nelson-Aalen cumulative hazard of the residuals
##########################################################################
# If the model fits, the residuals behave like censored draws from a unit exponential, so
# their estimated cumulative hazard follows the 45-degree line.
cs.cox <- coxSnellNA(res.cox$coxsnell, event)
cs.ln <- coxSnellNA(res.ln$coxsnell, event)
cs.cox.bmi <- lapply(0:1, function(g) coxSnellNA(res.cox$coxsnell[bmi==g], event[bmi==g]))

# how far from the line, up to the 99th percentile of the residuals
dev45 <- function(cs) { k <- cs$r <= quantile(cs$r, 0.99); max(abs(cs$H[k] - cs$r[k])) }
cat("max |H(r) - r| below the 99th percentile: Cox", round(dev45(cs.cox), 3), " log-normal", round(dev45(cs.ln), 3), "\n")


##########################################################################
# 3) Schoenfeld residuals and the test of proportional hazards
##########################################################################
sch <- schoenfeld(cx$cox.fit2$estimate, cx$cox.fit2$vcov, time, talpha, event, x)
print(round(cbind(sch$test, p.approx = sch$test.approx[, "p"]), 4))

# checks: raw residuals match residuals(coxph, "schoenfeld") (rows are exit times), and the
# score test matches cox.zph() with g(t) = log(t)
rs <- resid(coxph.fit2, type = "schoenfeld")
rs <- rs[order(as.numeric(rownames(rs))), ]
cat("check Schoenfeld vs residuals(coxph):", max(abs(sch$resid - rs)), "\n")
zph <- cox.zph(coxph.fit2, transform = "log", terms = F)
cat("check PH test vs cox.zph(): max |chisq diff| =", max(abs(sch$test[, "chisq"] - zph$table[, "chisq"])), "\n")


##########################################################################
# 4) martingale residuals against age, deviance residuals against the linear predictor
##########################################################################
# lowess smooth of the martingale residuals: flat if age enters linearly. iter = 0 turns off
# lowess's robustness weights, which would otherwise treat the (fewer, larger) residuals of
# people who left as outliers and pull the smooth down
mart.age <- lowess(age, res.cox$martingale, f = 0.5, iter = 0)
dev.lp <- data.frame(lp = xb.cox, deviance = res.cox$deviance, event = event)

saveRDS(list(res.cox = res.cox, res.ln = res.ln, cs.cox = cs.cox, cs.ln = cs.ln, cs.cox.bmi = cs.cox.bmi,
             sch = sch, mart.age = mart.age, dev.lp = dev.lp), "output/04_diagnostics.rds")
write.csv(data.frame(term = rownames(sch$test), sch$test, row.names = NULL), "output/ph_test_schoenfeld.csv", row.names = F)
