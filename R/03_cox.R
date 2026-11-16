##########################################################################
## 03_cox.R
## Cox proportional hazards models from our own partial likelihood
## (delayed entry, Breslow ties), maximised with optim, checked against coxph().
## Run from the repository root:  source("R/03_cox.R")
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

x1 <- makeX("bmi", d)    # BMI only
x <- makeX(covars, d)    # all covariates


##########################################################################
# 1) fit: BMI only and all covariates
##########################################################################
# coxPHScore() is the analytic gradient, so optim's Hessian (differences of the gradient)
# is accurate enough for the standard errors.
cox.fit1 <- fitML(coxPHLik, setNames(0, "bmi"), gr = coxPHScore,
                  time = time, talpha = talpha, event = event, x = x1)
cox.fit2 <- fitML(coxPHLik, setNames(rep(0, ncol(x)), colnames(x)), gr = coxPHScore,
                  time = time, talpha = talpha, event = event, x = x)

hr <- cbind(beta = cox.fit2$estimate, se = cox.fit2$se,
            HR = exp(cox.fit2$estimate),
            lo = exp(cox.fit2$estimate - 1.96*cox.fit2$se),
            hi = exp(cox.fit2$estimate + 1.96*cox.fit2$se))
print(round(hr, 5))
cat(sprintf("log partial likelihood: BMI only %.5f, all covariates %.5f\n", cox.fit1$loglik, cox.fit2$loglik))


##########################################################################
# 2) check against coxph() with the same data and ties
##########################################################################
coxph.fit2 <- coxph(Surv(talpha, time, event) ~ age + sex + black + other + wasmarried + nevermarried +
                      publicinsur + privateinsur + lowinc + midinc + highinc + bmi,
                    data = d, ties = "breslow")
cat("check vs coxph: max |beta diff| =", max(abs(cox.fit2$estimate - coef(coxph.fit2))),
    " max |se diff| =", max(abs(cox.fit2$se - sqrt(diag(vcov(coxph.fit2))))),
    " logLik", round(cox.fit2$loglik, 5), round(coxph.fit2$loglik[2], 5), "\n")

# Breslow baseline cumulative hazard at the estimates (used for the residuals in 04)
H0.cox <- coxBaseHaz(cox.fit2$estimate, time, talpha, event, x)
bh <- basehaz(coxph.fit2, centered = F)
cat("check baseline hazard vs basehaz(): max abs difference =", max(abs(H0.cox(bh$time) - bh$hazard)), "\n")

saveRDS(list(cox.fit1 = cox.fit1, cox.fit2 = cox.fit2, hr = hr, H0.cox = H0.cox), "output/03_cox.rds")
