##########################################################################
## 01_km.R
## Kaplan-Meier estimates of time to first exit from employment (overall and by
## BMI), the log-rank test, and the log(-log S) check of proportional hazards.
## Run from the repository root:  source("R/01_km.R")
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

cat("N =", length(time), " exits =", sum(event), " obese =", sum(bmi), "\n")


##########################################################################
# 1) Kaplan-Meier with left truncation
##########################################################################
kmovrll <- KMcalc(talpha, time, event)
normbmi <- KMcalc(talpha[bmi==0], time[bmi==0], event[bmi==0])
highbmi <- KMcalc(talpha[bmi==1], time[bmi==1], event[bmi==1])

# check: survfit() with counting-process data gives the same curve and risk sets
sf <- survfit(Surv(talpha, time, event) ~ bmi)
chk <- summary(sf, times = normbmi$time[normbmi$n.event > 0])$surv[1:sum(normbmi$n.event > 0)]
cat("check KM (BMI < 30) vs survfit, max abs difference:",
    max(abs(normbmi$surv[normbmi$n.event > 0] - chk)), "\n")


##########################################################################
# 2) log-rank test, with left-truncated risk sets
##########################################################################
lr <- logrank(talpha, time, event, bmi)
print(round(lr, 4))

# checks: (a) without the tie correction it is the Cox score test with Breslow ties
cat("check log-rank (Breslow variance) vs coxph score test:",
    round(lr["chisq.breslow"], 6), round(coxph(Surv(talpha, time, event) ~ bmi, ties = "breslow")$score, 6), "\n")
#         (b) ignoring left truncation it matches survdiff(), which cannot handle it
lr.notrunc <- logrank(rep(0, length(time)), time, event, bmi)
cat("check log-rank ignoring truncation vs survdiff:",
    round(lr.notrunc["chisq"], 4), round(survdiff(Surv(time, event) ~ bmi)$chisq, 4), "\n")
# The 2016 poster typed in chi2 = 8.97, p = 0.0027 by hand; the 2016 code has no log-rank call.


##########################################################################
# 3) proportional hazards check: log(-log S(t)) against log(t), by BMI
##########################################################################
# Under PH the two curves are parallel (a constant vertical gap of beta).
# Only times where 0 < S < 1, so the transform is finite.
phcheck <- list(
  normbmi = with(normbmi[normbmi$surv < 1 & normbmi$surv > 0, ], data.frame(time, lls = log(-log(surv)))),
  highbmi = with(highbmi[highbmi$surv < 1 & highbmi$surv > 0, ], data.frame(time, lls = log(-log(surv)))))

saveRDS(list(kmovrll = kmovrll, normbmi = normbmi, highbmi = highbmi, logrank = lr,
             logrank.notrunc = lr.notrunc, phcheck = phcheck), "output/01_km.rds")
