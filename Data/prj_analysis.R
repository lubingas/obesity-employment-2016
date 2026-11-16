rm(list=ls())
library(maxLik) ## a good optimizer
library(survival)
library(foreign)
library(tikzDevice)
library(xtable)
library(splines)
library(RColorBrewer)
col <- brewer.pal(3,"Dark2")

setwd("C:/Users/slubi/Dropbox/CSSS544/Project/Data/")

# prepare the data
d <- read.csv("mepsempcsv.csv")
talpha <- d$talpha
time <- d$time
event <- d$event
bmi <- d$bmi
sex <- d$sex
age <- d$age

white <- d$white
black <- d$black
other <- d$other

married <- d$married
nevermarried <- d$nevermarried
wasmarried <- d$wasmarried

uninsured <- d$uninsured
publicinsur <- d$publicinsur
privateinsur <- d$privateinsur

poor <- d$poor
lowinc <- d$lowinc
midinc <- d$midinc
highinc <- d$highinc

surv.obj <- Surv(talpha, time, event, type = "counting")

fit <- survfit(surv.obj ~ 1)
plot(fit)
summary(fit)

# write exponential likelihood function
expLik <- function(par, time, talpha, event) {
  if (par < 0)
    return(NA)
  ll <- event * (log(par) - par*time + par*talpha) + (1-event) * (par*talpha - par*time)
  return(sum(ll))
}

exp.fit <- maxLik(logLik = expLik, start = 0.1, time = time, talpha = talpha, event = event, method = "NR")
summary(exp.fit)

times <- seq(0,max(time), by = 1)

# write weibull likelihood, accounting for left truncation
weibLik <- function(par, time, talpha, event) {
  c <- par[1]
  b <- par[2]
  if (c < 0 || b < 0)
    return(NA)
  ll <- event * (log(c) + (c-1)*log(time) - c*log(b) - (time/b)^c + (talpha/b)^c) + (1-event) * ((talpha/b)^c - (time/b)^c)
  return(sum(ll))
}

weib.fit <- maxLik(logLik = weibLik, start = c(c = 1, b = 1),
                   time = time, talpha = talpha, event = event, method = "NR")
summary(weib.fit)

# write lognormal likelihood function
lnormLik <- function(par, time, talpha, event) {
  a <- par[1]
  b <- par[2]
  if (a < 0 || b < 0)
    return(NA)
  ll <- event * (dlnorm(time, a, b, log = T) - plnorm(talpha, a, b, lower.tail = F, log.p = T)) + (1 - event) * (plnorm(time, a, b, lower.tail = F, log.p = T) - plnorm(talpha, a, b, lower.tail = F, log.p = T))
  return(sum(ll))
}

lnorm.fit <- maxLik(logLik = lnormLik, start = c(a = 1, b = 1),
                    time = time, talpha = talpha, event = event, method = "NR")
summary(lnorm.fit)

# write gompertz likelihood function
gammaLik <- function(par, time, talpha, event) {
  a <- par[1]
  b <- par[2]
  if (a < 0 || b < 0)
    return(NA)
  ll <- event * (dgamma(time, a, b, log = T) - pgamma(talpha, a, b, lower.tail = F, log.p = T)) + (1 - event) * (pgamma(time, a, b, lower.tail = F, log.p = T) - pgamma(talpha, a, b, lower.tail = F, log.p = T))
  return(sum(ll))
}

gamma.fit <- maxLik(logLik = gammaLik, start = c(0.000001, 0.01), time = time, talpha = talpha, event = event, method = "BFGS")
summary(gamma.fit)

plot(x = times, y = exp(-exp.fit$estimate * times),type = "l")
points(x = times, y = pgamma(times, gamma.fit$estimate[1], gamma.fit$estimate[2], lower.tail = F, log.p = F), lty = 4, type = "l", col = "green")
points(x = times, y = plnorm(times, lnorm.fit$estimate[1], lnorm.fit$estimate[2], lower.tail = F, log.p = F), lty = 3, type = "l", col = "red")
points(x = times, y = exp(-(times/weib.fit$estimate[2])^weib.fit$estimate[1]), lty = 2, type = "l", col = "blue")

##### NOW TO ACTUALLY START MODELLING

# 1) first I look at empirical KM distribution
surv.obj <- Surv(talpha, time, event)

fit <- survfit(surv.obj ~ bmi)
plot(fit)
kmtab <- summary(fit)

## does the PH assumption hold for BMI?
with(fit, plot(log(-log(surv)) ~ log(time)))
plot(fit, fun = function(S) log(-log(S)), log = "x")


KMcalc <- function(talpha, time, failure, data = NULL) {
  if (!is.null(data)) {
    arguments <- as.list(match.call())
    time <- eval(arguments$time, data)
    failure <- eval(arguments$failure, data)
    talpha <- eval(arguments$talpha, data)
  }
  
  sorted <- data.frame(talpha, time, failure)
  
  # sort data
  sorted <- sorted[order(sorted$time),]
  
  # prepare variables
  time <- sorted$time
  talpha <- sorted$talpha
  failure <- sorted$failure
  censored <- as.numeric(sorted$failure==0)
  
  # initialize some variables we need
  riskset <- adj_riskset <- vector("integer", length=length(time))
  
  cumu.fail.no <- ave(failure, time, FUN = sum)
  cumu.cens.no <- ave(censored, time, FUN = sum)
  
  fal <- 0
  for (i in 1:length(time)) {
    riskset[i] <- adj_riskset[i] <- sum(talpha < time[i]) - fal
    fal <- fal + failure[i] + censored[i]
  }
  
  savetime <- -999.9
  saverisk <- 0.0
  
  for (i in 1:length(time)) {
    if (time[i]  == savetime) 
      adj_riskset[i] <- saverisk
    savetime <- time[i]
    saverisk <- adj_riskset[i]
  }
  
  time2 <- time[which(!duplicated(time))]
  adj_riskset <- adj_riskset[which(!duplicated(time))]
  cumu.fail.no <- cumu.fail.no[which(!duplicated(time))]
  cumu.cens.no <- cumu.cens.no[which(!duplicated(time))]
  
  # calculate conditional survival probability
  p.fail <- cumu.fail.no/adj_riskset
  p.survival <- 1-p.fail
  
  # calculate survival function
  survival <- vector("integer", length=length(p.fail))
  for (j in 2:length(p.fail)) {
    survival[1] <- p.survival[1]
    survival[j] <- survival[j-1] * p.survival[j]
  }
  return(list(time = time2, n.risk = adj_riskset, n.event = cumu.fail.no, n.censor = cumu.cens.no, survival = survival, events = sum(cumu.fail.no)))
}

# overall
kmovrll <- KMcalc(talpha = talpha, time = time, failure = event)

## by BMI
normbmi <- KMcalc(talpha = talpha[bmi==0], time = time[bmi==0], failure = event[bmi==0])
highbmi <- KMcalc(talpha = talpha[bmi==1], time = time[bmi==1], failure = event[bmi==1])

## check PH assumption for BMI
y1 <- log(-log(normbmi$survival))
y2 <- log(-log(highbmi$survival))
xx1 <- normbmi$time
xx2 <- highbmi$time

tikz(file = "C:/Users/slubi/Dropbox/CSSS544/Project/phcheck.tex", width = 5, height = 4)
par(mar = c(4,4,3,1)+0.1)
plot(x = xx1, y = y1, lwd = 2, col = col[2], ann = F,
     type = "l", xaxt="n", yaxt="n", log = "x", frame.plot = F)
points(x = xx2, y = y2, lwd = 2, type = "l", col=col[3])
axis(1, at=c(0.05, 0.5, 5, 50), cex.axis = 0.6)
axis(2, at=c(-6, -4, -2, 0), cex.axis = 0.6)
mtext(side = 1, text = paste("$t$"), line = 3, cex = 0.7)
mtext(side = 2, text = paste("$\\log(-\\log(S(t)))$"), line = 3, cex = 0.7)
dev.off()

# Now overlay the distributions
# 1) a log normal distribution
lnorm.lowbmi.fit <- maxLik(logLik = lnormLik, start = c(5, 2), time = time[bmi==0], talpha = talpha[bmi==0], event = event[bmi==0], method = "NR")
summary(lnorm.lowbmi.fit)

lnorm.highbmi.fit <- maxLik(logLik = lnormLik, start = c(5, 1), time = time[bmi==1], talpha = talpha[bmi==1], event = event[bmi==1], method = "NR")
summary(lnorm.highbmi.fit)

# 2) a weibull distribution
weib.lowbmi.fit <- maxLik(logLik = weibLik, start = c(0.5, 3000), time = time[bmi==0], talpha = talpha[bmi==0], event = event[bmi==0], method = "NR")
summary(weib.lowbmi.fit)

weib.highbmi.fit <- maxLik(logLik = weibLik, start = c(0.5, 3000), time = time[bmi==1], talpha = talpha[bmi==1], event = event[bmi==1], method = "NR")
summary(weib.highbmi.fit)

# Plot Kaplan Meier curves
#  options(tikzDefaultEngine = "xetex") #, tikzLatexPackages = c(getOption("tikzLatexPackages")))

y3 <- normbmi$survival
y4 <- normbmi$survival[normbmi$n.censor > 0]
y5 <- highbmi$survival
y6 <- highbmi$survival[highbmi$n.censor > 0]
y7 <- plnorm(times, lnorm.lowbmi.fit$estimate[1], lnorm.lowbmi.fit$estimate[2], lower.tail = F, log.p = F)
y8 <- plnorm(times, lnorm.highbmi.fit$estimate[1], lnorm.highbmi.fit$estimate[2], lower.tail = F, log.p = F)
y9 <- exp(-(times/weib.lowbmi.fit$estimate[2])^weib.lowbmi.fit$estimate[1])
y10 <- exp(-(times/weib.highbmi.fit$estimate[2])^weib.highbmi.fit$estimate[1])

xx3 <- normbmi$time
xx4 <- normbmi$time[normbmi$n.censor > 0]
xx5 <- highbmi$time
xx6 <- highbmi$time[highbmi$n.censor > 0]

tikz(file = "C:/Users/slubi/Dropbox/CSSS544/Project/kmlog.tex", width = 5, height = 4)
par(mar = c(4,4,3,1)+0.1)
plot(x = xx3, y = y3, ylim = c(0,1), yaxt="n", xaxt = "n", type = "s", col=col[1], 
     frame.plot = F, ann = F, cex.lab = 0.7)
points(x = xx4, y = y4, "p", pch = 3, col=col[1])
points(x = xx5, y = y5, "s", col = col[2])
points(x = xx6, y = y6, "p", pch = 3, col=col[2])
axis(2, at=c(0, 0.25, 0.5, 0.75, 1), cex.axis = 0.6)
axis(1, at=c(0, 10, 20, 30, 40), cex.axis = 0.6)
points(x = times, y = y7, lwd = 2, lty = 2, type = "l", col = col[1])
points(x = times, y = y8, lwd = 2, lty = 2, type = "l", col = col[2])
points(x = times, y = y9, lwd = 2, lty = 3, type = "l", col = col[1])
points(x = times, y = y10, lwd = 2, lty = 3, type = "l", col = col[2])
mtext(paste("Logrank test: $\\chi^2 = ", 8.97,", df = ", 1, ", p = ", 0.002738, "$"), cex = 0.6)
mtext(side = 1, text = paste("$t$"), line = 3, cex = 0.7)
mtext(side = 2, text = paste("$S(t)$"), line = 3, cex = 0.7)
legend(x=25, y = 0.95, , title = "Group", legend = c(paste("$<30$"), paste("$\\geq 30$")), bty = "n", lty = c(1, 1), col = c(col[1], col[2]), cex = 0.4, ncol = 2) 
legend(x=25, y = 0.75, title = "Dist.", legend = c("LN", "W"), bty = "n", lty = c(2, 3), cex = 0.4, ncol = 2)
dev.off()

## Now the cox proportional hazards model likelihood
coxPHLik <- function(par, time, talpha, failure, x) {
  
  # sorting by time
  talpha <- talpha[order(time)]
  failure <- failure[order(time)]
  censored <- as.numeric(failure==0)
  x <- x[order(time),]
  time <- sort(time)
  
  # compute partial likelihoods
  beta <- par[1:ncol(x)]
  xb <- x%*%beta
  pl <- exp(xb)
  
  ## now accumulate risk set and adjust for ties
  for (i in 1:length(time)) {
    # which people have failed by this time?
    riskset[i] <- adj_riskset[i] <- sum((talpha < time[i]) * pl) - sum((failure == 1 & time < time[i]) * pl) - sum((censored == 1 & time < time[i]) * pl)
  }
  
  ## breslow method to adjust for ties
  savetime <- -999.9
  saverisk <- 0.0
  
  for (i in 1:length(time)) {
    if (time[i]  == savetime) 
      adj_riskset[i] <- saverisk
    savetime <- time[i]
    saverisk <- adj_riskset[i]
  }
  
  ll <- log(ifelse(failure==1, pl/adj_riskset, 1))
  sum(ll)
}

## 1) prepare the data
x <- model.matrix(~ agey1x + sex + black + other + wasmarried + nevermarried + publicinsur + privateinsur + lowinc + midinc + highinc + bmi - 1, data = s)

## 2) run model
cox_mod2 <- maxLik(logLik = coxPHLik, start = c(Bage = 0, Bsex = 0, Bblack = 0, Bother = 0, Bwasmarried = 0,
                                                Bnevermarried = 0, Bpublicinsur = 0, Bprivateinsur = 0,
                                                Blowinc = 0, Bmidinc = 0, Bhighinc = 0, Bbmi = 0),
                   talpha = talpha, time = time, failure = event, x = x, method = "NR", iterlim = 10000)
summary(cox_mod2)

coxph.fit1 <- coxph(Surv(talpha, time, event) ~ bmi, method = "breslow")
summary(coxph.fit1)
kmfit1 <- survfit(coxph(Surv(talpha, time, event) ~ strata(bmi), method = "breslow"))
plot(kmfit1, fun = "cloglog")

coxph.fit2 <- coxph(Surv(talpha, time, event) ~ age + sex + black + other + wasmarried + nevermarried + publicinsur + privateinsur + lowinc + midinc + highinc + bmi, method = "breslow")
summary(coxph.fit2)

coxph.fit3 <- coxph(Surv(talpha, time, event) ~ age + sex + black + other + wasmarried + nevermarried + publicinsur + privateinsur + lowinc + midinc + highinc + strata(bmi), method = "breslow")

haz <- basehaz(coxph.fit3)
timhz <- haz$time; bhaz <- haz$hazard
tikz(file = "C:/Users/slubi/Dropbox/CSSS544/Project/hzlog.tex", width = 3, height = 3)
par(mar = c(4,4,3,1)+0.1)
plot(timhz, bhaz, ann = F, yaxt = "n", xaxt = "n", frame.plot = F, col = "darkgreen", cex= 0.5)
axis(2, at=c(0, 0.5, 1.0, 1.5, 2.0, 2.5), cex.axis = 0.5)
axis(1, at=c(0, 10, 20, 30, 40), cex.axis = 0.5)
mtext(side = 1, text = paste("$t$"), line = 3, cex = 0.5)
mtext(side = 2, text = paste("$H(t)$"), line = 3, cex = 0.5)
dev.off()

kmfit2 <- survfit(coxph.fit3)
plot(kmfit2, fun="cumhaz")


# let us run some diagnostics
# 1) Cox Snell residuals
coxsnellres <- event - resid(coxph.fit2, type = "martingale")
coxsnellfit <- survfit(Surv(coxsnellres, event) ~ 1)
Htilde <- cumsum(coxsnellfit$n.event/coxsnellfit$n.risk)
plot(coxsnellfit$time, log(Htilde), type='s', col=col[1])
abline(0, 1,lty = 2)

coxsnellfit0 <- survfit(Surv(coxsnellres[bmi==0], event[bmi==0]) ~ 1)
Htilde0 <- cumsum(coxsnellfit0$n.event/coxsnellfit0$n.risk)
plot(coxsnellfit0$time, log(Htilde0), type='s', col=col[2], lty = 2)

coxsnellfit1 <- survfit(Surv(coxsnellres[bmi==1], event[bmi==1]) ~ 1)
Htilde1 <- cumsum(coxsnellfit1$n.event/coxsnellfit1$n.risk)
points(coxsnellfit1$time, log(Htilde1), type='s', col=col[3], lty = 2)

# 2) Martingale residuals
marting <- resid(coxph.fit2,type='martingale')
plot(age, marting)
tikz(file = "C:/Users/slubi/Dropbox/CSSS544/Project/mres.tex", width = 3, height = 3)
par(mar = c(4,4,3,1)+0.1)
plot(bmi[b], marting[b], frame.plot = F, cex = 0.5, col = col[2], ann = F, yaxt = "n", xaxt = "n")
axis(1, at=c(0, 1), labels = c(paste("$< 30$"), paste("$\\geq 30$")), cex.axis = 0.5)
axis(2, at=c(-1.0, -0.5, 0, 0.5, 1.0), cex.axis = 0.5)
mtext(side = 2, text = "Martingale residuals", line = 2, cex = 0.5)
dev.off()


# deviance residuals
preds <- coxph.fit2$linear.predictors
dresids <- resid(coxph.fit2, type = "deviance")
# b <- sample(1:length(preds), length(seq(1, length(preds), 10)))
b <- seq(1, length(preds), 10)

tikz(file = "C:/Users/slubi/Dropbox/CSSS544/Project/devr.tex", width = 3, height = 3)
par(mar = c(4,4,3,1)+0.1)
plot(x = preds[b], y = dresids[b], col = col[3], cex = 0.5, cex.axis = 0.5, ann = F, frame.plot = F)
mtext(side = 1, text = "Linear Predictor", line = 2, cex = 0.5)
mtext(side = 2, text = "Residuals", line = 2, cex = 0.5)
abline(h=0, col = col[1], lty=2)
dev.off()

methods(plot)
getAnywhere(plot.cox.zph)

plot.cox.zph2 <- function (x, resid = TRUE, se = TRUE, df = 4, nsmo = 40, var, 
                           ...) 
{
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  if (qmat$rank < df) 
    stop("Spline fit is singular, try a smaller degrees of freedom")
  if (se) {
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
  }
  ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
  if (missing(var)) 
    var <- 1:nvar
  else {
    if (is.character(var)) 
      var <- match(var, dimnames(yy)[[2]])
    if (any(is.na(var)) || max(var) > nvar || min(var) < 
        1) 
      stop("Invalid variable requested")
  }
  if (x$transform == "log") {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  }
  else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])
    indx <- !duplicated(xx)
    apr1 <- approx(xx[indx], xtime[indx], seq(min(xx), max(xx), 
                                              length = 17)[2 * (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- approx(xtime[indx], xx[indx], temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  for (i in var) {
    y <- yy[, i]
    yhat <- pmat %*% qr.coef(qmat, y)
    if (resid) 
      yr <- range(yhat, y)
    else yr <- range(yhat)
    if (se) {
      temp <- 2 * sqrt(x$var[i, i] * seval)
      yup <- yhat + temp
      ylow <- yhat - temp
      yr <- range(yr, yup, ylow)
    }
    if (x$transform == "identity") 
      plot(range(xx), yr, type = "n", xlab = "Time", ylab = ylab[i], 
           ...)
    else if (x$transform == "log") 
      plot(range(xx), yr, type = "n", xlab = "Time", ylab = ylab[i], 
           log = "x", ...)
    else {
      plot(range(xx), yr, type = "n", xlab = "Time", ylab = ylab[i], 
           axes = FALSE, ...)
      axis(1, xaxisval, xaxislab)
      axis(2)
      box()
    }
    if (resid) 
      points(xx, y, ...)
    lines(pred.x, yhat)
    if (se) {
      lines(pred.x, yup, lty = 2)
      lines(pred.x, ylow, lty = 2)
    }
  }
}

# Shoenfeld
shresid <- resid(coxph.fit2, type = "schoenfeld")
failt <- time[event==1]
scatter.smooth(sort(failt), shresid[,12], log = "x")

temp <- cox.zph(coxph.fit2, transform = "log")
print(temp) # display the results

tikz(file = "C:/Users/slubi/Dropbox/CSSS544/Project/sho.tex", width = 3, height = 3)
par(mar = c(4,4,3,1)+0.1)
plot.cox.zph2(temp[12], frame.plot = F, ann = F, cex = 0.6, xaxt = "n", yaxt = "n")
axis(1, at=c(0, 0.05, 0.5, 5, 50), cex.axis = 0.5)
axis(2, at=c(-2, 0, 2, 4), cex.axis = 0.5)
abline(h=0, col=col[3], lty=2)
mtext(side = 1, text = paste("$t$"), line = 2, cex = 0.5)
mtext(side = 2, text = paste("$\\beta(t)$"), line = 2, cex = 0.5)
dev.off()

# estimate as a log-normal PH model
# the likelihood
lnormPHLik <- function(par, time, talpha, event, x) {
  a <- par[1]
  b <- par[2]
  if (a < 0 || b < 0)
    return(NA)
  beta <- par[3:(ncol(x)+2)]
  xb <- x%*%beta
  exfail <- ((dlnorm(time, a, b, log = F) * plnorm(time, a, b, lower.tail = F, log.p = F)^(exp(xb)- 1)  * exp(xb))^event) * ((plnorm(time, a, b, lower.tail = F, log.p = F)^exp(xb))^(1 - event))
  trunc <- ifelse(talpha > 0, (plnorm(talpha, a, b, lower.tail = F, log.p = F)^exp(xb)), 1)
  
  ll <- log(exfail/trunc)
  return(sum(ll))
} 

# BMI Only
x1 <- cbind(bmi)
lnormPH.fit1 <- maxLik(logLik = lnormPHLik, start = c(a = 1.0, b = 10.0, Bbmi = 0),
                      time = time, talpha = talpha, event = event, x = x1, 
                      method = "NR",iterlim = 10000)
summary(lnormPH.fit1)


# all covariates
x <- cbind(age, sex, black, other, wasmarried, nevermarried, publicinsur, privateinsur, lowinc, midinc, highinc, bmi)
lnormPH.fit2 <- maxLik(logLik = lnormPHLik, start = c(a = 1, b = 2, Bage = 0, Bsex = 0, Bblack = 0, Bother = 0,
                                                      Bwasmarried = 0, Bnevermarried = 0,
                                                      Bpublicinsur = 0, Bprivateinsur = 0,
                                                      Blowinc = 0, Bmidinc = 0, Bhighinc = 0, Bbmi = 0),
                       time = time, talpha = talpha, event = event, x = x, 
                       method = "NR",iterlim = 10000)
summary(lnormPH.fit2)
lnormPH.fit2.fisher_info <- solve(-lnormPH.fit2$hessian)
lnormPH.prop_sigma <- sqrt(diag(lnormPH.fit2.fisher_info))

se(coxph.fit2)

coxph.fit2.fisher_info <- solve(-coxph.fit2$var)
coxph.prop_sigma <- sqrt(diag(coxph.fit2.fisher_info))

lnorm2 <- c(0.000100000096, 0.572971048787, 1.289770412919,-0.01772858623,-0.07375788051,0.062376570548,
            -0.34527374794, 0.206011767883, 0.170343411284, 0.187202351611, -0.17886814871, -0.19415404664, -0.31718916119, -0.53303065189, 
            0.183981796822)

lnorm2se <- c(0.170576003708, 0.063757460710, 5.193965127669, 0.003461956775, 0.063813201392, 0.077194512005, 0.123477367187,
              0.096208890281, 0.079674936241, 0.100912760536, 0.078491845879, 0.092479327210, 0.088470920768, 0.104794552317,
              0.069246889104)

print(xtable(cbind(cox = c(sigma = NA, a = NA, b = NA, paste0(round(coxph.fit2$coefficients, 4), " (", round(unlist(summary(coxph.fit2)[7])[13:24], 4),")"), coxph.fit2$loglik[2]),
                   lnorm = c(NA,paste0(round(lnormPH.fit2$estimate, 4), " (",round(lnormPH.prop_sigma, 4),")"), lnormPH.fit2$maximum),
                   lnorm2 = c(paste0(round(lnorm2, 4), " (", round(lnorm2se,4), ")"), -3381.601)), 
             digits = c(0,5,5,5), caption = "Parameter estimates"), caption.placement = "top")

print(xtable(cbind(coxest = c(sigma = NA, a = NA, b = NA, coxph.fit2$coefficients, N = length(time), L = coxph.fit2$loglik[2]),
                   coxse = c(NA, NA, NA, unlist(summary(coxph.fit2)[7])[25:36], NA, NA),
                   lnormest = c(NA, lnormPH.fit2$estimate, length(time), lnormPH.fit2$maximum),
                   lnormse = c(NA, lnormPH.prop_sigma,NA, NA),
                   lnorm2est = c(lnorm2,length(time), -3381.601),
                   lnorm2se = c(lnorm2se, NA, NA)),
             digits = c(0,5,5,5,5,5,5), caption = "Parameter estimates"), caption.placement = "top")

paste0(lnormPH.fit2$estimate, "(",lnormPH.prop_sigma,")")

paste0(coxph.fit2$coefficients, "(", unlist(summary(coxph.fit2)[7])[13:24],")")
paste0(lnorm2, "(", lnorm2se, ")")

# how about a model with unobserved heterogeneity?
# the likelihood
lnormPHLik2 <- function(par, time, talpha, event, x) {
  a <- par[1]
  b <- par[2]
  sigma <- par[3]
  if (a < 0 || b < 0 || sigma < 0)
    return(NA)
  beta <- par[4:(ncol(x)+3)]
  newbeta <- c(beta, 1)
  i <- 1
  inll <- as.numeric()
  for (z in seq(-5, 5, 0.1)) {
    newx <- cbind(x, z)
    xb <- newx%*%newbeta
    p <- dnorm(z, 0, sigma)
    exfail <- ((dlnorm(time, a, b, log = F) * plnorm(time, a, b, lower.tail = F, log.p = F)^(exp(xb)- 1)  * exp(xb))^event) * ((plnorm(time, a, b, lower.tail = F, log.p = F)^exp(xb))^(1 - event))
    trunc <- ifelse(talpha > 0, (plnorm(talpha, a, b, lower.tail = F, log.p = F)^exp(xb)), 1)
    inll[i] <- p * sum(log(exfail/trunc))
    i <- i + 1
    }
  totll <- sum(inll)
  return(totll)
} 

lnormPH_mod2 <- maxLik(logLik = lnormPHLik2, start = c(a = 10, b = 1, sigma = 0.5, Bbmi = 0),
                       time = time, talpha = talpha, event = event, x = x1, 
                       method = "BFGS",iterlim = 10000)
summary(lnormPH_mod2)

lnormPH.fit3 <- maxLik(logLik = lnormPHLik2, start = c(a = 1, b = 2, sigma = 0.5, Bage = 0, Bsex = 0, Bblack = 0, Bother = 0,
                                                       Bwasmarried = 0, Bnevermarried = 0,
                                                       Bpublicinsur = 0, Bprivateinsur = 0,
                                                       Blowinc = 0, Bmidinc = 0, Bhighinc = 0, Bbmi = 0),
                       time = time, talpha = talpha, event = event, x = x, 
                       method = "BFGS",iterlim = 10000)

summary(lnormPH.fit3)


# how about a model with unobserved heterogeneity?

llk <- function(par, time, talpha, event, x) {
  ## integrate over z
  intfx <- function(z) {
    p <- dnorm(z, 0, sigma)
    exfail <- ((dlnorm(time, a, b, log = F) * plnorm(time, a, b, lower.tail = F, log.p = F)^(exp(xb+z)- 1)  * exp(xb+z))^event) * ((plnorm(time, a, b, lower.tail = F, log.p = F)^exp(xb+z))^(1 - event))
    trunc <- ifelse(talpha > 0, (plnorm(talpha, a, b, lower.tail = F, log.p = F)^exp(xb+z)), 1)
    sum(log(p * exfail/trunc))
  }
  return(llk)
}

intllk <- function(par, time, talpha, event, x) {
  a <- par[1]
  b <- par[2]
  sigma <- par[3]
  if (a < 0 || b < 0 || sigma < 0)
    return(NA)
  beta <- par[4:(ncol(x)+3)]
  xb <- x%*%beta
  integrate(llk(par, time, talpha, event, x), lower = -Inf, upper = Inf)$value
}


lnormPH.fit4 <- maxLik(logLik = intllk, start = c(a = 1, b = 1, sigma = 0.5, Bbmi = 0),
                       time = time, talpha = talpha, event = event, x = x1, 
                       method = "NR", iterlim = 10000)

summary(lnormPH.fit4)

exp <- data.frame(talpha,time,event,age, sex, black, other, wasmarried, nevermarried, publicinsur, privateinsur, lowinc, midinc, highinc, bmi)
write.table(exp, "anafil.dat",row.names = F)


beta <- 0.1
xb <- x1%*%beta
xb + 2

