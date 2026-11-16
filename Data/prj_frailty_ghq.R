rm(list=ls())
library(maxLik)   ## a good optimizer
library(survival)
library(statmod)  ## gauss.quad() gives the Gauss-Hermite nodes and weights

setwd("~/Projects/obesity-employment-2016/Data/")

##########################################################################
## Log-normal PH model with continuous (normal) unobserved heterogeneity,
## integrated with Gauss-Hermite quadrature instead of MLE's integrator.
## Starts from prj_analysis.R (lnormPHLik / lnormPHLik2) - see the
## "CHANGED" comments for what is different and why.
##
## Sections: 1) no heterogeneity  2) heterogeneity via Gauss-Hermite
##           3) fits + MLE comparison + LR tests  4) profile likelihood for sigma
##           5) other baselines (exponential, Weibull, Gompertz, log-logistic,
##              gamma, piecewise exponential) x (none, normal, gamma heterogeneity)
## Run time: about 1 minute for sections 1-3, 6-8 minutes each for 4 and 5.
##########################################################################

# prepare the data (same as prj_analysis.R)
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

x1 <- cbind(bmi)
x <- cbind(age, sex, black, other, wasmarried, nevermarried, publicinsur, privateinsur, lowinc, midinc, highinc, bmi)

# what MLE gave (hard-coded in prj_analysis.R as lnorm2): sigma, a, b, then the 12 betas
lnorm2 <- c(0.000100000096, 0.572971048787, 1.289770412919,-0.01772858623,-0.07375788051,0.062376570548,
            -0.34527374794, 0.206011767883, 0.170343411284, 0.187202351611, -0.17886814871, -0.19415404664, -0.31718916119, -0.53303065189,
            0.183981796822)
lnorm2.ll <- -3381.601


##########################################################################
# 1) log-normal PH WITHOUT heterogeneity (lnormPHLik from prj_analysis.R, rewritten on the log scale)
##########################################################################
# The 2016 version builds the likelihood on the probability scale and logs at the end:
#   [f(t) S(t)^(e^xb - 1) e^xb]^event [S(t)^e^xb]^(1-event) / S(talpha)^e^xb
# Since f(t)/S(t) = h0(t), this is the same as h0(t)^event e^(xb*event) S(t)^e^xb / S(talpha)^e^xb,
# so the log is
#   event*(log h0(t) + xb) + e^xb * (log S0(t) - log S0(talpha))
# CHANGED: work with log h0 and log S0 directly (plnorm(..., log.p = T)). Same answer, but no
#          underflow when S(t) is tiny, and S0(0) = 1 handles talpha = 0 without the ifelse().
# CHANGED: estimate log(b) instead of b, so b > 0 automatically; no "if (a < 0) return(NA)"
#          (a is the meanlog and is allowed to be negative).
lnormPHLik <- function(par, time, talpha, event, x) {
  a <- par[1]
  b <- exp(par[2])
  beta <- par[3:(ncol(x)+2)]
  xb <- drop(x%*%beta)
  logS0t <- plnorm(time, a, b, lower.tail = F, log.p = T)    # log S0(t)
  logS0a <- plnorm(talpha, a, b, lower.tail = F, log.p = T)  # log S0(talpha), = 0 when talpha = 0
  logh0t <- dlnorm(time, a, b, log = T) - logS0t             # log h0(t) = log f0(t) - log S0(t)
  ll <- event * (logh0t + xb) + exp(xb) * (logS0t - logS0a)
  return(sum(ll))
}


##########################################################################
# 2) the same model WITH heterogeneity: Gauss-Hermite quadrature
##########################################################################
# For each person the likelihood is an integral over z ~ N(0, sigma):
#   L_i = integral g(z|0,sigma) L_i(z) dz,   L_i(z) = the no-heterogeneity likelihood with xb -> xb + z
# and the log-likelihood is sum_i log(L_i). The integral sits INSIDE each person's likelihood.
#
# Why the earlier attempts didn't work:
#  - lnormPHLik2 (grid z = -5, -4.9, ..., 5) computed sum_z p(z) * sum_i log L_i(z): that weights the
#    log-likelihood of the whole sample by z, instead of sum_i log( sum_z p(z) L_i(z) dz ). It also left
#    out the grid step dz, and the grid was fixed at +/-5 whatever sigma was.
#  - the integrate() version gave integrate() a function that returns one value per person (6,075
#    numbers for each z); integrate() needs one number per z, so it would take 6,075 separate
#    integrate() calls per likelihood evaluation.
#
# Gauss-Hermite: integral exp(-u^2) f(u) du ~= sum_k w_k f(u_k). Substituting z = sqrt(2)*sigma*u
# turns the normal density into exp(-u^2)/sqrt(pi), so
#   L_i ~= sum_k [w_k / sqrt(pi)] * L_i(z_k),   z_k = sqrt(2)*sigma*u_k
# That is just a weighted sum over K fixed nodes, so we can do all people and all nodes at once
# as an N x K matrix.
#
# CHANGED: the sum over nodes is done on the log scale with the "log-sum-exp" trick
#          (log sum_k e^{v_k} = m + log sum_k e^{v_k - m}, m = max_k v_k) so nothing underflows.
# CHANGED: sigma is estimated as log(sigma) so it stays > 0.
#
# Left truncation, two versions (trunc = ):
#  "mle"   = the poster/MLE likelihood: divide by S(talpha|z) INSIDE the integral
#            L_i = integral g(z) L_i(z) / S(talpha|z) dz
#            (treats z as drawn when the person enters the survey)
#  "ratio" = standard delayed entry: divide the whole integral by the probability of still being
#            employed at entry, averaged over z
#            L_i = integral g(z) L_i(z) dz / integral g(z) S(talpha|z) dz
#            (z is fixed from the start of the job, so people still employed at entry are a
#             selected, lower-frailty group). The two agree when sigma = 0.

gh <- gauss.quad(30, kind = "hermite")  # 30 nodes; accurate to ~1e-6 per person for sigma <= 1

# log(sum_k exp(M[, k])) for every row of M
rowLogSumExp <- function(M) {
  m <- M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
  m + log(rowSums(exp(M - m)))
}

lnormPHLik.gh <- function(par, time, talpha, event, x, trunc = "mle", nodes = gh$nodes, weights = gh$weights) {
  a <- par[1]
  b <- exp(par[2])
  sigma <- exp(par[3])
  beta <- par[4:(ncol(x)+3)]
  xb <- drop(x%*%beta)
  logS0t <- plnorm(time, a, b, lower.tail = F, log.p = T)
  logS0a <- plnorm(talpha, a, b, lower.tail = F, log.p = T)
  logh0t <- dlnorm(time, a, b, log = T) - logS0t

  z <- sqrt(2) * sigma * nodes               # the K values of z
  lw <- log(weights) - 0.5*log(pi)           # log(w_k / sqrt(pi))
  eta <- outer(xb, z, "+")                   # N x K matrix of xb + z_k

  lnum <- event * (logh0t + eta) + exp(eta) * logS0t   # log of h(t|z)^event S(t|z)
  lden <- exp(eta) * logS0a                            # log of S(talpha|z)

  if (trunc == "mle") {
    ll <- rowLogSumExp(sweep(lnum - lden, 2, lw, "+"))
  } else {
    ll <- rowLogSumExp(sweep(lnum, 2, lw, "+")) - rowLogSumExp(sweep(lden, 2, lw, "+"))
  }
  return(sum(ll))
}

# quick check of the quadrature against integrate(), person by person (sigma = 0.5)
check.gh <- function(par, i, x) {
  a <- par[1]; b <- exp(par[2]); sigma <- exp(par[3]); beta <- par[4:(ncol(x)+3)]
  xb <- drop(x[i, , drop = F]%*%beta)
  Lz <- function(z) {
    lf <- event[i] * (dlnorm(time[i], a, b, log = T) - plnorm(time[i], a, b, lower.tail = F, log.p = T) + xb + z) +
      exp(xb + z) * (plnorm(time[i], a, b, lower.tail = F, log.p = T) - plnorm(talpha[i], a, b, lower.tail = F, log.p = T))
    out <- dnorm(z, 0, sigma) * exp(lf)
    out[!is.finite(out)] <- 0
    out
  }
  c(integrate = log(integrate(Lz, -12*sigma, 12*sigma, rel.tol = 1e-10)$value),
    gh = lnormPHLik.gh(par, time[i], talpha[i], event[i], x[i, , drop = F], trunc = "mle"))
}
ids <- c(which(event == 1)[1:3], which(event == 0)[1:3])
print(sapply(ids, check.gh, par = c(0.57, log(1.29), log(0.5), 0.18), x = x1))


##########################################################################
# 3) fit: BMI only and all covariates, no heterogeneity vs heterogeneity
##########################################################################
# Estimates come back as a, log(b), log(sigma); this reports b and sigma with delta-method SEs.
# If sigma has run to its lower bound (sigma -> 0, which is what MLE found), the Hessian is singular
# in log(sigma); SEs for the other parameters are then taken from the Hessian without that row/column.
showfit <- function(fit) {
  est <- coef(fit)
  H <- hessian(fit)
  logpars <- grep("^log", names(est))
  hetpar <- intersect(c("logsigma", "logtheta"), names(est))
  atbound <- length(hetpar) == 1 && exp(est[hetpar]) < 1e-2
  keep <- if (atbound) which(names(est) != hetpar) else seq_along(est)
  se <- rep(NA, length(est))
  se[keep] <- sqrt(diag(solve(-H[keep, keep])))
  out <- cbind(est = est, se = se)
  out[logpars, "est"] <- exp(est[logpars])
  out[logpars, "se"] <- exp(est[logpars]) * se[logpars]   # delta method
  rownames(out) <- sub("^log", "", names(est))
  out
}

# BMI only
lnormPH.fit1 <- maxLik(logLik = lnormPHLik, start = c(a = 1.0, logb = 0, Bbmi = 0),
                       time = time, talpha = talpha, event = event, x = x1,
                       method = "NR", iterlim = 10000)

# heterogeneity, start from the no-heterogeneity estimates with sigma = 0.5
st1 <- c(coef(lnormPH.fit1)[1:2], logsigma = log(0.5), coef(lnormPH.fit1)[-(1:2)])
lnormPH.fit1h.mle <- maxLik(logLik = lnormPHLik.gh, start = st1,
                            time = time, talpha = talpha, event = event, x = x1, trunc = "mle",
                            method = "BFGS", iterlim = 10000)
lnormPH.fit1h.ratio <- maxLik(logLik = lnormPHLik.gh, start = st1,
                              time = time, talpha = talpha, event = event, x = x1, trunc = "ratio",
                              method = "BFGS", iterlim = 10000)

# all covariates
lnormPH.fit2 <- maxLik(logLik = lnormPHLik, start = c(a = 1, logb = 0, Bage = 0, Bsex = 0, Bblack = 0, Bother = 0,
                                                      Bwasmarried = 0, Bnevermarried = 0,
                                                      Bpublicinsur = 0, Bprivateinsur = 0,
                                                      Blowinc = 0, Bmidinc = 0, Bhighinc = 0, Bbmi = 0),
                       time = time, talpha = talpha, event = event, x = x,
                       method = "NR", iterlim = 10000)

st2 <- c(coef(lnormPH.fit2)[1:2], logsigma = log(0.5), coef(lnormPH.fit2)[-(1:2)])
lnormPH.fit2h.mle <- maxLik(logLik = lnormPHLik.gh, start = st2,
                            time = time, talpha = talpha, event = event, x = x, trunc = "mle",
                            method = "BFGS", iterlim = 10000)
lnormPH.fit2h.ratio <- maxLik(logLik = lnormPHLik.gh, start = st2,
                              time = time, talpha = talpha, event = event, x = x, trunc = "ratio",
                              method = "BFGS", iterlim = 10000)

# side by side with MLE (lnorm2): sigma, a, b, betas
mle.compare <- cbind(MLE = lnorm2,
                     R.gh = showfit(lnormPH.fit2h.mle)[c(3, 1, 2, 4:15), "est"])
mle.compare <- rbind(mle.compare, logLik = c(lnorm2.ll, logLik(lnormPH.fit2h.mle)))
print(round(mle.compare, 5))

# LR test of sigma = 0. sigma = 0 is on the edge of the parameter space, so the LR statistic is a
# 50:50 mix of chi2(0) and chi2(1): halve the usual chi2(1) p-value.
lrtest.sigma <- function(fit0, fit1) {
  lr <- max(0, 2 * (logLik(fit1) - logLik(fit0)))
  c(LR = lr, p = 0.5 * pchisq(lr, df = 1, lower.tail = F))
}

het.summary <- rbind(
  "BMI only, MLE truncation"      = c(sigma = exp(coef(lnormPH.fit1h.mle)["logsigma"]),   ll0 = logLik(lnormPH.fit1), ll1 = logLik(lnormPH.fit1h.mle),   lrtest.sigma(lnormPH.fit1, lnormPH.fit1h.mle)),
  "BMI only, delayed entry"       = c(sigma = exp(coef(lnormPH.fit1h.ratio)["logsigma"]), ll0 = logLik(lnormPH.fit1), ll1 = logLik(lnormPH.fit1h.ratio), lrtest.sigma(lnormPH.fit1, lnormPH.fit1h.ratio)),
  "All covariates, MLE truncation" = c(sigma = exp(coef(lnormPH.fit2h.mle)["logsigma"]),   ll0 = logLik(lnormPH.fit2), ll1 = logLik(lnormPH.fit2h.mle),   lrtest.sigma(lnormPH.fit2, lnormPH.fit2h.mle)),
  "All covariates, delayed entry"  = c(sigma = exp(coef(lnormPH.fit2h.ratio)["logsigma"]), ll0 = logLik(lnormPH.fit2), ll1 = logLik(lnormPH.fit2h.ratio), lrtest.sigma(lnormPH.fit2, lnormPH.fit2h.ratio)))
colnames(het.summary)[1] <- "sigma"
print(round(het.summary, 4))

showfit(lnormPH.fit1h.ratio)
showfit(lnormPH.fit2h.ratio)


##########################################################################
# 4) profile likelihood for sigma (all covariates)
##########################################################################
# Fix sigma on a grid and maximise over everything else (maxLik's fixed = argument). If the profile
# is highest at sigma = 0, the boundary estimate is the global maximum, not a stuck optimiser.
sig.grid <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5)
gh60 <- gauss.quad(60, kind = "hermite")  # more nodes for the larger sigmas
prof <- matrix(NA, length(sig.grid), 2, dimnames = list(sig.grid, c("mle", "ratio")))
for (tr in c("mle", "ratio")) {
  st <- st2
  for (j in seq_along(sig.grid)) {
    st["logsigma"] <- log(sig.grid[j])
    f <- maxLik(logLik = lnormPHLik.gh, start = st, fixed = "logsigma",
                time = time, talpha = talpha, event = event, x = x, trunc = tr,
                nodes = gh60$nodes, weights = gh60$weights, method = "BFGS", iterlim = 10000)
    prof[j, tr] <- logLik(f)
    st <- coef(f)   # warm start the next sigma from this one
  }
}
prof <- rbind("0" = c(logLik(lnormPH.fit2), logLik(lnormPH.fit2)), prof)
print(round(prof, 3))

pdf("frailty_profile.pdf", width = 5, height = 4)
par(mar = c(4,4,1,1)+0.1)
plot(as.numeric(rownames(prof)), prof[, "ratio"], type = "b", pch = 19, col = "#0072B2",
     ylim = range(prof), xlab = expression(sigma), ylab = "Profile log-likelihood", frame.plot = F)
points(as.numeric(rownames(prof)), prof[, "mle"], type = "b", pch = 17, col = "#E69F00")
legend("bottomleft", legend = c("Delayed entry", "MLE truncation"), col = c("#0072B2", "#E69F00"),
       pch = c(19, 17), bty = "n")
dev.off()


##########################################################################
# 5) other baseline hazards, with and without heterogeneity
##########################################################################
# Same PH structure h(t|x,z) = h0(t) exp(xb + z) for every baseline, so the betas are always log
# hazard ratios and comparable across models. Each baseline only needs log h0(t) and the cumulative
# hazard H0(t) = -log S0(t). Parameters are on unconstrained scales (logs for positive ones).
#   exponential  h0 = lambda
#   weibull      S0 = exp(-(t/b)^c)                       (the weibLik parameterisation in prj_analysis.R)
#   gompertz     h0 = r exp(g t), H0 = (r/g)(exp(g t) - 1) (g can be negative)
#   loglogistic  S0 = 1/(1 + (t/b)^c)
#   lognormal    as above (a = meanlog, b = sdlog)
#   gamma        S0 = pgamma(t, k, r, lower.tail = F)
#   pwexp        piecewise-constant hazard, cuts at 0.5, 1, 2, 3, 5, 10, 20 years (closest to Cox)
cuts <- c(0, 0.5, 1, 2, 3, 5, 10, 20, Inf)

base <- list(
  exponential = list(names = "loglambda",
                     logh = function(p, t) rep(p[1], length(t)),
                     H = function(p, t) exp(p[1]) * t),
  weibull = list(names = c("logc", "logb"),
                 logh = function(p, t) { c <- exp(p[1]); b <- exp(p[2]); log(c) - log(b) + (c-1)*(log(t) - log(b)) },
                 H = function(p, t) (t/exp(p[2]))^exp(p[1])),
  gompertz = list(names = c("g", "logr"),
                  logh = function(p, t) p[2] + p[1]*t,
                  H = function(p, t) if (abs(p[1]) < 1e-8) exp(p[2])*t else exp(p[2]) * expm1(p[1]*t) / p[1]),
  loglogistic = list(names = c("logc", "logb"),
                     logh = function(p, t) { c <- exp(p[1]); b <- exp(p[2]); log(c) - log(b) + (c-1)*(log(t) - log(b)) - log1p((t/b)^c) },
                     H = function(p, t) log1p((t/exp(p[2]))^exp(p[1]))),
  lognormal = list(names = c("a", "logb"),
                   logh = function(p, t) dlnorm(t, p[1], exp(p[2]), log = T) - plnorm(t, p[1], exp(p[2]), lower.tail = F, log.p = T),
                   H = function(p, t) -plnorm(t, p[1], exp(p[2]), lower.tail = F, log.p = T)),
  gamma = list(names = c("logk", "logr"),
               logh = function(p, t) dgamma(t, exp(p[1]), exp(p[2]), log = T) - pgamma(t, exp(p[1]), exp(p[2]), lower.tail = F, log.p = T),
               H = function(p, t) -pgamma(t, exp(p[1]), exp(p[2]), lower.tail = F, log.p = T)),
  pwexp = list(names = paste0("loglambda", 1:(length(cuts)-1)),
               logh = function(p, t) p[findInterval(t, cuts, left.open = T, rightmost.closed = T)],
               H = function(p, t) {
                 # time spent in each piece: min(max(t - cut_j, 0), width_j), then sum of lambda_j * time
                 expo <- pmin(pmax(outer(t, cuts[-length(cuts)], "-"), 0),
                              matrix(diff(cuts), length(t), length(cuts)-1, byrow = T))
                 drop(expo %*% exp(p))
               })
)

# generic PH likelihood: het = "none", "normal" (GH quadrature, as above) or "gamma" (closed form)
PHLik <- function(par, time, talpha, event, x, dist, het = "normal", trunc = "ratio",
                  nodes = gh$nodes, weights = gh$weights) {
  k <- length(base[[dist]]$names)
  p <- par[1:k]
  if (het == "none") {
    beta <- par[(k+1):length(par)]
  } else {
    theta <- exp(par[k+1])   # sigma for "normal", variance for "gamma"
    beta <- par[(k+2):length(par)]
  }
  xb <- drop(x%*%beta)
  logh0t <- base[[dist]]$logh(p, time)
  H0t <- base[[dist]]$H(p, time)
  H0a <- base[[dist]]$H(p, talpha)

  if (het == "none") {
    ll <- event * (logh0t + xb) - exp(xb) * (H0t - H0a)
  } else if (het == "normal") {
    z <- sqrt(2) * theta * nodes
    lw <- log(weights) - 0.5*log(pi)
    eta <- outer(xb, z, "+")
    lnum <- event * (logh0t + eta) - exp(eta) * H0t
    lden <- -exp(eta) * H0a
    if (trunc == "mle") {
      ll <- rowLogSumExp(sweep(lnum - lden, 2, lw, "+"))
    } else {
      ll <- rowLogSumExp(sweep(lnum, 2, lw, "+")) - rowLogSumExp(sweep(lden, 2, lw, "+"))
    }
  } else if (het == "gamma") {
    # gamma frailty (mean 1, variance theta) integrates out exactly:
    #   S(t|x) = (1 + theta e^xb H0(t))^(-1/theta),  h(t|x) = h0(t) e^xb / (1 + theta e^xb H0(t))
    # (delayed-entry version: divide by the marginal S(talpha|x))
    ut <- log1p(theta * exp(xb) * H0t)
    ua <- log1p(theta * exp(xb) * H0a)
    ll <- event * (logh0t + xb - ut) - (ut - ua)/theta
  }
  if (any(!is.finite(ll))) return(NA)
  return(sum(ll))
}

# check: the generic function reproduces lnormPHLik and lnormPHLik.gh
p.chk <- c(0.6, log(1.3), log(0.5), 0.2)
c(lnormPHLik(p.chk[-3], time, talpha, event, x1), PHLik(p.chk[-3], time, talpha, event, x1, "lognormal", het = "none"))
c(lnormPHLik.gh(p.chk, time, talpha, event, x1, trunc = "ratio"), PHLik(p.chk, time, talpha, event, x1, "lognormal", het = "normal"))

# starting values for each baseline (rough, from the crude exponential rate)
lam <- sum(event) / sum(time - talpha)
base.start <- list(exponential = log(lam),
                   weibull = c(0, -log(lam)),
                   gompertz = c(-0.01, log(lam)),
                   loglogistic = c(0, -log(lam)),
                   lognormal = c(1, 0),
                   gamma = c(0, log(lam)),
                   pwexp = rep(log(lam), length(cuts)-1))

# no-heterogeneity fit for a baseline: BFGS from rough starts, then NR to polish
fit.none <- function(dist, xx) {
  st0 <- c(setNames(base.start[[dist]], base[[dist]]$names), setNames(rep(0, ncol(xx)), paste0("B", colnames(xx))))
  f0 <- maxLik(logLik = PHLik, start = st0, time = time, talpha = talpha, event = event, x = xx,
               dist = dist, het = "none", method = "BFGS", iterlim = 10000)
  maxLik(logLik = PHLik, start = coef(f0), time = time, talpha = talpha, event = event, x = xx,
         dist = dist, het = "none", method = "NR", iterlim = 10000)
}

# heterogeneity fit, starting from the no-heterogeneity estimates with sigma (or variance) = 0.5
fit.het <- function(f0, dist, xx, het, trunc = "ratio") {
  k <- length(base[[dist]]$names)
  st <- c(coef(f0)[1:k], logtheta = log(0.5), coef(f0)[-(1:k)])
  tryCatch(maxLik(logLik = PHLik, start = st, time = time, talpha = talpha, event = event, x = xx,
                  dist = dist, het = het, trunc = trunc, method = "BFGS", iterlim = 10000),
           error = function(e) NULL)
}

results <- list()
for (m in c("BMI only", "All covariates")) {
  xx <- if (m == "BMI only") x1 else x
  for (dist in names(base)) {
    cat(m, "-", dist, "\n")
    f0 <- fit.none(dist, xx)
    fits <- list("none" = f0,
                 "normal, MLE truncation" = fit.het(f0, dist, xx, "normal", "mle"),
                 "normal, delayed entry" = fit.het(f0, dist, xx, "normal", "ratio"),
                 "gamma, delayed entry" = fit.het(f0, dist, xx, "gamma", "ratio"))
    for (h in names(fits)) {
      f <- fits[[h]]
      if (is.null(f)) next
      est <- coef(f)
      theta <- if (h == "none") 0 else exp(est["logtheta"])
      # LL re-evaluated with 60 nodes, as a check on quadrature error
      ll60 <- if (grepl("^normal", h)) PHLik(est, time, talpha, event, xx, dist, het = "normal",
                                             trunc = ifelse(grepl("MLE", h), "mle", "ratio"),
                                             nodes = gh60$nodes, weights = gh60$weights) else as.numeric(logLik(f))
      se.bmi <- tryCatch(showfit(f)["Bbmi", "se"], error = function(e) NA)
      lr <- if (h == "none") c(LR = NA, p = NA) else lrtest.sigma(f0, f)
      results[[length(results)+1]] <- data.frame(
        model = m, baseline = dist, heterogeneity = h,
        npar = length(est), logLik = as.numeric(logLik(f)), logLik.gh60 = ll60,
        AIC = -2*as.numeric(logLik(f)) + 2*length(est),
        het.par = theta,   # sigma for "normal", variance for "gamma"
        HR.bmi = exp(est["Bbmi"]), lo = exp(est["Bbmi"] - 1.96*se.bmi), hi = exp(est["Bbmi"] + 1.96*se.bmi),
        LR = lr["LR"], p.LR = lr["p"], row.names = NULL)
    }
  }
}
results <- do.call(rbind, results)
print(results, digits = 4)
write.csv(results, "frailty_baseline_comparison.csv", row.names = F)

# Quadrature check: logLik.gh60 - logLik should be ~0. It drifts (by up to ~1 log-likelihood unit)
# only where sigma comes out large (> ~2, e.g. Weibull or gamma baselines with normal frailty).
# Plain Gauss-Hermite needs more nodes as sigma grows; refit with gh60 if one of those matters.
summary(results$logLik.gh60 - results$logLik)

# best model by AIC within each covariate set
results[order(results$model, results$AIC), ][!duplicated(results$model[order(results$model, results$AIC)]), ]
