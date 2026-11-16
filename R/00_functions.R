##########################################################################
## 00_functions.R
## Shared functions for the 2016 obesity and employment duration analysis,
## rebuilt in 2026. Every model is written out by hand and maximised with
## optim(); the survival package is only used to check the answers.
##
## Conventions used throughout:
##   time   = years since the job started (exit or censoring)
##   talpha = years since the job started when follow-up began (left truncation)
##   event  = 1 if the person left employment, 0 if censored
##   x      = design matrix (no intercept), beta = log hazard ratios
##   xb     = linear predictor, always built as  xb <- linpred(beta, x)
##   a person is at risk at time t if  talpha < t <= time
##########################################################################

library(statmod)   ## gauss.quad() for the Gauss-Hermite nodes and weights


##########################################################################
# 1) data and the linear predictor
##########################################################################

# the 2016 covariates, in the 2016 order. Reference groups: white, married, uninsured, poor.
# sex = 1 for men; bmi = 1 if the mean BMI over rounds 3 and 5 was above 30.
covars <- c("age", "sex", "black", "other", "wasmarried", "nevermarried",
            "publicinsur", "privateinsur", "lowinc", "midinc", "highinc", "bmi")

readData <- function(file = "Data/mepsempcsv.csv") {
  d <- read.csv(file)
  stopifnot(all(d$time > d$talpha))   # everyone enters before they leave
  d
}

# design matrix from a formula, without the intercept (the baseline hazard plays that role)
makeX <- function(vars, data) {
  f <- reformulate(vars)
  x <- model.matrix(f, data = data)[, -1, drop = F]
  x
}

# the linear predictor x %*% beta, as a plain vector
linpred <- function(beta, x) drop(x %*% beta)

# sum of w over everyone with v >= t, for each value of t (sorted cumulative sums, no loops).
# With v = time this is everyone still employed at t; with v = talpha it is everyone not yet
# in the survey at t. Their difference is the risk set {talpha < t <= time}.
sumAtOrAbove <- function(w, v, t) {
  w <- as.matrix(w)
  o <- order(v)
  cs <- apply(w[o, , drop = F], 2, function(col) rev(cumsum(rev(col))))
  cs <- rbind(matrix(cs, ncol = ncol(w)), 0)          # row n+1 = nobody left
  first <- findInterval(t, v[o], left.open = T) + 1    # first position with v >= t
  cs[first, , drop = F]
}

riskSum <- function(w, time, talpha, t) sumAtOrAbove(w, time, t) - sumAtOrAbove(w, talpha, t)


##########################################################################
# 2) optim wrapper: maximise a log-likelihood, return estimates, SEs, logLik
##########################################################################
# fnscale = -1 makes optim maximise. The Hessian optim returns is the Hessian of the
# log-likelihood itself, so the variance is solve(-H). drop = names of parameters to leave
# out of the variance (a frailty sigma sitting at its boundary, where the Hessian is singular).
fitML <- function(loglik, start, ..., gr = NULL, method = "BFGS", drop = NULL, maxit = 10000) {
  # ndeps: step for the numerical gradient. optim's default (1e-3) is too coarse for
  # coefficients like age (a step of 1e-3 moves xb by ~0.04) and stops BFGS early.
  fit <- optim(start, loglik, gr = gr, ..., method = method, hessian = T,
               control = list(fnscale = -1, maxit = maxit, reltol = 1e-14,
                              ndeps = rep(1e-5, length(start))))
  keep <- !(names(fit$par) %in% drop)
  vcov <- matrix(NA, length(fit$par), length(fit$par), dimnames = list(names(fit$par), names(fit$par)))
  vcov[keep, keep] <- solve(-fit$hessian[keep, keep])
  list(estimate = fit$par, se = sqrt(diag(vcov)), vcov = vcov, loglik = fit$value,
       convergence = fit$convergence, npar = length(fit$par), hessian = fit$hessian)
}

# estimates on the reporting scale: parameters estimated as log(b) come back as b, with a
# delta-method SE (se(b) = b * se(log b))
report <- function(fit) {
  est <- fit$estimate; se <- fit$se
  lg <- grepl("^log", names(est))
  est[lg] <- exp(fit$estimate[lg]); se[lg] <- est[lg] * fit$se[lg]
  names(est) <- names(se) <- sub("^log", "", names(est))
  cbind(est = est, se = se)
}

# likelihood ratio test of sigma = 0. sigma = 0 is on the edge of the parameter space, so the
# statistic is a 50:50 mix of chi2(0) and chi2(1): halve the usual chi2(1) p-value.
lrBoundary <- function(ll0, ll1) {
  lr <- max(0, 2 * (ll1 - ll0))
  c(LR = lr, p = 0.5 * pchisq(lr, df = 1, lower.tail = F))
}


##########################################################################
# 3) Kaplan-Meier with left truncation, and the log-rank test
##########################################################################
# n_k = #{j: talpha_j < t_k <= time_j},  S(t) = prod_{t_k <= t} (1 - d_k/n_k)
# The 2016 KMcalc() did this with loops; this version counts the risk set with sorted sums.
KMcalc <- function(talpha, time, event) {
  tk <- sort(unique(time))
  n.risk <- drop(riskSum(rep(1, length(time)), time, talpha, tk))
  n.event <- as.vector(tapply(event, factor(time, levels = tk), sum))
  n.censor <- as.vector(tapply(1 - event, factor(time, levels = tk), sum))
  surv <- cumprod(1 - n.event/n.risk)
  # Greenwood standard error
  se <- surv * sqrt(cumsum(n.event / (n.risk * (n.risk - n.event))))
  data.frame(time = tk, n.risk, n.event, n.censor, surv, se)
}

# log-rank test for two groups (group = 0/1), with left-truncated risk sets.
# At each event time: O - E for group 1 and the hypergeometric variance.
logrank <- function(talpha, time, event, group) {
  tk <- sort(unique(time[event == 1]))
  n <- drop(riskSum(rep(1, length(time)), time, talpha, tk))
  n1 <- drop(riskSum(as.numeric(group == 1), time, talpha, tk))
  d <- as.vector(tapply(event, factor(time, levels = tk), sum))
  d1 <- as.vector(tapply(event * (group == 1), factor(time, levels = tk), sum))
  d[is.na(d)] <- 0; d1[is.na(d1)] <- 0
  O1 <- sum(d1); E1 <- sum(d * n1/n)
  V <- sum(d * (n1/n) * (1 - n1/n) * ifelse(n > 1, (n - d)/(n - 1), 0))
  chisq <- (O1 - E1)^2 / V
  # without the (n - d)/(n - 1) tie correction this is the Cox score test with Breslow ties
  chisq.breslow <- (O1 - E1)^2 / sum(d * (n1/n) * (1 - n1/n))
  c(O1 = O1, E1 = E1, V = V, chisq = chisq, df = 1, p = pchisq(chisq, 1, lower.tail = F),
    chisq.breslow = chisq.breslow)
}


##########################################################################
# 4) Cox proportional hazards: partial likelihood with delayed entry, Breslow ties
##########################################################################
# PL = prod_i [ e^{x_i b} / sum_{j in R(t_i)} e^{x_j b} ]^{delta_i},  R(t) = {j: talpha_j < t <= time_j}
# Tied exits all share the same risk set (Breslow). This is the coxPHLik() from 2016, which
# never ran (riskset was not set up and x came from an undefined data frame).
coxPHLik <- function(par, time, talpha, event, x) {
  xb <- linpred(par, x)
  tf <- time[event == 1]
  S0 <- drop(riskSum(exp(xb), time, talpha, tf))
  sum(xb[event == 1] - log(S0))
}

# its gradient (the score): sum over exits of x_i minus the risk-set average of x
coxPHScore <- function(par, time, talpha, event, x) {
  w <- exp(linpred(par, x))
  tf <- time[event == 1]
  S0 <- drop(riskSum(w, time, talpha, tf))
  S1 <- riskSum(w * x, time, talpha, tf)
  colSums(x[event == 1, , drop = F] - S1/S0)
}

# Breslow estimate of the baseline cumulative hazard, H0(t) = sum_{t_i <= t} delta_i / S0(t_i),
# returned as a step function
coxBaseHaz <- function(beta, time, talpha, event, x) {
  w <- exp(linpred(beta, x))
  tf <- sort(time[event == 1])
  S0 <- drop(riskSum(w, time, talpha, tf))
  H <- cumsum(1/S0)
  keep <- !duplicated(tf, fromLast = T)
  stepfun(tf[keep], c(0, H[keep]))
}


##########################################################################
# 5) residuals (for any PH model with a cumulative hazard function)
##########################################################################
# Cox-Snell: r_i = e^{xb_i} [H0(time_i) - H0(talpha_i)], the hazard accumulated while observed.
#            If the model is right, these behave like censored unit exponentials.
# martingale: M_i = delta_i - r_i
# deviance:   sign(M_i) sqrt(-2 [M_i + delta_i log(delta_i - M_i)])
phResiduals <- function(xb, H0, time, talpha, event) {
  r <- exp(xb) * (H0(time) - H0(talpha))
  M <- event - r
  dev <- sign(M) * sqrt(-2 * (M + ifelse(event == 1, log(event - M), 0)))
  data.frame(coxsnell = r, martingale = M, deviance = dev)
}

# Nelson-Aalen cumulative hazard of the Cox-Snell residuals (should follow the 45-degree line)
coxSnellNA <- function(r, event) {
  km <- KMcalc(talpha = rep(-Inf, length(r)), time = r, event = event)   # everyone at risk from the start
  data.frame(r = km$time, H = cumsum(km$n.event / km$n.risk))
}

# Schoenfeld residuals for the Cox model: at each exit, r_i = x_i minus the risk-set average of x.
# Scaled (Grambsch-Therneau) residuals: beta + d * V r_i; smoothed against time they show beta(t).
#
# Test of proportional hazards: add x_k * g(t), g(t) = log(t), to the model and test its
# coefficient gamma_k = 0 with a score test at the fitted beta:
#   U_k = sum_i g(t_i) r_ik
#   Var(U) = I_gg - I_gb I_bb^-1 I_bg,  I_gg = sum_i g_i^2 V_i, I_gb = sum_i g_i V_i, I_bb = sum_i V_i
# where V_i is the covariance of x in the risk set at t_i. One covariate at a time, then all
# together (GLOBAL). This is the test cox.zph() uses since survival 3.0. The 2016-era cox.zph()
# used a quicker approximation (V_i replaced by its average); that is returned as test.approx.
schoenfeld <- function(beta, vcov, time, talpha, event, x) {
  w <- exp(linpred(beta, x))
  tf <- time[event == 1]
  S0 <- drop(riskSum(w, time, talpha, tf))
  S1 <- riskSum(w * x, time, talpha, tf)
  r <- x[event == 1, , drop = F] - S1/S0
  o <- order(tf)
  r <- r[o, , drop = F]; tf <- tf[o]; S0 <- S0[o]; S1 <- S1[o, , drop = F]
  d <- nrow(r); p <- ncol(x)
  g <- log(tf)

  # risk-set covariance of x at each exit: S2/S0 - (S1/S0)(S1/S0)'
  pairs <- which(upper.tri(diag(p), diag = T), arr.ind = T)
  S2 <- riskSum(w * x[, pairs[, 1]] * x[, pairs[, 2]], time, talpha, tf)   # tf is already sorted
  Igg <- Igb <- Ibb <- matrix(0, p, p)
  for (i in 1:d) {
    Vi <- matrix(0, p, p)
    Vi[pairs] <- S2[i, ] / S0[i]
    Vi[pairs[, 2:1]] <- S2[i, ] / S0[i]
    m <- S1[i, ] / S0[i]
    Vi <- Vi - tcrossprod(m)
    Ibb <- Ibb + Vi; Igb <- Igb + g[i] * Vi; Igg <- Igg + g[i]^2 * Vi
  }
  U <- drop(g %*% r)
  VU <- Igg - Igb %*% solve(Ibb, t(Igb))
  chisq <- U^2 / diag(VU)
  global <- drop(U %*% solve(VU, U))
  ph.test <- rbind(cbind(chisq = chisq, df = 1, p = pchisq(chisq, 1, lower.tail = F)),
                   GLOBAL = c(global, p, pchisq(global, p, lower.tail = F)))

  # scaled residuals and the older approximate test (g centred)
  scaled <- sweep(d * r %*% vcov, 2, beta, "+")
  gc <- g - mean(g)
  ch.a <- drop(gc %*% (d * r %*% vcov))^2 / (diag(vcov) * d * sum(gc^2))
  gl.a <- drop(U %*% vcov %*% U) * d / sum(gc^2)
  test.approx <- rbind(cbind(chisq = ch.a, df = 1, p = pchisq(ch.a, 1, lower.tail = F)),
                       GLOBAL = c(gl.a, p, pchisq(gl.a, p, lower.tail = F)))

  colnames(r) <- colnames(scaled) <- rownames(ph.test)[1:p] <- rownames(test.approx)[1:p] <- colnames(x)
  list(time = tf, resid = r, scaled = scaled, test = ph.test, test.approx = test.approx)
}


##########################################################################
# 6) parametric baselines: each needs log h0(t) and H0(t) = -log S0(t)
##########################################################################
# Parameters are on unconstrained scales (logs for the positive ones).
#   exponential  h0 = lambda
#   weibull      S0 = exp(-(t/b)^c)                        (the 2016 weibLik parameterisation)
#   lognormal    S0 = 1 - Phi((log t - a)/b)               (a = meanlog, b = sdlog, as in 2016)
#   gamma        S0 = 1 - pgamma(t, k, r)                  (shape k, rate r)
#   gompertz     h0 = r exp(g t), H0 = (r/g)(exp(g t) - 1)  (g can be negative)
#   loglogistic  S0 = 1/(1 + (t/b)^c)
base <- list(
  exponential = list(names = "loglambda",
                     logh = function(p, t) rep(p[1], length(t)),
                     H = function(p, t) exp(p[1]) * t),
  weibull = list(names = c("logc", "logb"),
                 logh = function(p, t) { c <- exp(p[1]); b <- exp(p[2]); log(c) - log(b) + (c-1)*(log(t) - log(b)) },
                 H = function(p, t) (t/exp(p[2]))^exp(p[1])),
  lognormal = list(names = c("a", "logb"),
                   logh = function(p, t) dlnorm(t, p[1], exp(p[2]), log = T) - plnorm(t, p[1], exp(p[2]), lower.tail = F, log.p = T),
                   H = function(p, t) -plnorm(t, p[1], exp(p[2]), lower.tail = F, log.p = T)),
  gamma = list(names = c("logk", "logr"),
               logh = function(p, t) dgamma(t, exp(p[1]), exp(p[2]), log = T) - pgamma(t, exp(p[1]), exp(p[2]), lower.tail = F, log.p = T),
               H = function(p, t) -pgamma(t, exp(p[1]), exp(p[2]), lower.tail = F, log.p = T)),
  gompertz = list(names = c("g", "logr"),
                  logh = function(p, t) p[2] + p[1]*t,
                  H = function(p, t) if (abs(p[1]) < 1e-8) exp(p[2])*t else exp(p[2]) * expm1(p[1]*t) / p[1]),
  loglogistic = list(names = c("logc", "logb"),
                     logh = function(p, t) { c <- exp(p[1]); b <- exp(p[2]); log(c) - log(b) + (c-1)*(log(t) - log(b)) - log1p((t/b)^c) },
                     H = function(p, t) log1p((t/exp(p[2]))^exp(p[1])))
)

# survival curve S0(t) for plotting
survCurve <- function(dist, p, t) exp(-base[[dist]]$H(p, t))


##########################################################################
# 7) parametric proportional hazards likelihood, with or without frailty
##########################################################################
# h(t | x, z) = h0(t) exp(xb + z). Left truncation: each person's likelihood is divided by the
# chance of still being employed when follow-up began.
#
# No frailty:
#   log L_i = delta_i [log h0(t_i) + xb_i] - e^{xb_i} [H0(t_i) - H0(talpha_i)]
#   (the 2016 lnormPHLik, written on the log scale so nothing underflows)
#
# Normal frailty z ~ N(0, sigma), integrated out with Gauss-Hermite quadrature:
#   integral g(z) f(u) dz ~= sum_k [w_k / sqrt(pi)] f(sqrt(2) sigma u_k)
#   trunc = "delayed": z is fixed from the start of the job, so the people still employed at
#           entry are a selected group. Divide the whole integral by S(talpha) averaged over z:
#           L_i = int g(z) h(t_i|z)^delta S(t_i|z) dz / int g(z) S(talpha_i|z) dz
#   trunc = "2016":    divide by S(talpha|z) inside the integral, as the 2016 poster and mle did
#           (as if z were drawn when the person entered the survey).
#   The two agree when sigma = 0.
#   sigma enters only through z = sigma * u with nodes symmetric about 0, so the likelihood is
#   the same for sigma and -sigma; it is estimated unconstrained and reported as |sigma|.
gh <- gauss.quad(30, kind = "hermite")

# log(sum_k exp(M[, k])) for every row of M, without underflow
rowLogSumExp <- function(M) {
  m <- M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
  m + log(rowSums(exp(M - m)))
}

PHLik <- function(par, time, talpha, event, x = NULL, dist = "lognormal", frailty = F,
                  trunc = "delayed", nodes = gh$nodes, weights = gh$weights, sum = T) {
  k <- length(base[[dist]]$names)
  p <- par[1:k]
  sigma <- if (frailty) par[k+1] else 0
  beta <- if (is.null(x)) numeric(0) else par[(k + frailty + 1):length(par)]
  xb <- if (is.null(x)) rep(0, length(time)) else linpred(beta, x)

  logh0t <- base[[dist]]$logh(p, time)
  H0t <- base[[dist]]$H(p, time)
  H0a <- base[[dist]]$H(p, talpha)

  if (!frailty) {
    ll <- event * (logh0t + xb) - exp(xb) * (H0t - H0a)
  } else {
    z <- sqrt(2) * sigma * nodes
    lw <- log(weights) - 0.5*log(pi)
    eta <- outer(xb, z, "+")                          # N x K matrix of xb + z_k
    lnum <- event * (logh0t + eta) - exp(eta) * H0t   # log of h(t|z)^delta S(t|z)
    lden <- -exp(eta) * H0a                           # log of S(talpha|z)
    if (trunc == "2016") {
      ll <- rowLogSumExp(sweep(lnum - lden, 2, lw, "+"))
    } else {
      ll <- rowLogSumExp(sweep(lnum, 2, lw, "+")) - rowLogSumExp(sweep(lden, 2, lw, "+"))
    }
  }
  if (!sum) return(ll)
  if (any(!is.finite(ll))) return(-1e10)   # keeps optim away from impossible values
  sum(ll)
}

# starting values for a baseline, from the crude exit rate
baseStart <- function(dist, time, talpha, event) {
  lam <- sum(event) / sum(time - talpha)
  st <- switch(dist,
               exponential = log(lam), weibull = c(0, -log(lam)), lognormal = c(1, 0),
               gamma = c(0, log(lam)), gompertz = c(-0.01, log(lam)), loglogistic = c(0, -log(lam)))
  setNames(st, base[[dist]]$names)
}

# fit a PH model: no frailty first, then (if asked) frailty starting from those estimates
fitPH <- function(dist, time, talpha, event, x = NULL, frailty = F, trunc = "delayed",
                  sigma.start = 0.5, start = NULL, ...) {
  if (is.null(start)) {
    st0 <- baseStart(dist, time, talpha, event)
    if (!is.null(x)) st0 <- c(st0, setNames(rep(0, ncol(x)), paste0("B", colnames(x))))
    f0 <- fitML(PHLik, st0, time = time, talpha = talpha, event = event, x = x, dist = dist, ...)
    if (!frailty) return(f0)
    k <- length(base[[dist]]$names)
    start <- c(f0$estimate[1:k], sigma = sigma.start, f0$estimate[-(1:k)])
  }
  f1 <- fitML(PHLik, start, time = time, talpha = talpha, event = event, x = x, dist = dist,
              frailty = T, trunc = trunc, ...)
  f1$estimate["sigma"] <- abs(f1$estimate["sigma"])
  # at the boundary (sigma ~ 0) the Hessian is singular in sigma: refit the variance without it
  if (f1$estimate["sigma"] < 0.01) {
    keep <- names(f1$estimate) != "sigma"
    f1$vcov[] <- NA
    f1$vcov[keep, keep] <- solve(-f1$hessian[keep, keep])
    f1$se <- sqrt(diag(f1$vcov))
  }
  f1
}
