##########################################################################
## 05_frailty.R
## Log-normal PH model with normal unobserved heterogeneity (frailty), integrated
## with Gauss-Hermite quadrature, and tests of whether the frailty is needed:
##  1) quadrature check against integrate(), person by person
##  2) fits: BMI only and all covariates; delayed entry (main) and the 2016/mle form
##  3) likelihood ratio tests of sigma = 0 (boundary corrected), and a 60-node check
##  4) profile likelihood for sigma
##  5) other baselines x (no, normal, gamma) frailty   [about 15 minutes; set run.baselines]
## Needs output/02_parametric.rds.
## Run from the repository root:  source("R/05_frailty.R")
##########################################################################
rm(list=ls())
source("R/00_functions.R")
run.baselines <- T

# prepare the data
d <- readData()
talpha <- d$talpha
time <- d$time
event <- d$event
x1 <- makeX("bmi", d)    # BMI only
x <- makeX(covars, d)    # all covariates

pm <- readRDS("output/02_parametric.rds")
gh60 <- gauss.quad(60, kind = "hermite")


##########################################################################
# 1) does 30-node Gauss-Hermite get each person's integral right?
##########################################################################
# delayed entry: log L_i = log int g(z) h(t|z)^delta S(t|z) dz - log int g(z) S(talpha|z) dz
check.gh <- function(par, i, x) {
  a <- par[1]; b <- exp(par[2]); sigma <- par[3]; beta <- par[4:length(par)]
  xb <- linpred(beta, x[i, , drop = F])
  logS0t <- plnorm(time[i], a, b, lower.tail = F, log.p = T)
  logS0a <- plnorm(talpha[i], a, b, lower.tail = F, log.p = T)
  logh0t <- dlnorm(time[i], a, b, log = T) - logS0t
  num <- function(z) dnorm(z, 0, sigma) * exp(event[i]*(logh0t + xb + z) + exp(xb + z)*logS0t)
  den <- function(z) dnorm(z, 0, sigma) * exp(exp(xb + z)*logS0a)
  c(integrate = log(integrate(num, -12*sigma, 12*sigma, rel.tol = 1e-10)$value) -
                log(integrate(den, -12*sigma, 12*sigma, rel.tol = 1e-10)$value),
    gh30 = PHLik(par, time[i], talpha[i], event[i], x[i, , drop = F], frailty = T, sum = F))
}
ids <- c(which(event == 1)[1:3], which(event == 0 & talpha > 5)[1:3])
for (sig in c(0.5, 1)) {
  chk <- sapply(ids, check.gh, par = c(0.57, log(1.29), sig, 0.18), x = x1)
  cat("check GH vs integrate(), sigma =", sig, ": max abs difference", max(abs(chk[1, ] - chk[2, ])), "\n")
}


##########################################################################
# 2) fits
##########################################################################
# no frailty (from 02) and normal frailty with each truncation form, starting at sigma = 0.5
fr <- list()
for (m in c("BMI only", "All covariates")) {
  xx <- if (m == "BMI only") x1 else x
  f0 <- if (m == "BMI only") pm$lnormPH.fit1 else pm$lnormPH.fit2
  st <- c(f0$estimate[1:2], sigma = 0.5, f0$estimate[-(1:2)])
  fr[[m]] <- list(none = f0,
                  delayed = fitPH("lognormal", time, talpha, event, xx, frailty = T, trunc = "delayed", start = st),
                  y2016 = fitPH("lognormal", time, talpha, event, xx, frailty = T, trunc = "2016", start = st))
}


##########################################################################
# 3) tests of sigma = 0
##########################################################################
frailty.tests <- do.call(rbind, lapply(names(fr), function(m) {
  do.call(rbind, lapply(c("delayed", "y2016"), function(tr) {
    f <- fr[[m]][[tr]]
    xx <- if (m == "BMI only") x1 else x
    ll60 <- PHLik(f$estimate, time, talpha, event, xx, frailty = T, trunc = ifelse(tr == "delayed", "delayed", "2016"),
                  nodes = gh60$nodes, weights = gh60$weights)
    data.frame(model = m, truncation = ifelse(tr == "delayed", "delayed entry", "2016 (inside the integral)"),
               sigma = f$estimate["sigma"], logLik.none = fr[[m]]$none$loglik, logLik.frailty = f$loglik,
               logLik.gh60 = ll60, t(lrBoundary(fr[[m]]$none$loglik, f$loglik)),
               HR.bmi = exp(f$estimate["Bbmi"]), se.bmi = f$se["Bbmi"], row.names = NULL)
  }))
}))
print(frailty.tests, digits = 6)


##########################################################################
# 4) profile likelihood for sigma (all covariates, delayed entry)
##########################################################################
# Fix sigma on a grid and maximise over everything else. If the profile is highest at
# sigma = 0, the boundary estimate is the maximum, not an optimiser that got stuck.
sig.grid <- c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5)
prof <- numeric(length(sig.grid))
st <- pm$lnormPH.fit2$estimate
for (j in seq_along(sig.grid)) {
  f <- optim(st, function(p) PHLik(c(p[1:2], sig.grid[j], p[-(1:2)]), time, talpha, event, x, frailty = T,
                                   nodes = gh60$nodes, weights = gh60$weights),
             method = "BFGS", control = list(fnscale = -1, maxit = 10000, reltol = 1e-12, ndeps = rep(1e-5, length(st))))
  prof[j] <- f$value
  st <- f$par   # warm start the next sigma from this one
}
profile <- data.frame(sigma = c(0, sig.grid), logLik = c(pm$lnormPH.fit2$loglik, prof))
print(profile, digits = 8)

saveRDS(list(fr = fr, frailty.tests = frailty.tests, profile = profile), "output/05_frailty.rds")
write.csv(frailty.tests, "output/frailty_tests.csv", row.names = F)
write.csv(profile, "output/frailty_profile.csv", row.names = F)


##########################################################################
# 5) other baselines, with no, normal or gamma frailty (delayed entry)
##########################################################################
# Same PH structure h(t|x,z) = h0(t) exp(xb + z) for every baseline. Gamma frailty
# (mean 1, variance theta) integrates out exactly:
#   S(t|x) = (1 + theta e^xb H0(t))^(-1/theta),  h(t|x) = h0(t) e^xb / (1 + theta e^xb H0(t))
# divided by the marginal S(talpha|x) for delayed entry. theta is estimated as log(theta).
gammaFrailtyLik <- function(par, time, talpha, event, x, dist) {
  k <- length(base[[dist]]$names)
  p <- par[1:k]
  theta <- exp(par[k+1])
  xb <- linpred(par[(k+2):length(par)], x)
  ut <- log1p(theta * exp(xb) * base[[dist]]$H(p, time))
  ua <- log1p(theta * exp(xb) * base[[dist]]$H(p, talpha))
  ll <- event * (base[[dist]]$logh(p, time) + xb - ut) - (ut - ua)/theta
  if (any(!is.finite(ll))) return(-1e10)
  sum(ll)
}

if (run.baselines) {
  t0 <- Sys.time()
  results <- list()
  for (m in c("BMI only", "All covariates")) {
    xx <- if (m == "BMI only") x1 else x
    for (dist in names(base)) {
      cat(m, "-", dist, "\n")
      f0 <- fitPH(dist, time, talpha, event, xx)
      k <- length(base[[dist]]$names)
      fn <- tryCatch(fitPH(dist, time, talpha, event, xx, frailty = T,
                           start = c(f0$estimate[1:k], sigma = 0.5, f0$estimate[-(1:k)])), error = function(e) NULL)
      fg <- tryCatch(fitML(gammaFrailtyLik, c(f0$estimate[1:k], logtheta = log(0.5), f0$estimate[-(1:k)]),
                           time = time, talpha = talpha, event = event, x = xx, dist = dist,
                           drop = "logtheta"), error = function(e) NULL)
      fits <- list("none" = f0, "normal" = fn, "gamma" = fg)
      for (h in names(fits)) {
        f <- fits[[h]]
        if (is.null(f)) next
        het <- switch(h, none = 0, normal = f$estimate["sigma"], gamma = exp(f$estimate["logtheta"]))
        lr <- if (h == "none") c(LR = NA, p = NA) else lrBoundary(f0$loglik, f$loglik)
        results[[length(results)+1]] <- data.frame(
          model = m, baseline = dist, frailty = h, npar = f$npar, logLik = f$loglik,
          AIC = -2*f$loglik + 2*f$npar,
          het.par = unname(het),   # sigma (normal) or variance theta (gamma)
          HR.bmi = exp(f$estimate["Bbmi"]),
          lo = exp(f$estimate["Bbmi"] - 1.96*f$se["Bbmi"]), hi = exp(f$estimate["Bbmi"] + 1.96*f$se["Bbmi"]),
          LR = lr["LR"], p.LR = lr["p"], row.names = NULL)
      }
    }
  }
  results <- do.call(rbind, results)
  print(results, digits = 4)
  write.csv(results, "output/frailty_baseline_comparison.csv", row.names = F)
  print(Sys.time() - t0)

  # best model by AIC within each covariate set
  print(results[order(results$model, results$AIC), ][!duplicated(results$model[order(results$model, results$AIC)]), ])
}
