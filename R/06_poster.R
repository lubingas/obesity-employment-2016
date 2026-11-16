##########################################################################
## 06_poster.R
## Figures (tikz, so the text is set in the poster's fonts) and Table 1 for
## prj_poster.tex, from the saved results of 01-05.
## Run from the repository root:  source("R/06_poster.R")
##########################################################################
rm(list=ls())
source("R/00_functions.R")
library(tikzDevice)
library(splines)

km <- readRDS("output/01_km.rds")
pm <- readRDS("output/02_parametric.rds")
cx <- readRDS("output/03_cox.rds")
dg <- readRDS("output/04_diagnostics.rds")
fr <- readRDS("output/05_frailty.rds")
d <- readData()

# Okabe-Ito colours
col <- c(normbmi = "#0072B2", highbmi = "#D55E00", cox = "#009E73", lnorm = "#E69F00",
         grey = "#7F7F7F", sky = "#56B4E9", purple = "#CC79A7")

# tikz set up to measure text with the poster's label font (Inter), so labels fit exactly
fontdir <- normalizePath("fonts")
options(tikzDefaultEngine = "xetex",
        tikzXelatexPackages = c("\\usepackage{tikz}",
                                "\\usepackage[active,tightpage,xetex]{preview}",
                                "\\usepackage{fontspec}",
                                paste0("\\setmainfont{Inter-Regular.otf}[Path=", fontdir, "/, ItalicFont=Inter-Italic.otf]"),
                                "\\usepackage{unicode-math}",
                                "\\setmathfont{Latin Modern Math}",
                                paste0("\\setmathfont{SourceSerif4-It.otf}[Path=", fontdir, "/, range=it/{latin,Latin,greek,Greek}]"),
                                "\\PreviewEnvironment{pgfpicture}",
                                "\\setlength\\PreviewBorder{0pt}"),
        tikzMetricsDictionary = file.path(tempdir(), "tikzMetrics"))

# common plot set up: no box, light axes, small margins
setpar <- function(mar = c(3.2, 3.6, 0.8, 0.6)) par(mar = mar, mgp = c(2.1, 0.55, 0), tcl = -0.25, las = 1,
                                                    cex.axis = 0.85, col.axis = "#3b3a35", fg = "#5c5b55")
times <- seq(0.01, 42, length.out = 400)

# Dense point clouds (thousands of residuals) would make the tikz files too big for TeX's memory,
# so the points are drawn into a 600 dpi PNG that fills the plot region exactly; axes, lines and
# text stay as tikz. Call after plot(NA, ...) has set up the axes.
rasterPoints <- function(x, y, col, cex = 0.25, pch = 16, res = 600) {
  usr <- par("usr"); pin <- par("pin"); lg <- par("xlog")
  xl <- if (lg) 10^usr[1:2] else usr[1:2]
  f <- tempfile(fileext = ".png")
  png(f, width = round(pin[1]*res), height = round(pin[2]*res), res = res, bg = "transparent")
  par(mar = c(0, 0, 0, 0))
  plot.new(); plot.window(xl, usr[3:4], xaxs = "i", yaxs = "i", log = ifelse(lg, "x", ""))
  points(x, y, pch = pch, cex = cex, col = col)
  dev.off()
  dev.set(dev.prev())   # back to the tikz device
  rasterImage(as.raster(png::readPNG(f)), xl[1], usr[3], xl[2], usr[4], interpolate = T)
}

# tikzDevice writes the PNG next to the .tex but refers to it without the folder; the poster
# compiles from the repository root, so add the folder to the reference
fixRaster <- function(file) {
  tx <- readLines(file)
  tx <- sub("]{fig_", "]{figures/fig_", tx, fixed = T)
  writeLines(tx, file)
}

# thin a step function to points at least 'by' apart in x (for long Cox-Snell curves)
thin <- function(x, y, by) { k <- !duplicated(floor(x/by)) | seq_along(x) == length(x); list(x = x[k], y = y[k]) }


##########################################################################
# Figure 1: Kaplan-Meier by BMI, with log-normal and Weibull curves fitted to each group
##########################################################################
tikz("figures/fig_km.tex", width = 7.05, height = 4.6, pointsize = 11)
setpar()
plot(NA, xlim = c(0, 42), ylim = c(0, 1), xaxs = "i", yaxs = "i", frame.plot = F, xaxt = "n",
     xlab = "Years since the job started", ylab = "$S(t)$, still employed")
axis(1, at = seq(0, 40, 10))
abline(h = seq(0.25, 1, 0.25), col = "#e2e1dc", lwd = 0.6)
for (g in c("normbmi", "highbmi")) {
  k <- km[[g]]
  lines(c(0, k$time), c(1, k$surv), type = "s", col = col[g], lwd = 1.8)
  i <- ifelse(g == "normbmi", 0, 1)
  lines(times, survCurve("lognormal", pm$bygroup[[paste0("lognormal.", i)]]$estimate, times), col = col[g], lwd = 1.4, lty = 2)
  lines(times, survCurve("weibull", pm$bygroup[[paste0("weibull.", i)]]$estimate, times), col = col[g], lwd = 1.4, lty = 3)
}
lr <- km$logrank
legend("topright", bty = "n", cex = 0.85, xpd = NA, lwd = c(1.8, 1.8, NA, 1.8, 1.4, 1.4), lty = c(1, 1, NA, 1, 2, 3),
       col = c(col["normbmi"], col["highbmi"], NA, "#3b3a35", "#3b3a35", "#3b3a35"),
       legend = c("BMI $\\leq 30$", "BMI $> 30$", "", "Kaplan--Meier", "Log-normal", "Weibull"))
text(42, 0.37, adj = 1, cex = 0.85, col = "#3b3a35",
     sprintf("Log-rank (left-truncated risk sets): $\\chi^2_1 = %.2f$, $p = %.3f$", lr["chisq"], lr["p"]))
dev.off()


##########################################################################
# Figure 2a: log(-log S(t)) against log(t), by BMI
##########################################################################
tikz("figures/fig_phcheck.tex", width = 2.33, height = 2.75, pointsize = 11)
setpar(mar = c(3.2, 3.3, 0.8, 0.4))
pc <- km$phcheck
plot(NA, xlim = c(0.05, 42), ylim = range(c(pc$normbmi$lls, pc$highbmi$lls)), log = "x", frame.plot = F,
     xaxt = "n", xlab = "Years (log scale)", ylab = "$\\log(-\\log \\hat{S}(t))$")
axis(1, at = c(0.05, 0.5, 5, 40), labels = c("0.05", "0.5", "5", "40"))
for (g in c("normbmi", "highbmi")) lines(pc[[g]]$time, pc[[g]]$lls, type = "s", col = col[g], lwd = 1.5)
legend("topleft", bty = "n", cex = 0.75, lwd = 1.5, col = col[c("normbmi", "highbmi")],
       legend = c("BMI $\\leq 30$", "BMI $> 30$"))
dev.off()


##########################################################################
# Figure 2b: scaled Schoenfeld residuals for BMI, against log(time)
##########################################################################
sch <- dg$sch
y <- sch$scaled[, "bmi"]; lt <- log(sch$time)
sm <- lm(y ~ ns(lt, df = 4))
pt <- seq(min(lt), max(lt), length.out = 200)
pr <- predict(sm, newdata = data.frame(lt = pt), se.fit = T)
tikz("figures/fig_schoenfeld.tex", width = 2.33, height = 2.75, pointsize = 11)
setpar(mar = c(3.2, 3.3, 0.8, 0.4))
plot(NA, xlim = range(sch$time), ylim = c(-3, 3.5), log = "x", frame.plot = F, xaxt = "n",
     xlab = "Years (log scale)", ylab = "$\\beta_{\\mathrm{BMI}}(t)$")
axis(1, at = c(0.05, 0.5, 5, 40), labels = c("0.05", "0.5", "5", "40"))
rasterPoints(sch$time, y, col = adjustcolor(col["grey"], 0.4), cex = 0.3)
polygon(exp(c(pt, rev(pt))), c(pr$fit + 2*pr$se.fit, rev(pr$fit - 2*pr$se.fit)), border = NA,
        col = adjustcolor(col["highbmi"], 0.18))
lines(exp(pt), pr$fit, col = col["highbmi"], lwd = 1.8)
abline(h = cx$cox.fit2$estimate["bmi"], lty = 2, col = "#3b3a35")
abline(h = 0, lty = 3, col = "#5c5b55")
dev.off()
fixRaster("figures/fig_schoenfeld.tex")


##########################################################################
# Figure 3a: Cox-Snell residuals, Cox and log-normal PH
##########################################################################
tikz("figures/fig_coxsnell.tex", width = 2.33, height = 2.75, pointsize = 11)
setpar(mar = c(3.2, 3.3, 0.8, 0.4))
lim <- c(0, 1.5)
plot(NA, xlim = lim, ylim = lim, frame.plot = F, xlab = "Cox--Snell residual $r$", ylab = "$\\hat{H}(r)$")
abline(0, 1, col = "#5c5b55", lty = 3)
# only the plotted range, thinned to steps at least 0.002 apart
cs1 <- with(dg$cs.cox[dg$cs.cox$r <= lim[2], ], thin(r, H, 0.002))
cs2 <- with(dg$cs.ln[dg$cs.ln$r <= lim[2], ], thin(r, H, 0.002))
lines(cs1$x, cs1$y, type = "s", col = col["cox"], lwd = 1.5)
lines(cs2$x, cs2$y, type = "s", col = col["lnorm"], lwd = 1.5, lty = 2)
legend("bottomright", bty = "n", cex = 0.75, xpd = NA, lwd = 1.5, lty = c(1, 2), col = col[c("cox", "lnorm")],
       legend = c("Cox", "Log-normal"))
dev.off()


##########################################################################
# Figure 3b: martingale residuals against age (functional form)
##########################################################################
tikz("figures/fig_martingale.tex", width = 2.33, height = 2.75, pointsize = 11)
setpar(mar = c(3.2, 3.3, 0.8, 0.4))
set.seed(544)
plot(NA, xlim = range(d$age) + c(-0.5, 0.5), ylim = c(-2.5, 1), frame.plot = F,
     xlab = "Age (years)", ylab = "Martingale residual")
rasterPoints(d$age + runif(nrow(d), -0.4, 0.4), dg$res.cox$martingale, col = adjustcolor(col["grey"], 0.3))
lines(dg$mart.age, col = col["highbmi"], lwd = 1.8)
abline(h = 0, lty = 3, col = "#5c5b55")
dev.off()
fixRaster("figures/fig_martingale.tex")


##########################################################################
# Figure 3c: deviance residuals against the linear predictor
##########################################################################
tikz("figures/fig_deviance.tex", width = 2.33, height = 2.75, pointsize = 11)
setpar(mar = c(3.2, 3.3, 0.8, 0.4))
dl <- dg$dev.lp
plot(NA, xlim = range(dl$lp), ylim = c(-2.5, 3.5), frame.plot = F,
     xlab = "Linear predictor $x\\hat{\\beta}$", ylab = "Deviance residual")
rasterPoints(dl$lp, dl$deviance,
             col = ifelse(dl$event == 1, adjustcolor(col["highbmi"], 0.4), adjustcolor(col["normbmi"], 0.25)))
abline(h = 0, lty = 3, col = "#5c5b55")
legend("topright", bty = "n", cex = 0.75, xpd = NA, pch = 16, col = col[c("highbmi", "normbmi")], legend = c("Exit", "Censored"))
dev.off()
fixRaster("figures/fig_deviance.tex")


##########################################################################
# Figure 4: profile log-likelihood for sigma (all covariates, delayed entry)
##########################################################################
# shows sigma up to 0.75; the full grid (to 1.5) keeps falling and is in output/frailty_profile.csv
pf <- fr$profile
pz <- pf[pf$sigma <= 0.75, ]
tikz("figures/fig_profile.tex", width = 2.33, height = 2.75, pointsize = 11)
setpar(mar = c(3.2, 3.3, 0.8, 0.4))
plot(pz$sigma, pz$logLik - pf$logLik[1], type = "b", pch = 16, cex = 0.8, col = col["normbmi"], lwd = 1.5,
     frame.plot = F, ylim = c(-13, 0.5), xlab = "Frailty SD $\\sigma$ (fixed)", ylab = "Profile $\\ell(\\sigma) - \\ell(0)$")
abline(h = -qchisq(0.90, 1)/2, lty = 3, col = "#5c5b55")
text(0.45, -qchisq(0.90, 1)/2, "95\\% limit", pos = 3, cex = 0.8, col = "#5c5b55", offset = 0.3)
dev.off()


##########################################################################
# Table 1: Cox, log-normal PH, log-normal PH + frailty
##########################################################################
labels <- c(age = "Age (years)", sex = "Male", black = "Black", other = "Other race",
            wasmarried = "Formerly married", nevermarried = "Never married",
            publicinsur = "Public insurance", privateinsur = "Private insurance",
            lowinc = "Low income", midinc = "Middle income", highinc = "High income",
            bmi = "\\textbf{Obese (BMI $>$ 30)}")
fmt <- function(v, dg = 3) ifelse(is.na(v), "", formatC(v, format = "f", digits = dg))
minus <- "{\\char\"2212}"   # a true minus sign (U+2212), not a hyphen
fmt.m <- function(v, dg = 3) gsub("-", minus, fmt(v, dg), fixed = T)

cox <- cx$cox.fit2
ln <- pm$lnormPH.fit2; ln.r <- report(ln)
lf <- fr$fr[["All covariates"]]$delayed; lf.r <- report(lf)
test <- fr$frailty.tests[fr$frailty.tests$model == "All covariates" & fr$frailty.tests$truncation == "delayed entry", ]

rows <- c(sprintf("$\\sigma$ & & & & & %s & --$^{a}$ \\\\", fmt(lf.r["sigma", "est"])),
          sprintf("$a$ & & & %s & %s & %s & %s \\\\", fmt.m(ln.r["a", "est"]), fmt(ln.r["a", "se"]), fmt.m(lf.r["a", "est"]), fmt(lf.r["a", "se"])),
          sprintf("$b$ & & & %s & %s & %s & %s \\\\", fmt(ln.r["b", "est"]), fmt(ln.r["b", "se"]), fmt(lf.r["b", "est"]), fmt(lf.r["b", "se"])),
          "\\midrule")
for (v in covars) {
  dgt <- ifelse(v == "age", 4, 3)
  rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", labels[v],
                          fmt.m(cox$estimate[v], dgt), fmt(cox$se[v], dgt),
                          fmt.m(ln.r[paste0("B", v), "est"], dgt), fmt(ln.r[paste0("B", v), "se"], dgt),
                          fmt.m(lf.r[paste0("B", v), "est"], dgt), fmt(lf.r[paste0("B", v), "se"], dgt)))
}
hr <- function(b, s) sprintf("%.2f (%.2f--%.2f)", exp(b), exp(b - 1.96*s), exp(b + 1.96*s))
fmt.ll <- function(v) gsub("-", minus, formatC(v, format = "f", digits = 2, big.mark = ","), fixed = T)
rows <- c(rows, "\\midrule",
          sprintf("HR obese, 95\\%% CI & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\",
                  hr(cox$estimate["bmi"], cox$se["bmi"]), hr(ln$estimate["Bbmi"], ln$se["Bbmi"]), hr(lf$estimate["Bbmi"], lf$se["Bbmi"])),
          sprintf("Log-likelihood & \\multicolumn{2}{c}{%s$^{b}$} & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\",
                  fmt.ll(cox$loglik), fmt.ll(ln$loglik), fmt.ll(lf$loglik)),
          sprintf("Parameters & \\multicolumn{2}{c}{%d} & \\multicolumn{2}{c}{%d} & \\multicolumn{2}{c}{%d} \\\\", cox$npar, ln$npar, lf$npar))
foot <- sprintf(paste0("$\\sigma$: frailty SD; $a$, $b$: mean and SD of log time in the log-normal baseline. $N$ = %s jobs, %s exits. Reference groups: white, married, uninsured, poor or near poor. ",
                       "$^{a}$At its boundary; LR test of $\\sigma = 0$: $\\chi^2 = %.2f$, $p = %.2f$ ",
                       "(\\textonehalf$\\chi^2_0$ + \\textonehalf$\\chi^2_1$). ",
                       "$^{b}$Log partial likelihood, not comparable with the full log-likelihoods."),
                formatC(nrow(d), big.mark = ","), formatC(sum(d$event), big.mark = ","), test$LR, test$p)
head <- c("\\begin{tabular}{@{}l rr rr rr@{}}", "\\toprule",
          " & \\multicolumn{2}{c}{Cox} & \\multicolumn{2}{c}{Log-normal PH} & \\multicolumn{2}{c}{+ frailty} \\\\",
          "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(l){6-7}",
          " & $\\beta$ & SE & $\\beta$ & SE & $\\beta$ & SE \\\\", "\\midrule")
writeLines(c("% Table 1, written by R/06_poster.R", head, rows, "\\bottomrule", "\\end{tabular}"), "figures/tab_estimates.tex")
writeLines(foot, "figures/tab_estimates_foot.tex")

# numbers quoted in the poster text
qt <- list(n = nrow(d), exits = sum(d$event), obese = sum(d$bmi), new.jobs = sum(d$talpha == 0),
           hr.ln = hr(ln$estimate["Bbmi"], ln$se["Bbmi"]),
           ph.bmi = sch$test["bmi", ], ph.global = sch$test["GLOBAL", ], ph.wasmarried = sch$test["wasmarried", ],
           bmi.only = fr$frailty.tests[fr$frailty.tests$model == "BMI only" & fr$frailty.tests$truncation == "delayed entry", ],
           aic = pm$aic0)
print(qt)
saveRDS(qt, "output/06_poster_numbers.rds")
