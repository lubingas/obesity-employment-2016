# Obesity and employment duration (2016)

Do obese workers leave employment sooner? This is my UW CSSS 544 (event history analysis) project from 2016, with one update in 2026: the frailty model, which I originally fitted in mle, is now fitted in R with Gauss–Hermite quadrature.

**Poster:** [`prj_poster.pdf`](prj_poster.pdf) (corrected November 2026). The poster as presented in March 2016 is [`prj_poster_2016.pdf`](prj_poster_2016.pdf).

## Data and methods

- **Data:** 6,075 workers aged 17–55 in Medical Expenditure Panel Survey (MEPS) Panel 14, 2009–2010. 1,131 left employment during follow-up. 1,682 (28%) were obese (BMI ≥ 30).
- **Clock:** time is years since the job started, with left truncation at the tenure reached when follow-up began in 2009.
- **Models:**
  - Kaplan–Meier estimates and a Cox model, both allowing for left truncation
  - parametric proportional-hazards models (exponential, Weibull, log-normal, Gompertz)
  - a log-normal proportional-hazards model with normal unobserved heterogeneity (frailty)
- **Covariates:** age, sex, race, marital status, insurance and income.

## Results

| | Hazard ratio, obese vs not (95% CI) |
|---|---|
| Cox model | 1.20 (1.06–1.37) |
| Log-normal PH | 1.20 (1.06–1.37) |
| Log-normal PH with frailty | 1.20 (1.06–1.37); σ at its boundary (0) |

Obese workers had about a 20% higher hazard of leaving employment. Allowing for unobserved heterogeneity changes nothing: the frailty variance is estimated at zero.

## The 2026 update: Gauss–Hermite quadrature in place of mle

In 2016 I fitted the frailty model with mle, a programming language for building likelihood models written by the course instructor, Darryl Holman (Holman 2003). It is free for non-commercial use and is distributed as a stand-alone executable, so I rewrote the model in R, where anyone can inspect and rerun it with widely used open tools.

[`Data/prj_frailty_ghq.R`](Data/prj_frailty_ghq.R) writes the frailty likelihood in R. Each person's likelihood is integrated over the frailty with 30-node Gauss–Hermite quadrature (checked against 60 nodes, and person by person against R's `integrate()`), on the log scale, and maximised with `maxLik`. It:

- **Reproduces the 2016 mle fit:** log-likelihood −3381.601, the same coefficients, and σ at its lower boundary.
- **Fixes the boundary standard errors:** mle's standard errors at the boundary were unreliable (an SE of 5.19 for the log-normal scale b). The R fit gives 0.064.
- **Tests σ = 0 at the boundary:** a likelihood-ratio test with the boundary correction (half of a χ²₁), p = 0.5. A profile likelihood for σ peaks at 0.
- **Compares two ways of handling left truncation:**
  - dividing by S(t_α | z) inside the integral, as the poster's likelihood does
  - standard delayed entry, which divides by S(t_α) averaged over the frailty

  With all covariates, both give σ ≈ 0.
- **Compares other baselines** (exponential, Weibull, Gompertz, log-logistic, gamma, piecewise exponential), each with no heterogeneity, normal frailty or gamma frailty. With all covariates, the log-normal model without frailty has the lowest AIC. The "frailty" that appears under the more rigid baselines, or with BMI alone, stands in for a misspecified baseline or omitted covariates. It is not evidence of unobserved heterogeneity.

Results are in [`Data/frailty_baseline_comparison.csv`](Data/frailty_baseline_comparison.csv) and [`Data/frailty_profile.pdf`](Data/frailty_profile.pdf).

## Poster corrections (November 2026)

- **Table 1, frailty column:** replaced with the R quadrature estimates and labelled "Log Normal + frailty".
- **Methods:** now say how the frailty integral is computed.
- **Outcome label:** "unemployment" is now "leaving employment". The outcome is any exit from work, including leaving the labour force, not only unemployment.
- **Frailty likelihood:** the hazard term read e^{xβ} where it should be e^{xβ+z}.
- **Confidence interval:** the lower bound rounds to 6%, not 5% (exp(0.184 − 1.96 × 0.066) = 1.056).
- **Dated note:** a note under the methods explains the move from mle to R and lists these corrections. The body text and equations are set slightly smaller to make room for it.

Apart from these, the Cox and log-normal columns, the figures and the conclusions are as in 2016.

## What comes next

Looking back at this analysis, four things limit it:

- **BMI timing.** BMI was measured partway through follow-up, so some of it comes after the exits it is meant to explain. Reverse causation is possible.
- **Exit dates.** Exit dates were reconstructed from MEPS round data and are not reliable.
- **One outcome.** All exits are treated alike, so being laid off, being fired, quitting for ill health and quitting for another job cannot be told apart.
- **One panel.** A single two-year MEPS panel limits the number of exits.

I am extending this work to fix these:

- **More data:**
  - nine MEPS panels (2009–2018)
  - the National Longitudinal Survey of Youth 1979 (NLSY79), which follows the same people from 1979 to 2022 and records job end months and reasons for leaving
- **BMI before follow-up:** BMI is measured before the follow-up it is used for.
- **Exits by reason:** exits are split by stated reason, with business closures as a negative control and moves to another job as a competing exit.

The results will be posted when ready.

## Files

| Path | What |
|---|---|
| `prj_poster.tex`, `prj_poster.pdf` | Corrected poster; the figures are the tikz files `kmlog.tex`, `phcheck.tex`, `sho.tex`, `devr.tex`, `mres.tex` |
| `prj_poster_2016.pdf` | The poster as presented in 2016 |
| `fonts/` | The poster's fonts (Encode Sans, Open Sans Condensed; SIL Open Font License), so it compiles anywhere with `xelatex prj_poster.tex` |
| `Data/prj_analysis.R` | The 2016 analysis: Kaplan–Meier, Cox, parametric PH likelihoods with left truncation, residual checks and the poster figures. Left as written in 2016, so it still points at the old Windows paths. |
| `Data/prj_frailty_ghq.R` | The 2026 R rewrite of the frailty model. Needs `maxLik`, `survival` and `statmod`. Sections 1–3 take about a minute; the profile likelihood and the baseline comparison take 6–8 minutes each. |
| `Data/mepsempcsv.csv` | The analysis file, one row per person |

## References

Holman DJ (2003) mle: A programming language for building likelihood models. Version 2.1. Volume 1, User's Manual, and Volume 2, Reference Manual. [faculty.washington.edu/djholman/mle](https://faculty.washington.edu/djholman/mle/index.html)

## Data source

Agency for Healthcare Research and Quality, Medical Expenditure Panel Survey, Panel 14 (2009–2010) public-use files, [meps.ahrq.gov](https://meps.ahrq.gov).

## Citation

Lubinga SJ. *Obesity and Employment Duration: Continuous-Time Hazards Analyses.* Poster, CSSS 544, University of Washington, 2016; corrected 2026. https://github.com/lubingas/obesity-employment-2016

Code: MIT License. Fonts: SIL Open Font License (see `fonts/`).
