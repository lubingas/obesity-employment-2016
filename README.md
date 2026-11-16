# Obesity and employment duration (2016)

Do obese workers leave employment sooner? This is my UW CSSS 544 (event history analysis) project from 2016, partly built in R at the time. In 2026 I fixed and streamlined the code, now fitting every model with `optim`, and rebuilt in R the frailty model I originally fitted in mle.

**Poster:** [`prj_poster.pdf`](prj_poster.pdf) (revised November 2026). The poster as presented in March 2016 is [`prj_poster_2016.pdf`](prj_poster_2016.pdf).

## Data and methods

- **Data:** 6,075 workers aged 17–55 in Medical Expenditure Panel Survey (MEPS) Panel 14, 2009–2010. 1,131 left employment during follow-up.
- **Exposure:** obese if the mean BMI over rounds 3 and 5 was above 30: 1,682 people (28%). Note that this is above 30, not 30 or more. 27 people with a mean of exactly 30.0 count as not obese.
- **Clock:** years since the job started. Jobs already under way when follow-up began are left-truncated at the tenure reached then; 1,227 jobs began during follow-up.
- **Covariates:** age, sex, race, marital status, insurance and income.
- **Models:**
  - Kaplan–Meier estimates and a log-rank test, with left-truncated risk sets
  - a Cox model (delayed entry, Breslow ties)
  - six parametric baselines (exponential, Weibull, log-normal, gamma, Gompertz, log-logistic)
  - a log-normal proportional-hazards (PH) model
  - the same model with normal unobserved heterogeneity (frailty)
- **Diagnostics:**
  - log(−log S) plots
  - a Schoenfeld-residual score test of proportional hazards
  - Cox–Snell, martingale and deviance residuals
  - a likelihood ratio test and profile likelihood for the frailty SD σ

## Results

| | Hazard ratio, obese vs not (95% CI) |
|---|---|
| Cox model | 1.20 (1.06–1.37) |
| Log-normal PH | 1.20 (1.06–1.37) |
| Log-normal PH with frailty | 1.20 (1.06–1.37); σ at its boundary (0) |

Obese workers had about a 20% higher hazard of leaving employment.

- **Log-rank test:** χ² = 8.98, p = 0.003, with left-truncated risk sets. Ignoring left truncation gives χ² = 0.67, so the truncation matters.
- **Proportional hazards:** the assumption holds for BMI (p = 0.57). The global test across all 12 covariates gives p = 0.064. Only being formerly married departs (p = 0.006).
- **Baseline:** the log-normal has the lowest AIC of the six baselines, with and without covariates.
- **No unobserved heterogeneity:**
  - With all covariates, σ is estimated at 0 (LR = 0, p = 0.5). The profile likelihood is highest at σ = 0.
  - With BMI alone and delayed entry, σ = 0.83 (p < 0.001). Here the frailty stands in for the omitted covariates.
  - Across the six baselines with normal or gamma frailty, a frailty appears only with too few covariates or a rigid baseline.

## The 2026 update

The 2016 code (`Data/prj_analysis.R`, left as it was) fitted the parametric models with `maxLik` and used `survfit()`, `coxph()` and `cox.zph()` for the rest. It also contained my own Kaplan–Meier and Cox partial-likelihood functions. The Cox one never ran: its risk-set vectors were not set up inside the function, and its design matrix came from a data frame that didn't exist.

The update in `R/` fixes and streamlines that code and finishes the job:

- **Shared design matrix:** one design matrix from `makeX()`, and the linear predictor is always `xb <- linpred(beta, x)`.
- **Kaplan–Meier and log-rank test:** written with left-truncated risk sets, counted with sorted sums instead of loops.
- **Cox model:** our own partial likelihood and its score, with delayed entry and Breslow ties, plus the Breslow baseline hazard.
- **Parametric likelihoods:** every one is written on the log scale with left truncation, for any baseline that supplies log h₀(t) and H₀(t).
- **Frailty:** normal frailty is integrated with 30-point Gauss–Hermite quadrature. It is checked against 60 points and against `integrate()`.
- **Residuals:** Cox–Snell, martingale, deviance and Schoenfeld residuals, plus the score test of proportional hazards.
- **Fitting:** everything is maximised with `optim` (BFGS), with standard errors from the Hessian.
- **Checks:** each estimator is checked against the `survival` package as it runs:
  - `survfit()`
  - `coxph()`, including its score test and residuals
  - `basehaz()`
  - `cox.zph()`

  They agree to within 10⁻⁶. The checks are printed when you run the scripts.

### What mle did, and the 2016 standard errors

In 2016 I fitted the frailty model with mle, a programming language for building likelihood models written by the course instructor, Darryl Holman (Holman 2003). It is free for non-commercial use and is distributed as a stand-alone executable.

The mle program file itself has not survived; only its data file, `anafil.dat`, has. But its output can be checked:

- **mle maximised the likelihood on the poster.** At mle's 2016 estimates, that likelihood gives −3381.60094 in R; mle reported −3381.601.
- **The 2016 table had the first three standard errors in the frailty column one row out of place.** It showed σ 0.17058, a 0.06376 and b 5.19397. Recomputed at mle's estimates, the SE of a is 0.163, the SE of b is 0.064, and the large value belongs to σ, which sits at its boundary. So the 5.19 was a copying slip in 2016. It does not mean mle's standard errors were unreliable, as an earlier version of this README said.
- **mle's other standard errors were 3–7% larger** than those from the exact Hessian. This is probably due to how mle approximates the Hessian, but without the program file I can't confirm it.

### Left truncation with frailty

The 2016 likelihood divides by S(t_α | z) inside the frailty integral, as if each person's frailty were drawn when the survey began. The poster now uses standard delayed entry instead. Frailty is fixed from the start of the job, so the people still employed at entry are a selected group, and the whole integral is divided by S(t_α) averaged over the frailty. The two agree when σ = 0, and both give σ = 0 with all covariates. `05_frailty.R` fits both.

### Numbers that moved

The Cox and log-normal coefficients and log-likelihoods match 2016 to every printed digit. In the frailty column, mle had stopped just short of the maximum (−3381.60100 against −3381.60085), so its 2016 estimates differ in the third or fourth decimal (a 0.573 against 0.571).

- **Standard errors, log-normal column:** these now come from `optim`'s Hessian, checked against `numDeriv`. They differ from the 2016 `maxLik` values by less than 1%; for example, the SE of a is 0.163 rather than 0.164.
- **Standard errors, frailty column:** these are now in the right rows. They come from the R refit, not mle.
- **Log-rank test:** it is now computed rather than typed in. It is χ² = 8.98, from the hypergeometric variance. The 2016 value, χ² = 8.97 with p = 0.002738, is exactly the Cox score test with Breslow ties, which also allows for truncation.
- **Proportional-hazards test:** now the exact score test that `cox.zph()` has used since survival 3.0. The 2016-era approximation is in `R/00_functions.R` as `test.approx`.

## Poster corrections (November 2026)

The poster was rebuilt in LaTeX with the fonts and colours of my website: Source Serif 4 and Inter, both under the SIL Open Font License, shipped in `fonts/`. All figures and Table 1 are produced by `R/06_poster.R`.

**Fixed from 2016:**

- **Cox partial likelihood:** the δᵢ exponent was missing, so as written censored people contributed. The risk set is now R(tᵢ).
- **Notation:** the Kaplan–Meier notation now uses d_k and strict entry (α_j < t_k). The frailty likelihood now has subscripts and dz.
- **Frailty column SEs:** the standard errors of σ, a and b are now in the right rows.
- **Obesity cut-off:** BMI above 30, as analysed.
- **Outcome:** "leaving employment", any exit from work, not "unemployment".
- **Log-likelihoods in the table:** the Cox log partial likelihood is now marked as not comparable with the full log-likelihoods.

**New:** the profile likelihood for σ, the Cox–Snell comparison of the Cox and log-normal models, and martingale residuals against age.

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

## Theoretical model for the extension

The extension treats leaving a job as the end of a worker–employer match. It asks whether obesity changes how fast matches end, and whose decision ends them. This section sets out the model that guides it; the results will be posted when ready.

### The separation model

A match lasts while both sides prefer it to their alternatives. Shocks arrive over time, and the match ends when one pushes either side past its threshold. For worker $`i`$ in job $`j`$ at tenure $`t`$:

```math
\text{employer keeps the match while}\quad y_{ij}(t) - w_{ij} - c_{ij} - d_{ij} \ge 0
```

```math
\text{worker stays while}\quad w_{ij} + a_{ij} - e_{ij}(t) \ge \max\{\, b_i(t),\ \text{best outside offer} \,\}
```

The terms are:

- $`y`$, productivity
- $`w`$, the wage
- $`c`$, non-wage cost, including the employer's share of health-insurance premiums
- $`d`$, any distaste the employer, coworkers or customers have for the worker
- $`a`$, job amenities
- $`e`$, the disutility of work
- $`b`$, the value of not working

Wages are set at hire and adjust slowly within a match, so shocks end matches rather than just changing pay. Whatever the wage absorbs at hire shows up as a wage gap, not as extra separations.

Each kind of shock produces its own exit, so the competing risks are defined by who ends the match and why:

| Exit (who ends it) | Shock |
|---|---|
| Layoff (employer) | Productivity or demand falls, so $`y < w + c + d`$ |
| Business closure (firm) | The whole business exits, whoever the worker is |
| Health quit (worker) | Illness raises $`e`$ or lowers $`y`$ |
| Other quit to non-work (worker) | $`b`$ rises: family, school, retirement |
| Job-to-job (worker) | An outside offer beats the current match |

Where the worker goes next (another job or none) is a separate outcome. It reflects outside offers.

Suppose shocks arrive at rates that depend on tenure, and covariates scale the chance that a shock crosses a threshold. Then each exit $`k`$ has a mixed proportional hazard:

```math
\theta_k(t \mid x_i, \mathrm{ob}_i, v_{ik}) = \theta_{0k}(t)\,\exp\left(x_i'\beta_k + \gamma_k\,\mathrm{ob}_i + v_{ik}\right)
```

Here $`t`$ is tenure, $`\mathrm{ob}_i`$ is obesity measured before the interval, and $`v_{ik}`$ is unobserved heterogeneity, which can be correlated across exits. With a single exit type, a log-normal $`\theta_0`$ and normal $`v`$, this is the 2016 frailty model.

The tenure shape also comes from theory. In Jovanovic's (1979) matching model, match quality is learned on the job. Poor matches end early and survivors build firm-specific skills, so the hazard rises briefly and then falls. That is the case for a hump-shaped baseline such as the log-normal.

### How obesity enters

| Channel | What moves in the model |
|---|---|
| Health and productivity | $`y`$ falls, $`e`$ rises, health shocks arrive faster |
| Taste-based discrimination (Becker 1957) | $`d > 0`$, so obese workers sit nearer the layoff threshold |
| Statistical discrimination (Altonji & Pierret 2001) | $`y`$ is judged lower early in a job; the gap fades as the employer learns |
| Employer health-insurance cost (Bhattacharya & Bundorf 2009) | $`c`$ rises, but only for workers with employer coverage |
| Job lock (Madrian 1994) | Workers who expect higher medical spending value employer coverage more, so quits fall |
| Outside offers (Rooth 2009) | Fewer and worse offers |
| Value of not working | $`b`$ shifts, through disability benefits or a spouse's income or coverage |

The employer-side channels raise separations only to the extent that the wage did not absorb them at hire.

### What each channel predicts

A + means obesity raises that hazard, − lowers it, and 0 leaves it unchanged.

| Channel | Layoff | Closure | Health quit | Other quit or retire | Job-to-job |
|---|---|---|---|---|---|
| Health and productivity | + | 0 | ++ | + | − |
| Taste-based discrimination | + | 0 | 0 | 0 | − |
| Statistical discrimination | + early in the job | 0 | 0 | 0 | − |
| Employer insurance cost | + if covered by the employer | 0 | 0 | 0 | 0 |
| Job lock | 0 | 0 | − | − | − |
| Outside offers | 0 | 0 | 0 | 0 | − |
| Unobserved confounders (time preference, earlier health, ability) | + | + if obese workers sort into fragile firms | + | + | + |

Two exits separate the channels from confounding:

- **Closures:** they happen whoever works at the firm, so no individual-level channel predicts an obesity effect. That makes them a negative-control outcome (Gibbons & Katz 1991; Lipsitch et al. 2010).
- **Job-to-job moves:** every causal channel predicts fewer moves or no change, while confounding by impatience predicts more.

### The likelihood

In MEPS, exits are seen only between interviews. For a job held at $`t_0`$, the chance of exit $`k`$ in the interval $`(t_0, t_1]`$ is:

```math
\Pr(\text{exit } k \text{ in } (t_0,t_1] \mid \text{in job at } t_0) \approx \frac{\Delta\Lambda_k}{\sum_m \Delta\Lambda_m}\left[1-\exp\left(-\sum_m \Delta\Lambda_m\right)\right],
\qquad \Delta\Lambda_k = e^{x'\beta_k+\gamma_k\mathrm{ob}}\left[\Lambda_{0k}(t_1)-\Lambda_{0k}(t_0)\right]
```

No exit has probability $`\exp(-\sum_m \Delta\Lambda_m)`$. The formula is exact when the exit-specific hazards stay proportional to each other within an interval (Narendranathan & Stewart 1993).

Most jobs are already under way when first observed. So frailty is conditioned on the job lasting to entry tenure $`t_{e,i}`$ (Guo 1993; van den Berg & Drepper 2016):

```math
L_i = \frac{\int \prod_{r} P_{ir}(v)\; S(t_{e,i}\mid v)\, g(v)\, dv}{\int S(t_{e,i}\mid v)\, g(v)\, dv}
```

Here $`P_{ir}(v)`$ is interval $`r`$'s contribution given frailty $`v`$, and $`S`$ is survival in the job across all exits. This is the delayed-entry form the revised poster uses for the 2016 frailty model.

## How to run it

From the repository root:

```
Rscript R/run_all.R          # about 10 minutes; writes output/ and figures/
xelatex prj_poster.tex       # run twice, for the header
```

You need these R packages: `statmod`, `tikzDevice`, `png` (with `survival` and `splines` from base R). `survival` is used only for the checks; `numDeriv` is optional.

## Files

| Path | What |
|---|---|
| `R/00_functions.R` | Every estimator: data and design matrix, `KMcalc`, `logrank`, `coxPHLik`/`coxPHScore`, `coxBaseHaz`, residuals and `schoenfeld`, parametric baselines, `PHLik` (with or without Gauss–Hermite frailty), and the `optim` wrapper `fitML` |
| `R/01_km.R` | Kaplan–Meier by BMI, log-rank test, log(−log S) check |
| `R/02_parametric.R` | Six baselines by AIC, log-normal and Weibull by BMI group, and log-normal PH models |
| `R/03_cox.R` | Cox models from our partial likelihood |
| `R/04_diagnostics.R` | Cox–Snell, martingale, deviance and Schoenfeld residuals, plus the PH test |
| `R/05_frailty.R` | Frailty fits, LR tests, the 60-node and `integrate()` checks, the profile likelihood, and other baselines × normal or gamma frailty |
| `R/06_poster.R` | Poster figures (tikz) and Table 1 |
| `R/run_all.R` | Runs 01–06 in order |
| `output/` | Results tables (`.csv`); the saved fits (`.rds`) are rebuilt by `run_all.R` and not tracked |
| `figures/` | Poster figures and Table 1, written by `06_poster.R` |
| `prj_poster.tex`, `prj_poster.pdf` | The revised poster |
| `prj_poster_2016.pdf` | The poster as presented in 2016 |
| `fonts/` | Source Serif 4 and Inter (SIL Open Font License) |
| `Data/prj_analysis.R` | The 2016 analysis, left as written in 2016, so it still points at the old Windows paths |
| `Data/mepsempcsv.csv` | The analysis file, one row per person |

## References

Altonji JG, Pierret CR (2001) Employer learning and statistical discrimination. *Quarterly Journal of Economics* 116(1):313–350.

Becker GS (1957) *The Economics of Discrimination.* University of Chicago Press.

Bhattacharya J, Bundorf MK (2009) The incidence of the healthcare costs of obesity. *Journal of Health Economics* 28(3):649–658. https://doi.org/10.1016/j.jhealeco.2009.02.009

Gibbons R, Katz LF (1991) Layoffs and lemons. *Journal of Labor Economics* 9(4):351–380.

Guo G (1993) Event-history analysis for left-truncated data. *Sociological Methodology* 23:217–243.

Holman DJ (2003) mle: A programming language for building likelihood models. Version 2.1. Volume 1, User's Manual, and Volume 2, Reference Manual. [faculty.washington.edu/djholman/mle](https://faculty.washington.edu/djholman/mle/index.html)

Jovanovic B (1979) Job matching and the theory of turnover. *Journal of Political Economy* 87(5):972–990.

Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010) Negative controls: a tool for detecting confounding and bias in observational studies. *Epidemiology* 21(3):383–388.

Madrian BC (1994) Employment-based health insurance and job mobility: is there evidence of job-lock? *Quarterly Journal of Economics* 109(1):27–54.

Narendranathan W, Stewart MB (1993) Modelling the probability of leaving unemployment: competing risks models with flexible base-line hazards. *Applied Statistics* 42(1):63–83.

Rooth D-O (2009) Obesity, attractiveness, and differential treatment in hiring: a field experiment. *Journal of Human Resources* 44(3):710–735.

van den Berg GJ, Drepper B (2016) Inference for shared-frailty survival models with left-truncated data. *Econometric Reviews* 35(6):1075–1098.

## Data source

Agency for Healthcare Research and Quality, Medical Expenditure Panel Survey, Panel 14 (2009–2010) public-use files, [meps.ahrq.gov](https://meps.ahrq.gov).

## Citation

Lubinga SJ. *Obesity and Employment Duration: Continuous-Time Hazards Analyses.* Poster, CSSS 544, University of Washington, 2016; revised 2026. https://github.com/lubingas/obesity-employment-2016

Code: MIT License. Fonts: SIL Open Font License (see `fonts/`).
