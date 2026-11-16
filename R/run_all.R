##########################################################################
## run_all.R
## Runs the whole analysis in order and rebuilds the poster's figures and table.
## From the repository root:  Rscript R/run_all.R   (or source("R/run_all.R") in R)
## Then compile the poster:   xelatex prj_poster.tex  (twice, for the header)
##
## Needs R packages: statmod, tikzDevice, png (and survival, splines, which come with R).
## Run time: about 10 minutes, most of it the frailty fits in 05 (set run.baselines <- F
## there to skip the 4-minute comparison of other baselines).
##########################################################################
dir.create("output", showWarnings = F)
dir.create("figures", showWarnings = F)
for (f in c("01_km.R", "02_parametric.R", "03_cox.R", "04_diagnostics.R", "05_frailty.R", "06_poster.R")) {
  cat("\n==========", f, "==========\n")
  source(file.path("R", f), echo = F)
}
