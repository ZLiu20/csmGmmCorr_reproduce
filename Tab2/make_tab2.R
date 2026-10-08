# Table 2: mean FDP at t = 1 under alternative estimators and data-generating processes.
# Uses the output of the folders Fig2, Fig3 and Tab2:  Rscript Tab2/make_tab2.R
here::i_am("Tab2/make_tab2.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))

T <- sources_means("FDP")
fwrite(T, here::here("Tab2", "Tab2_sources_fdp.csv"))
cap <- function(t0, powref) paste0("Mean FDP at $q=0.10$ and $t=", t0, "$ under alternative estimators and data-generating ",
  "processes, averaged over the first and fifth effect-size designs and one signal-proportion design of each target ",
  "($\\tau_1=0.02$ for pleiotropy and replication, $0.04$ for mediation); the last block uses the fifth effect-size design ",
  "with $\\rho=0.8$. Replicates: 500 (first row), 200 (estimation rows and $\\rho=0.8$; 50 for the fixed null correlation), ",
  "100 (independent test set and other data-generating rows; 200 for symmetric statistics). Power is in ", powref, ".")
write_tex(sources_tex(T, 1, "tab:sources", cap(1, "Supplementary Table~\\ref{supp-tab:sources-power}")), here::here("Tab2", "Tab2.tex"))
