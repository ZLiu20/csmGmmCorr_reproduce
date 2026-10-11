here::i_am("Fig4/make_tab3.R")
source(here::here("SupportingCode", "load_all.R"))
d <- rd_collect(here::here("Fig4", "output"))$summary
write_metric_table(d, "R_raw", here::here("Fig4", "Tab3.tex"),
                   "Discoveries in the real-data applications.", "tab:realdata",
                   c("analysis", "variant", "t"))
