here::i_am("SuppTab5/make_supptab5.R")
source(here::here("SupportingCode", "load_all.R"))
d <- design_means(read_sim("Fig2", "main"))
write_metric_table(d, "incon_score", here::here("SuppTab5", "SuppTab5.tex"),
                   "Average number of score-order incongruities on the common evaluation set.",
                   "tab:incongruity")
