here::i_am("Tab2/make_tab2.R")
source(here::here("SupportingCode", "load_all.R"))
fdp <- sources_means("FDP")
power <- sources_means("power")
lines <- c(sources_tex(fdp, 1L, "tab:source-fdp", "Mean FDP across simulation sources."),
           "\\clearpage",
           sources_tex(power, 1L, "tab:source-power", "Mean power across simulation sources."))
write_tex(lines, here::here("Tab2", "Tab2.tex"))
