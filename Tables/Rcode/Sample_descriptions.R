library(here)
library(tidyverse)
library(xtable)
library(lubridate)

samp.info <- read_csv("MetaData/SampleMetaDat.csv")

samp.info.short <- samp.info %>%
  select(`Sample ID`:L_digested) %>%
  arrange(Date_Collected) %>%
  rename(Date = Date_Collected,
         Sample = `Sample ID`,
         `Filter set up (\\textmu m)` = Filter_Set_Up,
         `Volume (L)` = L_Filtered,
         `Volume extracted for speciation (L)` = L_extracted,
         `Volume extracted for As digestion (L)` = L_digested) %>%
  mutate(Sample = Sample %>% str_replace_all("_PS", "-PS"))
         

table.latex <- xtable(samp.info.short, sanitize.text.function = identity)
caption(table.latex) <- "\\label{SampleDescriptions}Summary of samples collected and analyzed for total particulate arsenic and arsenic speciation of particles."
align(table.latex) <- c("p{0.1cm}","p{1.2cm}", "p{1cm}", "p{1cm}", "p{1cm}","p{1cm}","p{1cm}","p{2cm}","p{2cm}")

print(table.latex, type="latex", 
      file="Tables/ManuscriptReady/SampleDescriptions.tex",
      size="\\fontsize{7pt}{8pt}\\selectfont",
      sanitize.text.function = identity, include.rownames=FALSE)
