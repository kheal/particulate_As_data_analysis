library(here)
library(tidyverse)
library(xtable)
library(lubridate)

db_info <- read_csv("MetaData/DB_construction.csv")
db <- read_csv("MetaData/AsLipidDatabase/full_MS1db.csv")
db_summary <- db %>% group_by(LipidClass) %>%
  count() %>%
  rename(`unique lipids` = n,
         Lipid_Class = LipidClass)
db_summary_2 <- db %>% select(LipidClass, mz) %>% unique() %>%
  group_by(LipidClass) %>%
  count() %>%
  rename(`Unique empirical formulae` = n,
         Lipid_Class = LipidClass)


db_info_short <- db_info %>%
  left_join(db_summary, by = "Lipid_Class") %>%
  left_join(db_summary_2, by = "Lipid_Class") %>%
  mutate(`Arsenolipid class (short ID)` = paste(Long_name, " (", Lipid_Class, ")", sep = "")) %>%
  mutate(`Tail lengths` = paste(tail_min, "--", tail_max, sep = "")) %>%
  mutate(Saturations = paste(sat_min, "--", sat_max, sep = "")) %>%
  mutate(Saturations = ifelse(Saturations == "NA--NA", "", Saturations)) %>%
  mutate(`Tail lengths` = ifelse(`Tail lengths` == "NA--NA", "", `Tail lengths`)) %>%
  rename(References = manuscripts) %>%
  select(`Arsenolipid class (short ID)`, `Tail lengths`, Saturations, `Unique empirical formulae`,
         References)

table.latex <- xtable(db_info_short, sanitize.text.function = identity)
caption(table.latex) <- "\\label{DBSummary}Summary of lipids included in search database and associated literature. For the more rarely seen AsIsop, AsPE, and TMAsFOH lipids we used only the exact formulae seen before in literature values rather than a complete set of iterations possible lipids."
align(table.latex) <- c("p{0.1cm}", "p{4cm}","p{1.2cm}", "p{1.2cm}", "p{1.2cm}", "p{2cm}")

print(table.latex, type="latex", 
      file="Tables/ManuscriptReady/DBSummary.tex",
      size="\\fontsize{7pt}{8pt}\\selectfont",
      sanitize.text.function = identity, include.rownames=FALSE)
