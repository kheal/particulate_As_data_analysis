library(tidyverse)

# Load your files
db_info <- read_csv("MetaData/DB_construction.csv")
db <- read_csv("MetaData/AsLipidDatabase/full_MS1db.csv")

# Make it pretty
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
  rename(References = refs_numeric) %>%
  select(`Arsenolipid class (short ID)`, `Tail lengths`, Saturations, `Unique empirical formulae`,
         References)

# Write it out
write_csv(db_info_short, "Tables/ManuscriptReady/SuppTables/DB_summary.csv")
