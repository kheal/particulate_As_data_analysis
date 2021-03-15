library(tidyverse)
library(Rdisop)

db <- read_csv("MetaData/AsLipidDatabase/full_MS1db.csv")
db_ms2 <- read_csv("MetaData/AsLipidDatabase/MS2_Fragments.csv")

db_ms2_2 <- db_ms2 %>%
  filter(str_detect(Formula, "As")) %>%
  filter(LipidClass %in% db$LipidClass) %>%
  mutate(mz = sapply(Formula, 
                   function(x)getMolecule(x)$exactmass)) 


db_clean <- db_ms2_2 %>%
  select(LipidClass, Formula, mz) %>%
  rename(`Lipid Class (short)` = LipidClass,
         `Empirical Formula (in ionic state)` = Formula,
         `m/z` = mz)

write_csv(db_clean, "Tables/ManuscriptReady/SuppTables/DB_summary_full_MS2.csv")
