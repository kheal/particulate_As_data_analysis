library(tidyverse)

db <- read_csv("MetaData/AsLipidDatabase/full_MS1db.csv")
#db_ms2 <- read_csv("MetaData/AsLipidDatabase/MS2_Fragments.csv")

db_clean <- db %>%
  select(LipidClass:mz) %>%
  rename(`Lipid Class (short)` = LipidClass,
         `Lipid ID (ms1 only)` = Lipid_Name,
         `SN1 tail length` = tail_length1,
         `SN2 tail length` = tail_length2,
         `SN1 saturation` = saturation_num1,
         `SN2 saturation` = saturation_num2,
         `Empirical Formula (in ionic state)` = EmpiricalFormula,
         `m/z` = mz)

write_csv(db_clean, "Tables/ManuscriptReady/SuppTables/DB_summary_full_MS1.csv")
