#Load your libraries
library(tidyverse)
library(Rdisop)

#TO DO: finish incorporating the constructor db csv

#Load up the constructor csv
constructor <- read.csv("MetaData/DB_construction.csv")

# AsFAs-----
# define tail length and saturation possibilities 
lipid_class <- "AsFA"
constructor_sub <- constructor %>% filter(Lipid_Class == lipid_class)

tail_length1 <- c(constructor_sub$tail_min:constructor_sub$tail_max)
saturation_num1 <- c(constructor_sub$sat_min:constructor_sub$sat_max)

FA_db <- expand.grid(tail_length1, saturation_num1) 
colnames(FA_db) <- c("tail_length1", "saturation_num1")

FA_db <- FA_db %>%
  as_tibble() %>%
  mutate(LipidClass = lipid_class) %>%
  mutate(EmpiricalFormula = paste0("AsO3", "C", 
                                   tail_length1+3, "H", 
                                   2*tail_length1-2*saturation_num1+7)) %>%
  mutate(EmpiricalFormula_MH = paste0("AsO3", "C", 
                                      tail_length1+3, "H", 
                                    2*tail_length1-2*saturation_num1+8)) %>%
  mutate(mz = sapply(EmpiricalFormula_MH, 
                  function(x)getMolecule(x)$exactmass)) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_MH, 
                   function(x)getMolecule(x)$isotopes[[1]][1,2]))%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_MH, 
                         function(x)getMolecule(x)$isotopes[[1]][2,1]/getMolecule(x)$isotopes[[1]][2,2]))%>%
  mutate(Lipid_Name = paste0(LipidClass, floor(mz)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C,  tail_length1,  saturation_num1, ratio_12Cto13C)

# AsHCs----
# define tail length and saturation possibilities 
lipid_class <- "AsHC"
constructor_sub <- constructor %>% filter(Lipid_Class == lipid_class)
tail_length1 <- c(constructor_sub$tail_min:constructor_sub$tail_max)
saturation_num1 <- c(constructor_sub$sat_min:constructor_sub$sat_max)

HC_db <- expand.grid(tail_length1, saturation_num1) 
colnames(HC_db) <- c("tail_length1", "saturation_num1")

HC_db <- HC_db %>%
  as_tibble() %>%
  mutate(LipidClass = lipid_class) %>%
  mutate(EmpiricalFormula = paste0("AsO", "C", 
                                   tail_length1+2, "H", 
                                   2*tail_length1-2*saturation_num1+7)) %>%
  mutate(EmpiricalFormula_MH = paste0("AsO", "C", 
                                      tail_length1+2, "H", 
                                   2*tail_length1-2*saturation_num1+8)) %>%
  mutate(mz = sapply(EmpiricalFormula_MH, 
                     function(x)getMolecule(x)$exactmass)) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_MH, 
                         function(x)getMolecule(x)$isotopes[[1]][1,2]))%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_MH, 
                                 function(x)getMolecule(x)$isotopes[[1]][2,1]/getMolecule(x)$isotopes[[1]][2,2]))%>%
  mutate(Lipid_Name = paste0(LipidClass, floor(mz)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C,  tail_length1,  saturation_num1, ratio_12Cto13C)

# AsSugPLs  ----
# define tail length and saturation possibilities 
lipid_class <- "AsSugPL"
constructor_sub <- constructor %>% filter(Lipid_Class == lipid_class)

tail_length1 <- seq(constructor_sub$tail_min, constructor_sub$tail_max, 2)
tail_length2 <- seq(constructor_sub$tail_min, constructor_sub$tail_max, 2)
saturation_num1 <- c(constructor_sub$sat_min : constructor_sub$sat_max)
saturation_num2 <- c(constructor_sub$sat_min : constructor_sub$sat_max)

SugPL_db <- expand.grid(tail_length1, saturation_num1, 
                      tail_length2, saturation_num2) 
colnames(SugPL_db) <- c("tail_length1", "saturation_num1", 
                     "tail_length2", "saturation_num2")

SugPL_db <- SugPL_db %>%
  as_tibble() %>%
  mutate(LipidClass = lipid_class) %>%
  mutate(EmpiricalFormula = paste0("AsO14P", "C", 
                                   tail_length1+tail_length2+15, 
                                   "H", 
                                   2*(tail_length1+tail_length2)
                                   -2*(saturation_num1+saturation_num2)
                                   +28))%>%
  mutate(EmpiricalFormula_MH = paste0("AsO14P", "C", 
                                   tail_length1+tail_length2+15, 
                                   "H", 
                                   2*(tail_length1+tail_length2)
                                   -2*(saturation_num1+saturation_num2)
                                   +29)) %>%
  mutate(mz = sapply(EmpiricalFormula_MH, 
                     function(x)getMolecule(x)$exactmass)) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_MH, 
                         function(x)getMolecule(x)$isotopes[[1]][1,2]))%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_MH, 
                                 function(x)getMolecule(x)$isotopes[[1]][2,1]/getMolecule(x)$isotopes[[1]][2,2]))%>%
  mutate(Lipid_Name = paste0(LipidClass, floor(mz)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C, tail_length1, tail_length2, saturation_num1, saturation_num2) %>% 
  unique()



# lysoAsSugPLs ----
# define tail length and saturation possibilities 
lipid_class <- "lysoAsSugPL"
constructor_sub <- constructor %>% filter(Lipid_Class == lipid_class)

tail_length1 <- c(constructor_sub$tail_min:constructor_sub$tail_max)
saturation_num1 <- c(constructor_sub$sat_min:constructor_sub$sat_max)


lysoSugPL_db <- expand.grid(tail_length1, saturation_num1) 
colnames(lysoSugPL_db) <- c("tail_length1", "saturation_num1")

lysoSugPL_db <- lysoSugPL_db %>%
  as_tibble() %>%
  mutate(LipidClass = lipid_class) %>%
  mutate(EmpiricalFormula = paste0("AsO13P", "C", 
                                   tail_length1+14, 
                                   "H", 
                                   2*(tail_length1)
                                   -2*(saturation_num1)
                                   +28)) %>%
  mutate(EmpiricalFormula_MH = paste0("AsO13P", "C", 
                                   tail_length1+14, 
                                   "H", 
                                   2*(tail_length1)
                                   -2*(saturation_num1)
                                   +29)) %>%
  mutate(mz = sapply(EmpiricalFormula_MH, 
                     function(x)getMolecule(x)$exactmass)) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_MH, 
                         function(x)getMolecule(x)$isotopes[[1]][1,2]))%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_MH, 
                                 function(x)getMolecule(x)$isotopes[[1]][2,1]/getMolecule(x)$isotopes[[1]][2,2]))%>%
  mutate(Lipid_Name = paste0(LipidClass, floor(mz)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C,  tail_length1,  saturation_num1, ratio_12Cto13C)



# AsPCs  ----
# define tail length and saturation possibilities 
lipid_class <- "AsPC"
constructor_sub <- constructor %>% filter(Lipid_Class == lipid_class)

tail_length1 <- seq(constructor_sub$tail_min, constructor_sub$tail_max, 2)
tail_length2 <- seq(constructor_sub$tail_min, constructor_sub$tail_max, 2)
saturation_num1 <- c(constructor_sub$sat_min : constructor_sub$sat_max)
saturation_num2 <- c(constructor_sub$sat_min : constructor_sub$sat_max)
# construct the df for singly charged
AsPC_db <- expand.grid(tail_length1, saturation_num1, 
                        tail_length2, saturation_num2) 
colnames(AsPC_db) <- c("tail_length1", "saturation_num1", 
                        "tail_length2", "saturation_num2")


AsPC_db_singlycharged <- AsPC_db %>%
  as_tibble() %>%
  mutate(LipidClass = lipid_class) %>%
  mutate(EmpiricalFormula = paste0("AsO9PN", "C", 
                                   tail_length1+tail_length2+12, 
                                   "H", 
                                   2*(tail_length1+tail_length2)
                                   -2*(saturation_num1+saturation_num2)
                                   +25)) %>%
  mutate(EmpiricalFormula_MH = paste0("AsO9PN", "C", 
                                   tail_length1+tail_length2+12, 
                                   "H", 
                                   2*(tail_length1+tail_length2)
                                   -2*(saturation_num1+saturation_num2)
                                   +26)) %>%
  mutate(mz = sapply(EmpiricalFormula_MH, 
                     function(x)getMolecule(x)$exactmass)) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_MH, 
                         function(x)getMolecule(x)$isotopes[[1]][1,2]))%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_MH, 
                                 function(x)getMolecule(x)$isotopes[[1]][2,1]/getMolecule(x)$isotopes[[1]][2,2]))%>%
  mutate(Lipid_Name = paste0(LipidClass, floor(mz)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C, tail_length1, tail_length2, saturation_num1, saturation_num2) %>% 
  unique()

# construct the df for doubly charged
AsPC_db_doublycharged <- AsPC_db %>%
  as_tibble() %>%
  mutate(LipidClass = lipid_class) %>%
  mutate(EmpiricalFormula = paste0("AsO9PN", "C", 
                                   tail_length1+tail_length2+12, 
                                   "H", 
                                   2*(tail_length1+tail_length2)
                                   -2*(saturation_num1+saturation_num2)
                                   +25)) %>%
  mutate(EmpiricalFormula_M2H = paste0("AsO9PN", "C", 
                                      tail_length1+tail_length2+12, 
                                      "H", 
                                      2*(tail_length1+tail_length2)
                                      -2*(saturation_num1+saturation_num2)
                                      +27)) %>%
  mutate(mz = sapply(EmpiricalFormula_M2H, 
                     function(x)getMolecule(x)$exactmass)/2) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_M2H, 
                         function(x)getMolecule(x)$isotopes[[1]][1,2])/2)%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_M2H, 
                                 function(x)getMolecule(x)$isotopes[[1]][2,1]/
                                   getMolecule(x)$isotopes[[1]][2,2]))%>%
  mutate(Lipid_Name = paste0(LipidClass, floor(mz)*2-1, "_MH2")) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C, tail_length1, tail_length2, saturation_num1, saturation_num2) %>% 
  unique()


#Additional entries -----
extra_db <- read_csv("MetaData/AsLipidDatabase/Extra_entries.csv")

extra_db <- extra_db %>%
  mutate(mz = sapply(EmpiricalFormula_MH, 
                     function(x)getMolecule(x)$exactmass)) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_MH, 
                         function(x)getMolecule(x)$isotopes[[1]][1,2]))%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_MH, 
                                 function(x)getMolecule(x)$isotopes[[1]][2,1]/getMolecule(x)$isotopes[[1]][2,2])) %>%
  mutate(Lipid_Name = paste0(LipidClass, round(mz-1.5, digits = 0))) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C) %>% 
  unique()


# Combine and save out MS1 database to use as inclusion list ------
combine_db <- FA_db %>%
  bind_rows(HC_db) %>% bind_rows(SugPL_db) %>%
  bind_rows(lysoSugPL_db) %>% bind_rows(AsPC_db_singlycharged) %>%
  bind_rows(AsPC_db_doublycharged) %>%
  bind_rows(extra_db) %>%
  select(LipidClass, Lipid_Name, tail_length1, tail_length2, saturation_num1, saturation_num2, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C)

combine_db_lowermasses <- combine_db %>%
  filter(mz < 600) #range from 252.98 tp 547.33; 394 masses in total

combine_db_highermasses <- combine_db %>%
  filter(mz > 600) #range from 657.1657 tp 1071.6458; 236 masses in total

# Write out inclusion lists ------
write_csv(combine_db, "MetaData/AsLipidDatabase/full_MS1db.csv")
write_csv(combine_db_lowermasses, "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv")
write_csv(combine_db_highermasses, "MetaData/AsLipidDatabase/highmass_inclusionlist.csv")
