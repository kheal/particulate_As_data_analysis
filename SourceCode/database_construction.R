#Load your libraries
library(tidyverse)
library(Rdisop)

# AsFAs-----
# define tail length and saturation possibilities 
tail_length1 <- c(7:21)
saturation_num1 <- c(0:8)

FA_db <- expand.grid(tail_length1, saturation_num1) 
colnames(FA_db) <- c("tail_length1", "saturation_num1")

FA_db <- FA_db %>%
  as_tibble() %>%
  mutate(LipidClass = "AsFA") %>%
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
  mutate(Lipid_Name = paste0(LipidClass, round(mz, digits = 0)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C)

# AsHCs----
# define tail length and saturation possibilities 
tail_length1 <- c(13:30)
saturation_num1 <- c(0:6)

HC_db <- expand.grid(tail_length1, saturation_num1) 
colnames(HC_db) <- c("tail_length", "saturation_num1")

HC_db <- HC_db %>%
  as_tibble() %>%
  mutate(LipidClass = "AsHC") %>%
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
  mutate(Lipid_Name = paste0(LipidClass, round(mz, digits = 0)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C)


# AsPls  ----
# define tail length and saturation possibilities 
tail_length1 <- c(11:19)
tail_length2 <- c(11:19)
saturation_num1 <- c(0:4)
saturation_num2 <- c(0:4)

SugPL_db <- expand.grid(tail_length1, saturation_num1, 
                      tail_length2, saturation_num2) 
colnames(SugPL_db) <- c("tail_length1", "saturation_num1", 
                     "tail_length2", "saturation_num2")

SugPL_db <- SugPL_db %>%
  as_tibble() %>%
  mutate(LipidClass = "AsSugPL") %>%
  mutate(EmpiricalFormula = paste0("AsO14P", "C", 
                                   tail_length1+tail_length2+15, 
                                   "H", 
                                   2*(tail_length1+tail_length2)
                                   -2*(saturation_num1+saturation_num2)
                                   +28)) %>%
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
  mutate(Lipid_Name = paste0(LipidClass, round(mz, digits = 0)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C) %>% 
  unique()



# lysoAsPls ----
# define tail length and saturation possibilities 
tail_length1 <- c(11:19)
saturation_num1 <- c(0:4)

lysoSugPL_db <- expand.grid(tail_length1, saturation_num1) 
colnames(lysoSugPL_db) <- c("tail_length1", "saturation_num1")

lysoSugPL_db <- lysoSugPL_db %>%
  as_tibble() %>%
  mutate(LipidClass = "lysoAsSugPL") %>%
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
  mutate(Lipid_Name = paste0(LipidClass, round(mz, digits = 0)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C) %>% 
  unique()



# AsPCs  ----
# define tail length and saturation possibilities 
tail_length1 <- c(14:21)
tail_length2 <- c(14:21)
saturation_num1 <- c(0:6)
saturation_num2 <- c(0:6)

# construct the df
AsPC_db <- expand.grid(tail_length1, saturation_num1, 
                        tail_length2, saturation_num2) 
colnames(AsPC_db) <- c("tail_length1", "saturation_num1", 
                        "tail_length2", "saturation_num2")

AsPC_db <- AsPC_db %>%
  as_tibble() %>%
  mutate(LipidClass = "AsPC") %>%
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
  mutate(Lipid_Name = paste0(LipidClass, round(mz, digits = 0)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C) %>% 
  unique()



#Additional entries -----
extra_db <- read_csv("MetaData/AsLipidDatabase/Extra_entries.csv")

extra_db <- extra_db %>%
  mutate(mz = sapply(EmpiricalFormula_MH, 
                     function(x)getMolecule(x)$exactmass)) %>%
  mutate(mz_13C = sapply(EmpiricalFormula_MH, 
                         function(x)getMolecule(x)$isotopes[[1]][1,2]))%>%
  mutate(ratio_12Cto13C = sapply(EmpiricalFormula_MH, 
                                 function(x)getMolecule(x)$isotopes[[1]][2,1]/getMolecule(x)$isotopes[[1]][2,2]))%>%
  mutate(Lipid_Name = paste0(LipidClass, round(mz, digits = 0)-1)) %>%
  select(LipidClass, Lipid_Name, EmpiricalFormula, mz, mz_13C, ratio_12Cto13C) %>% 
  unique()


# Combine and save out MS1 database to use as inclusion list ------
combine_db <- FA_db %>%
  bind_rows(HC_db) %>% bind_rows(SugPL_db) %>%
  bind_rows(lysoSugPL_db) %>% bind_rows(AsPC_db) %>%
  bind_rows(extra_db) 

combine_db_lowermasses <- combine_db %>%
  filter(mz < 600) #range from 252.98 tp 547.33; 394 masses in total

combine_db_highermasses <- combine_db %>%
  filter(mz > 600) #range from 657.1657 tp 1071.6458; 236 masses in total

# Write out inclusion lists ------
write_csv(combine_db_lowermasses, "MetaData/AsLipidDatabase/lowmass_inclusionlist.csv")
write_csv(combine_db_highermasses, "MetaData/AsLipidDatabase/highmass_inclusionlist.csv")
