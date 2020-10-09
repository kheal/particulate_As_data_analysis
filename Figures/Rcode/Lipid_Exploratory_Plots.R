library(tidyverse)
library(here)
library(cowplot)
library(stats)

#tidy up ALOHA Crude
icapdat1 <- read_delim("2020_01_16_Heal_AsLipidsRP_3.csv", 
                       ";", escape_double = FALSE, trim_ws = TRUE, 
                       skip = 1)
ALOHA.crude.all <- icapdat1 %>% select(X2, X3, X4, ALOHA_crude) %>%
  rename(scan = X2, analyte  = X3, type = X4)
  
ALOHA.crude.as <- ALOHA.crude.all %>%
  filter(analyte == "75As" & type == "Y") %>%
  rename(As.intensity = ALOHA_crude) %>% select(-analyte, -type)

ALOHA.crude.time <- ALOHA.crude.all %>%
  filter(analyte == "75As" & type == "Time") %>%
  rename(As.time = ALOHA_crude) %>% select(-analyte, -type)

ALOHA.crude.time.as <- left_join(ALOHA.crude.as, ALOHA.crude.time )  %>%
  filter(As.time > 200 & As.time < 2000) %>% select(-scan)

ALOHA.to.plot <- tibble("Time" = runmed(ALOHA.crude.time.as$As.time, k = 7),
                       "Intensity" = runmed(ALOHA.crude.time.as$As.intensity, k = 7))


#tidy up BATS Crude
icapdat2 <- read_delim("2020_01_16_Heal_AsLipidsRP_4.csv", 
                       ";", escape_double = FALSE, trim_ws = TRUE, 
                       skip = 1)
BATS.crude.all <- icapdat2 %>% select(X2, X3, X4, BATS07_crude) %>%
  rename(scan = X2, analyte  = X3, type = X4)

BATS.crude.as <- BATS.crude.all %>%
  filter(analyte == "75As" & type == "Y") %>%
  rename(As.intensity = BATS07_crude) %>% select(-analyte, -type)

BATS.crude.time <- BATS.crude.all %>%
  filter(analyte == "75As" & type == "Time") %>%
  rename(As.time = BATS07_crude) %>% select(-analyte, -type)

BATS.crude.time.as <- left_join(BATS.crude.as, BATS.crude.time )  %>%
  filter(As.time > 200 & As.time < 2000) %>% select(-scan)

BATS.to.plot <- tibble("Time" = runmed(BATS.crude.time.as$As.time, k = 7),
                  "Intensity" = runmed(BATS.crude.time.as$As.intensity, k = 7))
  
  
#plot them up!
g <- ggplot(data = ALOHA.to.plot, aes(x = Time,  y = Intensity)) +
  geom_line() + 
  geom_vline(xintercept = 1400, linetype = "dashed", lwd = 1.5, alpha = 0.3) +
  geom_vline(xintercept = 875, linetype = "dashed", lwd = 1.5, alpha = 0.3) +
  ggtitle("Station ALOHA")
g2 <- ggplot(data = BATS.to.plot, aes(x = Time,  y = Intensity)) +
  geom_line() + 
  geom_vline(xintercept = 1400, linetype = "dashed", lwd = 1.5, alpha = 0.3) +
  geom_vline(xintercept = 875, linetype = "dashed", lwd = 1.5, alpha = 0.3) +
  ggtitle("Hydrostation S")

g3 <- plot_grid(g, g2, align = 'lr', ncol = 1)
