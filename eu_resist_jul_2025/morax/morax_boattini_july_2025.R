library(readxl)
library(tidyverse)
library(data.table)


# Look-up species-antibacterials combinations resistance cutoffs ------------

# MORAX_Complete_dataset <- read_xlsx(path="2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)
# names(MORAX_Complete_dataset)
# 
# Look_up <- MORAX_Complete_dataset %>% select(`Species identification`, `Benzylpenicillin MIC`:`Doxycycline inhib zone diam`) 
# 
# 
# Look_up <- Look_up %>% gather(Abx, Rx, `Benzylpenicillin MIC`:`Doxycycline inhib zone diam`) %>% drop_na() %>%
#   select(-Rx) %>%   distinct()
# 
# length(unique(Look_up$`Species identification`))
# 
# length(unique(Look_up$Abx))
# 
# Look_up <- Look_up %>% mutate(Abx=str_replace(Abx, " MIC", "")) %>%
#   mutate(Abx=str_replace(Abx, " inhib zone diam", "")) %>% distinct()
# 
# fwrite(Look_up, "Look_up_MORAX_Complete_dataset.csv")

# -----------
# Resistant/susceptible flags using MIC --------------------------------

Look_up <- fread("../data/Look_up_MORAX_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MORAX_Complete_dataset$`Species identification`))

sort(unique(MORAX_Complete_dataset$`Species identification`))

names(MORAX_Complete_dataset)

MIC_data <- MORAX_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MORAX_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Benzylpenicillin MIC`:`Doxycycline MIC`) 


MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="≤1", 0.1, MIC)) %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


summary_MIC_concent <- MIC_data %>% filter(!is.na(MIC)) %>%
  group_by(`Species identification`, Abx) %>% 
  summarise(mean=mean(MIC, na.rm=T), 
            sd=sd(MIC), 
            mean=mean(MIC), 
            median=median(MIC), 
            Q90=quantile(MIC, 0.90),
            n=n()) 


fwrite(summary_MIC_concent, "../out/summary_MIC_concent_Jul_16_MIC90.csv")



EUCAST_resist_counts <- MIC_data %>% filter(!is.na(EUCAST_mic_Resist )) %>%
  group_by(`Species identification`, Abx, EUCAST_mic_Resist) %>% count() %>%
  spread(key=EUCAST_mic_Resist, value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_eucast=`1`+`0`, perc_r_eucast=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

CLSI_m100_resist_counts <- MIC_data %>% filter(!is.na(CLIST_mic_m100_Resist )) %>%
  group_by(`Species identification`, Abx, CLIST_mic_m100_Resist  ) %>% count() %>%
  spread(key=CLIST_mic_m100_Resist  , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_clsi=`1`+`0`, perc_r_clsi=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

casfm_resist_counts <- MIC_data %>% filter(!is.na(CASFM_mic_Resist  )) %>%
  group_by(`Species identification`, Abx, CASFM_mic_Resist  ) %>% count() %>%
  spread(key=CASFM_mic_Resist  , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_casfm=`1`+`0`, perc_r_casfm=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))


fwrite(MIC_data, "../out/MIC_data_July_16.csv")

fwrite(EUCAST_resist_counts, "../out/EUCAST_resist_counts_Jul_16.csv")
fwrite(CLSI_m100_resist_counts, "../out/CLSI_m100_resist_counts_Jul_16.csv")
fwrite(casfm_resist_counts, "../out/casfm_resist_counts_Jul_16.csv")


# ---------------------------
# Resistant/susceptible flags using inhibitory zone diameter --------------------------------

Look_up <- fread("../data/Look_up_MORAX_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MORAX_Complete_dataset$`Species identification`))

sort(unique(MORAX_Complete_dataset$`Species identification`))

names(MORAX_Complete_dataset)

ZONE_data <- MORAX_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MORAX_Complete_dataset %>% select(contains(" inhib zone"))
  )

names(ZONE_data)

ZONE_data <- ZONE_data %>%
  gather(Abx, ZONE, `Benzylpenicillin inhib zone diam`:`Doxycycline inhib zone diam`) 

ZONE_data <- ZONE_data %>% mutate(Abx=str_replace(Abx, " inhib zone diam", "")) 

ZONE_data$ZONE <- as.numeric(ZONE_data$ZONE)

unique(ZONE_data$ZONE)

ZONE_data <- ZONE_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_mic_big , CLSI_mic_bigeq_m100 , CASFM_mic_big )))

ZONE_data <- ZONE_data %>% 
  mutate(EUCAST_Diam_Resist=ifelse(ZONE<EUCAST_disc_lower,1,0)) %>%
  mutate(CLIST_Diam_m100_Resist=ifelse(ZONE<CLSI_disc_lower_m100  ,1,0)) %>%
  mutate(CASFM_Diam_Resist=ifelse(ZONE<CASFM_disc_lower ,1,0)) 

summary_ZONE_concent <- ZONE_data %>% filter(!is.na(ZONE)) %>%
  group_by(`Species identification`, Abx) %>% 
  summarise(mean=mean(ZONE, na.rm=T), 
            sd=sd(ZONE), 
            mean=mean(ZONE), 
            median=median(ZONE), 
            Q1=quantile(ZONE, 0.25),
            Q3=quantile(ZONE, 0.75),
            n=n()) 


EUCAST_resist_counts_zone <- ZONE_data %>% filter(!is.na(EUCAST_Diam_Resist)) %>%
  group_by(`Species identification`, Abx, EUCAST_Diam_Resist) %>% count() %>%
  spread(key=EUCAST_Diam_Resist, value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_eucast=`1`+`0`, perc_r_eucast=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

casfm_resist_counts_zone <- ZONE_data %>% filter(!is.na(CASFM_Diam_Resist  )) %>%
  group_by(`Species identification`, Abx, CASFM_Diam_Resist ) %>% count() %>%
  spread(key=CASFM_Diam_Resist , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_casfm=`1`+`0`, perc_r_casfm=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))


fwrite(summary_ZONE_concent, "../out/summary_ZONE_concent_July_16.csv")

fwrite(ZONE_data, "../out/ZONE_data_July_16.csv")

fwrite(EUCAST_resist_counts_zone, "../out/EUCAST_resist_counts_zone_July_16.csv")
fwrite(casfm_resist_counts_zone, "../out/casfm_resist_counts_zone_July_16.csv")





# ---------------------


# Summary Table MIC and Diameters ----------------

MIC_workbook_jul_16 <- read_excel(path = "../out/mic_workbook_jul_16.xlsx",  sheet="Summary MIC Values", skip = 1)
unique(MIC_workbook_jul_16$`Species identification`)
TOP <-  MIC_workbook_jul_16 

MIC_workbook_jul_16 <- read_excel(path = "../out/mic_workbook_jul_16.xlsx",  sheet="EUCAST Resist Thresh", skip = 0)
MIC_workbook_jul_16 <- MIC_workbook_jul_16 %>% rename("N_EUCASAT_thre"="# Nr Isolates")

TOP <- TOP %>% rename("Abx"="Abx MIC tested") %>% left_join(MIC_workbook_jul_16)

MIC_workbook_jul_16 <- read_excel(path = "../out/mic_workbook_jul_16.xlsx",  sheet="CLSI 100 Resist Thresh", skip = 0)
MIC_workbook_jul_16 <- MIC_workbook_jul_16 %>% rename("N_CLSI_thre_m100"="# Nr Isolates") 

TOP <- TOP %>% left_join(MIC_workbook_jul_16)


MIC_workbook_jul_16 <- read_excel(path = "../out/mic_workbook_jul_16.xlsx",  sheet="CASFM Resist Thresh", skip = 0)
MIC_workbook_jul_16 <- MIC_workbook_jul_16 %>% rename("N_CASFM_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(MIC_workbook_jul_16 )

fwrite(TOP, "../out/MIC_Summary_All_Jul_16.csv")




diam_inhib_workbook_jul_16 <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="Summary Zone Diam Values", skip = 1)

TOP <- diam_inhib_workbook_jul_16

diam_inhib_workbook_jul_16  <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="EUCAST Resist Thresh", skip = 0)
diam_inhib_workbook_jul_16 <- diam_inhib_workbook_jul_16 %>% rename("N_EUCASAT_thre"="# Nr Isolates")

TOP <- TOP %>% rename("Abx"="Abx Zone tested") %>% left_join(diam_inhib_workbook_jul_16 )

diam_inhib_workbook_jul_16 <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="CASFM Resist Thresh", skip = 0)
diam_inhib_workbook_jul_16 <- diam_inhib_workbook_jul_16 %>% rename("N_CASFM_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(diam_inhib_workbook_jul_16)

data.frame(TOP)

fwrite(TOP, "../out/Diams_Summary_All_Jul_16.csv")



# -------------
# Overall figures ---------------------------

Look_up <- fread("../data/Look_up_MORAX_Complete_dataset.csv")

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

# isolates per year 

MORAX_Complete_dataset %>% select(`Code event`, `Species identification`, `2020=1`:`2024=1`) %>%
  filter(!is.na(`Species identification`)) %>%
  gather(Group, exp, `2020=1`:`2024=1`) %>%
  filter(exp==1) %>% group_by(`Code event`) %>% filter(Group==max(Group)) %>% ungroup() %>%
  group_by(Group) %>% count() %>% mutate(n=n/708)

#   Group      n
# 1 2020=1 0.120
# 2 2021=1 0.181
# 3 2022=1 0.216
# 4 2023=1 0.244
# 5 2024=1 0.239

# Isolates per country

data.frame(MORAX_Complete_dataset %>% select(`Code event`, `Species identification`, `Austria`:`Turkey`) %>%
             filter(!is.na(`Species identification`)) %>%
             gather(Group, exp, `Austria`:`Turkey`) %>%
             filter(exp==1) %>%
             group_by(Group) %>% count() %>% 
             arrange(n) %>% mutate(n=n/709) ) 

# 1        Slovakia 0.001410437
# 2          Turkey 0.001410437
# 3     Netherlands 0.002820874
# 4        Bulgaria 0.007052186
# 5     Switzerland 0.011283498
# 6         Croatia 0.015514810
# 7          Greece 0.016925247
# 8           Italy 0.018335684
# 9  Czech Republic 0.019746121
# 10         Poland 0.029619182
# 11         Norway 0.031029619
# 12        Denmark 0.033850494
# 13        Belgium 0.035260931
# 14        Austria 0.038081805
# 15       Portugal 0.039492243
# 16       Slovenia 0.040902680
# 17          Spain 0.060648801
# 18         Sweden 0.073342736
# 19        Germany 0.159379408
# 20        Ireland 0.174894217
# 21         France 0.188998590

# Medical ward

data.frame(MORAX_Complete_dataset %>% select(`Code event`, `Species identification`,  `Emergency=1`:`Neonatal ICU=1`) %>%
             filter(!is.na(`Species identification`)) %>%
             gather(Group, exp, `Emergency=1`:`Neonatal ICU=1`) %>%
             filter(exp==1) %>% distinct() %>% 
             group_by(Group) %>%
             count() %>%
             arrange(n) %>% mutate(n=n/683))


#                    Group           n
# 1         Neonatal ICU=1 0.007320644
# 2 Surgical paediatrics=1 0.008784773
# 3       Paediatric ICU=1 0.016105417
# 4        Surgical ward=1 0.048316252
# 5                  ICU=1 0.049780381
# 6  Medical Paediatrics=1 0.221083455
# 7            Emergency=1 0.320644217
# 8         Medical ward=1 0.327964861


# Most Common species

Look_up

TOP50 <- data.frame(MORAX_Complete_dataset %>%  filter(!is.na(`Species identification`)) %>% 
                      group_by(`Species identification`) %>% count() %>%
                      arrange(-n) %>% ungroup() %>% mutate(cum=cumsum(n)) %>% mutate(n=n/709)) 
  

top_20_species <- TOP50 %>% slice(1:20) %>% select(Species.identification) 

top_20_species <- list(top_20_species$Species.identification)[[1]]

# Most Common species by country TOP10 per country

temp <- data.frame(MORAX_Complete_dataset %>%  filter(!is.na(`Species identification`)) %>% 
                     gather(Group, exp, `Austria`:`Turkey`) %>%
                     filter(exp==1) %>%
                     mutate(`Species identification`=ifelse(`Species identification`%in%top_20_species, `Species identification`, "Other")) %>%
                     group_by(Group, `Species identification`) %>% count() %>%
                     rename("sub"="n") %>% ungroup() %>%
                     group_by(Group) %>% mutate(Tot=sum(sub)) %>%
                     mutate(perc=sub/Tot) %>%
                     arrange(Group, -perc ) %>%
                     group_by(Group) %>%  select(-c(sub, Tot)) %>%
                     spread(key=Group, value=perc))


fwrite(temp, "../out/temp.csv")

temp <- data.frame(MORAX_Complete_dataset %>%  filter(!is.na(`Species identification`)) %>%
              gather(Group, exp, `Emergency=1`:`Neonatal ICU=1`) %>%
             filter(exp==1) %>% distinct() %>% 
             mutate(`Species identification`=ifelse(`Species identification`%in%top_20_species, `Species identification`, "Other")) %>%
             group_by(Group,`Species identification`) %>%
             count() %>%
             rename("sub"="n") %>% ungroup() %>%
             group_by(Group) %>% mutate(Tot=sum(sub)) %>%
             mutate(perc=sub/Tot) %>%
             arrange(Group, -perc ) %>%
             group_by(Group) %>% select(-c(sub, Tot)) %>%
             spread(key=Group, value=perc))

fwrite(temp, "../out/temp.csv")



# --------------------

# EUCAST Plot species vs resistance rate MIC ---------------------

Look_up <- fread("../data/Look_up_MORAX_Complete_dataset.csv")

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

names(MORAX_Complete_dataset)


MIC_data <- MORAX_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MORAX_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Benzylpenicillin MIC`:`Doxycycline MIC`) 

MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="≤1",0.1, MIC)) %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


temp <- MIC_data %>% select(`Code event`, `Species identification`, Abx, EUCAST_mic_Resist)

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "EUCAST MIC \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_mic_eucast.svg", plot=plot, width=7, height=7)

resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, EUCAST_mic_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "EUCAST MIC \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_mic_eucast.svg", plot=plot, width=8, height=8)




# ----------------
# CLSI m100 Plot species vs resistance rate MIC ---------------------


Look_up <- fread("../data/Look_up_MORAX_Complete_dataset.csv")

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

names(MORAX_Complete_dataset)

MIC_data <- MORAX_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MORAX_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Benzylpenicillin MIC`:`Doxycycline MIC`) 

MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="≤1",0.1, MIC)) %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


temp <- MIC_data %>% select(`Code event`, `Species identification`, Abx, CLIST_mic_m100_Resist)

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CLIST_mic_m100_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "CLSI m100 MIC \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_mic_clsi_m100.svg", plot=plot, width=7, height=7)

resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, CLIST_mic_m100_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CLIST_mic_m100_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "CLSI m100 MIC \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_mic_clsi_m100.svg", plot=plot, width=8, height=8)




# ----------------
# CASFM Plot species vs resistance rate MIC ---------------------



Look_up <- fread("../data/Look_up_MORAX_Complete_dataset.csv")

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

names(MORAX_Complete_dataset)

MIC_data <- MORAX_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MORAX_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Benzylpenicillin MIC`:`Doxycycline MIC`) 


MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>% mutate(MIC=ifelse(MIC=="≤1",0.1, MIC)) %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


temp <- MIC_data %>% select(`Code event`, `Species identification`, Abx, CASFM_mic_Resist)

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CASFM_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "CASFM MIC \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_mic_casfm.svg", plot=plot, width=7, height=7)

resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, CASFM_mic_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CASFM_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "CASFM MIC \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_mic_casfm.svg", plot=plot, width=8, height=8)




# ----------------
# Compare % Resistance Across countries ---------------------------

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

MORAX_Complete_dataset <- MORAX_Complete_dataset  %>% 
  mutate(`Code event`=str_replace_all(`Code event`, "Ä°", "I"))  %>% 
  mutate(`Code event`=str_replace_all(`Code event`, "İ", "I"))

Countries_map <- MORAX_Complete_dataset %>% select(`Code event`, Austria:Turkey) %>%
  gather(Country, ON, Austria:Turkey) %>%
  filter(ON==1) %>% distinct() %>% select(-ON)


mic_workbook_jul_16 <- read_xlsx(path="../out/mic_workbook_jul_16.xlsx",sheet = "MIC Clean Data", skip=0, col_types = "text", trim_ws = TRUE)

mic_workbook_jul_16 <- mic_workbook_jul_16 %>% 
  mutate(`Code event`=str_replace_all(`Code event`, "Ä°", "I")) 

mic_workbook_jul_16 <- mic_workbook_jul_16 %>% select(`Code event`, `Species identification`, Abx, 
                                                      `EUCAST Resistant?`, `CLIST m100 Resistant?`, 
                                                      `CASFM Resistant?`) 

# Select the 4 columns
resistance_cols <- c("EUCAST Resistant?", "CLIST m100 Resistant?", "CASFM Resistant?")

# Convert to numeric: assume "Yes" = 1, "No" = 0, NA stays NA
mic_workbook_jul_16[resistance_cols] <- lapply(mic_workbook_jul_16[resistance_cols], function(x) {
  ifelse(x == "1", 1, ifelse(x == "0", 0, NA))
})

mic_workbook_jul_16$Resistance_Sum <- rowSums(mic_workbook_jul_16[resistance_cols], na.rm = TRUE)


mic_workbook_jul_16$All_Resistance_NA <- apply(
  mic_workbook_jul_16[resistance_cols],
  1,
  function(x) all(is.na(x))
)

mic_workbook_jul_16 <- mic_workbook_jul_16 %>% filter(All_Resistance_NA==FALSE) %>%
  select(-All_Resistance_NA)

mic_workbook_jul_16$Resistance_NonNA_Count <- rowSums(!is.na(mic_workbook_jul_16[resistance_cols]))


mic_workbook_jul_16 <- mic_workbook_jul_16 %>% mutate(perc=Resistance_Sum/Resistance_NonNA_Count) %>%
  select(`Code event`, `Species identification`, Abx,Resistance_Sum , Resistance_NonNA_Count , perc)

mic_workbook_jul_16 <- mic_workbook_jul_16 %>% left_join(Countries_map)


temp <- data.frame(mic_workbook_jul_16 %>%
  group_by(Country, `Species identification`, Abx) %>%
  summarise(den=sum(Resistance_NonNA_Count), num=sum(Resistance_Sum)) %>%
  ungroup() %>% mutate(perc=num/den) %>% filter(den>10) %>%
    select(-den, -num) %>%
  spread(key=Country, value=perc))

fwrite(temp, "temp.csv")

# ---------------
# Summary table >10 isolates overall ---------
mic_summary_jul_16 <- read_xlsx(path="../out/summaries_overall_jul_16.xlsx",sheet = "mic_summary_jul_16", skip=0, col_types = "text", trim_ws = TRUE)
diam_inhib_summary_jul_16 <- read_xlsx(path="../out/summaries_overall_jul_16.xlsx",sheet = "diam_inhib_summary_jul_16", skip=0, col_types = "text", trim_ws = TRUE)

temp <- mic_summary_jul_16 %>% 
  filter(N_EUCASAT_thre>10) %>%
  select(Species, Abx, Median, Q90, `% Resist EUCAST Threshold`) %>%
  full_join(
    mic_summary_jul_16 %>% 
    filter(N_CLSI_thre_m100 >10) %>%
    select(Species, Abx, Median, Q90, `% Resist CLSI Threshold m100`) 
  ) %>%
  full_join(
    mic_summary_jul_16 %>% 
    filter(N_CASFM_thre  >10) %>%
    select(Species, Abx, Median, Q90, `% Resist CASFM Threshold`) 
  ) %>%
  mutate(`% Resist EUCAST Threshold`=as.numeric(`% Resist EUCAST Threshold`)) %>%
  mutate(`% Resist CLSI Threshold m100`=as.numeric(`% Resist CLSI Threshold m100`)) %>%
  mutate(`% Resist CASFM Threshold`=as.numeric(`% Resist CASFM Threshold`))  %>%
  full_join(
    diam_inhib_summary_jul_16 %>% 
    filter(N_EUCASAT_thre>10) %>%
    select(Species, Abx, Median, Q1, Q3, `% Resist EUCAST Threshold`) %>% 
      rename("Diam_Median"="Median", "Diam_% Resist EUCAST Threshold"="% Resist EUCAST Threshold") %>%
      mutate(IQR=as.numeric(Q3)-as.numeric(Q1)) %>% select(-Q1, -Q3)
  ) %>%
  full_join(
    diam_inhib_summary_jul_16 %>% 
    filter(N_CASFM_thre  >10) %>%
    select(Species, Abx, Median, Q1, Q3, `% Resist CASFM Threshold`) %>% 
      rename("Diam_Median"="Median", "Diam_% Resist CASFM Threshold"="% Resist CASFM Threshold") %>%
      mutate(IQR=as.numeric(Q3)-as.numeric(Q1)) %>% select(-Q1, -Q3)
  )

fwrite(temp, "temp.csv")

# ------------


# Overall resistance -----------


MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

unique(MORAX_Complete_dataset$`Medical Paediatrics=1`)
unique(MORAX_Complete_dataset$`Surgical paediatrics=1`)
unique(MORAX_Complete_dataset$`Paediatric ICU=1`)


peds_codes <- MORAX_Complete_dataset %>% filter(`Medical Paediatrics=1`=="1"| 
                                    `Surgical paediatrics=1`=="1"|
                                    `Paediatric ICU=1`=="1") %>%
  select(`Code event`) %>% distinct()
  


MIC_Clean_Data <- read_xlsx(path="../out/mic_workbook_jul_16.xlsx",sheet = "MIC Clean Data", skip=0, col_types = "text", trim_ws = TRUE)


Inhib_Zone_Diam_Clean <- read_xlsx(path="../out/diam_inhib_workbook_jul_16.xlsx",sheet = "Inhib Zone Diam Clean Data", skip=0, col_types = "text", trim_ws = TRUE)


result <- MIC_Clean_Data %>% select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`, `CLIST m100 Resistant?`) %>%
  gather(Test, Result, , `EUCAST Resistant?`:`CLIST m100 Resistant?`) %>%
  bind_rows(
    Inhib_Zone_Diam_Clean %>% select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`, ) %>%
  gather(Test, Result, , `EUCAST Resistant?`:`EUCAST Resistant?`)
  ) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>% summarise(n=n(), res=sum(as.numeric(Result) )) %>%
  mutate(perc=res/n) %>% select(-res)

result <- result %>% arrange(`Species identification`, Abx)

result

fwrite(result, "overall_aeromonas_resistance_rates.csv")

# ------

# Compare adults vs pediatrics EUCAST ---------

MORAX_Complete_dataset <- read_xlsx(path="../data/2.MORAX-EU_dataset.xlsx",sheet = "1. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

unique(MORAX_Complete_dataset$`Medical Paediatrics=1`)
unique(MORAX_Complete_dataset$`Surgical paediatrics=1`)
unique(MORAX_Complete_dataset$`Paediatric ICU=1`)


peds_codes <- MORAX_Complete_dataset %>% filter(`Medical Paediatrics=1`=="1"| 
                                    `Surgical paediatrics=1`=="1"|
                                    `Paediatric ICU=1`=="1") %>%
  select(`Code event`) %>% distinct()
  


MIC_Clean_Data <- read_xlsx(path="../out/mic_workbook_jul_16.xlsx",sheet = "MIC Clean Data", skip=0, col_types = "text", trim_ws = TRUE)



Inhib_Zone_Diam_Clean <- read_xlsx(path="../out/diam_inhib_workbook_jul_16.xlsx",sheet = "Inhib Zone Diam Clean Data", skip=0, col_types = "text", trim_ws = TRUE)


result_peds <- MIC_Clean_Data %>% select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`) %>%
  bind_rows(
    Inhib_Zone_Diam_Clean %>% select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`) 
  ) %>% drop_na() %>% inner_join(peds_codes) %>%
  group_by(`Species identification`, Abx) %>% summarise(n=n(), res=sum(as.numeric(`EUCAST Resistant?`) )) %>%
  mutate(perc=res/n) %>% select(-res)

result_peds <- result_peds %>% arrange(`Species identification`, Abx)

result_peds



result_adults <- MIC_Clean_Data %>% select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`) %>%
  bind_rows(
    Inhib_Zone_Diam_Clean %>% select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`) 
  ) %>% drop_na() %>% anti_join(peds_codes) %>%
  group_by(`Species identification`, Abx) %>% summarise(n=n(), res=sum(as.numeric(`EUCAST Resistant?`) )) %>%
  mutate(perc=res/n) %>% select(-res)

result_adults <- result_adults %>% arrange(`Species identification`, Abx)

result_adults














library(dplyr)
library(tidyr)
library(purrr)

# --- PEDIATRICS ---
peds_summary <- MIC_Clean_Data %>%
  select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`) %>%
  bind_rows(
    Inhib_Zone_Diam_Clean %>%
      select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`)
  ) %>%
  drop_na() %>%
  inner_join(peds_codes, by = "Code event") %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_peds = n(),
    res_peds = sum(as.numeric(`EUCAST Resistant?`)),
    .groups = "drop"
  ) %>%
  mutate(susc_peds = n_peds - res_peds)

# --- ADULTS ---
adult_summary <- MIC_Clean_Data %>%
  select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`) %>%
  bind_rows(
    Inhib_Zone_Diam_Clean %>%
      select(`Code event`, `Species identification`, Abx, `EUCAST Resistant?`)
  ) %>%
  drop_na() %>%
  anti_join(peds_codes, by = "Code event") %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_adults = n(),
    res_adults = sum(as.numeric(`EUCAST Resistant?`)),
    .groups = "drop"
  ) %>%
  mutate(susc_adults = n_adults - res_adults)

# --- MERGE ---
comparison <- full_join(
  peds_summary,
  adult_summary,
  by = c("Species identification", "Abx")
)

# --- FISHER TEST ---
comparison <- comparison %>%
  rowwise() %>%
  mutate(
    fisher_p = {
      mat <- matrix(
        c(
          res_peds,
          susc_peds,
          res_adults,
          susc_adults
        ),
        nrow = 2
      )

      if(any(is.na(mat))) {
        NA_real_
      } else {
        fisher.test(mat)$p.value
      }
    }
  ) %>%
  ungroup()

comparison

fwrite(comparison, "morax_comparison_adults_peds_eucast.csv")


# -------------------