library(tidyverse)
library(rnaturalearth)
library(sf)
library(wesanderson)


world <- ne_countries(scale="large", returnclass = "sf") 

head(world)

world %>% filter(admin %in% c("Las Palmas")) %>% select(admin) %>% distinct()

unique(world$admin)

data.frame(world %>% filter(admin %in% c("Portugal", "Spain", "France", "Romania", "Netherlands", 
                              "Denmark", "Italy", "Germany", "Hungary", "Switzerland", "Croatia", 
                              "Greece", "Turkey", "Poland")) %>%
  select(admin) %>% distinct())

world <- world %>% filter(admin %in% c("Portugal", "Spain", "France", "Romania", "Netherlands", 
                              "Denmark", "Italy", "Germany", "Hungary", "Switzerland", "Croatia", "Greece", "Turkey", "Poland"))

target_crs <- "+proj=moll"

world_moll <- world %>% st_transform(crs = target_crs)

world_moll <- world_moll %>% mutate(beds=ifelse(admin=="Portugal", 1170+1150,
                                  ifelse(admin=="Spain", 1260+850+1400+660,
                                         ifelse(admin=="France", 1500+830+1960+2900,
                                                ifelse(admin=="Romania", 690+1170,
                                                       ifelse(admin=="Netherlands",1730,
                                                              ifelse(admin=="Denmark",1150,
                                                                     ifelse(admin=="Italy", 450+870+1900,
                                                                            ifelse(admin=="Germany", 830+2200,
                                                                                   ifelse(admin=="Hungary",2300,
                                                                                          ifelse(admin=="Switzerland", 950,
                                                                                                 ifelse(admin=="Croatia", 1420+2100,
                                                                                                        ifelse(admin=="Greece",770+800,
                                                                                                               ifelse(admin=="Turkey",370,
                                                                                                                      ifelse(admin=="Poland",910,NA)))))))))))))))



world_moll <- world_moll %>% mutate(cases=ifelse(admin=="Portugal", 874,
                                  ifelse(admin=="Spain", 1260,
                                         ifelse(admin=="France", 1718,
                                                ifelse(admin=="Romania", 143,
                                                       ifelse(admin=="Netherlands",234,
                                                              ifelse(admin=="Denmark",188,
                                                                     ifelse(admin=="Italy", 291,
                                                                            ifelse(admin=="Germany", 336,
                                                                                   ifelse(admin=="Hungary",291,
                                                                                          ifelse(admin=="Switzerland", 221,
                                                                                                 ifelse(admin=="Croatia", 345,
                                                                                                        ifelse(admin=="Greece",253,
                                                                                                               ifelse(admin=="Turkey",61,
                                                                                                                      ifelse(admin=="Poland",162,NA)))))))))))))))


sum(is.na(world_moll$beds))

world_moll %>%
  ggplot() +
  geom_sf(aes(fill=beds)) +
  scale_fill_viridis_c() +
  scale_x_continuous(labels=function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels=function(x) paste0(x, '\u00B0', "N")) +
  theme_minimal() +
  theme(panel.background=element_rect((fill="aliceblue"))) 


window_coord <- st_sfc(
  st_point(c(-21.45,24.05)),
  st_point(c(47.46,59.53)),
  crs=4326
)

window_coord_sf <- window_coord %>%
  st_transform(crs = target_crs) %>%
  st_coordinates()

world_moll <- world_moll %>% select(!contains("fclass"))

world_moll %>%
  ggplot() +
  geom_sf(aes(fill=log(beds))) +
  coord_sf(xlim=window_coord_sf[,"X"], ylim=window_coord_sf[,"Y"], expand=FALSE) +
  scale_fill_gradientn(colours = wes_palette("Zissou1", 10000, type = "continuous")) +
  #scale_fill_viridis_c(option="C") +
  scale_x_continuous(labels=function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels=function(x) paste0(x, '\u00B0', "N")) +
  theme_minimal() +
  theme(panel.background=element_rect((fill="white"))) 


world_moll %>%
  ggplot() +
  geom_sf(aes(fill=log(cases))) +
  coord_sf(xlim=window_coord_sf[,"X"], ylim=window_coord_sf[,"Y"], expand=FALSE) +
  scale_fill_gradientn(colours = wes_palette("Zissou1", 10000, type = "continuous")) +
  #scale_fill_viridis_c(option="C") +
  scale_x_continuous(labels=function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels=function(x) paste0(x, '\u00B0', "N")) +
  theme_minimal() +
  theme(panel.background=element_rect((fill="white"))) 

  
  
  
world_moll %>%
  ggplot() +
  geom_sf(aes(fill=cases/beds)) +
  coord_sf(xlim=window_coord_sf[,"X"], ylim=window_coord_sf[,"Y"], expand=FALSE) +
  scale_fill_gradientn(colours = wes_palette("Zissou1", 10000, type = "continuous")) +
  #scale_fill_viridis_c(option="C") +
  scale_x_continuous(labels=function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels=function(x) paste0(x, '\u00B0', "N")) +
  theme_minimal() +
  theme(panel.background=element_rect((fill="white"))) 

data.frame(world_moll %>% select(admin, cases, beds) %>% mutate(den=cases/beds))
