####################################################
####        Variabili climatiche Asromonte NP 2021.#
####        F. Cerini, Feb 2026                ####
###         Exploration of data by Daniele delle Monache####
####################################################

rm(list = ls())


#packages
require(tidyverse)
require(data.table)
require(openxlsx)
require(ggpubr)
require(corrplot)

## Load site environmental data ###
env <- read.xlsx("Data/Sites_variables_final_clc.xlsx") %>%
  mutate(
    across(where(is.character), as.factor),   # convert all character columns to factors
    across(where(is.logical), as.factor)      # (optional) also convert logicals to factors if any
  ) %>% 
  #add a vegetation code column so we don't have long description
  mutate(Code_veg = factor(Vegetation.type, 
                           levels = unique(Vegetation.type), 
                           labels = LETTERS[1:length(unique(Vegetation.type))])) 

str(env)
colnames(env)

ggscatter(env ,
          x = "T_media_Feb_Apr_90gg", 
          y = "T_Mean_Feb_Apr",
          add = "reg.line", size = 3, shape = 1)+
  stat_cor()+
  stat_regline_equation(label.y = 8)

ggplot(env %>% select(Site,
                      T_media_Feb_Apr_90gg,
                      T_Mean_Feb_Apr) %>% rename(T_local = T_media_Feb_Apr_90gg,
                                                 T_macro = T_Mean_Feb_Apr) %>% 
         pivot_longer(cols = T_local:T_macro,
                      names_to = "T_model"),
       aes(x = Site, y = value, fill = T_model))+
  geom_bar(stat = "identity", alpha = 0.8, 
           position = "identity")+
  theme_bw()+
  ggtitle("From Elisa")

hist(env$T_media_Feb_Apr_90gg, breaks = 10)
hist(env$T_Mean_Feb_Apr, breaks = 10)

ggscatter(env, x = "T_Mean_Feb_Apr", 
          y = "Soil_humidity",
          add = "reg.line", size = 3, shape = 1)+
  stat_cor()+
  stat_regline_equation(label.y = 8)

ggscatter(env, x = "Altitude", 
          y = "T_Mean_Feb_Apr",
          add = "reg.line", size = 3, shape = 1)+
  stat_cor()+
  stat_regline_equation(label.y = 8)
colnames(env)
ggscatter(env %>% 
            filter(!Site %in% c("S10","S14","S15","S18",
                                "S20","S6","S3")), x = "Altitude", 
          y = "T_media_Feb_Apr_90gg",
          add = "reg.line", size = 3, shape = 1)+
  stat_cor()+
  stat_regline_equation(label.y = 8)


# Load Rdata with raw temperatures made by DDM
load("Data/ws_pred.RData")
raw_t_ws<-as.data.frame(pred, row.names = NULL) %>% 
  rownames_to_column("Date") %>% 
  rename_with(~ paste0("S", .), -Date)

raw_t_ws_long<-raw_t_ws %>% pivot_longer(cols = S1:S22, names_to = "Site") %>% 
  mutate(Date = as.Date(Date))

ggplot(raw_t_ws_long %>% filter(between(Date,
                                    as.Date("2021-02-01"),
                                    as.Date("2021-04-30"))),
       aes(x = Date, 
           y = value,
           col = Site))+
  geom_point()+
  facet_wrap(~Site)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Weather stations")

raw_t_logg<- read.xlsx("Data/logger_dati_mean_2020_2021.xlsx", 
                       detectDates = T)

raw_t_logg_long<-raw_t_logg %>% pivot_longer(cols = S1:S22, names_to = "Site")

ggplot(raw_t_logg_long %>% 
         filter(between(date,
                        as.Date("2021-02-01"),
                        as.Date("2021-04-30"))),
       aes(x = date, 
           y = value,
           col = Site))+
  geom_point()+
  facet_wrap(~Site)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Data logger")


str(raw_t_logg)
colnames(raw_t_logg)
unique(raw_t_logg$date)

t_ws_feb_apr<-raw_t_ws %>% filter(between(Date,
                                                  as.Date("2021-02-01"),
                                                  as.Date("2021-04-30"))) %>% 
    pivot_longer(cols = S1:S22, names_to = "Site") %>% group_by(Site) %>% 
    summarise(Mean_t_ws_feb_apr = mean(value),
              sd_t_ws_feb_apr = sd(value)) %>% 
    filter(!Site %in% c("S13", "S17"))

t_ws_jul_aug<-raw_t_ws %>% filter(between(Date,
                                          as.Date("2021-07-01"),
                                          as.Date("2021-08-31"))) %>% 
  pivot_longer(cols = S1:S22, names_to = "Site") %>% group_by(Site) %>% 
  summarise(Mean_t_ws_jul_aug = mean(value),
            sd_t_ws_jul_aug = sd(value)) %>% 
  filter(!Site %in% c("S13", "S17"))

t_ws_corr<-left_join(t_ws_feb_apr,t_ws_jul_aug, by = "Site")

png("Results/Figures paper/Official/MACRO_SCALE_ONLY/temperatures_corr_supp.png",
    width = 3000, height = 2000, res = 300)

ggscatter(t_ws_corr, x = "Mean_t_ws_feb_apr", 
          y = "Mean_t_ws_jul_aug",
          add = "reg.line", size = 3, shape = 1)+
  stat_cor()

dev.off()


write.xlsx(t_ws_feb_apr, file = "Data/temperature_mean_sd_WS_feb_apr2021.xlsx")


colnames(raw_t_logg)
t_logg_feb_apr<-raw_t_logg %>% filter(between(date,
                                             as.Date("2021-02-01"),
                                             as.Date("2021-04-30"))) %>% 
  pivot_longer(cols = S1:S22, names_to = "Site") %>% group_by(Site) %>% 
  summarise(Mean_t_logg_feb_apr = mean(value)) %>% 
  filter(!Site %in% c("S13", "S17"))


t_df<-left_join(t_logg_feb_apr,
                t_ws_feb_apr,
                by = "Site") %>% 
  pivot_longer(Mean_t_logg_feb_apr:Mean_t_ws_feb_apr, names_to = "Estimate",
               values_to = "Mean_Temp")


ggplot(t_df, 
       aes(x = Site, y = Mean_Temp, fill = Estimate))+
  geom_bar(stat = "identity", 
           position = "identity",
           alpha = 0.9)+
  theme_bw()+
  ggtitle("From raw_data")




