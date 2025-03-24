#####HPLC data #####


#Load the libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)


#Read the data table
hplc <- read_delim("~/Documents/HPLC/Rade-Micro_HPLC_MLP_2012-2023/data-Tableau 1.csv",
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

hplc<- data_Tableau_1

#Basic inspection
summary(hplc)

hplc<- as.data.frame(hplc)


#Change the disposition of the data table
hplc <- hplc %>%
  mutate(across(
    chlc3:q_anthera,
    ~ as.numeric(gsub(",", ".", .))      # Replace comma by points and point in numeric format
  ))

head(hplc)

#####First: inspect quality of data by looking at chla vs TAP to look at the qc quality#####
hplc <- hplc %>%
  mutate(TAP = rowSums(across(c(chlc3, chlc2c1, tchlb, peri, fuco, hex, but, tcar, viola, diadino, diato, zea)), na.rm = TRUE),
         ind=chlda / tchla)

hplc<- hplc%>% mutate(
  TAP_log= log10(TAP),
  tchla_log=log10(tchla))

hplctest<-hplc%>%filter(depth %in% c(75))


#Plot the inspection plot
ggplot()+ geom_point(aes(x=tchla_log,y=TAP_log), data=hplc)+ theme_bw()


ggplot(hplc, aes(x = tchla_log, y = TAP_log)) +
  geom_point(aes(color = as.factor(q_tchla))) +  # points colored by tchla
  stat_smooth(method = "lm", color = "red", se = FALSE) +
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x) +
  stat_cor(aes(label = ..rr.label..), label.x.npc = "center", label.y.npc = "top") +
  theme_bw() +
  labs(color = "Q_tchla")

#check something
count_qc2 <- hplc %>%
  filter(q_tchla == 3) %>%
  summarize(count = n())

#With everything
hplc_long <- hplc %>%
  pivot_longer(
    cols = starts_with("q_"),
    names_to = "q_variable",
    values_to = "q_value"
  ) %>%
  mutate(q_value = as.factor(q_value))

#Some plots to look at the quality of data

ggplot(hplc_long, aes(x = tchla_log, y = TAP_log)) +
  geom_point(aes(color = q_value)) +
  stat_smooth(method = "lm", color = "red", se = FALSE) +
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x) +
  facet_wrap(~ q_variable, scales = "free_y") +
  theme(strip.background=element_rect(fill="lightblue"))+
  labs(color = "QC Value", title = "TAP vs Tchla with points colored by qc value of one pigment")

#Show only bad QC
ggplot(hplc_long, aes(x = tchla_log, y = TAP_log)) +
  geom_point(data = hplc_long %>% filter(q_value != "1"), aes(color = q_value)) +
  stat_smooth(method = "lm", color = "red", se = FALSE) +
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x) +
  facet_wrap(~ q_variable, scales = "free_y") +
  theme(strip.background = element_rect(fill = "lightblue")) +
  labs(color = "QC Value", title = "TAP vs Tchla with points colored by qc value of one pigment")

#Get how many qc==3 per variable
hplc_long <- hplc %>%
  pivot_longer(cols = starts_with("q_"),
               names_to = "q_variable",
               values_to = "q_value")

# Number of values where bad qc per variable
count_qc3_by_variable <- hplc_long %>%
  filter(q_value == 3) %>%
  group_by(q_variable,depth) %>%
  summarize(count = n(), .groups = "drop")

# Print the data_table
count_qc3_by_variable

# Number of values where bad qc per variable
count_qc0_by_variable <- hplc_long %>%
  filter(q_value == 0) %>%
  group_by(q_variable) %>%
  summarize(count = n(), .groups = "drop")

# Print the data_table
count_qc0_by_variable

#Check the degradation ratio
# Visualize the distribution of ratio per QC
ggplot(hplc_long, aes(x = q_value, y = ind, fill = as.factor(q_value))) +
  geom_boxplot(outlier.color = "red") +
  labs(
    x = "QC Value",
    y = "Degradation Ratio (chlorophyllide_a / TChla)",
    title = "Distribution of Degradation Ratio by QC Value"
  ) +
  theme_minimal()

#Look at the degradation ratio 
ggplot(hplc_long, aes(x = ind, fill = as.factor(q_value))) +
  geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7, position = "identity")  +
  facet_wrap(~ q_variable, scales = "free_y") +
  labs(
    x = "Degradation Ratio (chlorophyllide_a / TChla)",
    y = "Count",
    fill = "QC Value",
    title = "Histogram of Degradation Ratio by Variable and Colored by QC"
  ) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightblue"))

#####Second step: keep only good QC##### --> Remove QC=3 but keep QC=0 #####

hplcd <- hplc%>%
  select(-starts_with("q")) %>%
  pivot_longer(-c(date, latitude, longitude, file, depth, ...58:...66))%>%select(-c(...58:...66))%>%mutate(value=as.numeric(value))
hplcc <- hplc %>%
  select(date, depth, starts_with("q")) %>%
  pivot_longer(-c(date, depth)) %>%
  mutate(name=name %>% str_remove("^q") %>% str_remove("^_")) %>%
  rename(qcode=value)
# combine
#hplccomb <- inner_join(hplcd, hplcc) |>
#  mutate(
#   value=if_else(qcode %in% c(0, 3), NA, value)
#  ) %>%distinct()

hplccomb <- inner_join(hplcd, hplcc) %>%
  mutate(
    value=if_else(qcode==3, NA, value)
  ) %>%distinct()

hplccomb<-hplccomb%>%select(-qcode)%>%pivot_wider(names_from=name, values_from=value)

#####Pre-cleaning of data#####
hplccomb <- hplccomb %>% mutate(zea = ifelse(zea == 0, zea + lut, zea),
                                tchlb=ifelse(tchlb == 0 & (chlb != 0 | dvchlb != 0), chlb + dvchlb, tchlb))


##### Third step: Estimation of coeff: regression to look at diagnostic pigments #####

#####First way: Find my own coefficients #####
#hplccomb_depth<- hplccomb%>%filter(depth%in%c(0,1))%>%
#  mutate(date = as.Date(date, format = "%m/%d/%y"))
#hplccomb_depth %>%
#  mutate(date = as.Date(date, format = "%m/%d/%y"))%>%ggplot+geom_point(aes(x=date,y=tchla))+scale_y_log10()

#Test with tchla
#x <- hplccomb_depth%>%select(tchla,date)
#dec <- stlplus(x$tchla, x$date, n.p=26, s.window="periodic", t.window=10*26)
#par(cex.lab = 1.5,   
#    cex.axis = 1.3,  
#    cex.main = 1.7,  
#    cex.sub = 1.3)
#plot(dec)
# -> OK

#Try to integrate
hplccomb_int <- hplccomb %>%
  select(-depth) %>%
  group_by(date) %>%
  summarize(across(chlc3:anthera, ~ mean(., na.rm = TRUE)), .groups = "drop")


#Multiple regression
reg_multiple <- lm(tchla ~ fuco+peri+hex+but+allo+tchlb+zea, data = hplccomb)

summary(reg_multiple)

#When integrated
reg_multiple_integrate<-lm(tchla ~ fuco+peri+hex+but+allo+tchlb+zea, data = hplccomb_int)

summary(reg_multiple_integrate)

#Graph
ggplot(hplccomb, aes(x = fuco + peri + hex + but + allo + tchlb + zea, y = tchla))  + geom_point()+
  stat_smooth(method = "lm", color = "red", se = FALSE) +
  labs(x = "Explanatory variables", y = "Total chlorophyll a (tchla)", title = "Multiple regression") +
  theme_minimal() +
  theme(legend.position = "top") + scale_x_log10()+scale_y_log10()

#####Fourth step: calculation of size classes and taxonomic group with the three methods to have the coefficients#####


#With the coefficients of the multiple regression with all the points

hplccomb<- hplccomb%>%mutate(SDPw=1.87*fuco+ 1.70*peri+1.28*hex+0.53*but+3.44*allo+1.18*tchlb+1.46*zea,
                             fmicro=(1.87*fuco+1.70*peri)/SDPw,
                             fnano=(1.28*hex+0.53*but+3.44*allo)/SDPw,
                             fpico=(1.18*tchlb+1.46*zea)/SDPw,
                             micro_chla=fmicro*tchla,
                             nano_chla=fnano*tchla,
                             pico_chla=fpico*tchla
                             
)

######Taxonomic information######
hplccomb<- hplccomb%>% mutate(diatoms=(1.87*fuco)/SDPw,
                              dinoflagellates=(1.70*peri)/SDPw,
                              green_algae=(1.18*chlb)/SDPw,
                              prokaryotes=(1.46*zea)/SDPw,
                              prochlorococcus=(dvchla)/tchla, 
                              micro=(fuco+peri)/tchla,
                              pico=(tchlb+zea)/tchla
)
######Final data table#####


hplc_final<- hplccomb%>%select(date, file, depth, SDPw, fmicro, fnano, fpico, micro_chla, nano_chla, pico_chla, diatoms, dinoflagellates, green_algae, prokaryotes, prochlorococcus, tchla)#%>%mutate(tchla=log(tchla))

save(hplc_final, file="hplc_final.Rdata")