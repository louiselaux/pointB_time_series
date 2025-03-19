#Load the libraries
library(tidyverse)
library(dplyr)
library(ggplot2)

#Load the datasets
load(file="output/dl.Rdata")
load(file="output/dstl.Rdata")
load(file="output/stats_env.Rdata")


#Figure 1: Inspection plots of the time series#

#Inspection plot when depth is 1 meter
dl%>%filter(depth=="1")%>%filter(date>"1967-01-01")%>%filter(name%in%c("CHLA","T","S","O","NO2","NO3","Sigma","COP","NOP","PO4","SIOH4"))%>%ggplot() + facet_wrap(~name, scales="free_y") +
  geom_line(aes(x=date, y=value), size=0.2)+  theme_bw() +  
  theme(
    axis.title = element_text(size = 14),        
    axis.text = element_text(size = 12),         
    strip.text = element_text(size = 16),        
    legend.title = element_text(size = 14),      
    legend.text = element_text(size = 12),       
    plot.title = element_text(size = 18, hjust = 0.5),  
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "lightblue")
  ) +
  labs(x = "Date", y = "Value")


#Inspection plot when taking the mean of each variable
dl%>%group_by(date,name)%>%summarise(value=mean(value, na.rm=TRUE))%>%filter(name%in%c("CHLA","T","S","O","NO2","NO3","Sigma","COP","NOP","PO4","SIOH4"))%>%ggplot() + facet_wrap(~name, scales="free_y") +
  geom_line(aes(x=date, y=value), size=0.2)+  theme_bw() +  
  theme(
    axis.title = element_text(size = 14),        
    axis.text = element_text(size = 12),         
    strip.text = element_text(size = 16),        
    legend.title = element_text(size = 14),      
    legend.text = element_text(size = 12),      
    plot.title = element_text(size = 18, hjust = 0.5),  
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "lightblue")
  ) +
  labs(x = "Date", y = "Valeur")



#Figure 2: Long-term trend with environment time series #

#The plot
ggplot() +
  facet_wrap(name~., scales="free_y", ncol=2) +
  geom_path(aes(target_date, deseason), data=dstl %>% filter(name %in% c("CHLA","T","NO3","SIOH4","S")), colour="grey20") +
  geom_abline(aes(slope=gls_slope, intercept=gls_intercept), data=stats%>%filter(name %in% c("CHLA","T","NO3","SIOH4","S") & gls_signif %in% c("*", "**", "***")), colour="red", size=0.75, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) + #when plotting 2nd line for recent years
  xlab("Date") +
  ggtitle("Long term trends") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),        
    axis.text = element_text(size = 12),         
    strip.text = element_text(size = 16),        
    legend.title = element_text(size = 14),      
    legend.text = element_text(size = 12),      
    plot.title = element_text(size = 18, hjust = 0.5),  
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "lightblue")
  )



