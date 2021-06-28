setwd("/Users/katebuckeridge/Dropbox/R/UGrass/CNstab")


library(ggplot2)
library (ggExtra)
#library(ggpubr)
#library(plyr)

climate <- read.csv("field_climate.csv", header=TRUE, sep=",")

pltClim <- ggplot(climate, aes(x=Point,alpha=Data)) + 
  geom_bar(aes(y=Rain),
           fill = "grey",
           colour = "grey",
           position= position_dodge(width=0.9), 
           stat="identity") + 
  geom_line(aes(y=Max*5),
           colour=c("red"), 
            size = 2) +
  geom_line(aes(y=Min*5), 
            colour=c("blue"), 
            size = 2) +
  scale_alpha_manual(values = c(0.3,1)) +
  scale_x_discrete(name="", 
                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous(name = "Rain (mm)", 
                     sec.axis = sec_axis(~./5, name =expression(Temperature~(degree*C)))) + 
  theme(legend.position = "none",
        text = element_text(family = "Helvetica", 
                            color="black", 
                            size=12)) +
  theme_bw( ) 

pltClim

#### simpler version

climate2 <- read.csv("field_climate2.csv", header=TRUE, sep=",")
climate2$Date <- as.Date(climate2$Date,format='%d/%m/%Y')


pltClim2 <- ggplot(climate2, aes(x=Date)) + 
  geom_bar(aes(y=Rain/5),
           fill = "#99d8c9",
           #position= position_dodge(width=0.9), 
           stat="identity") + 
  geom_line(aes(y=Max),
            colour="#e34a33", 
            linetype="solid",
            show.legend=T,
            size = 2) +
  geom_line(aes(y=Min), 
            colour="#2b8cbe", 
            linetype="solid",
            size = 2) +
  scale_x_date(date_breaks="1 month",date_labels = "%m-%y") +
  scale_y_continuous(name = expression(Temperature~(degree*C)) , 
                     sec.axis = sec_axis(~.*5, name = "Rain (mm)")) + 
  theme(text = element_text(family = "Helvetica", 
                            color="black", 
                            size=12)) +
  theme_bw( ) 

pltClim2

png(file="MyerscoughClimate.png", units="in", width=5, height=3, res=300)
pltClim2
dev.off()
