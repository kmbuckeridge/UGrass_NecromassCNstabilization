setwd("/Users/katebuckeridge/Dropbox/R/UGrass/CNstab")


## note: all files changed when data transformed to atm% 2020 02 06
###### figs ####################################################

require(ggplot2)
library(ggpubr)
library(plyr)

######   summarySE function to define means and se for plots


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=F,
                      conf.interval=.95, .drop=TRUE) {
  
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=F) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

###### data for figs same as stats

field <- read.csv("field_depth.csv", header=TRUE, sep=",")


#################################################################################
### fine C (correct)

### first make table of means and se

Cf <- summarySE(data=field, measurevar="soilCf", groupvars=c("LUI","depth","time"), na.rm=F,
                conf.interval=.95)

### make new positions for stacked error bars
Cf$y_C = NA
Cf$y_C[Cf$depth == "subsurface"] = Cf$soilCf[Cf$depth == "subsurface"]
Cf$y_C[Cf$depth == "asurface"] = Cf$soilCf[Cf$depth == "asurface"] + 
  Cf$soilCf[Cf$depth == "subsurface"]

# define error bars
limits <- aes(ymax = y_C + se, ymin= y_C - se)
# make time a factor #may not be necessary if factor in field
Cf$time <- as.factor(Cf$time)

pltC <- ggplot(Cf, aes(fill=depth, y=soilCf, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Mineral-associated~C~from~added~necromass~(g~m^-2)),
                     limits=c(0,1.6)) +
  geom_errorbar(limits,
                position= "identity", 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltC

#### fine N 

Nf <- summarySE(data=field, measurevar="soilNf", groupvars=c("LUI","depth","time"), na.rm=F,
                conf.interval=.95)

#make new positions for stacked error bars
Nf$y_N = NA
Nf$y_N[Nf$depth == "subsurface"] = Nf$soilNf[Nf$depth == "subsurface"]
Nf$y_N[Nf$depth == "asurface"] = Nf$soilNf[Nf$depth == "asurface"] + 
  Nf$soilNf[Nf$depth == "subsurface"]

#define error bars
limitsNf <- aes(ymax = y_N + se, ymin= y_N - se)

#make time a factor
Nf$time <- as.factor(Nf$time)

pltN <- ggplot(Nf, aes(fill=depth, y=soilNf, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Mineral-associated~N~from~added~necromass~(g~m^-2)),
                     limits=c(0,1.6)) +
  geom_errorbar(limitsNf,
                position= "identity", 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) +  
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltN

######## fig

Fig4abc <- ggarrange(pltC, pltN,pltCN,
                   ncol = 2, nrow = 2, widths = c(1,1), common.legend = TRUE, align = "h")

png(file="Fig4abc.png", units="in", width=10, height=8, res=300)
Fig4abc
dev.off()

#### fine CN 

field$soilCNf <- (field$soilCf*12)/(field$soilNf*14)

CNf <- summarySE(data=field, measurevar="soilCNf", groupvars=c("LUI","depth","time"), na.rm=F,
                conf.interval=.95)

#define error bars
limitsCNf <- aes(ymax = soilCNf + se, ymin= soilCNf - se)
#make time a factor
CNf$time <- as.factor(CNf$time)


pltCN <- ggplot(CNf, aes(fill=depth, y=soilCNf, x=time)) +
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Carbon:nitrogen~ratio~from~added~necromass),
                     breaks=c(0,1,2,3),
                     labels=c(0,1,2,3),
                     limits=c(0,3)) +
  geom_errorbar(limitsCNf, 
                position= position_dodge(width=0.9), 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) +  
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 

pltCN


png(file="fieldCNratio.png", units="in", width=6, height=4, res=300)
pltCN
dev.off()





#######################################################################################

### bulk soil CN and C:N

### fine C (correct)

### first make table of means and se

Cb <- summarySE(data=field, measurevar="soilCb", groupvars=c("LUI","depth","time"), na.rm=F,
                conf.interval=.95)

### make new positions for stacked error bars
Cb$y_C = NA
Cb$y_C[Cb$depth == "subsurface"] = Cb$soilCb[Cb$depth == "subsurface"]
Cb$y_C[Cb$depth == "asurface"] = Cb$soilCb[Cb$depth == "asurface"] + 
  Cb$soilCb[Cb$depth == "subsurface"]

# define error bars
limits <- aes(ymax = y_C + se, ymin= y_C - se)
# make time a factor #may not be necessary if factor in field
Cb$time <- as.factor(Cb$time)

pltCb <- ggplot(Cb, aes(fill=depth, y=soilCb, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Bulk~soil~C~from~added~necromass~(g~m^-2)),
                     limits=c(0,2.5)) +
  geom_errorbar(limits,
                position= "identity", 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltCb

#### bulk N 

Nb <- summarySE(data=field, measurevar="soilNb", groupvars=c("LUI","depth","time"), na.rm=F,
                conf.interval=.95)

#make new positions for stacked error bars
Nb$y_N = NA
Nb$y_N[Nb$depth == "subsurface"] = Nb$soilNb[Nb$depth == "subsurface"]
Nb$y_N[Nb$depth == "asurface"] = Nb$soilNb[Nb$depth == "asurface"] + 
  Nb$soilNb[Nb$depth == "subsurface"]

#define error bars
limitsNb <- aes(ymax = y_N + se, ymin= y_N - se)

#make time a factor
Nb$time <- as.factor(Nb$time)

pltNb <- ggplot(Nb, aes(fill=depth, y=soilNb, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Bulk~soil~N~from~added~necromass~(g~m^-2)),
                     limits=c(0,2.5)) +
  geom_errorbar(limitsNb,
                position= "identity", 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) +  
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltNb


#### bulk C:N 

field$soilCNb <- (field$soilCb*12)/(field$soilNb*14)

CNb <- summarySE(data=field, measurevar="soilCNb", groupvars=c("LUI","depth","time"), na.rm=F,
                 conf.interval=.95)

#define error bars
limitsCNb <- aes(ymax = soilCNb + se, ymin= soilCNb - se)
#make time a factor
CNb$time <- as.factor(CNb$time)


pltCNb <- ggplot(CNb, aes(fill=depth, y=soilCNb, x=time)) +
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Carbon:nitrogen~ratio~from~added~necromass),
                     breaks=c(0,1,2,3),
                     labels=c(0,1,2,3),
                     limits=c(0,3)) +
  geom_errorbar(limitsCNb, 
                position= position_dodge(width=0.9), 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) +  
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 

pltCNb


png(file="fieldCNratioBulk.png", units="in", width=6, height=4, res=300)
pltCNb
dev.off()


######## fig for ms (supp)

Figsupp5abc <- ggarrange(pltCb, pltNb,pltCNb,
                     ncol = 2, nrow = 2, widths = c(1,1), common.legend = TRUE, align = "h")

png(file="Figsupp5abc.png", units="in", width=10, height=8, res=300)
Figsupp5abc
dev.off()





###################################################################################
#### MB figure  ###########

### mb13C

### first make table of means and se
mb13c <- summarySE(data=field, measurevar="MB13C", groupvars=c("LUI","depth","time"), na.rm=F,
                conf.interval=.95)

### make new positions for stacked error bars
mb13c$y_C = NA
mb13c$y_C[mb13c$depth == "subsurface"] = mb13c$MB13C[mb13c$depth == "subsurface"]
mb13c$y_C[mb13c$depth == "asurface"] = mb13c$MB13C[mb13c$depth == "asurface"] + 
  mb13c$MB13C[mb13c$depth == "subsurface"]

# define error bars
limitsmb13c <- aes(ymax = y_C + se, ymin= y_C - se)
#make time a factor
mb13c$time <- as.factor(mb13c$time)

pltmb13C <- ggplot(mb13c, aes(fill=depth, y=MB13C, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Microbial~biomass~C~from~added~necromass~(mg~m^-2)),
                     breaks=c(0,20,40,60,80),
                     labels=c(0,20,40,60,80),
                     limits=c(0,85)) +
  geom_errorbar(limitsmb13c,
                position= "identity", 
                width=0.6) +
   scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) + 
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltmb13C

# ### mb13C and mbc (dual axes) can't get it to work...
# 
# pltmbC2 <- ggplot(mapping = aes(x = fieldCN$time, y = fieldCN$mbC, fill = fieldCN$Depth)) +
#   geom_bar(position= "stack", 
#            stat="identity") + 
#   geom_point(x = fieldCN$time, y = fieldCN$bC*10, fill = fieldCN$Depth,
#              size = 2,
#              alpha = 0.5, 
#              position= "identity", 
#              stat="identity") +
#   scale_x_discrete(name="Time since necromass added (d)",
#                    labels=c("3","21","92","238")) +
#   scale_y_continuous(name=expression(Microbial~biomass~C~from~added~necromass~(mu*g~m^-2)),
#                      breaks=c(0,50,100,150,200,250),
#                      labels=c(0,50,100,150,200,250),
#                      limits=c(0,300),
#                      sec.axis = sec_axis(~./10, name =expression(Microbial~biomass~C~(g~m^-2) )) )+
#  # geom_errorbar(limitsmbC,
#  #               position= "identity", 
#  #               width=0.6) +
#   #geom_errorbar(limitsbC,
#      #           position = "identity",
#      #           width = 0.3,
#      #           color = gray) +
#   scale_fill_manual(name = "Depth",
#                     labels = c("Surface (0-5 cm)","Subsurface (5-10 cm)"),
#                     values = c("#f1a340","#998ec3")) + 
#   theme(text = element_text(family = "Helvetica", 
#                             color="black", 
#                             size=12),
#         axis.title.y.right = element_text (colour = gray) ) +
#   theme_bw( ) 
#   facet_grid(.~LUI)
# pltmbC2

### mbc 

### first make table of means and se
mbc <- summarySE(data=field, measurevar="mbc", groupvars=c("LUI","depth","time"), na.rm=F,
                   conf.interval=.95)

### make new positions for stacked error bars
mbc$y_C = NA
mbc$y_C[mbc$depth == "subsurface"] = mbc$mbc[mbc$depth == "subsurface"]
mbc$y_C[mbc$depth == "asurface"] = mbc$mbc[mbc$depth == "asurface"] + 
  mbc$mbc[mbc$depth == "subsurface"]

# define error bars
limitsmbc <- aes(ymax = y_C + se, ymin= y_C - se)
#make time a factor
mbc$time <- as.factor(mbc$time)

pltmbc <- ggplot(mbc, aes(fill=depth, y=mbc, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Microbial~biomass~C~(g~m^-2)),
                     breaks=c(0,25,50,75,100),
                     labels=c(0,25,50,75,100),
                     limits=c(0,120)) +
  geom_errorbar(limitsmbc,
                position= "identity", 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) + 
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltmbc


#### mb15N

### first make table of means and se
mb15n <- summarySE(data=field, measurevar="MB15N", groupvars=c("LUI","depth","time"), na.rm=F,
                   conf.interval=.95)

#make new positions for stacked error bars
mb15n$y_mbN = NA
mb15n$y_mbN[mb15n$depth == "subsurface"] = mb15n$MB15N[mb15n$depth == "subsurface"]
mb15n$y_mbN[mb15n$depth == "asurface"] = mb15n$MB15N[mb15n$depth == "asurface"] + 
  mb15n$MB15N[mb15n$depth == "subsurface"]

#define error bars
limitsmb15n <- aes(ymax = y_mbN + se, ymin= y_mbN - se)
#make time a factor
mb15n$time <- as.factor(mb15n$time)

pltmb15n <- ggplot(mb15n, aes(fill=depth, y=MB15N, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Microbial~biomass~N~from~added~necromass~(mg~m^-2)),
                     breaks=c(0,10,20,30,40),
                     labels=c(0,10,20,30,40),
                     limits=c(0,50)) +
  geom_errorbar(limitsmb15n,
                position= "identity", 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) + 
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltmb15n

### mbn 

### first make table of means and se
mbn <- summarySE(data=field, measurevar="mbn", groupvars=c("LUI","depth","time"), na.rm=F,
                 conf.interval=.95)

### make new positions for stacked error bars
mbn$y_n = NA
mbn$y_n[mbn$depth == "subsurface"] = mbn$mbn[mbn$depth == "subsurface"]
mbn$y_n[mbn$depth == "asurface"] = mbn$mbn[mbn$depth == "asurface"] + 
  mbn$mbn[mbn$depth == "subsurface"]

# define error bars
limitsmbn <- aes(ymax = y_n + se, ymin= y_n - se)
#make time a factor
mbn$time <- as.factor(mbn$time)

pltmbn <- ggplot(mbn, aes(fill=depth, y=mbn, x=time)) +
  geom_bar(position= "stack", 
           stat="identity", 
           aes(fill = depth)) + 
  scale_x_discrete(name="Time since necromass added (d)",
                   labels=c("3","21","92","238")) +
  scale_y_continuous(name=expression(Microbial~biomass~N~(g~m^-2)),
                     breaks=c(0,2.5,5,7.5,10),
                     labels=c(0,2.5,5,7.5,10),
                     limits=c(0,12)) +
  geom_errorbar(limitsmbn,
                position= "identity", 
                width=0.6) +
  scale_fill_manual(name = "Depth",
                    labels = c("0-5 cm","5-10 cm"),
                    values = c("#2c7fb8","#7fcdbb")) + 
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
pltmbn

#### ms figures #######


Fig3ab <- ggarrange(co2sub, n2osub,
                   ncol = 2, nrow = 1, widths = c(1,1), common.legend = TRUE, align = "h")
Fig3cd <- ggarrange(pltmb13C, pltmb15n,
                    ncol = 2, nrow = 1, widths = c(1,1), common.legend = TRUE, align = "h")

png(file="Fig3_gas.png", units="in", width=10, height=4, res=300)
Fig3ab
dev.off()
png(file="Fig3_mb.png", units="in", width=10, height=4, res=300)
Fig3cd
dev.off()

SuppFig3_gas <- ggarrange(co2, n2o,
                  ncol = 2, nrow = 1, widths = c(1,1), common.legend = TRUE, align = "h")

SuppFig3_mb <- ggarrange(pltmbc, pltmbn,
                          ncol = 2, nrow = 1, widths = c(1,1), common.legend = TRUE, align = "h")

png(file="Fig3Supp_gas.png", units="in", width=10, height=4, res=300)
SuppFig3_gas
dev.off()
png(file="Fig3Supp_mb.png", units="in", width=10, height=4, res=300)
SuppFig3_mb
dev.off()

######### deep soils  ############
######## fine fraction

deepFB <- read.csv("field_deep2.csv", header=TRUE, sep=",")
deepMB <- read.csv("field_deep.csv", header=TRUE, sep=",")

dfb <- summarySE(data=deepFB, measurevar="sub", groupvars=c("LUI","Pool","Element"), na.rm=F,
                 conf.interval=.95)
limits.dfb <- aes(ymax = sub + se, ymin= sub - se)
dFB <- ggplot(dfb, aes(x=Element, y=sub, fill=Pool)) +
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = Pool)) + 
  geom_errorbar(limits.dfb, 
                position= position_dodge(width=0.9), 
                width=0.8) +
  scale_x_discrete(name="",
                   labels=c("Carbon","Nitrogen")) +
  scale_y_continuous(name=expression(Recovery~of~added~necromass~(g~m^-2))) +
                    # limits=c(0,12),
                     #breaks = c(0,2,4,6,8,10)) +
  
  scale_fill_manual(values = c("#b2df8a","#a6cee3"), 
                    name = "Pool",
                    labels = c("Bulk soil","Mineral-associated")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
dFB

dmb <- summarySE(data=deepMB, measurevar="sub", groupvars=c("LUI","Pool","Element"), na.rm=F,
                 conf.interval=.95)
limits.dmb <- aes(ymax = sub + se, ymin= sub - se)
dMB <- ggplot(dmb, aes(x=Element, y=sub, fill=Pool)) +
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = Pool)) + 
  geom_errorbar(limits.dmb, 
                position= position_dodge(width=0.9), 
                width=0.8) +
  scale_x_discrete(name="",
                   labels=c("Carbon","Nitrogen")) +
  scale_y_continuous(name=expression(Recovery~of~added~necromass~(mg~m^-2))) +
  # limits=c(0,12),
  #breaks = c(0,2,4,6,8,10)) +
  
  scale_fill_manual(values = c("#1f78b4"), 
                    name = "Pool",
                    labels = c("Microbial biomass")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
dMB

FigS4ab <- ggarrange(dFB, dMB,
                    ncol = 2, nrow = 1, widths = c(1,0.75), align = "h")

png(file="deep.png", units="mm", width=267, height=100, res=300)
FigS4ab
dev.off()




#############################################################################
### GHG

# NOTE: excel screws up date format and reverses month and day - double-check csv before import

ghg <- read.csv("field_GHG.csv", header=TRUE, sep=",")

ghg$Date <- as.Date(ghg$Date,format='%d/%m/%Y')
ghg$sublui <- paste(ghg$substrate, ghg$LUI, sep="_") #this if I want to plot by substrate and LUI - aes can only handle one variable, so need to make a concatenated column for the two variables - then change the grouping and colour variable to 'sublui' - actually removed this from .csv fle, now just necro plots

ghgn <- summarySE(data=ghg, measurevar="mgN2ONd", groupvars=c("LUI","Date"), na.rm=F,
                   conf.interval=.95)
ghgn$Date <- as.Date(ghgn$Date,format='%d/%m/%Y')
ghgns <- summarySE(data=ghg, measurevar="cumN2Osub", groupvars=c("LUI","Date"), na.rm=F,
                   conf.interval=.95)
n2o_sub <- summarySE(data=ghg, measurevar="ngN2Osub", groupvars=c("LUI","Date"), na.rm=F,
                   conf.interval=.95)
n2o_sub$Date <- as.Date(n2o_sub$Date,format='%d/%m/%Y')

ghgc <- summarySE(data=ghg, measurevar="gCO2Cd", groupvars=c("LUI","Date"), na.rm=F,
                  conf.interval=.95)
ghgc

co2_sub <- summarySE(data=ghg, measurevar="mgCO2sub", groupvars=c("LUI","Date"), na.rm=F,
                  conf.interval=.95)
co2_sub$Date <- as.Date(co2_sub$Date,format='%d/%m/%Y')

ghgc$Date <- as.Date(ghgc$Date,format='%d/%m/%Y')
ghgcs <- summarySE(data=ghg, measurevar="cumCO2sub", groupvars=c("LUI","Date"), na.rm=F,
                   conf.interval=.95)
ghgcs
ghgcs$Date <- as.Date(ghgcs$Date,format='%d/%m/%Y')

# n2o
n2o <- ggplot(ghgn, aes(x=Date, y=mgN2ONd, group=LUI, colour=LUI, shape = LUI)) + 
     geom_errorbar(aes(ymin=mgN2ONd-se, ymax=mgN2ONd+se), 
                   width=.1, 
                   position=position_dodge(0.05)) +
    # geom_line(size=0.5) + 
     geom_point(aes(shape = factor(LUI)), size=3)+
     scale_x_date(date_breaks="1 month",date_labels = "%m-%y") +
     scale_y_continuous(name=expression(N[2]*O~(mg~N~m^-2~d^-1))) +
     scale_color_manual(values = c("black", "dark grey"), 
                        labels = c("High intensity", "Low intensity"),
                        name = "Land management") +
     scale_shape_manual(values = c(19,17), 
                        labels = c("High intensity", "Low intensity"),
                        name = "Land management") +
  theme_bw()+
  theme(axis.title.x = element_blank())
n2o

# n2o from substrate
n2osub<- ggplot(n2o_sub, aes(x=Date, y=ngN2Osub, group=LUI, colour=LUI, shape = LUI)) + 
  geom_errorbar(aes(ymin=ngN2Osub-se, ymax=ngN2Osub+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  #geom_line(size=0.5) + 
  geom_point(aes(shape = factor(LUI)), size=3)+
  scale_x_date(date_breaks="1 month",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(N[2]*O~from~necromass~(ng~N~m^-2))) +
  scale_color_manual(values = c("black", "dark grey"), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  scale_shape_manual(values = c(19,17), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  theme_bw()+
  theme(axis.title.x = element_blank())
n2osub


# n2o cumulative from substrate
n2os<- ggplot(ghgns, aes(x=Date, y=cumN2Osub, group=LUI, colour=LUI, shape = LUI)) + 
  geom_errorbar(aes(ymin=cumN2Osub-se, ymax=cumN2Osub+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  geom_line(size=0.5) + 
  geom_point(aes(shape = factor(LUI)), size=3)+
  scale_x_date(date_breaks="1 month",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(N[2]*O~from~necromass~(mu*g~N~m^-2))) +
  scale_color_manual(values = c("black", "dark grey"), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  scale_shape_manual(values = c(19,17), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  theme_bw()+
  theme(axis.title.x = element_blank())
n2os

# co2
co2<- ggplot(ghgc, aes(x=Date, y=gCO2Cd, group=LUI, colour=LUI, shape = LUI)) + 
    geom_errorbar(aes(ymin=gCO2Cd-se, ymax=gCO2Cd+se), 
                width=.1, 
                position=position_dodge(0.05)) +
   # geom_line(size=0.5) + 
    geom_point(aes(shape = factor(LUI)), size=3)+
    scale_x_date(date_breaks="1 month",date_labels = "%m-%y") +
    scale_y_continuous(name=expression(CO[2]~(g~C~m^-2~d^-1))) +
  scale_color_manual(values = c("black", "dark grey"), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  scale_shape_manual(values = c(19,17), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
    theme_bw()+
  theme(axis.title.x = element_blank())
co2

# co2 from substrate
co2sub<- ggplot(co2_sub, aes(x=Date, y=mgCO2sub, group=LUI, colour=LUI, shape = LUI)) + 
  geom_errorbar(aes(ymin=mgCO2sub-se, ymax=mgCO2sub+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  #geom_line(size=0.5) + 
  geom_point(aes(shape = factor(LUI)), size=3)+
  scale_x_date(date_breaks="1 month",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(CO[2]~from~necromass~(mg~C~m^-2))) +
  scale_color_manual(values = c("black", "dark grey"), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  scale_shape_manual(values = c(19,17), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  theme_bw()+
  theme(axis.title.x = element_blank())
co2sub

# co2 cumulative from substrate
co2s<- ggplot(ghgcs, aes(x=Date, y=cumCO2sub, group=LUI, colour=LUI, shape = LUI)) + 
  geom_errorbar(aes(ymin=cumCO2sub-se, ymax=cumCO2sub+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  geom_line(size=0.5) + 
  geom_point(aes(shape = factor(LUI)), size=3)+
  scale_x_date(date_breaks="1 month",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(CO[2]~from~necromass~(g~C~m^-2))) +
  scale_color_manual(values = c("black", "dark grey"), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  scale_shape_manual(values = c(19,17), 
                     labels = c("High intensity", "Low intensity"),
                     name = "Land management") +
  theme_bw()+
  theme(axis.title.x = element_blank())
co2s


Fig.ghg <- ggarrange(co2,n2o,co2s,n2os,
                     ncol = 2, nrow = 2, widths = c(1,1), common.legend = TRUE, align = "h")


png(file="fieldGHG.png", units="in", width=10, height=8, res=300)
Fig.ghg
dev.off()

