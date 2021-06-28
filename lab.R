setwd("/Users/katebuckeridge/Dropbox/R/UGrass/CNstab")


lab <- read.csv("lab.csv", header=TRUE, sep=",")
library(emmeans) #was lsmeans

##### stats

labC <- lab[lab$Element == "Carbon", ]
labN <- lab[lab$Element == "Nitrogen", ]

##### Bulk carbon immobilization of substrate higher with roots than leaves or microbial residue, no LUI diff

model.immob.1 <- lm(Bulk~LUI*Necro, data=labC)

opar = par (mfrow = c(2,2))
plot(model.immob.1)
par(opar)

anova(model.immob.1)
# Response: Bulk
#           Df  Sum Sq Mean Sq F value    Pr(>F)    
# LUI        1   93.27   93.27  1.1770 0.2887483    
# Necro      2 1707.96  853.98 10.7767 0.0004574 ***
# LUI:Necro  2  435.47  217.74  2.7477 0.0842322 .  
# Residuals 24 1901.84   79.24  

labC$bulk_log <- log(labC$Bulk)
model.immob.5 <- lm(bulk_log~LUI*Necro, data=labC)

opar = par (mfrow = c(2,2))
plot(model.immob.5)
par(opar)

anova(model.immob.5)
# Response: bulk_log
#           Df  Sum Sq Mean Sq F value   Pr(>F)    
# LUI        1 0.02178 0.02178  0.3717 0.547806    
# Necro      2 1.15955 0.57978  9.8964 0.000734 ***
# LUI:Necro  2 0.27897 0.13949  2.3809 0.113955    
# Residuals 24 1.40603 0.05858 

lsmeans(model.immob.5, pairwise ~ Necro)
# $lsmeans
# Necro     lsmean     SE df lower.CL upper.CL
# Leaves      3.47 0.0765 24     3.32     3.63
# Microbial   3.37 0.0765 24     3.21     3.53
# Roots       3.83 0.0765 24     3.67     3.99
# 
# Results are averaged over the levels of: LUI 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast           estimate    SE df t.ratio p.value
# Leaves - Microbial    0.104 0.108 24  0.957  0.6107 
# Leaves - Roots       -0.356 0.108 24 -3.285  0.0085 
# Microbial - Roots    -0.459 0.108 24 -4.241  0.0008 
# 
# Results are averaged over the levels of: LUI 
# P value adjustment: tukey method for comparing a family of 3 estimates  


####### Bulk N

model.immob.11 <- lm(Bulk~LUI*Necro, data=labN)

opar = par (mfrow = c(2,2))
plot(model.immob.11)
par(opar)

anova(model.immob.11)
# Response: Bulk
# Df Sum Sq Mean Sq F value  Pr(>F)  
# LUI        1   47.5   47.54  0.1889 0.66772  
# Necro      2  492.7  246.36  0.9788 0.39026  
# LUI:Necro  2 1325.0  662.52  2.6323 0.09256 .
# Residuals 24 6040.5  251.69 

labN$bulk_log <- log(labN$Bulk)
model.immob.15 <- lm(bulk_log~LUI*Necro, data=labN)

opar = par (mfrow = c(2,2))
plot(model.immob.15)
par(opar)

anova(model.immob.15)
# Response: bulk_log
# Df  Sum Sq  Mean Sq F value Pr(>F)
# LUI        1 0.00011 0.000115  0.0032 0.9554
# Necro      2 0.04984 0.024919  0.6953 0.5087
# LUI:Necro  2 0.17357 0.086785  2.4215 0.1102
# Residuals 24 0.86016 0.035840  



##### Fine mineral C immobilization higher in low LUI, Microbial residues retention trend higher than roots

model.immob.2 <- lm(Fine~LUI*Necro, data=labC)

opar = par (mfrow = c(2,2))
plot(model.immob.2)
par(opar)

anova(model.immob.2)
# Response: Fine
#           Df Sum Sq Mean Sq F value   Pr(>F)   
# LUI        1 278.43 278.432  8.8917 0.006478 **
# Necro      2 196.12  98.058  3.1315 0.061887 . 
# LUI:Necro  2 104.49  52.244  1.6684 0.209681   
# Residuals 24 751.53  31.314

labC$fine_log <- log(labC$Fine)
model.immob.3 <- lm(fine_log~LUI*Necro, data=labC)

opar = par (mfrow = c(2,2))
plot(model.immob.3)
par(opar)

anova(model.immob.3)
# Response: fine_log
#           Df  Sum Sq Mean Sq F value   Pr(>F)   
# LUI        1 0.43242 0.43242  9.2433 0.005638 **
# Necro      2 0.28425 0.14213  3.0380 0.066662 . 
# LUI:Necro  2 0.11165 0.05583  1.1933 0.320566   
# Residuals 24 1.12277 0.04678 

with(labC, tapply(Fine, list(LUI), FUN=mean, na.rm=TRUE))
# High intensity  Low intensity 
# 24.98088       31.07385 
with(labC, tapply(Fine, list(Necro), FUN=mean, na.rm=TRUE))
# Leaves Microbial     Roots 
# 27.56910  31.36267  25.15032 

lsmeans(model.immob.3, pairwise ~ LUI)
# $lsmeans
# LUI            lsmean     SE df lower.CL upper.CL
# High intensity   3.18 0.0558 24     3.07     3.30
# Low intensity    3.42 0.0558 24     3.31     3.54
# 
# Results are averaged over the levels of: Necro 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                       estimate    SE df t.ratio p.value
# High intensity - Low intensity    -0.24 0.079 24 -3.040  0.0056 
# 
# Results are averaged over the levels of: Necro 

lsmeans(model.immob.3, pairwise ~ Necro)
# $lsmeans
# Necro     lsmean     SE df lower.CL upper.CL
# Leaves      3.29 0.0684 24     3.15     3.43
# Microbial   3.43 0.0684 24     3.29     3.57
# Roots       3.19 0.0684 24     3.05     3.33
# 
# Results are averaged over the levels of: LUI 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast           estimate     SE df t.ratio p.value
# Leaves - Microbial   -0.137 0.0967 24 -1.411  0.3511 
# Leaves - Roots        0.101 0.0967 24  1.045  0.5568 
# Microbial - Roots     0.238 0.0967 24  2.456  0.0545
# 
# Results are averaged over the levels of: LUI 
# P value adjustment: tukey method for comparing a family of 3 estimates 

######### Fine N 

model.immob.12 <- lm(Fine~LUI*Necro, data=labN)

opar = par (mfrow = c(2,2))
plot(model.immob.12)
par(opar)

anova(model.immob.12)
#           Df Sum Sq Mean Sq F value    Pr(>F)    
# LUI        1 6105.2  6105.2 29.0272 1.562e-05 ***
# Necro      2 1552.9   776.5  3.6918    0.0400 *  
# LUI:Necro  2   24.7    12.3  0.0586    0.9432    
# Residuals 24 5047.8   210.3   

labN$fine_log <- log(labN$Fine)
model.immob.13 <- lm(fine_log~LUI*Necro, data=labN)

opar = par (mfrow = c(2,2))
plot(model.immob.13)
par(opar)

anova(model.immob.13)
# Response: fine_log
#           Df  Sum Sq Mean Sq F value    Pr(>F)    
# LUI        1 1.06028 1.06028 33.3548 5.934e-06 ***
# Necro      2 0.24008 0.12004  3.7764   0.03751 *  
# LUI:Necro  2 0.01480 0.00740  0.2329   0.79404    
# Residuals 24 0.76291 0.03179 

lsmeans(model.immob.13, pairwise ~ LUI)
# $lsmeans
# LUI            lsmean     SE df lower.CL upper.CL
#  High intensity   4.13 0.046 24     4.04     4.23
# Low intensity    4.51 0.046 24     4.41     4.60
# 
# Results are averaged over the levels of: Necro 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                       estimate     SE df t.ratio p.value
# High intensity - Low intensity   -0.376 0.0651 24 -5.775  <.0001 

# 
# Results are averaged over the levels of: Necro 

lsmeans(model.immob.13, pairwise ~ Necro)
# $lsmeans
# Necro     lsmean     SE df lower.CL upper.CL
# Leaves      4.23 0.0564 24     4.11     4.34
# Microbial   4.30 0.0564 24     4.18     4.41
# Roots       4.44 0.0564 24     4.32     4.56
# 
# Results are averaged over the levels of: LUI 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast           estimate     SE df t.ratio p.value
# Leaves - Microbial  -0.0723 0.0797 24 -0.907  0.6411 
# Leaves - Roots      -0.2153 0.0797 24 -2.700  0.0323 
# Microbial - Roots   -0.1430 0.0797 24 -1.793  0.1933 
# 
# Results are averaged over the levels of: LUI 

################

micro <- labC[11:20,]
model.immob.3 <- lm(fine_log~LUI, data=micro)
anova(model.immob.3)
# Response: fine_log
# Df  Sum Sq  Mean Sq F value   Pr(>F)   
# LUI        1 0.26138 0.261377  16.681 0.003513 **
# Residuals  8 0.12535 0.015669  
leaf <- labC[1:10,]
model.immob.3 <- lm(fine_log~LUI, data=leaf)
anova(model.immob.3)
#Response: fine_log
# Df  Sum Sq  Mean Sq F value Pr(>F)
# LUI        1 0.00105 0.001053  0.0146 0.9067
# Residuals  8 0.57552 0.071941 
root <- labC[21:30,]
model.immob.3 <- lm(fine_log~LUI, data=root)
anova(model.immob.3)
# Response: fine_log
# Df  Sum Sq Mean Sq F value  Pr(>F)  
# LUI        1 0.21738 0.21738  3.7904 0.08741 .
# Residuals  8 0.45880 0.05735 



########### GHG stats

lghg <- read.csv("lab_GHG.csv", header=TRUE, sep=",")

lghg2 <- lghg[lghg$Treatment != "W",] #remove NAs for cumulative CO2 and cum N2O

library(nlme)
library(gridExtra)
library(emmeans)

m.co2 = lme(cum13CO2 ~ LUI*Treatment,
            random = ~ 1 + Day|Plot,
            data = lghg2,
            method = "REML")

grid.arrange(plot(m.co2,type=c("p","smooth")),
             plot(m.co2,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.co2,resid(.,type="pearson")~Day,
                  type=c("p","smooth")),
             qqnorm(m.co2,abline=c(0,1)))
# residual vs fitted (non-linear patterns) should: bounce randomly around 0, form horizontal, no outliers
# scale-location (homoscedasticity = equal variance) should see: equal spread along line, line horizontal
# Residuals vs time looks for influential outliers
# normal Q-Q (normal distribution) should: be on the line or there is skew

summary(m.co2)
# Linear mixed-effects model fit by REML
# Data: lghg2 
# AIC      BIC    logLik
# 2893.618 2931.426 -1436.809
# 
# Random effects:
#   Formula: ~1 + Day | Plot
# Structure: General positive-definite, Log-Cholesky parametrization
#              StdDev     Corr  
# (Intercept)  0.2801321 (Intr)
# Day          0.1555713 -0.618
# Residual    18.4500008       
# 
# Fixed effects: cum13CO2 ~ LUI * Treatment 
#                                 Value Std.Error  DF   t-value p-value
# (Intercept)                  51.66852  3.049569 316 16.942893  0.0000
# LUILow intensity            -10.21430  4.312742   8 -2.368401  0.0454
# TreatmentN                  -17.53306  3.518277 316 -4.983421  0.0000
# TreatmentR                  -30.18540  3.518277 316 -8.579596  0.0000
# LUILow intensity:TreatmentN  15.26953  4.975595 316  3.068886  0.0023
# LUILow intensity:TreatmentR  11.03800  4.975595 316  2.218429  0.0272
# Correlation: 
#   (Intr) LUILwi TrtmnN TrtmnR LUILi:TN
# LUILow intensity            -0.707                              
# TreatmentN                  -0.577  0.408                       
# TreatmentR                  -0.577  0.408  0.500                
# LUILow intensity:TreatmentN  0.408 -0.577 -0.707 -0.354         
# LUILow intensity:TreatmentR  0.408 -0.577 -0.354 -0.707  0.500  
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.6683852 -0.8225108  0.1750345  0.6882192  1.9364194 
# 
# Number of Observations: 330
# Number of Groups: 10 


emmeans(m.co2, pairwise ~ LUI : Treatment)
# $emmeans
# LUI            Treatment emmean   SE df lower.CL upper.CL
# High intensity L           51.7 3.05  9     44.8     58.6
# Low intensity  L           41.5 3.05  8     34.4     48.5
# High intensity N           34.1 3.05  9     27.2     41.0
# Low intensity  N           39.2 3.05  8     32.2     46.2
# High intensity R           21.5 3.05  9     14.6     28.4
# Low intensity  R           22.3 3.05  8     15.3     29.3
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                            estimate   SE  df t.ratio p.value
# High intensity L - Low intensity L    10.214 4.31   8  2.368  0.2700 
# High intensity L - High intensity N   17.533 3.52 316  4.983  <.0001 
# High intensity L - Low intensity N    12.478 4.31   8  2.893  0.1375 
# High intensity L - High intensity R   30.185 3.52 316  8.580  <.0001 
# High intensity L - Low intensity R    29.362 4.31   8  6.808  0.0013 
# Low intensity L - High intensity N     7.319 4.31   8  1.697  0.5678 
# Low intensity L - Low intensity N      2.264 3.52 316  0.643  0.9876 
# Low intensity L - High intensity R    19.971 4.31   8  4.631  0.0143 
# Low intensity L - Low intensity R     19.147 3.52 316  5.442  <.0001 
# High intensity N - Low intensity N    -5.055 4.31   8 -1.172  0.8383 
# High intensity N - High intensity R   12.652 3.52 316  3.596  0.0050 
# High intensity N - Low intensity R    11.829 4.31   8  2.743  0.1676 
# Low intensity N - High intensity R    17.708 4.31   8  4.106  0.0277 
# Low intensity N - Low intensity R     16.884 3.52 316  4.799  <.0001 
# High intensity R - Low intensity R    -0.824 4.31   8 -0.191  0.9999 
# 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 6 estimates 


m.n2o = lme(cum15N2O ~ LUI*Treatment,
            random = ~ 1 + Day|Plot,
            data = lghg2,
            method = "REML")

grid.arrange(plot(m.n2o,type=c("p","smooth")),
             plot(m.n2o,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.n2o,resid(.,type="pearson")~Day,
                  type=c("p","smooth")),
             qqnorm(m.n2o,abline=c(0,1)))

summary(m.n2o)
# Linear mixed-effects model fit by REML
# Data: lghg2 
# AIC      BIC    logLik
# 1097.418 1135.225 -538.7089
# 
# Random effects:
#   Formula: ~1 + Day | Plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev      Corr  
# (Intercept) 0.887036417 (Intr)
# Day         0.004312339 0.972 
# Residual    1.160936316       
# 
# Fixed effects: cum15N2O ~ LUI * Treatment 
#                                 Value Std.Error  DF    t-value p-value
# (Intercept)                  3.538989 0.3141983 316  11.263555  0.0000
# LUILow intensity            -1.616827 0.4443434   8  -3.638688  0.0066
# TreatmentN                  -2.192706 0.2213819 316  -9.904631  0.0000
# TreatmentR                  -4.568995 0.2213819 316 -20.638523  0.0000
# LUILow intensity:TreatmentN  0.372666 0.3130812 316   1.190318  0.2348
# LUILow intensity:TreatmentR  2.095841 0.3130812 316   6.694241  0.0000
# Correlation: 
#   (Intr) LUILwi TrtmnN TrtmnR LUILi:TN
# LUILow intensity            -0.707                              
# TreatmentN                  -0.352  0.249                       
# TreatmentR                  -0.352  0.249  0.500                
# LUILow intensity:TreatmentN  0.249 -0.352 -0.707 -0.354         
# LUILow intensity:TreatmentR  0.249 -0.352 -0.354 -0.707  0.500  
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.32655762 -0.47951104  0.01115104  0.40399649  2.79345622 
# 
# Number of Observations: 330
# Number of Groups: 10 

emmeans(m.n2o, pairwise ~ LUI : Treatment)
# $emmeans
# LUI            Treatment emmean    SE df lower.CL upper.CL
# High intensity L          3.539 0.314  9    2.828    4.250
# Low intensity  L          1.922 0.314  8    1.198    2.647
# High intensity N          1.346 0.314  9    0.636    2.057
# Low intensity  N          0.102 0.314  8   -0.622    0.827
# High intensity R         -1.030 0.314  9   -1.741   -0.319
# Low intensity  R         -0.551 0.314  8   -1.276    0.174
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                            estimate    SE  df t.ratio p.value
# High intensity L - Low intensity L     1.617 0.444   8  3.639  0.0510 
# High intensity L - High intensity N    2.193 0.221 316  9.905  <.0001 
# High intensity L - Low intensity N     3.437 0.444   8  7.735  0.0005 
# High intensity L - High intensity R    4.569 0.221 316 20.639  <.0001 
# High intensity L - Low intensity R     4.090 0.444   8  9.205  0.0001 
# Low intensity L - High intensity N     0.576 0.444   8  1.296  0.7803 
# Low intensity L - Low intensity N      1.820 0.221 316  8.221  <.0001 
# Low intensity L - High intensity R     2.952 0.444   8  6.644  0.0015 
# Low intensity L - Low intensity R      2.473 0.221 316 11.171  <.0001 
# High intensity N - Low intensity N     1.244 0.444   8  2.800  0.1555 
# High intensity N - High intensity R    2.376 0.221 316 10.734  <.0001 
# High intensity N - Low intensity R     1.897 0.444   8  4.270  0.0224 
# Low intensity N - High intensity R     1.132 0.444   8  2.548  0.2155 
# Low intensity N - Low intensity R      0.653 0.221 316  2.950  0.0396 
# High intensity R - Low intensity R    -0.479 0.444   8 -1.078  0.8773 
# 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 6 estimates 

##################################################################
###### figs


require(ggplot2)
library(ggpubr)
library(plyr)

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



######## fine fraction

labCN <- read.csv("lab.csv", header=TRUE, sep=",")

lcn <- summarySE(data=labCN, measurevar="Fine", groupvars=c("LUI","Necro","Element"), na.rm=F,
                 conf.interval=.95)
limits.lcn <- aes(ymax = Fine + se, ymin=Fine - se)

CN <- ggplot(lcn, aes(x=Element, y=Fine, fill=Necro)) +
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = Necro)) + 
  geom_errorbar(limits.lcn, 
                position= position_dodge(width=0.9), 
                width=0.8) +
  scale_x_discrete(name="",
                   labels=c("Carbon","Nitrogen")) +
  scale_y_continuous(name=expression(Mineral-associated~({"%"}~recovered)),
                     limits=c(0,120),
                     breaks = c(0,20,40,60,80,100)) +
  scale_fill_manual(values = c("#b2df8a","#a6cee3", "#1f78b4"),
                    name = "Substrate",
                    labels = c("Leaf","Microbial necromass","Root","Control")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
CN

blcn <- summarySE(data=labCN, measurevar="Bulk", groupvars=c("LUI","Necro","Element"), na.rm=F,
                 conf.interval=.95)
limits.blcn <- aes(ymax = Bulk + se, ymin=Bulk - se)

bulkCN <- ggplot(blcn, aes(x=Element, y=Bulk, fill=Necro)) +
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = Necro)) + 
  geom_errorbar(limits.blcn, 
                position= position_dodge(width=0.9), 
                width=0.8) +
  scale_x_discrete(name="",
                   labels=c("Carbon","Nitrogen")) +
  scale_y_continuous(name=expression(Bulk~soil~({"%"}~recovered)),
                     limits=c(0,120),
                     breaks = c(0,20,40,60,80,100)) +
  scale_fill_manual(values = c("#b2df8a","#a6cee3", "#1f78b4"),
                    name = "Substrate",
                    labels = c("Leaf","Microbial necromass","Root","Control")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
bulkCN

png(file="labCNfine.png", units="in", width=5, height=4, res=300)
CN
dev.off()




#############################################################################
### GHG

# NOTE: excel screws up date format and reverses month and day - double-check csv before import


# web example for lines from discrete x axis:
#   ggplot(hist, aes(x=weekday, y=counts, group=1)) +
#   geom_point(stat='summary', fun.y=sum) +
#   stat_summary(fun.y=sum, geom="line")

# or from dates
# Standard deviation of the mean
# ggplot(df3, aes(x=dose, y=len, group=supp, color=supp)) + 
#   geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1) +
#   geom_line() + geom_point()+
#   scale_color_brewer(palette="Paired")+theme_minimal()
# # Use position_dodge to move overlapped errorbars horizontally
# ggplot(df3, aes(x=dose, y=len, group=supp, color=supp)) + 
#   geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1, 
#                 position=position_dodge(0.05)) +
#   geom_line() + geom_point()+
#   scale_color_brewer(palette="Paired")+theme_minimal()

lghg <- read.csv("lab_GHG.csv", header=TRUE, sep=",")

lghg$Date <- as.Date(lghg$Date,format='%d/%m/%Y')

#lghg$sublui <- paste(lghg$Treatment, lghg$LUI, sep="_") #this if I want to plot by substrate and LUI - aes can only handle one variable, so need to make a concatenated column for the two variables - then change the grouping and colour variable to 'sublui' - option is to facet wrap

lghgn <- summarySE(data=lghg, measurevar="N2O", groupvars=c("LUI","Treatment","Date"), na.rm=F,
                  conf.interval=.95)

lghg.2 <- lghg[lghg$Treatment != "W",]
lghgns <- summarySE(data=lghg.2, measurevar="cum15N2O", groupvars=c("LUI","Treatment","Date"), na.rm=F,
                   conf.interval=.95)

lghgc <- summarySE(data=lghg, measurevar="CO2", groupvars=c("LUI","Treatment","Date"), na.rm=F,
                  conf.interval=.95)

lghgcs <- summarySE(data=lghg.2, measurevar="cum13CO2", groupvars=c("LUI","Treatment","Date"), na.rm=F,
                   conf.interval=.95)



l.n2o <- ggplot(lghgn, aes(x=Date, y=N2O, group=Treatment, colour=Treatment, shape = Treatment)) + 
  geom_errorbar(aes(ymin=N2O-se, ymax=N2O+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  geom_line(size=0.5) + 
  geom_point(aes(shape = factor(Treatment)), size=3, alpha = 0.7)+
  scale_x_date(date_breaks="3 months",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(N[2]*O~(ng~N~g^-1~d^-1))) +
  scale_color_manual(values = c("#b2df8a","#a6cee3", "#1f78b4","black"),
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  scale_shape_manual(values = c(15,19,17, 20), 
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
l.n2o


ln2os<- ggplot(lghgns, aes(x=Date, y=cum15N2O, group=Treatment, colour=Treatment, shape = Treatment)) + 
  geom_errorbar(aes(ymin=cum15N2O-se, ymax=cum15N2O+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  geom_line(size=0.5) + 
  geom_point(aes(shape = factor(Treatment)), size=3, alpha = 0.7)+
  scale_x_date(date_breaks="3 months",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(N[2]*O-N~from~substrate~(cumulative~{"%"}~recovered))) +
  scale_color_manual(values = c("#b2df8a","#a6cee3", "#1f78b4","black"),
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  scale_shape_manual(values = c(15,19,17,20), 
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
ln2os


lco2<- ggplot(lghgc, aes(x=Date, y=CO2, group=Treatment, colour=Treatment, shape = Treatment)) + 
  geom_errorbar(aes(ymin=CO2-se, ymax=CO2+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  geom_line(size=0.5) + 
  geom_point(aes(shape = factor(Treatment)), size=3, alpha = 0.7)+
  scale_x_date(date_breaks="3 months",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(CO[2]~(mu*g~C~g^-1~d^-1))) +
  scale_color_manual(values = c("#b2df8a","#a6cee3", "#1f78b4","black"),
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  scale_shape_manual(values = c(15,19,17,20), 
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
lco2

lco2s<- ggplot(lghgcs, aes(x=Date, y=cum13CO2, group=Treatment, colour=Treatment, shape = Treatment)) + 
  geom_errorbar(aes(ymin=cum13CO2-se, ymax=cum13CO2+se), 
                width=.1, 
                position=position_dodge(0.05)) +
  geom_line(size=0.5) + 
  geom_point(aes(shape = factor(Treatment)), size=3, alpha = 0.7)+
  scale_x_date(date_breaks="3 months",date_labels = "%m-%y") +
  scale_y_continuous(name=expression(CO[2]-C~from~substrate~(cumulative~{"%"}~recovered))) +
  scale_color_manual(values = c("#b2df8a","#a6cee3", "#1f78b4","black"),
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  scale_shape_manual(values = c(15,19,17,20), 
                     name = "Substrate",
                     labels = c("Leaf","Microbial necromass","Root","Control")) +
  facet_wrap(~LUI) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=12)) 
lco2s

Fig.ghg <- ggarrange(lco2s,ln2os,
                        ncol = 2, nrow = 1, widths = c(1,1), common.legend = TRUE, align = "h")

Fig.lab <- ggarrange(bulkCN,CN,
                     ncol = 2, nrow = 1, widths = c(1,1), common.legend = TRUE, align = "h")

png(file="ghg.png", units="in", width=10, height=4, res=300)
Fig.ghg
dev.off()


png(file="lab.png", units="in", width=10, height=4, res=300)
Fig.lab
dev.off()


#############################################################################
## not used (sub on x-axis)

labCN <- read.csv("labfig.csv", header=TRUE, sep=",")
limits <- aes(ymax = Immob + se, ymin=Immob - se)

png(file="labCN2.png", units="in", width=5, height=4, res=300)
#png(file="labCN.pdf", units="in", width=5, height=4, res=300)
plt1 <- ggplot(labCN, aes(fill=LUI, y=Immob, x=Necro))
plt1 + 
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = LUI)) + 
  scale_x_discrete(name="Substrate",
                   labels=c("Roots","Leaves","Necromass")) +
  scale_y_continuous(name=expression(Carbon~immobilized~(mg~g^-1~soil))) +
  geom_errorbar(limits, 
                position= position_dodge(width=0.9), 
                width=0.8) +
  scale_fill_manual(name = "Grassland \n management\n intensity",
                    labels = c("High", "Abandoned"),
                    values = c("#edf8b1","#2c7fb8")) +  #"#7fcdbb"
  facet_grid(.~Pool) +
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)) +
  theme_bw( )
#text = element_text(family = "Helvetica", 
#color="black", 
#size=20),


dev.off()

### version with just fine soil

labCN <- read.csv("labfig.csv", header=TRUE, sep=",")
labCNf <- labCN[labCN$Pool == "Mineral-associated",]
limits.cnf <- aes(ymax = Immob + se, ymin=Immob - se)

#png(file="labCN2.png", units="in", width=5, height=4, res=300)
png(file="labCN.png", units="in", width=4, height=4, res=300)
plt1 <- ggplot(labCNf, aes(fill=Necro, y=Immob, x=LUI))
plt1 + 
  geom_bar(position= position_dodge(width=0.9), 
           stat="identity", 
           aes(fill = Necro)) + 
  scale_x_discrete(name="Grassland Management Intensity",
                   labels=c("High","Low")) +
  scale_y_continuous(name=expression(Mineral-associated~necromass-C~(mg~g^-1~soil))) +
  geom_errorbar(limits, 
                position= position_dodge(width=0.9), 
                width=0.8) +
  scale_fill_manual(name = "Necromass \n type",
                    labels = c("Roots", "Leaves", "Microbial"),
                    values = c("#edf8b1","#7fcdbb","#2c7fb8")) +  
  theme(text = element_text(family = "Helvetica", 
                            color="black", 
                            size=12)) +
  theme_bw( ) 
#facet_grid(.~Pool)

dev.off()

