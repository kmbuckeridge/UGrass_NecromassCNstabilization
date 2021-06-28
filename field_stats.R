setwd("/Users/katebuckeridge/OneDrive - University of Edinburgh/R/UGrass/CNstab")
setwd("C:/Users/kbuckeri/OneDrive - University of Edinburgh/R/UGrass/CNstab") ## at work

citation() # also citation("stats")
# R Core Team (2020). R: A language and environment for statistical computing. R
# Foundation for Statistical Computing, Vienna, Austria. URL
# https://www.R-project.org/.

citation("nlme")
# Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2020). _nlme: Linear and
# Nonlinear Mixed Effects Models_. R package version 3.1-144, <URL:
#   https://CRAN.R-project.org/package=nlme>.

citation("ggplot2")
# H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
# York, 2016.

# full data
field <- read.csv("field_depth.csv", header=TRUE, sep=",") 
fieldall <- read.csv("field.csv", header=TRUE, sep=",") #this is depths combined

#separate by depth (I am not interested in whether surface and sub-surface or different, just whether they have treatment diffs and time patterns)
field.surf <- field[field$depth == "asurface",]
field.sub <- field[field$depth == "subsurface",]


############################################################################################################################################ 

cor.test(field.surf$soilCf,field.surf$time) # -0.7023644, p-value = 4.386e-07
cor.test(field.sub$soilCf,field.sub$time) # -0.2624903, p-value = 0.1018
cor.test(field.surf$soilNf,field.surf$time) # -0.330566 , p-value = 0.03722
cor.test(field.sub$soilNf,field.sub$time) # 0.007942878, p-value = 0.9612
cor.test(field.surf$MB13C,field.surf$time) # -0.589175, p-value = 6.341e-05
cor.test(field.sub$MB13C,field.sub$time) # -0.5001698, p-value = 0.001014
cor.test(field.surf$MB15N,field.surf$time) # -0.5472689, p-value = 0.0002577
cor.test(field.sub$MB15N,field.sub$time) # -0.1446565, p-value = 0.3732
cor.test(fieldall$soilCf,fieldall$time) #-0.7033705, p-value = 4.154e-07
cor.test(fieldall$soilNf,fieldall$time) #-0.3508044 , p-value = 0.02646
cor.test(fieldall$MB13C,fieldall$time) # -0.720941, p-value = 1.547e-07
cor.test(fieldall$MB15N,fieldall$time) # -0.5785553, p-value = 9.212e-05

### We use nlme::lme because at present it is the only easy way to allow for temporal autocorrelation in a LMM in R.(Ben Bolker)

# http://bbolker.github.io/mixedmodels-misc/ecostats_chap.html
# Belshe (2013) Ecology Letters, tundra C over time, I use this model and rationale

# corCAR1, which implements a continuous-time first-order autocorrelation model (i.e. autocorrelation declines exponentially with time - this is essentially creating a flexible time covariate according to the auto-correlation structure), works because we have missing values in the data. The more standard discrete-time autocorrelation models (lme offers corAR1 for a first-order model and corARMA for a more general model) don't work with missing data.

################# Surface carbon on fines ###########################

library(nlme)

m.scf = lme(soilCf ~ LUI*time,
              random = ~ 1 + time|plot,
              correlation = corCAR1(form=~time|plot),
              data = field.surf,
              method = "REML")

#assess data shape
library(ggpubr)
library(gridExtra)

grid.arrange(plot(m.scf,type=c("p","smooth")),
             plot(m.scf,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.scf,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.scf,abline=c(0,1)))

# residual vs fitted (non-linear patterns) should: bounce randomly around 0, form horizontal, no outliers
# scale-location (homoscedasticity = equal variance) should see: equal spread along line, line horizontal
# Residuals vs time looks for influential outliers
# normal Q-Q (normal distribution) should: be on the line or there is skew

summary(m.scf)

# Linear mixed-effects model fit by REML
# Data: field.surf 
# AIC      BIC    logLik
# 301.7151 315.9668 -141.8576
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 1.054362e-03 (Intr)
# time        2.635317e-06 0     
# Residual    8.195659e+00       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: soilCf ~ LUI * time 
#                          Value Std.Error DF   t-value p-value
# (Intercept)           44.51221  2.536207 28 17.550702  0.0000
# LUILow_intensity      -9.33845  3.586738  8 -2.603606  0.0314
# time                  -0.10524  0.019811 28 -5.312353  0.0000
# LUILow_intensity:time  0.02682  0.028017 28  0.957377  0.3466
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.691  0.489       
# LUILow_intensity:time  0.489 -0.691 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6979510 -0.7259343 -0.1333524  0.4694279  1.9434010 
# 
# Number of Observations: 40
# Number of Groups: 10 

# this is a graphical version of the autocorrelation factor - bigger (absolute) numbers (approaching or >1) imply autocorrelation. The first lag (0) is always = 1.
grid.arrange(plot(ACF(m.scf),alpha=0.05),
             plot(ACF(m.scf,resType="normalized"),alpha=0.05),
             nrow=1)

# also, R handbook: https://rcompanion.org/handbook/I_09.html
# this is a data version of the autocorrelation factor - bigger (absolute) numbers (approching or >1) imply autocorrelation.  The first lag (0) is always = 1. 
# This output suggests not much AC
ACF(m.scf,
    form = ~ time | plot)

# lag         ACF
# 1   0  1.00000000
# 2   1 -0.02993632
# 3   2 -0.24289619
# 4   3  0.24241723

library(emmeans)
emmeans(m.scf, pairwise ~ LUI : time)
# $emmeans
# LUI            time emmean     SE df lower.CL upper.CL
# High intensity 88.5  0.753 0.0392  9    0.664    0.842
# Low intensity  88.5  0.604 0.0392  8    0.514    0.694
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                                 estimate     SE df t.ratio p.value
# High intensity,88.5 - Low intensity,88.5    0.149 0.0554  8 2.687   0.0276

# can also do full pairwise with all times:
emmeans(m.scf, pairwise ~ LUI : time,at=list(time=c(1,2,3,4)))

# $emmeans
# LUI            time emmean     SE df lower.CL upper.CL
# High intensity    1  0.950 0.0540  9    0.828    1.072
# Low intensity     1  0.751 0.0540  8    0.626    0.875
# High intensity    2  0.948 0.0537  9    0.826    1.069
# Low intensity     2  0.749 0.0537  8    0.625    0.873
# High intensity    3  0.945 0.0534  9    0.825    1.066
# Low intensity     3  0.747 0.0534  8    0.624    0.870
# High intensity    4  0.943 0.0531  9    0.823    1.063
# Low intensity     4  0.746 0.0531  8    0.623    0.868
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                            estimate       SE df t.ratio p.value
# High intensity,1 - Low intensity,1   0.19918 0.076309  8  2.610  0.2737 
# High intensity,1 - High intensity,2  0.00225 0.000424 28  5.312  0.0003 
# High intensity,1 - Low intensity,2   0.20086 0.076104  8  2.639  0.2644 
# High intensity,1 - High intensity,3  0.00450 0.000848 28  5.312  0.0003 
# High intensity,1 - Low intensity,3   0.20253 0.075901  8  2.668  0.2552 
# High intensity,1 - High intensity,4  0.00675 0.001271 28  5.312  0.0003 
# High intensity,1 - Low intensity,4   0.20421 0.075699  8  2.698  0.2464 
# Low intensity,1 - High intensity,2  -0.19693 0.076104  8 -2.588  0.2811 
# Low intensity,1 - Low intensity,2    0.00168 0.000424 28  3.958  0.0097 
# Low intensity,1 - High intensity,3  -0.19468 0.075901  8 -2.565  0.2888 
# Low intensity,1 - Low intensity,3    0.00335 0.000848 28  3.958  0.0097 
# Low intensity,1 - High intensity,4  -0.19243 0.075699  8 -2.542  0.2967 
# Low intensity,1 - Low intensity,4    0.00503 0.001271 28  3.958  0.0097 
# High intensity,2 - Low intensity,2   0.19861 0.075898  8  2.617  0.2716 
# High intensity,2 - High intensity,3  0.00225 0.000424 28  5.312  0.0003 
# High intensity,2 - Low intensity,3   0.20028 0.075694  8  2.646  0.2622 
# High intensity,2 - High intensity,4  0.00450 0.000848 28  5.312  0.0003 
# High intensity,2 - Low intensity,4   0.20196 0.075493  8  2.675  0.2531 
# Low intensity,2 - High intensity,3  -0.19635 0.075694  8 -2.594  0.2790 
# Low intensity,2 - Low intensity,3    0.00168 0.000424 28  3.958  0.0097 
# Low intensity,2 - High intensity,4  -0.19410 0.075493  8 -2.571  0.2867 
# Low intensity,2 - Low intensity,4    0.00335 0.000848 28  3.958  0.0097 
# High intensity,3 - Low intensity,3   0.19803 0.075490  8  2.623  0.2695 
# High intensity,3 - High intensity,4  0.00225 0.000424 28  5.312  0.0003 
# High intensity,3 - Low intensity,4   0.19971 0.075288  8  2.653  0.2601 
# Low intensity,3 - High intensity,4  -0.19578 0.075288  8 -2.600  0.2769 
# Low intensity,3 - Low intensity,4    0.00168 0.000424 28  3.958  0.0097 
# High intensity,4 - Low intensity,4   0.19746 0.075085  8  2.630  0.2674 
# 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 8 estimates 

################### Deep carbon on fines ################

m.dcf = lme(soilCf ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = field.sub,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

#assess data shape
grid.arrange(plot(m.dcf,type=c("p","smooth")),
             plot(m.dcf,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.dcf,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.dcf,abline=c(0,1)))

summary(m.dcf)
# same result as ben Bolker: The results generally look sensible: the only warning sign is that the among-site variation in baseline NEE ((Intercept)) and the among-site variation in slope are perfectly correlated (i.e., the -1 term in the Corr column under Random effects). We weren’t happy with this, but we kept the full random effects model anyway. The alternative, dropping the random effect of year, seemed unacceptable in our case. 
# also, added the control to allow for more iterations (it timed out otherwise)

# Linear mixed-effects model fit by REML
# Data: field.sub 
# AIC      BIC    logLik
# 274.1493 288.4009 -128.0746
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev     Corr  
# (Intercept) 4.26080901 (Intr)
# time        0.01447183 -1    
# Residual    5.00354446       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: soilCf ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)            7.151320  2.455278 28  2.912631  0.0070
# LUILow_intensity      14.478140  3.472287  8  4.169626  0.0031
# time                   0.005993  0.013717 28  0.436867  0.6656
# LUILow_intensity:time -0.057605  0.019399 28 -2.969440  0.0061
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.750  0.531       
# LUILow_intensity:time  0.531 -0.750 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.46560590 -0.31658296 -0.01004712  0.26173615  3.77313782 
# 
# Number of Observations: 40
# Number of Groups: 10 

library(emmeans)
emmeans(m.dcf, pairwise ~ LUI : time)
#$emmeans
# LUI            time emmean   SE df lower.CL upper.CL
# High_intensity 88.5   7.68 1.74  9     3.75     11.6
# Low_intensity  88.5  17.06 1.74  8    13.05     21.1
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                                 estimate   SE df t.ratio p.value
# High_intensity,88.5 - Low_intensity,88.5    -9.38 2.46  8 -3.811  0.0052 
# 
# Degrees-of-freedom method: containment  

# can also do full pairwise with all times:
emmeans(m.dcf, pairwise ~ LUI : time,at=list(time=c(1,2,3,4)))
# $emmeans
# LUI            time emmean   SE df lower.CL upper.CL
# High_intensity    1   7.16 2.45  9     1.63     12.7
# Low_intensity     1  21.58 2.45  8    15.94     27.2
# High_intensity    2   7.16 2.43  9     1.66     12.7
# Low_intensity     2  21.53 2.43  8    15.91     27.1
# High_intensity    3   7.17 2.42  9     1.68     12.7
# Low_intensity     3  21.47 2.42  8    15.88     27.1
# High_intensity    4   7.18 2.41  9     1.71     12.6
# Low_intensity     4  21.42 2.41  8    15.86     27.0
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                             estimate     SE df t.ratio p.value
# High_intensity,1 - Low_intensity,1  -14.42053 3.4578  8 -4.170  0.0382 
# High_intensity,1 - High_intensity,2  -0.00599 0.0137 28 -0.437  0.9998 
# High_intensity,1 - Low_intensity,2  -14.36892 3.4505  8 -4.164  0.0385 
# High_intensity,1 - High_intensity,3  -0.01199 0.0274 28 -0.437  0.9998 
# High_intensity,1 - Low_intensity,3  -14.31731 3.4433  8 -4.158  0.0388 
# High_intensity,1 - High_intensity,4  -0.01798 0.0412 28 -0.437  0.9998 
# High_intensity,1 - Low_intensity,4  -14.26570 3.4362  8 -4.152  0.0391 
# Low_intensity,1 - High_intensity,2   14.41454 3.4505  8  4.178  0.0379 
# Low_intensity,1 - Low_intensity,2     0.05161 0.0137 28  3.763  0.0157 
# Low_intensity,1 - High_intensity,3   14.40855 3.4433  8  4.184  0.0375 
# Low_intensity,1 - Low_intensity,3     0.10322 0.0274 28  3.763  0.0157 
# Low_intensity,1 - High_intensity,4   14.40256 3.4362  8  4.191  0.0372 
# Low_intensity,1 - Low_intensity,4     0.15484 0.0412 28  3.763  0.0157 
# High_intensity,2 - Low_intensity,2  -14.36293 3.4433  8 -4.171  0.0382 
# High_intensity,2 - High_intensity,3  -0.00599 0.0137 28 -0.437  0.9998 
# High_intensity,2 - Low_intensity,3  -14.31132 3.4361  8 -4.165  0.0385 
# High_intensity,2 - High_intensity,4  -0.01199 0.0274 28 -0.437  0.9998 
# High_intensity,2 - Low_intensity,4  -14.25970 3.4289  8 -4.159  0.0388 
# Low_intensity,2 - High_intensity,3   14.35694 3.4361  8  4.178  0.0378 
# Low_intensity,2 - Low_intensity,3     0.05161 0.0137 28  3.763  0.0157 
# Low_intensity,2 - High_intensity,4   14.35094 3.4289  8  4.185  0.0375 
# Low_intensity,2 - Low_intensity,4     0.10322 0.0274 28  3.763  0.0157 
# High_intensity,3 - Low_intensity,3  -14.30532 3.4288  8 -4.172  0.0381 
# High_intensity,3 - High_intensity,4  -0.00599 0.0137 28 -0.437  0.9998 
# High_intensity,3 - Low_intensity,4  -14.25371 3.4216  8 -4.166  0.0384 
# Low_intensity,3 - High_intensity,4   14.29933 3.4216  8  4.179  0.0378 
# Low_intensity,3 - Low_intensity,4     0.05161 0.0137 28  3.763  0.0157 
# High_intensity,4 - Low_intensity,4  -14.24772 3.4144  8 -4.173  0.0381 
# 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 8 estimates 


ACF(m.dcf,
    form = ~ time | plot)
# lag        ACF
# 1   0  1.0000000
# 2   1 -0.2049763
# 3   2 -0.6934181
# 4   3  0.7044272

######################## Total (surface + deep) C on fines #########################

m.cf = lme(soilCf ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = fieldall,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

summary(m.cf)
# Linear mixed-effects model fit by REML
# Data: fieldall 
# AIC      BIC   logLik
# 46.57661 60.82828 -14.2883
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 0.0451564274 (Intr)
# time        0.0001403291 -0.995
# Residual    0.2346610244       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: soilCf ~ LUI * time 
# Value  Std.Error DF   t-value p-value
# (Intercept)            1.1051047 0.07537328 28 14.661757  0.0000
# LUILow_intensity       0.1099401 0.10659392  8  1.031392  0.3325
# time                  -0.0021230 0.00057069 28 -3.720023  0.0009
# LUILow_intensity:time -0.0006585 0.00080708 28 -0.815855  0.4215
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.691  0.489       
# LUILow_intensity:time  0.489 -0.691 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.3215346 -0.5894007 -0.1734697  0.3232558  3.6953518 
# 
# Number of Observations: 40
# Number of Groups: 10 

m.cf.2 = lme(soilCf ~ LUI*time*depth,
           random = ~ 1 + time|plot,
           correlation = corCAR1(form=~time|plot),
           data = field,
           method = "REML",
           control=list(maxIter=10000, niterEM=10000))
## Error in Initialize.corCAR1(X[[i]], ...) : covariate must have unique values within groups for "corCAR1" objects

######################## Surface N on fines ##############################

m.snf = lme(soilNf ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = field.surf,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

#assess data shape
grid.arrange(plot(m.snf,type=c("p","smooth")),
             plot(m.snf,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.snf,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.snf,abline=c(0,1)))

summary(m.snf)
# Linear mixed-effects model fit by REML
# Data: field.surf 
# AIC     BIC    logLik
# 223.0164 237.268 -102.5082
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev     Corr  
# (Intercept) 1.11447356 (Intr)
# time        0.01170543 0.075 
# Residual    2.32403542       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: soilNf ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)           18.460985 0.8750109 28 21.098005  0.0000
# LUILow_intensity      -5.584568 1.2374523  8 -4.512956  0.0020
# time                  -0.016871 0.0076787 28 -2.197064  0.0365
# LUILow_intensity:time  0.004939 0.0108593 28  0.454784  0.6528
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.386  0.273       
# LUILow_intensity:time  0.273 -0.386 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6977644 -0.5320639 -0.0505189  0.5729580  2.4988920 
# 
# Number of Observations: 40
# Number of Groups: 10 


######################## Deep nitrogen on fines ##############################

m.dnf = lme(soilNf ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = field.sub,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

grid.arrange(plot(m.dnf,type=c("p","smooth")),
             plot(m.dnf,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.dnf,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.dnf,abline=c(0,1)))

summary(m.dnf)
# Linear mixed-effects model fit by REML
# Data: field.sub 
# AIC      BIC   logLik
# 187.0536 201.3053 -84.5268
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev      Corr  
# (Intercept) 0.899109884 (Intr)
# time        0.001321441 -0.997
# Residual    1.538426497       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: soilNf ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)            2.309325 0.6231608 28  3.705825   9e-04
# LUILow_intensity       5.365338 0.8812824  8  6.088103   3e-04
# time                   0.012820 0.0037654 28  3.404713   2e-03
# LUILow_intensity:time -0.025204 0.0053251 28 -4.733102   1e-04
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.623  0.440       
# LUILow_intensity:time  0.440 -0.623 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.59665743 -0.35933228 -0.04408618  0.34102358  2.85313817 
# 
# Number of Observations: 40
# Number of Groups: 10 

emmeans(m.dnf, pairwise ~ LUI : time)
#$emmeans
# LUI            time emmean    SE df lower.CL upper.CL
# High_intensity 88.5   3.35 0.477  9     2.27     4.43
# Low_intensity  88.5   6.40 0.477  8     5.30     7.50
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                                 estimate    SE df t.ratio p.value
# High_intensity,88.5 - Low_intensity,88.5    -3.05 0.675  8 -4.517  0.0020 
# 
# Degrees-of-freedom method: containment 

######################## Total (surface + deep) N on fines #########################

m.nf = lme(soilNf ~ LUI*time,
           random = ~ 1 + time|plot,
           correlation = corCAR1(form=~time|plot),
           data = fieldall,
           method = "REML",
           control=list(maxIter=10000, niterEM=10000))

summary(m.nf)

# Linear mixed-effects model fit by REML
# Data: fieldall 
# AIC        BIC   logLik
# -14.55644 -0.3047714 16.27822
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 0.0238347290 (Intr)
# time        0.0003293945 0.276 
# Residual    0.0919652703       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: soilNf ~ LUI * time 
# Value  Std.Error DF   t-value p-value
# (Intercept)            0.6431036 0.03039000 28 21.161686  0.0000
# LUILow_intensity      -0.0067879 0.04297795  8 -0.157940  0.8784
# time                  -0.0001254 0.00026668 28 -0.470278  0.6418
# LUILow_intensity:time -0.0006275 0.00037714 28 -1.663757  0.1073
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.486  0.344       
# LUILow_intensity:time  0.344 -0.486 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.9100877 -0.4849779 -0.2221044  0.4007441  2.3874370 
# 
# Number of Observations: 40
# Number of Groups: 10 


#################### 13C in surface microbial biomass ###############

m.smb13c = lme(MB13C ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = field.surf,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

#assess shape
grid.arrange(plot(m.smb13c,type=c("p","smooth")),
             plot(m.smb13c,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.smb13c,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.smb13c,abline=c(0,1)))

summary(m.smb13c)
# Linear mixed-effects model fit by REML
# Data: field.surf 
# AIC      BIC    logLik
# 341.6011 355.8528 -161.8005
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev      Corr  
# (Intercept)  4.50674494 (Intr)
# time         0.02322227 -1    
# Residual    13.95193874       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: MB13C ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)            46.23103  4.764789 28  9.702641  0.0000
# LUILow_intensity      -16.54587  6.738429  8 -2.455450  0.0396
# time                   -0.14187  0.035288 28 -4.020297  0.0004
# LUILow_intensity:time   0.04618  0.049905 28  0.925386  0.3627
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.723  0.511       
# LUILow_intensity:time  0.511 -0.723 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.0500625 -0.3042638 -0.1049641  0.2376921  2.7081903 
# 
# Number of Observations: 40
# Number of Groups: 10 

#################### 13C in deep microbial biomass ###############

m.dmb13c = lme(MB13C ~ LUI*time,
               random = ~ 1 + time|plot,
               correlation = corCAR1(form=~time|plot),
               data = field.sub,
               method = "REML",
               control=list(maxIter=10000, niterEM=10000))

#assess shape
grid.arrange(plot(m.dmb13c,type=c("p","smooth")),
             plot(m.dmb13c,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.dmb13c,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.dmb13c,abline=c(0,1)))

summary(m.dmb13c)
#Linear mixed-effects model fit by REML
# Data: field.sub 
# AIC      BIC    logLik
# 268.7435 282.9951 -125.3717
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev     Corr  
# (Intercept) 3.87934372 (Intr)
# time        0.01442398 -1    
# Residual    4.66810195       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: MB13C ~ LUI * time 
#                           Value Std.Error DF   t-value p-value
# (Intercept)           10.007945  2.257580 28  4.433041  0.0001
# LUILow_intensity       6.080554  3.192700  8  1.904517  0.0933
# time                  -0.018092  0.012998 28 -1.391974  0.1749
# LUILow_intensity:time -0.033666  0.018381 28 -1.831539  0.0777
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.765  0.541       
# LUILow_intensity:time  0.541 -0.765 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.82960875 -0.54596013 -0.04799051  0.51455793  2.43156076 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### 13C in total MBC ###################

m.mbc = lme(MB13C ~ LUI*time,
           random = ~ 1 + time|plot,
           correlation = corCAR1(form=~time|plot),
           data = fieldall,
           method = "REML",
           control=list(maxIter=10000, niterEM=10000))

summary(m.mbc)
# Linear mixed-effects model fit by REML
# Data: fieldall 
# AIC      BIC    logLik
# 337.7676 352.0192 -159.8838
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 4.148361e-04 (Intr)
# time        4.023518e-06 -0.005
# Residual    1.352221e+01       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: MB13C ~ LUI * time 
#                           Value Std.Error DF   t-value p-value
# (Intercept)            56.23898  4.184546 28 13.439684  0.0000
# LUILow_intensity      -10.46532  5.917842  8 -1.768435  0.1150
# time                   -0.15996  0.032686 28 -4.893793  0.0000
# LUILow_intensity:time   0.01251  0.046225 28  0.270737  0.7886
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.691  0.489       
# LUILow_intensity:time  0.489 -0.691 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.0869649 -0.3588657 -0.1282293  0.1849644  3.2222867 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### 15N in surface microbial biomass ###############

m.smb15n = lme(MB15N ~ LUI*time,
               random = ~ 1 + time|plot,
               correlation = corCAR1(form=~time|plot),
               data = field.surf,
               method = "REML",
               control=list(maxIter=10000, niterEM=10000))

#assess shape
grid.arrange(plot(m.smb15n,type=c("p","smooth")),
             plot(m.smb15n,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.smb15n,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.smb15n,abline=c(0,1)))

summary(m.smb15n)
# Linear mixed-effects model fit by REML
# Data: field.surf 
# AIC      BIC    logLik
# 312.0044 326.2561 -147.0022
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
#              StdDev       Corr  
# (Intercept) 1.415295e-04 (Intr)
# time        1.920959e-06 -0.004
# Residual    9.454694e+00       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: MB15N ~ LUI * time 
#                            Value Std.Error DF   t-value p-value
# (Intercept)            30.360161  2.925824 28 10.376618  0.0000
# LUILow_intensity      -11.351102  4.137740  8 -2.743309  0.0253
# time                   -0.085012  0.022854 28 -3.719773  0.0009
# LUILow_intensity:time   0.027229  0.032321 28  0.842474  0.4067
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.691  0.489       
# LUILow_intensity:time  0.489 -0.691 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.4861059 -0.4742320 -0.1863099  0.5692957  2.5736570 
# 
# Number of Observations: 40
# Number of Groups: 10 

#################### 15N in deep microbial biomass ###############

m.dmb15n = lme(MB15N ~ LUI*time,
               random = ~ 1 + time|plot,
               correlation = corCAR1(form=~time|plot),
               data = field.sub,
               method = "REML",
               control=list(maxIter=10000, niterEM=10000))

#assess shape
grid.arrange(plot(m.dmb15n,type=c("p","smooth")),
             plot(m.dmb15n,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.dmb15n,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.dmb15n,abline=c(0,1)))
#there is one strong outlier here on the first sample day

field.sub$log_mb15n <- log(field.sub$MB15N + 1)

m.dmb15n.1 = lme(log_mb15n ~ LUI*time,
               random = ~ 1 + time|plot,
               correlation = corCAR1(form=~time|plot),
               data = field.sub,
               method = "REML",
               control=list(maxIter=10000, niterEM=10000))
#assess shape
grid.arrange(plot(m.dmb15n.1,type=c("p","smooth")),
             plot(m.dmb15n.1,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.dmb15n.1,resid(.,type="pearson")~time,
                  type=c("p","smooth")),
             qqnorm(m.dmb15n.1,abline=c(0,1)))


summary(m.dmb15n.1)
# Linear mixed-effects model fit by REML
# Data: field.sub 
# AIC      BIC    logLik
# 239.6428 253.8945 -110.8214
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
#             StdDev     Corr  
# (Intercept) 3.09046669 (Intr)
# time        0.01576336 -1    
# Residual    3.04580810       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: MB15N ~ LUI * time 
#                           Value Std.Error DF    t-value p-value
# (Intercept)            1.825508 1.6728995 28  1.0912242  0.2845
# LUILow_intensity       3.977658 2.3658372  8  1.6812897  0.1312
# time                   0.005221 0.0101932 28  0.5122347  0.6125
# LUILow_intensity:time -0.022445 0.0144154 28 -1.5570021  0.1307
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.853  0.603       
# LUILow_intensity:time  0.603 -0.853 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6235809 -0.2866085 -0.1086815  0.2793336  4.2045374 
# 
# Number of Observations: 40
# Number of Groups: 10 

summary(m.dmb15n.1)  ## improved model (smaller residuals) but still not significant
# Linear mixed-effects model fit by REML
# Data: field.sub 
# AIC      BIC   logLik
# 103.0374 117.2891 -42.5187
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
#             StdDev      Corr  
# (Intercept) 0.694716404 (Intr)
# time        0.003477589 -1    
# Residual    0.419811780       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: log_mb15n ~ LUI * time 
#                            Value Std.Error DF   t-value p-value
# (Intercept)            0.8759236 0.3367548 28  2.601072  0.0147
# LUILow_intensity       0.6310100 0.4762432  8  1.324974  0.2218
# time                   0.0022782 0.0018570 28  1.226786  0.2301
# LUILow_intensity:time -0.0036806 0.0026262 28 -1.401501  0.1721
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.918  0.649       
# LUILow_intensity:time  0.649 -0.918 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6978336 -0.5505506 -0.1622284  0.5927323  2.2257057 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### 15N in total MBN ###################

m.mbn = lme(MB15N ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = fieldall,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

summary(m.mbn)
# Linear mixed-effects model fit by REML
# Data: fieldall 
# AIC      BIC    logLik
# 315.9472 330.1989 -148.9736
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev      Corr  
# (Intercept) 0.648986663 (Intr)
# time        0.003572002 -0.992
# Residual    9.976975525       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: MB15N ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)           32.18567  3.101060 28 10.378926  0.0000
# LUILow_intensity      -7.37344  4.385561  8 -1.681300  0.1312
# time                  -0.07979  0.024169 28 -3.301311  0.0026
# LUILow_intensity:time  0.00478  0.034181 28  0.139976  0.8897
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.693  0.490       
# LUILow_intensity:time  0.490 -0.693 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.30848311 -0.35568803 -0.08815416  0.51556774  2.52656084 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### total MBN ###################

m.mbn2 = lme(mbn ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = fieldall,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

summary(m.mbn2) #no diff between MBN in high and low LUI

# Linear mixed-effects model fit by REML
# Data: fieldall 
# AIC      BIC    logLik
# 206.0985 220.3502 -94.04926
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 4.119230e-05 (Intr)
# time        3.998474e-07 -0.002
# Residual    2.171887e+00       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: mbn ~ LUI * time 
#                        Value Std.Error DF   t-value p-value
# (Intercept)            6.080700 0.6721064 28  9.047227  0.0000
# LUILow_intensity      -0.621568 0.9505020  8 -0.653936  0.5315
# time                  -0.003053 0.0052500 28 -0.581521  0.5655
# LUILow_intensity:time  0.001036 0.0074246 28  0.139472  0.8901
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.691  0.489       
# LUILow_intensity:time  0.489 -0.691 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6030336 -0.6409396 -0.2096943  0.4430190  2.6363593 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### surface MBN ###################

m.smbn2 = lme(mbn ~ LUI*time,
              random = ~ 1 + time|plot,
              correlation = corCAR1(form=~time|plot),
              data = field.surf,
              method = "REML",
              control=list(maxIter=10000, niterEM=10000))

summary(m.smbn2)
# Linear mixed-effects model fit by REML
# Data: field.surf 
# AIC      BIC    logLik
# 204.4871 218.7388 -93.24356
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 3.549639e-05 (Intr)
# time        3.695342e-07 -0.003
# Residual    2.123819e+00       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: mbn ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)            4.442651 0.6572314 28  6.759646  0.0000
# LUILow_intensity      -0.825499 0.9294655  8 -0.888143  0.4004
# time                  -0.005905 0.0051338 28 -1.150243  0.2598
# LUILow_intensity:time  0.001761 0.0072602 28  0.242563  0.8101
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.691  0.489       
# LUILow_intensity:time  0.489 -0.691 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.8168808 -0.5890267 -0.2559443  0.6668077  2.0577759 
# 
# Number of Observations: 40
# Number of Groups: 10

################### deep MBN ###################

m.dmbn2 = lme(mbn ~ LUI*time,
              random = ~ 1 + time|plot,
              correlation = corCAR1(form=~time|plot),
              data = field.sub,
              method = "REML",
              control=list(maxIter=10000, niterEM=10000))

summary(m.dmbn2)
# Linear mixed-effects model fit by REML
# Data: field.sub 
# AIC      BIC    logLik
# 116.1202 130.3718 -49.06009
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev      Corr  
# (Intercept) 0.588336257 (Intr)
# time        0.001710343 -0.849
# Residual    0.518937876       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: mbn ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)            1.6380485 0.3082479 28  5.314062  0.0000
# LUILow_intensity       0.2039311 0.4359284  8  0.467809  0.6524
# time                   0.0028521 0.0014692 28  1.941265  0.0624
# LUILow_intensity:time -0.0007255 0.0020778 28 -0.349197  0.7296
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.685  0.484       
# LUILow_intensity:time  0.484 -0.685 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.65718003 -0.51420012 -0.05659008  0.48800398  1.87916552 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### total MBC ###################

m.mbc2 = lme(mbc ~ LUI*time,
             random = ~ 1 + time|plot,
             correlation = corCAR1(form=~time|plot),
             data = fieldall,
             method = "REML",
             control=list(maxIter=10000, niterEM=10000))

summary(m.mbc2) #MBC was lower in low LUI

# Linear mixed-effects model fit by REML
# Data: fieldall 
# AIC      BIC    logLik
# 357.0539 371.3056 -169.5269
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept)  6.575361676 (Intr)
# time         0.002355652 -0.888
# Residual    16.800373141       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: mbc ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)            72.86135  5.972994 28 12.198463  0.0000
# LUILow_intensity      -20.60779  8.447089  8 -2.439632  0.0406
# time                   -0.01800  0.040624 28 -0.443092  0.6611
# LUILow_intensity:time   0.00491  0.057451 28  0.085458  0.9325
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.613  0.433       
# LUILow_intensity:time  0.433 -0.613 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.51554263 -0.69714585 -0.08218904  0.48384406  2.21703259 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### surface MBC ###################

m.smbc2 = lme(mbc ~ LUI*time,
             random = ~ 1 + time|plot,
             correlation = corCAR1(form=~time|plot),
             data = field.surf,
             method = "REML",
             control=list(maxIter=10000, niterEM=10000))

summary(m.smbc2)
# Linear mixed-effects model fit by REML
# Data: field.surf 
# AIC      BIC    logLik
# 348.3247 362.5764 -165.1624
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 4.401902e-04 (Intr)
# time        4.032328e-06 -0.004
# Residual    1.565767e+01       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: mbc ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)            43.77645  4.845380 28  9.034679  0.0000
# LUILow_intensity      -15.02256  6.852401  8 -2.192306  0.0597
# time                   -0.03916  0.037848 28 -1.034641  0.3097
# LUILow_intensity:time   0.01281  0.053525 28  0.239243  0.8127
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.691  0.489       
# LUILow_intensity:time  0.489 -0.691 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.92049133 -0.63062683 -0.04874818  0.57738148  1.96179446 
# 
# Number of Observations: 40
# Number of Groups: 10 

################### deep MBC ###################

m.dmbc2 = lme(mbc ~ LUI*time,
              random = ~ 1 + time|plot,
              correlation = corCAR1(form=~time|plot),
              data = field.sub,
              method = "REML",
              control=list(maxIter=10000, niterEM=10000))

summary(m.dmbc2)

# Linear mixed-effects model fit by REML
# Data: field.sub 
# AIC      BIC    logLik
# 275.0331 289.2847 -128.5165
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev      Corr  
# (Intercept) 7.541179796 (Intr)
# time        0.005680621 -0.985
# Residual    4.299677627       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: mbc ~ LUI * time 
# Value Std.Error DF   t-value p-value
# (Intercept)           29.084903  3.625505 28  8.022304  0.0000
# LUILow_intensity      -5.585224  5.127238  8 -1.089324  0.3077
# time                   0.021159  0.010699 28  1.977618  0.0579
# LUILow_intensity:time -0.007896  0.015131 28 -0.521835  0.6059
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.464  0.328       
# LUILow_intensity:time  0.328 -0.464 -0.707
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.8585963 -0.7376387  0.1080798  0.6660063  1.9884672 
# 
# Number of Observations: 40
# Number of Groups: 10 

#######################################################################

############# GHG

ghg <- read.csv("field_GHG.csv", header=TRUE, sep=",")

m.co2 = lme(cumCO2sub ~ Day*LUI,
              random = ~ 1 + Day|Plot,
              correlation = corCAR1(form=~Day|Plot),
              data = ghg,
              method = "REML",
              control=list(maxIter=10000, niterEM=10000))

#assess data shape
grid.arrange(plot(m.co2,type=c("p","smooth")),
             plot(m.co2,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.co2,resid(.,type="pearson")~Day,
                  type=c("p","smooth")),
             qqnorm(m.co2,abline=c(0,1)))

summary(m.co2)
# Linear mixed-effects model fit by REML
# Data: ghg 
# AIC      BIC    logLik
# 165.5535 185.2604 -73.77677
# 
# Random effects:
#   Formula: ~1 + Day | Plot
# Structure: General positive-definite, Log-Cholesky parametrization
#             StdDev       Corr  
# (Intercept) 0.2152264102 (Intr)
# Day         0.0009206658 0.998 
# Residual    0.6304251529       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~Day | Plot 
# Parameter estimate(s):
#   Phi 
# 0.9380711 
# Fixed effects: cumCO2sub ~ Day * LUI 
#                  Value  Std.Error DF   t-value p-value
# (Intercept)  1.1844907 0.21866269 58  5.416977  0.0000
# Day          0.0060129 0.00153606 58  3.914508  0.0002
# LUILow      -0.1927135 0.30923574  8 -0.623193  0.5505
# Day:LUILow  -0.0014418 0.00217232 58 -0.663730  0.5095
# Correlation: 
#   (Intr) Day    LUILow
# Day        -0.551              
# LUILow     -0.707  0.390       
# Day:LUILow  0.390 -0.707 -0.551
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.04087725 -1.13019215  0.06089181  0.66606493  1.49209575 
# 
# Number of Observations: 70
# Number of Groups: 10 

m.co2mg = lme(mgCO2sub ~ Day*LUI,
            random = ~ 1 + Day|Plot,
            correlation = corCAR1(form=~Day|Plot),
            data = ghg,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))
# Error in lme.formula(mgCO2sub ~ Day * LUI, random = ~1 + Day | Plot, correlation = corCAR1(form = ~Day |  : 
#                                                                                              nlminb problem, convergence error code = 1
#                                                                                            message = singular convergence (7)

m.n2o = lme(cumN2Osub ~ Day*LUI,
            random = ~ 1 + Day|Plot,
            correlation = corCAR1(form=~Day|Plot),
            data = ghg,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))
#assess data shape
grid.arrange(plot(m.n2o,type=c("p","smooth")),
             plot(m.n2o,sqrt(abs(resid(.)))~fitted(.),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(m.n2o,resid(.,type="pearson")~Day,
                  type=c("p","smooth")),
             qqnorm(m.n2o,abline=c(0,1)))

summary(m.n2o)
# Linear mixed-effects model fit by REML
# Data: ghg 
# AIC      BIC    logLik
# 627.8842 647.5911 -304.9421
# 
# Random effects:
#   Formula: ~1 + Day | Plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev      Corr  
# (Intercept) 14.01526189 (Intr)
# Day          0.07401003 1     
# Residual    18.61427632       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~Day | Plot 
# Parameter estimate(s):
#   Phi 
# 0.9091755 
# Fixed effects: cumN2Osub ~ Day * LUI 
# Value Std.Error DF   t-value p-value
# (Intercept)  0.95102  8.279526 58 0.1148637  0.9089
# Day          0.03934  0.053770 58 0.7316553  0.4673
# LUILow      33.67180 11.709019  8 2.8757151  0.0206
# Day:LUILow   0.09217  0.076042 58 1.2121492  0.2304
# Correlation: 
#   (Intr) Day    LUILow
# Day         0.074              
# LUILow     -0.707 -0.052       
# Day:LUILow -0.052 -0.707  0.074
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.97795328 -0.19936276 -0.01597218  0.08448789  4.67883475 
# 
# Number of Observations: 70
# Number of Groups: 10 

########################################################################
########################   Decay rates ################################

######## example from online:
#data from: https://stackoverflow.com/questions/19453861/how-to-fit-and-plot-exponential-decay-function-using-ggplot2-and-linear-approxim

summary(lm(log(y) ~ time, data = x, subset = Factor)) # I need the summary statistics to compare models
g1 <- ggplot(x, aes(x = time, y = y, color = Factor)) + geom_point() +expand_limits(y=0)
g1 + geom_smooth(method = "glm", family = gaussian(link="log"), start=c(5,0))

glm(y ~ Factor*time-1, data = x, family = gaussian(link = "log"), start = c(5,5,5,5,0,0,0,0))
# Coefficients:
#   FactorA       FactorB       FactorC       FactorD          time  FactorB:time  
# 3.800792      3.907178      3.884973      3.964836     -0.010518     -0.028596  
# FactorC:time  FactorD:time  
# 0.004722     -0.029597  
# 
# Degrees of Freedom: 158 Total (i.e. Null);  150 Residual
# Null Deviance:	    224800 
# Residual Deviance: 13930 	AIC: 1174

fit2 = lm(y ~ time + time:Factor - 1, data = x)
summary(fit2)
# 
# Call:
#   lm(formula = y ~ time + time:Factor - 1, data = x)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -30.60  13.13  31.94  48.49  74.94 
# 
# Coefficients: (1 not defined because of singularities)
#               Estimate Std. Error t value Pr(>|t|)
# time          0.29718    0.22266   1.335    0.184
# time:FactorA -0.22577    0.23009  -0.981    0.328
# time:FactorB  0.01546    0.34869   0.044    0.965
# time:FactorC -0.16693    0.22970  -0.727    0.468
# time:FactorD       NA         NA      NA       NA
# 
# Residual standard error: 37.89 on 154 degrees of freedom
# Multiple R-squared:  0.06086,	Adjusted R-squared:  0.03647 
# F-statistic: 2.495 on 4 and 154 DF,  p-value: 0.04519

fit1 = lm(y ~ time*Factor, data = x)
summary(fit1)
(x$pred = exp(predict(fit1)) )

############ data (same as above) ##############

# full data
field <- read.csv("field_depth.csv", header=TRUE, sep=",") 
fieldall <- read.csv("field.csv", header=TRUE, sep=",") #this is depths combined

#separate by depth (I am not interested in whether surface and sub-surface or different, just whether they have treatment diffs and time patterns)
field.surf <- field[field$depth == "asurface",]
field.sub <- field[field$depth == "subsurface",]

############ fine soil C (total) ################


### can do this with a lm, or with my model above (lme) but with log-transformed data. I think it is more correct to use the glm instead of the mixed model. Does it make sense to incorporate an autocorrelation correction when assessing decay rate...I don't think so. In other words, use the 2nd version of below:

fieldall$lnsfc <- log(fieldall$soilCf)

m.scf.ln = lme(lnsfc ~ LUI*time,
            random = ~ 1 + time|plot,
            correlation = corCAR1(form=~time|plot),
            data = field.surf,
            method = "REML",
            control=list(maxIter=10000, niterEM=10000))

summary(m.scf.ln)
# Linear mixed-effects model fit by REML
# Data: field.surf 
# AIC      BIC    logLik
# 45.18745 59.43912 -13.59372
# 
# Random effects:
#   Formula: ~1 + time | plot
# Structure: General positive-definite, Log-Cholesky parametrization
# StdDev       Corr  
# (Intercept) 0.0359571837 (Intr)
# time        0.0006964042 0.995 
# Residual    0.2128466744       
# 
# Correlation Structure: Continuous AR(1)
# Formula: ~time | plot 
# Parameter estimate(s):
#   Phi 
# 0.2 
# Fixed effects: lnsfc ~ LUI * time 
# Value  Std.Error DF  t-value p-value
# (Intercept)            3.383451 0.06780148 28 49.90232  0.0000
# LUILow_intensity      -0.243194 0.09588577  8 -2.53629  0.0349
# time                  -0.003236 0.00060142 28 -5.38129  0.0000
# LUILow_intensity:time  0.000216 0.00085054 28  0.25442  0.8010
# Correlation: 
#   (Intr) LUILw_ time  
# LUILow_intensity      -0.707              
# time                  -0.452  0.320       
# LUILow_intensity:time  0.320 -0.452 -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.66801943 -0.67432103 -0.07003746  0.62569810  1.83806648 
# 
# Number of Observations: 40
# Number of Groups: 10

#MRT for fine soil necro-C (high LUI): Adding 1 day of time decreases soil fine 13C by -0.003236 xg in High LUI
1/-0.003236
# [1] -309.0235 days
#low LUI: the effect of time on soil fine C is higher by 0.000216 xg under Low LUI compared to High LUI. In other words the slope between time and soil fine C is:
((-0.003236)+(0.000216))
# = -0.00302 xg/d in Low LUI.

# Therefore, Low LUI MRT:
1/((-0.0023700)+(0.000216))
# [1] -464.2526 days

fit.sfc = lm(lnsfc ~ time*LUI, data = fieldall)
summary(fit.sfc)
#(field$pred.sfc = exp(predict(fit.sfc)) ) # if I want to back-transform and save the predicted values
#field$conf.sfc = exp(predict(fit.sfc, interval = "confidence")) # if I want to save the confodence interval

# The output of this model will show the estimated intercept (estimate) for the reference level of LUI (high) , the estimated slope for the reference level ("time" estimate = slope for High LUI), and the difference in intercepts ("LUILow_intensity" estimate) and slopes ("time:LUILow_intensity" estimates) between the reference level and all other levels.

# Call:
# lm(formula = lnsfc ~ time * LUI, data = fieldall)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.37867 -0.14013 -0.01279  0.10668  0.58085 
# 
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            3.5294421  0.0648297  54.442  < 2e-16 ***
# time                  -0.0023700  0.0005064  -4.680 3.98e-05 ***
# LUILow_intensity       0.0889182  0.0916830   0.970    0.339    
# time:LUILow_intensity -0.0006606  0.0007162  -0.922    0.362    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2095 on 36 degrees of freedom
# Multiple R-squared:  0.6167,	Adjusted R-squared:  0.5848 
# F-statistic: 19.31 on 3 and 36 DF,  p-value: 1.241e-07

#MRT for fine soil necro-C (high LUI): Adding 1 day of time decreases soil fine 13C by -0.0023700 xg in High LUI
1/-0.0023700
# [1] -421.9409 days
#low LUI: the effect of time on soil fine C is lower by 0.001831 xg under Low LUI compared to High LUI. In other words the slope between time and soil fine C is:
((-0.0023700)+(-0.0006606))
# = -0.0030306 xg/d in Low LUI.

# Therefore, Low LUI MRT:
1/((-0.0023700)+(-0.0006606))
# [1] -329.9677 days

g1 <- ggplot(field, aes(x = time, y = soilCf, color = LUI)) + geom_point() +expand_limits(y=0)
g1 + geom_smooth(method = "glm", method.args= list (family = gaussian(link="log"), start=c(5,0)))

############ fine soil C (surface) ################

field.surf$lnsfc <- log(field.surf$soilCf)
fit.ssfc = lm(lnsfc ~ time*LUI, data = field.surf)
summary(fit.ssfc)

# Call:
#   lm(formula = lnsfc ~ time * LUI, data = field.surf)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.36708 -0.18977 -0.00964  0.14249  0.48383 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.0638640  0.0737676  -0.866   0.3924    
# time                  -0.0032364  0.0005762  -5.617 2.26e-06 ***
# LUILow intensity      -0.2376858  0.1043231  -2.278   0.0287 *  
# time:LUILow intensity  0.0001908  0.0008149   0.234   0.8162    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2384 on 36 degrees of freedom
# Multiple R-squared:  0.6541,	Adjusted R-squared:  0.6252 
# F-statistic: 22.69 on 3 and 36 DF,  p-value: 2.017e-08

#MRT for surface fine soil necro-C (high LUI): 
1/-0.0032364
# [1] -308.9853 days

#Low LUI MRT:
1/((-0.0032364)+(0.0001908))
# [1] -328.3425days

gC5 <- ggplot(field.surf, 
              aes(x = time, y = soilCf, color = LUI)) + 
  geom_point() +
  expand_limits(y=0) +
  geom_smooth(method = "glm", 
              method.args= list (family = gaussian(link="log"), 
                                 start=c(5,0)))
gC5
############ fine soil C (deep) ################

field.sub$lnsfc <- log(field.sub$soilCf)
fit.dsfc = lm(lnsfc ~ time*LUI, data = field.sub)
summary(fit.dsfc)

# Call:
#   lm(formula = lnsfc ~ time * LUI, data = field.sub)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.65781 -0.23349 -0.00716  0.19155  0.87099 
# 
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            1.4945101  0.1126025  13.272 1.96e-15 ***
# time                   0.0009062  0.0008796   1.030  0.30977    
# LUILow_intensity       1.1060881  0.1592440   6.946 3.87e-08 ***
# time:LUILow_intensity -0.0038775  0.0012439  -3.117  0.00358 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3639 on 36 degrees of freedom
# Multiple R-squared:  0.6105,	Adjusted R-squared:  0.5781 
# F-statistic: 18.81 on 3 and 36 DF,  p-value: 1.649e-07

#MRT for deep fine soil necro-C (high LUI): 
1/0.0009062
# [1] 1103.509 days (slight increase)

#Low LUI MRT:
1/((0.0009062)+(-0.0038775))
# [1] -336.553 days

gC10 <- ggplot(field.sub, 
               aes(x = time, y = soilCf, color = LUI)) + 
  geom_point() +
  expand_limits(y=0) +
  geom_smooth(method = "glm", 
              method.args= list (family = gaussian(link="log"), 
                                 start=c(5,0)))
gC10
############ fine soil N (total) ################

fieldall$lnsfn <- log(fieldall$soilNf)
fit.sfn = lm(lnsfn ~ time*LUI, data = fieldall)
summary(fit.sfn)
# Call:
#   lm(formula = lnsfn ~ time * LUI, data = fieldall)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.30099 -0.12186 -0.00968  0.07853  0.33199 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -0.4459295  0.0528648  -8.435 4.76e-10 ***
# time                  -0.0002396  0.0004129  -0.580   0.5654    
# LUILow_intensity      -0.0215692  0.0747621  -0.289   0.7746    
# time:LUILow_intensity -0.0011354  0.0005840  -1.944   0.0597 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1708 on 36 degrees of freedom
# Multiple R-squared:  0.3147,	Adjusted R-squared:  0.2575 
# F-statistic: 5.509 on 3 and 36 DF,  p-value: 0.003217

#MRT for fine soil necro-N (high LUI): 
1/-0.0002396
# [1] -4173.623 days 

#Low LUI MRT:
1/((-0.0002396)+(-0.0011354))
# [1] -727.2727 days

g1 <- ggplot(fieldall, aes(x = time, y = soilNf, color = LUI)) + geom_point() +expand_limits(y=0)
g1 + geom_smooth(method = "glm", method.args= list (family = gaussian(link="log"), start=c(5,0)))

############ fine soil N (surface) ################

field.surf$lnsfn <- log(field.surf$soilNf)
fit.ssfn = lm(lnsfn ~ time*LUI, data = field.surf)
summary(fit.ssfn)

# Call:
#   lm(formula = lnsfn ~ time * LUI, data = field.surf)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.3699 -0.1445 -0.0025  0.1273  0.4458 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            2.890e+00  6.528e-02  44.266  < 2e-16 ***
#   time                  -1.121e-03  5.100e-04  -2.199 0.034371 *  
#   LUILow_intensity      -3.800e-01  9.233e-02  -4.115 0.000215 ***
#   time:LUILow_intensity  4.210e-06  7.212e-04   0.006 0.995374    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.211 on 36 degrees of freedom
# Multiple R-squared:  0.5385,	Adjusted R-squared:  0.5001 
# F-statistic:    14 on 3 and 36 DF,  p-value: 3.302e-06

#MRT for surface fine soil necro-N (high LUI): 
1/-1.121e-03
# [1] -892.0607 days 

#Low LUI MRT:
1/((-1.121e-03)+(4.210e-06))
# [1] -895.4235 days

gN5 <- ggplot(field.surf, 
              aes(x = time, y = soilNf, color = LUI)) + 
  geom_point() +
  expand_limits(y=0) +
  geom_smooth(method = "glm", 
              method.args= list (family = gaussian(link="log"), 
                                 start=c(5,0)))
gN5

############ fine soil N (deep) ################

field.sub$lnsfn <- log(field.sub$soilNf)
fit.dsfn = lm(lnsfn ~ time*LUI, data = field.sub)
summary(fit.dsfn)

# Call:
#   lm(formula = lnsfn ~ time * LUI, data = field.sub)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.62175 -0.23726  0.02415  0.15383  0.60532 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.7885639  0.0987816   7.983 1.76e-09 ***
#   time                   0.0036064  0.0007716   4.674 4.05e-05 ***
#   LUILow_intensity       1.1726479  0.1396983   8.394 5.35e-10 ***
#   time:LUILow_intensity -0.0054515  0.0010912  -4.996 1.52e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3192 on 36 degrees of freedom
# Multiple R-squared:  0.6737,	Adjusted R-squared:  0.6465 
# F-statistic: 24.77 on 3 and 36 DF,  p-value: 7.156e-09

#MRT for deep fine soil necro-N (high LUI): 
1/0.0036064
# [1] 277.2848 days ( MRT does not make sense, b/c increases with time)

#Low LUI MRT:
1/((0.0036064)+(-0.0054515))
# [1] -541.976 days

gN10 <- ggplot(field.sub, 
               aes(x = time, y = soilNf, color = LUI)) + 
  geom_point() +
  expand_limits(y=0) +
  geom_smooth(method = "glm", 
              method.args= list (family = gaussian(link="log"), 
                                 start=c(5,0)))
gN10

############ Biomass soil C (total) ################

fieldall$lnmbc <- log(fieldall$MB13C)
fit.mbc = lm(lnmbc ~ time*LUI, data = fieldall)
summary(fit.mbc)

# Call:
#   lm(formula = lnmbc ~ time * LUI, data = fieldall)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.69845 -0.15209  0.03224  0.13206  0.64083 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            4.0283418  0.0939915  42.859  < 2e-16 ***
# time                  -0.0047678  0.0007342  -6.494 1.53e-07 ***
# LUILow_intensity      -0.1818900  0.1329241  -1.368    0.180    
# time:LUILow_intensity -0.0011197  0.0010383  -1.078    0.288    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3037 on 36 degrees of freedom
# Multiple R-squared:  0.7616,	Adjusted R-squared:  0.7418 
# F-statistic: 38.35 on 3 and 36 DF,  p-value: 2.651e-11

#MRT for MBC from necro-C (high LUI): 
1/-0.0047678
# -209.7403 days 

#Low LUI MRT:
1/((-0.0047678)+(-0.0011197))
# -169.8514 days

############ Biomass soil C (surface) ################

field.surf$lnmbc <- log(field.surf$MB13C)
fit.smbc = lm(lnmbc ~ time*LUI, data = field.surf)
summary(fit.smbc)
# Call:
#   lm(formula = lnmbc ~ time * LUI, data = field.surf)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.16062 -0.24373  0.02923  0.23393  1.00410 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            3.8154817  0.1482755  25.732  < 2e-16 ***
#   time                  -0.0056419  0.0011582  -4.871 2.23e-05 ***
#   LUILow_intensity      -0.4510364  0.2096932  -2.151   0.0383 *  
#   time:LUILow_intensity -0.0008508  0.0016380  -0.519   0.6067    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4791 on 36 degrees of freedom
# Multiple R-squared:  0.6512,	Adjusted R-squared:  0.6222 
# F-statistic: 22.41 on 3 and 36 DF,  p-value: 2.331e-08

#MRT for surface MBC from necro-C (high LUI): 
1/-0.0056419
# -177.2453 days 

#Low LUI MRT:
1/((-0.0056419)+(-0.0008508))
# -154.0191 days

############ Biomass soil C (deep) ################

field.sub$lnmbc <- log(field.sub$MB13C)
fit.dmbc = lm(lnmbc ~ time*LUI, data = field.sub)
summary(fit.dmbc)
# Call:
# lm(formula = lnmbc ~ time * LUI, data = field.sub)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.8499 -0.3232 -0.0182  0.2874  0.8954 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            2.211576   0.141030  15.682   <2e-16 ***
# time                  -0.001936   0.001102  -1.757   0.0874 .  
# LUILow_intensity       0.418812   0.199446   2.100   0.0428 *  
# time:LUILow_intensity -0.002629   0.001558  -1.688   0.1001    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4557 on 36 degrees of freedom
# Multiple R-squared:  0.3785,	Adjusted R-squared:  0.3268 
# F-statistic:  7.31 on 3 and 36 DF,  p-value: 0.0005985

#MRT for surface MBC from necro-C (high LUI): 
1/-0.001936
# -516.5289 days 

#Low LUI MRT:
1/((-0.001936)+(-0.002629))
# -219.0581 days

############ Biomass soil N (total) ################

fieldall$lnmbn <- log(fieldall$MB15N)
fit.mbn = lm(lnmbn ~ time*LUI, data = fieldall)
summary(fit.mbn)

# Call:
# lm(formula = lnmbn ~ time * LUI, data = fieldall)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.26494 -0.24809  0.05675  0.33441  0.81902 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            3.4575205  0.1568857  22.038  < 2e-16 ***
#   time                  -0.0042316  0.0012255  -3.453  0.00143 ** 
#   LUILow_intensity      -0.3523606  0.2218698  -1.588  0.12100    
# time:LUILow_intensity -0.0008015  0.0017331  -0.462  0.64654    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.507 on 36 degrees of freedom
# Multiple R-squared:  0.4983,	Adjusted R-squared:  0.4565 
# F-statistic: 11.92 on 3 and 36 DF,  p-value: 1.433e-05

#MRT for MBN from necro-N (high LUI): 
1/-0.0042316
# -236.3172 days 

#Low LUI MRT:
1/((-0.0042316)+(-0.0008015))
# -198.6847 days

############ Biomass soil N (surface) ################

field.surf$lnmbn <- log(field.surf$MB15N)
fit.smbn = lm(lnmbn ~ time*LUI, data = field.surf)
summary(fit.smbn)
# Call:
# lm(formula = lnmbn ~ time * LUI, data = field.surf)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.50771 -0.31456  0.04845  0.41751  1.08801 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            3.4110265  0.1818589  18.756  < 2e-16 ***
#   time                  -0.0054272  0.0014205  -3.821 0.000507 ***
#   LUILow_intensity      -0.5570985  0.2571874  -2.166 0.037002 *  
#   time:LUILow_intensity -0.0007365  0.0020089  -0.367 0.716073    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5877 on 36 degrees of freedom
# Multiple R-squared:  0.5535,	Adjusted R-squared:  0.5163 
# F-statistic: 14.88 on 3 and 36 DF,  p-value: 1.843e-06

#MRT for surface MBN from necroN (high LUI): 
1/-0.0054272
# -184.2571 days 

#Low LUI MRT:
1/((-0.0054272)+(-0.0007365))
# -162.2402 days

############ Biomass soil N (deep) ################

field.sub$lnmbn <- log(field.sub$MB15N)
library(dplyr)
field.sub.2 <- filter(field.sub, lnmbn != "-Inf")
fit.dmbn = lm(lnmbn ~ time*LUI, data = field.sub.2)
summary(fit.dmbn)
# Call:
# lm(formula = lnmbn ~ time * LUI, data = field.sub.2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.04288 -0.18887 -0.01127  0.32362  1.59725 
# 
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.8174854  0.2022894   4.041 0.000326 ***
# time                   0.0009165  0.0014578   0.629 0.534151    
# LUILow_intensity       0.7777639  0.2779943   2.798 0.008766 ** 
# time:LUILow_intensity -0.0045240  0.0020309  -2.228 0.033304 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.562 on 31 degrees of freedom
# Multiple R-squared:  0.2437,	Adjusted R-squared:  0.1705 
# F-statistic:  3.33 on 3 and 31 DF,  p-value: 0.03216


#MRT for deep MBN from necro-N (high LUI): 
1/0.0009165
# 1091.107 days (increase)

#Low LUI MRT:
1/((0.0009165)+(-0.0045240))
# -277.2003 days  




###############################################################################
# microbial community analysis - 16S and ITS

otu <- read.csv("otu_Myerscough.csv", header=TRUE, sep=",", row.names = 1) ##this keeps a header row, empty cells as a comma, and OTU numbers as row names
otuits <- read.csv("otu_ITS_Myerscough.csv", header=TRUE, sep=",", row.names = 1)

library(data.table)

### 16S
#getrid of all of the all-zero rows
otu1 <- subset(otu, sum != 0)
otu1 <- otu1[-grep('sum',colnames(otu1))] #remove sum column
head(otu1)

tax <- otu1$taxonomy # makes a taxonomy only object

otu2 <- otu1[,1:10] # makes a file without taxonomy column
head(otu2, 3)

### ITS
#getrid of all of the all-zero rows
otuits1 <- subset(otuits, sum != 0)
otuits1 <- otuits1[-grep('sum',colnames(otuits1))]
head(otuits1)

taxits <- otuits1$taxonomy # makes a taxonomy only object

otuits2 <- otuits1[,1:10] # makes a file without taxonomy column
head(otuits2, 3)

##transpose tables
otu3 <- t(otu2)
otuits3 <- as.matrix(t(otuits2))

### taxonomy
# write a function to extract name, from: http://statweb.stanford.edu/~susan/Summer11/Labs/Lab1ngs_data_manipulation.pdf
extract.name.level= function(x, level){
  a=c(unlist(strsplit(x,';')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}
#a = as.character(otu$taxonomy[1])
a = as.character(otuits$taxonomy[1])
a
extract.name.level(a, 2)
extract.name.level(a, 3)
extract.name.level(a, 5)

# same source, write a function to summarize data at different taxonomic levels
otu3taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa)){
    taxa=colnames(x)
  }
  if(length(taxa)!=dim(x)[2]){
    print("ERROR: taxonomy should have the same length as the number of columns in OTU table")
    return;
  }
  level.names=sapply(as.character(taxa), 
                     function(x)
                       extract.name.level(x,level=level))
  t(apply(x, 1, 
          function(y) 
            tapply(y,level.names,sum)))
}

d.phylum = otu3taxonomy(otu3,level=2,taxa=tax)
d.class = otu3taxonomy(otu3,level=3,taxa=tax)
d.order = otu3taxonomy(otu3,level=4,taxa=tax)
d.family = otu3taxonomy(otu3,level=5,taxa=tax)
d.genus = otu3taxonomy(otu3,level=6,taxa=tax)
d.species = otu3taxonomy(otu3,level=7,taxa=tax)

d.phylumits = otu3taxonomy(otuits3,level=2,taxa=taxits)
d.classits = otu3taxonomy(otuits3,level=3,taxa=taxits)
d.orderits = otu3taxonomy(otuits3,level=4,taxa=taxits)
d.familyits = otu3taxonomy(otuits3,level=5,taxa=taxits)
d.genusits = otu3taxonomy(otuits3,level=6,taxa=taxits)
d.speciesits = otu3taxonomy(otuits3,level=7,taxa=taxits)

# summarize the OTU table at different taxonomic levels

### 16S

phylum2 <-t(d.phylum)
class2 <-t(d.class)
order2 <-t(d.order)
family2 <-t(d.family)
genus2 <-t(d.genus)
species2 <-t(d.species)
write.table(phylum2, file="output/phyla.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE) 
write.table(class2, file="output/classes.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(order2, file="output/orders.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(family2, file="output/families.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(genus2, file="output/genera.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(species2, file="output/species.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)

### ITS

phylumits2 <-t(d.phylumits)
classits2 <-t(d.classits)
orderits2 <-t(d.orderits)
familyits2 <-t(d.familyits)
genusits2 <-t(d.genusits)
speciesits2 <-t(d.speciesits)
write.table(phylumits2, file="output/phylaits.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE) 
write.table(classits2, file="output/classesits.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(orderits2, file="output/ordersits.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(familyits2, file="output/familiesits.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(genusits2, file="output/generaits.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(speciesits2, file="output/speciesits.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)

#### beta diversity

# vegan
library(vegan)
library(permute)
library(mgcv)

meta <- read.csv("metaCN.csv", header=TRUE, row.names=1, sep=",", na.strings=c("na", "NA"))

## making a distance matrix for pcoa ####

vegotu <- vegdist(otu3, binary=FALSE)#Bray-Curtis matrix, use this, Binary performs p/a standardization first
vegotuITS <- vegdist(otuits3, binary=FALSE)

### pcoa 
### This method is also known as MDS (Metric Multidimensional Scaling)
### PCoA provides Euclidean representation of a set of objects whose relationship is measured by any similarity or distance measure chosen by the user
### PCoA returns a set of orthogonal axes whose importance is measured by eigenvalues. This means that calculating PCoA on Euclidean distances among samples yields the same results as PCA calculated on covariance matrix of the same dataset (if scaling 1 is used). 

pcoa <- cmdscale(vegotu)
pcoa
pcoaits <- cmdscale(vegotuITS)
pcoaits

library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
# install.packages('devtools')
# library(devtools)
# install_github('fawda123/ggord')
# library(ggord)

#### PCOA figure ########

df <- data.frame(PCOA1 = pcoa[,1], PCOA2 = pcoa[,2])
p1 <- ggplot(data = df, aes(PCOA1, PCOA2)) + 
  geom_point(shape = 19,
             size=4, 
             aes(color = meta$LUI)) +
  scale_color_manual(values = c("#253494","#41b6c4"),
                    name="Land \nmanagement") +
  geom_text(x = -0.15, y = -0.35,
            label = "16S \nP=0.0076",
            size=5) +
  theme_bw() +
  theme(axis.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black")) 
p1

dfits <- data.frame(PCOA1 = pcoaits[,1], PCOA2 = pcoaits[,2])
p2 <- ggplot(data = dfits, aes(PCOA1, PCOA2)) + 
  geom_point(shape = 19,
             size=4, 
             aes(color = meta$LUI)) +
  scale_color_manual(values = c("#253494","#41b6c4"),
                     name="Land \nmanagement") +
  geom_text(x = -0.2, y = -0.35,
            label = "ITS \nP=0.0083",
            size=5) +
  theme_bw() +
  theme(axis.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black")) 
p2

x <- ggarrange(p1, p2, 
                  ncol = 2, nrow=1, common.legend=TRUE, legend="right")
png(file="CNpcoa.png", width=240, height=100, units='mm', res=600)
x
dev.off() 

#### hypothesis testing - PERMANOVA

ano = anosim(vegotu, meta$LUI, distance = "bray", permutations = 9999)
ano
# Call:
#   anosim(x = vegotu, grouping = meta$LUI, permutations = 9999,      distance = "bray") 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.556 
# Significance: 0.0076 
# 
# Permutation: free
# Number of permutations: 9999

anoits = anosim(vegotuITS, meta$LUI, distance = "bray", permutations = 9999)
anoits
# Call:
#   anosim(x = vegotuITS, grouping = meta$LUI, permutations = 9999,      distance = "bray") 
# Dissimilarity: bray 
# 
# ANOSIM statistic R:  0.94 
# Significance: 0.0083 
# 
# Permutation: free
# Number of permutations: 9999

#note adonis is for multivariate
#adonis (vegotu ~ LUI, meta)
