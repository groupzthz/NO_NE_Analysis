#### Title: Analysis - Novel object test - tracking system data ####
# compare the day before, during and after the NO test #
# version 5 - final version #

rm(list=ls())
dev.off()


#### packages ####
library("data.table")
library("tidyverse")
library("lme4")
library("DHARMa")
library("effects")
library("plotly")
library("MuMIn")
library("performance")
library("car")
library("pbkrtest")


library("lmerTest")
library("broom")
library("multcomp")
library("broom.mixed")
library("lemon")
library("gridExtra")




#### LOAD DATA #####################################################################

# set working directory
setwd("//nas-vetsuisse.campus.unibe.ch/VETSUISSE/Benutzer/candelotto/E3_laying hens_NO_NE/6_Statistical_Evaluation/Analysis_tracking/NO_tracking_analysis")

### load data
Litter <- read.table ('NO_trackingVariables_Litter.csv', header= TRUE, sep= ',')



#### Preparations ####

### format column data type
str(Litter)
Litter$Event <- as.factor(Litter$Event)
Litter$Pen <- as.factor(Litter$Pen)
Litter$BirdID <- as.factor(Litter$BirdID)
Litter$Date <- as.Date(Litter$Date)

### sort Events
Litter$Event <- relevel(Litter$Event, "post-test")
Litter$Event <- relevel(Litter$Event, "test")
Litter$Event <- relevel(Litter$Event, "pre-test")





####_______________________#########################
#### TRANSITIONS INTO ZONE #################################

### quick plot
plot(Litter$Event, Litter$intoZoneTrans, ylab = "transitions")
# individual variation in transitions
plot(Litter$BirdID, Litter$intoZoneTrans)



### RI MODEL ####################################################
# model with random intercept


### initial model 
intoTrans.lmer <- lmer(intoZoneTrans ~ Event +
                         (1|Pen/BirdID) + (1|day), 
                       data = Litter, 
                       REML= FALSE)
summary(intoTrans.lmer)



#### Residuals 
qqnorm (resid (intoTrans.lmer))
qqline (resid (intoTrans.lmer))
qqnorm (ranef (intoTrans.lmer) [['BirdID:Pen']] [, 1])
scatter.smooth (fitted (intoTrans.lmer), resid (intoTrans.lmer))

resid.test<- simulateResiduals(intoTrans.lmer, 1000)
plot(resid.test) 
plotResiduals(resid.test, form = Litter$Event) #good
testDispersion(intoTrans.lmer)
simulationOutput.fin <- simulateResiduals(fittedModel = intoTrans.lmer, plot = F,quantreg=T)
plot(simulationOutput.fin)



#### -> Group effects ######################
summary(intoTrans.lmer)
anova(intoTrans.lmer)

## estimates
as.data.frame(effect("Event",intoTrans.lmer))

## model conditional and marginal r^2
r.squaredGLMM(intoTrans.lmer)



#### -> Repeatability #########################
# (Variation explained by individual differences in behavioural type)
print(VarCorr(intoTrans.lmer),comp = c("Variance","Std.Dev."))

VarCorr(intoTrans.lmer)$"BirdID:Pen"[1] /
  (VarCorr(intoTrans.lmer)$"BirdID:Pen"[1] +
     VarCorr(intoTrans.lmer)$"day"[1] +
     VarCorr(intoTrans.lmer)$"Pen"[1] +
     attr(VarCorr(intoTrans.lmer), "sc")^2)






#### RIS MODEL ####################################################
# model with random intercept and random slope

### ris model
intoTrans.lmer.ris <- lmer(intoZoneTrans ~ Event +
                         (Event|BirdID), 
                       data = Litter, 
                       REML= FALSE)
summary(intoTrans.lmer.ris)
# full model with random slope is overfitting, 
# when removing day and Pen, the model is still overfitting 
# but the Std Error of fixed factor looks better


#### model comparison according to AIC
anova(intoTrans.lmer, intoTrans.lmer.ris)
# Δ = 24.3
# -> RI model is better 







####_______________________________#########################
#### TOTAL ZONE DURATION #################################

### quick plots
plot(Litter$Event, Litter$totalZoneDur/60/60, ylab = "tot. duration Litter [h]")
# individual variation in total duration
plot(Litter$BirdID, Litter$totalZoneDur/60/60)




### RI MODEL ##########################################
totDur.lmer <- lmer(totalZoneDur ~ Event +
                      #  (1|Pen/BirdID) + (1|day), 
                      (1|BirdID),
                    data = Litter, 
                    REML= FALSE)
summary(totDur.lmer)
# overfitting, because day and pen don't explain any of the variance 
#             therefore remove Pen and day



#### Residuals 
qqnorm (resid (totDur.lmer))
qqline (resid (totDur.lmer))
qqnorm (ranef (totDur.lmer) [['BirdID']] [, 1])
scatter.smooth (fitted (totDur.lmer), resid (totDur.lmer))

resid.test<- simulateResiduals(totDur.lmer, 1000)
plot(resid.test) 
plotResiduals(resid.test, form = Litter$Event) 
testDispersion(totDur.lmer) 




#### -> Group effects ####
anova(totDur.lmer)

## model conditional and marginal r^2
r.squaredGLMM(totDur.lmer)


## effect size 
effectSize <- as.data.frame(effect("Event",totDur.lmer))
# change from sec into min
sec_to_min <- function(x){
  return ((x)/60) }
# applying the custom function to every value except for first col
cbind(Event = effectSize$Event, data.frame(lapply(effectSize[,-1],sec_to_min)))
# difference preTest - test = 19.77min



#### -> Repeatability ####
# (Variation explained by individual differences in behavioural type)
print(VarCorr(totDur.lmer),comp = c("Variance","Std.Dev."))

VarCorr(totDur.lmer)$"BirdID"[1] /
  (VarCorr(totDur.lmer)$"BirdID"[1] +
     attr(VarCorr(totDur.lmer), "sc")^2)





#### RIS MODEL ################################################

totDur.lmer.ris <- lmer(totalZoneDur ~ Event +
                          (Event|BirdID),
                        data = Litter, 
                        REML= FALSE)
summary(totDur.lmer.ris)
# singular fit
# slope doesn't seem to explain much of the variance
# but standard error is ok


#### model comparison
anova(totDur.lmer, totDur.lmer.ris)
# model without random slope better




#### PLOT MODEL PREDICTIONS #########################
#### plot RI model predictions 


# get predicted values
RI <-augment(totDur.lmer)%>%
  dplyr::select(totalZoneDur,Event,BirdID,.fitted) %>%
  dplyr::group_by(BirdID,Event) %>%
  dplyr::summarise(Pred=mean(.fitted))


# sort event
RI$Event <- factor(RI$Event, levels = c("pre-test", "test", "post-test"))


# use plotly to make plot smoother
plotly::ggplotly(
  ggplot(RI,aes(x=Event,y=(Pred/60/60), group = BirdID))+
    geom_line(size=0.6, col = "darkgray") +
    theme_classic() +
    labs(y="", x="") +
    ylab("Duration in Litter [h]") +
    xlab ("") +
    # scale_y_continuous(
    #   breaks = c(-400, -200, 0, 200)) +
    theme(text = element_text(size=12),
          axis.text.x = element_text(size=14)) +
    theme(legend.position="none")
)










####_______________________________##################################
#### LATENCY TO FIRST LITTER VISIT #################################



#### DATA PREPARATIONS ####################################
### special cases:
## NA's 
# Birds that never went into the litter on the particular day
# giving them the maximum recording time (15h = 54’000s) would introduce outliers into the dataset
# Whether birds go into the litter or not
#        is accounted for in the other models i.e. duration of 0, transitions of 0
#        -> therefore exclude tests in which the bird didn't visit the litter for latency analysis 
### Zeros
# Birds that were already present in the litter when the recording started at 02:00
# In those cases it would probably be best to check each of them individually
# and adjust the tracking data based on whether they moved in just before 2:00 
# or based on whether it was only flickering or whether they actually slept in the litter at least for NE birds
# for now, exclude tests with zeros


# exclude all tests with at least one day 
# of a zero or a NA (to later analyze only complete tests)

#### remove zero rows and NA rows
Litter_adj <- subset(Litter, firstVisit != 0 &
                       !is.na(firstVisit))
# remove the full test round for all birds that now have less than 3 rows per test round
# in order to later use complete cases.
Litter_adj <- Litter_adj %>%
  group_by(BirdID, Test_count) %>%
  filter(n() == 3)
Litter_adj$BirdID <- droplevels(Litter_adj$BirdID)
# check if all went well
table(Litter_adj$BirdID, Litter_adj$Test_count)
unique(Litter_adj$BirdID)
# reduces the number of birds to 101







#### RI MODEL #######################################################

### quick plots
plot(Litter_adj$Event, (Litter_adj$firstVisit/60/60), ylab = "Latency [h]")
plot(Litter_adj$Test_count:Litter_adj$Event, Litter_adj$firstVisit)



#### -> model selection #########

### initial model 
latency.lmer <- lmer(firstVisit ~ Event +
                       #(1|Pen/BirdID) + (1|day), 
                       (1|BirdID),
                     data = Litter_adj, 
                     REML= FALSE)
summary(latency.lmer)
## singular fit, because day and Pen both don't explain any of the variance
# therefore remove day and Pen


#### Residuals 
qqnorm (resid (latency.lmer))
qqline (resid (latency.lmer)) 
qqnorm (ranef (latency.lmer) [['BirdID:Pen']] [, 1])
scatter.smooth (fitted (latency.lmer), resid (latency.lmer))

resid.test<- simulateResiduals(latency.lmer, 1000)
plot(resid.test) 
plotResiduals(resid.test, form = Litter_adj$Event) 
testDispersion(latency.lmer) 
simulationOutput.fin <- simulateResiduals(fittedModel = latency.lmer, plot = F,quantreg=T)
plot(simulationOutput.fin) 



#### try glmer 

# with poisson
latency.glmer.pois <- glmer(firstVisit ~ Event +
                            #(1|Pen/BirdID) + (1|day), 
                            (1|BirdID) + (1|day), 
                            data = Litter_adj,
                       family = "poisson")
summary(latency.glmer.pois)
# singular fit, pen explains almost nothing (^-9), 
# excluding pen works

### Residuals
scatter.smooth (fitted (latency.glmer.pois), resid (latency.glmer.pois))
qqnorm (resid (latency.glmer.pois)); qqline (resid (latency.glmer.pois)) 
testDispersion(latency.glmer.pois) 
check_overdispersion(latency.glmer.pois)


# with negative binomial
latency.glmer.nb <- glmer.nb(firstVisit ~ Event +
                              #(1|Pen/BirdID) + (1|day), 
                              (1|BirdID) + (1|day),
                          data = Litter_adj)
summary(latency.glmer.nb)
# singular fit, pen explains almost nothing (^-10), 
#  excluding pen works

### Residuals 
scatter.smooth (fitted (latency.glmer.nb), resid (latency.glmer.nb))
qqnorm (resid (latency.glmer.nb)); qqline (resid (latency.glmer.nb)) 
testDispersion(latency.glmer.nb) 
check_overdispersion(latency.glmer.nb)






#### compare the models with each other
anova(latency.glmer.pois, latency.glmer.nb) # nb ist besser als poisson
anova(latency.glmer.nb, latency.lmer)
anova(latency.glmer.pois,latency.lmer)

# model conditional and marginal r^2
r.squaredGLMM(latency.lmer)
r.squaredGLMM(latency.glmer.pois)
r.squaredGLMM(latency.glmer.nb)






#### -> final model ######################

latency.glmer.nb <- glmer.nb(firstVisit ~ Event +
                               #(1|Pen/BirdID) + (1|day), 
                               (1|BirdID) + (1|day),
                             data = Litter_adj)



#### -> group level results ##############
summary(latency.glmer.nb)
anova(latency.glmer.nb)


#### R^2 
# Pseudo-R-squared
# https://search.r-project.org/CRAN/refmans/MuMIn/html/r.squaredGLMM.html
# preferably use the trigamma approximation
r.squaredGLMM(latency.glmer.nb)
## with different package (performance) Nakagawa's R2
# https://easystats.github.io/performance/reference/r2_nakagawa.html
r2_nakagawa(latency.glmer.nb)

 

#### p - value

### Wald-Z approximation, though not very good
Anova(latency.glmer.nb, type=3) 

### Parametric bootstrap
# https://www.rdocumentation.org/packages/pbkrtest/versions/0.3-5/topics/PBmodcomp
# create null model (since there is only one fixed factor)
latency.glmer.nb.NULL <- update(latency.glmer.nb, .~.-Event)
# anova model comparison
anova(latency.glmer.nb, latency.glmer.nb.NULL) # AICΔ = 7, p=0.004
# parametric bootstrap approach
PBmodcomp(latency.glmer.nb, latency.glmer.nb.NULL, nsim=1000)





### calculate backtransformed estimates
# first adding the coefficients then backtransforming by using the exponent.
estimates <- summary(latency.glmer.nb)$coefficients[,1]
# pre-test in min
exp(estimates[1]) /60
# test day in min
exp(estimates[1] + estimates[2]) /60
# post test day
exp(estimates[1] + estimates[3]) /60

## effect size 
effectSize <- as.data.frame(effect("Event",latency.glmer.nb))
# back transform and change from sec into min
secTOmin <- function(x){
  return (x/60)}
cbind(Event = effectSize$Event, data.frame(lapply(effectSize[,-1],secTOmin)))
# (same as calculated above by hand so seems to be fine)




#### -> Repeatability ######################

print(VarCorr(latency.glmer.nb),comp = c("Variance","Std.Dev."))
# negative binomial doesn't estimate residual variance
# because no gaussian distribution is assumed and thus 
# mean and variance of those distributions are defined by the same parameter.
# see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2016q2/024823.html
#     and https://stats.stackexchange.com/questions/153611/interpreting-random-effect-variance-in-glmer



## use icc from performance package
# https://easystats.github.io/performance/reference/icc.html
icc(latency.glmer.nb, by_group = TRUE)






#### RIS MODEL ###########################################

### final ri model
ri <- glmer.nb(firstVisit ~ Event +
                 (1|BirdID) + (1|day),
               data = Litter_adj)



### ris model
ris <- glmer.nb(firstVisit ~ Event +
                  (Event|BirdID) + (1|day),
                data = Litter_adj)
summary(ris)
# singular fit but removing day will still produce a singular fit




#### PLOT MODEL PREDICTIONS ####
# plot RI vs. RIS


### extract predictions for ri and ris
RI <-augment(ri)%>%
  dplyr::select(firstVisit,Event,BirdID,.fitted)

RIS <-augment(ris)%>%
  dplyr::select(firstVisit,Event,BirdID,.fitted)

RI$.fittedRIS <-RIS$.fitted

df <- RI %>%
  dplyr::group_by(BirdID,Event) %>%
  dplyr::summarise(RI=mean(.fitted),
                   RIS=mean(.fittedRIS))%>%
  gather(type,Value,"RI":"RIS")



## highlight interesting birds
df$animal_id <- as.character(df$BirdID)

df$ID <- ifelse(df$animal_id %in%
                  c("10_sp","5_wg", "3_sb", "12_sws"), #4ws auch gut #4_wp  4_gg
                df$animal_id,"Other individuals")


### sort event
df$Event <- factor(df$Event, levels = c("pre-test", "test", "post-test"))



### prepare ri and ris plots
# use plotly to make plot smoother


plot_ri <- plotly::ggplotly(
  ggplot(df[df$type =="RI",],
         # backtransform predictions
         aes(x=Event,y=(exp(Value)/60/60),group=animal_id,color=relevel(as.factor(ID), "Other individuals")))+
    geom_line(size=0.8) +theme_classic() +
    labs(y="latency [h]", x="") +
    scale_color_manual(values = c("gray", "blue", "darkorange" , "red" ,"#7CAE00")) + ## #FFCC00","#00BFC4","red","gray")) +
    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    ## for a closer look at lower values
    #coord_cartesian(ylim=c(0,1)) +
    theme(text = element_text(size=12),
          axis.text.x = element_text(size=14)) +
    theme(legend.position="none")
)


plot_ris <- plotly::ggplotly(
  ggplot(df[df$type =="RIS",],
         # backtransform predictions
         aes(x=Event,y=(exp(Value)/60/60),group=animal_id,color=relevel(as.factor(ID), "Other individuals")))+
         geom_line(size=0.8) +
    theme_classic() +
    labs(y="latency [h]", x="") +
    scale_color_manual(values = c("gray", "blue", "darkorange" , "red" ,"#7CAE00"))+ #FFCC00","#00BFC4","red","gray")) +
    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    ## for a closer look at lower values
    #coord_cartesian(ylim=c(0,1)) +
    theme(text = element_text(size=12),
          axis.text.x = element_text(size=14)) +
    theme(legend.position="none")
)


### plot without plotly
# mylegend <- g_legend(plot_ri)
# grid.arrange(arrangeGrob(plot_ri + theme(legend.position="none"),
#                          plot_ris + theme(legend.position="none"),
#                          ncol=2, mylegend) )


### plot with plotly
subplot(plot_ri, plot_ris, nrows=1)







