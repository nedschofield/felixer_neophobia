#reproducible R code

#Evasive invasives? Implications of neophobia for feral cat (Felis catus) control.

#Authors
#Ned L. Ryan-Schofield *, John L. Read, Hugh W. McGregor, Todd J. McWhorter, Katherine E. Moseby
# *corresponding author nedscho555@gmail.com

#load packages
library(dplyr)
library(MuMIn)
library(glmmTMB)
library(lmtest)
library(ggplot2)
library(ggiraph)
library(ggiraphExtra)


##########################################################################
######### inside vs outside Arid Recovery experiment analysis ############
##########################################################################

#download data off of github https://github.com/nedschofield/felixer_neophobia
#set wd to folder containing "cat_neophobia_data_collated_v4.csv"

##import data
ind.data.all <- read.csv("~/University/2020/Masters/April Arid recovery trial/Felixer experiment/cat_neophobia_data_collated_v4.csv")


#subset to just columns of interest
catpasses <-subset(ind.data.all, select = c(1:37))

#create new variable passage for all the types of passes that could happen
catpasses <- catpasses %>%
  mutate(passage = if_else(wing1_detected == 0 & centre_detected == 0 & wing2_detected == 0, 000, 
                           if_else(wing1_detected == 0 & centre_detected == 0 & wing2_detected == 1, 001,
                                   if_else(wing1_detected == 0 & centre_detected == 1 & wing2_detected == 1, 011,
                                           if_else(wing1_detected == 1 & centre_detected == 1 & wing2_detected == 1, 111,
                                                   if_else(wing1_detected == 0 & centre_detected == 1 & wing2_detected == 0, 010,
                                                           if_else(wing1_detected == 1 & centre_detected == 0 & wing2_detected == 0, 100,
                                                                   if_else(wing1_detected == 1 & centre_detected == 0 & wing2_detected == 1, 101,
                                                                           if_else(wing1_detected == 1 & centre_detected == 1 & wing2_detected == 0, 110, 2)))))))))

#subset to passes that were detected on 2 or more cameras. i.e. those that were walking along the fence for 20+ metres
catpasses <- subset(catpasses, passage == 111 | passage == 101 | passage == 110 | passage == 011)


####################################
####### summary statistics #########
####################################

# number of passes in experimental paddock and control area
all.pass.exp <- nrow(subset(ind.data.all, location == "dingo"))
all.pass.exp
all.pass.cont <- nrow(subset(ind.data.all, location == "wildWest"))
all.pass.cont

#spread of data across sites
sites <- ind.data.all %>% 
  group_by(site) %>%
  summarise(n = n()) %>%
  add_row(site = "CR03", n = 0)
mean(sites$n)
min(sites$n)
max(sites$n)
sd(sites$n)

#number of passes detected on 2 or more cameras
two.cam.all <- nrow(catpasses)
two.cam.all
two.cam.exp <- nrow(subset(catpasses, location == "dingo"))
two.cam.exp
two.cam.cont <- nrow(subset(catpasses, location == "wildWest"))
two.cam.cont

#how many cats avoided cameras for all passes detected on 2 or more cameras
avoided.two.cam.all <- nrow(subset(catpasses, passage == "101"))
avoided.two.cam.all
#percentage of passes avoided
(avoided.two.cam.all/two.cam.all)*100

#avoidance in experimental paddock over first 37 days
#number of passes
exp.37.all <- nrow(subset(catpasses, location == "dingo" & days_since_site_set <= 37))
exp.37.all
#number avoided 
exp.37.avoid <- nrow(subset(catpasses, location == "dingo" & passage == "101" & days_since_site_set <= 37)) 
exp.37.avoid
#percentage of passes avoided
(exp.37.avoid/exp.37.all)*100

#avoidance in control area over first 37 days
#number of passes
cont.37.all <- nrow(subset(catpasses, location == "wildWest" & days_since_site_set <= 37))
cont.37.all
#number avoided 
cont.37.avoid <- nrow(subset(catpasses, location == "wildWest" & passage == "101" & days_since_site_set <= 37)) 
cont.37.avoid
#percentage of passes avoided
(cont.37.avoid/cont.37.all)*100


######################################
##### model avoidance behaviour ######
######################################

##create models
m0 <- glm(avoided ~ 1, family = binomial, catpasses)
m1 <- glm(avoided ~ treatment, family = binomial, catpasses)
m2 <- glm(avoided ~ days_since_site_set, family = binomial, catpasses)
m3 <- glm(avoided ~ days_since_site_set + treatment, family = binomial, catpasses)
m4 <- glm(avoided ~ days_since_site_set * treatment, family = binomial, catpasses)

#include models with site type (crate vs felixer)
m0.1 <- glm(avoided ~ site_type, family = binomial, catpasses)
m1.1 <- glm(avoided ~ treatment + site_type, family = binomial, catpasses)
m2.1 <- glm(avoided ~ days_since_site_set + site_type, family = binomial, catpasses)
m3.1 <- glm(avoided ~ days_since_site_set + treatment + site_type, family = binomial, catpasses)
m4.1 <- glm(avoided ~ days_since_site_set * treatment * site_type, family = binomial, catpasses)


#compare in an information theory framework
model.sel(m0, m1, m2, m3, m4, m0.1, m1.1, m2.1, m3.1, m4.1)

#likelihood ratio test shows no added value from interaction term (m4) or site_type (m3.1) which are both not significant
#logLik is pretty much the same so dismiss based on Arnold (2010)
lrtest(m3, m4)
lrtest(m3, m3.1)

#get statistics for the top candidate model
summary(m3)
#get odds ratio and conf int for top model
exp(coef(m3))
exp(cbind(OR = coef(m3), confint(m3)))


##visualise predicted effects

#create predicted effects with ggPredict
pred.m3 <- ggPredict(m3, se=T, jitter = F)

#plot predicted effects (colour blind friendly)
plot(pred.m3) + labs(x="Days since Felixer set", 
                     y="Probability of avoiding Felixer")+
  theme_bw()+
  theme(legend.position = c(0.8, 0.8), panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_viridis_d()+
  scale_colour_viridis_d()




#################################################################
###### analysis of individual variation in collared cats ########
#################################################################

#subset to just cats that were individually identified
ind.data <-subset(catpasses, ID_dataset == "catID")

#number of passes, and number of avoidances
nrow(ind.data)
nrow(subset(ind.data, passage == "101"))


#explore using MuMin for model selection and averaging
options(na.action = "na.fail") # Required for dredge to run

#create a saturated model (needs to be realistic with maximum of 31 non-fixed predictors - not to many interaction terms)
sat <- glmmTMB(formula = walkpast ~ days_cat_exposed + days_cat_in  + captured_no + days_since_last_capture + weight_g + sex + (1 | catID), family = binomial, data = ind.data)

#dredge all nested models of saturated model
sat.dredged <- dredge(sat, beta = F, evaluate = T, rank = AICc)

options(na.action = "na.omit") # set back to default


#rank models by AICc
model.sel(sat.dredged)

#get top candidate models and average parameters
summary(model.avg(sat.dredged, subset = delta <= 2))

#get importance and number of models
sw(model.avg(sat.dredged, subset = delta <=2))
nrow(sat.dredged)

#test null model against null model with catID as random effect
no <- glmmTMB(avoided ~ 1, family = binomial, ind.data)
no1 <- glmmTMB(avoided ~ 1 + (1|catID), family = binomial, ind.data)
lrtest(no, no1)

