library(readr)
library(ggplot2)
library(dplyr)
library(tibble)
library(ipw)
library(survey)
library(ggthemes)
library(lme4)

# Preparing the preDIVA data set for analysis
{
# Loading the data
data = read_csv("dataAlc.csv")

# Create T variables (denoted as DT)
data["DT12"] = as.numeric(as.Date(data$VisitDt_V2, format = "%m/%d/%Y") - as.Date(data$RandomisatieDatum, format = "%Y/%m/%d"))
data["DT23"] = as.numeric(as.Date(data$VisitDt_V3, format = "%m/%d/%Y") - as.Date(data$VisitDt_V2, format = "%m/%d/%Y"))

# Constructing adherence per year
data["Adherence12"] = data$aantalbezoekenvaatgroep_V2_2/data$DT12*365.25
data["Adherence23"] = data$aantalbezoekenvaatgroep_V3_2/data$DT23*365.25

# Binary adherence (Border is set as attending at least 66% of the scheduled visits (2 per year))
data["Adherence12yn"] = as.numeric(data$Adherence12 >=2)
data["Adherence23yn"] = as.numeric(data$Adherence23 >=2)

# Outcome variable
data["WHOchange"] = data$WHOpa_V3-data$WHOpa_V1

# Create final dataset
t = na.omit(data[,c("VisitDt_V1","VisitDt_V2","VisitDt_V3","Randomisation","Adherence12yn","Adherence23yn","AgeVisit1","Gender","WHOchange")])
nrow(t[t$Randomisation==1,])
nomis = na.omit(data[,c("StudieNr","Adherence12yn","AgeVisit1","Gender","DT12","Adherence23yn","DT23","WHOchange","GPNr","Randomisation")])
nrow(nomis) # Print amount of data points
predivdat = nomis[nomis["Randomisation"] == 1,]
}
# Plot time between baseline and first follow-up
{
ggplot(data,aes(x=DT12)) + geom_histogram(binwidth=30,fill='white',color='black') +theme_minimal() +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 0, l = 0))) +
  labs(
    x = "Days between randomization and first follow-up measurement",
    y = "Frequency"  )

{
hist(data$DT12,breaks=100)
}
}
#Inspect association T and Adherence
{
predivdat["GPNr"] = as.factor(predivdat$GPNr)
predivdat["DT12norm"] = (predivdat$DT12- mean(predivdat$DT12))/sd(predivdat$DT12)
m = glmer(Adherence23yn ~ DT12norm + AgeVisit1 + Gender +(1| GPNr), family = "binomial", data = predivdat)
summary(m)
}
# IPTW calculation
{
PropModel12 = ipwpoint(exposure = Adherence12yn,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ 1,
                       denominator = ~ AgeVisit1+Gender+GPNr,
                       trunc = 0.0125, data = as.data.frame(predivdat))
Weights12 = PropModel12$ipw.weights

ipwplot(weights = PropModel12$ipw.weights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(PropModel12$weights.trunc)))

PropModel23 = ipwpoint(exposure = Adherence23yn,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ Adherence12yn,
                       denominator = ~ Adherence12yn + AgeVisit1 + Gender + GPNr,
                       trunc = 0.0125, data = as.data.frame(predivdat))
Weights23 = PropModel23$ipw.weights

ipwplot(weights = PropModel23$ipw.weights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(PropModel12$weights.trunc)))
IPTWweights = PropModel12$ipw.weights*PropModel23$weights.trunc
ipwplot(weights = IPTWweights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(IPTWweights)))
}
# TAC calculation
{
library(ipw)
PropModel12 = ipwpoint(exposure = Adherence12yn,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ 1,
                       denominator = ~ AgeVisit1+Gender+GPNr,
                       trunc = 0.0125, data = as.data.frame(predivdat))
Weights12 = PropModel12$ipw.weights

ipwplot(weights = PropModel12$ipw.weights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(PropModel12$weights.trunc)))

PropModel23 = ipwpoint(exposure = Adherence23yn,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ Adherence12yn,
                       denominator = ~ Adherence12yn + AgeVisit1 + Gender + GPNr + DT12,
                       trunc = 0.0125, data = as.data.frame(predivdat))
Weights23 = PropModel23$ipw.weights

ipwplot(weights = PropModel23$ipw.weights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(PropModel12$weights.trunc)))
TACweights = PropModel12$ipw.weights*PropModel23$weights.trunc
ipwplot(weights = TACweights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(TACweights)))
}
# RMT calculation
{
library(ipw)
PropModel12 = ipwpoint(exposure = Adherence12yn,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ 1,
                       denominator = ~ AgeVisit1+Gender+GPNr,
                       trunc = 0.0125, data = as.data.frame(predivdat))
Weights12 = PropModel12$ipw.weights

ipwplot(weights = PropModel12$ipw.weights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(PropModel12$weights.trunc)))

PropModel23 = ipwpoint(exposure = Adherence23yn,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ Adherence12yn + DT12,
                       denominator = ~ Adherence12yn + AgeVisit1 + Gender + GPNr + DT12,
                       trunc = 0.0125, data = as.data.frame(predivdat))
Weights23 = PropModel23$ipw.weights

ipwplot(weights = PropModel23$ipw.weights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(PropModel12$weights.trunc)))
RMTweights= PropModel12$ipw.weights*PropModel23$weights.trunc
ipwplot(weights = RMTweights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(RMTweights)))

predivdat["Binarytime"] = predivdat$DT12<mean(predivdat$DT12)
MeasTimeModel = ipwpoint(exposure = Binarytime,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ Adherence12yn,
                       denominator = ~ Adherence12yn + AgeVisit1 + Gender + GPNr,
                       trunc = 0.0125, data = as.data.frame(predivdat))
MeasTimeWeights = MeasTimeModel$ipw.weights
ipwplot(weights = MeasTimeWeights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(MeasTimeWeights)))

}
# IPCW calculation
{
#Probability of dropout
Censdat = data[,c("StudieNr","Adherence12yn","AgeVisit1","Gender","DT12","Adherence23yn","DT23","WHOchange","GPNr","Randomisation")][data$Randomisation==1,]

Censdat["CensFU1"] = is.na(Censdat$Adherence12yn)|is.na(Censdat$DT12)
Censdat["CensFU2"] = is.na(Censdat$Adherence12yn)|is.na(Censdat$DT12)|is.na(Censdat$Adherence23yn)|is.na(Censdat$DT23)|is.na(Censdat$WHOchange)

CensModel12 = ipwpoint(exposure = CensFU1,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ 1,
                       denominator = ~ AgeVisit1+Gender+GPNr,
                       trunc = 0.0125, data = as.data.frame(Censdat))
Censweights12 = CensModel12$ipw.weights

{
CensModel12 = ipwpoint(exposure = CensFU1,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ 1,
                       denominator = ~ AgeVisit1+Gender+GPNr,
                       trunc = 0.0125, data = as.data.frame(Censdat))
Censweights12 = CensModel12$ipw.weights
ipwplot(weights = Weights12, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(Weights12)))

CensModel23 = ipwpoint(exposure = CensFU2,
                       family = 'binomial',
                       link = "logit",
                       numerator = ~ Adherence12yn,
                       denominator = ~ Adherence12yn + AgeVisit1 + Gender + GPNr + DT12,
                       trunc = 0.0125, data = as.data.frame(Censdat[Censdat$CensFU1 == 0,]))
Censweights23 = CensModel23$ipw.weights
ipwplot(weights = Weights12, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(Weights12)))

Censdat["Censweight12"] = Censweights12

Censdat23 = Censdat[Censdat$CensFU1 == 0,]
Censdat23["CensWeight23"] = Censweights23

#Match the censoring weights to the corresponding people
predivdat["Censweight12"] = numeric(nrow(predivdat))
for(i in 1:nrow(predivdat)){
  predivdat[i,"Censweight12"] = Censdat$Censweight12[which(Censdat$StudieNr == c(predivdat[i,"StudieNr"]))]
}

predivdat["Censweight23"] = numeric(nrow(predivdat))
for(i in 1:nrow(predivdat)){
  predivdat[i,"Censweight23"] = Censdat23$CensWeight23[which(Censdat23$StudieNr == c(predivdat[i,"StudieNr"]))]
}

}
}
# Final estimation of Causal effect
{
predivdat["IPTWweights"] = IPTWweights*predivdat$Censweight12*predivdat$Censweight23
predivdat["TACweights"] = TACweights*predivdat$Censweight12*predivdat$Censweight23
predivdat["RMTweights"] = RMTweights*MeasTimeWeights*predivdat$Censweight12*predivdat$Censweight23

ipwplot(weights = predivdat$RMTweights, logscale = FALSE,
        main = "Stabilized weights", xlim = c(0, max(predivdat$RMTweights)))

predivdat["AlwaysNever"] = NA

predivdat["AlwaysNever"][predivdat$Adherence12yn==1 & predivdat$Adherence23yn == 1,"AlwaysNever"] = 1
predivdat["AlwaysNever"][predivdat$Adherence12yn==0 & predivdat$Adherence23yn == 0,"AlwaysNever"] = 0

ANdat = predivdat[!is.na(predivdat$AlwaysNever),]

msmUnw = svyglm(WHOchange ~ AlwaysNever,design = svydesign(~1,data = ANdat))
summary(msmUnw)
msmIPTW = svyglm(WHOchange ~ AlwaysNever,design = svydesign(~1,weights = ~ANdat$IPTWweights,data = ANdat))
summary(msmIPTW)
msmTAC = svyglm(WHOchange ~ AlwaysNever,design = svydesign(~1,weights = ~ANdat$TACweights,data = ANdat))
summary(msmTAC)
msmRMT = svyglm(WHOchange ~ AlwaysNever,design = svydesign(~1,weights = ~ANdat$RMTweights,data = ANdat))
summary(msmRMT)
}
