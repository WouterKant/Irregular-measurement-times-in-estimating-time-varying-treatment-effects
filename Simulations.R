### Function definitions and packages

library(ipw)
library(ggplot2)
library(lme4)
library(gtools)
library(dplyr)
library(ggpubr)
library(jtools)
library(writexl)

cfYgenerator = function(n, cf, XO,XL1,XY,bO,bOA0,bOU,bOL0,bL1,bL1A0,bL1L0,bL1O,bL1U,bY,bYA0,bYA1,bYL1,
                        bYL0,bYU,RY11,RY00,U,L0){
  tA0 = rep(cf,n)
  tborderO = inv.logit(bO+bOA0*tA0+bOU*U+bOL0*L0)
  tO = XO < tborderO
  tborderL1 = inv.logit(bL1+bL1A0*tA0+ bL1L0*L0+bL1O*tO+bL1U*U)
  tL1 = XL1 < tborderL1
  tA1 = rep(cf,n)
  cf0Y = bY+bYA0*tA0+bYA1*tA1+bYL1*tL1+bYL0*L0+bYU*U+bYO*tO + ((tA0+tA1) == 2)*RY11 + ((tA0+tA1) == 0)*RY00
  return(cf0Y)
}
Experiment = function(B,n){
  ngrid = c(n)
  true = numeric(B)
  truerec = numeric(length(ngrid))
  
  results = matrix(data = NA, nrow = B, ncol = 6)
  cresults = matrix(data = NA, nrow = B, ncol = 6)
  
  resultsContrast = matrix(data = NA, nrow = B, ncol = 6)
  
  multNresults1 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNresults2 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNresults3 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNresults4 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNresults5 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNresults6 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  
  multNcresults1 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNcresults2 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNcresults3 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNcresults4 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNcresults5 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  multNcresults6 = matrix(data = NA, nrow = length(ngrid),ncol = 2)
  
  for(R in 1:length(ngrid)){
    n = ngrid[R]
    for(b in 1:B){
      ### Data Generating Process
      XU = runif(n)
      XL0 = runif(n)
      XA0 = runif(n)
      XO = runif(n)
      XL1 = runif(n)
      XA1 = runif(n)
      RY11 = rnorm(n,0,0.5)
      RY00 = rnorm(n,0,0.5)
      RY10 = rnorm(n,0,0.5)
      RY01 = rnorm(n,0,0.5)
      
      ## Generation of data
      borderU = inv.logit(bU)
      U = XU<borderU
      borderL0 = inv.logit(bL0+bL0U*U)
      L0 = XL0<borderL0
      borderA0 = inv.logit(bA0+bA0L0*L0)
      A0 = XA0 < borderA0
      borderO = inv.logit(bO+bOA0*A0+bOU*U+bOL0*L0)
      O = XO < borderO
      borderL1 = inv.logit(bL1+bL1A0*A0+bL1L0*L0 + bL1O*O+bL1U*U)
      L1 = XL1 < borderL1
      borderA1 = inv.logit(bA1+bA1A0*A0+bA1L1*L1+bA1L0*L0+bA1O*O+bA1U*U)
      A1 = XA1 < borderA1
      Y = bY+bYA0*A0+bYA1*A1+bYL1*L1+bYL0*L0+bYU*U+bYO*O + ((A0+A1) == 2)*RY11 + ((A0+A1) == 0)*RY00
      
      dat = as.data.frame(cbind(U,L0,A0,O,L1,A1,Y))
      
      
      
      # Counterfactual if A0, A1 = 0, no dropout
      cf0Y = cfYgenerator(n,0,XO,XL1,XY,bO,bOA0,bOU,bOL0,bL1,bL1A0,bL1L0,bL1O,bL1U,bY,bYA0,bYA1,bYL1,
                          bYL0,bYU,RY11,RY00,U,L0)
      cf1Y = cfYgenerator(n,1,XO,XL1,XY,bO,bOA0,bOU,bOL0,bL1,bL1A0,bL1L0,bL1O,bL1U,bY,bYA0,bYA1,bYL1,
                          bYL0,bYU,RY11,RY00,U,L0)
      true[b] = mean(cf1Y-cf0Y)
      
      ### Modelling
      #Naive
      naiveres = mean(dat$Y[dat$A0 == 1 & dat$A1 == 1]) - mean(dat$Y[dat$A0 == 0 & dat$A1 == 0])
      
      #ALAx # Normal IPTW
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit", numerator = ~A0, denominator = ~ A0+L1+L0, data = dat) 
      W = t0$ipw.weights * t1$ipw.weights
      ALAx = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #AOLAOx # Normal IPTW with stabilized O
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0+ O, denominator = ~ A0+L1+L0 + O, data = dat) 
      W = t0$ipw.weights * t1$ipw.weights
      AOLAOx = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #ALAOx # Time-as-confounder
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0, denominator = ~ A0+L1+L0 + O, data = dat) 
      W = t0$ipw.weights * t1$ipw.weights
      ALAOx = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #AOLAOv # Reweighting by N
      tn1 = ipwpoint(O,"binomial","logit",numerator = ~A0, denominator = ~ L0 + A0, data = dat)
      t0 = ipwpoint(A0,"binomial","logit", denominator = ~ L0, data = dat)
      t1 = ipwpoint(A1,"binomial","logit",numerator = ~A0 + O, denominator = ~ A0 + L1 + L0 + O, data = dat)
      W = tn1$ipw.weights * t0$ipw.weights * t1$ipw.weights
      AOLAOv = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #TACROT
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0, denominator = ~ A0+L1+L0 + O, data = dat)
      tn1 = ipwpoint(O,"binomial","logit",numerator = ~A0, denominator = ~ L0 + A0, data = dat)
      
      WTACROT = t0$ipw.weights * t1$ipw.weights * tn1$ipw.weights
      TACROT = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],WTACROT[dat$A0 == 1 & dat$A1 == 1]) -
        weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],WTACROT[dat$A0 == 0 & dat$A1 == 0])
      
      results[b,] = c(naiveres, ALAx,AOLAOx,ALAOx,AOLAOv,TACROT)
      cresults[b,] = results[b,]-true[b]
      print(c(R,b))
    }
    multNresults1[R,] = c(mean(results[,1]), sd(results[,1]))
    multNresults2[R,] = c(mean(results[,2]), sd(results[,2]))
    multNresults3[R,] = c(mean(results[,3]), sd(results[,3]))
    multNresults4[R,] = c(mean(results[,4]), sd(results[,4]))
    multNresults5[R,] = c(mean(results[,5]), sd(results[,5]))
    multNresults6[R,] = c(mean(results[,6]), sd(results[,6]))
#    multNresults6[R,] = c(mean(results[,6],na.rm = T), sd(results[,6],na.rm =T))
    truerec[R] = mean(true)
  }
  #fuldat = cbind(as.data.frame(rbind(multNresults1,multNresults2,multNresults3,multNresults4,multNresults5)), rep(c("Naive", "IPTW","IPTW-O","TAC","ROT"), each = length(ngrid)))
  fuldat = cbind(as.data.frame(rbind(multNresults2,multNresults3,multNresults4,multNresults5,multNresults6)), rep(truerec,each=5), rep(c("IPTW","IPTW-T","TAC","RMT","TAC+RMT"), each = length(ngrid)))
  fuldatCI = cbind(fuldat, fuldat[,1] - 1.96*fuldat[,2]/sqrt(B),fuldat[,1] + 1.96*fuldat[,2]/sqrt(B), ngrid)
  colnames(fuldatCI) = c('Mean', 'Standdev', 'True','Name', 'LowCI', 'HighCI','n')
  fuldatCI['Mean bias'] = fuldatCI$Mean-fuldatCI$True
  fuldatCI['MSE'] = fuldatCI$`Mean bias`^2+fuldatCI$Standdev^2
  fuldatCI['% bias'] = fuldatCI$`Mean bias`/fuldatCI$True*100
  return(fuldatCI)
}
PlotPerStrengthTimeBiasExperiment = function(B,n){
  Ogrid = c(0,0.5,1,1.5,2)
  true = numeric(B)
  truerec = numeric(length(Ogrid))
  results = matrix(data = NA, nrow = B, ncol = 5)
  resultsContrast = matrix(data = NA, nrow = B, ncol = 5)
  
  multNresults1 = matrix(data = NA, nrow = length(Ogrid),ncol = 2)
  multNresults2 = matrix(data = NA, nrow = length(Ogrid),ncol = 2)
  multNresults3 = matrix(data = NA, nrow = length(Ogrid),ncol = 2)
  multNresults4 = matrix(data = NA, nrow = length(Ogrid),ncol = 2)
  multNresults5 = matrix(data = NA, nrow = length(Ogrid),ncol = 2)
  
  for(R in 1:length(Ogrid)){
    bA1O = Ogrid[R]
    for(b in 1:B){
      ### Data Generating Process
      XU = runif(n)
      XL0 = runif(n)
      XA0 = runif(n)
      XO = runif(n)
      XL1 = runif(n)
      XA1 = runif(n)
      RY11 = rnorm(n,0,1)
      RY00 = rnorm(n,0,1)
      RY10 = rnorm(n,0,1)
      RY01 = rnorm(n,0,1)
      
      ## Generation of data
      borderU = inv.logit(bU)
      U = XU<borderU
      borderL0 = inv.logit(bL0+bL0U*U)
      L0 = XL0<borderL0
      borderA0 = inv.logit(bA0+bA0L0*L0)
      A0 = XA0 < borderA0
      borderO = inv.logit(bO+bOA0*A0+bOU*U+bOL0*L0)
      O = XO < borderO
      borderL1 = inv.logit(bL1+bL1A0*A0+bL1L0*L0 + bL1O*O+bL1U*U)
      L1 = XL1 < borderL1
      borderA1 = inv.logit(bA1+bA1A0*A0+bA1L1*L1+bA1L0*L0+bA1O*O+bA1U*U)
      A1 = XA1 < borderA1
      Y = bY+bYA0*A0+bYA1*A1+bYL1*L1+bYL0*L0+bYU*U+bYO*O + ((A0+A1) == 2)*RY11 + ((A0+A1) == 0)*RY00
      
      dat = as.data.frame(cbind(U,L0,A0,O,L1,A1,Y))
      
      # Counterfactual if A0, A1 = 0, no dropout
      cf0Y = cfYgenerator(n,0,XO,XL1,XY,bO,bOA0,bOU,bOL0,bL1,bL1A0,bL1L0,bL1O,bL1U,bY,bYA0,bYA1,bYL1,
                          bYL0,bYU,RY11,RY00,U,L0)
      cf1Y = cfYgenerator(n,1,XO,XL1,XY,bO,bOA0,bOU,bOL0,bL1,bL1A0,bL1L0,bL1O,bL1U,bY,bYA0,bYA1,bYL1,
                          bYL0,bYU,RY11,RY00,U,L0)
      true[b] = mean(cf1Y-cf0Y)
      
      ### Modelling
      #Naive
      naiveres = mean(dat$Y[dat$A0 == 1 & dat$A1 == 1]) - mean(dat$Y[dat$A0 == 0 & dat$A1 == 0])
      
      #ALAx # Normal IPTW
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit", numerator = ~A0, denominator = ~ A0+L1+L0, data = dat) 
      W = t0$ipw.weights * t1$ipw.weights
      ALAx = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #AOLAOx # Normal IPTW with stabilized O
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0+ O, denominator = ~ A0+L1+L0 + O, data = dat) 
      W = t0$ipw.weights * t1$ipw.weights
      AOLAOx = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #ALAOx # Time-as-confounder
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0, denominator = ~ A0+L1+L0 + O, data = dat) 
      W = t0$ipw.weights * t1$ipw.weights
      ALAOx = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #AOLAOv # Reweighting by N
      tn1 = ipwpoint(O,"binomial","logit",numerator = ~A0, denominator = ~ L0 + A0, data = dat)
      t0 = ipwpoint(A0,"binomial","logit", denominator = ~ L0, data = dat)
      t1 = ipwpoint(A1,"binomial","logit",numerator = ~A0 + O, denominator = ~ A0 + L1 + L0 + O, data = dat)
      W = tn1$ipw.weights * t0$ipw.weights * t1$ipw.weights
      AOLAOv = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      #TACROT
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0, denominator = ~ A0+L1+L0 + O, data = dat)
      tn1 = ipwpoint(O,"binomial","logit",numerator = ~A0, denominator = ~ L0 + A0, data = dat)
      
      WTACROT = t0$ipw.weights * t1$ipw.weights * tn1$ipw.weights
      TACROT = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],WTACROT[dat$A0 == 1 & dat$A1 == 1]) -
        weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],WTACROT[dat$A0 == 0 & dat$A1 == 0])
      
      
      results[b,] = c(naiveres, ALAx,ALAOx,AOLAOv,TACROT)
      print(c(R,b))
    }
    #contrastresults = results-mean(true)
    #multNresults1[R,] = c(mean(contrastresults[,1]), sd(contrastresults[,1]))
    #multNresults2[R,] = c(mean(contrastresults[,2]), sd(contrastresults[,2]))
    #multNresults3[R,] = c(mean(contrastresults[,3]), sd(contrastresults[,3]))
    #multNresults4[R,] = c(mean(contrastresults[,4]), sd(contrastresults[,4]))
    #multNresults5[R,] = c(mean(contrastresults[,5]), sd(contrastresults[,5]))
    multNresults1[R,] = c(mean(results[,1]), sd(results[,1]))
    multNresults2[R,] = c(mean(results[,2]), sd(results[,2]))
    multNresults3[R,] = c(mean(results[,3]), sd(results[,3]))
    multNresults4[R,] = c(mean(results[,4]), sd(results[,4]))
    multNresults5[R,] = c(mean(results[,5]), sd(results[,5]))
    truerec[R] = mean(true)
  }
  
  #fuldat = cbind(as.data.frame(rbind(multNresults1,multNresults2,multNresults3,multNresults4,multNresults5)), rep(c("Naive", "IPTW","IPTW-O","TAC","ROT"), each = length(ngrid)))
  fuldat = cbind(as.data.frame(rbind(multNresults2,multNresults3,multNresults4,multNresults5)),rep(truerec,each=4), rep(c("IPTW","TAC","RMT","TAC+RMT"), each = length(Ogrid)))
  fuldatCI = cbind(fuldat, fuldat[,1] - 1.96*fuldat[,2]/sqrt(B),fuldat[,1] + 1.96*fuldat[,2]/sqrt(B), Ogrid)
  colnames(fuldatCI) = c('Mean', 'Standdev', 'True','Name', 'LowCI', 'HighCI','O')
  fuldatCI['Mean bias'] = fuldatCI$Mean-fuldatCI$True
  fuldatCI['MSE'] = fuldatCI$`Mean bias`^2+fuldatCI$Standdev^2
  fuldatCI['% bias'] = fuldatCI$`Mean bias`/fuldatCI$True*100
  
#  fuldatCI = cbind(fuldat, fuldat[,1] - 1.96*fuldat[,2]/sqrt(B),fuldat[,1] + 1.96*fuldat[,2]/sqrt(B), Ogrid)
#  colnames(fuldatCI) = c('Mean', 'Standdev', 'Name', 'LowCI', 'HighCI','bAO')
#  pd <- position_dodge(0.15) # move the dots to the left and right
#  fuldatCI['MSE'] = fuldatCI$Mean^2+fuldatCI$Standdev^2
  return(fuldatCI)
}
SelectionBiasExperiment = function(B,n){
  Selgrid = c(0,1,2,3,4)
  true = numeric(B)
  truerec = numeric(length(Selgrid))
  results = matrix(data = NA, nrow = B, ncol = 5)
  resultsContrast = matrix(data = NA, nrow = B, ncol = 5)
  
  multNresults1 = matrix(data = NA, nrow = length(Selgrid),ncol = 2)
  multNresults2 = matrix(data = NA, nrow = length(Selgrid),ncol = 2)
  multNresults3 = matrix(data = NA, nrow = length(Selgrid),ncol = 2)
  multNresults4 = matrix(data = NA, nrow = length(Selgrid),ncol = 2)
  multNresults5 = matrix(data = NA, nrow = length(Selgrid),ncol = 2)
  
  for(R in 1:length(Selgrid)){
    for(b in 1:B){
      ### Data Generating Process
      XU = runif(n)
      XL0 = runif(n)
      XA0 = runif(n)
      XO = runif(n)
      XL1 = runif(n)
      XA1 = runif(n)
      RY11 = rnorm(n,0,0.1)
      RY00 = rnorm(n,0,0.1)
      RY10 = rnorm(n,0,0.1)
      RY01 = rnorm(n,0,0.1)
      
      ## Generation of data
      borderU = inv.logit(bU)
      U = XU<borderU
      borderL0 = inv.logit(bL0+bL0U*U)
      L0 = XL0<borderL0
      borderA0 = inv.logit(bA0+bA0L0*L0)
      A0 = XA0 < borderA0
      borderO = inv.logit(bO+bOA0*A0+bOU*U+bOL0*L0)
      O = XO < borderO
      borderL1 = inv.logit(bL1+bL1A0*A0+bL1L0*L0 + bL1O*O+bL1U*U)
      L1 = XL1 < borderL1
      borderA1 = inv.logit(bA1+bA1A0*A0+bA1L1*L1+bA1L0*L0+bA1O*O+bA1U*U)
      A1 = XA1 < borderA1
      Y = bY+bYA0*A0+bYA1*A1+bYL1*L1+bYL0*L0+bYU*U+bYO*O + ((A0+A1) == 2)*RY11 + ((A0+A1) == 0)*RY00
      
      dropout = runif(n)>=inv.logit(2-Selgrid[R]*A1 + Selgrid[R]*L1)
      datincldropout = as.data.frame(cbind(U,L0,A0,O,L1,A1,Y,dropout))
      dat = datincldropout[!datincldropout$dropout,]
      
      #dat = as.data.frame(cbind(U,L0,A0,O,L1,A1,Y,dropout))
      # Counterfactual if A0, A1 = 0, no dropout
      cf0Y = cfYgenerator(n,0,XO,XL1,XY,bO,bOA0,bOU,bOL0,bL1,bL1A0,bL1L0,bL1O,bL1U,bY,bYA0,bYA1,bYL1,
                          bYL0,bYU,RY11,RY00,U,L0)
      cf1Y = cfYgenerator(n,1,XO,XL1,XY,bO,bOA0,bOU,bOL0,bL1,bL1A0,bL1L0,bL1O,bL1U,bY,bYA0,bYA1,bYL1,
                          bYL0,bYU,RY11,RY00,U,L0)
      true[b] = mean(cf1Y-cf0Y)
      
      ### Modelling
      #Naive
      naiveres = mean(dat$Y[dat$A0 == 1 & dat$A1 == 1]) - mean(dat$Y[dat$A0 == 0 & dat$A1 == 0])
      
      
      # ROT
      tn1 = ipwpoint(O,"binomial","logit",numerator = ~A0, denominator = ~ L0 + A0, data = dat)
      t0 = ipwpoint(A0,"binomial","logit", denominator = ~ L0, data = dat)
      t1 = ipwpoint(A1,"binomial","logit",numerator = ~A0 + O, denominator = ~ A0 + L1 + L0 + O, data = dat)
      W = tn1$ipw.weights * t0$ipw.weights * t1$ipw.weights
      ROT = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      # TAC 
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = dat) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0, denominator = ~ A0+L1+L0 + O, data = dat) 
      W = t0$ipw.weights * t1$ipw.weights
      TAC = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      # ROT + ROT
      tn1 = ipwpoint(O,"binomial","logit",numerator = ~A0, denominator = ~ L0 + A0, data = datincldropout)
      tn2 = ipwpoint(dropout,"binomial","logit", denominator = ~ L1 + A1, data = datincldropout)
      t0 = ipwpoint(A0,"binomial","logit", denominator = ~ L0, data = datincldropout)
      t1 = ipwpoint(A1,"binomial","logit",numerator = ~A0 + O, denominator = ~ A0 + L1 + L0 + O, data = datincldropout)
      W = tn1$ipw.weights[!datincldropout$dropout] * tn2$ipw.weights[!datincldropout$dropout] * t0$ipw.weights[!datincldropout$dropout] * t1$ipw.weights[!datincldropout$dropout]
      ROTROT = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      # TAC + ROT
      tn2 = ipwpoint(dropout,"binomial","logit", denominator = ~ L1 + A1, data = datincldropout)
      t0 = ipwpoint(A0, "binomial", "logit", denominator = ~ L0, data = datincldropout) 
      t1 = ipwpoint(A1, "binomial", "logit",numerator = ~A0, denominator = ~ A0+L1+L0 + O, data = datincldropout) 
      W = tn2$ipw.weights[!datincldropout$dropout] * t0$ipw.weights[!datincldropout$dropout] * t1$ipw.weights[!datincldropout$dropout]
      TACROT = weighted.mean(dat$Y[dat$A0 == 1 & dat$A1 == 1],W[dat$A0 == 1 & dat$A1 == 1]) - weighted.mean(dat$Y[dat$A0 == 0 & dat$A1 == 0],W[dat$A0 == 0 & dat$A1 == 0])
      
      
      
      results[b,] = c(naiveres, ROT, TAC, ROTROT, TACROT)
      print(c(R,b))
    }
    contrastresults = results-mean(true)
    multNresults1[R,] = c(mean(results[,1]), sd(results[,1]))
    multNresults2[R,] = c(mean(results[,2]), sd(results[,2]))
    multNresults3[R,] = c(mean(results[,3]), sd(results[,3]))
    multNresults4[R,] = c(mean(results[,4]), sd(results[,4]))
    multNresults5[R,] = c(mean(results[,5]), sd(results[,5]))
    truerec[R] = mean(true)
    
  }
  
  
  #fuldat = cbind(as.data.frame(rbind(multNresults1,multNresults2,multNresults3,multNresults4,multNresults5)), rep(c("Naive", "IPTW","IPTW-O","TAC","ROT"), each = length(ngrid)))
  fuldat = cbind(as.data.frame(rbind(multNresults3,multNresults4,multNresults5)),rep(truerec,each=3), rep(c("TAC","RMT","TAC+RMT"), each = length(Selgrid)))
  
  fuldatCI = cbind(fuldat, fuldat[,1] - 1.96*fuldat[,2]/sqrt(B),fuldat[,1] + 1.96*fuldat[,2]/sqrt(B), Selgrid)
  colnames(fuldatCI) = c('Mean', 'Standdev', 'True','Name', 'LowCI', 'HighCI','O')
  fuldatCI['Mean bias'] = fuldatCI$Mean-fuldatCI$True
  fuldatCI['MSE'] = fuldatCI$`Mean bias`^2+fuldatCI$Standdev^2
  fuldatCI['% bias'] = fuldatCI$`Mean bias`/fuldatCI$True*100
  
  return(fuldatCI)
}

# SCENARIO 1: Direct confounding (DC)
{
  set.seed(123)
  bU = 0
  bL0 = 0
  bL0U = 0
  bA0 = 0
  bA0L0 = -0.5
  bA0U = 0
  bO = 0
  bOA0 = 2
  bOL0 = 0
  bOU = 0
  bL1 = 0
  bL1A0 = 1
  bL1L0 = 1
  bL1O = -2
  bL1U = 0
  bA1 = 0
  bA1A0 = 0
  bA1L1 = -1
  bA1L0 = -1
  bA1O = 1
  bA1U = 0
  bY = 0
  bYA0 = 1
  bYA1 = 1
  bYL0 = 1
  bYL1 = 2
  bYU = 0
  bYO = 2
  
  DCresults = Experiment(10000,2000)
  ggplot(DCresults, aes(x=Name, y=`% bias`)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = 'red') +
    geom_errorbar(aes(ymin=((LowCI-True)/True*100), ymax=((HighCI-True)/True*100)), width=0.6, position = pd) +
    theme_minimal()+
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 0, l = 0))) +
    ylab("Mean percentage bias") + xlab("Reweighting method") #+ ggtitle("Mean bias under DC")
}
write_xlsx(DCresults,"DCresults.xlsx")
DCresults %>% 
  mutate_if(is.numeric, round,digits=4)

# Scenario 2: Confounding through Measured Variables (CMV)
{
  set.seed(123)
  bU = 0
  bL0 = 0
  bL0U = 2
  bA0 = 0
  bA0L0 = -1
  bA0U = 0
  bO = 0
  bOA0 = 2
  bOL0 = -2
  bOU = 0
  bL1 = 0
  bL1A0 = 1
  bL1L0 = 1
  bL1O = 0
  bL1U = 0
  bA1 = 0
  bA1A0 = 0
  bA1L1 = -1
  bA1L0 = -1
  bA1O = 2
  bA1U = 0
  bY = 0
  bYA0 = 1
  bYA1 = 1
  bYL0 = 2
  bYL1 = 1
  bYU = 0
  bYO = 0
  
  CMVresults = Experiment(10000,2000)
  
  ggplot(CMVresults, aes(x=Name, y=`% bias`)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = 'red') +
    geom_errorbar(aes(ymin=((LowCI-True)/True*100), ymax=((HighCI-True)/True*100)), width=0.6, position = pd) +
    theme_minimal()+
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 0, l = 0))) +
    ylab("Mean percentage bias") + xlab("Reweighting method") #+ ggtitle("Mean bias under DC") 
  }
write_xlsx(CMVresults,"CMVresults.xlsx")
CMVresults %>% 
  mutate_if(is.numeric, round,digits=4)

# SCENARIO 3: Confounding through Unmeasured Variables (CUV)
{
  set.seed(123)
  bU = 0
  bL0 = 0
  bL0U = 2
  bA0 = 0
  bA0L0 = -1
  bA0U = 0
  bO = 0
  bOA0 = 2
  bOL0 = 0
  bOU = 2
  bL1 = 0
  bL1A0 = 1
  bL1L0 = 1
  bL1O = -2
  bL1U = -2
  bA1 = 0
  bA1A0 = 0
  bA1L1 = -1
  bA1L0 = -1
  bA1O = 2
  bA1U = 0
  bY = 0
  bYA0 = 1
  bYA1 = 1
  bYL0 = 2
  bYL1 = 2
  bYU = 2
  bYO = 0
  CUVresults = Experiment(10000,2000)
  ggplot(CUVresults, aes(x=Name, y=`% bias`)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = 'red') +
    geom_errorbar(aes(ymin=((LowCI-True)/True*100), ymax=((HighCI-True)/True*100)), width=0.6, position = pd) +
    theme_minimal()+
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 0, l = 0))) +
    ylab("Mean percentage bias") + xlab("Reweighting method") #+ ggtitle("Mean bias under DC") 
}
write_xlsx(CUVresults,"CUVresults.xlsx")
CUVresults %>% 
  mutate_if(is.numeric, round,digits=4)

# SENARIO 4: Selection bias and CUV
{
  set.seed(123)
  bU = 0
  bL0 = 0
  bL0U = 2
  bA0 = 0
  bA0L0 = -1
  bA0U = 0
  bO = 0
  bOA0 = 2
  bOL0 = 0
  bOU = 2
  bL1 = 0
  bL1A0 = 1
  bL1L0 = 1
  bL1O = -2
  bL1U = -2
  bA1 = 0
  bA1A0 = 0
  bA1L1 = -1
  bA1L0 = -1
  bA1O = 3
  bA1U = 0
  bY = 0
  bYA0 = 1
  bYA1 = 1
  bYL0 = 2
  bYL1 = 2
  bYU = 2
  bYO = 0
  SelectBiasResults = SelectionBiasExperiment(10000,2000)
  ggplot(SelectBiasResults, aes(x=O, y=`Mean bias`, colour=Name, shape = Name)) + 
    # geom_errorbar(aes(ymin=LowCI, ymax=HighCI), width=1, position = pd) +
    geom_line() +
    geom_hline(yintercept = 0, color = 'black') +
    geom_point(size = 2.5) +
    theme_minimal() +
    xlab("Strength of selection bias") + ylab("Mean percentage bias") #+ ggtitle("All confounding biases \n by measurement time")
  }
write_xlsx(SelectBiasResults,"SelectBiasResults.xlsx")

# EXTRA: Bias per strength of confounding bias for plot
{
  set.seed(123)
  bU = 0
  bL0 = 0
  bL0U = 2
  bA0 = 0
  bA0L0 = -1
  bA0U = 0
  bO = 0
  bOA0 = 2
  bOL0 = 0
  bOU = 2
  bL1 = 0
  bL1A0 = 1
  bL1L0 = 1
  bL1O = -2
  bL1U = -2
  bA1 = 0
  bA1A0 = 0
  bA1L1 = -1
  bA1L0 = -1
  bA1O = 2
  bA1U = 0
  bY = 0
  bYA0 = 1
  bYA1 = 1
  bYL0 = 2
  bYL1 = 2
  bYU = 2
  bYO = 0
  BiasPerO = PlotPerStrengthTimeBiasExperiment(10000,2000)
  pd=position_dodge(0)
  ggplot(BiasPerO, aes(x=O, y=`Mean bias`, colour=Name)) + 
    # geom_errorbar(aes(ymin=LowCI, ymax=HighCI), width=1, position = pd) +
    geom_line() +
    geom_hline(yintercept = 0, color = 'black') +
    geom_point() +
    theme_minimal() +
    xlab("Causal effect of measurement time \n on treatment assignment") + ylab("Mean percentage bias") #+ ggtitle("All confounding biases \n by measurement time")
  
}
write_xlsx(BiasPerO,"BiasPerO.xlsx")


