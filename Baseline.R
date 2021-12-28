#source: https://www.leonardoegidi.com/covid-19

############
# Covid 19 #
############

# load packages
library(rstanarm)
library(ggplot2)
library(bayesplot)
library(tidyverse)
library(lme4)
library(arm)
library(loo) 
library(kableExtra)
#_______________________________________________________________________________
#setwd("C:/Users/m/Desktop/Project/tables")
# import daily data from github
regioni <- read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv")
head(regioni)
#write.csv(regioni,"C:/Users/m/Desktop/Project/tables//Dataset.csv", row.names=FALSE)
regione <- filter(regioni, codice_regione==6) #fvg
write.csv(regione,"C://Users//m//Desktop//Project//Dataset.csv", row.names=FALSE)
#_______________________________________________________________________________
# store variables
n_osp_giornalieri <- length(regioni$totale_ospedalizzati) #check having 21 rows (regions)

l_times <- length(unique(regioni$data))#see n° of distinct times (1 means ok, uniform times)
times <- 1:l_times

n_reg <- length(unique(regioni$denominazione_regione)) #number of regions (check==21)
n <- dim(regioni)[1] #(check==21)
#_______________________________________________________________________________
# lag source data (here not lagged indeed)
osp_lag <- c() # lagged intensive care units
initialize <- rep(0, n_reg) #init 1x21 empty vector
osp_lag <- initialize
for (j in 0:n){
  osp_lag[j] <- regioni$totale_ospedalizzati[j] #right counter?
}

icu_lag <- c() # lagged intensive care units
initialize <- rep(0, n_reg)
icu_lag <- initialize
for (j in 0:n){
  icu_lag[j] <- regioni$terapia_intensiva[j] #right counter?
}

tamp_lag <- c() # lagged swabs
initialize <- rep(0, n_reg)
tamp_lag <- initialize
for (j in 0:n){
  tamp_lag[j] <- regioni$tamponi[j] #right counter?
}
#_______________________________________________________________________________
# re-codifying distinct times
times <- as.factor(regioni$data)
levels(times) <- c(1:l_times)
times <- as.numeric(times)
regioni$codice_regione <- as.factor(regioni$codice_regione) #levels 1

# if launched before 0.00 AM
# yesterday_day <- substr(format(Sys.Date(), format="%d%b-%Y"),1,5)
# exact_day <- substr(format(Sys.Date()+1, format="%d%b-%Y"),1,5)

yesterday_day_long <- Sys.Date()
exact_day_long <- Sys.Date()+1
#_______________________________________________________________________________
# plot observed total hospitalized & ICU collected until today in each region (response variables):
df <- data.frame(times = times,
                 reg = unlist(list(regioni$denominazione_regione)),
                 y = regioni$totale_ospedalizzati)
ggplot(df, aes(x=times, y =y))+geom_point(aes(x=times, y =y))+xlab("Days")+
  ylab("Total Hospitalized")+ scale_x_discrete(limit = c( 8, 15, 22, 29, 34 ), 
                                               labels = c( "2-3", "9-3", "16-3", 
                                               "23-3", "28-3") )+
  facet_wrap("reg", scales ="free")+ theme(strip.text.x = element_text(size = 12, 
                                    colour = "black"),axis.text.x = element_text(face="bold", 
                                                                                 color="#993333", 
                                                                                 angle=45, size =9),
                                    axis.title.x = element_text(size=22),
                                    axis.title.y = element_text(size = 22))
#just 1 data point!
#_______________________________________________________________________________
resp <- regioni$totale_ospedalizzati #response variable Y
n.iter <- 100 #default iterations number

# A) model 1: GLM Poisson canonical link 1 covariate (times)

mod.1 <- stan_glm(resp ~ times, family=poisson, data=regioni, iter=n.iter) # bayesian
mod.1.cl <- glm(resp ~ times, family=poisson, data=regioni) # classical

loo.glm.1 <- loo(mod.1)$estimates[3,1] # better bayesian models, lower LOOIC values
summary.glm.1 <- summary(mod.1.cl) 


# B) model 2: GLM, canonical link Poisson 2 covariates (times, times^2)
mod.2 <- stan_glm(resp ~ times + I(times^2), family=poisson, data=regioni, iter=n.iter)
mod.2.cl <- glm(resp ~ times + I(times^2), family=poisson, data=regioni)

loo.glm.2 <- loo(mod.2)$estimates[3,1]
summary.glm.2 <- summary(mod.2.cl)

# c(loo.glm.1, loo.glm.2) #[1] 926215.3 903766.7 (2nd model looks better)

# C) model 3: GLM, canonical link Poisson 3 covariates (times, times^2, swabs)
mod.3 <- stan_glm(resp ~  times + I (times^2) + tamponi, family = poisson, data = regioni, iter = n.iter)
mod.3.cl <- glm(resp ~ times + I (times^2) + tamponi, family = poisson, data = regioni)

loo.glm.3 <- loo(mod.3)$estimates[3,1]
summary.glm.3 <- summary(mod.3.cl)


# D) model 4: GLM, canonical link Poisson 4 covariates (times, times^2, swabs, lagged hospitalizations)
mod.4 <- stan_glm(resp ~  times + I(times^2) + tamponi + osp_lag, family=poisson, data=regioni, iter=n.iter)
mod.4.cl <- glm(resp ~  times + I (times^2) + tamponi + osp_lag, family=poisson, data=regioni)

loo.glm.4 <- loo(mod.4)$estimates[3,1]
summary.glm.4 <- summary(mod.4.cl)

# c(loo.glm.3, loo.glm.4) #[1] 384315.7 300293.7 (4th model looks better)
#_______________________________________________________________________________
#Hierarchical Poisson models

# E) model 5: GLMM, hierarchical Poisson 3 covariates (regions, times, times^2)

mod.hier.1 <- stan_glmer(resp ~ (1|codice_regione)+ times + I(times^2), family=poisson, data=regioni, iter=n.iter)
mod.hier.1.cl <- glmer(resp ~ (1|codice_regione)+ times + I(times^2),family=poisson, data=regioni)

loo.hier.1 <- loo(mod.hier.1)$estimates[3,1]
d1 <- display(mod.hier.1.cl)
dic.hier.1 <- as.double(d1$DIC)


# F) model 6: GLMM, hierarchical Poisson 4 covariates (regions, times, times^2, swabs)
mod.hier.2 <- stan_glmer(resp ~ (1|codice_regione)+ times + I(times^2) + tamponi, data=regioni, family=poisson, iter=n.iter)
mod.hier.2.cl <- glmer(resp ~ (1|codice_regione)+ times + I(times^2) + tamponi, data=regioni, family=poisson)

loo.hier.2 <- loo(mod.hier.2)$estimates[3,1]
d2 <- display(mod.hier.2.cl)
dic.hier.2 <- as.double(d2$DIC)

# c(loo.hier.1, loo.hier.2) # [1] 21594.55 12976.74 (6th model looks better)
#_______________________________________________________________________________
# redefine variables
regioni_ristr <- cbind(times, osp_lag, regioni)
regioni_tibble <- as_tibble(regioni_ristr)
regioni_ristr<- regioni_tibble %>% filter( 
          (times > 3 & (denominazione_regione == "Campania"|
                        denominazione_regione == "Puglia" |
                        denominazione_regione == "Basilicata"|
                        denominazione_regione == "Molise" |
                        denominazione_regione == "Calabria")) |
          (times > 1 & denominazione_regione == "Sicilia")|
          (times > 8 & denominazione_regione == "Sardegna")|
          (times > 0 & (denominazione_regione == "Piemonte"|
                       denominazione_regione == "Valle d'Aosta" |
                       denominazione_regione == "Abruzzo"|
                       denominazione_regione == "P.A. Bolzano" |
                       denominazione_regione == "Emilia Romagna"|
                       denominazione_regione == "Friuli Venezia Giulia"|
                       denominazione_regione == "Lazio" |
                       denominazione_regione == "Liguria"|
                       denominazione_regione == "Lombardia" |
                       denominazione_regione == "Marche"|
                       denominazione_regione == "P.A. Trento"|
                       denominazione_regione == "Toscana"|
                       denominazione_regione == "Umbria"|
                       denominazione_regione == "Veneto")))
# lockdown variable
lockdown <- c()
eff<-function(xx) { vec= rep(0, length(xx))
for (i in 1:length(xx)) {
  if (xx[i] <= 15)
  {vec[i]=xx[i]^2/225}
  else
  {vec[i]=1}}
vec
}

# defined to be zero before march 8 2020
lockdown <- c(rep(0, 14), eff(1:(l_times-14)))
resp <- regioni_ristr$totale_ospedalizzati
times <- regioni_ristr$times
lockdown_rep <- c()
for (j in 1: length(resp)){
  lockdown_rep[j] <- lockdown[times[j]]
}

n_reg <- length(unique(regioni_ristr$denominazione_regione))
times_scaled <- times - mean(times)

dummy_reg <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0)
dummy_reg_rep <- c()
veri_codici <- unique(regioni_ristr$codice_regione)

for (j in 1: dim(regioni_ristr)[1]){
  dummy_reg_rep[j] <- dummy_reg[
    (1:(n_reg))[regioni_ristr$codice_regione[j]] ]
}

# G) model 7: GLMM, hierarchical Poisson 4 covariates (regions, times, times^2, lockdown dummy)
# N.B. distinct starting points & times are scaled

mod.hier.3 <- stan_glmer(resp ~ (1|codice_regione)+ times_scaled + I(times_scaled^2) + lockdown_rep, 
                         data=regioni_ristr, family=poisson, iter=n.iter)
mod.hier.3.cl <- glmer(resp ~ (1|codice_regione)+ times_scaled + I(times_scaled^2) + lockdown_rep, 
                       data=regioni_ristr, family=poisson)

loo.hier.3 <- loo(mod.hier.3)$estimates[3,1]
d3 <- display(mod.hier.3.cl)
dic.hier.3 <- as.double(d3$DIC)

# H) model 8: GLMM, hierarchical Poisson 6 covariates (regions, times, times^2, lockdown dummy, swabs)
# N.B. distinct starting points & times are scaled

mod.hier.4 <- stan_glmer(resp ~ (1|codice_regione)+ times_scaled + I(times_scaled^2) + tamponi + lockdown_rep, 
                         data=regioni_ristr, family=poisson, iter=n.iter)
mod.hier.4.cl <- glmer(resp ~ (1|codice_regione)+ times_scaled + I(times_scaled^2) + tamponi + lockdown_rep, 
                       data=regioni_ristr, family=poisson)

loo.hier.4 <- loo(mod.hier.4)$estimates[3,1]
d4 <- display(mod.hier.4.cl)
dic.hier.4 <- as.double(d4$DIC)

# I) model 9: GLMM, hierarchical Poisson 4 covariates (regions, times^2, swabs, interaction lockdown-time)
# N.B. distinct starting points & times are scaled

mod.hier.4.int <- stan_glmer(resp ~ (1|codice_regione) + lockdown_rep*times_scaled + I(times_scaled^2) + tamponi,
                             data=regioni_ristr, family=poisson, iter=n.iter)
mod.hier.4.cl.int <- glmer(resp ~ (1|codice_regione)+ lockdown_rep*times_scaled + I(times_scaled^2) + tamponi,
                           data=regioni_ristr, family=poisson)

loo.hier.4.int <- loo(mod.hier.4.int)$estimates[3,1]
d4.int <- display(mod.hier.4.cl.int)
dic.hier.4.int <- as.double(d4.int$DIC)

# J) model 10: GLMM, hierarchical Poisson 4 covariates (regions, interaction lockdown-time, interaction lockdown-times^2, swabs)
# N.B. distinct starting points & times are scaled

mod.hier.4.int.2 <- stan_glmer(resp ~ (1|codice_regione) + lockdown_rep*times_scaled + lockdown_rep*I(times_scaled^2) + tamponi,
                               data=regioni_ristr, family=poisson, iter=n.iter)
mod.hier.4.cl.int.2 <- glmer(resp ~ (1|codice_regione) + lockdown_rep*times_scaled + lockdown_rep*I(times_scaled^2) + tamponi ,
                             data=regioni_ristr, family=poisson)

loo.hier.4.int.2 <- loo(mod.hier.4.int.2)$estimates[3,1]
d4.int.2 <- display(mod.hier.4.cl.int.2)
dic.hier.4.int.2 <- as.double(d4.int.2$DIC)

# K) model 11: GLMM, hierarchical Poisson 6 covariates (regions, times, times^2, swabs, lockdown, lagged hosp.)
# N.B. distinct starting points & times are scaled

mod.hier.5 <- stan_glmer(resp ~ (1|codice_regione) + times_scaled + I(times_scaled^2) + tamponi + lockdown_rep + osp_lag,
                         data=regioni_ristr, family=poisson, iter=n.iter)
mod.hier.5.cl <- glmer(resp ~ (1|codice_regione) + times_scaled + I(times_scaled^2) + tamponi + lockdown_rep + osp_lag,
                       data=regioni_ristr, family=poisson)

loo.hier.5 <- loo(mod.hier.5)$estimates[3,1]
d5 <- display(mod.hier.5.cl)
dic.hier.5 <- as.double(d5$DIC)

# L) model 12: GLMM, hierarchical Poisson 5 covariates (regions, interaction lockdown-time, times^2, swabs, lagged hosp.)
# N.B. distinct starting points & times are scaled

mod.hier.5.int <- stan_glmer(resp ~ (1|codice_regione) + lockdown_rep*times_scaled + I(times_scaled^2) + tamponi + osp_lag,
                             data=regioni_ristr, family=poisson, iter=n.iter)
mod.hier.5.cl.int <- glmer(resp ~ (1|codice_regione) + lockdown_rep*times_scaled + I(times_scaled^2) + tamponi + osp_lag,
                           data=regioni_ristr, family=poisson)

loo.hier.5.int <- loo(mod.hier.5.int)$estimates[3,1]
d5.int <- display(mod.hier.5.cl.int)
dic.hier.5.int <- as.double(d5.int$DIC)

# M) model 13: GLMM, hierarchical Poisson 5 covariates (regions, interaction lockdown-time, interaction lockdown-times^2, swabs, lagged hosp.)
# N.B. distinct starting points & times are scaled

mod.hier.5.int.2 <- stan_glmer(resp ~ (1|codice_regione) + lockdown_rep*times_scaled + lockdown_rep*I(times_scaled^2) + tamponi + osp_lag,
                               data=regioni_ristr, family=poisson, iter=n.iter)
mod.hier.5.cl.int.2 <- glmer(resp ~ (1|codice_regione) + lockdown_rep*times_scaled + lockdown_rep*I(times_scaled^2) + tamponi + osp_lag,
                             data=regioni_ristr, family=poisson)

loo.hier.5.int.2 <- loo(mod.hier.5.int.2)$estimates[3,1]
d5.int.2 <- display(mod.hier.5.cl.int.2)
dic.hier.5.int.2 <- as.double(d5.int.2$DIC)

# c(loo.hier.3, loo.hier.4, loo.hier.4.int, loo.hier.4.int.2, loo.hier.5, loo.hier.5.int, loo.hier.5.int.2)
#[1] 21368.48 12752.79 12635.67 12658.37 11267.50 11239.85 11244.17
# c(dic.hier.3, dic.hier.4, dic.hier.4.int, dic.hier.4.int.2, dic.hier.5, dic.hier.5.int, dic.hier.5.int.2)
#[1] 12505.717  4159.200  4025.308  4024.579  2767.834  2722.902  2722.899

# N) model 14: GLMM, hierarchical Poisson 10 covariates (regions, region dummy, interaction region dummy-swab,interaction region dummy-lagged hosp., interaction lockdown-time, interaction lockdown-times^2, interaction lockdown-times^3, swabs, lagged swabs, lagged hosp.)
# N.B. distinct starting points & times are scaled

mod.hier.5.int.3 <- stan_glmer(resp ~ (1|codice_regione)+ as.factor(dummy_reg_rep)*tamp_lag +
                               as.factor(dummy_reg_rep)*tamponi + as.factor(dummy_reg_rep)*osp_lag +
                               #as.factor(dummy_reg_rep) 
                               + lockdown_rep*I(times_scaled) + lockdown_rep*I(times_scaled^2) + 
                                lockdown_rep*I(times_scaled^3) + tamp_lag + osp_lag + tamponi,
                               data=regioni_ristr, family=poisson, iter=n.iter, cores=4)

loo.hier.5.int.3 <- loo(mod.hier.5.int.3)$estimates[3,1]
loo.hier.5.int.3
#_______________________________________________________________________________
# single-region GLMs

# MARCHE
regioni_marche <- regioni_c%>%filter(denominazione_regione=="Marche")
resp_marche <- regioni_marche$totale_ospedalizzati
times_marche <- c(1:l_times)
tamp_lag_marche <- tamp_lag[regioni$denominazione_regione=="Marche"]
osp_lag_marche <- osp_lag[regioni$denominazione_regione=="Marche"]
tamponi_marche <- regioni$tamponi[regioni$denominazione_regione=="Marche"]

mod.glm.marche <- stan_glm(resp_marche ~ #as.factor(dummy_reg_rep) +
                           lockdown*I(times_marche) +lockdown*I(times_marche^2) + 
                           lockdown*I(times_marche^3) + tamp_lag_marche + osp_lag_marche,
                           data=regioni_marche, family=poisson, iter=n.iter, cores=4)

# FRIULI VENEZIA GIULIA
regioni_fvg <- regioni_c%>%filter(denominazione_regione=="Friuli Venezia Giulia")
resp_fvg <- regioni_fvg$totale_ospedalizzati
times_fvg <- c(1:l_times)
tamp_lag_fvg <- tamp_lag[regioni$denominazione_regione=="Friuli Venezia Giulia"]
osp_lag_fvg <- osp_lag[regioni$denominazione_regione=="Friuli Venezia Giulia"]
tamponi_fvg <- regioni$tamponi[regioni$denominazione_regione=="Friuli Venezia Giulia"]

mod.glm.fvg <- stan_glm(resp_fvg ~ #as.factor(dummy_reg_rep) +
                        lockdown*I(times_fvg) + lockdown*I(times_fvg^2) + 
                        lockdown*I(times_fvg^3) + tamp_lag_fvg + osp_lag_fvg,
                        data=regioni_fvg, family=poisson, iter=n.iter, cores=4)
#_______________________________________________________________________________
# graphical model comparisons (LOOIC, DIC) for the hierarchical models

loo.v <- c(loo.hier.1, loo.hier.2, loo.hier.3, loo.hier.4, loo.hier.4.int, loo.hier.4.int.2, loo.hier.5, loo.hier.5.int, loo.hier.5.int.2)
dic.v <- c(dic.hier.1, dic.hier.2, dic.hier.3, dic.hier.4, dic.hier.4.int, dic.hier.4.int.2, dic.hier.5, dic.hier.5.int, dic.hier.5.int.2)
sort.loo.v <- sort.int(loo.v, index.return = TRUE)$x
sort.dic.v <- sort.int(dic.v, index.return = TRUE)$x

par(xaxt="n", mfrow=c(1,2))
plot(sort.loo.v, type="b", xlab="", ylab="LOOIC")
par(xaxt="s")
axis(1, c(1:9), c("hier1", "hier2", "hier3", "hier4", "hier4int", "hier4int2", "hier5", "hier5int", "hier5int2")[sort.int(loo.v, index.return = TRUE)$ix], las=2)
par(xaxt="n")
plot(sort.dic.v, type="b", xlab="", ylab="DIC")
par(xaxt="s")
axis(1, c(1:9), c("hier1", "hier2", "hier3", "hier4", "hier4int", "hier4int2", "hier5", "hier5int", "hier5int2")[sort.int(dic.v, index.return = TRUE)$ix], las=2)
#_______________________________________________________________________________
# model posterior estimates just for the best model
# labeling
alpha_names <- paste0("alpha[", n_reg:1, "]")
beta_names <- paste0("beta[", 5:1, "]")
new_names <- c(expression(sigma[alpha]), alpha_names, beta_names, expression(mu))
posterior <- as.array(mod.hier.5.int)
mcmc_intervals(posterior)+ xaxis_text(on =TRUE, size=rel(1.9)) + yaxis_text(on =TRUE, size=rel(1.9))+ scale_y_discrete(labels = rev((parse(text = new_names))))+ ggtitle("Parameter estimation")+
theme(plot.title = element_text(hjust = 0.5, size =rel(2)))
#_______________________________________________________________________________
#compare the empirical distribution of true data with empirical distribution of replications
color_scheme_set("blue")
ppc_ecdf_overlay(y=resp, yrep=posterior_predict(mod.hier.5.int, draws=200)) # ecdf
ppc_stat(y=resp, yrep=posterior_predict(mod.hier.5.int, draws=200), stat="mean") # mean
ppc_stat(y=resp, yrep=posterior_predict(mod.hier.5.int, draws=200), stat="sd") # sd
#_______________________________________________________________________________
# make predictions for the 4 next days (via local polynomial fitting on time and squared time)

# data manipulation for out-of-sample
mod <- mod.hier.4.int
data <- regioni_ristr
prev <- 4 # we are going to predict 4 data points
times <- as.factor(data$data)
levels(times) <- c(1:l_times)
times <- as.numeric(times)
new_times <- c((l_times):(l_times+prev))
pred <- posterior_predict(mod, type ="response")
starting <- rep(l_times, n_reg)

for (j in 1:prev){
  times_out <- c(starting, rep(l_times+j, n_reg)) 
  starting <- times_out
}

times_tot <- c(times, times_out)
centered_times_tot <- times_tot-mean(times)
stringa_regioni <- data$denominazione_regione[502:522]
stringa_codici <- data$codice_regione[502:522]
tamponi <- regioni_ristr$tamponi
new_tamponi <- matrix(NA, n_reg, prev+1)
se_new <- matrix(NA, n_reg, prev+1)

#loess prediction for the swabs
par(mfrow=c(2,3))
for (j in 1:n_reg){
  obj <- loess(tamponi[denominazione_regione==stringa_regioni[j]] ~
               times[denominazione_regione==stringa_regioni[j]]+
               I(times[denominazione_regione==stringa_regioni[j]])^2,
               data = regioni_ristr, control = loess.control(surface = "direct"))
  new_tamponi[j, 1] <- obj$fitted[length(obj$fitted)]
  
  for (t in (new_times-1)){
    fit.loess <- predict(obj, newdata=data.frame(times = t+1, denominazione_regione=stringa_regioni[j]), se=TRUE)
    new_tamponi[j, t-l_times+2] <- fit.loess$fit
    se_new[j, 1] <- predict(obj, newdata=data.frame(times=l_times, denominazione_regione=stringa_regioni[j]), se=TRUE)$se
    se_new[j,t-l_times+2] <- fit.loess$se
  }
  plot(c(times[regioni_ristr$denominazione_regione==stringa_regioni[j]], new_times[2]:(l_times + prev)),  
       c(tamponi[regioni_ristr$denominazione_regione==stringa_regioni[j]], rep(NA, length(2: (prev+1)))),
       xlab="Days", ylab="swabs", main=stringa_regioni[j], ylim=c(0, new_tamponi[j, prev+1]+0.1*new_tamponi[j, prev+1]))
  points(new_times[1]:(l_times + prev), new_tamponi[j, 1:(prev+1)], col="red", lwd=2)
}

# prediction of n° of total hospitalized (21 regions x next 4 days) 
newdata <- data.frame(times_scaled = centered_times_tot[(length(times)+1):(length(times)+length(times_out))],
                      denominazione_regione = rep(stringa_regioni, prev+1), codice_regione = rep(stringa_codici, prev+1),
                      lockdown_rep = rep(1, n_reg*(prev+1)), tamponi = as.vector(new_tamponi))

# set up function for plotting predictions
regioni_plot <- function(mod, data, prev, newdata, new_points){
  pred.out <- posterior_predict(mod, newdata=newdata, type="response")
  pred_0025 <- apply(pred, 2, function(x) quantile(x,0.025))
  pred_025 <- apply(pred, 2, function(x) quantile(x,0.25))
  pred_med <- apply(pred, 2, function(x) quantile(x,0.5))
  pred_075 <- apply(pred, 2, function(x) quantile(x,0.75))
  pred_0975 <- apply(pred, 2, function(x) quantile(x,0.975))
  
  pred_0025_out <- apply(pred.out, 2, function(x) quantile(x,0.025))
  pred_025_out <- apply(pred.out, 2, function(x) quantile(x,0.25))
  pred_med_out <- apply(pred.out, 2, function(x) quantile(x,0.5))
  pred_075_out <- apply(pred.out, 2, function(x) quantile(x,0.75))
  pred_0975_out <- apply(pred.out, 2, function(x) quantile(x,0.975))
  
  if (missing(new_points)){new_points <- rep(NA, length(new_times)*(n_reg))
    }else{ new_points <- c(rep(NA, n_reg), new_points)
      }
  
  df.regioni <- data.frame(y = c(resp, new_points), times = c(times, times_out), 
  med = c(pred_med, rep(NA, length(new_times)*n_reg)),
  lo = c(pred_0025, rep(NA, length(new_times)*n_reg)),
  lo.2 = c(pred_025, rep(NA, length(new_times)*n_reg)),
  hi = c(pred_0975, rep(NA, length(new_times)*n_reg)), 
  hi.2 = c(pred_075, rep(NA, length(new_times)*n_reg)),
  reg = unlist(list(data$denominazione_regione, rep(stringa_regioni, prev+1))),
  med.out = c(rep(NA, length(times)), pred_med_out),
  lo.out = c(rep(NA, length(times)), pred_0025_out),
  lo.2.out = c(rep(NA, length(times)), pred_025_out),
  hi.out= c(rep(NA, length(times)), pred_0975_out),
  hi.2.out = c(rep(NA, length(times)),pred_075_out)
  )
  
  if (prev <=7){
    table_reg <- matrix(NA, n_reg, prev)
    str_days <- c("30 mar", "31 mar", "1 apr", "2 apr", "3 apr", "4 apr", "5 apr", "6 apr")
    options(knitr.table.format = "html")
    for (j in 1:prev){
      table_reg[,j] <-paste(round(pred_med_out[(n_reg*j+1):(n_reg*(j+1))])," (",
              round(pred_0025_out[(n_reg*j+1):(n_reg*(j+1))]), ", ",
              round(pred_0975_out[(n_reg*j+1):(n_reg*(j+1))]), ")", sep="")
    }
    dimnames(table_reg) <- list(stringa_regioni, str_days[1:prev])
    # kable(as.data.frame(table_reg), "html",
    #       caption =
    #         "Predicted number of total hospitalized (in
    #         parentheses, the lower and the upper bound)")%>%
    #   kable_styling(bootstrap_options = "striped",
    #                 full_width = F)
    ospt <- write.csv(table_reg, "osp30mar.csv") #output 1
  }
  # plot
  ggplot(df.regioni, aes(times, y))+geom_ribbon(aes(x=times, ymin=lo, ymax=hi, group=1),
                data=df.regioni, fill = color_scheme_get("blue")[[2]])+
    geom_ribbon(aes(x=times, ymin=lo.out, ymax=hi.out, group=1),
                data=df.regioni, fill = color_scheme_get("red")[[2]])+
    geom_line(aes(x= times, y= med), data=df.regioni,
              color = color_scheme_get("blue")[[4]], size =1.1)+
    geom_line(aes(x= times, y= med.out), data=df.regioni,
              color = color_scheme_get("red")[[4]], size =1.1)+
    geom_point(aes(x = times, y = y))+ xlab("Days")+ ylab("Total Hospitalized")+
    scale_x_discrete(limit = c( 8, 15, 22, 29, 36 ), 
                     labels = c( "2-3", "9-3", "16-3", "23-3", "30-3"))+
    facet_wrap("reg", scales ="free")+
    theme(strip.text.x = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(face="bold", color="#993333", angle=45, size =9),
          axis.title.x = element_text(size=22),
          axis.title.y = element_text(size = 22))

  # return(list(pred_0025_out = pred_0025_out,
  #             pred_med_out = pred_med_out,
  #             pred_0975_out = pred_0975_out))
}
# 4-days predictions plot
regioni_plot(mod, data, prev, newdata) #call the function
#_______________________________________________________________________________
# make predictions of total hospitalized for the 30 next days
# data for predictions
newdata <- data.frame(times_scaled=centered_times_tot[(length(times)+1):(length(times)+length(times_out))],
                      denominazione_regione=rep(stringa_regioni, prev+1),
                      codice_regione = rep(stringa_codici, prev+1),
                      lockdown_rep = rep(1, n_reg*(prev+1)),
                      tamponi = as.vector(new_tamponi))

regioni_plot(mod, data, prev, newdata)
#_______________________________________________________________________________
# modeling intensive care units (ICU)
icu <- regioni_ristr$terapia_intensiva

# A) model 1: GLMM, hierarchical Poisson 3 covariates (regions, interaction lockdown-time, times^2)
# N.B. distinct starting points & times are scaled
mod.hier.icu.2.int <- stan_glmer(icu ~ lockdown_rep*times_scaled + (1|codice_regione) + I(times_scaled^2),
                                  data=regioni_ristr, family=poisson, iter=n.iter)
print(mod.hier.icu.2.int)
#_______________________________________________________________________________
# model posterior estimates just for the best model (ICUs)
# labeling
alpha_names <- paste0("alpha[", n_reg:1, "]")
beta_names <- paste0("beta[", 4:1, "]")
new_names <- c(expression(sigma[alpha]), alpha_names, beta_names, expression(mu))
posterior2 <- as.array(mod.hier.icu.2.int)
mcmc_intervals(posterior2)+ xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+ scale_y_discrete(labels = rev((parse(text = new_names))))+
  ggtitle("Parameter estimation")+ theme(plot.title = element_text(hjust = 0.5, size =rel(2)))

###END######END######END######END######END######END######END######END######END######END###