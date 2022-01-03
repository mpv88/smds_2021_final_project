#lsf.str("package:dplyr") #to list all functions in package
#environmentName(environment(select)) #to get package name from function
############
# Covid 19 #
############

# required packages
library(ggplot2)
library(splines)
library(tidyverse)
library(PerformanceAnalytics)
library(vars)

#_______________________________________________________________________________
setwd("C:/Users/m/Desktop/Project/smds_2021_final_project")

# import data from github
regioni <- read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv")
regioni$data <- as.Date(regioni$data) #get rid of timestamps
veneto_data <- dplyr::filter(regioni, codice_regione==5) #filters region
#write.csv(veneto_data,"veneto.csv", row.names=FALSE)

vaccini <- read.csv("https://raw.githubusercontent.com/italia/covid19-opendata-vaccini/master/dati/somministrazioni-vaccini-summary-latest.csv")
vaccini$data_somministrazione <- as.Date(vaccini$data_somministrazione)
vaccini_veneto <- dplyr::filter(vaccini, codice_regione_ISTAT==5)
#write.csv(vaccini_veneto,"vaccini.csv", row.names=FALSE)

#add residents
veneto_data$residenti <- 4879133 # data at 01/01/2020 (pre-Covid) http://dati.istat.it/Index.aspx?QueryId=18460
#_______________________________________________________________________________
# data manipulation

#join, filter, drop, merge, clean, lag, log-transform etc..
cleaned_dataset <- dplyr::left_join(x=veneto_data, y=vaccini_veneto, by=c("data"="data_somministrazione"), suffix=c(".x",".y"))
cleaned_dataset <- select(dataset, -c(stato, codice_regione, denominazione_regione, lat, long, note, note_test, note_casi,
                                      codice_nuts_1, codice_nuts_2, area, dose_addizionale_booster, codice_NUTS1, codice_NUTS2,
                                      codice_regione_ISTAT, nome_area)) # drop useless cols

# transform predictors in log-scale and add
cleaned_dataset$log_ICU<-log(cleaned_dataset$terapia_intensiva)
cleaned_dataset$log_hospitalized<-log(cleaned_dataset$totale_ospedalizzati)
cleaned_dataset$log_swabs<-log(cleaned_dataset$tamponi)

# split train/test dataset
cleaned_dataset_train <- dplyr::filter(cleaned_dataset, between(data, as.Date("2020-10-01"), as.Date("2021-02-01")))
cleaned_dataset_test <- dplyr::filter(cleaned_dataset, between(data, as.Date("2020-02-01"), as.Date("2021-02-15")))

#add time as increment
cleaned_dataset_train$days <- seq(1, length(unique(cleaned_dataset_train$data))) %>% relocate(days, .after=data)
cleaned_dataset_train$log_days<-log(cleaned_dataset_train$days)
cleaned_dataset_test$days <- seq(1, length(unique(cleaned_dataset_test$data))) %>% relocate(days, .after=data)
cleaned_dataset_test$log_days<-log(cleaned_dataset_test$days)

#discover best lag
best_lag <- VARselect(cleaned_dataset$terapia_intensiva)#find best lag for Y variable
n <- best_lag$selection[1]
#auto-regressive distributed lag (ADL) 
#_______________________________________________________________________________
# data visualization
str(cleaned_dataset) #discover data-types
summary(cleaned_dataset) #descriptive stats for each column

# observe distribution of Y and some covariate candidates
par(mfrow=c(1,3))
hist(cleaned_dataset$terapia_intensiva, probability=TRUE, breaks=15)
lines(density(car_price$price), col=4)
hist(cleaned_dataset$totale_ospedalizzati, probability=TRUE, breaks=15)
hist(cleaned_dataset$totale_casi, probability=TRUE, breaks=15)

# plot response variable vs time
ggplot(cleaned_dataset, aes(x=data, y=terapia_intensiva))+
  geom_point(aes(x=data, y=terapia_intensiva))+
  labs(title="ICUs time evolution",
       x="time (days)", y="intensive care units (ICUs")

# correlation matrix
round(cor(cleaned_dataset), 4)

# scatterplot matrix
chart.Correlation(cleaned_dataset[,c("terapia_intensiva", "totale_ospedalizzati", "tamponi", "casi_testati")])

#_______________________________________________________________________________
#MODEL 1: LM with log-covariates
fit1 <- lm(terapia_intensiva ~ ns(days, 3))


#check residuals
par(mfrow = c(2,2))
plot(fit1)
#_______________________________________________________________________________
#MODEL 2: GLM Poisson with intercept + cubic spline on time
fit2 <- glm(terapia_intensiva ~ ns(days, 3), data=cleaned_dataset, family=poisson(link="log"), offset=log(residenti))
summary(fit2) #residual deviance[168.17]>df[120] signals overdispersion!
an2 <- anova(fit2, test="Chisq") #check which covariates to keep

predict2 <- predict(fit2, list(Days=seq(15,1)))
#_______________________________________________________________________________
#MODEL 3: GLM Quasi-Poisson with intercept + cubic spline on time
fit3 <- glm(terapia_intensiva ~ ns(days, 3), data=cleaned_dataset, family=quasipoisson(link="log"), offset=log(residenti))
summary(fit3)
an3 <- anova(fit3, test="F")


#_______________________________________________________________________________
#MODEL 4: ARIMA with best lag
#_______________________________________________________________________________
#MODEL 5: ARIMA with best lag and cross correlation best lags

#Evaluate MAE/MSE on predictions

###END######END######END######END######END######END######END######END######END######END###