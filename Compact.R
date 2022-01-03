# __CASE: COVID-19 - Veneto's ICUs__ ------------------------------------------
# Authors: Renzo Testa, Lorenzo Cabriel, Mattia Pividori
# Preliminary Notes: LM (with log-transform), GLM (Poisson, Quasi-Poisson with offset)
                  # GAM?, Bayes? etc..

# A) LIBRARIES ----------------------------------------------------------------
library(caret)
library(DataExplorer)
library(ggplot2)
library(performance)
library(PerformanceAnalytics)
library(splines)
library(tibble)
library(tidymodels)
library(tidyverse)
library(vars)
library(vip)

# B) HELPER FUNCTIONS ---------------------------------------------------------

# 1.Generates prediction tables (true data & predictions)
prediction_table <- function(fitted_model, input_data, is_tidy) {
  if (is_tidy == TRUE) {
    result <- fitted_model %>%
              stats::predict(input_data) %>%
              dplyr::rename(pred = .pred) %>%
              dplyr::mutate(actual = input_data$ICU,
                            pred_real = pred^2, # use if we decide to transform Y (see F.1)
                            actual_real = actual^2)  # use if we decide to transform Y (see F.1)
  } else {result <- fitted_model %>% 
                    stats::predict(input_data) %>%
                    tibble::as_tibble_col(column_name = "pred") %>%
                    dplyr::mutate(actual = input_data$ICU,
                                  pred_real = pred^2, # use if we decide to transform Y (see F.1)
                                  actual_real = actual^2) # use if we decide to transform Y (see F.1)
  }
}

# 2.Extracts performance metrics from models
pull_metrics <- function(input_tibble) {
  metrics <- yardstick::metric_set(rmse, rsq, mae)
  
  tranformed_result <- metrics(input_tibble, truth=actual, estimate=pred) %>%
    dplyr::pull(.estimate)
  
  real_result <- metrics(input_tibble, truth=actual_real, estimate=pred_real) %>%
    dplyr::pull(.estimate)
  
  exit <- cbind(transf_metrics = tranformed_result, real_metrics = real_result) %>%
    as.data.frame() %>% magrittr::set_rownames(c("rmse", "rsq", "mae"))
}

# 3.Save model coefficients for a fitted model object to re-apply to k-folds
get_model <- function(x) {
  hardhat::extract_fit_parsnip(x) %>% generics::tidy()
}

# C) IMPORT RAW DATA ------------------------------------------------------
setwd("C:/Users/m/Desktop/Project/smds_2021_final_project")

# import data from different sources
regioni <- readr::read_csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv")

vaccini <- readr::read_csv("https://raw.githubusercontent.com/italia/covid19-opendata-vaccini/master/dati/somministrazioni-vaccini-summary-latest.csv")


# D) DATA WRANGLING -------------------------------------------------------
regioni %>% dim()
vaccini %>% dim()

veneto_covid <- regioni %>%
                dplyr::select(-c(stato, denominazione_regione, 
                                 lat, long, note, note_test, note_casi,
                                 codice_nuts_1, codice_nuts_2)) %>%
                dplyr::mutate(data = as.Date(data)) %>%
                dplyr::rename(date = data,
                              ICU = terapia_intensiva,
                              hospitalized_total = totale_ospedalizzati,
                              swabs = tamponi,
                              total_cases = totale_casi,
                              deaths = deceduti,
                              self_isolation = isolamento_domiciliare) %>%
                dplyr::filter(codice_regione=="05")

veneto_vaccini <- vaccini %>%                
                  dplyr::select(-c(area, dose_addizionale_booster, codice_NUTS1,
                                   codice_NUTS2, nome_area)) %>%
                  dplyr::rename(date = data_somministrazione,
                                vaccinated_total = totale,
                                vaccinated_males = sesso_maschile,
                                vaccinated_females = sesso_femminile,
                                shot_first = prima_dose,
                                shot_second = seconda_dose,
                                already_infected = pregressa_infezione) %>%
                  dplyr::filter(codice_regione_ISTAT=="5")

# merge all data & add residents
veneto_alldata <- veneto_covid %>%
                  dplyr::left_join(veneto_vaccini, by=c("date"="date")) %>%
                  dplyr::select(-c(codice_regione, codice_regione_ISTAT)) %>%
                  tibble::add_column(residents = 4879133) #pre-Covid residents at 01/01/2020


#TODO: discuss is worth filtering on dates
dplyr::filter(cleaned_dataset, between(data, as.Date("2020-10-01"), as.Date("2021-02-01")))
dplyr::filter(cleaned_dataset, between(data, as.Date("2020-02-01"), as.Date("2021-02-15")))


# E) EXPLORATORY DATA ANALYSIS (EDA) --------------------------------------
str(veneto_alldata) #discover data-types, see if R already sees some categorical columns as Factors
length(unique(bike_all$holiday)) # check which columns are labeled as 'num' even if are categorical ('Factor')->convert them!

summary(veneto_alldata) #descriptive stats for each column

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

# F) DATA PREPARATION (all data) -------------------------------------------

## F.1) Checking transforms ------------------------------------------------
# we wish a normal response variable (ICUs) so we test which transform to use:
# original response 
PerformanceAnalytics::skewness(veneto_alldata$ICU) #[1] 0.9354436 positive skew
PerformanceAnalytics::kurtosis(veneto_alldata$ICU) #[1] -0.4690721 (excess) platykurtic

# try natural log
PerformanceAnalytics::skewness(log(veneto_alldata$ICU)) #[1] Nan
PerformanceAnalytics::kurtosis(log(veneto_alldata$ICU)) #[1] Nan

# try log10
PerformanceAnalytics::skewness(log10(veneto_alldata$ICU)) #[1] Nan
PerformanceAnalytics::kurtosis(log10(veneto_alldata$ICU)) #[1] Nan

# try log(1+x)
PerformanceAnalytics::skewness(log1p(veneto_alldata$ICU)) #[1] -0.5695241
PerformanceAnalytics::kurtosis(log1p(veneto_alldata$ICU)) #[1] -0.5856917

# try square root
PerformanceAnalytics::skewness(sqrt(veneto_alldata$ICU)) #[1] 0.3493079
PerformanceAnalytics::kurtosis(sqrt(veneto_alldata$ICU)) #[1] -1.135296

# try cubic root
PerformanceAnalytics::skewness(veneto_alldata$ICU^(1/3)) #[1] 0.01831247
PerformanceAnalytics::kurtosis(veneto_alldata$ICU^(1/3)) #[1] -1.070609

## F.2) Categorical variables as factors ----------------------------------
# categorical variables are converted to factors and ordered when needed
# bike_all$season <- factor(bike_all$season,
#                           levels = c(1, 2, 3, 4),
#                           labels = c("spring", "summer", "autumn", "winter"))
# 
# bike_all$workingday <- factor(bike_all$workingday,
#                               levels = c(0, 1), labels = c(FALSE, TRUE))
# 
# bike_all$weather <- factor(bike_all$weather,
#                            levels = c(1, 2, 3, 4),
#                            labels = c("clear", "cloudy", "rainy", "heavy rain"),
#                            ordered = TRUE)

# transform predictors in log-scale and add
cleaned_dataset$log_ICU<-log(cleaned_dataset$terapia_intensiva)
cleaned_dataset$log_hospitalized<-log(cleaned_dataset$totale_ospedalizzati)
cleaned_dataset$log_swabs<-log(cleaned_dataset$tamponi)




## F.3) Set up a train/test split -----------------------------------------

set.seed(25) # fix seed to make reproducible results

split <- veneto_alldata %>% rsample::initial_split(prop = 0.8) # we can change % here
train_data <- split %>% rsample::training()
test_data <- split %>% rsample::testing()

train_data %>% dim() # check if is 80% of data points [1] 0.79941
test_data %>% dim() # checks if is 20% of data points [1] 0.20059

## F.4) Set up k-fold(s) cross validation ---------------------------------

train_cv <- train_data %>% rsample::vfold_cv(v = 10, repeats = 2) # change k-folds here and repetitions


# G) FEATURE ENGINEERING (separate train/test data) -----------------------
# using recipe library to keep scalable steps

recipe_lm <- recipes::recipe(ICU ~ swabs, data = train_data) %>%
             #recipes::step_rm(year, month, weekday) %>%
             #recipes::step_date(date) %>%
             #recipes::step_corr(all_numeric(), threshold = 0.8) %>%
             #recipes::step_poly(all_predictors()) %>% 
             #recipes::step_interact(~ all_predictors():all_predictors())
             #recipes::step_ns(days, deg_free = tune::tune("days df")) %>% 
             recipes::step_log(swabs, base = 10) %>% 
             recipes::step_dummy(all_nominal_predictors())

#add time as increment --> inserisci qui sopra come step
cleaned_dataset_train$days <- seq(1, length(unique(cleaned_dataset_train$data)))



#discover best lag
best_lag <- VARselect(cleaned_dataset$terapia_intensiva)#find best lag for Y variable
n <- best_lag$selection[1]
#auto-regressive distributed lag (ADL) 


# H) DATA MODELING ---------------------------------------------------------
# we used tidymodels & caret for all modeling section

cores <- parallel::detectCores() # get cores for parallelization

## H.1) Model 1: LM with log-transforms ------------------------------------
lm_engine <- parsnip::linear_reg() %>% 
             parsnip::set_engine("lm") %>%
             parsnip::set_mode("regression")

lm_workflow <- workflows::workflow() %>%
               #workflows::remove_variables() %>% 
               workflows::add_recipe(recipe_lm) %>%
               workflows::add_model(lm_engine)

# optional: create a grid and tune spline params with cv
# https://tune.tidymodels.org/articles/getting_started.html
# https://stackoverflow.com/questions/65274224/tidymodels-tunable-models-involving-10-fold-cross-validation-using-the-function
#recap which are the tuning params in recipe
lm_param <- recipe_lm %>% 
            dials::parameters() %>% 
            recipes::update(`days df` = spline_degree())

spline_grid <- grid_max_entropy(lm_param, size = 10)

lm_tune <- tune::tune_grid(lm_workflow,
                           resamples = train_cv,
                           grid = spline_grid,
                           metrics = yardstick::metric_set(rmse))
#_______________________________________________________________________________
# fit just the best model 
lm_train_fit <- fit(lm_workflow, train_data)

# store results for evaluation
lm_results <- lm_train_fit %>% 
              hardhat::extract_fit_engine()

# quick peek at results
lm_results %>% summary()
generics::glance(lm_train_fit)

# check if results are in line with tidymodels: OK!
fit1 <- stats::lm(train_data$ICU ~ log10(train_data$swabs))
summary(fit1)

## H.2) Model 2: GLM Poisson with intercept + cubic spline on time ---------

fit2 <- glm(terapia_intensiva ~ ns(days, 3), data=cleaned_dataset, family=poisson(link="log"), offset=log(residenti))
summary(fit2) #residual deviance[168.17]>df[120] signals overdispersion!
an2 <- anova(fit2, test="Chisq") #check which covariates to keep

predict2 <- predict(fit2, list(Days=seq(15,1)))

## H.3) Model 3: GLM Quasi-Poisson with intercept + cubic spline on time ----
fit3 <- glm(terapia_intensiva ~ ns(days, 3), data=cleaned_dataset, family=quasipoisson(link="log"), offset=log(residenti))
summary(fit3)
an3 <- anova(fit3, test="F")


## H.4) Model 4: ARIMA (1,1,1) ---------------------------------------------
## H.5) Model 5: ARIMA with best lag and cross correlation best lags -------


# I) DATA PREDICTION & EVALUATION ---------------------------------------
# here we valuate MAE/MSE/RMSE/R_2 on predictions for each model to choose best

## I.1) Model 1: LM with log-transforms ------------------------------------

# a) predict on train set
lm_train_res <- prediction_table(lm_train_fit, train_data, TRUE)

# b) get metrics on train set
pull_metrics(lm_train_res) %>% print()

# c) predict on validation set
lm_validation_res <- lm_workflow %>%  
                     tune::fit_resamples(resamples = train_cv,
                     control = tune::control_resamples(extract = get_model),
                     metrics = yardstick::metric_set(rmse, rsq, mae))

# d) get metrics on validation set (resampling-estimates average error over K folds x R repetitions)
tune::collect_metrics(lm_validation_res, summarize = TRUE)

# e) predict on test set
lm_test_res <- prediction_table(lm_train_fit, test_data, TRUE)

# f) get metrics on test set   
pull_metrics(lm_test_res) %>% print()

# f.bis) alternative test error check via workflow
lm_workflow %>% 
  last_fit(split = split) %>% 
  collect_metrics()

# g) check predictors importance
vip::vip(fit1) #TODO: fare un primo lm con tutti i predittori e le interazioni e poi scremo
vip::vip(lm_results)

# h) plot model vs true data


# i) plot residuals
lm_results %>% performance::check_model()
graphics::par(mfrow = c(2,2))
lm_results %>% plot(pch = 16,
              col = '#006EA1')


## I.2) Model 2: GLM Poisson with intercept + cubic spline on time ---------

# J) COMPARISONS & CONCLUSIONS -------------------------------------------
rbind(base_train_rmse, base_test_rmse,
      tree_tidy_train_rmse, tree_tidy_test_rmse,
      tree_caret_train_rmse, tree_caret_test_rmse)
#lsf.str("package:dplyr") #to list all functions in package
#environmentName(environment(select)) #to get package name from function
