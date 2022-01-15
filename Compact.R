# __CASE: COVID-19 - Veneto's ICUs__ ------------------------------------------
# Authors: Renzo Testa, Lorenzo Cabriel, Mattia Pividori
# Preliminary Notes: LM (with transforms), GLM (Poisson & Quasi-Poisson with offset),
#                       ARIMA etc..
#
# A) LIBRARIES ----------------------------------------------------------------
library(DataExplorer)
library(ggplot2)
library(lubridate)
library(MASS)
library(modeltime)
library(performance)
library(PerformanceAnalytics)
library(splines)
library(tibble)
library(tidymodels)
library(tidyverse)
library(vars)
library(vip)

# B) HELPER FUNCTIONS ---------------------------------------------------------

# 1.Generates prediction tables (true data & predictions for LM & GLM)
prediction_table <- function(fitted_model, input_data, is_delta) {
  if (is_delta == TRUE) {
    result <- fitted_model %>%
              stats::predict(input_data) %>%
              dplyr::rename(pred = .pred) %>%
              dplyr::mutate(actual = input_data$delta_ICU,
                            pred_real = pred, # use if we decide to transform Y (see F.1)
                            actual_real = actual)  # use if we decide to transform Y (see F.1)
  } else {result <- fitted_model %>% 
                    stats::predict(input_data) %>% # type = "response"
                    tibble::as_tibble_col(column_name = "pred") %>%
                    dplyr::mutate(actual = log(input_data$ICU),
                                  pred_real = exp(pred), # use if we decide to transform Y (see F.1)
                                  actual_real = exp(actual)) # use if we decide to transform Y (see F.1)
  }
}

# 2.Generates calibration tables (true data & predictions for ARIMAs)
calibration_table <- function(fitted_model, input_data) {
    result <- fitted_model %>%
              modeltime::modeltime_calibrate(input_data) %>%
              purrr::pluck(.,".calibration_data") %>% 
              purrr::pluck(., 1) %>%
              dplyr::mutate(.residuals = NULL,
                            actual_real = .actual,
                            pred_real = .prediction) %>% 
              dplyr::rename(actual = .actual,
                            pred = .prediction)
 }

# 3.Extracts performance metrics from models
pull_metrics <- function(input_tibble) {
  metrics <- yardstick::metric_set(yardstick::rmse, yardstick::rsq, yardstick::mae,
                                   yardstick::mape, yardstick::mase, yardstick::smape)
  
  tranformed_result <- metrics(input_tibble, truth = actual, estimate = pred) %>%
                       dplyr::pull(.estimate)
  
  real_result <- metrics(input_tibble, truth = actual_real, estimate = pred_real) %>%
                 dplyr::pull(.estimate)
  
  exit <- cbind(transf_metrics = tranformed_result, real_metrics = real_result) %>%
          as.data.frame() %>% magrittr::set_rownames(c("rmse", "rsq", "mae", 
                                                       "mape", "mase", "smape"))
}

# 4.Save model coefficients for a fitted model object to re-apply to k-folds
get_model <- function(x) {
  hardhat::extract_fit_parsnip(x) %>% generics::tidy()
}

# 5.Plot fitted model estimates vs. real-world data
plot_fitted <- function(prediction_table) {
  ggplot2::ggplot(prediction_table, ggplot2::aes(x = actual, y = pred)) + 
  ggplot2::geom_abline(lty = 2) + 
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::labs(x = "Official ICUs registered", y = "Predicted ICUs") +
  tune::coord_obs_pred()
}

# 6.Plot fitted estimates series vs. real-world series
plot_series <- function(prediction_table) {
  df <- prediction_table %>% dplyr::mutate(residuals_real = actual_real-pred_real)
                             
  df1 <- df %>% dplyr::select(date, actual_real, pred_real) %>%
                tidyr::gather(key = "variable", value = "value", -date)
     
  p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = date, y = value)) + 
        ggplot2::geom_line(aes(color = variable, linetype = variable)) + 
        ggplot2::scale_color_manual(values = c("darkred", "steelblue")) +
        ggplot2::labs(x = NULL, y = "true vs. predicted ICUs") +
        ggplot2::theme(legend.position = "top") 
  
  df2 <- df %>% dplyr::select(date, residuals_real) %>%
                tidyr::gather(key = "variable", value = "value", -date)
  
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = date, y = value)) + 
        ggplot2::geom_bar(stat = "identity") + ggplot2::theme_minimal() + 
        ggplot2::theme(axis.title.x = element_blank(),
                       axis.text.x = element_text(angle=90)) +
        ggplot2::labs(y = "residuals")
  
  gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(2, 1),
                          top = grid::textGrob("ICUs series",
                          gp = grid::gpar(fontsize=20,font=3)),
                          bottom = grid::textGrob("time (days)"))
}

# 7.Plot residuals of fitted model
plot_residuals <- function(model_results) {
  graphics::par(mfrow = c(2,2))
  tryCatch(
  plot1 <- model_results %>% plot(pch = 16, col = "#006EA1")
  )
  plot2 <- model_results %>% performance::check_model()
  list(plot1, plot2)
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
                              self_isolation = isolamento_domiciliare,
                              positives_total = totale_positivi) %>%
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
                  #tibble::column_to_rownames("date") %>%
                  dplyr::select(-c(codice_regione, codice_regione_ISTAT,
                                   already_infected,shot_second,shot_first,
                                   vaccinated_females, vaccinated_males, tamponi_test_antigenico_rapido,
                                   tamponi_test_molecolare, totale_positivi_test_antigenico_rapido,
                                   totale_positivi_test_molecolare, ingressi_terapia_intensiva,
                                   casi_testati,casi_da_screening, casi_da_sospetto_diagnostico,
                                   dimessi_guariti, variazione_totale_positivi, nuovi_positivi,
                                   ricoverati_con_sintomi)) %>%
                  tibble::add_column(residents = 4879133) # pre-Covid residents at 01/01/2020

veneto_alldata %>% DataExplorer::introduce()
veneto_alldata %>% DataExplorer::plot_intro()

# E) EXPLORATORY DATA ANALYSIS (EDA) --------------------------------------

# 0.Visual dataset structure
veneto_alldata %>% DataExplorer::plot_str() # tree based
veneto_alldata %>% DataExplorer::plot_str(type = "r") # radial
veneto_alldata %>% DataExplorer::introduce() # introductory data table
veneto_alldata %>% DataExplorer::plot_intro() # introductory data plot
veneto_alldata %>% DataExplorer::plot_missing() # visual missing data
veneto_alldata %>% DataExplorer::profile_missing() # missing data profile for additional analysis

# 1.Data Inspection, a collection of methods
veneto_alldata %>% structure() # peek at tibble details + head rows
veneto_alldata %>% glimpse() # glimpse at all columns
veneto_alldata %>% str() # discover data-types, is R already seeing cat cols as Factors?
veneto_alldata %>% summary() # descriptive stats for each col

# unique occurrences of each col, categorical data to be later converted to factors
veneto_alldata %>%
  dplyr::summarise(dplyr::across(tidyselect::everything(), n_distinct))

# 2.Observing Y distribution
# plot response variable vs time
ggplot(veneto_alldata, aes(x = date, y=ICU))+
  geom_point(aes(x=date, y = ICU))+
  labs(title="ICUs time evolution",
       x="time (days)", y = "intensive care units (ICUs)")

# plot distribution, do we need rescaling?
ggplot2::ggplot(veneto_alldata, ggplot2::aes(x = ICU)) + 
  ggplot2::geom_histogram(bins = 50) 
#+ scale_x_log10()
#TODO: here you can try visually which is the best transform 

# 3.Observe covariates distributions

# FOR CATEGORICAL/DISCRETE VARIBLES
veneto_alldata %>% DataExplorer::plot_bar() # bar-charts to get frequency distribution
veneto_alldata %>% DataExplorer::plot_bar(with = "ICU") # bivariate frequency distribution of each discrete feature, by X or Y
#TODO do we have categoricals or we want to add them?

veneto_alldata %>% # visualize just discrete var correlations heatmap
  stats::na.omit() %>% DataExplorer::plot_correlation(type = "d")

veneto_alldata %>% # boxplots, distributions of some numerical X conditional on Y
  dplyr::select(c("ICU", "hospitalized_total", "swabs", "deaths")) %>%
  DataExplorer::plot_boxplot(by = "ICU")

# FOR NUMERICAL/CONTINOUS VARIBLES
n_obs <- 680L # n° obs we want to consider for analysis purposes

veneto_alldata %>% DataExplorer::plot_histogram() # histograms with all X distributions 

veneto_alldata %>% # qq-plot of selected X, on selected n° obs/rows
  dplyr::select(c("ICU", "hospitalized_total", "swabs", "deaths")) %>% 
  DataExplorer::plot_qq(sampled_rows = n_obs)

veneto_alldata %>% # qq-plot of all vars conditioned on 1 variable (better if categorical)
  dplyr::select(c("ICU", "residents")) %>% 
  DataExplorer::plot_qq(by = "ICU", sampled_rows = n_obs)

#TODO do we want to re-scale some vars? here you can try visually which is the best transform
veneto_alldata %>% DataExplorer::update_columns(c("ICU", "swabs"), function(x) log(x)) %>%
                   dplyr::select(c("ICU", "swabs")) %>% 
                   DataExplorer::plot_qq(sampled_rows = n_obs) # for example here we log(ICU, swabs)

veneto_alldata %>% DataExplorer::plot_density() # plot density distribution all vars

veneto_alldata %>% # correlation matrix of selected variables
  dplyr::select(c("ICU", "hospitalized_total", "swabs", "deaths")) %>% 
  stats::cor() %>% round(digits = 4)

veneto_alldata %>% # visualize var correlations heatmap
  dplyr::select(c("ICU", "hospitalized_total", "swabs", "deaths")) %>%
  stats::na.omit() %>% DataExplorer::plot_correlation(maxcat = 5L)

veneto_alldata %>% # scatterplot matrix with distributions
  dplyr::select(c("ICU", "hospitalized_total", "swabs", "deaths")) %>%
  PerformanceAnalytics::chart.Correlation()

veneto_alldata %>% # scatterplots of each single X vs Y
  dplyr::select(c("ICU", "hospitalized_total", "swabs", "deaths")) %>%
  DataExplorer::plot_scatterplot(by = "ICU", sampled_rows = n_obs)

veneto_alldata %>% # perform/visualize PCA on selected features to reduce dimension of problem
  dplyr::select(c("ICU", "hospitalized_total", "swabs", "deaths")) %>%
  stats::na.omit() %>% DataExplorer::plot_prcomp(variance_cap = 0.9, nrow = 2L, ncol = 2L)


# F) DATA PREPARATION (perform transform on all data) ----------------------

## F.1) Checking transforms ------------------------------------------------
# if we wish a normal response (ICUs) we may want to test which transform is best:
# original response
veneto_alldata %>% dplyr::summarise(dplyr::across(ICU, mean)) #[1] 109
veneto_alldata %>% dplyr::summarise(dplyr::across(ICU, sd)) #[1] 113
veneto_alldata %>% dplyr::select(ICU) %>% PerformanceAnalytics::skewness() #[1] 0.9264596 positive skew
veneto_alldata %>% dplyr::select(ICU) %>% PerformanceAnalytics::kurtosis() #[1] -0.4800243 (excess) platykurtic

# try natural log
veneto_alldata %>% dplyr::mutate(transf = log(ICU)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::skewness() #[1] Nan
veneto_alldata %>% dplyr::mutate(transf = log(ICU)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::kurtosis() #[1] Nan
# try log10
veneto_alldata %>% dplyr::mutate(transf = log10(ICU)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::skewness() #[1] Nan
veneto_alldata %>% dplyr::mutate(transf = log10(ICU)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::kurtosis() #[1] Nan
# try log(1+x)
veneto_alldata %>% dplyr::mutate(transf = log1p(ICU)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::skewness() #[1] -0.5764408
veneto_alldata %>% dplyr::mutate(transf = log1p(ICU)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::kurtosis() #[1] -0.5802109
# try square root
veneto_alldata %>% dplyr::mutate(transf = (ICU)^(1/2)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::skewness() #[1] 0.339607
veneto_alldata %>% dplyr::mutate(transf = (ICU)^(1/2)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::kurtosis() #[1] -1.142681
# try cubic root
veneto_alldata %>% dplyr::mutate(transf = (ICU)^(1/3)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::skewness() #[1] 0.009550655
veneto_alldata %>% dplyr::mutate(transf = (ICU)^(1/3)) %>% 
  dplyr::select(transf) %>% PerformanceAnalytics::kurtosis() #[1] -1.073347


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

# other useful columns-wise transformations to predictors

## F.3) Numerical variables transformations ------------------------------------------
veneto_alldata <- veneto_alldata %>% 
  dplyr::mutate(vaccinated_total = tidyr::replace_na(vaccinated_total, 0),
                delta_ICU = ICU - dplyr::lag(ICU, n = 1L, default = 0),
                days = dplyr::row_number())


## F.4) Set up a train/test split -----------------------------------------

set.seed(25) # fix seed to make reproducible results

# on required time series
veneto_alldata <- veneto_alldata %>%
  dplyr::filter(dplyr::between(date, lubridate::as_date("2020-10-01"), lubridate::as_date("2021-02-15"))) %>% 
  dplyr::mutate(days = row_number())

split <- veneto_alldata %>% rsample::initial_time_split(prop = 0.90, lag = 0) #TODO: change % here
train_data <- split %>% rsample::training() #%>% tibble::rownames_to_column(var = "date")
test_data <- split %>% rsample::testing() #%>% tibble::rownames_to_column(var = "date")

train_data %>% dim() # check if is 80% of data points [1] 0.8985507
test_data %>% dim() # checks if is 20% of data points [1] 0.1014493

## F.5) Set up k-fold(s) cross validation ----------------------------------

train_cv <- train_data %>% rsample::loo_cv() # leave one out cv
#train_cv <- train_data %>% rsample::vfold_cv(v = 10, repeats = 2) # NOT FOR TS (k-folds + repetitions) cv
train_cv <- train_data %>% rsample::sliding_index(lookback = 0L, # FOR TS ONLY
                                                  assess_start = 1L,
                                                  assess_stop = 1L,
                                                  complete = TRUE,
                                                  step = 1L,
                                                  skip = 0L)

# G) FEATURE ENGINEERING (separate transforms on train/test data) ----------
# using recipe library see https://recipes.tidymodels.org/articles/Ordering.html

# preditors_to_exclude <- c("date", "ICU", "residents", "deaths", "hospitalized_total", 
#                            "self_isolation", "positives_total", "total_cases", "swabs", "vaccinated_total")
preditors_to_exclude <- c("hospitalized_total", "self_isolation")

preditors_to_include <- c("hospitalized_total", "self_isolation", "positives_total",
                          "total_cases", "swabs", "vaccinated_total")
## G.1) Recipes LM -------------------------------------------------------
recipe_lm <- recipes::recipe(delta_ICU ~ days+hospitalized_total+self_isolation+swabs, data = split) %>%
             #recipes::step_mutate(days = row_number()) %>%   
             #recipes::step_select(ICU, days) %>%  #ERROR table predicts!!!
             #recipes::update_role(all_of(preditors_to_exclude), new_role = "excluded_predictors") %>%  
             recipes::step_lag(swabs, hospitalized_total, self_isolation, lag = 1) %>% 
             recipes::step_ns(days, deg_free = 3) %>% 
             #recipes::step_rm("deaths") %>%
             #recipes::step_date(date) %>%
             #recipes::step_poly(days, degree = 3)
             #recipes::step_corr(all_numeric(), threshold = 0.80) %>%
             #recipes::step_center(all_predictors()) %>%
             #recipes::step_scale(all_predictors()) %>%
             #recipes::step_ns(days, deg_free = tune::tune("days df")) %>%
             #recipes::step_dummy(all_nominal_predictors())
             recipes::step_interact(~ all_predictors():all_predictors())           
             #recipes::step_normalize(all_numeric()) %>%
             #recipes::step_log(swabs, base = 10) %>%
             #recipes::step_naomit(all_predictors())

summary(recipe_lm)
# apply recipe to train data & test data
veneto_train_preprocessed_lm <- recipe_lm %>%
  recipes::prep(train_data) %>% 
  recipes::juice() %>% view()

veneto_test_preprocessed_lm <- recipe_lm %>%
  recipes::prep(test_data) %>% 
  recipes::juice() %>% view()

## G.2) Recipes GLM -------------------------------------------------------
recipe_glm_pois <- recipes::recipe(ICU ~ days + residents, case_weight = residents, data = split) #%>%
                   #recipes::step_mutate(days = row_number()) %>%
                   #recipes::step_ns(days, deg_free = 3)

summary(recipe_glm_pois)
# apply recipe to train data & test data
veneto_train_preprocessed_glm_pois <- recipe_glm_pois %>%
  recipes::prep(train_data) %>% 
  recipes::juice() %>% view()

veneto_test_preprocessed_glm_pois <- recipe_glm_pois %>%
  recipes::prep(test_data) %>% 
  recipes::juice() %>% view()
                
## G.3) Recipes ARIMA -------------------------------------------------------
recipe_arima <- recipes::recipe(ICU ~ date, data = split) %>%
                recipes::step_select(date, ICU)

summary(recipe_arima)
# apply recipe to train data & test data
veneto_train_preprocessed_arima <- recipe_arima %>%
  recipes::prep(train_data) %>% 
  recipes::juice() %>% view()

# apply recipe to train data & test data
veneto_test_preprocessed_arima <- recipe_arima %>%
  recipes::prep(test_data) %>% 
  recipes::juice() %>% view()

# ts plot
veneto_alldata %>% timetk::plot_time_series(date, ICU, .interactive = FALSE)
veneto_train_preprocessed_arima %>% timetk::plot_time_series(date, ICU, .interactive = FALSE)
veneto_test_preprocessed_arima %>% timetk::plot_time_series(date, ICU, .interactive = FALSE)

# format transform tibble to ts 
veneto_train_preprocessed_arima <- veneto_train_preprocessed_arima %>% 
  tsibble::as_tsibble() %>% stats::as.ts()
veneto_test_preprocessed_arima <- veneto_test_preprocessed_arima %>% 
  tsibble::as_tsibble() %>% stats::as.ts()


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
# recap which are the tuning params in recipe
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
fit1 <- stats::lm(delta_ICU ~ stats::poly(days, degree = 3), data = veneto_train_preprocessed_lm)
fit1 <- stats::lm(delta_ICU ~ splines::bs(days, df = 3), data=veneto_train_preprocessed_lm)
fit1 <- stats::lm(delta_ICU ~ lag_1_swabs+lag_1_hospitalized_total+lag_1_self_isolation, data=veneto_train_preprocessed_lm)
summary(fit1) # check with tidy
stats::extractAIC(fit1) # starting AIC + check tidy
stepAIC(fit1, list(upper = ~ lag_1_swabs+lag_1_hospitalized_total+lag_1_self_isolation, lower = ~ 1), direction = "backward")


## H.2) Model 2: GLM Poisson with intercept + cubic spline on time ---------
glm_pois_engine <- poissonreg::poisson_reg() %>% 
                   parsnip::set_engine("glm", family=stats::poisson(link="log")) %>% #, offset = log(residents)
                   parsnip::set_mode("regression")

glm_pois_workflow <- workflows::workflow() %>%
                     #workflows::remove_variables() %>% 
                     workflows::add_recipe(recipe_glm_pois) %>%
                     workflows::add_model(glm_pois_engine)

# fit just the best model 
glm_pois_train_fit <- fit(glm_pois_workflow, train_data)

# store results for evaluation
glm_pois_results <- glm_pois_train_fit %>% hardhat::extract_fit_engine()

# quick peek at results
glm_pois_results %>% summary()
generics::glance(glm_pois_train_fit)

# check if results are in line with tidymodels: OK!
fit2.1 <- stats::glm(ICU ~ splines::ns(days, 3), data = veneto_train_preprocessed_glm_pois,
                     family = stats::poisson(link = "log"))#, offset = log(residents))
summary(fit2.1)
fit2.2 <- stats::glm(ICU/residents ~ splines::ns(days, 3), data = veneto_train_preprocessed_glm_pois,
                     family = stats::poisson(link = "log"), weights = residents)
summary(fit2.2)
anova2 <- anova(fit2.1, test="Chisq") # check which predictors to keep
# residual deviance[168.17] > df[120] signals over-dispersion!

## H.3) Model 3: GLM Quasi-Poisson with intercept + cubic spline on time ----



# check if results are in line with tidymodels: OK!
fit3 <- stats::glm(ICU ~ ns(days, 3), data = veneto_train_preprocessed_glm_pois,
                   family = stats::quasipoisson(link = "log"), offset = log(residents))
summary(fit3)
an3 <- anova(fit3, test="F")


## H.4) Model 4: auto-ARIMA (1,2,2) -----------------------------------------
arima_engine <- modeltime::arima_reg() %>% 
  parsnip::set_engine("auto_arima") %>%
  parsnip::set_mode("regression")

arima_workflow <- workflows::workflow() %>%
  #workflows::remove_variables() %>% 
  workflows::add_recipe(recipe_arima) %>%
  workflows::add_model(arima_engine)

# fit just the best model 
arima_train_fit <- fit(arima_workflow, train_data)

# store results for evaluation
arima_results <- arima_train_fit %>% hardhat::extract_fit_engine()

# quick peek at results
arima_results

# check if results are in line with forecast: OK!
fit4 <- forecast::auto.arima(veneto_train_preprocessed_arima, trace = TRUE)
summary(fit4)


## H.5) Model 5: ARIMA with best lag and cross correlation best lags -------
# discover best lag for Y variable
veneto_train_preprocessed_arima %>% diff() %>% stats::decompose() %>% plot()
veneto_train_preprocessed_arima %>% diff() %>% stats::acf() # auto correlation to get q
veneto_train_preprocessed_arima %>% diff() %>% stats::pacf() # partial auto correlation to get p
veneto_train_preprocessed_arima %>% diff() %>% stats::Box.test(lag = 10, type = "Ljung-Box") # daily change is random, uncorrelated with previous days
#veneto_train_preprocessed_arima %>% diff() %>% tseries::adf.test(alternative="stationary", k=0)

best_lag <-  VARselect(diff(veneto_train_preprocessed_arima), lag.max = 10)
n <- best_lag$selection[1] # best AR lag is 8 
#auto-regressive distributed lag (ADL) 

arima_engine2 <- modeltime::arima_reg(seasonal_period = "auto", # periodic nature of seasonality, default = "auto" 
                                      non_seasonal_ar = 0, # p-order of NSAR terms
                                      non_seasonal_differences = 2, # d-order of NS differencing 
                                      non_seasonal_ma = 1, # q-order of NSMA terms
                                      seasonal_ar = 0, # P-order of SAR terms
                                      seasonal_differences = 0, # D-order of S differencing 
                                      seasonal_ma = 1) %>% # Q-order of SMA terms
                 parsnip::set_engine("arima") %>%
                 parsnip::set_mode("regression")

arima_workflow2 <- workflows::workflow() %>%
  #workflows::remove_variables() %>% 
  workflows::add_recipe(recipe_arima) %>%
  workflows::add_model(arima_engine2)

# fit just the best model (default ARIMA(0,0,0)
arima_train_fit2 <- fit(arima_workflow2, train_data)

# store results for evaluation
arima_results2 <- arima_train_fit2 %>% 
  hardhat::extract_fit_engine()

# quick peek at results
arima_results2


# I) DATA PREDICTION & EVALUATION ---------------------------------------
# here we valuate predictions for each model to choose the best

## I.1) LM with time spline ---------------------------------------------

# a) predict on train set
lm_train_res <- lm_train_fit %>% prediction_table(train_data, is_delta = TRUE)

# b) plot train results
plot_fitted(lm_train_res)
lm_train_res %>% dplyr::mutate(date = row_number()) %>% plot_series()
lm_results %>% vip::vip()
lm_results %>% plot_residuals()

# c) get metrics on train set
lm_train_res %>% pull_metrics() %>% print()

# d) predict on validation set
lm_validation_res <- lm_workflow %>%  
                     tune::fit_resamples(resamples = train_cv,
                     #control = tune::control_resamples(extract = get_model),
                     metrics = yardstick::metric_set(rmse, rsq, mae))

# e) get metrics on validation set
lm_validation_res %>% tune::collect_metrics(summarize = TRUE)

# f) predict on test set
lm_test_res <- lm_train_fit %>% prediction_table(test_data, is_delta = TRUE)

# g) plot test results (best practice is to keep transformed scale)
lm_test_res %>% plot_fitted()
lm_test_res %>% dplyr::mutate(date = row_number()) %>% plot_series()
fit1 %>% forecast::forecast(newdata = veneto_test_preprocessed_lm %>%
                                      mutate(across(everything(), .fns = ~replace_na(.,0))), # na_omit
                            level = c(95)) %>% plot() #see table with confidence bands

# h) get metrics on test set   
lm_test_res %>% pull_metrics() %>% print()

# h.bis) alternative test error check via workflow
lm_workflow %>% 
  tune::last_fit(split = split) %>% 
  tune::collect_metrics()


## I.2) GLM Poisson --------------------------------------------------------

# a) predict on train set
glm_pois_train_res1 <- glm_pois_train_fit %>% prediction_table(train_data, is_delta = FALSE)
glm_pois_train_res <- fit2.1 %>% prediction_table(train_data, is_delta = FALSE)
glm_pois_train_fit %>% stats::predict(train_data)# type = "response"
fit2.1 %>% stats::predict(train_data,  type = "response")

# b) plot train results
glm_pois_train_res %>% plot_fitted()
glm_pois_train_res %>% dplyr::mutate(date = row_number()) %>% plot_series()
glm_pois_results %>% vip::vip()
glm_pois_results %>% plot_residuals()

# c) get metrics on train set
glm_pois_train_res %>% pull_metrics() %>% print()

# d) predict on validation set
glm_pois_validation_res <- glm_pois_workflow %>%  
                           tune::fit_resamples(resamples = train_cv,
                           #control = tune::control_resamples(extract = get_model),
                           metrics = yardstick::metric_set(rmse, rsq, mae))

# e) get metrics on validation set
glm_pois_validation_res %>% tune::collect_metrics(summarize = TRUE)

# f) predict on test set
glm_pois_test_res1 <- glm_pois_train_fit %>% prediction_table(test_data, is_delta = FALSE)
glm_pois_test_res <- fit2.1 %>% prediction_table(test_data, is_delta = FALSE)

# g) plot test results (best practice is to keep transformed scale)
glm_pois_test_res %>% plot_fitted()
glm_pois_test_res %>% dplyr::mutate(date = row_number()) %>% plot_series()
fit2.1 %>% forecast::forecast(h = 14, level = c(95)) %>% plot() #see table with confidence bands

# h) get metrics on test set   
glm_pois_test_res %>% pull_metrics() %>% print()

# h.bis) alternative test error check via workflow
glm_pois_workflow %>% 
  tune::last_fit(split = split) %>% 
  tune::collect_metrics()


## I.3) ARIMA ----------------------------------------------------

# a) predict on train set
arima_train_res <- arima_train_fit %>% calibration_table(train_data)

# b) plot train results
arima_train_res %>% plot_fitted()
arima_train_res %>% plot_series()
fit4 %>% forecast::checkresiduals() # also (arima_results$data$.residuals)

# c) get metrics on train set
arima_train_res %>% pull_metrics() %>% print()

# d) predict on validation set
arima_validation_res <- arima_workflow %>%  
                        tune::fit_resamples(resamples = train_cv,
                        #control = tune::control_resamples(extract = get_model),
                        metrics = yardstick::metric_set(rmse, rsq, mae))

# e) get metrics on validation set
arima_validation_res %>% tune::collect_metrics(summarize = TRUE)

# f) predict on test set
arima_test_res <- arima_train_fit %>% calibration_table(test_data)

# g) plot test results (best practice is to keep transformed scale)
arima_test_res %>% plot_fitted()
arima_test_res %>% plot_series()
fit4 %>% forecast::forecast(h = 14, level = c(95)) %>% plot() #see table with confidence bands

# h) get metrics on test set   
arima_test_res %>% pull_metrics() %>% print()

# h.bis) alternative test error check via workflow
arima_workflow %>% 
  tune::last_fit(split = split) %>% 
  tune::collect_metrics()


# J) COMPARISONS & CONCLUSIONS -------------------------------------------
rbind(base_train_rmse, base_test_rmse,
      tree_tidy_train_rmse, tree_tidy_test_rmse,
      tree_caret_train_rmse, tree_caret_test_rmse)
#lsf.str("package:dplyr") #to list all functions in package
#environmentName(environment(select)) #to get package name from function
