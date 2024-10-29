rm(list=(ls()))

# Setting up the workflow
## Loading the necessary packages:
library(tidyverse)
library(tidyr)
library(forecast)
library(vars)
library(tseries)
library(AER)
library(urca)
library(lmtest)
library(fpp3)
library(MASS)

## Importing the dataset:
data <- read.csv("Global Dataset of Inflation.csv")

## Selecting the analyzed countries:
selected_countries <- data %>% filter(Country %in% c("Hungary", "Italy", "Spain"))


## Transforming the data into longitudinal
data_long <- selected_countries %>%
  gather(key = "year", value = "inflation", X1970:X2022) %>%
  mutate(year = as.numeric(gsub("X", "", year)))

## EDA
ggplot(data_long, aes(x=year, y=inflation, color=Country)) +
  facet_wrap(~Series.Name, scales = "free_y") +
  geom_line() +
  labs(title="Inflation Trends: 1970-2022", y="Inflation Rate (%)", x="Year") +
  theme_minimal()


#  ARIMA

## The following function runs the ARIMA model for a given country and a given inflation type:
ARIMA_function <- function(first_year, ARIMA_country, ARIMA_inflation){

  print(ARIMA_country)
  print(ARIMA_inflation)
  
  country_inflation_values <- data_long %>% 
    filter(Country == ARIMA_country & Series.Name == ARIMA_inflation & year >= first_year) %>%
    arrange(year) %>%
    pull(inflation)
  
  adf_test_ts <- ts(country_inflation_values, start=first_year, frequency = 1)
  print("Seeing if the values are stationary:")
  adf_test <- adf.test(adf_test_ts, alternative = "stationary")
  acf(adf_test_ts, main=paste0("ACF - ",ARIMA_country," , ", ARIMA_inflation))
  pacf(adf_test_ts, main=paste0("PACF - ",ARIMA_country," , ", ARIMA_inflation))
  
  if (adf_test$p.value < 0.05) {
    d <- 0
  } else {
    p_value_diff1 <- adf.test(diff(adf_test_ts, differences = 1))$p.value
    
    if (p_value_diff1 < 0.05) {
      d <- 1
    } else {
      p_value_diff2 <- adf.test(diff(adf_test_ts, differences = 2))$p.value
      
      if (p_value_diff2 < 0.05) {
        d <- 2
      } else {
        d <- 3
      }
    }
  }
  
  
  ts_as_tsibble <- ts(country_inflation_values, start=first_year, frequency = 1) %>% as_tsibble
  arima_model <- ts_as_tsibble %>%
    model(ARIMA(value ~ pdq(c(NA, d, NA)), ic = "aic", stepwise = FALSE, greedy = FALSE))
  report(arima_model)
  
  
  forecast_arima <- arima_model %>% forecast(h=8)
  last_20_years <- tail(ts_as_tsibble, 20)
  
  
  autoplot(forecast_arima) +
    labs(
      title = paste0("Forecast of ", ARIMA_country,"'s ", ARIMA_inflation, " using ARIMA"),
      subtitle = "With 95% Confidence Interval",
      x = "Year",
      y = "Inflation"
    ) +
    theme_minimal() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) + autolayer(last_20_years)

  # Checking how accurate the forecast is:
  train_length <- floor(0.8 * length(country_inflation_values))
  train_data <- head(country_inflation_values, train_length)
  test_data <- tail(country_inflation_values, -train_length)
  
  training_ts <- ts(train_data, start=first_year, frequency = 1)
  arima_fit <- forecast::auto.arima(training_ts)
  
  forecast_obj <- forecast::forecast(arima_fit, h=length(test_data))
  
  accuracy_vals <- forecast::accuracy(forecast_obj, test_data)
  print(accuracy_vals)
  
  
}



## Running the function for all combinations of countries and inflation types:
ARIMA_function(1970,"Hungary","Headline Consumer Price Inflation")
ARIMA_function(1970,"Italy","Headline Consumer Price Inflation")
ARIMA_function(1970,"Spain","Headline Consumer Price Inflation")

ARIMA_function(1970,"Hungary","Energy Consumer Price Inflation")
ARIMA_function(1970,"Italy","Energy Consumer Price Inflation")
ARIMA_function(1970,"Spain","Energy Consumer Price Inflation")

ARIMA_function(1970,"Hungary","Food Consumer Price Inflation")
ARIMA_function(1970,"Italy","Food Consumer Price Inflation")
ARIMA_function(1970,"Spain","Food Consumer Price Inflation")

ARIMA_function(1991,"Hungary","Official Core Consumer Price Inflation")
ARIMA_function(1970,"Italy","Official Core Consumer Price Inflation")
ARIMA_function(1976,"Spain","Official Core Consumer Price Inflation")

ARIMA_function(1979,"Hungary","Producer Price Inflation")
ARIMA_function(1975,"Italy","Producer Price Inflation")
ARIMA_function(1973,"Spain","Producer Price Inflation")




# VAR


## VAR function:
VAR_function <- function(VAR_inflation){
  var_data_long <- data_long %>% 
    filter(Series.Name == VAR_inflation) %>%
    spread(key = Country, value = inflation)
  
  var_data <- var_data_long %>%
    group_by(year) %>%
    summarise(Hungary = sum(Hungary, na.rm = TRUE),
              Italy = sum(Italy, na.rm = TRUE),
              Spain = sum(Spain, na.rm = TRUE))
  
  var_data_tsi <- var_data %>% 
    as_tsibble(index = year)
  
  
  var_diff <- var_data_tsi %>% 
    mutate(Hungary = difference(Hungary), 
           Italy = difference(Italy), 
           Spain = difference(Spain)) %>%
    drop_na()
  
  models <- var_diff %>%
    model(var_aic = VAR(vars(Hungary, Italy, Spain), ic = "aic"),
          var_bic = VAR(vars(Hungary, Italy, Spain), ic = "bic"))
  
  ## The results of both models
  aic_model <- models$var_aic[[1]]
  report(aic_model)
  bic_model <- models$var_bic[[1]]
  report(bic_model)
  glance(models)
  
  aic_residuals <- residuals(models$var_aic[[1]])
  aic_residuals %>%
    features(Hungary, features = ljung_box, lag = 20)
  
  
  print(grangertest(Hungary ~ Italy, order = 2, data = var_diff))
  print(grangertest(Hungary ~ Spain, order = 2, data = var_diff))
  print(grangertest(Italy ~ Hungary, order = 2, data = var_diff))
  print(grangertest(Italy ~ Spain, order = 2, data = var_diff))
  print(grangertest(Spain ~ Hungary, order = 2, data = var_diff))
  print(grangertest(Spain ~ Italy, order = 2, data = var_diff))
  
  model.bic <- vars::VAR(as.data.frame(var_diff), p = 1)
  plot(irf(model.bic, impulse = "Hungary", response = "Italy"))
  plot(irf(model.bic, impulse = "Italy", response = "Spain"))
  plot(irf(model.bic, impulse = "Spain", response = "Hungary"))
  
  plot(fevd(model.bic, n.ahead = 8))
  
  
  ## Evaluating the model:
  train_length <- floor(0.8 * nrow(var_diff))
  train_data <- head(var_diff, train_length)
  train_data_ts <- ts(train_data, frequency = 1) 
  test_data <- tail(var_diff, -train_length)
  
  var_fit <- vars::VAR(train_data_ts, ic = "aic")
  
  forecast::accuracy(var_fit$varresult[[2]])
}
## Running the function for all inflation types:
VAR_function("Headline Consumer Price Inflation")
VAR_function("Energy Consumer Price Inflation")
VAR_function("Food Consumer Price Inflation")
VAR_function("Official Core Consumer Price Inflation")
VAR_function("Producer Price Inflation")


