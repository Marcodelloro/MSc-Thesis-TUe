import numpy as np
import warnings
import itertools
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.tsa.holtwinters import ExponentialSmoothing
from prophet import Prophet

from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from pmdarima import auto_arima
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller

params = pd.read_csv('table_parAUS.csv')
params_to_forecast = ['alpha', 'gamma', 'delta', 'sigma', 'tau']
df = params[['alpha', 'gamma', 'delta', 'sigma', 'tau']]
weeksIn = [0, 7, 14, 21, 28, 35]
days_ahead = 21   # Estimation horizon (3 weeks)
# horizon = 4 * days_ahead # --> 'ForecastsAUS_CI.xlsx'
# horizon = 8 * days_ahead # --> 'ForecastsAUS_CI2.xlsx'
horizon = 10 * days_ahead # --> 'ForecastsAUS_CI3.xlsx'

with pd.ExcelWriter('ForecastsAUS_CI.xlsx', engine='xlsxwriter') as writer:
    for k in range(len(weeksIn)):
        df_sub = df.iloc[0:horizon + weeksIn[k]]
        forecasts = pd.DataFrame(index=np.arange(days_ahead))  # Creating an empty DataFrame with the right index

        date_range_prophet = pd.date_range(start=df_sub.index[-35], periods=35, freq='D')
        combined_forecast = pd.DataFrame()
        last_values = {}

        for param in params_to_forecast:
            train_holt = df_sub[param][-35:]  # Training data for the current parameter
            test_holt = df[param][horizon:horizon + weeksIn[k] + days_ahead]
            last_value_holt = train_holt.iloc[-1]
            hwmodel_damped = ExponentialSmoothing(train_holt, trend='mul', damped_trend=True, seasonal=None).fit(smoothing_trend=0.5)
            test_pred = hwmodel_damped.forecast(days_ahead)
            num_simulations = 1000
            simulated_paths = hwmodel_damped.simulate(nsimulations=days_ahead, repetitions=num_simulations, error='mul')
            lower_quantiles = np.percentile(simulated_paths, 10, axis=1)
            upper_quantiles = np.percentile(simulated_paths, 90, axis=1)

            test_pred_reset_index = test_pred.reset_index(drop=True)

            forecasts[f'{param}_yhat'] = test_pred_reset_index
            forecasts[f'{param}_yhat_lower'] = lower_quantiles
            forecasts[f'{param}_yhat_upper'] = upper_quantiles

    # Plotting the results of the simulation with Holm

            plt.figure(figsize=(10, 6))
            plt.plot(train_holt.index, train_holt, label='Training Data')
            plt.plot(test_pred.index, test_pred, label='Forecast', color='red')
            plt.plot(test_holt.index, test_holt, label='Test Data', color='gray')
            plt.fill_between(test_pred.index, lower_quantiles, upper_quantiles, color='red', alpha=0.3, label='95% Confidence Interval')
            plt.xlabel('Time')
            plt.ylabel(param)
            plt.title(f'Holt-Winters Forecast with Confidence Interval for {param}')
            plt.legend()
            plt.show()

    # Simulation with Prophet
            prhData= df_sub[param].tail(35)
            fitProphet = pd.DataFrame({
                'ds': date_range_prophet,
                'y': prhData.values
            })

            # Fit the Prophet model
            m = Prophet()
            m.fit(fitProphet)
            future = m.make_future_dataframe(periods=days_ahead)
            forecast = m.predict(future)

            # Extract and rename the last 21 values from the Prophet forecast
            last_21_forecast = forecast[['yhat', 'yhat_lower', 'yhat_upper']].tail(21).reset_index(drop=True)
            last_21_forecast.columns = [f'{param}_Prophet_yhat', f'{param}_Prophet_yhat_lower',
                                        f'{param}_Prophet_yhat_upper']
            forecasts = pd.concat([forecasts, last_21_forecast], axis=1)

            last_value = df_sub[param].iloc[-1]  # Get the last value of the training set for the current parameter
            last_values[f'{param}_yhat'] = [last_value]
            last_values[f'{param}_yhat_lower'] = [last_value]
            last_values[f'{param}_yhat_upper'] = [last_value]
            last_values[f'{param}_Prophet_yhat'] = [last_value]
            last_values[f'{param}_Prophet_yhat_lower'] = [last_value]
            last_values[f'{param}_Prophet_yhat_upper'] = [last_value]

        last_values_df = pd.DataFrame(last_values, index=[-1])  # Create a DataFrame with the last values

        forecasts = pd.concat([last_values_df, forecasts]).reset_index(drop=True)

        forecasts.to_excel(writer, sheet_name=f'Horizon_{weeksIn[k]}_days', index=False)