import numpy as np 
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

'''
Unless you want to pip install all the libraries/modules
you should use Jupyter for this code
'''

# You probably have a different path
allData = pd.read_csv('./AllData.csv')

# Risk Free
rf = 0.15

# Risk free subtracted from stocks
monthly_rf = allData.drop('S&P500', axis=1).copy() - rf

# Linear Regression

# Risk free subtracted from 
# S&P500 returns used as market returns 
market_rf = allData['S&P500'] - rf

# smf.OLS doesn't provide a constant automatically
# so we have to add it using this function
market_rf = sm.add_constant(market_rf)

# creates a dictionary of fits where key is the stock name
# and value is the fit
results = {x : smf.OLS(monthly_rf[x], market_rf).fit() for x in monthly_rf}

# look at all the information
# run through the list of stock names
for x in monthly_rf:
	print(results[x].summary())

# Prints the Alpha and Beta estimates for each stock
for x in monthly_rf:
	print('Stock: %s   Alpha: %8.4f   Beta: %8.4f' % (x, results[x].params[0], results[x].params[1]))

# creates a list of underperforming stocks
under_performers = [key for key in results.keys() if results[key].tvalues[0] < -1.697]

# Prints out the underperformers
for x in under_performers:
	print('%s' % x)

# creates a list of overperforming stocks
over_performers = [key for key in results.keys() if results[key].tvalues[0] > 1.697]

# Prints out the overperformers
for x in over_performers:
	print('%s' % x)

# Sample market return
market_return = allData['S&P500'].mean()

# Sample return from 30 DJI stocks
monthly_return = allData.drop('S&P500', axis=1).copy().mean()

# Sample covariance from 30 DJI stocks
monthly_covariance = allData.drop('S&P500', axis=1).copy().cov()

# calculates returns through CAPM
returns = [(rf + results[x].params[1]*(market_return - rf)) for x in results]

# Array of 1s
u = np.array([1]*len(monthly_return))

# Calculates market weights with CAPM returns
market_weights = np.dot(np.linalg.inv(monthly_covariance),
	returns - rf*u)/(np.dot(np.dot(u.T, np.linalg.inv(monthly_covariance)),
		returns - rf*u))

# Calculates market weights with sample returns
og_market_weights = np.dot(np.linalg.inv(monthly_covariance),
	monthly_return - rf*u)/(np.dot(np.dot(u.T, np.linalg.inv(monthly_covariance)),
		monthly_return - rf*u))
