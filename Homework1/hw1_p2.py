import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DJIA = pd.read_csv('DJIA-Data.csv') # read in the file into a DataFrame

noa = len(DJIA.columns) # number of assets

mu = DJIA.mean() # sample mean vector for each asset

cov_mat = DJIA.cov() # covariance matrix for DataFrame

std_dev = DJIA.std() # standard deviation for each asset

'''
This function computes the sample mean of the portfolio
'''
def sampleMean(weights, mu):
	return np.sum(mu * weights)

'''
This function computes the volatility of the portfolio
'''
def sampleVol(weights, cov):
	return np.sqrt(np.dot(np.dot(cov, weights),
		weights.T))



u = np.array([1] * len(mu)) # initialize a unit vector

'''
The following is an implementation to find components for
the weights of the minimum variance portfolio
found in our textbook, Statistical Models and Methods 
for Financial Markets
'''
A = np.dot(np.dot(mu.T, np.linalg.inv(cov_mat)), u)
B = np.dot(np.dot(mu.T, np.linalg.inv(cov_mat)), mu)
C = np.dot(np.dot(u.T, np.linalg.inv(cov_mat)), u)
D = B*C - A*A
p1 = np.dot(B * np.linalg.inv(cov_mat), u)
p2 = np.dot(A * np.linalg.inv(cov_mat), mu)
p3 = np.dot(C*np.linalg.inv(cov_mat), mu) - np.dot(A*np.linalg.inv(cov_mat),u)
    
eff_mu = [] # initalize lists to store the efficient means
eff_vol = [] # initialize lists to store the efficient volalities

'''
Calculates the weights of the minimum variance portfolio 
where (A/C) is chosen to minimize the variance
'''
mvp_weights = ((p1 - p2 + (A/C) *(p3)) / D)
    
mvp_mu = np.dot(mu, mvp_weights.T) # calculates mean of MVP
mvp_vol = np.sqrt((B - 2 * mvp_mu*A + (mvp_mu**2)*C) /D) # calculates volatility of MVP
print(mvp_vol)

'''
Next, we want to loop through possible target values for the mean
from which we can gather the efficient portfolios.
'''
for x in np.arange(mvp_mu, 12, 0.1):
    eff_mu.append(sampleMean(((p1 - p2 + x *(p3)) / D), mu))
    eff_vol.append(sampleVol(((p1 - p2 + x *(p3)) / D), cov_mat))
    
eff_mu = np.array(eff_mu)
eff_vol = np.array(eff_vol)



specific_weights = ((p1 - p2 + 0.3 *(p3)) / D) # risk-free rate of 0.3

'''
Question 2(a)
'''
print(mvp_weights) # weights of the MVP



'''
Question 2(b)
'''

market_weights = np.dot(np.linalg.inv(cov_mat), mu - 0.15*u)/(np.dot(np.dot(u.T, np.linalg.inv(cov_mat)), mu - 0.15*u))
market_mu = np.dot(mu, market_weights.T)
market_volatility = np.sqrt((B - 2 * market_mu*A + (market_mu**2)*C) /D)

print(market_volatility)

'''
Question 2(c)
'''
print(specific_weights) # weights for expected return of 0.3%


capm_mu = []
capm_vol = [x for x in np.arange(0, 20, 0.1)]
for x in capm_vol:
    capm_mu.append(0.15 + (market_mu - 0.15)/market_volatility*x)


'''
This function plots the efficient frontier
along with the 30 points and the minimum variance portfolio
'''

plt.figure(figsize=(8,6))
plt.title('Efficient Frontier')
plt.xlabel('Standard Deviation of the Portfolio')
plt.ylabel('Return of the Portfolio')
plt.scatter(eff_vol, eff_mu, color='black', marker='x')
plt.plot(mvp_vol, mvp_mu, color='green', marker='o')
plt.scatter(std_dev, mu)
plt.plot(market_volatility, market_mu, color='red', marker='^')
plt.plot(0, 0.15, color='blue', marker='o')
plt.plot(capm_vol, capm_mu, color='gray', marker='*')
plt.show()