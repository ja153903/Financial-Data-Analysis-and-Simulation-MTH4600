import numpy as np 
import matplotlib.pyplot as plt 

S0 = 100.0 # stock price
K = 110.0 # strike price
T = 0.25 # time at maturity
r = 0.03 # short rate
M = 10000 # number of simulations
I = 252 # number of stock periods
dt = 1.0 / I # number of steps

discount_factor = np.exp(-r * T) # discount factor for option evaluation

op1_price = {} # dictionary to hold the values of the rho-options price

'''
Initialize these three volatilities
'''
s1_vol = 0.2
s2_vol = 0.25
s3_vol = 0.3

'''
Calculate drift terms once and for all
'''
s1_drift = (r - s1_vol*s1_vol*0.5) * dt 
s2_drift = (r - s2_vol*s2_vol*0.5) * dt 
s3_drift = (r - s3_vol*s3_vol*0.5) * dt 

'''
Calculate np.sqrt()*vol once and for all
'''
s1_rand = np.sqrt(dt) * s1_vol
s2_rand = np.sqrt(dt) * s2_vol
s3_rand = np.sqrt(dt) * s3_vol

rhos = [x for x in np.arange(0.0, 1.1, 0.1)] # generates the rhos

op1_price = {}
op2_price = {}

vbar = v2bar = avgbar = avg2bar = 0
print('rho', 'vbar', 'errorv', 'avgbar', 'errora')
for rho in rhos:
	s1 = s2 = s3 = np.zeros(I+1)
	covariance_matrix = [[1, rho, rho],[rho, 1, rho], [rho, rho, 1]]
	op1_payoff_sum = 0

	if rho != 1.0:
		cholesky_matrix = np.linalg.cholesky(covariance_matrix) # cholesky decomp
	else:
		cholesky_matrix = [[1, 0, 0], [1, 0, 0], [1, 0, 0]]


	for i in range(1, M+1):

		uncorrelated_matrix = np.random.standard_normal((3, 63)) # generate uncorrelated normal
		correlated_matrix = np.dot(cholesky_matrix, uncorrelated_matrix) # generate correlated normal

		s1cs = np.cumsum(correlated_matrix[0])
		s2cs = np.cumsum(correlated_matrix[1])
		s3cs = np.cumsum(correlated_matrix[2])

		s1 = S0 * np.exp(s1_drift + s1_rand * s1cs)
		s2 = S0 * np.exp(s2_drift + s2_rand * s2cs)
		s3 = S0 * np.exp(s3_drift + s3_rand * s3cs)

		smax = max(np.max(s1), np.max(s2), np.max(s3))
		savgmax = max((s1[len(s1)-1] + s2[len(s2)-1] + s3[len(s3)-1])/3, 0)

		payoff = max(smax - K, 0) * discount_factor
		avgpayoff = max(savgmax - K, 0) * discount_factor
		
		vbar = ((i-1) * vbar + payoff)/i
		avgbar = ((i-1) * avgbar + avgpayoff)/i

		v2bar = ((i-1) * v2bar + payoff*payoff)/i
		avg2bar = ((i-1) * avg2bar + avgpayoff*avgpayoff)/i

		if i % 10000 == 0:

			std_dev_v = np.sqrt(v2bar - vbar * vbar)
			std_dev_a = np.sqrt(avg2bar - avgbar * avgbar)

			error_v = 1.96 * std_dev_v / np.sqrt(i)
			error_a = 1.96 * std_dev_a / np.sqrt(i)

			print(rho, vbar, error_v, avgbar, error_a)
      
        op1_price[rho] = vbar
	op2_price[rho] = avgbar





