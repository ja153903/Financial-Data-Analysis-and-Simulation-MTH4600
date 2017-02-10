import numpy as np
import time

k = 110
T = 0.25
r = 0.03
M = 1000
I = 1000
dt = T/I
discount_factor = np.exp(-r * T)

lookback_price = {}
asian_price = {}

rhos = [x for x in np.arange(0, 1, 0.1)] + [0.99999999]

cov_mat = np.zeros((3,3))

xdrift = (r - 0.5 * (0.2 * 0.2)) * dt 
xrand = np.sqrt(dt) * 0.2

ydrift = (r - 0.5 * (0.25 * 0.25)) * dt
yrand = np.sqrt(dt) * 0.25

zdrift = (r - 0.5 * (0.3 * 0.3)) * dt
zrand = np.sqrt(dt) * 0.3

t0 = time.clock()

for rho in rhos:
	
	cov_mat[0, :] = [1, rho, rho]
	cov_mat[1, :] = [rho, 1, rho]
	cov_mat[2, :] = [rho, rho, 1]
	
	cho_mat = np.linalg.cholesky(cov_mat)
	
	lookback_payoff_sum = 0
	asian_payoff_sum = 0
	
	for i in range(0, M):
		uncorr_normal = np.random.standard_normal((3, I))
		chol_normal = np.dot(cho_mat, uncorr_normal)
		
		xnorm = np.cumsum(chol_normal[0])
		ynorm = np.cumsum(chol_normal[1])
		znorm = np.cumsum(chol_normal[2])
	
		sx = 100 * np.exp(xdrift + xrand*xnorm)
		sy = 100 * np.exp(ydrift + yrand*ynorm)
		sz = 100 * np.exp(zdrift + zrand*znorm)
		
		lookback_payoff_sum += max(max(np.max(sx) - 110, np.max(sy)-110, np.max(sz)-110), 0)/M
		asian_payoff_sum += max(np.mean([sx[len(sx)-1], sy[len(sy)-1], sz[len(sy)-1]])-110, 0)/M
	
	lookback_price[rho] = lookback_payoff_sum * discount_factor
	asian_price[rho] = asian_payoff_sum *discount_factor
	

print('Lookback Option: ', lookback_price)
print('\n')
print('Asian Option: ', asian_price)
