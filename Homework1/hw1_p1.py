'''
MTH 4600: Data Analysis and Simulation
Homework 1
Question 1(a) and (b)
Group: Jaime Abbariao, Bell Chen, Jonnathan Romero
'''
import numpy as np
import time
import matplotlib.pyplot as plt

s0 = 100 # initial price
k = 110 # strike price
T = 0.25 # time to expiration (3 months)
r = 0.03 # interest rate
M = 1000 # number of simulations
I = 1000 # number of samples
dt = T/I # increments of steps
discount_factor = np.exp(-r * T) # discount factor for option pricing

lookback_price = {} # dictionary for prices by the rho value
asian_price = {} # dictionary for prices by the rho value

# generate a list of rho values from 0 to 1 by steps of 0.1
rhos = [x for x in np.arange(0, 1, 0.1)] + [0.99999999]

'''
To ease the calculation of the stock pricing,
we partition parts of the Geometric Brownian Motion calculation here
'''
xdrift = (r - 0.5 * (0.2 * 0.2)) * dt  
xrand = np.sqrt(dt) * 0.2

ydrift = (r - 0.5 * (0.25 * 0.25)) * dt
yrand = np.sqrt(dt) * 0.25

zdrift = (r - 0.5 * (0.3 * 0.3)) * dt
zrand = np.sqrt(dt) * 0.3

'''
Modular version of the path simulation
to make code more readable
This function returns an array of the values 
of the price simulation.
'''
def simulatePath(drift, rand, norm):
	return s0 * np.exp(drift + rand*norm)

'''
This function returns the payoff of the lookback option
in 1(a) divided by the number of simulations.
'''
def lookbackPayoff(s0, s1, s2):
	lp = max(s0.max(), s1.max(), s2.max()) - 110
	return lp/M if lp > 0 else 0

'''
This function returns the payoff of the asian option
in 1(b) divided by the number of simulations
'''
def asianPayoff(s0, s1, s2):
	ap = np.average([s0[len(s0)-1], s1[len(s1)-1], s2[len(s2)-1]])-110
	return ap/M if ap > 0 else 0


for rho in rhos: # we loop through the values of rho
	
	'''
	This calculates the Cholesky Decomposed matrix of the covariance matrix
	'''
	cho_mat = np.linalg.cholesky([[1, rho, rho], [rho, 1, rho], [rho, rho, 1]])
	
	lookback_payoff_sum = asian_payoff_sum = 0 # initializing the running sum
	
	
	for i in range(0, M): # loop through simulations
		
		uncorr_normal = np.random.standard_normal((3, I)) # generate 3 random Gaussian vectors
		chol_normal = np.dot(cho_mat, uncorr_normal) # dot product produces 3 correlated Gaussian vectors
		
		'''
		To speed up the simulation of the path, we take the cumulative sum
		of the correlated Gaussian vectors.
		'''
		xnorm = np.cumsum(chol_normal[0])
		ynorm = np.cumsum(chol_normal[1])
		znorm = np.cumsum(chol_normal[2])
		
		
		'''
		We simulate the 3 paths here.
		'''
		sx = simulatePath(xdrift, xrand, xnorm)
		sy = simulatePath(ydrift, yrand, ynorm)
		sz = simulatePath(zdrift, zrand, znorm)
		
		
		lookback_payoff_sum += lookbackPayoff(sx, sy, sz) # collect the sum of the payoffs 
		
		asian_payoff_sum += asianPayoff(sx, sy, sz) # collect the sum of the payoffs
	
	'''
	The price of the option can be calculated discounting the average payoff
	After that, we enter the rho and price into the dictionary as a key-value pair
	'''
	lookback_price[rho] = lookback_payoff_sum * discount_factor
	asian_price[rho] = asian_payoff_sum *discount_factor


print('Lookback Option: ')
print([lookback_price[x] for x in rhos])
print('\n')
print('Asian Option: ')
print([asian_price[x] for x in rhos])

plt.scatter(list(lookback_price.keys()), list(lookback_price.values()), label='Lookback')
plt.title('1(a) and 1(b) Option')
plt.xlabel('Rho')
plt.ylabel('Price of the Option')
plt.xlim([-0.1, 1.1])
plt.scatter(list(asian_price.keys()), list(asian_price.values()), color='red', label='Asian')
plt.legend()
plt.show()





