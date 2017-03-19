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
I = 252 # number of samples
dt = T/I # increments of steps
discount_factor = np.exp(-r * T) # discount factor for option pricing

op1_price = {} # dictionary for prices by the rho value
op2_price = {} # dictionary for prices by the rho value

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
This function returns the payoff of the op1 option
in 1(a) divided by the number of simulations.
'''
def op1Payoff(s0, s1, s2):
	lp = max(s0.max(), s1.max(), s2.max()) - k
	return lp/M if lp > 0 else 0

'''
This function returns the payoff of the op2 option
in 1(b) divided by the number of simulations
'''
def op2Payoff(s0, s1, s2):
	ap = np.average([s0[len(s0)-1], s1[len(s1)-1], s2[len(s2)-1]])-k
	return ap/M if ap > 0 else 0


for rho in rhos: # we loop through the values of rho
	
	'''
	This calculates the Cholesky Decomposed matrix of the covariance matrix
	'''
	cho_mat = np.linalg.cholesky([[1, rho, rho], [rho, 1, rho], [rho, rho, 1]])
	
	op1_payoff_sum = op2_payoff_sum = 0 # initializing the running sum
	
	
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
		
		
		op1_payoff_sum += op1Payoff(sx, sy, sz) # collect the sum of the payoffs 
		
		op2_payoff_sum += op2Payoff(sx, sy, sz) # collect the sum of the payoffs
	
	'''
	The price of the option can be calculated discounting the average payoff
	After that, we enter the rho and price into the dictionary as a key-value pair
	'''
	op1_price[rho] = op1_payoff_sum * discount_factor
	op2_price[rho] = op2_payoff_sum * discount_factor


print('op1 Option: ')
print([op1_price[x] for x in rhos])
print('\n')
print('op2 Option: ')
print([op2_price[x] for x in rhos])

plt.scatter(list(op1_price.keys()), list(op1_price.values()), label='op1')
plt.title('1(a) and 1(b) Option')
plt.xlabel('Rho')
plt.ylabel('Price of the Option')
plt.xlim([-0.1, 1.1])
plt.scatter(list(op2_price.keys()), list(op2_price.values()), color='red', label='op2')
plt.legend()
plt.show()





