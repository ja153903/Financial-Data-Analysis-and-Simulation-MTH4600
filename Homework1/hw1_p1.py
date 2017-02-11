'''
MTH 4600 - Data Analysis and Simulation
Homework 1
Group Members: Jaime Abbariao, Bell Chen, Jonnathan Romero
'''

import numpy as np
from scipy import linalg
'''
Question 1

Price Path Simulation
'''

k = 110
T = .25
r = 0.03
N = 1000
M = 50
dt = T/N

xdrift = (r - 0.5 * (0.2*0.2))*dt
xrand = (np.sqrt(dt) * 0.2)

ydrift = (r - 0.5  * (0.25*.25))*dt
yrand = np.sqrt(dt) * 0.25

zdrift = (r - 0.5 * (0.3 * 0.3))*dt
zrand = np.sqrt(dt) * 0.3

rhos = [x for x in np.arange(0, 1, 0.1)] + [0.999999]

cov_mat = np.zeros((3,3))

for rho in rhos:

	cov_mat[0, :] = [1, rho, rho]
	cov_mat[1, :] = [rho, 1, rho]
	cov_mat[2, :] = [rho, rho, 1]

	chol_mat = linalg.cholesky(cov_mat)

	for i in range(0, M):

		corr_normal = np.random.standard_normal((3, N))
		chol_normal = np.dot(chol_mat, corr_normal)

		sx = np.zeros_like(corr_normal[0])
		sy = np.zeros_like(corr_normal[0])
		sz = np.zeros_like(corr_normal[0])

		xnorm = np.cumsum(chol_normal[0])
		ynorm = np.cumsum(chol_normal[1])
		znorm = np.cumsum(chol_normal[2])
		
		sx = 100 * np.exp(xdrift+ xrand * xnorm)
		sy = 100 * np.exp(ydrift + yrand * ynorm)
		sz = 100 * np.exp(zdrift + zrand * znorm)





