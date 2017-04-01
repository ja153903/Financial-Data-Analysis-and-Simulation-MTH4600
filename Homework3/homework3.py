import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import norm


class GARCH11:

    def __init__(self):
        self.S0 = 100 # initial stock price
        self.r = 0.05 # risk-free rate
        self.T = 0.5 # maturity
        self.K = [50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150] # list of strikes

        self.long_run_vol_squared = 0.3 * 0.3 # long-run volatility
        self.stochastic_vol_squared = 0.35 * 0.35 # initial stochastic volatility

        self.alpha = 0.080 # weight for long-run vol
        self.beta = 0.885 # weight for previous volatility
        self.gamma = 0.035 # weight for current return * return
        self.steps = 252 # number of steps for gbm
        self.dt = self.T / self.steps
        self.square_root_dt = np.sqrt(self.dt)

    def generatePathMatrix(self, sim):
        '''
		This function returns an array of paths for the stock
		and the stochastic volatility squared

		To decrease the time complexity of the Monte Carlo simulation,
		we can take advantage of Numpy's vectorization and simulatenously 
		run the simulations in a multidimensional array.
		'''

        # initialize 1000 arrays with N periods
        St = np.zeros((sim, self.steps))

        # intialize the period to contain the first stock price
        St[:,0] = self.S0

        # initialize 1000 arrays with N periods
        stochastic_vol = np.zeros((sim, self.steps))

        # initialize the first volatility to be 0.35 * 0.35
        stochastic_vol[:, 0] = self.stochastic_vol_squared

        # initialize 1000 arrays with N periods to 0
        B = np.zeros((sim, self.steps))

        # initialize 1000 arrays with N standard normals
        std_normal = np.random.standard_normal((sim, self.steps))

        # Now we just need to generate the path
        for i in range(1, self.steps):
            mu = (self.r - 0.5 * stochastic_vol[:, i-1]) * self.dt
            B[:, i] = self.square_root_dt * np.sqrt(stochastic_vol[:, i-1]) * std_normal[:, i-1]
            St[:, i] = St[:, i-1] * np.exp(mu + B[:, i])
            stochastic_vol[:, i] = self.gamma * self.long_run_vol_squared + self.beta * stochastic_vol[:, i-1] + self.alpha * B[:, i] * B[:, i] / self.dt

        return St, stochastic_vol

    def priceOption(self, price_path, K):
        '''
        This function returns an array of call option prices
        depending on different strikes
        '''
        payoff_sum = np.sum(np.exp(-self.r * self.T) * np.maximum(price_path[:,self.steps - 1] - K, 0.0))

        return payoff_sum / len(price_path[:,self.steps-1])

    def blackScholesFormula(self, T, S0, K, sigma, r):
        '''
        Returns option price using Black-Scholes Formula
		'''
        v = r * np.sqrt(T)

        if v < 0.0000001:
            S0 = S0 * np.exp(r * T)
            if S0 > K:
               value = S0 - K
            else:
               value = 0
        else:
            d0 = (np.log(S0/K) + r * T) / v
            d_plus = d0 + 0.5 * v
            d_minus = d0 - 0.5 * v

            # find an alternative to this to calculate the standard normal CDF
            value = S0 * norm.cdf(d_plus) - K * np.exp(-r * T) * norm.cdf(d_minus)

        return value if value > 0 else 0

    def impliedVolatility(self, K, call_price):

        sigma = 0
        step = 0.1

        init_value = self.blackScholesFormula(self.T, self.S0, K, 0.0, self.r)

        if init_value > call_price:
            return -1

        if call_price < init_value + 0.0000001:
            return 0.0

        value = init_value

        for i in range(1, 6):

            while value <= call_price:

                sigma = sigma + step
                value = self.blackScholesFormula(self.T, self.S0, K, sigma, self.r)

            sigma = sigma - step
            value = self.blackScholesFormula(self.T, self.S0, K, sigma, self.r)

            step = step / 10

        return sigma

if __name__ == '__main__':

    print("Hello World")

    test1 = GARCH11()

    # test with 1000 simulations
    stock_price_path, stochastic_vol_squared = test1.generatePathMatrix(1000)

    call_prices = {}
    for x in test1.K:
        call_prices[x] = test1.priceOption(stock_price_path, x)

    for x in test1.K:
        print(test1.impliedVolatility(x, call_prices[x]))