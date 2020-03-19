import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from mpl_toolkits import mplot3d


def diagonal(arr, k: int):
    n = len(arr) + abs(k)
    a = np.zeros((n, n))
    if k >= 0:
        for i, e in enumerate(arr):
            a[i][i + k] = e
    else:
        k = abs(k)
        for i, e in enumerate(arr):
            a[i + k][i] = e
    return a


class BSMSolver:
    def __init__(self, option_type='CALL'):
        if option_type not in ['CALL', 'PUT']:
            raise Exception("Type of option should be either CALL or PUT.")
        self.r = 0.05  # interest rate
        self.sigma = 0.2  # volatility
        self.T = 1.0  # maturity time in years
        self.Smin = 0  # min value of the underlying asset
        self.Smax = 100  # max value of the underlying asset
        self.K = 60  # strike price of the underlying asset
        self.N = 100  # number of points in time
        self.M = 200  # number of points in price

        self.ds = (self.Smax - self.Smin) / self.M  # price step
        self.dt = self.T / self.N  # time step
        self.is_call = option_type == "CALL"

    def initial_condition(self, s):
        if self.is_call:
            return max(s - self.K, 0)
        else:
            return max(self.K - s, 0)

    def solve(self) -> np.array:
        np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

        N = self.N  # N is time
        M = self.M  # M is price
        prices = np.zeros((M + 2, N + 2))
        r = self.r
        dt = self.dt

        times = np.array([i * self.dt for i in range(0, N + 1)])
        space_prices = np.array([self.Smin + i * self.ds for i in range(0, M + 1)])

        # Initial Condition
        for i in range(1, M + 2):
            prices[i][-1] = self.initial_condition((i - 1) * self.ds)

        # print(prices[:, N])
        # Boundary Condition
        for i in range(1, N + 2):
            common_factor = (np.e ** (-r * (self.T - (i - 1) * dt)))
            if self.is_call:
                prices[1][i] = 0
                prices[-1][i] = (self.Smax - self.K) * common_factor
            else:
                prices[1][i] = (self.K - self.Smin) * common_factor
                prices[-1][i] = 0

        # coefficients
        sig2 = self.sigma * self.sigma

        a = np.zeros(M + 1)
        b = np.zeros(M + 1)
        c = np.zeros(M + 1)
        for i in range(0, M + 1):
            i2 = i * i
            a[i] = (dt / 4) * (sig2 * i2 - r * i)
            b[i] = -(dt / 2) * (sig2 * i2 + r)
            c[i] = (dt / 4) * (sig2 * i2 + r * i)

        # Matrix formulation
        C = -diagonal(a[2:M], -1) + diagonal(1 - b[1:M], 0) - diagonal(c[1:M - 1], 1)
        D = diagonal(a[2:M], -1) + diagonal(1 + b[1:M], 0) + diagonal(c[1:M - 1], 1)
        (P, L, U) = la.lu(C)
        L = P @ L

        L_INV = np.linalg.inv(L)
        U_INV = np.linalg.inv(U)

        # main loop
        offset = np.zeros(D.shape[1])
        # print("offset shape", offset.shape)
        # print("D.shape", D.shape)

        # print("prices:")
        # print(prices)
        for i in range(N, 0, -1):
            if len(offset) == 1:
                offset = a[2] * (prices[1, i] + prices[1, i + 1]) + \
                         c[-1] * (prices[-1, i] + prices[-1, i + 1])
            else:
                offset[0] = a[2] * (prices[1, i] + prices[1, i + 1])
                offset[-1] = c[-1] * (prices[-1, i] + prices[-1, i + 1])
            px = prices[2:M + 1, i + 1]

            x = D @ px
            x += offset
            x2 = L_INV @ x
            z = U_INV  @ x2

            prices[2:M + 1, i] = z

        # Visualize
        useful_prices = prices[1:, 1:]
        # print(useful_prices)
        # 3D
        plt.figure()
        ax = plt.axes(projection='3d')
        print(times.shape, space_prices.shape, useful_prices.shape)
        X, Y = np.meshgrid(times, space_prices)
        ax.contour3D(X,
                     Y,
                     useful_prices,
                     150,
                     cmap='hot')
        ax.set_title('European Call option price')
        ax.set_xlabel('t')
        ax.set_ylabel('S')
        ax.set_zlabel('C')
        plt.show()
        # 2D
        plt.figure()
        plt.plot(space_prices, useful_prices[:, 0])
        plt.xlabel("St")
        plt.ylabel("Call price")
        plt.show()
        print("ds:", self.ds)
        print("dt:", self.dt)


if __name__ == "__main__":
    print("Start")
    BSMSolver().solve()
    print("End")