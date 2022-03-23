import numpy as np
from scipy.interpolate import lagrange
from scipy.special import roots_chebyt
from scipy.special import eval_chebyt

import matplotlib.pyplot as plot

n = 7
min_x = -2
max_x = 2

n = n - 1


def f(x):
    return 2 * x ** 3 - x ** 2 + 3 * x + 6


X = np.arange(min_x, max_x + 0.01, 0.1)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / n)
print(X)
print(x_pos)


def drawFunction():
    plot.plot(X, f(X), label='Function')
    plot.legend()
    plot.show()


def drawExpected():
    def SciPyLagrange(x):
        return lagrange(x, f(x))

    scipyLagrange = SciPyLagrange(x_pos)

    plot.scatter(x_pos, f(x_pos), label='data')
    plot.plot(X, scipyLagrange(X), label='Lagrange using SciPy')
    plot.legend()
    plot.show()


def Lagrange(x, x_pos):
    def L(k, x):
        prod = 1
        for i in range(0, n + 1):
            if i != k:
                prod *= (x - x_pos[i]) / (x_pos[k] - x_pos[i])
        return prod

    sum = 0
    for k in range(0, n + 1, 1):
        sum += (f(x_pos[k]) * L(k, x))
    return sum


def drawLagrangeEquidistant():
    plot.scatter(x_pos, f(x_pos), label='data')
    plot.plot(X, Lagrange(X, x_pos), label='Lagrange interpolated with equidistants')
    plot.legend()
    plot.show()


def drawLagrangeCheby():
    x, _ = roots_chebyt(n + 1)
    x += ((max_x + min_x) / 2)
    x *= ((max_x - min_x) / 2)

    plot.scatter(x, f(x), label='data')
    plot.plot(X, Lagrange(X, x), label='Lagrange interpolated with Chebyshev')
    plot.legend()
    plot.show()


def Newton(x, x_pos):
    # def fnot(b, e):
    #     if b == e:
    #         return f(x_pos[b])
    #     if b + 1 == e:
    #         return (f(x_pos[e] - f(x_pos[b]))) / (x_pos[e] - x_pos[b])
    #     return (fnot(b + 1, e) - fnot(b, e - 1)) / (x_pos[e] - x_pos[b])
    #
    # sum = fnot(0, 0)
    # for k in range(1, n + 1):
    #     prod = fnot(0, k)
    #     for i in range(0, k):
    #         prod *= (x - x_pos[i])
    #     sum += prod
    # return sum
    def b(k):
        suma = 0
        for i in range(0, k + 1):
            prod = 1
            for j in range(0, k + 1):
                if j != i:
                    prod *= (x_pos[i] - x_pos[j])
            suma += (f(x_pos[i]) / prod)
        return suma

    suma = 0
    for k in range(0, n + 1):
        prod = b(k)
        for i in range(0, k):
            prod *= (x - x_pos[i])
        suma += prod
    return suma


def drawNewtonEquidistant():
    plot.scatter(x_pos, f(x_pos), label='data')
    plot.plot(X, Newton(X, x_pos), label='Newton interpolated with equidistants')
    plot.legend()
    plot.show()


def drawNewtonCheby():
    x, _ = roots_chebyt(n + 1)
    x += ((max_x + min_x) / 2)
    x *= ((max_x - min_x) / 2)

    plot.scatter(x, f(x), label='data')
    plot.plot(X, Newton(X, x), label='Newton interpolated with Chebyshev')
    plot.legend()
    plot.show()


drawFunction()
drawExpected()
drawLagrangeEquidistant()
drawLagrangeCheby()
drawNewtonEquidistant()
drawNewtonCheby()
