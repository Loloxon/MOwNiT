from math import pi, cos, e

import numpy as np
from scipy.special import roots_chebyt

import matplotlib.pyplot as plot

# 4. Korohoda Nikodem f4, k=4, m=2, [-pi, 3pi]

k = 4
m = 2
# N = 1000

n = 19
min_x = -pi
max_x = 3 * pi


def f(x):
    if not isinstance(x, float):
        return [e**(k*cos(m*i)) for i in x]
    return e**(k*cos(m*x))

X = np.arange(min_x, max_x + 0.01, 0.01)
N = len(X)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n-1))

def drawFunction():
    # f(X) = np.vectorize(f(X))
    plot.plot(X, f(X), label='Function')
    plot.legend()
    plot.show()

def Lagrange(x, x_pos):
    def L(k, x):
        prod = 1
        for i in range(0, n):
            if i != k:
                prod *= (x - x_pos[i]) / (x_pos[k] - x_pos[i])
        return prod

    sum = 0
    for k in range(0, n, 1):
        sum += (f(x_pos[k]) * L(k, x))
    return sum


def drawLagrangeEquidistant(show):
    if show:
        axis[0,0].scatter(x_pos, f(x_pos), label='data')
        axis[0,0].plot(X, Lagrange(X, x_pos), label='Lagrange interpolated with equidistants')
        axis[0,0].plot(X, f(X), label='Function')
        axis[0,0].set_title("L, EQ")
        # axis[,].show()
    return Lagrange(X, x_pos)


def drawLagrangeCheby(show):
    x, _ = roots_chebyt(n)
    x *= ((max_x - min_x) / 2)
    x += ((max_x + min_x) / 2)

    if show:
        axis[0,1].scatter(x, f(x), label='data')
        axis[0,1].plot(X, Lagrange(X, x), label='Lagrange interpolated with Chebyshev')
        axis[0,1].plot(X, f(X), label='Function')
        axis[0,1].set_title("L, Cheb")
        # axis[,].show()
    return Lagrange(X, x)


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
    for k in range(0, n):
        prod = b(k)
        for i in range(0, k):
            prod *= (x - x_pos[i])
        suma += prod
    return suma


def drawNewtonEquidistant(show):
    if show:
        axis[1,0].scatter(x_pos, f(x_pos), label='data')
        axis[1,0].plot(X, Newton(X, x_pos), label='Newton interpolated with equidistants')
        axis[1,0].plot(X, f(X), label='Function')
        axis[1,0].set_title("N, EQ")
        # axis[,].show()
    return Newton(X, x_pos)


def drawNewtonCheby(show):
    x, _ = roots_chebyt(n)
    x *= ((max_x - min_x) / 2)
    x += ((max_x + min_x) / 2)

    if show:
        axis[1,1].scatter(x, f(x), label='data')
        axis[1,1].plot(X, Newton(X, x), label='Newton interpolated with Chebyshev')
        axis[1,1].plot(X, f(X), label='Function')
        axis[1,1].set_title("N, Cheb")
        # axis[,].show()
    return Newton(X,x)

def comp_abs(Y):
    ans = 0
    for i in range(0, N):
        idx = min_x + (max_x - min_x)*i / N
        ans = max(ans, abs(Y[i]-f(idx)))
    return ans

def comp_sqr(Y):
    ans = 0
    for i in range(0, N):
        idx = min_x + (max_x - min_x)*i / N
        # ans = max(ans, abs(Y[i]-f(idx)))
        ans += (Y[i] - f(idx))**2
    return ans

ctr = 30

min1=min2=min3=min4=100000000
i1=i2=i3=i4=-1

for i in range(3, ctr+1):
    n=i
    X = np.arange(min_x, max_x + 0.01, 0.01)
    x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n-1))
    
    d1 = comp_abs(drawLagrangeEquidistant(0))
    d2 = comp_abs(drawLagrangeCheby(0))
    d3 = comp_abs(drawNewtonEquidistant(0))
    d4 = comp_abs(drawNewtonCheby(0))
    
    print(n, end=": ")
    print(d1, end="; ")
    print(d2, end="; ")
    print(d3, end="; ")
    print(d4)
    
    if min1>d1:
        min1=d1
        i1 = n
    if min2>d2:
        min2=d2
        i2 = n
    if min3>d3:
        min3=d3
        i3 = n
    if min4>d4:
        min4=d4
        i4 = n

print()
print(i1,":", min1, ", ",i2,":", min2,", ",i3,":", min3,", ",i4,":", min4, sep="")
print("---------------------")

figure, axis = plot.subplots(2, 2)

n = i1
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n-1))

drawLagrangeEquidistant(1)
drawNewtonEquidistant(1)

n = i2
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n-1))

drawLagrangeCheby(1)
drawNewtonCheby(1)
plot.show()



min1=min2=min3=min4=100000000
i1=i2=i3=i4=-1

for i in range(4, ctr + 2):
    n = i
    X = np.arange(min_x, max_x + 0.01, 0.01)
    x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n-1))

    d1 = comp_sqr(drawLagrangeEquidistant(0))
    d2 = comp_sqr(drawLagrangeCheby(0))
    d3 = comp_sqr(drawNewtonEquidistant(0))
    d4 = comp_sqr(drawNewtonCheby(0))

    print(n+1, end=": ")
    print(d1, end="; ")
    print(d2, end="; ")
    print(d3, end="; ")
    print(d4)

    if min1 > d1:
        min1 = d1
        i1 = n
    if min2 > d2:
        min2 = d2
        i2 = n
    if min3 > d3:
        min3 = d3
        i3 = n
    if min4 > d4:
        min4 = d4
        i4 = n

print()
print(i1, ":", min1,", ", i2, ":", min2,", ", i3, ":", min3,", ", i4, ":", min4, sep="")
print("---------------------")



# drawFunction()

figure, axis = plot.subplots(2, 2)

n = i1
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n-1))

drawLagrangeEquidistant(1)
drawNewtonEquidistant(1)

n = i2
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n-1))

drawLagrangeCheby(1)
drawNewtonCheby(1)
plot.show()