from math import pi, cos, e, sin

import numpy as np
from scipy.special import roots_chebyt

import matplotlib.pyplot as plot

k = 4
m = 2

n = 19
min_x = -pi
max_x = 3 * pi


def f(x):
    if not isinstance(x, float):
        return [e ** (k * cos(m * i)) for i in x]
    return e ** (k * cos(m * x))


def df(x):
    if not isinstance(x, float):
        return [-k * m * (e ** (k * cos(m * i))) * sin(m * i) for i in x]
    return -k * m * (e ** (k * cos(m * x))) * sin(m * x)


# print(f(0.0))
# print(df(0.0))
# print(f(0.1))
# print(df(0.1))
# exit(0)

X = np.arange(min_x, max_x + 0.01, 0.01)
N = len(X)
print("N: ", N)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))


def drawFunction():
    plot.plot(X, f(X), label='Function')
    plot.legend()
    plot.show()


def Hermit(x, x_pos):
    tab = [[0 for i in range(n * 2 + 1)] for j in range(n * 2 + 1)]
    for i in range(1, n * 2 + 1, 2):
        tab[i][0] = x_pos[i // 2]
        tab[i + 1][0] = x_pos[i // 2]
        tab[i][1] = f(x_pos[(i - 1) // 2])
        tab[i + 1][1] = f(x_pos[(i - 1) // 2])
        tab[i + 1][2] = df(x_pos[i // 2])

    for i in range(1, n * 2 + 1):
        for j in range(2, i + 1):
            if i % 2 == 1 or j > 2:
                tab[i][j] = (tab[i][j - 1] - tab[i - 1][j - 1]) / (x_pos[(i - 1) // 2] - x_pos[(i - j) // 2])
    ans = 0
    prod = 1
    for i in range(1, n * 2 + 1):
        ans += prod * tab[i][i]
        prod *= (x - tab[i][0])
    return ans


def drawHermitEquidistant(N):
    if N:
        axis[0].scatter(x_pos, f(x_pos), label='data')
        axis[0].plot(X, Hermit(X, x_pos), label='Hermit interpolated with equidistants')
        axis[0].plot(X, f(X), label='Function')
        axis[0].set_title("Metoda Hermita, punkty równoodległe," + str(N) + " węzłów", fontweight="bold")
        axis[0].set_xlabel("x")
        axis[0].set_ylabel("y")
    return Hermit(X, x_pos)


def drawHermitCheby(N):
    x, _ = roots_chebyt(n)
    x *= ((max_x - min_x) / 2)
    x += ((max_x + min_x) / 2)

    if N:
        axis[1].scatter(x, f(x), label='data')
        axis[1].plot(X, Hermit(X, x), label='Hermit interpolated with Chebyshev')
        axis[1].plot(X, f(X), label='Function')
        axis[1].set_title("Metoda Hermita, punkty Chebysheva," + str(N) + " węzłów", fontweight="bold")
        axis[1].set_xlabel("x")
        axis[1].set_ylabel("y")
    return Hermit(X, x)


def comp_abs(Y):
    ans = 0
    for i in range(0, N):
        idx = min_x + (max_x - min_x) * i / N
        ans = max(ans, abs(Y[i] - f(idx)))
    return ans


def comp_sqr(Y):
    ans = 0
    for i in range(0, N):
        idx = min_x + (max_x - min_x) * i / N
        ans += (Y[i] - f(idx)) ** 2
    return ans


ctr = 30

min1 = min2 = 100000000
i1 = i2 = -1
print("Obliczane różnice pierwszym sposobem")
print("Węzły; Hermit Equidistant; Hermit Chebyshev")
for i in range(3, ctr + 1):
    n = i
    X = np.arange(min_x, max_x + 0.01, 0.01)
    x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

    d1 = comp_abs(drawHermitEquidistant(0))
    d2 = comp_abs(drawHermitCheby(0))

    print(n, end=": ")
    print(d1, end="; ")
    print(d2)

    if min1 > d1:
        min1 = d1
        i1 = n
    if min2 > d2:
        min2 = d2
        i2 = n

print()
print("Najmniejsze niedokładności otrzymano dla:")
print(i1, ": ", min1, "; ", i2, ": ", min2, sep="")
print("---------------------")

figure, axis = plot.subplots(2)

figure.suptitle("Highest difference")

n = i1
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

drawHermitEquidistant(n)

n = i2
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

drawHermitCheby(n)

plot.show()

min1 = min2 = 100000000
i1 = i2 = -1
print("Obliczane różnice drugim sposobem")
print("Węzły; Hermit Equidistant; Hermit Chebyshev")
for i in range(3, ctr + 1):
    n = i
    X = np.arange(min_x, max_x + 0.01, 0.01)
    x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

    d1 = comp_abs(drawHermitEquidistant(0))
    d2 = comp_abs(drawHermitCheby(0))

    print(n, end=": ")
    print(d1, end="; ")
    print(d2)

    if min1 > d1:
        min1 = d1
        i1 = n
    if min2 > d2:
        min2 = d2
        i2 = n

print()
print("Najmniejsze niedokładności otrzymano dla:")
print(i1, ": ", min1, "; ", i2, ": ", min2, sep="")
print("---------------------")

figure, axis = plot.subplots(2)

figure.suptitle("Squared sum")

n = i1
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

drawHermitEquidistant(n)

n = i2
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

drawHermitCheby(n)

plot.show()
