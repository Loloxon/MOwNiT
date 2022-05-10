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

X = np.arange(min_x, max_x + 0.01, 0.01)
N = len(X)
print("N: ", N)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

wykres_licznik=8

def drawFunction():
    plot.plot(X, f(X), label='Funkcja')
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
        axis[0].scatter(x_pos, f(x_pos), label='Węzły')
        axis[0].plot(X, f(X), label='Funkcja')
        axis[0].plot(X, Hermit(X, x_pos), label='Metoda Hermita z punktami równoodległymi')
        axis[0].plot(X, Lagrange(X, x_pos), label='Metoda Lagrange z punktami równoodległymi')
        axis[0].set_title("Wyk. "+str(wykres_licznik)+"., Metoda Hermita, punkty równoodległe, " + str(N) + " węzłów", fontweight="bold")
        axis[0].set_xlabel("x")
        axis[0].set_ylabel("y")
        axis[0].legend()
    return Hermit(X, x_pos)


def drawHermitCheby(N):
    x, _ = roots_chebyt(n)
    x *= ((max_x - min_x) / 2)
    x += ((max_x + min_x) / 2)

    if N:
        axis[1].scatter(x, f(x), label='Węzły')
        axis[1].plot(X, f(X), label='Funkcja')
        axis[1].plot(X, Hermit(X, x), label='Metoda Hermita z punktami Chebysheva')
        axis[1].plot(X, Lagrange(X, x), label='Metoda Lagrange z punktami Chebysheva')
        axis[1].set_title("Wyk. "+str(wykres_licznik)+"., Metoda Hermita, punkty Chebysheva, " + str(N) + " węzłów", fontweight="bold")
        axis[1].set_xlabel("x")
        axis[1].set_ylabel("y")
        axis[1].legend()
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


ctr = 80

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

n = 23
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

drawHermitEquidistant(n)
wykres_licznik += 1

n = 23
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

drawHermitCheby(n)
wykres_licznik += 1

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
wykres_licznik += 1

n = i2
X = np.arange(min_x, max_x + 0.01, 0.01)
x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

drawHermitCheby(n)
wykres_licznik += 1

plot.show()
