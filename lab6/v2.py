import math
from copy import deepcopy
from math import pi, e

import matplotlib.pyplot as plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

pd.set_option("display.precision", 3)

graph_counter = 15


def f(samples):
    k = 4
    m = 2
    return e ** (k * math.cos(m * samples))

f = np.vectorize(f)
a = -1 * pi
b = 3 * pi
x_space = np.arange(a, b + 0.01, 0.01)

def plott(space, *functions, points=None, title=None):
    plt.rcParams['figure.figsize'] = [9, 6]

    if points != None:
        plt.scatter(points[0], points[1], label="Punkty")
    for foo, lbl, line in functions:
        plt.plot(space, foo(space), line, label=lbl)

    if title:
        plt.title(title, fontweight="bold")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend(bbox_to_anchor=(0.85, 0.23), borderaxespad=0)
    plt.grid()
    plt.show()


def get_norm(y1, y2, mode):
    n = len(y1)
    if mode == 'max':
        return max([abs(y1[i] - y2[i]) for i in range(n)])
    elif mode == 'sse':
        return sum([(y1[i] - y2[i]) ** 2 for i in range(n)])


class TrigonometricApproximation:
    def __init__(self, m, n, start, end):
        self.__m = m
        self.__n = n
        self.__xs = np.linspace(start, end, n)
        self.__ys = f(xs)
        self.range_start = start
        self.range_end = end
        self.rng_len = end - start
        self.__a = np.zeros(self.__m + 1)
        self.__b = np.zeros(self.__m + 1)
        self.__scale_xs_to_2pi()
        self.__fill_a_b()
        self.__scale_xs_from_2pi()

    def __scale_xs_to_2pi(self):
        for i in range(self.__n):
            self.__xs[i] = self.__scale_to_2pi(self.__xs[i])

    def __scale_xs_from_2pi(self):
        for i in range(self.__n):
            self.__xs[i] = self.__scale_from_2pi(self.__xs[i])

    def __scale_to_2pi(self, x):
        return (x - pi) / 2

    def __scale_from_2pi(self, x):
        return x * 2 + pi



    def __fill_a_b(self):
        for j in range(self.__m + 1):
            self.__a[j] = 2 * sum(self.__ys[i] * np.cos(j * self.__xs[i]) for i in range(self.__n)) / self.__n
            self.__b[j] = 2 * sum(self.__ys[i] * np.sin(j * self.__xs[i]) for i in range(self.__n)) / self.__n

    def approximate(self, argX):
        x = deepcopy(argX)
        x = self.__scale_to_2pi(x)
        return (self.__a[0]/2) + sum(
            [self.__a[j] * np.cos(j * x) + self.__b[j] * np.sin(j * x) for j in range(1, self.__m + 1)]
        )


AR = []
R = []
C = []
c_flag = 0

best = -1
for m in range(30, 1, -1):  # stopien
    AR.append([])
    R.append(m)
    for n in range(3, 81):  # l. punktów
        if c_flag == 0:
            C.append(n)
        if (n <= m * 2):
            AR[len(AR) - 1].append(0)
        else:
            xs = np.linspace(a, b, n)

            ys = f(xs)

            ta = TrigonometricApproximation(m, n, a, b)
            F = ta.approximate

            diff = get_norm(F(x_space), f(x_space), 'max')
            AR[len(AR) - 1].append(diff)
            # plot(x_space, [f, "f", "g-"], [F, "F", "b-"], points=[xs, f(xs)], title=f"n={n}  m={m}")
            if best == -1 or best > diff:
                best = diff
                m1 = m
                n1 = n
            print(m, n)

            # if (m==10 or m==15) and (n%40==0):
            #     plott(x_space, [f, "Funkcja", "g-"], [F, "Przybliżenie", "r-"], points=[xs, f(xs)],
            #           title=f"Wyk. " + str(graph_counter)
            #                 + "., funkcja aproksymująca " + str(m) + " stopnia dla " + str(n) +
            #                 " punktów")
            #     graph_counter += 1
    c_flag = 1

    # if m==2:
xs = np.linspace(a, b, n1)

ys = f(xs)

ta = TrigonometricApproximation(m1, n1, a, b)
F = ta.approximate
print(m1, n1)
print(best)
print(get_norm(F(x_space), f(x_space), 'max'))
plott(x_space, [f, "Funkcja", "g-"], [F, "Przybliżenie", "r-"], points=[xs, f(xs)], title=f"Wyk. " + str(graph_counter)
    + "., funkcja aproksymująca " + str(m1) + " stopnia dla " + str(n1) + " punktów")
graph_counter+=1


print(AR)
df = pd.DataFrame(AR, R, C)

figure, axis = plot.subplots(1)

# Default heatmap
# cmap = sns.cm.rocket_r
hm = sns.heatmap(df, square=True, cmap="Blues", cbar_kws=dict(use_gridspec=False, location="bottom"))
hm.set(xlabel='Liczba punktów', ylabel='Stopień wielomianu')

# plot.show()
# plott(hm)
