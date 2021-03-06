import math
from math import pi, cos, e, sin
import numpy as np
import matplotlib.pyplot as plot
import numpy.linalg

import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm

min_x = -pi
max_x = 3 * pi

graph_counter = 1

samples = np.arange(min_x, max_x + 0.01, 0.01)
samples_n = len(samples)


def f(samples):
    k = 4
    m = 2
    if not isinstance(samples, float):
        return [e ** (k * cos(m * i)) for i in samples]
    return e ** (k * cos(m * samples))


def genEquidistant(nodes_n):
    return np.arange(min_x, max_x + 0.01, (max_x - min_x) / (nodes_n - 1))


def aprox(nodes, values, degree):
    S = []
    for i in range(degree * 2 + 1):
        tmp = 0
        for n in nodes:
            tmp += (n ** i)
        S.append(tmp)
    T = []
    for i in range(degree + 1):
        tmp = 0
        for n in range(len(nodes)):
            tmp += (values[n] * nodes[n] ** i)
        T.append(tmp)
    L = []
    for i in range(degree + 1):
        L.append([])
        for j in range(degree + 1):
            L[i].append(S[i + j])
    M = np.linalg.solve(L, T)
    # for i in L:
    #     print(i)
    # print(T)
    # print(M)
    ANS = []
    for s in samples:
        tmp = 0
        for j in range(degree + 1):
            tmp += M[j] * s ** j
        ANS.append(tmp)
    return ANS


def aproxtry(nodes, values, degree):
    # global a0
    # degree = len(nodes) // 2
    nodes_n = len(nodes)
    # nodes_n2 = len(nodes)//2
    a0 = 0
    for i in values:
        a0 += i
    a0 /= nodes_n
    A = []
    for j in range(degree):
        tmp = 0
        for i in range(nodes_n):
            # tmp += values[i] * cos(2 * pi * i * j / nodes_n)
            tmp += values[i] * cos((j+1)*nodes[i])
        tmp *= (2 / nodes_n)
        A.append(tmp)
    B = []
    for j in range(degree):
        tmp = 0
        for i in range(nodes_n):
            # tmp += values[i] * sin(2 * pi * i * j / nodes_n)
            tmp += values[i] * sin((j+1)*nodes[i])
        tmp *= (2 / nodes_n)
        B.append(tmp)
    ANS = []
    for s in samples:
        tmp = a0
        for i in range(degree):
            tmp += A[i] * cos((i + 1) * s) #(s-min_x)*pi*(2-1/nodes_n)/(max_x-min_x)-pi
            # tmp += A[i] * cos((i + 1) * ((s-min_x)*pi*(2-2/nodes_n)/(max_x-min_x)-pi)) #(s-min_x)*pi*(2-1/nodes_n)/(max_x-min_x)-pi
            # tmp += A[i] * cos(s)**(i + 1)
        for i in range(degree):
            tmp += B[i] * sin((i + 1) * s)
            # tmp += B[i] * sin((i + 1) * (s-min_x)*pi*(2-2/nodes_n)/(max_x-min_x)-pi)
            # tmp += B[i] * sin(s)**(i + 1)
        ANS.append(tmp)
    return ANS



def showGraph(x, y, nodes_n, degree, graph_counter):
    nodes = genEquidistant(nodes_n)
    # axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, aprox(nodes, f(nodes), degree), '-.')
    axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, aproxtry(nodes, f(nodes), degree), '-.')

    axis[x][y].legend(['Punkty', 'Funkcja', "Przybli??enie"], loc='best')
    axis[x][y].set_xlabel("x")
    axis[x][y].set_ylabel("y")
    figure.suptitle(
        "Wyk. " + str(graph_counter) + "., funkcja aproksymuj??ca " + str(degree) + " stopnia dla " + str(nodes_n) +
        " punkt??w", fontweight="bold")
    # axis[x][y].set_title(
    #     "Wyk. " + str(graph_counter) + "., funkcja aproksymuj??ca " + str(degree) + " stopnia dla " + str(nodes_n) +
    #     " punkt??w", fontweight="bold")
    graph_counter += 1
    return graph_counter



# figure, axis = plot.subplots(2,2)
# figure.suptitle("Wykresy")
# # nodes = genEquidistant(6)
# # aprox(nodes, f(nodes), 2)
# m = 35   # stopien
# n = 50  # wezly
# if(m<=n):
# showGraph(0, 0, 12, 27, graph_counter)
# plot.show()

def comp_abs(Y):
    ans = 0
    for i in range(0, samples_n):
        idx = min_x + (max_x - min_x) * i / samples_n
        ans = max(ans, abs(Y[i] - f(idx)))
    return ans


min1 = -1

AR = []
R = []
C = []
c_flag = 0
for m in range(30, 1, -1):  # stopie??
    AR.append([])
    R.append(m)
    # for n in range(m+1, m + 16):  # stopie??       stopie??<punkty/2
    for n in range(3, 80):  # punkty       stopie??<punkty/2
        if c_flag == 0:
            C.append(n)
        if (n <= m*2):
            AR[len(AR) - 1].append(0)
        else:
            nodes = genEquidistant(n)
            d1 = comp_abs(aproxtry(nodes, f(nodes), m))

            AR[len(AR) - 1].append(d1)
            if min1 > d1 or min1 == -1:
                min1 = d1
                n1 = n
                m1 = m

            print("punkty: ", n, ", stopie??: ", m, "; ", d1)
            # if n%10==0 and m==60:
    if m==2:
    #     print()
    #     print("Najmniejsze niedok??adno??ci otrzymano dla:")
    #     print(n1, ": ", m1, ": ", min1, "; ")
    #     print("---------------------")
        figure, axis = plot.subplots(2, 2)
        figure.suptitle("Wykresy")
        graph_counter = showGraph(0, 0, n1, m1, graph_counter)
    # graph_counter = showGraph(0, 1, n+10, m, graph_counter)
    # graph_counter = showGraph(1, 0, n+20, m, graph_counter)
    # graph_counter = showGraph(1, 1, n+30, m, graph_counter)
    # graph_counter = showGraph(i3, 1, 0, graph_counter)
    # graph_counter = showGraph(i5, 0, 1, graph_counter)
    # graph_counter = showGraph(i7, 1, 1, graph_counter)
        plot.show()
    c_flag = 1

# Create a dataset
# ar = np.array([[1, 2], [3, 4]])
# d = {'col1': [1, 2], 'col2': [3, 4], 'col3': [3,0]}
# df = pd.DataFrame(np.random.random((2,2)))
# df = pd.DataFrame(ar)
print(AR)
df = pd.DataFrame(AR, R, C)

# Default heatmap
# cmap = sns.cm.rocket_r
hm = sns.heatmap(df, square=True, cmap="Blues", cbar_kws = dict(use_gridspec=False,location="bottom"))
hm.set(xlabel='Liczba punkt??w', ylabel='Stopie?? wielomianu')

plot.show()
