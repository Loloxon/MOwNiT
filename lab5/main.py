from math import pi, cos, e, sin
import numpy as np
import matplotlib.pyplot as plot
import numpy.linalg

min_x = -pi
max_x = 3 * pi

graph_counter = 20

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


def showGraph(x, y, nodes_n, degree, graph_counter):
    nodes = genEquidistant(nodes_n)
    axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, aprox(nodes, f(nodes), degree), '-.')
    # axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, aproxtry(nodes, 1), '-.')

    axis[x][y].legend(['Punkty', 'Funkcja', "Przybliżenie"], loc='best')
    axis[x][y].set_xlabel("x")
    axis[x][y].set_ylabel("y")
    axis[x][y].set_title(
        "Wyk. " + str(graph_counter) + "., funkcja aproksymująca " + str(degree) + " stopnia dla " + str(nodes_n) +
        " punktów", fontweight="bold")
    graph_counter += 1
    return graph_counter


# figure, axis = plot.subplots(2)
# figure.suptitle("Wykresy")
# # nodes = genEquidistant(6)
# # aprox(nodes, f(nodes), 2)
# m = 35   # stopien
# n = 50  # wezly
# if(m<=n):
#     showGraph(n, m, graph_counter)
#     plot.show()

def comp_abs(Y):
    ans = 0
    for i in range(0, samples_n):
        idx = min_x + (max_x - min_x) * i / samples_n
        ans = max(ans, abs(Y[i] - f(idx)))
    return ans


min1 = -1

for n in range(1, 50):  # punkty
    for m in range(2, n+1):  # stopień       stopień<=punkty
        nodes = genEquidistant(n)
        d1 = comp_abs(aprox(nodes, f(nodes), m))

        if min1 > d1 or min1 == -1:
            min1 = d1
            n1 = n
            m1 = m

        print("punkty: ", n, ", stopień: ", m, "; ", d1)
        if n%10==0 and m==60:
            print()
            print("Najmniejsze niedokładności otrzymano dla:")
            print(n1, ": ", m1, ": ", min1, "; ")
            print("---------------------")
            figure, axis = plot.subplots(2,2)
            figure.suptitle("Wykresy")
            graph_counter = showGraph(0, 0, n, m, graph_counter)
            graph_counter = showGraph(0, 1, n+10, m, graph_counter)
            graph_counter = showGraph(1, 0, n+20, m, graph_counter)
            graph_counter = showGraph(1, 1, n+30, m, graph_counter)
            # graph_counter = showGraph(i3, 1, 0, graph_counter)
            # graph_counter = showGraph(i5, 0, 1, graph_counter)
            # graph_counter = showGraph(i7, 1, 1, graph_counter)
            plot.show()
