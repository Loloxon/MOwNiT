from math import pi, cos, e, sin
import numpy as np
import matplotlib.pyplot as plot
import numpy.linalg
from scipy.special import roots_chebyt

# nodes_n = 10
min_x = -pi
max_x = 3 * pi

graph_counter = 10

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


def genChebyshev(nodes_n):
    nodes, _ = roots_chebyt(nodes_n)
    nodes *= ((max_x - min_x) / 2)
    nodes += ((max_x + min_x) / 2)
    return nodes


def cubicSpline(nodes, mode):
    y = f(nodes)
    nodes_n = len(nodes)

    def h(i):
        return nodes[i] - nodes[i - 1]

    def delt(i):
        return (y[i] - y[i - 1]) / (nodes[i] - nodes[i - 1])

    def delt2(i):
        return (delt(i + 1) - delt(i)) / (nodes[i + 1] - nodes[i - 1])

    def delt3(i):
        return (delt2(i + 1) - delt2(i)) / (nodes[i + 2] - nodes[i - 1])

    A = [[0 for i in range(nodes_n)] for j in range(nodes_n)]
    P = [0 for _ in range(nodes_n)]
    if mode == 0:
        A[0][0] = 2
        A[0][1] = 1
    else:
        A[0][0] = 1
        A[0][1] = 0
    for i in range(1, nodes_n - 1):  # 1 brzegi3
        A[i][i - 1] = h(i)
        A[i][i] = 2 * (h(i) + h(i + 1))
        A[i][i + 1] = h(i + 1)

    if mode == 0:
        A[nodes_n - 1][nodes_n - 2] = 2
        A[nodes_n - 1][nodes_n - 1] = 1
    else:
        A[nodes_n - 1][nodes_n - 2] = 0
        A[nodes_n - 1][nodes_n - 1] = 1
    if mode==0:
        P[0] = 0
    else:
        P[0] = 0
    for i in range(1, nodes_n - 1):
        P[i] = delt(i + 1) - delt(i)
    if mode==0:
        P[nodes_n - 1] = 0
    else:
        P[nodes_n-1] = 0
    ro = np.linalg.solve(A, P)

    def b(i):
        return (y[i] - y[i - 1]) / h(i) - h(i) * (ro[i] + 2 * ro[i - 1])

    def c(i):
        return 3 * ro[i - 1]

    def d(i):
        return (ro[i] - ro[i - 1]) / h(i)

    ans = []
    for j in range(samples_n):
        i = 0
        while i < len(nodes) - 1 and nodes[i] < samples[j]:
            i += 1
        ans.append(y[i - 1] + b(i) * (samples[j - 1] - nodes[i - 1]) + c(i) * (samples[j - 1] - nodes[i - 1]) ** 2 +
                   d(i) * (samples[j - 1] - nodes[i - 1]) ** 3)
    return ans


def quadriSpline(nodes, mode):
    y = f(nodes)
    nodes_n = len(nodes)

    def h(i):
        return nodes[i] - nodes[i - 1]

    def delt(i):
        return (y[i] - y[i - 1]) / (nodes[i] - nodes[i - 1])

    B2 = []

    def b2(i):
        return B2[i - 1]

    def c2(i):
        return (b2(i + 1) - b2(i)) / (2 * h(i))

    L = [[0 for i in range(nodes_n)] for j in range(nodes_n)]
    R = [0 for _ in range(nodes_n)]
    if mode == 0:
        L[0][0] = 1
        L[0][1] = 0
    else:
        L[0][0] = 1
        L[0][1] = 0
    for i in range(1, nodes_n - 1):  # 1 brzegi2
        L[i][i - 1] = 1
        L[i][i] = 1

    L[nodes_n - 1][nodes_n - 2] = 1
    L[nodes_n - 1][nodes_n - 1] = 1
    if mode==0:
        R[0] = delt(1)
    else:
        R[0] = 0
    for i in range(1, nodes_n):
        R[i] = 2 * delt(i)
    B2 = np.linalg.solve(L, R)

    ans = []
    for j in range(samples_n):
        i = 0
        while i < len(nodes) - 1 and nodes[i] < samples[j]:
            i += 1
        ans.append(y[i - 1] + b2(i) * (samples[j - 1] - nodes[i - 1]) + c2(i) * (samples[j - 1] - nodes[i - 1]) ** 2)
    return ans


def showGraph(nodes_n, x, y, graph_counter):
    nodes = genEquidistant(nodes_n)
    if y == 0:
        if x == 0:
            axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, quadriSpline(nodes, 0), '-.')
        else:
            axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, quadriSpline(nodes, 1), '-.')
        name = "f. 2. rzedu"
    else:
        if x == 0:
            axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, cubicSpline(nodes, 0), '-.')
        else:
            axis[x][y].plot(nodes, f(nodes), 'o', samples, f(samples), samples, cubicSpline(nodes, 1), '-.')
        name = "f. 3. rzedu"

    axis[x][y].legend(['Wezly', 'Funkcja', name], loc='best')
    axis[x][y].set_xlabel("x")
    axis[x][y].set_ylabel("y")
    axis[x][y].set_title(
        "Wyk. " + str(graph_counter) + "., funkcja sklejana, " + str(x+1) + ". Warunek, dla " + str(nodes_n) +
        " wezlow", fontweight="bold")
    graph_counter += 1
    return graph_counter


def comp_abs(Y):
    ans = 0
    for i in range(0, samples_n):
        idx = min_x + (max_x - min_x) * i / samples_n
        ans = max(ans, abs(Y[i] - f(idx)))
    return ans

print("n ",
      "quadri equidist 0; ",
      # "quadri chebys 0; ",
      "quadri equidist 1; ",
      # "quadri chebys 1;   ",
      "cubic equidist 0;   ",
      # "cubic chebys 0;   ",
      "cubic equidist 1;   ",
      # "cubic chebys 1   "
      )
min1 = min3 = min5 = min7 = -1
# '''
for i in range(5, 81):
    d1 = comp_abs(quadriSpline(genEquidistant(i), 0))
    d3 = comp_abs(quadriSpline(genEquidistant(i), 1))
    d5 = comp_abs(cubicSpline(genEquidistant(i), 0))
    d7 = comp_abs(cubicSpline(genEquidistant(i), 1))

    print(i, end=": ")
    print(d1, end="; ")
    print(d3, end="; ")
    print(d5, end="; ")
    print(d7)

    # figure, axis = plot.subplots(2, 2)
    # figure.suptitle("Highest difference")
    # graph_counter = showGraph(10, 0, 0, graph_counter, 0)
    # # graph_counter = showGraph(10, 1, 0, graph_counter, 0)
    # graph_counter = showGraph(10, 0, 1, graph_counter, 0)
    # # graph_counter = showGraph(10, 1, 1, graph_counter, 0)
    # plot.show()
    if min1 > d1 or min1 == -1:
        min1 = d1
        i1 = i
    if min3 > d3 or min3 == -1:
        min3 = d3
        i3 = i
    if min5 > d5 or min5 == -1:
        min5 = d5
        i5 = i
    if min7 > d7 or min7 == -1:
        min7 = d7
        i7 = i

    if i == 81 - 1:
        print()
        print("Najmniejsze niedokładności otrzymano dla:")
        print(i1, ": ", min1, "; ", i3, ": ", min3, "; ", i5, ": ", min5, "; ", i7, ": ", min7, sep="")
        print("---------------------")
        figure, axis = plot.subplots(2, 2)
        figure.suptitle("Wykresy")
        graph_counter = showGraph(i1, 0, 0, graph_counter)
        graph_counter = showGraph(i5, 0, 1, graph_counter)
        graph_counter = showGraph(i3, 1, 0, graph_counter)
        graph_counter = showGraph(i7, 1, 1, graph_counter)
        plot.show()

'''
figure, axis = plot.subplots(2,2)
figure.suptitle("Wykresy")
graph_counter = showGraph(10, 0, 0, graph_counter)
graph_counter = showGraph(10, 1, 0, graph_counter)
graph_counter = showGraph(10, 0, 1, graph_counter)
graph_counter = showGraph(10, 1, 1, graph_counter)
plot.show()
figure, axis = plot.subplots(2,2)
figure.suptitle("Wykresy")
graph_counter = showGraph(25, 0, 0, graph_counter)
graph_counter = showGraph(25, 1, 0, graph_counter)
graph_counter = showGraph(25, 0, 1, graph_counter)
graph_counter = showGraph(25, 1, 1, graph_counter)
plot.show()
figure, axis = plot.subplots(2,2)
figure.suptitle("Wykresy")
graph_counter = showGraph(55, 0, 0, graph_counter)
graph_counter = showGraph(55, 1, 0, graph_counter)
graph_counter = showGraph(55, 0, 1, graph_counter)
graph_counter = showGraph(55, 1, 1, graph_counter)
plot.show()'''
