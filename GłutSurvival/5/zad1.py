import sys
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame


def save(filename, results):
    filename += '.xlsx'
    df = DataFrame(data=results)
    print(df)
    df.to_excel(filename, sheet_name='sheet1', index=False, header=False)


def f(x):
    return np.sin(x) * np.sin(x ** 2 / np.pi)


def get_ys(xs):
    return [f(x) for x in xs]


def spline3(x_points, y_points, xs, boundary_cond):
    size = len(x_points) - 2
    matrix = np.zeros((size, size))

    for i in range(size):
        for j in range(size):
            if i == j:
                matrix[i][j] = 4
            if j == i+1 or j == i-1:
                matrix[i][j] = 1

    g = np.zeros(size); h = []
    for i in range(size):
        h.append(x_points[i+1] - x_points[i])
        g[i] = 6 / (h[i]**2) * (y_points[i] - 2*y_points[i+1] + y_points[i+2])
    h.append(x_points[-1] - x_points[-2])
    z = np.linalg.solve(matrix, g)

    # the boundary condition is set here
    z = list(z)
    if boundary_cond == 1:
        z = [0] + z + [0]
    else:
        z = [z[0]] + z + [z[-1]]

    a = []; b = []; c = []; d = []
    for i in range(size+1):
        a.append((z[i+1] - z[i]) / (6 * h[i]))
        b.append(0.5 * z[i])
        c.append((y_points[i+1] - y_points[i]) / h[i] - (z[i+1] + 2 * z[i]) / 6 * h[i])
        d.append(y_points[i])

    nr_fun = 0
    ys = []
    for i in range(len(xs)):
        while x_points[nr_fun + 1] < xs[i] < x_points[-1]:
            nr_fun += 1
        ys.append(get_val([d[nr_fun], c[nr_fun], b[nr_fun], a[nr_fun]], x_points[nr_fun], xs[i]))

    return ys


def spline2(x_points, y_points, xs, boundary_cond):
    size = len(x_points) - 1
    matrix = np.zeros((size, size))

    for i in range(size):
        for j in range(size):
            if i == j:
                matrix[i][j] = 1
            if j == i-1:
                matrix[i][j] = 1

    g = np.zeros(size); h = []
    for i in range(size):
        h.append(x_points[i+1] - x_points[i])
        g[i] = 2 / h[i] * (y_points[i+1] - y_points[i])
    b = np.linalg.solve(matrix, g)

    # the boundary condition is set here
    b = list(b)
    #if boundary_cond == 1:
    b = [0] + b
    #else:
    #   b = [b[-1]] + b
        # b = [b[0]] + b

    a = []; c = []
    for i in range(size):
        a.append((b[i+1] - b[i]) / (2 * h[i]))
        c.append(y_points[i])
    if boundary_cond == 2:
        b[0] = (y_points[1] - y_points[0]) / (x_points[1] - x_points[0])
        a[0] = 0

    nr_fun = 0
    ys = []
    for i in range(len(xs)):
        while x_points[nr_fun + 1] < xs[i] < x_points[-1]:
            nr_fun += 1
        ys.append(get_val([c[nr_fun], b[nr_fun], a[nr_fun]], x_points[nr_fun], xs[i]))

    return ys


def get_val(coeff, xi, x):
    val = 0
    for i, elem in enumerate(coeff):
        val += elem * (x - xi) ** i
    return val


def get_norm(y1, y2, which):
    if which == "max":
        return np.linalg.norm(np.subtract(y2, y1), np.inf)
    if which == "eu":
        return np.linalg.norm(np.subtract(y2, y1))


def chebyshew(a, b, no):
    ret = sorted(np.array([np.cos(((2 * j - 1) / (2 * no)) * np.pi) for j in range(1, no + 1)], dtype=np.float64))
    for i in range(len(ret)):
        ret[i] = 0.5 * (a + b) + 0.5 * (b - a) * ret[i]
    return ret


def get_xs(a, b, n):
    step = (b-a)/(n-1)
    ret = []
    for i in range(n):
        ret.append(a)
        a += step
    return ret


start = -np.pi
end = 2*np.pi
n_points = 8
n_draw = 1000
spline = 2

xs = get_xs(start, end, n_draw)
ys_or = get_ys(xs)
res = [['Liczba węzłów', 'spl2, nat', 'spl2, lin', 'spl3, nat', 'spl3, par']]

for i in range(3, 21):
    print(i)
    plt.figure(figsize=(9, 6))
        #xp = get_xs(start, end, i)
    xp = get_xs(start, end, i)
    yp = get_ys(xp)
    plt.plot(xp, yp, 'y.', markersize=10)
    plt.plot(xs, ys_or, 'y', label='Interpolated')
    r = [i]
    #if True: #spline == 2 or spline == 0:
    ys = spline2(xp, yp, xs, 1)
    plt.plot(xs, ys, 'g', label='2nd degree natural spline')
    r.append(get_norm(ys, ys_or, 'eu'))
    ys = spline2(xp, yp, xs, 2)
    plt.plot(xs, ys, 'm', label='2nd degree first function is linear')
    r.append(get_norm(ys, ys_or, 'eu'))
    #if True: #spline == 3 or spline == 0:
    ys = spline3(xp, yp, xs, 1)
    plt.plot(xs, ys, 'grey', label='3rd degree natural spline')
    r.append(get_norm(ys, ys_or, 'eu'))
    ys = spline3(xp, yp, xs, 2)
    plt.plot(xs, ys, 'b', label='3rd degree parabolic spline')
    r.append(get_norm(ys, ys_or, 'eu'))
    res.append(r)
    s = "plot" + str(i) + ".pdf"
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.savefig(s)


save('res_eq', res)

