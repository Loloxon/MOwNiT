import numpy as np
import math
import matplotlib.pyplot as plt
import sys


def coef(nodes, values):
    n = len(nodes)
    a = []
    for i in range(n):
        a.append(values[i])

    for j in range(1, n):

        for i in range(n-1, j-1, -1):
            a[i] = float(a[i]-a[i-1])/float(nodes[i] - nodes[i - j])

    return np.array(a)


def newton(f, x, nodes):
    a = coef(nodes, list(map(f, nodes)))
    n = len(a) - 1
    temp = a[n]
    for i in range(n - 1, -1, -1):
        temp = temp * (x - nodes[i]) + a[i]
    return temp


def fun(x):
    return math.sin(x) * math.sin(x**2/math.pi)


def chebyshew(x0, x1, n):
    result = []
    for i in range(1, n+1, 1):
        result.append(1/2*(x0 + x1) + 1/2*(x1 - x0)*math.cos((2*i - 1)*math.pi/(2*n)))
    return result


def equidistant(x0, x1, n):
    result = []
    for i in range(n):
        result.append(x0 + i*(x1-x0)/(n-1))
    return result


amount = 1000
x0 = -math.pi
x1 = 2 * math.pi
numbers = range(amount)
points = list(map(lambda x: (x0 + x*(x1-x0)/amount), numbers))
values = list(map(fun, points))
interpolated = []
n = int(sys.argv[1])
sys.argv[2] = sys.argv[2] + sys.argv[1] + ".pdf"
for i in range(amount):
    interpolated.append(newton(fun, points[i], chebyshew(x0, x1, n)))
'''
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Interpolation')
plt.plot(points, values, 'b-', points, interpolated, 'r-',
         equidistant(x0, x1, n), list(map(fun, equidistant(x0, x1, n))), 'y.')
plt.savefig(sys.argv[2])
'''
diff = []
for i in range(amount):
    diff.append(interpolated[i] - values[i])
print(n, end=' ')
print(np.linalg.norm(diff))