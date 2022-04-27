import numpy as np
import math
import matplotlib.pyplot as plt
import sys


def fun(x):
    return math.sin(x) * math.sin(x**2/math.pi)


def lagrange(f, x, n, nodes):
    values = list(map(f, nodes))
    sum = 0
    for i in range(n):
        term = 1
        for j in range(n):
            if i != j:
                term = term * (x - nodes[j]) / (nodes[i] - nodes[j])
        term = term*values[i]
        sum += term
    return sum


def chebyshew(x0, x1, n):
    result = []
    for i in range(1, n+1, 1):
        result.append(1/2*(x0 + x1) + 1/2*(x1 - x0)*math.cos((2*i - 1)*math.pi/(2*n)))
    return result


def equadistant(x0, x1, n):
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
    interpolated.append(lagrange(fun, points[i], n, equadistant(x0, x1, n)))
'''
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Interpolation')
plt.plot(points, values, 'b-', points, interpolated, 'r-',
         chebyshew(x0, x1, n), list(map(fun, chebyshew(x0, x1, n))), 'y.')
plt.savefig(sys.argv[2])
'''
diff = []
for i in range(amount):
    diff.append(interpolated[i] - values[i])
print(n, end=' ')
print(np.linalg.norm(diff))
