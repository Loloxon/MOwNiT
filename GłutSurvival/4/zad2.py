import numpy as np
import math
import matplotlib.pyplot as plt
import sys


def fun(x):
    return math.sin(x) * math.sin(x**2/math.pi)


def dfun(x):
    return 2*x*math.sin(x)*math.cos(x**2/math.pi)/math.pi + math.sin(x**2/math.pi)*math.cos(x)


def hermie(nodes, f, df, x):
    n = len(nodes)
    z = []
    for i in range(n):
        z.append(nodes[i])
        z.append(nodes[i])
    n2 = 2*n
    matrix = np.zeros((n2, n2))
    for i in range(n2):
        for j in range(i+1):
            if j == 0:
                matrix[i][j] = f(z[i])
            elif j == 1 & i % 2 == 1:
                matrix[i][j] = df(z[i])
            else:
                matrix[i][j] = matrix[i][j-1] - matrix[i-1][j-1]
                matrix[i][j] = matrix[i][j] / (z[i] - z[i-j])

    result = 0
    helper = 1
    for i in range(n2):
        result = result + matrix[i][i] * helper
        helper = helper * (x - z[i])
    return result


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
x1 = 2*math.pi
numbers = range(amount)
points = list(map(lambda x: (x0 + x*(x1-x0)/amount), numbers))
values = list(map(fun, points))
interpolated = []
n = int(sys.argv[1])
sys.argv[2] = sys.argv[2] + sys.argv[1] + ".pdf"
for i in range(amount):
    interpolated.append(hermie(chebyshew(x0, x1, n), fun, dfun, points[i]))
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
