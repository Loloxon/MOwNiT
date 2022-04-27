import numpy as np
import math
import matplotlib.pyplot as plt
import sys


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


def get_error(x, y):
    err = 0
    for i in range(len(x)):
        err +=(x[i] - y[i])**2
    err = err / len(x)
    return np.sqrt(err)


def fun(x):
    return math.sin(x) * math.sin(x**2/math.pi)


def poly(vec, x):
    result = 0
    for i in range(len(vec)):
        result += vec[i]*x**i
    return result


def regress(f, nodes, degree):
    values = np.array(list(map(f, nodes)))
    n = len(nodes)
    matrix = np.zeros((n, degree+1))
    for i in range(n):
        for j in range(degree+1):
            matrix[i][j] += nodes[i]**j
    return np.linalg.inv(np.transpose(matrix).dot(matrix)).dot(np.transpose(matrix)).dot(values)


n = int(sys.argv[1])
degree = int(sys.argv[2])
if n<=degree:
    sys.exit()
sys.argv[3] = sys.argv[3]+sys.argv[1]+"_"+sys.argv[2]+".pdf"
amount = 1000
x0 = -math.pi
x1 = 2 * math.pi
numbers = range(amount)
points = list(map(lambda x: (x0 + x*(x1-x0)/amount), numbers))
values = list(map(fun, points))
nodes = equadistant(-np.pi, 2*np.pi, n)
regressed = []

for i in range(amount):
    regressed.append(poly(regress(fun, nodes, degree), points[i]))

plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Regression')
plt.plot(points, values, 'b-', points, regressed, 'r-',
         nodes, list(map(fun, nodes)), 'y.')
plt.savefig(sys.argv[3])

norm = get_error(regressed, values)
z = "%i %i %f" % (n, degree, norm)
print(z)
