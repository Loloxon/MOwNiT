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


def get_trig_coeff(m, k, nodes):
    akl = []
    bkl = []
    n = len(nodes)
    for i in range(n):
        akl.append(fun(nodes[i]*3/2) * np.cos(k*nodes[i]))
        bkl.append(fun(nodes[i]*3/2) * np.sin(k*nodes[i]))
    return (sum(akl) / (n/2)), (sum(bkl) / (n/2))


def approximate_trig(x, nodes, m):  # N is for degree of the approximation
    a0 = get_trig_coeff(m, 0, nodes)[0]
    res = 0
    for k in range(1, m):
        ak, bk = get_trig_coeff(m, k, nodes)
        res += ak * np.cos(2/3 * k * x)
        res += bk * np.sin(2/3 * k * x)
    res += a0 / 2
    return res


n = int(sys.argv[1])
degree = int(sys.argv[2])
if n <= degree:
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
    regressed.append(approximate_trig(points[i], nodes, degree))

plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Regression')
plt.plot(points, values, 'b-', points, regressed, 'r-',
         nodes, list(map(fun, nodes)), 'y.')
plt.savefig(sys.argv[3])

norm = get_error(regressed, values)
z = "%i %i %f" % (n, degree, norm)
print(z)
