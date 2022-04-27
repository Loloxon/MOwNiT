import numpy as np
import matplotlib.pyplot as plt
import sys
k = 1.0
m = 2.0
x0 = 0
xn = (2*np.pi + k)/m


def u(x):
    return m**2*k*x


def w(x):
    return 0


def v(x):
    return -m**2


def rfun(x):
    return -k*np.sin(m*x) + k*x


def max_error(x_in, x_calc):
    for i in range(len(x_in)):
        x_in[i] = abs(x_in[i] - x_calc[i])
    return max(np.array(x_in))


ddf_eq = {'u': u, 'v': v, 'w': w}


# ddf_eq = {'u': u(x), 'v': v(x), 'w': w(x)} \ - lambda
def finite_difference_method(x0, y0, xn, yn, n, ddf_eq):
    h = (xn - x0) / n

    x = np.linspace(x0, xn, n)
    a, b, d, c = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])
    for i in range(1, n-1):
        a[i] = -(1 + ddf_eq['w'](x[i]) * h / 2)
        b[i] = -ddf_eq['u'](x[i]) * h**2
        c[i] = -(1 - ddf_eq['w'](x[i]) * h / 2)
        d[i] = (2 + ddf_eq['v'](x[i]) * h**2)

    A = np.zeros([n, n])
    for i in range(2, n-2):
        A[i][i] = d[i]
        A[i][i-1] = a[i]
        A[i][i+1] = c[i]
    A[1][1], A[n-2][n-2] = d[1], d[n-2]
    A[1][2] = c[1]
    A[n-2][n-3] = a[n-2]

    b[1] = b[1] - a[1] * y0
    b[n-2] = b[n-2] - c[n-2] * yn
    b[0], b[n-1] = y0, yn
    A[0][0] = 1
    A[n-1][n-1] = 1

    y = np.linalg.solve(A, b)   # Ay = b

    return x, y


amount = 1000
numbers = np.linspace(x0, xn, 1000)
sys.argv[2] = sys.argv[2]+sys.argv[1]+".pdf"
n = int(sys.argv[1])
points, calculated = finite_difference_method(x0, 0, xn, rfun(xn), n, ddf_eq)
print(max_error(list(map(rfun, points)), calculated))
print(rfun(xn))
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('MRS')
plt.plot(numbers, list(map(rfun, numbers)), 'b-')
#plt.show()
plt.savefig(sys.argv[2])

