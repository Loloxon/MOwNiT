import numpy as np
import matplotlib.pyplot as plt
import sys
k = 0.5
m = 2.0
x0 = np.pi/4
xn = 3*np.pi


def fun(x, y):
    return k**2 * m * np.sin(m*x) * np.cos(m*x) + k * m * y * np.sin(m*x)


def rfun(x):
    return np.e**(-k*np.cos(m*x)) - k*np.cos(m*x) + 1


def euler_method(x0, y0, n, xn, f):
    h = (xn - x0) / (n - 1)

    x = np.linspace(x0, xn, n)
    y = np.zeros([n])
    y[0] = y0

    for i in range(1, n):
        y[i] = h * f(x[i-1], y[i-1]) + y[i-1]

    return x, y


def max_error(x_in, x_calc):
    for i in range(len(x_in)):
        x_in[i] = abs(x_in[i] - x_calc[i])
    return max(np.array(x_in))


amount = 1000
numbers = np.linspace(x0, xn, 1000)
sys.argv[2] = sys.argv[2]+sys.argv[1]+".pdf"
n = int(sys.argv[1])
points, calculated = euler_method(x0, rfun(x0), n, xn, fun)
print(max_error(list(map(rfun, points)), calculated))

plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Euler')
plt.plot(numbers, list(map(rfun, numbers)), 'b-', points, calculated, 'r-')
plt.savefig(sys.argv[2])


