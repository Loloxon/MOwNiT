from typing import Callable, List

import mpmath as mp
import numpy as np
from mpmath import *


# zad 1 – funkcja d) n=15, m=10, [-1.5, 0.3];
# zad 2 – d)

# initial data
a = -1.5
b = 0.3
def f(x):
    n=15
    m=10
    return x**n+x**m

def df(x):
    n=15
    m=10
    return n*x**(n-1)+m*x**(m-1)

epsilon = 10e-24
max_iteration = 100000
mp.dps = 8
mp.prec = 8

def newton_method(x0, epsilon, mode):
    x1 = x0 - f(x0)/df(x0)
    for i in range(1,max_iteration):
        if mode==1:
            if abs(x0-x1)<epsilon:
                return "moduł różnicy: {}".format(abs(x0-x1)), i, x1
        else:
            if abs(f(x1))<epsilon:
                return "moduł funkcji: {}".format(abs(f(x1))), i, x1
        x0, x1 = x1, x1 - f(x1) / df(x1)
    return -1, -1, x1

def secant_method(a, b, epsilon, mode):
    x0, x1 = a, b
    # if f(x0) * f(x1) > 0:
    #     return "-", "-", "-"
    for i in range(0,max_iteration):
        if mode==1:
            if abs(x0-x1)<epsilon:
                return "moduł różnicy: {}".format(abs(x0-x1)), i, x1
        else:
            if abs(f(x1))<epsilon:
                return "moduł funkcji: {}".format(abs(f(x1))), i, x1
        if f(x1)-f(x0)==0:
            return "-", "-", "-"
        x0, x1 = x1, x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        # if(f(x1)*f(tmp)<=0):
        #     x0, x1 = x, x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        # else:
        #     x0, x1 = x1, x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
    return -1, -1, x1

def test_newton(epsilon, mode):
    print("Od lewego krańca")
    for i in range(int(a*10), int(b*10)+1):
        print(i/10, end=" & ")
        if i!=0:
            for epsilon in [0.01, 0.0001]:
                diff, iters1, x1 = newton_method(i/10, epsilon, mode)
                if epsilon != 0.0001:
                    print(iters1,",",x1, end=" & ")
                else:
                    print(iters1,",",x1, end=" \\\\ \\hline\n")
        else:
            print("-", ",", "-", end=" & ")
            print("-", ",", "-", end=" \\\\ \\hline\n")
    for i in range(int(a*10), int(b*10)+1):
        print(i/10, end=" & ")
        if i!=0:
            for epsilon in [0.000001, 0.00000001]:
                diff, iters1, x1 = newton_method(i/10, epsilon, mode)
                if epsilon != 0.00000001:
                    print(iters1,",",x1, end=" & ")
                else:
                    print(iters1,",",x1, end=" \\\\ \\hline\n")
        else:
            print("-", ",", "-", end=" & ")
            print("-", ",", "-", end=" \\\\ \\hline\n")

    print("Od prawego krańca")
    for i in range(int(b*10), int(a*10)-1, -1):
        print(i/10, end=" & ")
        if i!=0:
            for epsilon in [0.01, 0.0001]:
                diff, iters1, x1 = newton_method(i/10, epsilon, mode)
                if epsilon != 0.0001:
                    print(iters1,",",x1, end=" & ")
                else:
                    print(iters1,",",x1, end=" \\\\ \\hline\n")
        else:
            print("-", ",", "-", end=" & ")
            print("-", ",", "-", end=" \\\\ \\hline\n")
    for i in range(int(b*10), int(a*10)-1, -1):
        print(i/10, end=" & ")
        if i!=0:
            for epsilon in [0.000001, 0.00000001]:
                diff, iters1, x1 = newton_method(i/10, epsilon, mode)
                if epsilon != 0.00000001:
                    print(iters1,",",x1, end=" & ")
                else:
                    print(iters1,",",x1, end=" \\\\ \\hline\n")
        else:
            print("-", ",", "-", end=" & ")
            print("-", ",", "-", end=" \\\\ \\hline\n")

def test_secant(epsilon, mode):
    print("Od lewego krańca")
    for i in range(int(b*10), int(a*10), -1):
        # print(a, i/10, end='''         ''')
        # print(secant_method(a, i/10, epsilon, mode))

        print([a, i/10], end=" & ")
        for epsilon in [0.01, 0.0001]:
            diff, iters1, x1 = secant_method(a, i/10, epsilon, mode)
            if epsilon != 0.0001:
                print(iters1,",",x1, end=" & ")
            else:
                print(iters1,",",x1, end=" \\\\ \\hline\n")
    for i in range(int(b*10), int(a*10), -1):
        # print(a, i/10, end='''         ''')
        # print(secant_method(a, i/10, epsilon, mode))

        print([a, i/10], end=" & ")
        for epsilon in [0.000001, 0.00000001]:
            diff, iters1, x1 = secant_method(a, i/10, epsilon, mode)
            if epsilon != 0.00000001:
                print(iters1,",",x1, end=" & ")
            else:
                print(iters1,",",x1, end=" \\\\ \\hline\n")
    print("Do prawego krańca")
        # print(i/10, b, end='''         ''')
        # print(secant_method(i/10, b, epsilon, mode))

        # print([i/10, b], end=" & ")
        # for epsilon in [0.01, 0.0001, 0.000001, 0.00000001]:
        #     diff, iters1, x1 = secant_method(i/10, b, epsilon, mode)
        #     if epsilon != 0.00000001:
        #         print(iters1,",",x1, end=" & ")
        #     else:
        #         print(iters1,",",x1, end=" \\\\ \\hline\n")

    for i in range(int(a*10), int(b*10)):
        # print(a, i/10, end='''         ''')
        # print(secant_method(a, i/10, epsilon, mode))

        print([i/10, b], end=" & ")
        for epsilon in [0.01, 0.0001]:
            diff, iters1, x1 = secant_method(i/10, b, epsilon, mode)
            if epsilon != 0.0001:
                print(iters1,",",x1, end=" & ")
            else:
                print(iters1,",",x1, end=" \\\\ \\hline\n")
    for i in range(int(a*10), int(b*10)):
        # print(a, i/10, end='''         ''')
        # print(secant_method(a, i/10, epsilon, mode))

        print([i/10, b], end=" & ")
        if f(i/10)!=f(b):
            for epsilon in [0.000001, 0.00000001]:
                diff, iters1, x1 = secant_method(i/10, b, epsilon, mode)
                if epsilon != 0.00000001:
                    print(iters1,",",x1, end=" & ")
                else:
                    print(iters1,",",x1, end=" \\\\ \\hline\n")
        else:
            print("-", ",", "-", end=" & ")
            print("-", ",", "-", end=" \\\\ \\hline\n")
epsilon = 0.00000001
test_newton(epsilon, 0)
print("=====================================")
# test_secant(epsilon, 0)

'''
def F(X):
    ret = [0,0,0]
    ret[0] = X[0]**2 + X[1]**2 + X[2] - 1
    ret[1] = 2*X[0]**2 + X[1]**2 + X[2]**3 - 2
    ret[2] = 3*X[0] - 2*X[1]**3 - 2*X[2]**2 - 3
    return ret

def J(X):
    ret = [[2*X[0],     2*X[1],     1],
           [4*X[0],     2*X[1],     3*X[2]**2],
           [3,          -6*X[1]**2, -4*X[2]]]
    return ret

def newton_matrix(F: Callable[[List[float]], List[float]], J: Callable[[List[float]], List[float]], X: List[float], epsilon: float, exit: int):
    X = np.array(X)
    iters = 0
    while True:
        A = np.copy(X)
        try:
            S = np.linalg.solve(J(X), F(X))
        except np.linalg.LinAlgError:
            return ['x', 'x', 'x'], 'x'
        X = X - S
        iters += 1
        if exit == 1:
            if np.linalg.norm(X-A) < epsilon: return X, iters
        elif exit == 2:
            if np.linalg.norm(F(X)) < epsilon: return X, iters
        if iters > 10000:
            return ['x','x','x'], 'x'

for end1 in [-1, -0.3, 0.3, 1]:
    for end2 in [-1, -0.3, 0.3, 1]:
        for end3 in [-1, -0.3, 0.3, 1]:
            print(f"[{end1:.1f}, {end2:.1f}, {end3:.1f}]", end=" & ")
            for prec in [0.000001]:
                x1, iters2 = newton_matrix(F, J, [end1, end2, end3], prec, 1)
                x2, iters1 = newton_matrix(F, J, [end1, end2, end3], prec, 2)
                if prec != 0.000001:
                    print(iters2,",",x1, end=" & ")
                else:
                    print(iters2,",",x1, end=" \\\\ \\hline\n")
'''