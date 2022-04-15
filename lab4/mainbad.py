from math import pi, cos, e, sin

import numpy as np
# from scipy.special import roots_chebyt

import matplotlib.pyplot as plot

from scipy.interpolate import interp1d
# from scipy.interpolate import BSpline, make_interp_spline

import numpy.linalg
from scipy.special import roots_chebyt

k = 4
m = 2

n = 10
min_x = -pi
max_x = 3 * pi


wykres_licznik=1

def f(x):
    if not isinstance(x, float):
        return [e ** (k * cos(m * i)) for i in x]
    return e ** (k * cos(m * x))


def df(x):
    if not isinstance(x, float):
        return [-k * m * (e ** (k * cos(m * i))) * sin(m * i) for i in x]
    return -k * m * (e ** (k * cos(m * x))) * sin(m * x)


def ddf(x):
    if not isinstance(x, float):
        return [16 * (e ** (4 * cos(2 * i))) * (4 * (sin(2 * i) ** 2) - cos(2 * i)) for i in x]
    return 16 * (e ** (4 * cos(2 * x))) * (4 * (sin(2 * x) ** 2) - cos(2 * x))


X = np.arange(min_x, max_x + 0.01, 0.01)
N = len(X)
# x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

x_pos, _ = roots_chebyt(n)
x_pos *= ((max_x - min_x) / 2)
x_pos += ((max_x + min_x) / 2)
y = f(x_pos)

# x = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))
x1 = np.linspace(min_x, max_x, num=30, endpoint=True)



for i in range(5,50,5):
    figure, axis = plot.subplots(2,2)
    figure.suptitle("Highest difference")
    n=i
    x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))

    # x_pos, _ = roots_chebyt(n)
    # x_pos *= ((max_x - min_x) / 2)
    # x_pos += ((max_x + min_x) / 2)
    y = f(x_pos)
    def h(i):
        return x_pos[i] - x_pos[i - 1]


    def delt(i):
        return (y[i] - y[i - 1]) / (x_pos[i] - x_pos[i - 1])


    def delt2(i):
        return (delt(i + 1) - delt(i)) / (x_pos[i + 1] - x_pos[i - 1])


    def delt3(i):
        return (delt2(i + 1) - delt2(i)) / (x_pos[i + 2] - x_pos[i - 1])


    A = [[0 for i in range(n)] for j in range(n)]
    # B = [0 for j in range(n)]
    P = [0 for _ in range(n)]
    A[0][0] = -h(1)
    A[0][1] = h(1)
    for i in range(1, n - 1): #1 brzegi3
        A[i][i - 1] = h(i)
        A[i][i] = 2 * (h(i) + h(i + 1))
        A[i][i + 1] = h(i + 1)

    A[n - 1][n - 2] = h(n - 1)
    A[n - 1][n - 1] = -h(n - 1)
    P[0] = h(1) ** 2 * delt3(1)
    for i in range(1, n - 1):
        P[i] = delt(i + 1) - delt(i)
    P[n - 1] = -(h(n - 1) ** 2) * delt3(n - 3)

    # for i in A:
    #     print(i)
    # print("P:")
    # print(P)
    ro = np.linalg.solve(A, P)
    # print(ro)



    def s(x, degree):
        # def ro(i1):
        #
        #
        #     # return ddf(x_pos[i]) / 6
        #     return roo[i1]

        def b(i):
            return (y[i] - y[i - 1]) / h(i) - h(i) * (ro[i] + 2 * ro[i - 1])

        def c(i):
            return 3 * ro[i - 1]

        def d(i):
            return (ro[i] - ro[i - 1]) / h(i)

        B2 = []

        def b2(i):
            # print(i)
            return B2[i-1]

        def c2(i):
            return (b2(i+1)-b2(i))/(2*h(i))

        L = [[0 for i in range(n)] for j in range(n)]
        # B = [0 for j in range(n)]
        R = [0 for _ in range(n)]
        L[0][0] = 1
        L[0][1] = 0
        for i in range(1, n - 1):#1 brzegi2
            L[i][i - 1] = 1
            L[i][i] = 1
            # L[i][i + 1] = h(i + 1)

        L[n - 1][n - 2] = 1
        L[n - 1][n - 1] = 1
        R[0] = 0
        for i in range(1, n):
            R[i] = 2*delt(i)
        # R[n - 1] = 2*delt(n-1)
        B2 = np.linalg.solve(L, R)
        # print(B2)



        # def c2(i):
        #     return (b2(i+1)-b2(i))/(2*h(i))

        ans = []
        # for k in range(len(x_pos)):
        for j in range(len(x)):
            # if x[j]
            i = 0
            while i < len(x_pos)-1 and x_pos[i] < x[j]:
                i += 1
            # i -= 1
            if degree == 3:
                ans.append(y[i - 1] + b(i) * (x[j - 1] - x_pos[i - 1]) + c(i) * (x[j - 1] - x_pos[i - 1]) ** 2 + d(i) * (
                        x[j - 1] - x_pos[i - 1]) ** 3)
            elif degree == 2:
                ans.append(y[i - 1] + b2(i) * (x[j - 1] - x_pos[i - 1]) + c2(i) * (x[j - 1] - x_pos[i - 1]) ** 2)
        return ans

    axis[0][0].plot(x_pos, y, 'o', X, f(X),
              # x_pos, f1(x_pos), '-', x_pos, f2(x_pos), '--', x_pos, f3(x_pos), '--',
              X, s(X, 2), '--',
              X, s(X, 3), '-.'
              )
    axis[0][0].legend(['Wezly', 'Funkcja',
                 # '1 stopnia', '2 stopnia', '3 stopnia',
                 '2rd?',
                 '3rd?'
                 ], loc='best')
    axis[0][0].set_xlabel("x")
    axis[0][0].set_ylabel("y")
    axis[0][0].set_title("Wyk. "+str(wykres_licznik)+"., funkcja sklejana " + str(3) + " rzedu dla " + str(n) +
                      " wezlow rownomiernych", fontweight="bold")
    wykres_licznik+=1
    # plot.show()

    x_pos, _ = roots_chebyt(n)
    x_pos *= ((max_x - min_x) / 2)
    x_pos += ((max_x + min_x) / 2)
    y = f(x_pos)

    A = [[0 for i in range(n)] for j in range(n)]
    # B = [0 for j in range(n)]
    P = [0 for _ in range(n)]
    A[0][0] = -h(1)
    A[0][1] = h(1)
    for i in range(1, n - 1):#1 brzegi3
        A[i][i - 1] = h(i)
        A[i][i] = 2 * (h(i) + h(i + 1))
        A[i][i + 1] = h(i + 1)

    A[n - 1][n - 2] = h(n - 1)
    A[n - 1][n - 1] = -h(n - 1)
    P[0] = h(1) ** 2 * delt3(1)
    for i in range(1, n - 1):
        P[i] = delt(i + 1) - delt(i)
    P[n - 1] = -(h(n - 1) ** 2) * delt3(n - 3)

    # for i in A:
    #     print(i)
    # print("P:")
    # print(P)
    ro = np.linalg.solve(A, P)
    # print(ro)



    def s(x, degree):
        # def ro(i1):
        #
        #
        #     # return ddf(x_pos[i]) / 6
        #     return roo[i1]

        def b(i):
            return (y[i] - y[i - 1]) / h(i) - h(i) * (ro[i] + 2 * ro[i - 1])

        def c(i):
            return 3 * ro[i - 1]

        def d(i):
            return (ro[i] - ro[i - 1]) / h(i)

        B2 = []

        def b2(i):
            # print(i)
            return B2[i-1]

        def c2(i):
            return (b2(i+1)-b2(i))/(2*h(i))

        L = [[0 for i in range(n)] for j in range(n)]
        # B = [0 for j in range(n)]
        R = [0 for _ in range(n)]
        L[0][0] = 1
        L[0][1] = 0
        for i in range(1, n - 1):#1 brzegi2
            L[i][i - 1] = 1
            L[i][i] = 1
            # L[i][i + 1] = h(i + 1)

        L[n - 1][n - 2] = 1
        L[n - 1][n - 1] = 1
        R[0] = 0
        for i in range(1, n):
            R[i] = 2*delt(i)
        # R[n - 1] = 2*delt(n-1)
        B2 = np.linalg.solve(L, R)
        # print(B2)



        # def c2(i):
        #     return (b2(i+1)-b2(i))/(2*h(i))

        ans = []
        # for k in range(len(x_pos)):
        for j in range(len(x)):
            # if x[j]
            i = 0
            while i < len(x_pos)-1 and x_pos[i] < x[j]:
                i += 1
            # i -= 1
            if degree == 3:
                ans.append(y[i - 1] + b(i) * (x[j - 1] - x_pos[i - 1]) + c(i) * (x[j - 1] - x_pos[i - 1]) ** 2 + d(i) * (
                        x[j - 1] - x_pos[i - 1]) ** 3)
            elif degree == 2:
                ans.append(y[i - 1] + b2(i) * (x[j - 1] - x_pos[i - 1]) + c2(i) * (x[j - 1] - x_pos[i - 1]) ** 2)
        return ans

    axis[1][0].plot(x_pos, y, 'o', X, f(X),
              # x_pos, f1(x_pos), '-', x_pos, f2(x_pos), '--', x_pos, f3(x_pos), '--',
              X, s(X, 2), '--',
              X, s(X, 3), '-.'
              )
    axis[1][0].legend(['Wezly', 'Funkcja',
                 # '1 stopnia', '2 stopnia', '3 stopnia',
                 '2rd?',
                 '3rd?'
                 ], loc='best')
    axis[1][0].set_xlabel("x")
    axis[1][0].set_ylabel("y")
    axis[1][0].set_title("Wyk. "+str(wykres_licznik)+"., funkcja sklejana " + str(3) + " rzedu dla " + str(n) +
                      " wezlow Chebysheva", fontweight="bold")
    wykres_licznik+=1
    
    x_pos = np.arange(min_x, max_x + 0.01, (max_x - min_x) / (n - 1))
    
    # x_pos, _ = roots_chebyt(n)
    # x_pos *= ((max_x - min_x) / 2)
    # x_pos += ((max_x + min_x) / 2)
    y = f(x_pos)
    def h(i):
        return x_pos[i] - x_pos[i - 1]
    
    
    def delt(i):
        return (y[i] - y[i - 1]) / (x_pos[i] - x_pos[i - 1])
    
    
    def delt2(i):
        return (delt(i + 1) - delt(i)) / (x_pos[i + 1] - x_pos[i - 1])
    
    
    def delt3(i):
        return (delt2(i + 1) - delt2(i)) / (x_pos[i + 2] - x_pos[i - 1])
    
    
    A = [[0 for i in range(n)] for j in range(n)]
    # B = [0 for j in range(n)]
    P = [0 for _ in range(n)]
    A[0][0] = -h(1)
    A[0][1] = h(1)
    for i in range(1, n - 1): #2 brzegi3
        A[i][i - 1] = h(i)
        A[i][i] = 2 * (h(i) + h(i + 1))
        A[i][i + 1] = h(i + 1)
    
    A[n - 1][n - 2] = h(n - 1)
    A[n - 1][n - 1] = -h(n - 1)
    P[0] = 0
    for i in range(1, n - 1):
        P[i] = delt(i + 1) - delt(i)
    P[n - 1] = 0
    
    # for i in A:
    #     print(i)
    # print("P:")
    # print(P)
    ro = np.linalg.solve(A, P)
    # print(ro)
    
    
    
    def s(x, degree):
        # def ro(i1):
        #
        #
        #     # return ddf(x_pos[i]) / 6
        #     return roo[i1]
    
        def b(i):
            return (y[i] - y[i - 1]) / h(i) - h(i) * (ro[i] + 2 * ro[i - 1])
    
        def c(i):
            return 3 * ro[i - 1]
    
        def d(i):
            return (ro[i] - ro[i - 1]) / h(i)
    
        B2 = []
    
        def b2(i):
            # print(i)
            return B2[i-1]
    
        def c2(i):
            return (b2(i+1)-b2(i))/(2*h(i))
    
        L = [[0 for i in range(n)] for j in range(n)]
        # B = [0 for j in range(n)]
        R = [0 for _ in range(n)]
        L[0][0] = 1
        L[0][1] = 0
        for i in range(1, n - 1):#2 brzegi2
            L[i][i - 1] = 1
            L[i][i] = 1
            # L[i][i + 1] = h(i + 1)
    
        L[n - 1][n - 2] = 1
        L[n - 1][n - 1] = 1
        R[0] = 0
        for i in range(1, n):
            R[i] = 2*delt(i)
        # R[n - 1] = 2*delt(n-1)
        B2 = np.linalg.solve(L, R)
        # print(B2)
    
    
    
        # def c2(i):
        #     return (b2(i+1)-b2(i))/(2*h(i))
    
        ans = []
        # for k in range(len(x_pos)):
        for j in range(len(x)):
            # if x[j]
            i = 0
            while i < len(x_pos)-1 and x_pos[i] < x[j]:
                i += 1
            # i -= 1
            if degree == 3:
                ans.append(y[i - 1] + b(i) * (x[j - 1] - x_pos[i - 1]) + c(i) * (x[j - 1] - x_pos[i - 1]) ** 2 + d(i) * (
                        x[j - 1] - x_pos[i - 1]) ** 3)
            elif degree == 2:
                ans.append(y[i - 1] + b2(i) * (x[j - 1] - x_pos[i - 1]) + c2(i) * (x[j - 1] - x_pos[i - 1]) ** 2)
        return ans
    
    axis[0][1].plot(x_pos, y, 'o', X, f(X),
                 # x_pos, f1(x_pos), '-', x_pos, f2(x_pos), '--', x_pos, f3(x_pos), '--',
                 X, s(X, 2), '--',
                 X, s(X, 3), '-.'
                 )
    axis[0][1].legend(['Wezly', 'Funkcja',
                    # '1 stopnia', '2 stopnia', '3 stopnia',
                    '2rd?',
                    '3rd?'
                    ], loc='best')
    axis[0][1].set_xlabel("x")
    axis[0][1].set_ylabel("y")
    axis[0][1].set_title("Wyk. "+str(wykres_licznik)+"., funkcja sklejana " + str(3) + " rzedu dla " + str(n) +
                      " wezlow rownomiernych", fontweight="bold")
    wykres_licznik+=1
    # plot.show()
    
    x_pos, _ = roots_chebyt(n)
    x_pos *= ((max_x - min_x) / 2)
    x_pos += ((max_x + min_x) / 2)
    y = f(x_pos)
    
    A = [[0 for i in range(n)] for j in range(n)]
    # B = [0 for j in range(n)]
    P = [0 for _ in range(n)]
    A[0][0] = -h(1)
    A[0][1] = h(1)
    for i in range(1, n - 1):#2 brzegi3
        A[i][i - 1] = h(i)
        A[i][i] = 2 * (h(i) + h(i + 1))
        A[i][i + 1] = h(i + 1)
    
    A[n - 1][n - 2] = h(n - 1)
    A[n - 1][n - 1] = -h(n - 1)
    P[0] = 0
    for i in range(1, n - 1):
        P[i] = delt(i + 1) - delt(i)
    P[n - 1] = 0
    
    # for i in A:
    #     print(i)
    # print("P:")
    # print(P)
    ro = np.linalg.solve(A, P)
    # print(ro)
    
    
    
    def s(x, degree):
        # def ro(i1):
        #
        #
        #     # return ddf(x_pos[i]) / 6
        #     return roo[i1]
    
        def b(i):
            return (y[i] - y[i - 1]) / h(i) - h(i) * (ro[i] + 2 * ro[i - 1])
    
        def c(i):
            return 3 * ro[i - 1]
    
        def d(i):
            return (ro[i] - ro[i - 1]) / h(i)
    
        B2 = []
    
        def b2(i):
            # print(i)
            return B2[i-1]
    
        def c2(i):
            return (b2(i+1)-b2(i))/(2*h(i))
    
        L = [[0 for i in range(n)] for j in range(n)]
        # B = [0 for j in range(n)]
        R = [0 for _ in range(n)]
        L[0][0] = 1
        L[0][1] = 0
        for i in range(1, n - 1):#2 brzegi2
            L[i][i - 1] = 1
            L[i][i] = 1
            # L[i][i + 1] = h(i + 1)
    
        L[n - 1][n - 2] = 1
        L[n - 1][n - 1] = 1
        R[0] = 0
        for i in range(1, n):
            R[i] = 2*delt(i)
        # R[n - 1] = 2*delt(n-1)
        B2 = np.linalg.solve(L, R)
        # print(B2)
    
    
    
        # def c2(i):
        #     return (b2(i+1)-b2(i))/(2*h(i))
    
        ans = []
        # for k in range(len(x_pos)):
        for j in range(len(x)):
            # if x[j]
            i = 0
            while i < len(x_pos)-1 and x_pos[i] < x[j]:
                i += 1
            # i -= 1
            if degree == 3:
                ans.append(y[i - 1] + b(i) * (x[j - 1] - x_pos[i - 1]) + c(i) * (x[j - 1] - x_pos[i - 1]) ** 2 + d(i) * (
                        x[j - 1] - x_pos[i - 1]) ** 3)
            elif degree == 2:
                ans.append(y[i - 1] + b2(i) * (x[j - 1] - x_pos[i - 1]) + c2(i) * (x[j - 1] - x_pos[i - 1]) ** 2)
        return ans
    
    axis[1][1].plot(x_pos, y, 'o', X, f(X),
                 # x_pos, f1(x_pos), '-', x_pos, f2(x_pos), '--', x_pos, f3(x_pos), '--',
                 X, s(X, 2), '--',
                 X, s(X, 3), '-.'
                 )
    axis[1][1].legend(['Wezly', 'Funkcja',
                    # '1 stopnia', '2 stopnia', '3 stopnia',
                    '2rd?',
                    '3rd?'
                    ], loc='best')
    axis[1][1].set_xlabel("x")
    axis[1][1].set_ylabel("y")
    axis[1][1].set_title("Wyk. "+str(wykres_licznik)+"., funkcja sklejana " + str(3) + " rzedu dla " + str(n) +
                      " wezlow Chebysheva", fontweight="bold")
    wykres_licznik+=1
    
    plot.show()
'''
    def comp_abs(Y):
        ans = 0
        for i in range(0, N):
            idx = min_x + (max_x - min_x) * i / N
            ans = max(ans, abs(Y[i] - f(idx)))
        return ans


    def comp_sqr(Y):
        ans = 0
        for i in range(0, N):
            idx = min_x + (max_x - min_x) * i / N
            ans += (Y[i] - f(idx)) ** 2
        return ans

    d1 = comp_abs(s(X,2))
    d2 = comp_abs(drawHermitCheby(0))

    print(n, end=": ")
    print(d1, end="; ")
    print(d2)

    if min1 > d1:
        min1 = d1
        i1 = n
    if min2 > d2:
        min2 = d2
        i2 = n

print()
print("Najmniejsze niedokładności otrzymano dla:")
print(i1, ": ", min1, "; ", i2, ": ", min2, sep="")
print("---------------------")'''

