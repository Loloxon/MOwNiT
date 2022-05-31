import numpy
import numpy as np
import time
from decimal import Decimal
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plot
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize


def create_X(n):
    return np.array([float(-1) ** i for i in range(n)])


def create_matrix_A(n, k=11, m=3):
    A = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = k
            else:
                A[i, j] = m / (n - i - j - 2 + 0.5)
    return A


def create_vector_B(A, n):
    X = create_X(n)
    B = np.array([float(0) for _ in range(n)])
    for i in range(n):
        for j in range(n):
            B[i] += A[i, j] * X[j]
    return B


def f1(x, x_new):
    return max([abs(x[i] - x_new[i]) for i in range(len(x))])


def f2(x, A, B):
    return max([abs(np.dot(A, x[i]) - B[i]) for i in range(len(x))])


def jacoby(A, x, B, criterion, eps):
    n = len(x)
    L = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    D = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    U = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = A[i, j]
            if i == j:
                D[i, j] = A[i, j]
            if i < j:
                U[i, j] = A[i, j]

    N = numpy.linalg.inv(D)
    M = np.dot(-N, (L + U))

    iterations = 1
    x_new = np.dot(M, x) + np.dot(N, B)
    if criterion == 0:
        while f1(x, x_new) > eps:
            x, x_new = x_new, np.dot(M, x) + np.dot(N, B)
    elif criterion == 1:
        while f2(x, A, B) > eps:
            x, x_new = x_new, np.dot(M, x) + np.dot(N, B)
    return iterations, x
    # print(L, "L")
    # print()
    # print(D, "D")
    # print()
    # print(U, "U")
    # print()
    # print(N, "D^-1 (N)")
    # print()
    # print(L+U, "L+U")
    # print()
    # print(M, "-N*(L+U)")
    # print()
    # print(np.dot(M,x)+np.dot(N,B))
    # for i in range(15):
    #     x = np.dot(M, x) + np.dot(N, B)
    #     print(np.dot(M, x) + np.dot(N, B))
    # print()
