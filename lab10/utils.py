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
                A[i, j] = m/(n-i-j-2+0.5)
    return A

def create_vector_B(A, n):
    X = create_X(n)
    B = np.array([float(0) for _ in range(n)])
    for i in range(n):
        for j in range(n):
            B[i] += A[i, j] * X[j]
    return B

def jacoby(A, x):
    n = len(x)
    L = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = A[i, j]
    D = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i == j:
                D[i, j] = A[i, j]
    U = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i < j:
                U[i, j] = A[i, j]
    # N = D.transpose()
    N = numpy.linalg.inv(D)
    M = -N*(L+U)
    print(L, "l")
    print(D, "d")
    print(U, "u")
    print(N, "n")
    print(L+U, "n")
    print(M, "m")


