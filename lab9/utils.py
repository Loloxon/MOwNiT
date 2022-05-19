import numpy as np
from decimal import Decimal


def gauss_decimal(A, b, n):
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            L = A[i, k] / A[k, k]
            for j in range(k + 1, n):
                A[i, j] = A[i, j] - L * A[k, j]
            b[i] = b[i] - L * b[k]

    x = b.copy()
    for i in range(n - 1, -1, -1):
        S = Decimal(0.0)
        for j in range(i + 1, n):
            S = S + A[i, j] * x[j]
        x[i] = Decimal((x[i] - S) / A[i, i])
    return x


def create_X_decimal(n):
    return np.array([Decimal(-1) ** i for i in range(n)])


def create_matrix_A_decimal(n):
    A = np.array([[Decimal(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i == 0:
                A[i, j] = Decimal(1)
            else:
                A[i, j] = Decimal(1 / (i + j + 1))
    return A


def create_vector_B_decimal(A, n):
    X = create_X_decimal(n)
    B = np.array([Decimal(0) for _ in range(n)])
    for i in range(n):
        for j in range(n):
            B[i] += A[i, j] * X[j]
    return B


def find_diff_decimal(newX, n):
    oldX = create_X_decimal(n)
    diff = 0
    for i in range(n):
        diff = max(diff, abs(oldX[i] - newX[i]))
    return diff


######################################################################

def thomas(A, B, n):
    beta = np.array([float(0) for _ in range(n)])
    gamma = np.array([float(0) for _ in range(n)])
    beta[0] = -A[0,1]/A[0,0]
    gamma[0] = B[0]/A[0,0]
    for i in range(1,n):
        denominator = (A[i, i-1]*beta[i-1]+A[i, i])
        if i == n-1:
            beta[i] = 0
        else:
            beta[i] = -A[i, i+1]/denominator
        gamma[i] = (B[i]-A[i, i-1]*gamma[i-1])/denominator

    x = B.copy()
    x[n-1] = gamma[n-1]
    for i in range(n-2, -1, -1):
        x[i] = beta[i]*x[i+1] + gamma[i]
    return x


def gauss(A, B, n):
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            L = A[i, k] / A[k, k]
            for j in range(k + 1, n):
                A[i, j] = A[i, j] - L * A[k, j]
            B[i] = B[i] - L * B[k]

    x = B.copy()
    for i in range(n - 1, -1, -1):
        S = 0.0
        for j in range(i + 1, n):
            S = S + A[i, j] * x[j]
        x[i] = (x[i] - S) / A[i, i]
    return x


def create_X(n):
    return np.array([float(-1) ** i for i in range(n)])


def create_matrix_A(n):
    A = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i == 0:
                A[i, j] = 1
            else:
                A[i, j] = (1 / (i + j + 1))
    return A


def create_matrix_A2(n):
    A = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if j >= i:
                A[i, j] = 2 * (i + 1) / (j + 1)
            else:
                A[i, j] = A[j, i]
    return A


def create_matrix_A3(n, k, m):
    A = np.array([[float(0) for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = k
            elif j == i + 1:
                A[i, j] = 1 / (i + 1 + m)
            elif j == i - 1:
                A[i, j] = k / (i + 1 + m + 1)
            else:
                A[i, j] = 0
    return A


def create_vector_B(A, n):
    X = create_X(n)
    B = np.array([float(0) for _ in range(n)])
    for i in range(n):
        for j in range(n):
            B[i] += A[i, j] * X[j]
    return B


def find_diff(newX, n):
    oldX = create_X(n)
    diff = 0
    for i in range(n):
        diff = max(diff, abs(oldX[i] - newX[i]))
    return diff
