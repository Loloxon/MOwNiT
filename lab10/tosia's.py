import math
import numpy
import numpy as np
import random
import time
import xlsxwriter

e = 0.00000001


def jacob(A, b):
    X = [0 for _ in range(len(A))]
    L = [[0 for _ in range(len(A))] for _ in range(len(A))]
    D = [[0 for _ in range(len(A))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A)):
            if i != j:
                L[i][j] = A[i][j]
            if i == j:
                val = 1.0/A[i][j]
                D[i][j] = val
    X_n = np.dot(D, numpy.subtract(b, np.dot(L, X)))
    iterations = 1
    while max([math.fabs(X[i] - X_n[i]) for i in range(len(X))]) > e:
        iterations += 1
        X, X_n = np.dot(D, numpy.subtract(b, np.dot(L, X_n))), X
    return X, iterations


def jacob2(A, b):
    X = [0 for _ in range(len(A))]
    L = [[0 for _ in range(len(A))] for _ in range(len(A))]
    D = [[0 for _ in range(len(A))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A)):
            if i != j:
                L[i][j] = A[i][j]
            if i == j:
                val = 1.0/A[i][j]
                D[i][j] = val
    X_n = np.dot(D, numpy.subtract(b, np.dot(L, X)))
    iterations = 1
    while max(np.subtract(np.dot(A, X), b)) > e:
        iterations += 1
        X, X_n = np.dot(D, numpy.subtract(b, np.dot(L, X_n))), X
    return X, iterations


def matrix_1(n, m, k):
    A = []
    for i in range(n):
        temp = []
        for j in range(n):
            if i == j:
                temp.append(k)
            else:
                temp.append(m/(n - i - j + 0.5))
        A.append(temp)
    return A


workbook = xlsxwriter.Workbook('results.xlsx')
sheet = workbook.add_worksheet()
s = [3, 4, 5, 7, 10, 13, 15, 18, 20, 25, 30, 50, 100, 150, 200, 300, 500]
r = 0
for size in s:
    print(size)
    r += 1
    matrix = matrix_1(size, 1, 10)
    xs = []
    for i in range(size):
        xs.append(random.choice([-1.0, 1.0]))
    bs = np.dot(matrix, xs)
    t_start = time.time()
    ans, k = jacob(matrix, bs)
    t_end = time.time()
    max_diff = max([math.fabs(ans[i] - xs[i]) for i in range(size)])
    sheet.write(r, 1, str(size))
    sheet.write(r, 2, str(k))
    sheet.write(r, 3, str(t_end - t_start))
    sheet.write(r, 4, str(round(max_diff, 17)))
    t_start = time.time()
    ans, k = jacob2(matrix, bs)
    t_end = time.time()
    max_diff = max([math.fabs(ans[i] - xs[i]) for i in range(size)])
    sheet.write(r, 6, str(size))
    sheet.write(r, 7, str(k))
    sheet.write(r, 8, str(t_end - t_start))
    sheet.write(r, 9, str(round(max_diff, 17)))
workbook.close()