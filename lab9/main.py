from utils import *

def zad1(max, min=2):
    for n in range(min, max+1):
        A = create_matrix_A(n)
        B = create_vector_B(A, n)
        X = gauss(A, B, n)
        print(find_diff(X, n), end=" & ")
        A = create_matrix_A_decimal(n)
        B = create_vector_B_decimal(A, n)
        X = gauss_decimal(A, B, n)
        print(find_diff_decimal(X, n), end=" \\\\ \\hline\n")

def zad2(max, min=2):
    for n in range(min, max+1):
        A = create_matrix_A2(n)
        # print(A)
        B = create_vector_B(A, n)
        X = gauss(A, B, n)
        print(find_diff(X, n), end=" \\\\ \\hline\n")

def zad3(max, min = 2, k=4, m=4):
    for n in range(min, max+1):
        A = create_matrix_A3(n, k, m)
        # print(A)
        B = create_vector_B(A, n)
        X = gauss(A, B, n)
        # print(X)
        print(find_diff(X, n), end=" & ")
        B = create_vector_B(A, n)
        X = thomas(A, B, n)
        # print(X)
        print(find_diff(X, n), end=" \\\\ \\hline\n")


# zad1(20)
# zad2(20)
zad3(200,200)

