from utils import *


def zad1(max, criterion, min=2):
    for n in range(min, max+1):
        for eps in 0.1, 0.001, 0.00001:
            # print(n, end=" & ")
            A = create_matrix_A(n)
            B = create_vector_B(A, n)
            print(create_X(n))
            print()
            print(A)
            print()
            print(B)
            print()
            iterations, X = jacoby(A, np.array([float(0) for _ in range(n)]), B, criterion, eps)
            print(X)
            print()
            print(iterations)
            print("====")
            # start = calc_time()
            # X = gauss(A, B, n)
            # time1 = calc_time(start)
            # times[0].append(time1)
            # print(find_diff(X, n), end=" & ")
            # print(time1, end=" & ")
            #
            # B = create_vector_B(A, n)
            # start = calc_time()
            # X = thomas(A, B, n)
            # time2 = calc_time(start)
            # times[1].append(time2)
            # print(find_diff(X, n), end=" & ")
            # print(time2, end=" \\\\ \\hline\n")
            # print(find_diff(X, n), end=" \\\\ \\hline\n")

zad1(4, 0, 4)