from utils import *


def zad1(max, min=2):
    for n in range(min, max+1):
        # print(n, end=" & ")
        A = create_matrix_A(n)
        B = create_vector_B(A, n)
        jacoby(A, np.array([float(0) for _ in range(n)]))
        print(create_X(n))
        print(A)
        print(B)
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

zad1(4)