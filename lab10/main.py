from utils import *


def zad1(ranges, x, time):
    print("wielkość macierzy", end=" & ")
    print("ro=1.00e-01", end=" & ")
    print("ro=1.00e-03", end=" & ")
    print("ro=1.00e-05", end=" & ")
    print("ro=1.00e-07", end=" \\\\ \\hline\n")
    for n in ranges:
        X_init = create_X(n)
        X = np.array([float(x) for _ in range(n)])
        A = create_matrix_A(n)
        B = create_vector_B(A, n)
        N, M = jacoby(A)
        print(n, end=" & ")
        for eps in 0.1, 0.001, 0.00001, 0.0000001:
            iterations, X_ans, T = jacoby_iter(A, M, N, X, B, 0, eps)
            if time==0:
                if eps!=0.0000001:
                    print(iterations, f'{f1(X_init, X_ans):.2e}', sep=", ", end=" & ")
                else:
                    print(iterations, f'{f1(X_init, X_ans):.2e}', sep=", ", end=" \\\\ \\hline\n")
            else:
                if eps!=0.0000001:
                    print(f'{T:.2e}', sep=", ", end=" & ")
                else:
                    print(f'{T:.2e}',  sep=", ", end=" \\\\ \\hline\n")
        print(end=" & ")
        for eps in 0.1, 0.001, 0.00001, 0.0000001:
            iterations, X_ans, T = jacoby_iter(A, M, N, X, B, 1, eps)
            if time==0:
                if eps!=0.0000001:
                    print(iterations, f'{f1(X_init, X_ans):.2e}', sep=", ", end=" & ")
                else:
                    print(iterations, f'{f1(X_init, X_ans):.2e}', sep=", ", end=" \\\\ \\hline\n")
            else:
                if eps!=0.0000001:
                    print(f'{T:.2e}', sep=", ", end=" & ")
                else:
                    print(f'{T:.2e}',  sep=", ", end=" \\\\ \\hline\n")

def zad2(ranges):
    print("wielkość macierzy", end=" & ")
    print("promień spektralny", end=" \\\\ \\hline\n")
    for n in ranges:
        print(n, end=" & ")
        A = create_matrix_A(n)
        N, M = jacoby(A)
        rey = max(np.linalg.eigvals(M))
        print(f'{rey:.5e}', end=" \\\\ \\hline\n")

# zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 0, 0)
# zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 0, 1)
# zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 10, 0)
# zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 10, 1)
# zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 100, 0)
zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 100, 1)

# zad2([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500])
# np.linalg.eigvals()

