from utils import *


def zad1(ranges, x, time):
    print("wielkość macierzy", end=" & ")
    print("ro=0.1", end=" & ")
    print("ro=0.001", end=" & ")
    print("ro=0.00001", end=" & ")
    print("ro=0.0000001", end=" \\\\ \\hline\n")
    for n in ranges:
        X = np.array([float(x) for _ in range(n)])
        print(n, end=" & ")
        for eps in 0.1, 0.001, 0.00001, 0.0000001:
            A = create_matrix_A(n)
            B = create_vector_B(A, n)

            iterations, X, T= jacoby(A, X, B, 0, eps)
            iterations2, X2, T2 = jacoby(A, X, B, 1, eps)

            if time==0:
                if eps!=0.0000001:
                    print(iterations, iterations2, sep=", ", end=" & ")
                else:
                    print(iterations, iterations2, sep=", ", end=" \\\\ \\hline\n")
            else:
                if eps!=0.0000001:
                    print(f'{T:.2e}', f'{T2:.2e}', sep=", ", end=" & ")
                else:
                    print(f'{T:.2e}', f'{T2:.2e}', sep=", ", end=" \\\\ \\hline\n")

def zad2(ranges):
    print("wielkość macierzy", end=" & ")
    print("promień spektralny", end=" \\\\ \\hline\n")
    for n in ranges:
        X = np.array([float(0) for _ in range(n)])
        print(n, end=" & ")
        A = create_matrix_A(n)
        M = getM(A, X)
        rey = max(np.linalg.eigvals(M))
        print(f'{rey:.5e}', end=" \\\\ \\hline\n")

# zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 100, 0)
# zad1([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500], 100, 1)

zad2([3, 4, 5, 7, 9, 11, 15, 19, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500])
# np.linalg.eigvals()

