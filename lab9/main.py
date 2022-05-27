from utils import *

def zad1(max, min=2):
    print("wielkość macierzy", end=" & ")
    print("float", end=" & ")
    print("decimal", end=" \\\\ \\hline\n")
    for n in range(min, max+1):
        print(n, end=" & ")
        A = create_matrix_A(n)
        B = create_vector_B(A, n)
        X = gauss(A, B, n)
        print(f'{find_diff(X, n):.4e}', end=" & ")
        A = create_matrix_A_decimal(n)
        B = create_vector_B_decimal(A, n)
        X = gauss_decimal(A, B, n)
        print(f'{find_diff_decimal(X, n):.4e}', end=" \\\\ \\hline\n")

def zad2(max, min=2):
    # print("wielkość macierzy", end=" & ")
    # print("błąd zad2.", end=" & ")
    # print("błąd zad1.", end=" \\\\ \\hline\n")

    for n in range(min, max+1):
        print(n, end=" & ")
        A = create_matrix_A2(n)
        B = create_vector_B(A, n)
        X = gauss(A, B, n)
        print(f'{find_diff(X, n):.4e}', end=" & ")

        A = create_matrix_A(n)
        B = create_vector_B(A, n)
        X = gauss(A, B, n)
        print(f'{find_diff(X, n):.4e}', end=" \\\\ \\hline\n")

def zad3(max, min = 2, k=4, m=4):
    # print("wielkość macierzy", end=" & ")
    # print("algorytm gaussa", end=" & ")
    # print("czas obliczeń", end=" & ")
    # print("algorytm thomasa", end=" & ")
    # print("czas obliczeń", end=" \\\\ \\hline\n")
    # print("algorytm thomasa", end=" \\\\ \\hline\n")
    times=[[],[]]
    for n in range(min, max+1):
        print(n, end=" & ")
        A = create_matrix_A3(n, k, m)
        B = create_vector_B(A, n)
        start = calc_time()
        X = gauss(A, B, n)
        time1 = calc_time(start)
        times[0].append(time1)
        # print(f'{find_diff(X, n):.4e}', end=" & ")
        print(f'{time1:.4e}', end=" & ")

        B = create_vector_B(A, n)
        start = calc_time()
        X = thomas(A, B, n)
        time2 = calc_time(start)
        times[1].append(time2)
        # print(f'{find_diff(X, n):.4e}', end=" & ")
        print(f'{time2:.4e}', end=" \\\\ \\hline\n")
        # print(f'{find_diff(X, n):.4e}', end=" \\\\ \\hline\n")

    df = pd.DataFrame(times,["gauss", "thomas"],[i for i in range(min,max+1)])
    plot.subplots(1)
    hm = sns.heatmap(df, center=0.1, cbar_kws=dict(use_gridspec=False, location="bottom"))
    hm.set(xlabel='Wielkość macierzy', ylabel='Rodzaj algorytmu')
    # plot.show()


# zad1(20)
# zad2(200)
for i in range(300,1001, 100):
    zad2(i, i)

