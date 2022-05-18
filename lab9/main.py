n = 10
A = [[0 for _ in range(n)]for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i==0:
            A[i][j] = 1
        else:
            A[i][j] = 1/(i+j+1)

for i in A:
    print(i)
