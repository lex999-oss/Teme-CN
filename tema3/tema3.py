import numpy as np

A_file_handle = open('a.txt', 'r')
n = int(A_file_handle.readline())

A = np.zeros(shape=(n, n))

A_file_handle.readline()

for line in A_file_handle.readlines():
    line = line.split(", ")
    num = float(line[0])
    i_1 = int(line[1])
    i_2 = int(line[2])
    A[i_1][i_2] = num

A_file_handle.close()

a, b, c = [], [], []

B_file_handle = open('b.txt', 'r')

n = int(B_file_handle.readline())
p = int(B_file_handle.readline())
q = int(B_file_handle.readline())

B_file_handle.readline()

i = 0
while i < n:  # diagonala "a"
    x = float(B_file_handle.readline())
    a.append(x)
    i += 1

B_file_handle.readline()

i = 0
while i < n - p:  # diagonala "b"
    x = float(B_file_handle.readline())
    b.append(x)
    i += 1

B_file_handle.readline()

i = 0
while i < n - q:  # diagonala "c"
    x = float(B_file_handle.readline())
    c.append(x)
    i += 1

B_file_handle.close()

# if n < len(A):
n = len(A)

C = np.array(A)

for i in range(n): # adunarea se face foarte rapid folosind aceasta metoda
    for j in range(n):
        if i == j:
            C[i][j] += a[i]
        elif j - i == q:
            C[i][j] += b[i]
        elif i - j == p:
            C[i][j] += c[i - 1]

A_plus_B = []

for i in range(n):
    for j in range(n):
        if C[i][j] != 0:
            z = [C[i][j], i, j]
            A_plus_B.append(z)

result_file_handler = open("my_aplusb.txt", "w")

for elem in A_plus_B:
    res_str = f'{elem} \n'
    result_file_handler.write(res_str)
result_file_handler.close()

B = np.zeros(shape=(n, n))
for i in range(n):
    for j in range(n):
        if i == j:
            B[i][j] = a[i]
        elif j - i == q:
            B[i][j] = b[i]
        elif i - j == p:
            B[i][j] = c[i - 1]

V = A @ B

A_ori_B = []

for i in range(n):
    for j in range(n):
        if V[i][j] != 0:
            z = [V[i][j], i, j]
            A_ori_B.append(z)

result_file_handler = open("my_aorib.txt", "w")

for elem in A_ori_B:
    res_str = f'{elem} \n'
    result_file_handler.write(res_str)
result_file_handler.close()
