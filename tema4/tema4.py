import math

EPSILON = 10 ** -10

N = 0
P = 0
Q = 0


def read_matrix(path1, path2, a, b, c, B):
    global N
    global P
    global Q

    # open file for reading a
    with open(path1, 'r') as f1:
        n1 = int(f1.readline())
        P = int(f1.readline())
        Q = int(f1.readline())
        f1.readline()

        i = 0
        while i < n1:  # diagonala "a"
            x = float(f1.readline())
            a.append(x)
            i += 1

        f1.readline()

        i = 0
        while i < n1 - Q:  # diagonala "c"
            x = float(f1.readline())
            c.append(x)
            i += 1

        f1.readline()

        i = 0
        while i < n1 - P:  # diagonala "b"
            x = float(f1.readline())
            b.append(x)
            i += 1

    # open file for reading b
    f2 = open(path2, 'r')
    n2 = int(f2.readline())
    if n1 == n2:
        N = n1
    f2.readline()

    for i in range(n2):
        x = float(f2.readline())
        B.append(x)
    f1.close()
    f2.close()


def verify_diagonal(a):
    global N
    ok = 0

    for i in range(N):
        if a[i] == 0:
            ok = 1
    if ok == 1:
        return False
    else:
        return True


def gaussSeidel(x, a, b, c, B):
    global P
    global Q
    global N

    maxIterations = 10**2
    for i in range(maxIterations):
        dx = 0
        for j in range(N):
            xprev = x[j]
            summ1 = 0
            summ2 = 0
            if j == 0:
                summ2 = b[j] * x[j + 1]
            elif j == N - 1:
                summ1 = c[j - 1] * x[j - 1]
            else:
                summ1 = c[j - 1] * x[j - 1]
                summ2 = b[j] * x[j + 1]
            x[j] = (B[j] - summ1 - summ2) / a[j]
            dx += pow((x[j] - xprev), 2)

        dx = math.sqrt(dx)

        if (dx < EPSILON) and i != 0:
            print("Sequence converges to [", end="")
            for j in range(N - 1):
                print(x[j], ",\n", end="")
            print(x[N - 1], "].\nTook", i + 1, "iterations.")
            return
    print("Doesn't converge.")


def main():
    global N
    a = []
    b = []
    c = []
    B = []

    index = int(input("Selecteaza matricea si vectorul aferent(1..5): "))

    read_matrix(f'a{index}.txt', f'f{index}.txt', a, b, c, B)
    x = [i + 1 for i in range(0, N)]
    if not verify_diagonal(a):
        print("The matrix has 0 on diagonal a")
        return
    else:
        gaussSeidel(x, a, b, c, B)

    Ax = [0] * N
    for i in range(N):  # A*x
        if i == 0:
            Ax[i] = a[i] * x[i] + b[i] * x[i+1]
        elif i == N - 1:
            Ax[i] = c[i-1] * x[i-1] + a[i] * x[i]
        else:
            Ax[i] = c[i - 1] * x[i - 1] + a[i] * x[i] + b[i] * x[i + 1]

    Ax_minus_f = []
    for i in range(N):  # A*x - f
        v = math.fabs(Ax[i] - B[i])
        Ax_minus_f.append(v)
    inf_norm = max(Ax_minus_f)

    print(f'Norma infinita pentru Ax-f este: {inf_norm} \n')


if __name__ == '__main__':
    main()
