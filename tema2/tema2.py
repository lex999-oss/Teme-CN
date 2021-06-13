import numpy.linalg  # NumPy Linear Algebra Library
from math import sqrt


def cholesky(A):
    """Performs a Cholesky decomposition of A, which must
    be a symmetric and positive definite matrix. The function
    returns the lower variant triangular matrix, L."""
    n = len(A)

    # Create zero matrix for L
    L = [[0.0] * n for i in range(n)]

    # Perform the Cholesky decomposition
    for i in range(n):
        for k in range(i + 1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in range(k))

            if i == k:  # Diagonal elements
                # LaTeX: l_{kk} = \sqrt{ a_{kk} - \sum^{k-1}_{j=1} l^2_{kj}}
                L[i][k] = sqrt(A[i][i] - tmp_sum)
            else:
                # LaTeX: l_{ik} = \frac{1}{l_{kk}} \left( a_{ik} - \sum^{k-1}_{j=1} l_{ij} l_{kj} \right)
                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L


def forward_substitution(L, b):
    # Get number of rows
    n = L.shape[0]

    # Allocating space for the solution vector
    y = numpy.zeros_like(b, dtype=numpy.double)

    # Here we perform the forward-substitution.
    # Initializing  with the first row.
    y[0] = b[0] / L[0, 0]

    # Looping over rows in reverse (from the bottom  up),
    # starting with the second to last row, because  the
    # last row solve was completed in the last step.
    for i in range(1, n):
        y[i] = (b[i] - numpy.dot(L[i, :i], y[:i])) / L[i, i]

    return y


def back_substitution(U, y):
    # Number of rows
    n = U.shape[0]

    # Allocating space for the solution vector
    x = numpy.zeros_like(y, dtype=numpy.double)

    # Here we perform the back-substitution.
    # Initializing with the last row.
    x[-1] = y[-1] / U[-1, -1]

    # Looping over rows in reverse (from the bottom up),
    # starting with the second to last row, because the
    # last row solve was completed in the last step.
    for i in range(n - 2, -1, -1):
        x[i] = (y[i] - numpy.dot(U[i, i:], x[i:])) / U[i, i]

    return x


def partea1(A, B):
    """

    :param A: matricea A
    :param B: matricea B

    Rezolvarea primelor 4 subpuncte ale temei folosind metode proprii atat pentru descompunerea matricei A,
    cat si pentru substitutia directa si inversa
    """
    L = numpy.array(cholesky(A))

    LT = numpy.transpose(L)

    detA = numpy.linalg.det(A)  # det(A)
    detA_x = numpy.dot(numpy.linalg.det(L),
                       numpy.linalg.det(LT))  # calculam det(A) folosind matricea inf. triunghiulara
    # si transpusa acesteia

    print("Determinantul lui A folosind desc. Cholesky: %d" % detA_x)
    print("Determinantul lui A: %d" % detA)

    # A * x = B

    '''
    L * LT = A
    => y = L * B si x = LT * B
    '''

    y = forward_substitution(L, B)  # substitutie directa
    x = back_substitution(LT, y)  # substitutie inversa

    print("Solutia aproximata pentru A * x = B: ")
    print(x)
    check_x = numpy.linalg.norm(numpy.subtract(numpy.dot(A, x), B), ord=None, axis=None,
                                keepdims=False)  # euclid_norm(A * x - B)
    if check_x < 10 ** (-8):
        print("Solutia X este buna!")


def partea2(A, B):
    """

    :param A: matricea A
    :param B: matricea B

    Rezolvarea subpunctelor care implica folosirea bibliotecii NumPy
    """
    L = numpy.linalg.cholesky(A)

    LT = L.T.conj()

    # A * x = B

    '''
    L * LT = A
    => y = L * B si x = LT * B
    '''

    y = numpy.linalg.solve(L, B)  # substitutie directa
    x = numpy.linalg.solve(LT, y)  # substitutie inversa

    print("Solutia cu NumPy pentru A * x = B: ")
    print(x)

    '''
    A = L*LT => A_inv = L_inv * LT_inv
    '''

    L_inv = numpy.linalg.inv(L)
    LT_inv = numpy.linalg.inv(LT)

    A_I = numpy.dot(LT_inv, L_inv)
    A_I_direct = numpy.linalg.inv(A)

    print("A inversat:")
    print(A_I)

    check_x = numpy.linalg.norm(numpy.subtract(A_I, A_I_direct), ord=None, axis=None, keepdims=None)

    if check_x < 10 ** (-8):
        print("Diferenta este neglijabila -- euclid_norm( (A_Inv_Chol - A_Inv_bibl) )")


def main():
    """
    Functia main a programului.
    Preia informatii de la linia de comanda si apeleaza functiile partea1() si partea(2) pentru a rezolva subpunctele
    din problema.
    :return:
    """
    dataA = []
    dataB = []

    n = int(input("Matrice N*N, N = "))
    if n < 0:
        print("Valorile negative nu sunt permise!")
        exit(0)

    print("Introduceti datele pentru matricea A!")
    for i in range(n):
        for j in range(n):
            v = int(input())
            if v < 0:
                print("Valorile negative nu sunt permise!")
                exit(0)
            else:
                dataA.append(v)

    A = numpy.array(dataA).reshape(n, n)  # matricea A
    # A = numpy.array([[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]])
    print("Introduceti datele pentru matricea B!")
    for i in range(n):
        v = int(input())
        if v < 0:
            print("Valorile negative nu sunt permise!")
            exit(0)
        else:
            dataB.append(v)

    B = numpy.array(dataB)  # matricea B
    # B = numpy.array([1,2,3,4])

    try:
        L = numpy.linalg.cholesky(A)  # matricea inferior triunghiulara
    except numpy.linalg.LinAlgError:
        print("Matricea nu este pozitiv definita!")
        exit(0)

    partea1(A, B)
    partea2(A, B)


if __name__ == "__main__":
    main()
