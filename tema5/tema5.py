import numpy as np
import platform


def jacobi_eigenvalue(n, a, it_max):
    """
    :param n: order of matrix A
    :param a: matrix A(N, N) - square and symmetric
    :param it_max: Maximum number of iterations
    :return: U(N, N) - matrix of eigenvectors
             LAMBDA - vector of eigenvalues
             IT_NUM - the total number of iterations
             ROT_NUM - the total number of rotations
    """
    u = np.zeros([n, n])
    lam = np.zeros(n)

    for j in range(0, n):
        for i in range(0, n):
            u[i, j] = 0.0
        u[j, j] = 1.0

    for i in range(0, n):
        lam[i] = a[i, i]

    bw = np.zeros(n)
    zw = np.zeros(n)

    for i in range(0, n):
        bw[i] = lam[i]

    it_num = 0
    rot_num = 0

    while it_num < it_max:

        it_num = it_num + 1
        #
        #  The convergence threshold is based on the size of the elements in
        #  the strict upper triangle of the matrix.
        #
        thresh = 0.0
        for j in range(0, n):
            for i in range(0, j):
                thresh = thresh + a[i, j] ** 2

        thresh = np.sqrt(thresh) / float(4 * n)

        if thresh == 0.0:
            break

        for p in range(0, n):
            for q in range(p + 1, n):

                gapq = 10.0 * abs(a[p, q])
                termp = gapq + abs(lam[p])
                termq = gapq + abs(lam[q])
                #
                #  Annihilate tiny offdiagonal elements.
                #
                if 4 < it_num and termp == abs(lam[p]) and termq == abs(lam[q]):

                    a[p, q] = 0.0
                #
                #  Otherwise, apply a rotation.
                #
                elif thresh <= abs(a[p, q]):

                    h = lam[q] - lam[p]
                    term = abs(h) + gapq

                    if term == abs(h):
                        t = a[p, q] / h
                    else:
                        theta = 0.5 * h / a[p, q]
                        t = 1.0 / (abs(theta) + np.sqrt(1.0 + theta * theta))
                        if theta < 0.0:
                            t = - t

                    c = 1.0 / np.sqrt(1.0 + t * t)
                    s = t * c
                    tau = s / (1.0 + c)
                    h = t * a[p, q]
                    #
                    #  Accumulate corrections to diagonal elements.
                    #
                    zw[p] = zw[p] - h
                    zw[q] = zw[q] + h
                    lam[p] = lam[p] - h
                    lam[q] = lam[q] + h

                    a[p, q] = 0.0
                    #
                    #  Rotate, using information from the upper triangle of A only.
                    #
                    for j in range(0, p):
                        g = a[j, p]
                        h = a[j, q]
                        a[j, p] = g - s * (h + g * tau)
                        a[j, q] = h + s * (g - h * tau)

                    for j in range(p + 1, q):
                        g = a[p, j]
                        h = a[j, q]
                        a[p, j] = g - s * (h + g * tau)
                        a[j, q] = h + s * (g - h * tau)

                    for j in range(q + 1, n):
                        g = a[p, j]
                        h = a[q, j]
                        a[p, j] = g - s * (h + g * tau)
                        a[q, j] = h + s * (g - h * tau)
                    #
                    #  Accumulate information in the eigenvector matrix.
                    #
                    for j in range(0, n):
                        g = u[j, p]
                        h = u[j, q]
                        u[j, p] = g - s * (h + g * tau)
                        u[j, q] = h + s * (g - h * tau)

                    rot_num = rot_num + 1

        for i in range(0, n):
            bw[i] = bw[i] + zw[i]
            lam[i] = bw[i]
            zw[i] = 0.0
    #
    #  Restore upper triangle of input matrix.
    #
    for j in range(0, n):
        for i in range(0, j):
            a[i, j] = a[j, i]

    return u, lam, it_num, rot_num


def is_eigen_right(n, k, a, u, lam):
    """
    :param n: order of matrix A
    :param k: number of eigenvectors
    :param a: matrix A(N, N)
    :param u: matrix X(N, K) of eigenvectors
    :param lam: vector LAM(K) of eigenvalues
    :return: The Fobenius Norm of A * U - U * LAMBDA
    """
    c = mat_mul(n, n, k, a, u)

    for j in range(0, k):
        for i in range(0, n):
            c[i, j] = c[i, j] - lam[j] * u[i, j]

    value = norm_fro(c)

    return value


def mat_mul(n1, n2, n3, a, b):
    """
    :param a: matrix A(N1, N2)
    :param b: matrix B(N2, N3)
    :return: product matrix C(N1, N3) = A * B
    """
    c = np.zeros((n1, n3))

    for j in range(0, n3):
        for i in range(0, n1):
            for k in range(0, n2):
                c[i, j] = c[i, j] + a[i, k] * b[k, j]

    return c


def norm_fro(a):
    """
    :param a: Matrix A
    :return: Frobenius Norm of Matrix A
    """
    return np.sqrt(sum(sum(a ** 2)))


def cholesky_contig_decomp(a):
    """
    :param a: matrix A
    :return: A**(k) cholesky decomp.
    A**(k) = transpose(L**(k-1)) * L**(k-1)
    """
    it_max = 100
    it_num = 0
    eps = 10 ** (-4)
    try:
        L = np.linalg.cholesky(a)
    except np.linalg.LinAlgError:
        return "Matrix is not positive definite"
    LT = L.T.conj()
    a1 = np.dot(LT, L)
    temp = a1
    while it_num < it_max:
        temp1 = temp
        L = np.linalg.cholesky(temp1)
        LT = L.T.conj()
        temp = np.dot(LT, L)
        d = np.linalg.norm(np.subtract(temp, temp1), 2)
        it_num += 1
        if d < eps:
            break
    return temp


def get_max_sv(s, n):
    """
    :param n: number of elements for matrix S
    :param s: SVD array of matrix A
    :return: the rank of matrix A(the highest positive singular value)
    """
    max_sv = 0
    for i in range(n):
        if s[i] > max_sv:
            max_sv = s[i]
    return max_sv


def get_matrix_cond(s, n):
    """
    :param n: number of elements for matrix S
    :param s: SVD array of matrix A
    :return: the cond. value of matrix A(the highest positive singular value/the lowest positive singular value)
    """
    min_sv = 999999999999
    max_sv = get_max_sv(s, n)
    for j in range(n):
        if min_sv > s[j] > 0.0:
            min_sv = s[j]
    return max_sv / min_sv


def get_mat_si(s, n, p):
    si = np.zeros(shape=(n, p))
    for i in range(n):
        si[i, i] = 1 / s[i]
    return si


if __name__ == '__main__':
    # ------------------------ FIRST PART OF THE EXERCISE ---------------------------
    n = 4
    a = np.array([
        [4.0, -30.0, 60.0, -35.0, ],
        [-30.0, 300.0, -675.0, 420.0, ],
        [60.0, -675.0, 1620.0, -1050.0, ],
        [-35.0, 420.0, -1050.0, 700.0]])

    print('Python version: %s' % (platform.python_version()))
    print('JACOBI_EIGENVALUE computes the eigenvalues LAMBDA')
    print('and eigenvectors U of a symmetric matrix A so that A * U = U * V.')

    print('Input matrix A:')
    print(a)

    it_max = 100

    u, lam, it_num, rot_num = jacobi_eigenvalue(n, a, it_max)

    print('')
    print('Number of iterations = %d' % it_num)
    print('Number of rotations  = %d' % rot_num)

    print('Eigenvalues LAMBDA:')
    print(lam)

    print('Eigenvector matrix U:')
    print(u)
    #
    #  Compute eigentest.
    #
    error_fro = is_eigen_right(n, n, a, u, lam)
    print('')
    print('Frobenius norm error of A*U-U*V = %g' % error_fro)
    # ------------------------ SECOND PART OF THE EXERCISE ---------------------------
    print('')
    print("A**(k) matrix derived from A, following the formula:")
    print("A**(k) = L**(k-1) * transpose(L**(k-1))")
    print('')
    # compute cholesky decomposition
    a_k = cholesky_contig_decomp(a)
    print(a_k)
    # ------------------------ THIRD PART OF THE EXERCISE ---------------------------
    p = 5
    a = np.array([
        [4.0, -30.0, 60.0, -35.0],
        [-30.0, 300.0, -675.0, 420.0],
        [60.0, -675.0, 1620.0, -1050.0],
        [-35.0, 420.0, -1050.0, 700.0],
        [35.0, -300.0, 675.0, -420.0]])
    # get U, S, V from SVD of A
    u_svd, s_svd, v_svd = np.linalg.svd(a)
    # compute rank of A
    rank_A = len(s_svd)
    # compute conditional value of A
    cond_A = get_matrix_cond(s_svd, n)
    # compute inverse of matrix S
    si_svd = get_mat_si(s_svd, n, p)
    # transpose of matrix U
    u_svd_trans = np.transpose(u_svd)
    # compute the Moore-Penrose inverse of A
    v_svd = np.transpose(v_svd)
    inverse_a_mp = np.dot(np.dot(v_svd, si_svd), u_svd_trans)
    # compute the pseudo-inverse matrix of A in the sense of the smallest squares
    a_t = a.T.conj()  # transpose of A
    a_j_1 = np.dot(a_t, a)  # transpose(A) * A
    a_j_2 = np.linalg.inv(a_j_1)  # (transpose(A) * A)**(-1)
    a_j = np.dot(a_j_2, a_t)  # ((transpose(A) * A)**(-1)) * transpose(A)
    # compute the norm of A_I - A_J
    norm = np.linalg.norm(np.subtract(inverse_a_mp, a_j), ord=1)
    print('')
    print("Input matrix A:")
    print(a)
    print(f'Rank of matrix A: {rank_A}')
    print(f'Conditional value of matrix A: {cond_A}')
    print("Singular Value Decomposition of A:")
    print("matrix U:")
    print(u_svd)
    print("array S:")
    print(s_svd)
    print("matrix V:")
    print(v_svd)
    print("Moore-Penrose inverse of A(A_I):")
    print(inverse_a_mp)
    print("Pseudo-inverse matrix of A in the sense of the smallest squares(A_J):")
    print(a_j)
    print("norm(A_I - A_J):")
    print(norm)
