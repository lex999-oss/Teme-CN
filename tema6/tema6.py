import numpy


def f(x):
    return x ** 2 - 12 * x + 30


def least_squares(n, m, x0, xn, x):
    h = (xn - x0) / n
    points = [(x0, f(x0))]

    for i in range(1, n):
        points.append((x0 + i * h, f(x0 + i * h)))

    points.append((xn, f(xn)))
    print(points)
    a = numpy.zeros(shape=(m, m))
    for i in range(0, m):
        for j in range(0, m):
            sum = 0
            for point in points:
                sum += point[0] ** (i + j)
            a[i, j] = sum

    b = numpy.zeros(shape=(m,))
    for i in range(0, m):
        sum = 0
        for point in points:
            sum += point[1] * (point[0] ** i)
        b[i] = sum
    c = numpy.linalg.solve(a, b)
    c1 = numpy.array(c)
    for i in range(0, len(c)):
        c1[i] = c[len(c) - 1 - i]
    print(a)
    print('x =', x)
    print('f(x) =', f(x))
    print('Pm(x) =', schema_horner(c1, x))
    print('|Pm(x) - f(x)| =', abs(schema_horner(c1, x) - f(x)))


def schema_horner(c, x0):
    d = c[0]
    for coeficient in c[1:]:
        d = coeficient + d * x0
    return d


x0 = int(input("x0: "))
xn = int(input("xn: "))
least_squares(5, 3, x0, xn, 3)
