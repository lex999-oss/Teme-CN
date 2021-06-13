import numpy

def roots(coefs, x0):
    x = 0
    k = 0
    dx = 0
    max_it = 100
    eps = 10**(-8)
    p = numpy.poly1d(coefs)
    while eps <= abs(dx) <= 10 ** 8 and k < max_it:
        if abs(p.deriv(1)) <= eps:
            break
        c = ((p**2) * (p.deriv(2))) / (p.deriv(1) ** 3)
        dx = p / p.deriv(1) + c/2
        x -= dx
        k += 1
    if abs(dx) < eps:
        print(f'Convergenta pentru x* = {x}')
    else:
        print('Divergenta')

def schema_horner(c, x0):
    d = c[0]
    for coeficient in c[1:]:
        d = coeficient + d * x0
    return d

if __name__ == '__main__':
    coefs = [1.0, -6.0, 11.0, -6.0]
    x0 = float(input('Punct de start:'))
    roots(coefs, x0)