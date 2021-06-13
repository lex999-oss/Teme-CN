from math import tan, pi
import random
import time


def ex1():
    u = 1
    while True:
        if 1.0 + u == 1.0:
            break
        else:
            u = u / 10
    return u


def ex2():
    u = ex1()
    a = 1.0
    b = u
    c = b
    if (a + b) + c != a + (b + c):
        print("Operatia + nu este asociativa!")
    else:
        print("Operatia + este asociativa!")

    a = 1.11
    b = 1.01
    c = b
    if (a * b) * c != a * (b * c):
        print("Operatia * nu este asociativa!")
    else:
        print("Operatia * este asociativa!")


p1 = 1/3
p2 = 2/15
p3 = 17/315
p4 = 62/2853


def ex3_tan_poly(number):
    n = number ** 2
    n1 = number ** 3
    n2 = n1 * n
    n3 = n2 * n
    n4 = n3 * n
    return number + p1 * n1 + p2 * n2 + p3 * n3 + p4 * n4


ex2()
# generez 10000 numere, apoi contorizez doar tangenta pe acelasi vector
delta = 0
for i in range(1, 10000):
    x = random.uniform(-pi / 4, pi / 4)
    answer1 = ex3_tan_poly(x)
    answer2 = tan(x)
    delta += abs(answer1 - answer2)
avg_delta = delta / 10000
print("Media diferentei pentru 10000 de operatii este: %f" % avg_delta)

x = random.uniform(-pi / 4, pi / 4)
start_time = time.perf_counter_ns()
ex3_tan_poly(x)
my_tan_time = time.perf_counter_ns() - start_time

start_time = time.perf_counter_ns()
tan(x)
python_tan_time = time.perf_counter_ns() - start_time

time_delta = abs(my_tan_time - python_tan_time)
print("Timpii de rulare")
print("Pentru My_Tan: %.2f ns" % my_tan_time)
print("Pentru python_tan: %.2f ns" % python_tan_time)
print("Diferenta de timp: %.2f ns" % time_delta)
