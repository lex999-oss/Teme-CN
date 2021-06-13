#include <stdio.h>

double ex1()
{
    double x = 1.0;
    double u = 1;
    while(1)
    {
        if (x + u == x)
            break;
        else {
            u = u / 10.0;
        }
    }
    printf("u: %2.21lf\n", u);
    return u;
}

int main() {
    double u = ex1();
    double a = 1.0;
    double b = u / 10;
    double c = b;

    if ((a + b) + c != a + (b + c))
        printf("Operatia + nu este asociativa!");
}
