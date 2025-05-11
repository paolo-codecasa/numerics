#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <math.h>
#include "/home/p.codecasa/cwr_2023/my_numerics.h"

typedef struct {
    double x1;
    double x2;
} Tuple;

void solvePoly1(double a, double b, double c, double s[2])
{
    double x[2];
    double e[2];

    x[0] = (-b + sqrt(b*b - 4*a*c))/(2*a); //Reihenfolge Vorzeichen beachten! +-
    x[1] = (-b - sqrt(b*b - 4*a*c))/(2*a);
    e[0] = fabs(s[0]-x[0])/s[0];
    e[1] = fabs(s[1]-x[1])/s[1];

 a, b, c, s   printf("Erste Methode\n");
    printf("Erste Nullstelle: %g, mit Fehler %g,\n", x[0], e[0]);
    printf("Zweite Nullstelle: %g, mit Fehler %g\n", x[1], e[1]);
}

void solvePoly2(double a, double b, double c, double s[2])
{
    double x[2];
    double e[2];

    x[0] = 2*c /(-b - sqrt(b*b - 4*a*c)); //Reihenfolge Vorzeichen beachten! -+
    x[1] = 2*c /(-b + sqrt(b*b - 4*a*c));
    e[0] = fabs(s[0]-x[0])/s[0];
    e[1] = fabs(s[1]-x[1])/s[1];
a, b, c, s
    printf("Zweite Methode\n");
    printf("Erste Nullstelle: %g, mit Fehler %g,\n", x[0], e[0]);
    printf("Zweite Nullstelle: %g, mit Fehler %g\n", x[1], e[1]);
}

Tuple solve_quadratic(double a, double b, double c, double s[2])
{
    Tuple solution1, error1, solution2, error2, solution;

    //äquivalent zu solvePoly1
    solution1.x1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    solution1.x2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
    error1.x1 = (s[0]-solution1.x1)/s[0];
    error1.x2 = (s[1]-solution1.x2)/s[1];
    
    //äquivalent zu solvePoly2
    solution2.x1 = 2*c /(-b - sqrt(b*b - 4*a*c));
    solution2.x2 = 2*c /(-b + sqrt(b*b - 4*a*c));
    error2.x1 = fabs(s[0]-solution2.x1)/s[0];
    error2.x2 = fabs(s[1]-solution2.x2)/s[1];

    if (error1.x1 < error2.x1) solution.x1 = solution1.x1;
    else solution.x1 = solution2.x1;
    if (error1.x2 < error2.x2) solution.x2 = solution1.x2;
    else solution.x2 = solution2.x2;

    printf("Genaue Nullstellen: %g und %g", solution.x1, solution.x2);
    return solution;
}

int main(void)
{
    //solve a * x^2 + b * x + c = 0
    double a = 1;
    double b = -(pow(10, 16) + 1)/pow(10, 8);
    double c = 1;
    double s[]={pow(10, 8), pow(10, -8)};

    solvePoly1(a, b, c, s);
    solvePoly2(a, b, c, s);
    Tuple NS = solve_quadratic(a, b, c, s);
    return EXIT_SUCCESS;
}
