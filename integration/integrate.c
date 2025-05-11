#include <stdio.h>
#include <stdlib.h>

double p1(double y)
{
return y*y*y - y/2;
}

double integrate(double left, double right,
int N, double integrand(double))
{
double sum = 0;
double delta_y = (right - left) / N;
for(int k = 0; k < N; ++k) {
double y = left + k * delta_y;
double f = integrand(y);
double A = f * delta_y;
sum += A;
}
return sum;
}

int main(void)
{
double right = 20.0;   //Eingabewert (rechte Integrationsgrenze)
double left =0.0;  //linke Integrationsgrenze
int N = 10;  //hier ändern für resolution
double num = integrate(left, right, N, p1);
double ana = right*right*right*right/4 - right*right/4 ;
printf("Integralwert numerisch: %f\n", num);
printf("Integralwert analytisch: %f\n", ana);
double abweichung;
printf("Abweichung: %f\n", (num-ana)/ana);
return EXIT_SUCCESS;
}
