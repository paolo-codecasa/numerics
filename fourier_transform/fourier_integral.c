#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include "../my_numerics.h"



double fourier_integrand_cosine(double x, double k)
{
    double omega = 4;
    return cos(omega * x) * cos(k * x);
}

void do_stuff(double M, char *file_name)
{
    double right_M = M;   //rechte Integrationsgrenze
    double left_M = -M;  //linke Integrationsgrenze
    double dx = 0.001; //Resolution Näherung Integral
    //int N_M = (int)((right_M - left_M) / dx);  //Anzahl zu integrierenden Intervalle
    double right_k = 10.0;   //Anfang Intervall Funktionwerte (k ist Variable der Fourier-Transformierten)
    double left_k = -10.0;  //Ende Intervall Funktionwerte
    int N_k = 1000;  //Anzahl Messwerte
    double step_k = (right_k - left_k) / N_k; //Schritt von Messwert zu Messwert
    FILE* myFile = fopen(file_name, "w");
    for (int i=0; i < N_k; i++){
        double k = left_k + i * step_k;
        double integral = mn_integrate_simpson(left_M, right_M, dx, fourier_integrand_cosine, k);
        fprintf(myFile, "%g, %g\n", k, integral);
    }
    fclose(myFile);
}

int main(void)
{
    //double M = 10;   //Integrationsgrenze hier einstellen

    do_stuff(10, "data10.csv");
    do_stuff(30, "data30.csv");
    do_stuff(100, "data100.csv");
    do_stuff(500, "data500.csv");

    return EXIT_SUCCESS;
}
