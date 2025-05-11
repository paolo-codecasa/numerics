#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <math.h>
#include "../../my_numerics.h"

int main(void)
{
    double N = 83 * pow(10, 6); //gesamtpopulation (überprüfe ob erhalten)
    double E_0 = 30000; //exposed
    double J_0 = 9000; //infectious
    double R_0 = 0; //removed
    double S_0 = N - (E_0 + J_0 + R_0);
    double alpha = 1/6.1; //Kehrwert präinfektiöser zeit [tage^(-1)]
    double gamma = 1/3.7; //Kehrwert infektiöser Zeit [tage^(-1)]
    double r_0[] = {1.25, 1.5, 2}; //Reproduktionsrate
    int n=sizeof(r_0)/sizeof(double); //anzahl an getesteten Reproduktionsraten
    double beta[n];
    for(int i=0; i<n; i++) beta[i]= r_0[i] * gamma;
    double t_stop = 100; //tage
    //double t_red[3] = {15, 20, 35}; //Zeitpunkt der Reduktionsmaßnahmen [tage]
    //double r_red = 0.9; //Reduktionsrate nach Maßnahmen

    FILE* myFile = fopen("dataEpidemy.csv", "w");
    fprintf(myFile, "Zeit, Susceptible, Exposed, Infectious, Removed"); //alle abhängig von beta -->13 spalten
    double S[n]; //susceptible
    double E[n]; //exposed
    double J[n]; //infectious
    double R[n]; //removed
    for(int i=0; i<n; i++)
    {
        S[i]= S_0;
        E[i]= E_0;
        J[i]= J_0;
        R[i]= R_0;
    }
    int samples = 100; //samples pro tag
    int dt = 1/samples;
    for (int i=0; i < t_stop * samples; i++)
    {
        for (int j=0; j<3; j++)
        {
            //Reihenfolge ist wichtig !!!!!!!
            S[j] = -beta[j] * J[j] * S[j] * dt / N + S[j];
            E[j] = beta[j] * J[j] * S[j] * dt / N - alpha*E[j] + E[j]; //E ist bei dieser implementation S immer ein schritt dt voraus, ist aber egal
            J[j] = alpha * E[j] * dt -gamma * J[j] * dt + J[j];
            R[j] = gamma * J[j] * dt + R[j];
        }
        fprintf(myFile, "%d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", i, S[0], S[1], S[2], E[0], E[1], E[2], J[0], J[1], J[2], R[0], R[1], R[2]);
    }

    fclose(myFile);

    return EXIT_SUCCESS;
}
