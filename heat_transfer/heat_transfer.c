/*---------------------------------------------*\
|  CWR 2023                                     |
|  Author  p.codecasa@stud.uni-goettingen.de    |
\*---------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../../my_numerics.h"

// Physikalische Konstanten

double L = 100.0;
double delta_x = 0.1;
double D = 0.1;

double T0_L = 1.0;
double T0_R = - 1.0;

size_t N = 100;

// wärmeleitungsgleichung ftcs

void heatFTCS(double t, double y[], double dt)
{
    double* y_new = (double*) malloc(N * sizeof(double));
    y_new[0] = T0_L; //Anfangstemperatur links
    y_new[N-1] = T0_R; //Anfangstemperatur rechts

    for(size_t i = 0; i < N; ++i) //size_t damit modulo pos; Achte auf Ränder
    {
        //Modulo: Zugriffe auf Array in Grenzen --> periodische RandBed
        double T = y[i % N];
        double T_W = y[(i - 1) % N];
        double T_O = y[(i + 1) % N];

        //nur für periodisch
        /*y_new[0] = T0_L;
        y_new[N-1] = T0_R;*/

        y_new[i] = T + dt * D * (T_W - 2 * T + T_O) / (delta_x);
    }

    for(size_t i = 0; i<N; i++) y[i] = y_new[i];
    
    free(y_new);
}

void write_to_file(double t, double y[], FILE* file){

    fprintf(file,"%g", t);
    for(size_t i = 0; i<N; i++){
        fprintf(file, ",%g", y[i]);
    }
    fprintf(file,"\n");
}

int main() {
    // Initialisierung
    double y[N];
    for(size_t i = 0; i<N; i++){
        y[i] = 0;
    }
    y[0] = T0_L;
    y[N-1] = T0_R;

    // Datei
    FILE* file = fopen("temp.csv", "w");
    write_to_file(0,y,file);

    // Zeitschritt       
    double delta_t = delta_x/(2.0*D);

    // Hauptsimulation
    double t = delta_t;
    while(t<100.0){

        printf("Generating %f out of %f\n.", t, 2.0);

        t+=delta_t;
        heatFTCS(t, y, delta_t);
        write_to_file(t,y,file);
        
    }

    printf("Finished\n");

    fclose(file);
}
