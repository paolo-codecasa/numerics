/*---------------------------------------------*\
|  CWR 2023                                     |
|  Author  p.codecasa@stud.uni-goettingen.de    |
\*---------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "my_numerics.h"
#include<gsl/gsl_rng.h> //random generator
#include<time.h> //seed


double density_func_unitSphere(int D, double *x){
    double sum=0;
    for (int i=0; i<D; i++){ sum += x[i]*x[i];}
    if (sqrt(sum)<=1) { return 1;}
    else { return 0;}
}

double domain_volume(int D, double *R){
    double prod = 1; //not zero, because beginnig value will be multiplied to the next ones
    for (int i=0; i<D*2; i+=2){
        prod = prod * fabs(R[i]-R[i+1]); //Formel Volumen Hyper-Würfel
    }
    return prod;
}

double mc_integrate(gsl_rng* generator, int D, double *R, double integrand(int, double*), int N){
    //double a = gsl_rng_uniform(generator) //zufällige Zahl in [0,1)
    double x[D];
    double sum=0;

    for (int i=0; i<N; i++){
        for (int j=0; j<D; j++){
            
            x[j] = gsl_rng_uniform(generator)*fabs(R[j*2]-R[j*2+1]) + R[j*2]; //Zufallszahl für jede Dimension und auf Intervall Strecken
        }
        sum += integrand(D, x);
    }
    double V=domain_volume(D, R);
    return V/N *sum;
}

int main() {
    //random stuff
    gsl_rng *myGenerator; //pointer to GSL_RNG structure
    myGenerator = gsl_rng_alloc(gsl_rng_randu); //build wanted generator
    gsl_rng_set(myGenerator, time(NULL)); //planting seed
    //double a = gsl_rng_uniform(myGenerator) //random number in [0,1)


    ///////////////////////////////////////////// 
    // Simuliere Volumenintegral mit montecarlo für Hyperkugeln unterschiedlicher Dimensionen und vergleiche sie mit den analytischen Werten

    //params
    int N=100000; //Anzahl ausgewerteten Stellen
    double D_from=1; // Dimension Körper
    double D_to=7;
    int sims=100; //Anzahl simulationen
    FILE* file;

    ///////////////////////////////////////////// numerical
    fprintf(stderr,"Aufgabe 4\n");
    file = fopen("sphere_num.csv", "w"); //open file to write to

    for (int D=D_from; D<=D_to; D++){ // Gehe alle zu simulierenden Dimensionen durch
        //Umgebung des Volumen über den integriert wird
        /////////////////Hyperwürfel um die Hyper-Einheitskugel
        fprintf(stderr, "Durchgang mit D=%d\n", D);
        double R[D*2];
        for (int i=0; i<D*2; i+=2)
        {
            R[i] = -1; //untere Grenze Hyper-Würfel (entspricht der der Hyperkugel)
            R[i+1] = 1; //obere Grenze Hyper-Würfel
        }

        //simuliere so oft wie nötig
        for (int i=0; i<sims; i++){
            double I = mc_integrate(myGenerator, D, R, density_func_unitSphere, N); //simulation, I ist das Volumenintegral
            fprintf(file, "%d, %f \n", D, I);
        }
    }
    fclose(file);

    ///////////////////////////////////////////// analytical
    file = fopen("sphere_ana.csv", "w"); 
    for (double D=D_from; D<D_to+1; D++){ // "für alle benötigten x-Werte"
        double r=1; //Radius hyperkugel, hier kein Array nötig
        
            double gamma = tgamma(D/2.0 +1);
            /*if ((D/2.0 + 1.0) % 1 == 0) { //even input
                double n = D/2.0 + 1;
                gamma = faculty(n-1);
                printf("even");
            }
            else if ((D/2.0 + 1.0) % 1 == 0.5) { //odd input
                double n = D/2.0 + 0.5;
                gamma = faculty(2*n) * sqrt(M_PI) / (faculty(n) * pow(4, n));
                printf("odd");
            } else { printf("is D=%f a natural number?", D); }*/
            double vol = pow(r, D) * (pow(M_PI, (D/2.0)))/gamma;
            fprintf(file, "%f, %f \n", D, vol);
    }
    fclose(file);





    /////////////////////////////////////////////////////
    // Simuliere Integral für Hyperkugel bestimmter Dimension mit mehreren Anzahlen an Stützstellen und vergleiche die Abweichungen zur analytischen Lösung
    int D=5; // Dimension Körper
    double R[D*2];
    for (int i=0; i<D*2; i+=2)
    {
        R[i] = -1; //untere Grenze Hyper-Würfel (entspricht der der Hyperkugel)
        R[i+1] = 1; //obere Grenze Hyper-Würfel
    }

    // analytical
    double r=1; //Radius hyperkugel
    double gamma = tgamma(D/2.0 +1);
    double vol = pow(r, D) * (pow(M_PI, (D/2.0)))/gamma;


    file = fopen("logError_montecarlo.csv", "w");
    
    for(int i = 0; i < 100; i++)
    {
        N = pow(10, 2 + 0.06 * i); //Anzahl ausgewerteten Stellen
        //(0.06 is the length of on of the 100 equally distant logaritmic intervalls between 10^2 and 10^8)

        double I = mc_integrate(myGenerator, D, R, density_func_unitSphere, N); //simulation, I ist das Volumenintegral
        fprintf(file, "%d, %f \n", N, fabs(I - vol));
    }
    fclose(file);

    ///////////////////////////////////////////////////// 
    // Simuliere Integral für Hyperkugel bestimmter Dimension mit mehreren Anzahlen an Stützstellen und unterteile die jeweilige Simulation
    // in mehreren Untersimulationen, um die Standardabweichung zu bestimmen
    fprintf(stderr,"Aufgabe 6\n");

    D=5; // Dimension Körper
    //double r=1; //Radius hyperkugel (wurde schon deklariert)
    int M=100; //Anzahl Untersimulationen jeweils
    double res[M]; //Array mit M mögliche Ergebnisse für die jeweilige Simulation

    // analytical
    gamma = tgamma(D/2.0 +1); // (wurde schon deklariert)
    vol = pow(r, D) * (pow(M_PI, (D/2.0)))/gamma; // (wurde schon deklariert)
    //double R[D*2]; // schon allgemein deklariert
    
    
    file = fopen("standardDeviation_montecarlo.csv", "w");

    for(int i = 0; i < 100; i++)
    {
        N = pow(10, 2 + 0.06 * i); //Anzahl ausgewerteten Stellen
        fprintf(stderr,"Durchgang %d, mit N = %d \n", i, N);
        //(0.06 is the length of on of the 100 equally distant logaritmic intervalls between 10^2 and 10^8)

        for(int j = 0; j < M; j++) {
            double I = mc_integrate(myGenerator, D, R, density_func_unitSphere, N/M); //simulation, I ist das Volumenintegral
            res[j] = I;
        }
        
        double standDev = standardDeviation(res, M);
        
        fprintf(file, "%d, %f \n", N, standDev);

    }
    fclose(file);



    //////////////////////////////////////////////////////////////////////// end
    gsl_rng_free(myGenerator); //free up Space of random numbers generator
}

