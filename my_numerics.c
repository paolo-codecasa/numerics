/*---------------------------------------------*\
|  CWR 2023                                     |
|  Numerik-Bibliothek                           |
|  Author  p.codecasa@stud.uni-goettingen.de    |
\*---------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<tgmath.h>
#include<math.h>
#include"my_numerics.h"
#include<gsl/gsl_rng.h> //random generator
#include<lapacke.h>


//git commit -m "Abgabe jetzt fertig" --date="25 May 2023 20:00 +0100"

////////////////////////////////////////////////////////////AUFGABE35

double solve_linear_system(lapack_int n){

    clock_t begin = clock();

    // solving of linear system


    lapack_int m = n;

    FILE *file = fopen("../../matrices.dat", "r");

    // matrix A
    double *A = (double *) malloc(m*m*sizeof(double));
    // vector B
    double *b = (double *) malloc(m*sizeof(double));

    // fill with random values from matrices.dat
    for(int i = 0; i<m; i++){
        for(int j = 0; j<m; j++){
            fscanf(file, "%lf", &A[i * m + j]);
        }
        fscanf(file, "%lf", &b[i]);
    }

    // print A and b


    lapack_int info, lda=m, ldb=1, nrhs=1;

    // lda: The leading dimension of the matrix A
    // ldb: The leading dimension of the array b
    // nrhs: The number of columns of the right hand size
    // LAPACK_ROW_MAJOR: Inserted values in A and b going throw the rows
    // U: Store the matrix A in the upper triangular form.

    lapack_int *ipiv = (lapack_int*) malloc(m*sizeof(lapack_int));

    info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', m, nrhs, A, lda, ipiv, b, ldb);


    // free allocated space
    free(A);
    free(b);
    free(ipiv);

    fclose(file);

    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;


    return time_spent;
}

////////////////////////////////////////////////////////////AUFGABE31

/*
double mc_integrate(gsl_rng* generator, int D, double *R, double integrand(int, double*), int N){
    //double a = gsl_rng_uniform(generator) //zufällige Zahl in [0,1)
    double x[D];
    double sum;

    for (int i=0; i<N; i++){
        for (int j=0; j<D; j++){
            x[j] = gsl_rng_uniform(generator)*fabs(R[j][0]-R[j][1]) + R[j][0]; //Zufallszahl für jede Dimension und auf Intervall Strecken
        }
        sum += integrand(D, *x);
    }
    double V=domain_volume(D, *R);
    return V/N *sum;
}

double domain_volume(int D, double *R){
    double prod;
    for (int i=0; i<D; i++){
        prod = prod * fabs(R[i][0]-R[i][1]); //Formel Volumen Hyper-Würfel
    }
    return prod;
}
*/

double standardDeviation(double x[],int N) {
    double sum = 0.0;
    double meanValue;
    double standDev = 0.0;

    for (int i = 0; i < N; ++i) 
    {
        sum += x[i];
    }

    meanValue = sum / N;
    
    for (int i = 0; i < N; ++i)
    {
        standDev += (x[i] - meanValue) * (x[i] - meanValue);
    }

    standDev = sqrt(standDev / N);

    return standDev;
}

////////////////////////////////////////////////////////////AUFGABE30
//double mn_middle_point_integration(double lim[][], double dx, double integrand(double *)){
    /*double delta_y = (double) (right-left) / (double) N;

double sum = 0.0; 

for (int k=0; k < N; k++){

    double y = left + k * delta_y;
    double f = integrand(y);
    double A = f*delta_y;
    sum += A;
}

return sum;
}*/

double faculty(double n){
    if (n==0) {return 1; }
    else { return n * faculty(n-1); }
}

////////////////////////////////////////////////////////////AUFGABE23

void FTCS_Step(double *y, void *params){
    //printf("FTCS stepping \n");
    
    double* param_array = (double*) params;
    double dt = param_array[0];
    double dy = param_array[1];
    double D = param_array[2];
    double h = param_array[3];
    int N = (int)(h/dy); //Anzahl steps
    

    
    double* Y_out = (double*) malloc(N * sizeof(double)); //double Y_out[N]; //output

    for(unsigned int i = 1 ; i < (unsigned) N-1 ; i++)
    {
        double here = y[i%N];
        double west = y[(i-1)%N];
        double east = y[(i+1)%N];
        //für Dirichlet RB %N hinzufügen (Zugriff auf Array in Speichergrenzen, periodische Randbed erzeugt)

        //CALC NEXT STEP
        Y_out[i] = here + dt * D * (west - 2*here + east) / (dy * dy);

        //if (i%(N/10) < (i-1)%(N/10)){ printf("%d. tenth of FTCS done \n", (int)i/(N/10));}
    }

    //////////////////////////////////////////RandBed
    //Dirichlet
    //("..." ersetzen mit expliziten RandBed)
    /*
       double left = ... ;
       double righLim = ... ;
    */

    //Von Neumann (keine Änderung am Rand: w=0)
    //in Aufgabe23 verwendet
    Y_out[0] = y[1]; //far west RB
    Y_out[N-1] = y[N-2]; //far east RB

    for(size_t j = 0; j< (unsigned) N; j++){
        y[j] = Y_out[j];
    } //update original array

    free(Y_out);
    //printf("FTCS stepped \n");
}

static double mn_gaussian(double x){
    return exp(-x*x);
}

double mn_integrateSquare(double left, double right, int N, double integrand(double)){

double delta_y = (double) (right-left) / (double) N;

double sum = 0.0; 

for (int k=0; k < N; k++){

    double y = left + k * delta_y;
    double f = integrand(y);
    double A = f*delta_y;
    sum += A;
}

return sum;
}

double mn_erf_midpoint(double x, double delta_x){

    double left = 0;
    int N = abs((int) (x / delta_x));



    double sum = 0.0; 

    for (int k=0; k < N; k++){

        double y = left + k * delta_x;
        double f = mn_gaussian(y);
        double A = f*delta_x;
        sum += A;
    }

    return (2.0/sqrt(3.141593))*sum;

}

double mn_erf_simpson(double x, double delta_x){
    int sign = 1;		// Vorzeichen Grenzen
	double left = 0; //linke Integrationsgrenze
    if(x < left){	//falls intervall umgekehrt (unten>oben)
		sign = -1;
		x = -x;
	}
    
    int N = abs((int)(x/delta_x));

    double sum = 0.0; 

    for (int k=0; k < N; k++){
        double y = left + k * delta_x;
        double f =  mn_gaussian(y);
        sum += f;
    }    

    double A = (2/sqrt(M_PI))* delta_x * (0.5 * mn_gaussian(left) + 0.5 * mn_gaussian(x) + sum);

    return A*sign;
}

double mn_integrate_simpson(double left, double right, double dx, double integrand(double, double), double k){


    double sum = 0.0; 

    double N = (right-left)/dx;

    for (int i = 0; i<N; i++){

        double x_i = left + i * dx;
        double x_i_1 = left + (i+1) * dx;
        double m_i = x_i + 0.5 * dx;

        sum += (dx/6) * (integrand(x_i, k) + 4 * integrand(m_i, k) + integrand(x_i_1, k));

    }    

    return sum;

}

double diff(double x, double delta, double func(double), void* params){
    return (func(x + delta) - func(x-delta))/(2*delta);
}

Tuple mn_pq (double a, double b, double c){
    Tuple sol;
    sol.x1 = (-b + sqrt(b*b - 4.0 * a * c))/ (2.0*a); // big number divided by big number
    sol.x2 = (-b - sqrt(b*b - 4.0 * a * c))/ (2.0*a); // big number divided by small number -> high error
return sol;
}

Tuple mn_pq_2 (double a, double b, double c){
    Tuple sol;
    sol.x1 = 2.0*c/(-b - sqrt(b*b - 4.0*a*c)); // big number divided by small number -> high error
    sol.x2 = 2.0*c/(-b + sqrt(b*b - 4.0*a*c)); // big number divided by big number
return sol;
}

Tuple mn_solve_quadratic (double a, double b, double c){
    Tuple sol;
    if(b<=0){
    sol.x1 = (-b + sqrt(b*b - 4.0 * a * c))/ (2.0*a);
    sol.x2 = 2.0*c/(-b + sqrt(b*b - 4.0*a*c));
    }else{
        sol.x1 = 2.0*c/(-b - sqrt(b*b - 4.0*a*c));
        sol.x2 = (-b - sqrt(b*b - 4.0 * a * c))/ (2.0*a);
    }
    return sol;
}

void euler_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params){
    // malloc space for F
    double *f = (double*) malloc(dimension * sizeof(double));
    if (f == NULL ) {
		printf("Speicherallokation fehlgeschlagen in euler_step\n");
		exit(EXIT_FAILURE); // oder `return -1
	}
    // calc derivatives v and accelerations (F/m)
    func(t,y,f,params);
    // update state vector
    for(int i = 0; i<dimension; i++){
        y[i] += delta_t * (*(f+i));
    }
    free(f);
}

void rk4_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params){
    // k vectors
    double *k1 = (double*) malloc(dimension * sizeof(double));
    double *k2 = (double*) malloc(dimension * sizeof(double));
    double *k3 = (double*) malloc(dimension * sizeof(double));
    double *k4 = (double*) malloc(dimension * sizeof(double));
    // support intermediate states
    double *sup = (double*) malloc(dimension * sizeof(double));

    if (k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL || sup == NULL) {
		printf("Speicherallokation fehlgeschlagen in rk4_step.\n");
		exit(EXIT_FAILURE); // oder `return -1;`
    	}


    //  calc k vectors 1 to 4
    func(t,y,k1,params);
    for(int i = 0;i < dimension;i++) {
        sup[i] = y[i] +0.5*delta_t *k1[i];
    }

    func(t+0.5*delta_t,sup,k2,params);
    for(int i = 0;i < dimension;i++) {
        sup[i] = y[i] +0.5*delta_t *k2[i];
    }

   func(t+0.5*delta_t,sup,k3,params);
   for(int i = 0;i < dimension;i++) {
        sup[i] = y[i] +delta_t *k3[i];
    }
    
    func(t+delta_t,sup,k4,params);
    for(int i = 0;i < dimension;i++) {
        y[i] += delta_t/6 *(k1[i]+2*k2[i] +2 *k3[i] + k4[i]); // update state vector (factor delta_t and weighting was missing)
    }

    //free space
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(sup);
}

void rk2_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params){

    //support (in rk4_step wird es y_ genannt)
    double *support = (double*) malloc(dimension * sizeof(double));
    
    double *f = (double*) malloc(dimension * sizeof(double));
    
    double *k1 = (double*) malloc(dimension * sizeof(double));
    double *k2 = (double*) malloc(dimension * sizeof(double));
   
    func(t,y,f,params);
    for(int i = 0; i<dimension; i++){
        k1[i] = delta_t * f[i];
        support[i] = y[i] + k1[i]*0.5;
    }
    func(t + 0.5*delta_t, support, f, params);
    for(int i = 0; i<dimension;i++){
        k2[i] = delta_t*f[i];
        support[i] = y[i] + k2[i]*0.5;
    }

    // update state vector
    for(int i = 0; i<dimension; i++){
        y[i] = y[i] + k2[i];
    }

    // free space
    free(f);
    free(k1);
    free(k2);
    free(support);
  
}

/*double find_root(double func(double,void*), double x0, double delta, const double rel_tol, const int max_iter,void *params) {

    double x = x0;
    int iteration  = 0;

    while(iteration < max_iter) {
        
        double f = func(x,params);
        double df = diff(x,delta,func,params);
        double x_old = x;
        
        x -= f/df;
        
        if(fabs(x_old-x)/fabs(x) < rel_tol) break;
        printf("%g\n", df);
        
        iteration++;

    }

    return x;

}

double fourietransformationsberechnung(double M,double k, double f(double, void*),void* params) {
        double dx = ((double*)params)[0];
        double neueparams[] = {k,((double*)params)[1]};

        double x = integrate_simpson(-M,M,dx,f,neueparams);

        return x;
}*/

