#ifndef MYNUMERICS_H
#define MYNUMERICS_H
#include<stdio.h>
#include<math.h>
#include<time.h>

/*double mc_integrate(gsl_rng* generator, int D, double *R, double integrand(int, double*), int N);

double domain_volume(int D, double *R);
*/

double standardDeviation(double x[],int N);

//double mn_middle_point_integration(double lim[][], double dx, double integrand(double *))

double faculty(double n);

void FTCS_Step(double *y, void *params);

/* Implementations of the error function */
double mn_erf_midpoint(double x, double delta_x);
double mn_erf_simpson(double x, double delta_x);

/* Numerical Integration of scalar 1D function */
double mn_integrateSquare(double left, double right, int N, double integrand(double));

/* Numerical Integration of scalar 2D function with Simpson rule */
double mn_integrate_simpson(double left, double right, double dx, double integrand(double, double), double k);


/* Numerical derivation*/
double diff(double x, double delta, double func(double), void* params);

/* Tuple */

typedef struct {
    double x1;
    double x2;
}Tuple;

//... do not use because of bad accuracy:///////////////////////// 
Tuple mn_pq(double a, double b, double c);                      //
                                                                //
Tuple mn_pq_2 (double a, double b, double c);                   //
//...                                                           //
//////////////////////////////////////////////////////////////////      

/* quadratic solver */
Tuple mn_solve_quadratic(double a, double b, double c);        


/* ODE */
typedef
int ode_func(double, const double[], double[], void *);

/***** ODE SOLVERS *****/
/* Euler step */
void euler_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params);

/* Runge Kutta 4 */
void rk4_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params);

/* Runge Kutta 2 */
void rk2_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params);

//double find_root(double func(double,void*), double x0, double delta, const double rel_tol, const int max_iter,void *params);

//double fourietransformationsberechnung(double M,double k, double f(double, void*),void* params);

#endif


