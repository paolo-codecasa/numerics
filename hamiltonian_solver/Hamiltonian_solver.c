/*---------------------------------------------*\
|  CWR 2023                                     |
|  Autor  p.codecasa@stud.uni-goettingen.de    |
\*---------------------------------------------*/

// imports
#include<stdio.h>
#include<math.h>
#include<lapacke.h>
#include<time.h>

//Konstanten
double h = 1; //Plank'sches Wirkungsquantum
double m = 1; //Masse
double k_spring = 1; //Federkonstante

const double L = 10; //Intervalllänge
int N = 100; //Anzahl Schritte
int n=20;
    

double noPot(double x)
{
    return 0 * x;
}

double harmPot(double x)
{
    return (k_spring / 2.0) * x * x; //vorzeichen?
}

void init_hamiltonian(double *H, int N, double dx, double V(double))
{
    printf("H initialized\n");
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            int index = j*N + i;
            H[index]=0;
            if (i==j) { H[index] = (1/(dx*dx)) + V (i*dx - L/2.0); }
            else if (i==j+1 || i==j-1) { H[index] = -(1/(2*dx*dx)); }
        }
    }

    int j=0;
    for (int i=0; i<N*N; i++)
    {
        if (j>=N)
        {
            printf("\n");
            j=0;
        }
        printf("%g ", H[i]);
        j++;
    }
    printf("\n");
}

lapack_int lapacke_diagonalize(double *A, double *W, int N)
{
    /*
    LAPACKE_dsyev()
    lapack_int LAPACKE_dsyev
    ( 	int  	matrix_layout,      LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR
		char  	jobz,               'N': eigenvalues only; 'V': eigenvalues and eigenvectors
		char  	uplo,               'U':  Upper triangle of A is stored; 'L': Lower triangle of A is stored.
		lapack_int  	n,          The order of the matrix A.  N >= 0.
		double *  	a,              
		lapack_int  	lda,        The leading dimension of the array A.  LDA >= max(1,N).
		double *  	w 
	) 		
    */
    lapack_int info=LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, A, N, W);
    int j=0;
    printf("H: \n");
    for (int i=0; i<N*N; i++)
    {
        if (j>=N)
        {
            printf("\n");
            j=0;
        }
        printf("%g ", A[i]);
        j++;
    }
    printf("\nW: \n");
    for (int i=0; i<N; i++)
    {
        printf("%g ", W[i]);
    }
    printf("\n");
    return info;
}

int main(){
    double dx = (double) L/N; //Schrittweite

    double* H = (double *) malloc(N * N * sizeof(double));
    double* W = (double *) malloc(N * sizeof(double));

    ///////////////////////////////////////////////// 

    init_hamiltonian(H, N, dx, noPot); //Initialisierung kein Potential
    double omega=H[0], t_h=H[1];
    lapacke_diagonalize(H, W, N);



    FILE* myFile1 = fopen("numerisch.csv", "w");
    fprintf(myFile1,"# n, lambda\n");
    for(int i=0; i<n; i++)
    {
        fprintf(myFile1,"%d,%.16g\n", i+1, W[i]);
    }
    fclose(myFile1);

    FILE* myFile2 = fopen("analytisch.csv", "w");
    fprintf(myFile2,"# n, EigenWert\n");
    for(int i=1; i<n+1; i++)
    {
        //int index = N-i; //Reihenfolge Output lapacke aufsteigend, analytisch absteigend
        double EW = omega + 2*t_h*cos(i*M_PI/(N+1));
        fprintf(myFile2,"%d,%.16g\n", i, EW);
    }
    fclose(myFile2);

    FILE* myFile3 = fopen("EigenVektoren.csv", "w");
    fprintf(myFile3,"# x, psi(x, n=1), psi(x, n=2), psi(x, n=3)\n");
    for(int k=0; k<N; k++){
        fprintf(myFile3,"%g, %g, %g, %g\n", k*dx-L/2.0, H[k], H[k+N*1], H[k+N*2]);
    }
    fclose(myFile3);


    ///////////////////////////////////////////////// 

    init_hamiltonian(H, N, dx, harmPot); //Initialisierung kein Potential
    //for (int i=0; i<N; i++) W[i]=0;
    lapacke_diagonalize(H, W, N);

    FILE* myFile4 = fopen("numerisch_HARM.csv", "w");
    fprintf(myFile4,"# n, lambda\n");
    for(int i=0; i<n; i++)
    {
        fprintf(myFile4,"%d,%.16g\n", i+1, W[i]);
    }
    fclose(myFile4);

    FILE* myFile5 = fopen("analytisch_HARM.csv", "w");
    fprintf(myFile5,"# n, EigenWert\n");
    omega = 1;
    for(int i=0; i<n+1; i++) //***************** i=0 statt 1
    {
        //int index = N-i; //Reihenfolge Output lapacke aufsteigend, analytisch absteigend
        double EW = h * omega * (i + 0.5);
        fprintf(myFile5,"%d,%.16g\n", i, EW);
    }
    fclose(myFile5);

    free(H);
    free(W);
}
