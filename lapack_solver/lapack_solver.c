/*------------------------------------------------------------------------------*\
|  CWR 2023                                                                      |
|  Nutzung der LAPACK-Bibliothek                         |
|  Author  p.codecasa@stud.uni-goettingen.de                                     |
|          mit g.hunder@stud.uni-goettingen.de                                  |
\*------------------------------------------------------------------------------*/


// imports
#include<stdio.h>
#include<math.h>
#include<lapacke.h>
#include<time.h>



void print_matrix_rowmajor(char *desc, lapack_int m, lapack_int n, double *mat,
                           lapack_int ldm)	
{
  lapack_int i, j;
  printf( "\n %s\n", desc ); 
  for( i = 0; i < m; i++ ) 
  {
    for( j = 0; j < n; j++ ) 
      printf( " %6.2f", mat[i*ldm+j] );
    printf( "\n" );
  }
}

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


int main(){


    FILE* cpu_time_solve_file = fopen("cpu_time_solve.txt", "w");
    fprintf(cpu_time_solve_file, "# n, cpu-time");
    for(int i = 0; i<500; i+=50){

        double cpu_time = solve_linear_system(i);
        fprintf(cpu_time_solve_file,"%d,%g\n", i, cpu_time);

    }
    fclose(cpu_time_solve_file);

    
}

