#include <stdio.h>
#include <mpi.h>

#include <stdlib.h>

#include "utils.h"

/* a matrix multiplication without locality (column-first)*/
void parallelMatrixMultiplication(int rank, int P, int dimension, double *A, double *B, double *C)
{
    MPI_Status status;

    int r, n;
    n = dimension;
    r = n / P;    
    double tempS[r][n], tempR[r][n], distA[r][n], distC[n][n];

    // Distribute the A rows to process
    MPI_Scatter( A, r * n, MPI_DOUBLE, distA, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    
    //Distribute the C rows to process
    MPI_Scatter( C, r * n, MPI_DOUBLE, distC, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

    // Distribute the B rows to process
    MPI_Scatter( B, r * n, MPI_DOUBLE, tempS, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

    int step = 0;
    for(step = 0; step < P; step++){
        //Send the part of B to next process
        MPI_Send(tempS, r * n, MPI_DOUBLE, (rank+1) % P, rank, MPI_COMM_WORLD);
        //Recive the part of B to previus process
        MPI_Recv(tempR, r * n, MPI_DOUBLE, (rank-1) % P, (rank-1) % P, MPI_COMM_WORLD, &status);
    
        int block = abs((rank - step) % P);

        //locl computation C
        int l = 0;
        int j = 0;
        int k = 0;
        int i = 0;
        for( l = 0; l < P; l++){
            for( i = 0; i < r; i++){
                for( j = 0; j < r; j++ ) {
                    for(k = 0; k < r; k++ ){ 
                        distC[i][ l * r + j] = distC[i][ l * r + j] + distA[i][block * r + k] * tempS[k][l*r+j];
                    }
                }
            }
        }

        int t, y;
        for(t=0; t < r; t++){
            for(y= 0; y < n; y++) {
                tempS[t][y] = tempR[t][y];
            }
        }
    }
    
    MPI_Gather(distC, r * n, MPI_DOUBLE, C, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


int main(int argc, char *argv[])
{
    unsigned int exp ;
    double *A, *B ,*C;
    double *A_check, *B_check ,*C_check;

    unsigned int mat_size=0;

    int my_rank;
    int w_size;

    double start=0, av=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &w_size);

    if(argc != 2){
        printf("usage: %s matrix_size\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    else{
        mat_size = atoi(argv[1]);
    }
    
    if(my_rank == 0){

        printf("test with a matrix of size %u x %u\n", mat_size, mat_size);
        // printf("blyath %d %d ", my_rank, w_size);

        A = allocMatrix(mat_size);
        B = allocMatrix(mat_size);
        C = allocMatrix(mat_size);

        
    }

#ifdef PERF_EVAL
    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        if(my_rank == 0){
            
            initMatrix(mat_size, A);
            initMatrix(mat_size, B);
            initMatrixZero(mat_size, C);
            
            start = MPI_Wtime();
            
            sequentialMatrixMultiplication_REF(mat_size, A, B , C);
            
            experiments [exp] = MPI_Wtime() - start;
        }
    }

    if(my_rank == 0){
        av = average_time() ;  
        
        printf ("\n REF sequential time \t\t\t %.3lf seconds\n\n", av) ;
    }
    
    
    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        if(my_rank == 0){
            initMatrix(mat_size, A);
            initMatrix(mat_size, B);
            initMatrixZero(mat_size, C);
            
            start = MPI_Wtime();
            parallelMatrixMultiplication(my_rank, w_size, mat_size, A, B , C);
            
            experiments [exp] = MPI_Wtime() - start;
        }
        
    }

    if(my_rank == 0){
        av = average_time() ;  
        
        printf ("\n my mat_mult \t\t\t %.3lf seconds\n\n", av) ;
    }
    

#endif /*PERF_EVAL*/

#ifdef CHECK_CORRECTNESS
    /* running my sequential implementation of the matrix
       multiplication */
    if(my_rank == 0){
        initMatrix(mat_size, A);
        initMatrix(mat_size, B);
        initMatrixZero(mat_size, C);

        A_check = createMatrixCopy(mat_size, A);
        B_check = createMatrixCopy(mat_size, B);
        C_check = allocMatrix(mat_size);

        initMatrixZero(mat_size, C_check);
    }

    /* check for correctness */
    // if(my_rank == 0){
    parallelMatrixMultiplication(my_rank, w_size, mat_size, A, B , C);    
    
    if(my_rank == 0){
        sequentialMatrixMultiplication_REF(mat_size, A_check, B_check , C_check);

        if(checkMatricesEquality(mat_size, C, C_check)){
            printf("\t CORRECT matrix multiplication result \n");
        }
        else{
            printf("\t FAILED matrix multiplication !!! \n");
        }

        /* printMatrix(mat_size, C); */
        /* printMatrix(mat_size, C_check); */
        
        free(A_check);
        free(B_check);
        free(C_check);
    }

#endif /* CHECK_CORRECTNESS */

    if(my_rank == 0){
        free(A);
        free(B);
        free(C);
    }

    MPI_Finalize();

    return 0;
}
