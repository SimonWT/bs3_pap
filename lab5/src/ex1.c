#include <stdio.h>
#include <mpi.h>

#include <stdlib.h>

#include "utils.h"

/* a matrix multiplication without locality (column-first)*/
void parallelMatrixMultiplication(int rank, int P, int dimension, double *A, double *B, double *C)
{
 
    // int rank = my_rank(); 
    // int P = num_procs();
    MPI_Status status;
    
    if( rank == 0 ){
        printMatrix(dimension, A);

        printMatrix(dimension, B);

    }

    
    printf("%d %d \n", rank, P);

    int r, n;
    n = dimension;
    r = n / P;    
    double tempS[r][n], tempR[r][n], distA[r][n], distC[n][n];

    // tempS = B; /* copy of the values */ 
    // memcpy(B, tempS, dimension * dimension)
    printf("rank %d \n", rank);
    // Distribute the A rows to process
    // if(rank == 0){
        MPI_Scatter( A, r * n, MPI_DOUBLE, distA, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        // printf("A");
        MPI_Scatter( C, r * n, MPI_DOUBLE, distC, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

        // Distribute the B rows to process
        MPI_Scatter( B, r * n, MPI_DOUBLE, tempS, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    // }
    
    // MPI_Send(tempS, r * n, MPI_DOUBLE, (rank+1) % P, rank, MPI_COMM_WORLD);
    // MPI_Recv(tempR, r * n, MPI_DOUBLE, (rank-1) % P, (rank-1) % P, MPI_COMM_WORLD, &status);

    int step = 0;
    for(step = 0; step < P; step++){
        printf("iteration: %d \n\n", step);
        MPI_Send(tempS, r * n, MPI_DOUBLE, (rank+1) % P, rank, MPI_COMM_WORLD);
        printf(" temps S: \n");
        int t, y;
        for(t= 0; t < r; t++){
            for(y= 0; y < n; y++) {
                printf("\t %.1lf",tempS[t][y]);
            }
            printf("\n");
        }

        MPI_Recv(tempR, r * n, MPI_DOUBLE, (rank-1) % P, (rank-1) % P, MPI_COMM_WORLD, &status);
        printf(" temps R: \n");
        for(t= 0; t < r; t++){
            for(y= 0; y < n; y++) {
                printf("\t %.1lf",tempR[t][y]);
            }
            printf("\n");
        }
        // Recv(tempR, (rank−1) % P)
        int block = abs((rank - step) % P);

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
        printf("step distC, block:%d \n", block);
        //memcpy(tempS, tempR, r * n);
        for(t= 0; t < r; t++){
            for(y= 0; y < n; y++) {
                printf("\t %.1lf",distC[t][y]);
                tempS[t][y] = tempR[t][y];
            }
            printf("\n");
        }
    }
    // MPI_( A, r * n, MPI_DOUBLE, distA, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    printf("\n");
    printf(" dist C: \n");
    int t, y;
    for(t= 0; t < r; t++){
        for(y= 0; y < n; y++) {
            printf("\t %.1lf",distC[t][y]);
        }
        printf("\n");
    }
    // printMatrix(dimension, distC);
    MPI_Gather(distC, r * n, MPI_DOUBLE, C, r * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // printMatrix(dimension, C);
    if(rank == 0){
        printf("FINAL:  ");
        printMatrix(dimension, C);
    }
    // }
    printf("end\n");
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
            // print("suka")
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
