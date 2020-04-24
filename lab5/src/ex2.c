#include <stdio.h>
#include <mpi.h>

#include <stdlib.h>

#include "utils.h"

/* a matrix multiplication without locality (column-first)*/
void foxMatrixMultiplication(int rank, int P, int dimension, double *A, double *B, double *C){

    int dim = (int) sqrt(P);
    int dims[2] = { dim, dim };
    int periods[2] = { 1, 1 };
    int left, right, up, down;

    MPI_Status status;
    MPI_Comm cart_comm;

    MPI_Dims_create(P, 2, dims); //Creates a division of processors in a cartesian grid
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);//Makes a new communicator ordered in cartesian grid

    MPI_Cart_shift(cart_comm, 0, 1, &left, &right);
    MPI_Cart_shift(cart_comm, 1, 1, &up, &down); 

    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    int n = dimension;
    int r = n / dim;    

    double *locA, *locB, *locC;

    locA = allocMatrix(r);
    locB = allocMatrix(r);
    locC = allocMatrix(r);

    initMatrixZero(r, locA);
    initMatrix(r, locB);
    initMatrixZero(r, locC);

    // int a,b;
    // for (a = 0; a < r; a++){
    //     for(b=0; b<r; b++){
    //         locB[a * r + b] = B[(r*coords[0] + a) * r + (r * coords[1] + b)];
    //     }
    // }

    double diagonal[n];

    // printf("%d. locB:\n",rank);
    // printMatrix(r, locB);

    int step;
    int y, w;
    int u,g;
    // int s =0;
    // for(a =0; a < n; a++)
    //   for(s = 0; s< n; s++)
    //     printf("%f ", A[a * n + s]);
    // printf("\n");

    // printf("before for loop\n");
    for(step = 0; step < n; step++){
        //calculate diagonal
        if(rank == 0){
            // printf("%d. diag:\n",rank);
            for(y = 0; y < n; y++){
                diagonal[y] = A[y*r + (y + step) % n];
                // printf("choose from %d %d = %f \n", y,(y + step) % n , A[y][(y + step) % n]);
                // printf("%f ", diagonal[y]);
            }
            // printf("\n");
            
            // printMatrix(r, locB);
        }

        // printf("after calc diag\n");
        //broadcast diagonal
        MPI_Bcast(diagonal, P * n, MPI_INT, 0, MPI_COMM_WORLD);
        // printf("after bcast \n");
        //set corresponding value
        
        // printf("gets coordinat %d %d \n", coords[0], coords[1]);

        // locA = diagonal[coords[0]];
        double row;
        for(u=0; u< r; u++){
            row = diagonal[u + coords[0]*r];
            // printf("%d diagonal:\t", rank);
            // printf("%f ", row);
            for(g=0; g<r; g++){
                locA[u * r + g] = row;
            }
        }
        printf("\n");
        printf("%d my locB\n", rank);
        printMatrix(r, locB);
        //calculate local C
        for(u=0; u< r; u++){
            for(g=0; g<r; g++){
                locC[u * r + g] = locC[u*r + g] + locA[u*r + g] * locB[u*r + g];
            }
        }
        // locC = locC + locA * locB;
        double row1[r];
        int o;
        for(o = 0; o < r; o++)
            row1[o] = locB[o];
        // printf("my locC= %f \n", locC);
        //shift B
        MPI_Sendrecv_replace(&row1, r, MPI_INT, up, 0, down, 0, cart_comm, MPI_STATUS_IGNORE); //Sends and receives using a single buffer
        
        for(u=0; u< r; u++){
            for(g=0; g < r; g++){
                if(u == r - 1)
                     locB[u * r + g] = row1[g];
                else
                    locB[u * r + g] = locB[(u + 1)* r + g ];
            }
        }

	    MPI_Barrier(cart_comm); // Blocks until all processes in the communicator have reached this routine.

    }

    printf("%d . my locC: \n", rank);
    printMatrix(r, locC);


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
            foxMatrixMultiplication(my_rank, w_size, mat_size, A, B , C);
            
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
// if(my_rank == 0){
        printMatrix(mat_size, A);
        printMatrix(mat_size, B);
// }
        initMatrixZero(mat_size, C_check);
    }

    /* check for correctness */
    // if(my_rank == 0){
    foxMatrixMultiplication(my_rank, w_size, mat_size, A, B , C);    
    
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
