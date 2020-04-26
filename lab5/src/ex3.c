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
    initMatrixZero(r, locC);

    //some strange part
    int starts[2] = { 0,0 };
    int localSize[2] = { r, r };
    int globalSize[2] = { n, n };
    MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(2, globalSize, localSize, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, r * sizeof(double), &subarrtype);
	MPI_Type_commit(&subarrtype);

    // Scatter the array to all processors
	int* sendCounts = (int*)malloc(sizeof(int) * P);
	int* displacements = (int*)malloc(sizeof(int) * P);

    int i, j = 0;
	if (rank == 0) {
		for (i = 0; i < P; i++) {
			sendCounts[i] = 1;
		}
		int disp = 0;
        // printf("DIPLACEMENT:\n");
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				displacements[i * dim + j] = disp;
                disp += 1;
			}
            disp += (r - 1) * dim;
            // printf("\n");
		}
	}

    //distribute B
    MPI_Scatterv(B, sendCounts, displacements, subarrtype, locB,
		n * n / (P), MPI_DOUBLE,
		0, MPI_COMM_WORLD);

    // printf("(%d). Local B: \n", rank);
    // printMatrix(r, locB);

    double diagonal[n];

    int step;
    int y, w;
    int u,g;

    // printf("before for loop\n");
    for(step = 0; step < n; step++){
        //calculate diagonal
        if(rank == 0){
            for(y = 0; y < n; y++){
                diagonal[y] = A[y * n + (y + step) % n];
                // printf("(%d,%d)", y, (y + step) % n);
                // printf(" %.1lf", diagonal[y]);
            }
            // printf("\n");
        }

       
        //broadcast diagonal
        MPI_Bcast(diagonal, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        //set corresponding value to each row
        double row;
        for(u=0; u< r; u++){
            row = diagonal[u + coords[0] * r];
            // printf("(%d). row = %.1f\n", rank, row);
            for(g=0; g<r; g++){
                locA[u * r + g] = row;
            }
        }
        // printf("%d. locA:\n", rank);
        // printMatrix(r, locA);

        //calculate local C
        for(u=0; u< r; u++){
            for(g=0; g<r; g++){
                locC[u * r + g] = locC[u * r + g] + locA[u * r + g] * locB[u * r + g];
            }
        }
        // printf("%d. locC:\n", rank);
        // printMatrix(r, locC);
        
        //allocate row which be shift to upper row of processes
        double row1[r];
        int o;
        for(o = 0; o < r; o++)
            row1[o] = locB[o];

        //shift B
        MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
        MPI_Sendrecv_replace(&(row1[0]), r, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE); //Sends and receives using a single buffer
        
        //set income row and shift old
        for(u=0; u< r; u++){
            for(g=0; g < r; g++){
                if(u == r - 1)
                     locB[u * r + g] = row1[g];
                else
                    locB[u * r + g] = locB[(u + 1)* r + g ];
            }
        }

	     // Blocks until all processes in the communicator have reached this routine.
         MPI_Barrier(cart_comm);
    }
    

    //Gather the local C to C in root process
    // MPI_Gather(locC, r * r, MPI_DOUBLE, C, r * r, MPI_DOUBLE, 0, cart_comm);
    // MPI_Scatterv(locC, sendCounts, displacements, subarrtype, locB,
	// 	n * n / (P), MPI_DOUBLE,
	// 	0, MPI_COMM_WORLD);
    MPI_Gatherv(locC, n * n / (P), MPI_DOUBLE,
		C, sendCounts, displacements, subarrtype,
		0, MPI_COMM_WORLD);
    
    //Print result
    // if(rank == 0 ){
    //     printf("FINAL C:\t");
    //     printMatrix(n, C);
    // }

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

        printMatrix(mat_size, A);
        printMatrix(mat_size, B);

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
            printMatrix(mat_size, C_check);
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
