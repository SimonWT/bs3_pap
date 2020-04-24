#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"

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

    int n = dimension;
    int r = n / dim;    

    double *locA, *locB, *locC;

    locA = allocMatrix(r);
    locB = allocMatrix(r);
    locC = allocMatrix(r);

    initMatrixZero(r, locA);
    initMatrix(r, locB);
    initMatrixZero(r, locC);

    double diagonal[n];

    int coords[2];

    int step;
    int y, w;

    // int a,s =0;
    // for(a =0; a < n; a++)
    //   for(s = 0; s<n; s++)
    //     printf("%f ", A[a][s]);
    // printf("\n");

    // printf("before for loop\n");
    for(step = 0; step < n; step++){
        //calculate diagonal
        if(rank == 0){
            for(y = 0; y < n; y++){
                diagonal[y] = A[y*r + (y + step) % n];
                // printf("choose from %d %d = %f \n", y,(y + step) % n , A[y][(y + step) % n]);
                // printf("%f ", diagonal[y]);
            }
            // printf("\n");
        }
        // printf("after calc diag\n");
        //broadcast diagonal
        MPI_Bcast(diagonal, P, MPI_INT, 0, MPI_COMM_WORLD);
        // printf("after bcast \n");
        //set corresponding value
        MPI_Cart_coords(cart_comm, rank, 2, coords);
        // printf("gets coordinat %d %d \n", coords[0], coords[1]);

        // locA = diagonal[coords[0]];
        double row;
        int u,g;
        for(u=0; u< r; u++){
            row = diagonal[u + coords[0]*r];
            for(g=0; g<r; g++){
                locA[u * r + g] = row;
            }
        }

        // printf("my locA= %f \n", locA);
        //calculate local C
        for(u=0; u< r; u++){
            for(g=0; g<r; g++){
                locC[u*r + g] = locC[u*r + g] + locA[u*r + g] * locB[u*r + g];
            }
        }
        // locC = locC + locA * locB;
        // printf("my locC= %f \n", locC);
        //shift B
        MPI_Sendrecv_replace(&locB, n, MPI_INT, up, 0, down, 0, cart_comm, MPI_STATUS_IGNORE); //Sends and receives using a single buffer
	    MPI_Barrier(cart_comm); // Blocks until all processes in the communicator have reached this routine.

    }


}

int main(int argc, char** argv){

    double *A, *B ,*C;
    unsigned int mat_size = 0;

    if(argc != 2){
        printf("usage: %s matrix_size\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    // else{
    mat_size = atoi(argv[1]);
    // }

    A = allocMatrix(mat_size);
    B = allocMatrix(mat_size);
    C = allocMatrix(mat_size);

    initMatrix(mat_size, A);
    initMatrix(mat_size, B);
    initMatrixZero(mat_size, C);

    int rank, size;

    MPI_Status status;
    MPI_Comm cart_comm;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    foxMatrixMultiplication(rank, size, mat_size, A, B, C);

    // free(A);
    // free(B);
    // free(C);
    MPI_Finalize();
    return 0;
}

