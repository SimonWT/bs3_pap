#include <stdio.h>
#include <math.h>
#include "mpi.h"

void shiftUp(double *Array, int mat_size, int up, int down, MPI_Comm communicator){
    MPI_Sendrecv_replace(Array, mat_size, MPI_INT, up, 0, down, 0, communicator, MPI_STATUS_IGNORE); //Sends and receives using a single buffer
	MPI_Barrier(&communicator); // Blocks until all processes in the communicator have reached this routine.
}

int main(int argc, char** argv){

    int rank, size;
    // int msg[10]; 
    // int rbuf[10];    

    MPI_Status status;
    MPI_Comm cart_comm;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int dims[2] = {2, 2};
    int periods[2] = {1,1};
    int left, right, up, down;

    MPI_Dims_create(size, 2, dims); //Creates a division of processors in a cartesian grid
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);//Makes a new communicator ordered in cartesian grid

    MPI_Cart_shift(cart_comm,0,1,&left,&right); //
    MPI_Cart_shift(cart_comm,1,1,&up,&down); // 
    

	double A[2][2] = {{1, 2},{1, 2}};
    double locA, locB, locC;
    locC = 0;
    double B[2][2] = {{1, 1},{0, 0}};

    double diagonal[2];

    if ( rank %2 == 0){
        locB = 1;
    }else{
        locB = 0;
    }

    int coords[2];



    int step;
    int y, w;
    int n = 2;

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
                diagonal[y] = A[y][(y + step) % n];
                // printf("choose from %d %d = %f \n", y,(y + step) % n , A[y][(y + step) % n]);
                // printf("%f ", diagonal[y]);
            }
            // printf("\n");
        }
        // printf("after calc diag\n");
        //broadcast diagonal
        MPI_Bcast(diagonal, 4, MPI_INT, 0, MPI_COMM_WORLD);
        // printf("after bcast \n");
        //set corresponding value
        MPI_Cart_coords(cart_comm, rank, 2, coords);
        // printf("gets coordinat %d %d \n", coords[0], coords[1]);
        locA = diagonal[coords[0]];
        printf("my locA= %f \n", locA);
        //calculate local C
        locC = locC + locA * locB;
        printf("my locC= %f \n", locC);
        //shift B
        MPI_Sendrecv_replace(&locB, 1, MPI_INT, up, 0, down, 0, cart_comm, MPI_STATUS_IGNORE); //Sends and receives using a single buffer
	    MPI_Barrier(cart_comm); // Blocks until all processes in the communicator have reached this routine.

    }

	printf("Hello␣i'm  %d/%d, (%d, %d)\n", rank, size, coords[0], coords[1]);
    
    // if (rank == 3){
        printf("my  neighbors: %d %d %d %d\n", left, up, right, down);
    // }
    
    // printf("rank␣%d:␣I␣received\n", rank);
    // int i;
	// for(i = 0; i < 2; i++){
	// 	printf("%d ",rbuf[i]);
	// }
    // printf("\n");
    MPI_Finalize();
    return 0;
}

