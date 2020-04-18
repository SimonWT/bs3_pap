#include <stdio.h>
#include "mpi.h"

int main(int argc, char** argv){

    int rank, size;
    float msg[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};     

    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   

	if (rank == 0) {
	    printf("Hello␣i'm sender! %d/%d\n", rank, size);
	    MPI_Send(msg, 10, MPI_FLOAT, 1, 99, MPI_COMM_WORLD);
	}
	else {
	    MPI_Recv(msg, 10, MPI_FLOAT, 0, 99, MPI_COMM_WORLD, &status);
            printf("I␣received␣%s!\n", msg);
        }

    MPI_Finalize();
    return 0;
}

