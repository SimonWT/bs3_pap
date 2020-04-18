#include <stdio.h>
#include "mpi.h"

int main(int argc, char** argv){

    int rank, size;
    int msg[10];     

    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
	int pid = getpid();
        int localsum[1] = {0};
	int globalsum[1] = {0};
	
	printf("Hello␣i'm  %d/%d, pid: %d\n", rank, size, pid);
	localsum[0] = localsum[0] + pid;
	MPI_Reduce(localsum, globalsum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	//if(rank == 0)
	  MPI_Bcast(globalsum, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("rank␣%d:␣I␣received sum: %d\n", rank, globalsum[0]);
       
    MPI_Finalize();
    return 0;
}

