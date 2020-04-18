#include <stdio.h>
#include "mpi.h"

int main(int argc, char** argv){

    int rank, size;
    int msg[10]; 
    int rbuf[10];    

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   

	if (rank == 0){
		int newmsg[10] = {1, 2, 3, 4, 5,6, 7, 8, 9, 0};
		memcpy(msg, newmsg, sizeof(newmsg));
	}
    MPI_Scatter( msg, 2, MPI_INT, rbuf, 2, MPI_INT, 0, MPI_COMM_WORLD); 

	printf("Hello␣i'm  %d/%d\n", rank, size);
    printf("rank␣%d:␣I␣received\n", rank);
    int i;
	for(i = 0; i < 2; i++){
		printf("%d ",rbuf[i]);
	}
    printf("\n");
    MPI_Finalize();
    return 0;
}

