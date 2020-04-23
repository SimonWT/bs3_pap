#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"


MPI_Comm cannon_comm;
int my_rank;
int w_size;
double *A, *B ,*C;
unsigned int mat_size=0;
int i,k,j,t;
int left,right,up,down;
int dims[2];
int periods[2];

double start=0, av=0;

void init_grid(){

    dims[0]=0; dims[1]=0;
    periods[0]=1; periods[1]=1;
    MPI_Dims_create(w_size,2,dims); //Creates a division of processors in a cartesian grid
    MPI_Cart_create(MPI_COMM_WORLD,2, dims, periods, 1, &cannon_comm);//Makes a new communicator ordered in cartesian grid
    MPI_Cart_shift(cannon_comm,0,1,&left,&right);
    MPI_Cart_shift(cannon_comm,1,1,&up,&down);

    /*
    MPI_Comm_rank(grid_comm, *rank_grid);  // Determines the rank in the new comunicatior
	MPI_Cart_coords(grid_comm, rank_grid, 2, *coord_grid); // Translate the rank to coordinates 

    // Let's subdivide the communicator into subgroups which form lower-dimensional cartesian subgrids, In our case lines and columnes
	MPI_Cart_sub(grid_comm,{0,1}, MPI_Comm *row_comm); // Lines
	MPI_Comm_rank(row_comm,*rank_row);
	MPI_Cart_coords(row_comm, rank_row, 2, *coord_row); 
	
	MPI_Cart_sub(grid_com,{1,0}, MPI_Comm *col_comm); // Columnes
	MPI_Comm_rank(col_comm, *rank_col);
	MPI_Cart_coords(col_comm, rank_col, 2, *coord_col); 
    */

}



// Performs the horizontal shift in A 
void shiftA(){

    MPI_Sendrecv_replace(A,mat_size*mat_size,MPI_DOUBLE,left,0,right,0,cannon_comm,MPI_STATUS_IGNORE); //Sends and receives using a single buffer
	// if (coord_grid[0] >= 1){
	// 	int destination = rank_row - 1;
	// 	int source = (rank_row +1)%sqrt(w_size);
	// 	if (destination < 0){
    //         destination = sqrt(w_size) - 1;
    //     }

	// 	MPI_Sendrecv_replace(&Me.a, 1,MPI_INT,destination,0,source,0,row_comm,MPI_STATUS_IGNORE); //Sends and receives using a single buffer
	// }

	MPI_Barrier(cannon_comm); // Blocks until all processes in the communicator have reached this routine.
    
}

// Performs a Vertical shift in B 
void shiftB(){

    MPI_Sendrecv_replace(B,mat_size,MPI_INT,up,0,down,0,cannon_comm,MPI_STATUS_IGNORE); //Sends and receives using a single buffer
	// if (coord_grid[1] >= 1){
	// 	int destination = rank_col - 1;
	// 	int source = (rank_col +1)%sqrt(w_size);
	// 	if (destination < 0){
    //         destination = sqrt(w_size) - 1;
    //     } 

	// 	MPI_Sendrecv_replace(&Me.b, 1,MPI_INT,destination,0,source,0,col_comm,MPI_STATUS_IGNORE); //Sends and receives using a single buffer
	// }

	MPI_Barrier(cannon_comm); // Blocks until all processes in the communicator have reached this routine.

}

// Depending of the position of the processors inside the grid, it shifts them horizontaly
void preskewingA(){

    for(i = 0; i < mat_size; i++)
	{
        for(j = 0; j < i ; j++){
		    shiftA();
        }
	}

}

// Depending of the position of the processors inside the grid, it shifts them Vertically
void preskewingB(){

    for(i = 0; i < mat_size; i++)
	{
        for(j = 0; j < i ; j++){
		    shiftB();
        }
	}

}


int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // Determines the rank of the calling process in the communicator
    MPI_Comm_size(MPI_COMM_WORLD, &w_size); // Determines the size of the group associated with a communicator

    if(argc != 2){
        printf("usage: %s matrix_size\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    else{
        mat_size = atoi(argv[1]);
    }

    //if(my_rank == 0){

    A = allocMatrix(mat_size);
    B = allocMatrix(mat_size);
    C = allocMatrix(mat_size);

    initMatrix(mat_size, A);
    initMatrix(mat_size, B);
    initMatrixZero(mat_size, C);

    //buf=(double*)malloc(mat_size*mat_size*sizeof(double));

    init_grid();

    start = MPI_Wtime();

    preskewingA();
    preskewingB();

    for(t = 0; t < dims[0]; t++){
        for(i=0;i<mat_size;i++)
            for(k=0;k<mat_size;k++)
                for(j=0;j<mat_size;j++)
                    C[i*mat_size+j] += A[i*mat_size+k] * B[k*mat_size+j];
        shiftA();
        shiftB();
    }

    preskewingA();
    preskewingB();

    printMatrix(mat_size, A);
    printMatrix(mat_size, B);
    printMatrix(mat_size, C);

    int time = MPI_Wtime() - start;
    printf("Time: %.4fs\n",time);
    
    free(A);
    free(B);
    free(C);
    //}

    MPI_Finalize();
    return 0;

}