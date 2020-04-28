#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"


MPI_Comm cannon_comm;
int my_rank, cannon_rank;
int mycoords[2];
int w_size;
double *A, *B ,*C;
double *locA, *locB, *locC;
unsigned int mat_size=0;
unsigned int local_size=0;
int i,k,j;
int left,right,up,down;
int dims[2], periods[2];
double buf[1];
int shifts;


double start=0, av=0;

void init_grid(){
    dims[0]=0; dims[1]=0;
    periods[0]=1; periods[1]=1;
    MPI_Dims_create(w_size,2,dims); //Creates a division of processors in a cartesian grid
    MPI_Cart_create(MPI_COMM_WORLD,2, dims, periods, 1, &cannon_comm);//Makes a new communicator ordered in cartesian grid
    MPI_Comm_rank(cannon_comm, &cannon_rank);
    MPI_Cart_coords(cannon_comm, cannon_rank, 2,  mycoords);
}

// Performs the horizontal shift in A 
void shiftA(){
    MPI_Cart_shift(cannon_comm, 1, 1, &left, &right);
    MPI_Sendrecv_replace(locA, local_size * local_size, MPI_DOUBLE, left, 1, right, 1, cannon_comm, MPI_STATUS_IGNORE);
}

// Performs a Vertical shift in B 
void shiftB(){
    MPI_Cart_shift(cannon_comm, 0, 1, &up, &down);
    MPI_Sendrecv_replace(locB, local_size * local_size, MPI_DOUBLE, up, 1, down, 1, cannon_comm, MPI_STATUS_IGNORE);
}

// Depending of the position of the processors inside the grid, it shifts them horizontaly
void preskewingA(){

    for (i = 0; i < local_size; i++){

        shifts = mycoords[0] * local_size + i;

        for(j = 0; j < shifts; j++){
            buf[0] = locA[i*local_size];
            MPI_Cart_shift(cannon_comm, 1, 1, &left, &right);
            MPI_Sendrecv_replace(&buf, 1, MPI_DOUBLE, left, 0, right, 0, cannon_comm, MPI_STATUS_IGNORE);
            for( k = 0; k < local_size; k++){
                locA[i * local_size + k] = locA[i *  local_size + k  + 1];
            } 

            locA[local_size * (i+1) - 1] = buf[0];
           
        }
    }
}

// Depending of the position of the processors inside the grid, it shifts them Vertically
void preskewingB(){

    for (i = 0; i < local_size; i++){

        shifts = mycoords[1] * local_size + i;
        
        for(j = 0; j < shifts; j++){
            buf[0] = locB[i];
            MPI_Cart_shift(cannon_comm, 0, 1, &up, &down);
            MPI_Sendrecv_replace(&buf, 1, MPI_DOUBLE, up, 1, down, 1, cannon_comm, MPI_STATUS_IGNORE);
            for( k = 0; k < local_size-1; k++){
                locB[k *local_size + i] = locB[(k+1) *  local_size + i];
            } 

            locB[local_size * (local_size-1) + i] = buf[0];
           
        }
    }
}

// Depending of the position of the processors inside the grid, it shifts them horizontaly
void postskewingA(){
   
    for (i = 0; i < local_size; i++){

        shifts = mycoords[0] * local_size + i;
        shifts = (mat_size - shifts)%mat_size;

        for(j = 0; j < shifts; j++){
            buf[0] = locA[i*local_size];
            MPI_Cart_shift(cannon_comm, 1, 1, &left, &right);
            MPI_Sendrecv_replace(&buf, 1, MPI_DOUBLE, left, 0, right, 0, cannon_comm, MPI_STATUS_IGNORE);
            for( k = 0; k < local_size; k++){
                locA[i * local_size + k] = locA[i *  local_size + k  + 1];
            } 

            locA[local_size * (i+1) - 1] = buf[0];
           
        }
    }
}

// Depending of the position of the processors inside the grid, it shifts them Vertically
void postskewingB(){

    for (i = 0; i < local_size; i++){

        shifts = mycoords[1] * local_size + i;
        shifts = (mat_size - shifts)%mat_size;

        for(j = 0; j < shifts; j++){
            buf[0] = locB[i];
            MPI_Cart_shift(cannon_comm, 0, 1, &up, &down);
            MPI_Sendrecv_replace(&buf, 1, MPI_DOUBLE, up, 1, down, 1, cannon_comm, MPI_STATUS_IGNORE);
            for( k = 0; k < local_size-1; k++){
                locB[k *local_size + i] = locB[(k+1) *  local_size + i];
            } 

            locB[local_size * (local_size-1) + i] = buf[0];
           
        }
    }    
}


int main(int argc, char *argv[]){
    
    MPI_Init(&argc, &argv); // Initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // Determines the rank of the calling process in the communicator
    MPI_Comm_size(MPI_COMM_WORLD, &w_size); // Determines the w_size of the group associated with a communicator
     
    // Gather the matrix size
    if(argc != 2){
        printf("usage: %s matrix_size\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    else{
        mat_size = atoi(argv[1]);
    }

    // Initialize the matrices in rank 0
    if (my_rank == 0){
        A = allocMatrix(mat_size);
        B = allocMatrix(mat_size);
        C = allocMatrix(mat_size);

        initMatrix(mat_size, A);
        initMatrix(mat_size, B);
        initMatrixZero(mat_size, C);
    }
   
    start = MPI_Wtime();

    // Initialize the cartesian grid and creates the new communicator
    init_grid();

    // Computes the local_size of each processor
    local_size = mat_size/dims[0]; 

    // Allocate and initialize the block matrices 
    locA = allocMatrix(local_size);
    locB = allocMatrix(local_size);
    locC = allocMatrix(local_size);

    initMatrixZero(local_size, locA);
    initMatrixZero(local_size, locB);
    initMatrixZero(local_size, locC);

    // Create datatype to describe the submatrices of the global matrix
	int globalSize[2] = { mat_size, mat_size };
	int localSize[2] = { local_size, local_size };
	int starts[2] = { 0,0 };
	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(2, globalSize, localSize, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, local_size * sizeof(double), &subarrtype);
	MPI_Type_commit(&subarrtype);

    // Scatter the array to all processors
	int* sendCounts = (int*)malloc(sizeof(int) * w_size);
	int* displacements = (int*)malloc(sizeof(int) * w_size);

	if (my_rank == 0) {
        int i,j;
		for (i = 0; i < w_size; i++) {
			sendCounts[i] = 1;
		}
		int disp = 0;
		for (i = 0; i < dims[0]; i++) {
			for (j = 0; j < dims[0]; j++) {
				displacements[i * dims[0] + j] = disp;
				disp += 1;
			}
			disp += (local_size-1) * dims[0];
		}

	}

    MPI_Scatterv(A, sendCounts, displacements, subarrtype, locA, mat_size * mat_size / (w_size), MPI_DOUBLE,0, MPI_COMM_WORLD);
	MPI_Scatterv(B, sendCounts, displacements, subarrtype, locB, mat_size * mat_size / (w_size), MPI_DOUBLE,0, MPI_COMM_WORLD);
   
    // Perform Preskewing A & B
    preskewingA();
    preskewingB();

    // Perform the multiplication of the matrices
    for(i = 0; i < dims[0]; i++){
        sequentialMatrixMultiplication_REF(local_size,locA,locB,locC);
        shiftA();
        shiftB();
    }

    // Perform Postskewing A & B
    postskewingA();
    postskewingB();

    // Gather local matrices
    MPI_Gatherv(locC, mat_size * mat_size / w_size, MPI_DOUBLE, C, sendCounts, displacements, subarrtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(locA, mat_size * mat_size / w_size, MPI_DOUBLE, A, sendCounts, displacements, subarrtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(locB, mat_size * mat_size / w_size, MPI_DOUBLE, B, sendCounts, displacements, subarrtype, 0, MPI_COMM_WORLD);

    // Print time and the three final matrices 
    if(my_rank == 0){
    int time = MPI_Wtime() - start;
    printf("Time: %.4fs\n",time);

    printMatrix(mat_size, A);
    printMatrix(mat_size, B);
    printMatrix(mat_size, C);
    
    free(A);
    free(B);
    free(C);
    }


    MPI_Finalize();
    return 0;

}