#include <mpi.h>
#include <stdio.h>


double f(double x){

	return x*x;
}


int main(int argc, char** argv) {

	float a=0,b=1;
	float mesh = 40000;
	float h = (b-a)/mesh;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int i;
    float I=0;
    for(i=world_rank*10000;i<world_rank*10000+10000;i++){
    	I+=f(i*h);
    }
    if(world_rank!=0){
    	MPI_Send(&I,1,MPI_INT,0,1,MPI_COMM_WORLD);
    }

    if(world_rank==0){
    	float I1,I2,I3;
    	MPI_Recv(&I1,1,MPI_INT,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	MPI_Recv(&I2,1,MPI_INT,2,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	MPI_Recv(&I3,1,MPI_INT,3,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	I +=I1+I2+I3;
    	I = I*h;
    	printf("Risultato: %f\n",I);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
} 
