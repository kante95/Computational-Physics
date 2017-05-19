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
    int npoints_core = mesh/world_size;
    int nstart = npoints_core*world_rank;
    int nend = npoints_core*(world_rank+1);
    if(world_rank == world_size-1){
    		nend += mesh-npoints_core*(world_size);
    }
    for(i=nstart;i<nend;i++){
    	I+=f(i*h);
    }
    if(world_rank!=0){
    	MPI_Send(&I,1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
    }

    if(world_rank==0){
    	
    	int i;
    	float I1;

    	for(i=1;i<world_size;i++){
    		MPI_Recv(&I1,1,MPI_FLOAT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    		I+=I1;

    	}

    	I = I*h;
    	printf("Risultato: %f\n",I);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
} 
