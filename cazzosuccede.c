#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

#define min(a,b) ((a > b)?b:a)

#define N 125
#define first_temp 1.5
#define last_temp 3


void printmatrix(double matrix[N][3]){
	int i;
	for(i=0;i<N;i++){
		printf("Particella %d: [%lf] [%lf] [%lf]\n",i+1,matrix[i][0],matrix[i][1],matrix[i][2]);

	}

}

int find_maximum(double a[], int n) {
  int c, max, index;
 
  max = a[0];
  index = 0;
 
  for (c = 1; c < n; c++) {
    if (a[c] > max) {
       index = c;
       max = a[c];
    }
  }
 
  return index;
}


void copy_matrix(double source[N][3],double destination[N][3]){
	int i,j;
	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			destination[i][j] = source[i][j];
}


void initialize(double L,double positions[N][3]){
    float i = 5;
    int index = 0;
    float j,k,m;
    for(j =0;j<i;j++){
    	for(k =0;k<i;k++){
    		for(m =0;m<i;m++){
    			positions[index][0]= j*(L/i) - L/2.0;
                positions[index][1] = k*(L/i) - L/2.0;
                positions[index][2] = m*(L/i) - L/2.0;
                //printf("x: %lf y: %lf z: %lf",positions[index][0],positions[index][1],positions[index][2]);
                index++;
			}
    	}
    }
}

double myrandom(double from,double to){
	return ((double)rand()/(double)(RAND_MAX))*(to-from) + from;
}

double lennard_jones(double r){
	double r2 = r*r;
	double r6 = r2*r2*r2;
	double r12 = r6*r6;
	return 4*(1.0/r12-1.0/r6);
}

double calculate_potential(double positions[N][3],double L){

    double potential = 0;
    double r;
    int i,j,x;
    double d;
    for(i=0;i<N;i++){
    	for(j=i+1;j<N;j++){
    		r=0;
    		for(x = 0; x<3;x++){
    			d = positions[i][x] - positions[j][x]-L*rint((positions[i][x] - positions[j][x])/L);
    			r+=d*d;
    		}
    		r = sqrt(r);
    		if(r<L/2.0 && r!=0.0){
    			potential+=2*lennard_jones(r);
    		}
    	}
    }

    return potential;
}

double find_delta(double temperature,double L, double delta){
    int steps = 0;
    int acceptance = 0;
    float acceptance_ratio=0;
    double beta = 1.0/temperature;
    int i,j;
    double positions[N][3];
    double potential;
    double new_potential;
    double p;
    double xi;
    double new_positions[N][3];

    initialize(L,positions);
    potential = calculate_potential(positions,L);
    while(steps<1000){
    	//printmatrix(positions);
        for(j=0; j<N;j++){
        	for(i=0; i<3;i++){
        		new_positions[j][i] = positions[j][i]+myrandom(-delta/2,delta/2);
                new_positions[j][i] -= L*rint(new_positions[j][i]/L);
        	}
        }  
        //printmatrix(new_positions);
        new_potential = calculate_potential(new_positions,L);
        p = min(1,exp(-beta*(new_potential - potential)));
        xi = ((double)rand()/(double)(RAND_MAX));
        if (xi < p){
        	copy_matrix(new_positions,positions);
            potential = new_potential;
            acceptance+=1;
        }
        steps+=1;

        
    }
    acceptance_ratio = acceptance*100/(steps-1);
    //printf("Delta: %lf, Acceptance ratio: %lf\n",delta,acceptance_ratio);
    if(acceptance_ratio<60 && acceptance_ratio>40){
        return delta;
    }
    else if(acceptance_ratio == 0.0){
        return find_delta(temperature,L, delta/1.5);
    }
    else{
        return find_delta(temperature,L, delta*(acceptance_ratio/50)); //???? più graduale delta + (acceptance_ratio-50)%(L/2)
        //delta*(acceptance_ratio/50)
    }
}


double simulation(float density,int M,int rank){

    float temp_simul = 2;
 
 	int num_temp = 50;
 	double temp[num_temp];
	int i,t,k,j;
	double step = (last_temp-first_temp)/num_temp;
	for (i = 0; i < num_temp; i++) {
		temp[i] = step*i + first_temp;
	}

    double Vt[num_temp];
    double Vt2[num_temp];
    double betas[num_temp];
    for (i = 0; i < num_temp; i++) {
		betas[i] = 1.0/temp[i];
	}
    double normalization[num_temp];
    
    double V = N/density;
    double L =pow(V,1.0/3.0);
    //printf("lato: %lf",L);

    //printf("Attendi, trovo il migliore delta.... Processore %d\n",rank);
   
    double delta = find_delta(temp_simul,L,0.01/(pow(density,1.0/3.0)));//find_delta(temp_simul,L,0.02/(pow(density,1.0/3.0)));//0.05/(pow(density,1.0/3.0));//find_delta(temp_simul,L,0.01);
  
    //printf("Processore %d Delta migliore trovato: %lf, adesso inizio la simulazione\n",rank,delta);

    int acceptance = 0;

    double beta = 1.0/temp_simul;
    //inizializzazione
    double positions[N][3];
    initialize(L,positions); 
    //printmatrix(positions);
    double potential = calculate_potential(positions,L);
    //printf("Potenziale %lf\n",potential);
    double Vt_simul = potential;
    double Vt2_simul =  potential*potential;

    double new_positions[N][3];
    double new_potential;
    double p;
    double xi;
    double Cv[num_temp];
    int cv_max;
    double temp_max;

    double prova;
    char filename[80];
    FILE *f;

    for(t=0;t<num_temp;t++)
    {
    	Vt[t] =  potential*exp((beta - betas[t])*Vt_simul);
        Vt2[t] =  potential*potential*exp((beta - betas[t])*Vt_simul);
        normalization[t] = exp((beta - betas[t])*Vt_simul);

    }
    //loop principale della catena
    for(k=0;k<M;k++)
    {
        for(j=0 ; j<N;j++){
        	for(i=0 ; i<3;i++){
                prova = myrandom(-delta/2,delta/2);
                //printf("Random: %lf\n",prova);
        		new_positions[j][i] = positions[j][i]+prova;
                new_positions[j][i] -= L*rint(new_positions[j][i]/L);
        	}
        }
        //printmatrix(new_positions);
        new_potential = calculate_potential(new_positions,L);
        //printf("Potenziale %lf\n",new_potential);
        //sleep(10);
        p = min(1,exp(-beta*(new_potential - potential)));
        xi = ((double)rand()/(double)(RAND_MAX));
        if (xi < p){
        	copy_matrix(new_positions,positions);
            potential = new_potential;
            acceptance+=1;
            //printf("accettata!!\n");
        }
        Vt_simul += potential;
        Vt2_simul+= potential*potential;

        for(t=0;t<num_temp;t++)
    	{
    		Vt[t] +=  potential*exp((beta - betas[t])*potential);
        	Vt2[t] +=  potential*potential*exp((beta - betas[t])*potential);
        	normalization[t] += exp((beta - betas[t])*potential);
    	}
    	//printf("temp 1.5 V: %lf V2: %lf \n",Vt[0],Vt2[0]);
    }

    
    for(t=0;t<num_temp;t++)
    {
        Vt[t] /=  normalization[t];
        Vt2[t] /=  normalization[t];       
    }
    Vt_simul/=M;
    Vt2_simul/=M;

    //printf("Processore: %d, Acceptance ratio: %lf\n",rank,acceptance*100/(double)M);
    
    for(i=0;i<num_temp;i++)
    {
        Cv[i] = (Vt2[i]-Vt[i]*Vt[i])/(temp[i]*temp[i]);
    }
    
    cv_max = find_maximum(Cv,num_temp);
    temp_max = temp[cv_max];

    sprintf(filename,"density%f.csv",density);
    f = fopen(filename, "w+");

    for(i=0;i<num_temp;i++){
        //printf("Cv: %lf temp: %lf\n",Cv[i],temp[i]);
        //printf("Processore %d: %lf,%lf\n",rank,Cv[i],temp[i]);
        fprintf(f,"%lf,%lf\n",Cv[i],temp[i]);
    }
    //printf("Temperatura massima: %lf\n",temp_max);

    fclose(f);
    printf("Processore %d Densità: %f Temperatura massima: %lf acceptance ratio: %lf\n",rank,density,temp_max,acceptance*100/(double)M);

    return temp_max;
}

int main(){
    time_t t;
    int num_density=28;
    srand((unsigned) time(&t));

    int steps = 50000;
    float density[num_density];
    float temp[num_density];

    int i;
    double step = (1.0-0.01)/num_density;
    for (i = 0; i < num_density; i++) {
        density[i] = step*i + 0.01;
    }

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    int nstart = world_rank*num_density/world_size;
    int nend = nstart + num_density/world_size;

    for(i=nstart;i<nend;i++){
        temp[i] = simulation(density[i],steps,world_rank);
        //printf("Processore: %d, densità: %lf, temperatura: %lf\n",world_rank,density[i],temp[i]);
    }




    MPI_Finalize();


    //printf("Densità |Temperatura\n");
    //for(i=0;i<12;i++){
     //   printf("%lf, %lf\n",density[i],temp[i]);
    //}
    return 0;
}

/*int main(){
	time_t t;
	float density = 0.3;
    double temp;
	int steps = 50000;
	srand((unsigned) time(&t));
	temp = simulation(density,steps,1);
	printf("temp: %lf densità: %f\n",temp,density);
	return 0;
} */
