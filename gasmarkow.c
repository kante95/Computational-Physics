#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

#define min(a,b) ((a > b)?b:a)

#define N 125


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


void initialize(float L,double positions[N][3]){
    float i = pow(N,1.0/3.0);
    int index = 0;
    for(float j =0;j<L-L/i;j+=L/i){
    	for(float k =0;k<L-L/i;k+=L/i){
    		for(float m =0;m<L-L/i;m+=L/i){
    			positions[index][0]= j - L/2.0;
                positions[index][1] = k - L/2.0;
                positions[index][2] = m - L/2.0;
                //printf("x: %lf y: %lf z: %lf",positions[index][0],positions[index][1],positions[index][2]);
                index++;
			}
    	}
    }
}

double myrandom(double from,double to){
	return ((float)rand()/(float)(RAND_MAX))*(to-from) + from;
}

double lennard_jones(double r){
	return 4*(pow(r,-12.0)-pow(r,-6.0));
}

double calculate_potential(double positions[N][3],float L){

    double potential = 0;
    double r=0;
    int i,j,x;
    double d;
    for(i=0;i<N;i++){
    	for(j=i+1;j<N;j++){
    		for(x = 0; x<3;x++){
    			d = positions[i][x] - positions[j][x]-rint((positions[i][x] - positions[j][x])/(double)L);
    			//printf("i:%d j:%d x:%d positions[i,x]:%lf positions[j,x]%lf r: %lf\n",i,j,x,positions[i][x],positions[j][x],d);
    			//sleep(1);
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

double find_delta(double temperature,float L){
    double delta = 0.1;
    int steps = 0;
    int acceptance = 0;
    float acceptance_ratio;
    double beta = 1.0/temperature;
    int i,j;
    double positions[N][3];
    initialize(L,positions);
    double potential = calculate_potential(positions,L);
    double new_potential;
    double p;
    double xi;
    double new_positions[N][3];

    while(acceptance_ratio<45 || acceptance_ratio > 55){
        for(j=0; j<N;j++){
        	for(i=0; i<3;i++){
        		new_positions[j][i] = positions[j][i]+myrandom(-delta/2,delta/2);
                new_positions[j][i] = new_positions[j][i]-rint(new_positions[j][i]/L);
        	}
        }  
        new_potential = calculate_potential(new_positions,L);
        //printf("potential: %lf\n",potential);
        p = min(1,pow(2.718,-beta*(new_potential - potential)));
        xi = ((double)rand()/(double)(RAND_MAX));
        if (xi < p){
        	copy_matrix(new_positions,positions);
            potential = new_potential;
            acceptance+=1;
        }
        steps+=1;
        if (steps%1000==0){
            acceptance_ratio = acceptance/10;
            //printf("Delta: %lf, Acceptance ratio: %f\n",delta,acceptance_ratio);
            if(acceptance==0.0){
                delta = delta/10;
                acceptance = 0;
            }
            else{
                acceptance = 0;
                delta = delta*(acceptance_ratio/50);
            }
        }
    }
    return delta;
}


double simulation(float density,int M){

    float temp_simul = 2;
 
 	int num_temp = 50;
 	double temp[num_temp];
	int i,t,k,j;
	double step = 3.0/num_temp;
	for (i = 0; i < num_temp; i++) {
		temp[i] = step*i + 1.5;
	}

    double Vt[num_temp];
    double Vt2[num_temp];
    double betas[num_temp];
    for (i = 0; i < num_temp; i++) {
		betas[i] = 1.0/temp[i];
	}
    double normalization[num_temp];
    
    float V = N/density;
    float L =pow(V,1.0/3.0);

    printf("Attendi, trovo il migliore delta....\n");
   
    double delta = find_delta(temp_simul,L);
  
    printf("Delta migliore trovato: %lf, adesso inizio la simulazione\n",delta);

    int acceptance = 0;

    double beta = 1/temp_simul;
    //inizializzazione
    double positions[N][3];
    initialize(L,positions); 
    double potential = calculate_potential(positions,L);
    
    double Vt_simul = potential;
    double Vt2_simul =  potential*potential;

    double new_positions[N][3];
    double new_potential;
    double p;
    double xi;
    double Cv[num_temp];
    int cv_max;
    double temp_max;

    for(t=0;t<num_temp;t++)
    {
    	Vt[t] =  potential*pow(2.718,(beta - betas[t])*Vt_simul);
        Vt2[t] =  potential*potential*pow(2.718,(beta - betas[t])*Vt_simul);
        normalization[t] = pow(2.7281,(beta - betas[t])*Vt_simul);

    }
    //loop principale della catena
    for(k=0;k<M;k++)
    {
        for(j=0 ; j<N;j++){
        	for(i=0 ; i<3;i++){
        		new_positions[j][i] = positions[j][i]+myrandom(-delta/2,delta/2);
                new_positions[j][i] -= rint(new_positions[j][i]/L);
        	}
        }

        new_potential = calculate_potential(new_positions,L);
        p = min(1,pow(2.718,-beta*(new_potential - potential)));
        xi = ((double)rand()/(double)(RAND_MAX));
        if (xi < p){
        	copy_matrix(new_positions,positions);
            potential = new_potential;
            acceptance+=1;
        }
        Vt_simul += potential;
        Vt2_simul+= potential*potential;

        for(t=0;t<num_temp;t++)
    	{
    		Vt[t] +=  potential*pow(2.718,(beta - betas[t])*Vt_simul);
        	Vt2[t] +=  potential*potential*pow(2.718,(beta - betas[t])*Vt_simul);
        	normalization[t] += pow(2.7281,(beta - betas[t])*Vt_simul);
    	}
    }

    
    for(t=0;t<num_temp;t++)
    {
        Vt[t] /=  normalization[t];
        Vt2[t] /=  normalization[t];       
    }
    Vt_simul/=M;
    Vt2_simul/=M;

    printf("temp 1.5 V: %lf V2: %lf acceptance ratio: %lf\n",Vt[0],Vt2[0],acceptance*100/(double)M);
    
    for(i=0;i<num_temp;i++)
    {
        Cv[i] = (Vt2[i]-Vt[i]*Vt[i])/(temp[i]*temp[i]);
    }
    
    cv_max = find_maximum(Cv,num_temp);
    temp_max = temp[cv_max];

    for(i=0;i<num_temp;i++){
    	printf("Cv: %lf temp: %lf\n",Cv[i],temp[i]);
    }
    printf("Temperatura massima: %lf\n",temp_max);

    return temp_max;
}


int main(){
	time_t t;
	srand((unsigned) time(&t));
	float density = 0.01;
	int step = 50000;

	simulation(density,step);
	return 0;
}