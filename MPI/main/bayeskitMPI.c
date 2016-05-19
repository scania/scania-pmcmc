/*
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "gsl-sprng.h"
#include <mpi.h>

const int SAMPLENUM = 1000;
int obsnum;
gsl_rng * r;

typedef struct {
   	int prey;
   	int predator;
}lvstate;

typedef struct {
	double time;
	double rawdata;
}obsdata;

void simPrior(lvstate *simData);
void rowsample(unsigned int *rows, double *w);
void stepLV(lvstate *state, double *t0p, double *dtp, double *lvParam);
void peturb(double *lvParam);
void ran_gen();

double obsLik(lvstate mystate, double obs)
{
	ran_gen();
	const double SIGMA=10.0;
	return log(gsl_ran_gaussian_pdf((obs-mystate.prey), SIGMA));
}

//particle filter: output new estimated lvstate and its likehood
void pfPropPar(lvstate ppstate[16], obsdata *myobsdata, double *lvParam, double *ll)
{
	typedef struct {
		lvstate vpptuple[SAMPLENUM];
	}stateVec;

	int i,j,k,l;
	stateVec myStateVec[obsnum];
	stateVec tempStateVec[obsnum];
	lvstate simData[SAMPLENUM];
	lvstate tmpSimData;
	double timeUnit = 2.0;
	int deltaTimeUnit = 1;
	double delTime = timeUnit * deltaTimeUnit;
	double curTime = 0.0;

	double lw[SAMPLENUM];
	double w[SAMPLENUM];
	double wSum;
	double maxLw;
	double likehood = 0.0;
	unsigned int rows[SAMPLENUM];

	//initiate state: get simulated prey-predator tuple
	simPrior(simData);

	*ll = 0.0;
	for (j=0;j<SAMPLENUM;j++) {
		tempStateVec[0].vpptuple[j] = simData[j];
	}

	for (i=1;i<obsnum;i++) {
		curTime = myobsdata[i-1].time;

		wSum = 0.0;

		for (j=0;j<SAMPLENUM;j++) {
			stepLV(&simData[j], &curTime, &delTime, lvParam);
			tempStateVec[i].vpptuple[j] = simData[j];
			lw[j] = obsLik(tempStateVec[i].vpptuple[j], myobsdata[i].rawdata);

			if(j==0) {
				maxLw = lw[0];
			}
			else{
				if(lw[j] > maxLw) {
					maxLw = lw[j];
				}
			}
		}

		for(l=0;l<SAMPLENUM;l++){
			w[l] = exp(lw[l] - maxLw);
			wSum = wSum + w[l];
		}
		rowsample(rows, w);

		//reorganize the vector with only keeping the rows in 'rows'
		for (k=0;k<obsnum;k++) {
			for (j=0;j<SAMPLENUM;j++) {
				myStateVec[k].vpptuple[j] = tempStateVec[k].vpptuple[rows[j]];
			}
		}
		likehood = likehood + maxLw + log(wSum/SAMPLENUM);
		*ll = likehood;
	}

	for (i=0;i<obsnum;i++) {
		ppstate[i].prey = myStateVec[i].vpptuple[0].prey;
		ppstate[i].predator = myStateVec[i].vpptuple[0].predator;
	}
}

void stepLV(lvstate *state, double *t0p, double *dtp, double *lvParam)
{
	ran_gen();
	//printf("Starting stepLV...\n");
	double t0=*t0p, dt=*dtp, t;
	double h0,h1,h2,h3,u;

	while (dt>0) {
		h1 = lvParam[0] * state->prey;
		h2 = lvParam[1] * state->prey * state->predator;
		h3 = lvParam[2] * state->predator;
		h0 = h1 + h2 + h3;

		if ((h0<(1e-10))||(state->prey>=1000000))
			t=1e99;
		else{
			t= gsl_ran_exponential(r, 1/h0);
		}
		if (t > dt) {
			//PutRNGstate();
			//printf("out state %i, %i\n", state->prey, state->predator);
			return;
		}
		else {
			u=gsl_rng_uniform(r);
			if (u<(h1/h0)){
				state->prey = state->prey + 1;
				t0 = t0 + t;
				dt = dt - t;
			}
			else if (u<((h1+h2)/h0)) {
				state->prey = state->prey - 1;
				state->predator = state->predator + 1;
				t0 = t0 + t;
				dt = dt - t;
			}
			else {
				state->predator = state->predator - 1;
				t0 = t0 + t;
				dt = dt - t;
			}
		}
	}

	return;
}

//mcmc process
void runPmmhPath(int its, double *lvParam, double *obslik, lvstate ppstate[16], obsdata *myobsdata)
{
	ran_gen();
	const double SIGMA = 0.01;
	int i, j;
	double propParam[3];
	double curParam[3];
	double ll;
	double propMll, curMll=-1e99;
	lvstate curPath[obsnum];
	lvstate propPath[obsnum];
	const char *text;

    FILE *fp = fopen("mcmc-out-mpi.txt", "wab+");
	if (fp == NULL) {
		printf("file is null.\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fp, "th1,th2,th3,");
    for(i=0;i<=30;i++)
    {
        if((i%2) == 0){
        	if (i==30) {
        		fprintf(fp, "x%i,y%i", i, i);
        	} else {
        		fprintf(fp, "x%i,y%i,", i, i);
        	}
        }
    }
    fprintf(fp, "\n");



    memcpy(curPath, ppstate, sizeof(lvstate) * obsnum);
	memcpy(curParam, lvParam, sizeof(double) * 3);
	for (i=its;i>0;i--)
	{
		//printf("---rng---%i\n", r->type->get);
		//printf("---rng---%i\n", r->state);
		//printf("---rng---%i\n", r->type->min);
		//calculate measurements and decide moving to next state or not
		memcpy(propParam, curParam, sizeof(double) * 3);
		peturb(propParam);

		pfPropPar(propPath, myobsdata, lvParam, &ll);
		propMll = ll;
		if (log(gsl_ran_flat(r, 0.0, 0.1)) < (propMll - curMll)){
			curMll = propMll;
			memcpy(curParam, propParam, sizeof(double) * 3);
			memcpy(curPath, propPath, sizeof(lvstate) * obsnum);
		}
		fprintf(fp, "%f,%f,%f,", curParam[0], curParam[1], curParam[2]);

		for (j=0;j<obsnum;j++) {
			if (j==(obsnum-1)) {
				fprintf(fp, "%i,%i", curPath[j].prey, curPath[j].predator);
			} else {
				fprintf(fp, "%i,%i,", curPath[j].prey, curPath[j].predator);
			}
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void rowsample(unsigned int *rows, double *w)
{
  	ran_gen();
  	gsl_ran_discrete_t * grdp;
  	double row;

  	int i;

  	for(i=0;i<SAMPLENUM;i++){
		grdp = gsl_ran_discrete_preproc(SAMPLENUM, w);
		row = gsl_ran_discrete(r, grdp);
		rows[i] = (int)row;
	}

  	return;
}

void peturb(double *lvParam)
{
	ran_gen();
	const double SIGMA = 0.01;
	//printf("lvParam_before %f, %f, %f\n", lvParam[0], lvParam[1], lvParam[2]);
  	lvParam[0] = lvParam[0] * exp(gsl_ran_gaussian(r, SIGMA));
  	lvParam[1] = lvParam[1] * exp(gsl_ran_gaussian(r, SIGMA));
  	lvParam[2] = lvParam[2] * exp(gsl_ran_gaussian(r, SIGMA));
  	//printf("lvParam_after  %f, %f, %f\n", lvParam[0], lvParam[1], lvParam[2]);
}

//to get the prior simulate value of prey and predator
void simPrior(lvstate *simData)
{
	ran_gen();
	lvstate sim[SAMPLENUM];
	int i; 

	const double PREY_MEAN = 50.0;
	const double PRED_MEAN = 100.0;
  	
	for(i=0;i<SAMPLENUM;i++){
		simData[i].prey = gsl_ran_poisson(r, PREY_MEAN);
		simData[i].predator = gsl_ran_poisson(r, PRED_MEAN);
	}

}

void runModel(int its)
{
	obsnum = 16;
	obsdata myobsdata[obsnum];
	lvstate simData[SAMPLENUM], mylvstate[obsnum];
	double ll = 0.0;

	double lvParam[3] = {1.0, 0.005, 0.6};

	//produce time list(with even number in)
    int i;
    int j = 0;

    for(i=0;i<=30;i++)
    {
        if((i%2) == 0){
        	myobsdata[j].time = (double)i;
            j++;
        }
    }
    
    //read external txt file
    FILE *file = fopen("LVpreyNoise10.txt", "r");
    char line[1024];
	if (file == NULL) {
		printf("file is null.\n");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
    if (file != NULL) {
    	i = 0;
		while(fgets(line,1024,file)) {
			myobsdata[i].rawdata = atof(line);
        	i++;
    	}
    	fclose(file);
    }

    runPmmhPath(its, lvParam, &ll, mylvstate, myobsdata);
    return;
}

void ran_gen()
{
	const gsl_rng_type * T;
	long seed;
	gsl_rng_env_setup();
  	T = gsl_rng_default;
  	r = gsl_rng_alloc (T);
    seed = time (NULL) * getpid();    // set seed
  	gsl_rng_set (r, seed);

}

void tets1()
{
	printf("poisson test1 %i\n", gsl_ran_poisson(r, 50));
	printf("poisson test1 %i\n", gsl_ran_poisson(r, 50));
	printf("poisson test1 %i\n", gsl_ran_poisson(r, 50));
	printf("poisson test1 %i\n", gsl_ran_poisson(r, 50));
	printf("poisson test1 %i\n", gsl_ran_poisson(r, 50));
	printf("poisson test1 %i\n", gsl_ran_poisson(r, 50));
}

int main(int argc,char *argv[])
{
	ran_gen();
	int its, i, k;
	double x, can, a, alpha, Nsum, mytime, maxtime, mintime;

	printf("Starting lv-pmcmc ...\n");
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&k);
	mytime = MPI_Wtime();

	if (argc == 1) {
		its = 40;
	}
	else {
		its = atoi(argv[1]);
	}
	runModel(its);
	
	gsl_rng_free (r);
    printf("Running for %i", its);
    printf(" iterations, Done.\n");

    mytime = MPI_Wtime() - mytime;
    MPI_Reduce(&mytime,&Nsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&mytime, &maxtime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    if (k==0){
      printf("pmcmc time in total is %f\n", Nsum);
      printf("pmcmc maxtime is %f\n", maxtime);
      printf("pmcmc mintime is %f\n", mintime);
    }
    MPI_Finalize();
    return(EXIT_SUCCESS);
    //return 0;
}
*/
