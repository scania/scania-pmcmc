#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <errno.h>
#include <string.h>

const int SAMPLENUM = 1000;
int obsnum;
gsl_rng * r;
long seed;

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

double obsLik(lvstate *mystate, double obs)
{
	//ran_gen();
	const double SIGMA=10.0;
	return log(gsl_ran_gaussian_pdf((obs- mystate->prey), SIGMA));
}

//particle filter: output new estimated lvstate and its likehood
void pfPropPar(lvstate ppstate[16], obsdata *myobsdata, double *lvParam, double *ll)
{
	printf("Starting pfPropPar...\n");
	typedef struct {
		lvstate *vpptuple;
	}stateVec;

	int i,j,k,l;

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

	printf("Starting pfPropPar 1 ...\n");
	stateVec *myStateVec = malloc(obsnum * sizeof(stateVec));
	if (myStateVec == NULL)
	{
		printf("Fail to initiate myStateVec...Exit");
		goto freeMem;
	}
	memset(myStateVec,0,obsnum * sizeof(stateVec));

	stateVec *tempStateVec = malloc(obsnum * sizeof(stateVec));
	if (tempStateVec == NULL)
	{
		printf("Fail to initiate tempStateVec...Exit");
		goto freeMem;
	}
	memset(tempStateVec,0,obsnum * sizeof(stateVec));

	lvstate *simData = malloc(SAMPLENUM * sizeof(lvstate));
	if (simData == NULL)
	{
		printf("Fail to initiate simData...EXIT!");
		goto freeMem;
	}
	memset(simData,0,SAMPLENUM * sizeof(lvstate));
	//initiate state: get simulated prey-predator tuple
	simPrior(simData);

	*ll = 0.0;
	printf("Starting pfPropPar 2..\n");
	for (j=0;j<SAMPLENUM;j++) {
		tempStateVec->vpptuple = (lvstate*)  malloc(SAMPLENUM * sizeof(lvstate));
		if (tempStateVec->vpptuple == NULL)
		{
			printf("Fail to initiate myStateVec->vpptuple...Exit");
			goto freeMem;
			//free(tempStateVec);
			//goto error;
		}
		memset(tempStateVec->vpptuple,0,SAMPLENUM * sizeof(lvstate));

		//printf("simdata j= %i\n", j);
		tempStateVec->vpptuple[j].predator = simData[j].predator;
		tempStateVec->vpptuple[j].prey = simData[j].prey;
		//printf("simData %i\n", tempStateVec[0].vpptuple[j].prey);
		//simData[j].prey;
	}

	for (i=1;i<obsnum;i++) {
		curTime = myobsdata[i-1].time;

		wSum = 0.0;

		(tempStateVec+i)->vpptuple = (lvstate*)  malloc(SAMPLENUM * sizeof(lvstate));
		if ((tempStateVec+i)->vpptuple == NULL)
		{
			printf("Fail to initiate myStateVec->vpptuple...Exit");
			goto freeMem;
		}
		memset((tempStateVec+i)->vpptuple,0,SAMPLENUM * sizeof(lvstate));

		for (j=0;j<SAMPLENUM;j++) {
			stepLV(simData+j, &curTime, &delTime, lvParam);
			((tempStateVec+i)->vpptuple+j)->predator = (simData+j)->predator;
			((tempStateVec+i)->vpptuple+j)->prey = (simData+j)->prey;
			//printf("tempState %i vpptuple %i , %p\n", i,j,(tempStateVec+i)->vpptuple+j);
			lw[j] = obsLik(((tempStateVec+i)->vpptuple+j), myobsdata[i].rawdata);

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
//		rows[SAMPLENUM] = {2,3,4,1,2,6,7,4,1,8};
		/*
		for (i=0;i<3;i++)
		{
			printf("row %d,%i\n", i,rows[i]);
		}

		for (i=3;i<6;i++)
		{
			printf("row %d,%i\n", i,rows[i]);
		}
		*/
		//reorganize the vector with only keeping the rows in 'rows'
		for (k=0;k<=i;k++) {
			(myStateVec+k)->vpptuple = (lvstate*) malloc(SAMPLENUM * sizeof(lvstate));
//			memset(myStateVec[k].vpptuple,0,SAMPLENUM * sizeof(lvstate));
			if ((myStateVec+k)->vpptuple == NULL)
			{
				printf("Fail to initiate myStateVec->vpptuple...Exit");
				goto freeMem;
			}
			//printf("reorganize, k %i\n",k);
			for (j=0;j<SAMPLENUM;j++) {
				//printf("row %i %i,%i,%i,%p,%p\n", i,k,j,rows[j],(myStateVec+k)->vpptuple+j,(tempStateVec+k)->vpptuple+rows[j]);
				//printf("prey  %i\n", (tempStateVec[k].vpptuple+rows[j])->prey);
				((myStateVec+k)->vpptuple+j)->predator = ((tempStateVec+k)->vpptuple+rows[j])->predator;
				((myStateVec+k)->vpptuple+j)->prey = ((tempStateVec+k)->vpptuple+rows[j])->prey;
				//(myStateVec+k)->(vpptuple+j)->prey = (tempStateVec+k)->(vpptuple+rows[j])->prey;
			}
			for (j=0;j<SAMPLENUM;j++) {
				((tempStateVec+k)->vpptuple+j)->predator = ((myStateVec+k)->vpptuple+j)->predator;
				((tempStateVec+k)->vpptuple+j)->prey = ((myStateVec+k)->vpptuple+j)->prey;
			}
			for (j=0;j<SAMPLENUM;j++) {
				simData[j].predator = ((myStateVec+i)->vpptuple+j)->predator;
				simData[j].prey = ((myStateVec+i)->vpptuple+j)->prey;
			}
		}
		likehood = likehood + maxLw + log(wSum/SAMPLENUM);
		*ll = likehood;
	}

	for (i=0;i<obsnum;i++) {
		ppstate[i].prey = ((myStateVec+i)->vpptuple)->prey;
		ppstate[i].predator = ((myStateVec+i)->vpptuple)->predator;
	}

	printf("Starting pfPropPar 3...\n");

freeMem:
	//printf("free mem...");
	if (myStateVec != NULL) {
		//printf("free mem 2...");
		for (i=0;NULL != (myStateVec+i)->vpptuple && i<obsnum;i++) {
			//printf("free mem...");
			free((myStateVec+i)->vpptuple);
			(myStateVec+i)->vpptuple = NULL;
		}
		free(myStateVec);
		myStateVec = NULL;
	}
	if (tempStateVec != NULL) {
		for (i=0;NULL != (tempStateVec+i)->vpptuple && i<obsnum;i++) {
			free((tempStateVec+i)->vpptuple);
			(tempStateVec+i)->vpptuple = NULL;
		}
		free(tempStateVec);
		tempStateVec = NULL;
	}
	if (simData != NULL) {
		free(simData);
		simData = NULL;
	}

}

void stepLV(lvstate *state, double *t0p, double *dtp, double *lvParam)
{
	//printf("Starting stepLV...\n");
	//ran_gen();
	double t0=*t0p, dt=*dtp;
	double t;
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
	printf("Starting runPmmhPath...\n");
	//ran_gen();
	int i, j;
	double propParam[3];
	double curParam[3];
	double ll;
	double propMll, curMll=-1e99;
	lvstate curPath[obsnum];
	lvstate propPath[obsnum];

    FILE *fp = fopen("mcmc-out.txt", "wab+");
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

		pfPropPar(propPath, myobsdata, propParam, &ll);
		propMll = ll;
		if (log(gsl_ran_flat(r, 0.0, 1.0)) < (propMll - curMll)){
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
	//printf("Starting rowsample...\n");
  	gsl_ran_discrete_t * grdp;
  	double row;

  	int i;

  	for(i=0;i<SAMPLENUM;i++){
		grdp = gsl_ran_discrete_preproc(SAMPLENUM, w);
		row = gsl_ran_discrete(r, grdp);
		rows[i] = (int)row;
	}

  	free(grdp);
  	return;
}

void peturb(double *lvParam)
{
	//printf("Starting peturb...\n");
	const double SIGMA = 0.01;
  	lvParam[0] = lvParam[0] * exp(gsl_ran_gaussian(r, SIGMA));
  	lvParam[1] = lvParam[1] * exp(gsl_ran_gaussian(r, SIGMA));
  	lvParam[2] = lvParam[2] * exp(gsl_ran_gaussian(r, SIGMA));
}

//to get the prior simulate value of prey and predator
void simPrior(lvstate *simData)
{
	//printf("Starting simPrior...\n");
	int i; 

	const double PREY_MEAN = 50.0;
	const double PRED_MEAN = 100.0;
  	
	for(i=0;i<SAMPLENUM;i++){
		simData[i].prey = gsl_ran_poisson(r, PREY_MEAN);
		simData[i].predator = gsl_ran_poisson(r, PRED_MEAN);
	}
	//printf("End simPrior...\n");

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
	const gsl_rng_type * T;
	gsl_rng_env_setup();
  	T = gsl_rng_default;
  	r = gsl_rng_alloc (T);
  	ran_gen();

	int its;
	printf("Starting MPI main...\n");

	if (argc == 1) {
		its = 20;
	}
	else {
		its = atoi(argv[1]);
	}
	runModel(its);
	
	/*
	printf("poisson %i\n", gsl_ran_poisson(r, 50));
	printf("poisson %i\n", gsl_ran_poisson(r, 50));
	printf("poisson %i\n", gsl_ran_poisson(r, 50));
	printf("poisson %i\n", gsl_ran_poisson(r, 50));
	printf("poisson %i\n", gsl_ran_poisson(r, 50));
	printf("poisson %i\n", gsl_ran_poisson(r, 50));

	tets1();
	*/

	gsl_rng_free (r);
    printf("Running for %i", its);
    printf(" iterations, Done.");
    return 0;
}
