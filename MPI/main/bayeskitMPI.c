#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "gsl-sprng.h"
#include <mpi.h>

const int SAMPLENUM = 500;
int obsnum;
gsl_rng * r;
int procnum;

typedef struct {
	int prey;
	int predator;
} lvstate;

typedef struct {
	double time;
	double rawdata;
} obsdata;

void simPrior(lvstate *simData, MPI_Datatype datatype, MPI_Comm comm);
void rowsample(int *rows, double *w);
void stepLV(lvstate *state, double *t0p, double *dtp, double *lvParam, MPI_Datatype datatype, MPI_Comm comm);
void peturb(double *lvParam, MPI_Datatype datatype, MPI_Comm comm);

double obsLik(lvstate *mystate, double obs)
{

		const double SIGMA = 10.0;
		return log(gsl_ran_gaussian_pdf((obs - mystate->prey), SIGMA));

}

void mpiStepLV(lvstate *simData, double *t0p, double *dtp, double *lvParam, MPI_Datatype datatype, MPI_Comm comm)
{
	//printf("...callStepLV...");
	int i, j, offset, destNode, sourceNode, rank;
	int procnum, avgsize, samplemod;
	double t0, dt;
	MPI_Status status;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	//const unsigned int MAX_MPI_BUFFER = 524300000;
	//printf("rank in mpiStepLV %i \n", rank);
	//for (j = 0; j < SAMPLENUM; j++) {
		//MPI_Bcast(&simData[j].predator, 1, MPI_INT, 0, comm);
		//MPI_Bcast(&simData[j].prey, 1, MPI_INT, 0, comm);
	//}
	MPI_Bcast(lvParam, 3, MPI_DOUBLE, 0, comm);
	MPI_Bcast(t0p, 1, MPI_DOUBLE, 0, comm);
	MPI_Bcast(dtp, 1, MPI_DOUBLE, 0, comm);

	avgsize = SAMPLENUM/(procnum-1);
	samplemod = SAMPLENUM%(procnum-1);

	t0 = *t0p;
	dt = *dtp;
	if (rank==0)
	{
		int prey[SAMPLENUM];
		int pred[SAMPLENUM];
		for(i=0;i<SAMPLENUM;i++)
		{
			prey[i] = simData[i].prey;
			pred[i] = simData[i].predator;
		}
	}

	if (rank==0)
	{
		//printf("simdata 0 in call step lv: %i \n", simData[10].prey);
		//printf("---2---");

		//send tasks to other nodes

		offset = 0;

		for (destNode=1; destNode<procnum; destNode++){
			//MPI_Send(&lvParam, 3, MPI_DOUBLE, destNode, 0, MPI_COMM_WORLD);
			MPI_Send(&offset, 1, MPI_INT, destNode, 0, MPI_COMM_WORLD);
			for (i=0;i<avgsize;i++)
			{
				MPI_Send(&simData[offset+i].prey, 1, MPI_INT, destNode, 0, MPI_COMM_WORLD);
				MPI_Send(&simData[offset+j].predator, 1, MPI_INT, destNode, 0, MPI_COMM_WORLD);
			}
			offset = offset + avgsize;
		}


		for (j=0;j<samplemod;j++){
			stepLV(simData+offset+j, &t0, &dt, lvParam, MPI_FLOAT, MPI_COMM_WORLD);
		}


		// Wait to receive results from each node
		for (sourceNode=1; sourceNode<procnum; sourceNode++){

			MPI_Recv(&offset, 1, MPI_INT, sourceNode, 0, MPI_COMM_WORLD, &status);
			for (i=0;i<avgsize;i++)
			{
				MPI_Recv(&simData[offset+i].prey, 1, MPI_INT, sourceNode, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&simData[offset+i].predator, 1, MPI_INT, sourceNode, 0, MPI_COMM_WORLD, &status);
			}
		}
	}
	else if(rank>0 && rank<procnum)
	{
		//Receive data distributed from master node
		MPI_Recv(&offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		for (j=0;j<avgsize;j++)
		{
			MPI_Recv(&simData[offset+j].prey, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&simData[offset+j].predator, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}

		//doing stepLV for every slave processor/node
		for (j=0;j<avgsize;j++)
		{
			stepLV(simData+offset+j, &t0, &dt, lvParam, MPI_FLOAT, MPI_COMM_WORLD);
		}
		//printf("---5---");
		MPI_Send(&offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		for (j=0;j<avgsize;j++)
		{
			MPI_Send(&simData[offset+j].prey, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&simData[offset+j].predator, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}

		//printf("after On node rank %i offset %i simdata %i\n", rank, offset, simData[offset].prey);
		//printf("---6---");
	}


	//MPI_Buffer_detach( &buff, &b_size );
	return;
}

//particle filter: output new estimated lvstate and its likehood
void pfPropPar(lvstate ppstate[16], obsdata *myobsdata, double *lvParam, double *ll, MPI_Datatype datatype, MPI_Comm comm)
{
	//printf("...pfproppar...");
	int rank, procnum;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	//printf("rank in pfPropPar %i \n", rank);

	typedef struct {
		lvstate *vpptuple;
	} stateVec;

	int i, j, k, l;

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

	lvstate *simData = malloc(SAMPLENUM * sizeof(lvstate));
	if (simData == NULL) {
		printf("Fail to initiate simData...EXIT!");
		goto freeMem;
	}
	else {
		//printf("memset to simData to rank %i...", rank);
		memset(simData, 0, SAMPLENUM * sizeof(lvstate));
	}
	//initiate state: get simulated prey-predator tuple


		stateVec *myStateVec = malloc(obsnum * sizeof(stateVec));
		if (myStateVec == NULL) {
			printf("Fail to initiate myStateVec...Exit");
			goto freeMem;
		}
		memset(myStateVec, 0, obsnum * sizeof(stateVec));

		stateVec *tempStateVec = malloc(obsnum * sizeof(stateVec));
		if (tempStateVec == NULL) {
			printf("Fail to initiate tempStateVec...Exit");
			goto freeMem;
		}
		memset(tempStateVec, 0, obsnum * sizeof(stateVec));

	if (rank==0)
	{
		simPrior(simData, MPI_FLOAT, MPI_COMM_WORLD);
		for (j=0;j<SAMPLENUM;j++)
		{
			//printf("simdata in pfprop %i\n", (simData + j)->predator);
		}
		*ll = 0.0;

		if ((tempStateVec + 0)->vpptuple == NULL) {
			(tempStateVec + 0)->vpptuple = (lvstate*) malloc(
					SAMPLENUM * sizeof(lvstate));
			if ((tempStateVec + 0)->vpptuple == NULL) {
				printf("Fail to initiate myStateVec->vpptuple...Exit");
				goto freeMem;
			}
		}

		memset(tempStateVec[0].vpptuple, 0, SAMPLENUM * sizeof(lvstate));
		for (j = 0; j < SAMPLENUM; j++) {
			((tempStateVec + 0)->vpptuple + j)->predator = (simData + j)->predator;
			((tempStateVec + 0)->vpptuple + j)->prey = (simData + j)->prey;
		}
	}

	for (i = 1; i < obsnum; i++)
	{
		curTime = myobsdata[i-1].time;
		wSum = 0.0;

		//printf("-----b-----\n");
		if (rank==0)
		{
			if ((tempStateVec + i)->vpptuple == NULL) {
				//printf("init temp \n");
				(tempStateVec + i)->vpptuple = (lvstate*) malloc(
						SAMPLENUM * sizeof(lvstate));
				if ((tempStateVec + i)->vpptuple == NULL) {
					printf("Fail to initiate myStateVec->vpptuple...Exit");
					goto freeMem;
				}
				memset((tempStateVec + i)->vpptuple, 0,
						SAMPLENUM * sizeof(lvstate));
			}
		}


		//stepLV(simData + j, &curTime, &delTime, lvParam); put this into a function
		//printf("-----a-----\n");
		mpiStepLV(simData, &curTime, &delTime, lvParam, MPI_FLOAT, MPI_COMM_WORLD);

		if (rank==0)
		{
			for (j = 0; j < SAMPLENUM; j++)
			{
				((tempStateVec + i)->vpptuple + j)->predator =
						(simData + j)->predator;
				((tempStateVec + i)->vpptuple + j)->prey = (simData + j)->prey;
				//printf("tempState %i vpptuple %i , %p\n", i,j,(tempStateVec+i)->vpptuple+j);
				lw[j] = obsLik(((tempStateVec + i)->vpptuple + j), myobsdata[i].rawdata);
			}

			//get the max line weight
			for (j = 0; j < SAMPLENUM; j++) {
				if (j == 0) {
					maxLw = lw[0];
				} else {
					if (lw[j] > maxLw) {
						maxLw = lw[j];
					}
				}
			}
			//normalize the line wight and output list w[]
			for (l = 0; l < SAMPLENUM; l++) {
				w[l] = exp(lw[l] - maxLw);
				wSum = wSum + w[l];
			}
			rowsample(rows, w);
//		}

			//reorganize the vector with only keeping the rows in 'rows'
			for (k = 0; k <= i; k++)
			{
				if ((myStateVec + k)->vpptuple == NULL) {
					//printf("init mystatevec  \n");
					(myStateVec + k)->vpptuple = (lvstate*) malloc(
							SAMPLENUM * sizeof(lvstate));
					if ((myStateVec + k)->vpptuple == NULL) {
						printf("Fail to initiate myStateVec->vpptuple...Exit");
						goto freeMem;
					}
				}
	//			memset(myStateVec[k].vpptuple,0,SAMPLENUM * sizeof(lvstate));

				for (j = 0; j < SAMPLENUM; j++) {
					((myStateVec + k)->vpptuple + j)->predator =
							((tempStateVec + k)->vpptuple + rows[j])->predator;
					((myStateVec + k)->vpptuple + j)->prey =
							((tempStateVec + k)->vpptuple + rows[j])->prey;
					//printf("row %i %i,%i,%i,%i,%i\n", i,k,j,rows[j],((myStateVec+k)->vpptuple+j)->prey,((tempStateVec+k)->vpptuple+rows[j])->prey);
				}
				for (j = 0; j < SAMPLENUM; j++) {
					((tempStateVec + k)->vpptuple + j)->predator =
							((myStateVec + k)->vpptuple + j)->predator;
					((tempStateVec + k)->vpptuple + j)->prey =
							((myStateVec + k)->vpptuple + j)->prey;
				}
			}
			for (j = 0; j < SAMPLENUM; j++) {
				simData[j].predator = ((myStateVec + i)->vpptuple + j)->predator;
				simData[j].prey = ((myStateVec + i)->vpptuple + j)->prey;
				//printf("simdata in pfprop i%\n", simData[j].predator);
			}
			likehood = likehood + maxLw + log(wSum/SAMPLENUM);
			*ll = likehood;
			rows[SAMPLENUM] = NULL;
		}

	}
	if (rank==0)
	{
		for (i = 0; i < obsnum; i++) {
			ppstate[i].prey = ((myStateVec + i)->vpptuple + 0)->prey;
			ppstate[i].predator = ((myStateVec + i)->vpptuple + 0)->predator;
		}
	}

freeMem:
	//printf("free mem...");


	if (myStateVec != NULL) {
		//printf("free myStateVec.. rank %i\n.", rank);
		for (i = 0; ((myStateVec + i)->vpptuple != NULL) && i < obsnum; i++) {
			//printf("free myStateVec[vpptuple].. rank %i\n.", rank);
			free((myStateVec + i)->vpptuple);
			(myStateVec + i)->vpptuple = NULL;
		}
		free(myStateVec);
		myStateVec = NULL;
	}
	if (tempStateVec != NULL) {
		for (i = 0; NULL != (tempStateVec + i)->vpptuple && i < obsnum; i++) {
			//printf("free tempStateVec...i %i\n", i);
			free((tempStateVec + i)->vpptuple);
			(tempStateVec + i)->vpptuple = NULL;
		}
		free(tempStateVec);
		tempStateVec = NULL;
	}

	if (simData != NULL) {
		//printf("free simdata.. rank %i\n.", rank);
		free(simData);
		simData = NULL;
	}

	//printf("\n");
}

void arrStepLV(int *prey, int *pred, double *t0p, double *dtp, double *lvParam, MPI_Datatype datatype, MPI_Comm comm) {
	//printf("Starting stepLV...\n");
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	double t0 = *t0p, dt = *dtp, t;
	double h0, h1, h2, h3, u;

	while (dt > 0) {
		h1 = lvParam[0] * (*prey);
		h2 = lvParam[1] * (*prey) * (*pred);
		h3 = lvParam[2] * (*pred);
		h0 = h1 + h2 + h3;

		if ((h0 < (1e-10)) || ((*prey) >= 1000000))
			t = 1e99;
		else {
			t = gsl_ran_exponential(r, 1 / h0);
		}
		if (t > dt) {
			//PutRNGstate();
			//printf("out state %i, %i\n", state->prey, state->predator);
			return;
		} else {
			u = gsl_rng_uniform(r);
			if (u < (h1 / h0)) {
				(*prey) = (*prey) + 1;
				t0 = t0 + t;
				dt = dt - t;
			} else if (u < ((h1 + h2) / h0)) {
				(*prey) = (*prey) - 1;
				(*pred) = (*pred) + 1;
				t0 = t0 + t;
				dt = dt - t;
			} else {
				(*pred) = (*pred) - 1;
				t0 = t0 + t;
				dt = dt - t;
			}
		}
	}
	return;
}

void stepLV(lvstate *state, double *t0p, double *dtp, double *lvParam, MPI_Datatype datatype, MPI_Comm comm) {
	//printf("Starting stepLV...\n");
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	double t0 = *t0p, dt = *dtp, t;
	double h0, h1, h2, h3, u;

	while (dt > 0) {
		h1 = lvParam[0] * state->prey;
		h2 = lvParam[1] * state->prey * state->predator;
		h3 = lvParam[2] * state->predator;
		h0 = h1 + h2 + h3;

		if ((h0 < (1e-10)) || (state->prey >= 1000000))
			t = 1e99;
		else {
			t = gsl_ran_exponential(r, 1 / h0);
		}
		if (t > dt) {
			//PutRNGstate();
			//printf("out state %i, %i\n", state->prey, state->predator);
			return;
		} else {
			u = gsl_rng_uniform(r);
			if (u < (h1 / h0)) {
				state->prey = state->prey + 1;
				t0 = t0 + t;
				dt = dt - t;
			} else if (u < ((h1 + h2) / h0)) {
				state->prey = state->prey - 1;
				state->predator = state->predator + 1;
				t0 = t0 + t;
				dt = dt - t;
			} else {
				state->predator = state->predator - 1;
				t0 = t0 + t;
				dt = dt - t;
			}
		}
	}
	return;
}

//mcmc process
void runPmmhPath(int its, double *lvParam, double *obslik, lvstate ppstate[16],
		obsdata *myobsdata, MPI_Datatype datatype, MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	//printf("rank in runPmmhPath %i \n", rank);

	int i, j;
	double propParam[3];
	double curParam[3];
	double ll;
	double propMll, curMll = -1e99;
	lvstate curPath[obsnum];
	lvstate propPath[obsnum];
	FILE *fp = fopen("mcmc-out.txt", "wab+");
	if (fp == NULL) {
		printf("file is null.\n");
		return;
	}

	if (rank==0)
	{

		fprintf(fp, "th1,th2,th3,");
		for (i = 0; i <= 30; i++) {
			if ((i % 2) == 0) {
				if (i == 30) {
					fprintf(fp, "x%i,y%i", i, i);
				} else {
					fprintf(fp, "x%i,y%i,", i, i);
				}
			}
		}
		fprintf(fp, "\n");
		memcpy(curPath, ppstate, sizeof(lvstate) * obsnum);
		memcpy(curParam, lvParam, sizeof(double) * 3);
	}

	for (i = its; i > 0; i--) {
		//calculate measurements and decide moving to next state or not

		if (rank==0)
		{
			memcpy(propParam, curParam, sizeof(double) * 3);
			peturb(propParam, MPI_FLOAT, MPI_COMM_WORLD);
		}

		pfPropPar(propPath, myobsdata, propParam, &ll, MPI_FLOAT, MPI_COMM_WORLD);

		if (rank==0)
		{
			propMll = ll;
			if (log(gsl_ran_flat(r, 0.0, 1.0)) < (propMll - curMll)) {
				curMll = propMll;
				memcpy(curParam, propParam, sizeof(double) * 3);
				memcpy(curPath, propPath, sizeof(lvstate) * obsnum);
			}
		}

		if (rank==0)
		{
			printf("run for the %ith iteration.\n", i);
			//write current parameters and state sequence into output
			fprintf(fp, "%f,%f,%f,", curParam[0], curParam[1], curParam[2]);

			for (j = 0; j < obsnum; j++) {
				if (j == (obsnum - 1)) {
					fprintf(fp, "%i,%i", curPath[j].prey, curPath[j].predator);
				} else {
					fprintf(fp, "%i,%i,", curPath[j].prey, curPath[j].predator);
				}
			}
			fprintf(fp, "\n");
			printf("end for the %ith iteration.\n", i);
		}

	}
	fclose(fp);

}

void rowsample(int *rows, double *w)
{

	gsl_ran_discrete_t * grdp;
	int row;
	int i;

	for (i = 0; i < SAMPLENUM; i++) {
		grdp = gsl_ran_discrete_preproc(SAMPLENUM, w);
		row = (int) gsl_ran_discrete(r, grdp);
		rows[i] = row;
		gsl_ran_discrete_free(grdp);
	}

}

void peturb(double *lvParam, MPI_Datatype datatype, MPI_Comm comm) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//printf("Starting peturb...\n");
	if (rank==0)
	{
		const double SIGMA = 0.01;
		//printf("lvParam_before %f, %f, %f\n", lvParam[0], lvParam[1], lvParam[2]);
		lvParam[0] = lvParam[0] * exp(gsl_ran_gaussian(r, SIGMA));
		lvParam[1] = lvParam[1] * exp(gsl_ran_gaussian(r, SIGMA));
		lvParam[2] = lvParam[2] * exp(gsl_ran_gaussian(r, SIGMA));
		//printf("lvParam_after  %f, %f, %f\n", lvParam[0], lvParam[1], lvParam[2]);
	}
	return;
}

//to get the prior simulate value of prey and predator
void simPrior(lvstate *simData, MPI_Datatype datatype, MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//printf("Starting simPrior...\n");

	int i;

	const double PREY_MEAN = 50.0;
	const double PRED_MEAN = 100.0;

	for (i = 0; i < SAMPLENUM; i++) {
		simData[i].prey = gsl_ran_poisson(r, PREY_MEAN);
		simData[i].predator = gsl_ran_poisson(r, PRED_MEAN);
		//printf("simPrior %i\n", simData[i].predator);
	}
	//printf("End simPrior...\n");

	return;

}

void runModel(int its, MPI_Datatype datatype, MPI_Comm comm) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	obsnum = 16;
	obsdata myobsdata[obsnum];
	lvstate mylvstate[obsnum];
	double ll = 0.0;
	double lvParam[3] = { 1.0, 0.005, 0.6 };

	//produce time list(with even number in)
	if (rank==0) {
		int i;
		int j = 0;

		for (i = 0; i <= 30; i++) {
			if ((i % 2) == 0) {
				myobsdata[j].time = (double) i;
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
			while (fgets(line, 1024, file)) {
				myobsdata[i].rawdata = atof(line);
				i++;
			}
			fclose(file);
		}
		//printf("obsdata %f\n", myobsdata[5].rawdata);

	}

	runPmmhPath(its, lvParam, &ll, mylvstate, myobsdata, MPI_FLOAT, MPI_COMM_WORLD);
	return;
}

int main(int argc,char *argv[])
{
	int its, i, rank;
	double x, can, a, alpha, Nsum, mytime, maxtime, mintime;

	const gsl_rng_type * T;
	long seed;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	seed = time(NULL) * getpid();    // set seed
	gsl_rng_set(r, seed);

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);

	//printf("rank in main %i \n", rank);

	if (argc == 1) {
		its = 1000;
	}
	else {
		its = atoi(argv[1]);
	}

	runModel(its, MPI_FLOAT, MPI_COMM_WORLD);

	if (rank==0) {
		printf("Starting lv-pmcmc ...\n");
		printf("Running for %i", its);
		printf(" iterations, Done.\n");
	}

	gsl_rng_free (r);

	MPI_Finalize();
    return(EXIT_SUCCESS);
}
