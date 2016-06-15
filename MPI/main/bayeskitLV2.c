/*
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>

const int SAMPLENUM = 200;
int obsnum;
gsl_rng * r;

typedef struct {
	int prey;
	int predator0;
	int predator1;
} lvstate;

typedef struct {
	double time;
	double rawdata;
} obsdata;

void simPrior(lvstate *simData);
void rowsample(int *rows, double *w);
void stepLV(lvstate *state, double *t0p, double *dtp, double *lvParam);
void peturb(double *lvParam);

double obsLik(lvstate *mystate, double obs) {
	const double SIGMA = 20.0;
	double ll = 0.0;
	ll = log(gsl_ran_gaussian_pdf((obs - mystate->prey), SIGMA));
	return ll;
}

void callStepLV(lvstate *simData, double *t0p, double *dtp, double *lvParam) {
	int j;
	for (j=0;j<SAMPLENUM;j++) {
		stepLV(simData + j, t0p, dtp, lvParam);
	}

}

//particle filter: output new estimated lvstate and its likehood
void pfPropPar(lvstate ppstate[16], obsdata *myobsdata, double *lvParam,
		double *ll) {
	//printf("Starting pfPropPar...\n");
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

	//printf("Starting pfPropPar 1 ...\n");
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

	lvstate *simData = malloc(SAMPLENUM * sizeof(lvstate));
	if (simData == NULL) {
		printf("Fail to initiate simData...EXIT!");
		goto freeMem;
	}
	memset(simData, 0, SAMPLENUM * sizeof(lvstate));
	//initiate state: get simulated prey-predator tuple
	simPrior(simData);

	*ll = 0.0;

	if ((tempStateVec + 0)->vpptuple == NULL) {
		(tempStateVec + 0)->vpptuple = (lvstate*) malloc(
				SAMPLENUM * sizeof(lvstate));
		if ((tempStateVec + 0)->vpptuple == NULL) {
			printf("Fail to initiate myStateVec->vpptuple...Exit");
			goto freeMem;
			//free(tempStateVec);
			//goto error;
		}
	}

	memset(tempStateVec[0].vpptuple, 0, SAMPLENUM * sizeof(lvstate));
	for (j = 0; j < SAMPLENUM; j++) {
		((tempStateVec + 0)->vpptuple + j)->predator0 = (simData + j)->predator0;
		((tempStateVec + 0)->vpptuple + j)->predator1 = (simData + j)->predator1;
		((tempStateVec + 0)->vpptuple + j)->prey = (simData + j)->prey;
	}

	for (i = 1; i < obsnum; i++) {
		curTime = myobsdata[i - 1].time;

		wSum = 0.0;

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

		//printf("tempvec 1 %i\n", ((tempStateVec+0)->vpptuple+1)->prey);
		callStepLV(simData, &curTime, &delTime, lvParam);

		for (j = 0; j < SAMPLENUM; j++) {
			//stepLV(simData + j, &curTime, &delTime, lvParam);
			((tempStateVec + i)->vpptuple + j)->predator0 =
					(simData + j)->predator0;
			((tempStateVec + i)->vpptuple + j)->predator1 =
					(simData + j)->predator1;
			((tempStateVec + i)->vpptuple + j)->prey = (simData + j)->prey;
			//printf("tempState %i vpptuple %i , %p\n", i,j,(tempStateVec+i)->vpptuple+j);
			lw[j] = obsLik(((tempStateVec + i)->vpptuple + j),
					myobsdata[i].rawdata);

			if (j == 0) {
				maxLw = lw[0];
			} else {
				if (lw[j] > maxLw) {
					maxLw = lw[j];
				}
			}
		}

		for (l = 0; l < SAMPLENUM; l++) {
			w[l] = exp(lw[l] - maxLw);
			wSum = wSum + w[l];
		}
		rowsample(rows, w);

		//reorganize the vector with only keeping the rows in 'rows'
		for (k = 0; k <= i; k++) {
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
				((myStateVec + k)->vpptuple + j)->predator0 =
						((tempStateVec + k)->vpptuple + rows[j])->predator0;
				((myStateVec + k)->vpptuple + j)->predator1 =
						((tempStateVec + k)->vpptuple + rows[j])->predator1;
				((myStateVec + k)->vpptuple + j)->prey =
						((tempStateVec + k)->vpptuple + rows[j])->prey;
				//printf("row %i %i,%i,%i,%i,%i\n", i,k,j,rows[j],((myStateVec+k)->vpptuple+j)->prey,((tempStateVec+k)->vpptuple+rows[j])->prey);
			}
			for (j = 0; j < SAMPLENUM; j++) {
				((tempStateVec + k)->vpptuple + j)->predator0 =
						((myStateVec + k)->vpptuple + j)->predator0;
				((tempStateVec + k)->vpptuple + j)->predator1 =
						((myStateVec + k)->vpptuple + j)->predator1;
				((tempStateVec + k)->vpptuple + j)->prey =
						((myStateVec + k)->vpptuple + j)->prey;
			}
		}
		for (j = 0; j < SAMPLENUM; j++) {
			simData[j].predator0 = ((myStateVec + i)->vpptuple + j)->predator0;
			simData[j].predator1 = ((myStateVec + i)->vpptuple + j)->predator1;
			simData[j].prey = ((myStateVec + i)->vpptuple + j)->prey;
		}
		likehood = likehood + maxLw + log(wSum / SAMPLENUM);
		//printf("likehood: %f, maxLW: %f, mean: %f, wSum: %f\n", likehood, maxLw, log(wSum/SAMPLENUM), wSum);
		*ll = likehood;
		rows[SAMPLENUM] = NULL;
	}

	for (i = 0; i < obsnum; i++) {
		ppstate[i].prey = ((myStateVec + i)->vpptuple + 0)->prey;
		ppstate[i].predator0 = ((myStateVec + i)->vpptuple + 0)->predator0;
		ppstate[i].predator1 = ((myStateVec + i)->vpptuple + 0)->predator1;
	}

	freeMem:
	//printf("free mem...");
	if (myStateVec != NULL) {
		//printf("free mem 2...");
		for (i = 0; ((myStateVec + i)->vpptuple != NULL) && i < obsnum; i++) {
			//printf("free myStateVec...i %i\n", i);
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
		free(simData);
		simData = NULL;
	}

}

void stepLV(lvstate *state, double *t0p, double *dtp, double *lvParam) {
	//printf("Starting stepLV...\n");
	register double t0 = *t0p, dt = *dtp, t;
	register double h0, h1, h2, h3, h4, h5, u, l0, l1, l2, l3, l4;
	register int prey, predator0, predator1;
	prey = state->prey;
	predator0 = state->predator0;
	predator1 = state->predator1;
	l0 = lvParam[0];
	l1 = lvParam[1];
	l2 = lvParam[2];
	l3 = lvParam[3];
	l4 = lvParam[4];

	while (dt > 0) {
		h1 = l0 * prey;
		h2 = l1 * prey * predator0;
		h3 = l2 * prey * predator1;
		h4 = l3 * predator0;
		h5 = l4 * predator1;
		h0 = h1 + h2 + h3 + h4 + h5;

		if ((h0 < (1e-10)) || (prey >= 1000000))
			t = 1e99;
		else {
			t = gsl_ran_exponential(r, 1/h0);
		}
		if (t > dt) {
			state->prey = prey;
			state->predator0 = predator0;
			state->predator1 = predator1;
			return;
		} else {
			u = gsl_rng_uniform(r);
			if (u < (h1 / h0)) {
				prey = prey + 1;
			} else if (u < ((h1 + h2) / h0)) {
				prey = prey - 1;
				predator0 = predator0 + 1;
			} else if (u < ((h1 + h2 + h3) / h0)) {
				prey = prey - 1;
				predator1 = predator1 + 1;
			} else if (u < ((h1 + h2 + h3 + h4) / h0)) {
				predator0 = predator0 - 1;
			} else {
				predator1 = predator1 - 1;
			}
			dt = dt - t;

		}
	}

	return;
}

//mcmc process
void runPmmhPath(int its, double *lvParam, double *obslik, lvstate ppstate[16],
		obsdata *myobsdata) {
	//printf("Starting runPmmhPath...\n");
	int i, j;
	double propParam[5];
	double curParam[5];
	double ll;
	double propMll, curMll = -1e99;
	lvstate curPath[obsnum];
	lvstate propPath[obsnum];

	FILE *fp = fopen("mcmc-out.txt", "wab+");
	if (fp == NULL) {
		printf("file is null.\n");
		return;
	}
	fprintf(fp, "th1,th2,th3,th4,th5,");
	for (i = 0; i <= 30; i++) {
		if ((i % 2) == 0) {
			if (i == 30) {
				fprintf(fp, "x%i,y0_%i,y1_%i", i, i, i);
			} else {
				fprintf(fp, "x%i,y0_%i,y1_%i,", i, i, i);
			}
		}
	}
	fprintf(fp, "\n");
	memcpy(curPath, ppstate, sizeof(lvstate) * obsnum);
	memcpy(curParam, lvParam, sizeof(double) * 5);

	for (i = its; i > 0; i--) {
		//calculate measurements and decide moving to next state or not
		memcpy(propParam, curParam, sizeof(double) * 5);
		peturb(propParam);

		pfPropPar(propPath, myobsdata, propParam, &ll);
		propMll = ll;
		if (log(gsl_ran_flat(r, 0.0, 1.0)) < (propMll - curMll)) {
			curMll = propMll;
			memcpy(curParam, propParam, sizeof(double) * 5);
			memcpy(curPath, propPath, sizeof(lvstate) * obsnum);
		}
		//write current parameters and state sequence into output
		printf("run for the %ith iteration.\n", i);
		//write current parameters and state sequence into output
		fprintf(fp, "%f,%f,%f,%f,%f,", curParam[0], curParam[1], curParam[2], curParam[3], curParam[4]);

		for (j = 0; j < obsnum; j++) {
			if (j == (obsnum - 1)) {
				fprintf(fp, "%i,%i,%i", curPath[j].prey, curPath[j].predator0, curPath[j].predator1);
			} else {
				fprintf(fp, "%i,%i,%i,", curPath[j].prey, curPath[j].predator0, curPath[j].predator1);
			}
		}
		fprintf(fp, "\n");
		printf("end for the %ith iteration.\n", i);
	}
	fclose(fp);
}

void rowsample(int *rows, double *w) {
	//printf("Starting rowsample...\n");
	gsl_ran_discrete_t * grdp;
	int row;
	int i;
	grdp = gsl_ran_discrete_preproc(SAMPLENUM, w);
	for (i = 0; i < SAMPLENUM; i++) {

		row = (int) gsl_ran_discrete(r, grdp);
		rows[i] = row;

	}
	gsl_ran_discrete_free(grdp);
	return;
}

void peturb(double *lvParam) {
	//printf("Starting peturb...\n");

	const double SIGMA = 0.035;
	//printf("lvParam_before %f, %f, %f\n", lvParam[0], lvParam[1], lvParam[2]);
	lvParam[0] = lvParam[0] * exp(gsl_ran_gaussian(r, SIGMA));
	lvParam[1] = lvParam[1] * exp(gsl_ran_gaussian(r, SIGMA));
	lvParam[2] = lvParam[2] * exp(gsl_ran_gaussian(r, SIGMA));
	lvParam[3] = lvParam[3] * exp(gsl_ran_gaussian(r, SIGMA));
	lvParam[4] = lvParam[4] * exp(gsl_ran_gaussian(r, SIGMA));
	//printf("lvParam_after  %f, %f, %f\n", lvParam[0], lvParam[1], lvParam[2]);
}

//to get the prior simulate value of prey and predator
void simPrior(lvstate *simData) {
	//printf("Starting simPrior...\n");
	int i;

	const double PREY_MEAN = 800.0;
	const double PRED0_MEAN = 500.0;
	const double PRED1_MEAN = 600.0;

	for (i = 0; i < SAMPLENUM; i++) {
		simData[i].prey = gsl_ran_poisson(r, PREY_MEAN);
		simData[i].predator0 = gsl_ran_poisson(r, PRED0_MEAN);
		simData[i].predator1 = gsl_ran_poisson(r, PRED1_MEAN);
		//printf("simPrior %i\n", simData[i].predator);
	}
	//printf("End simPrior...\n");
	return;

}

void runModel(int its) {
	obsnum = 16;
	obsdata myobsdata[obsnum];
	lvstate mylvstate[obsnum];
	double ll = 0.0;

	double lvParam[5] = { 10.0, 0.005, 0.0025, 6.0, 3.0 };

	//produce time list(with even number in)
	int i;
	int j = 0;

	for (i = 0; i <= 30; i++) {
		if ((i % 2) == 0) {
			myobsdata[j].time = (double) i;
			j++;
		}
	}

	//read external txt file
	FILE *file = fopen("LV12.txt", "r");
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
	runPmmhPath(its, lvParam, &ll, mylvstate, myobsdata);
	return;
}

int main(int argc, char *argv[]) {
	const gsl_rng_type * T;
	long seed;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	seed = time(NULL) * getpid();    // set seed
	gsl_rng_set(r, seed);
	clock_t begin, end;
	double time_spent;

	begin = clock();
	int its;
	printf("Starting main...\n");

	if (argc == 1) {
		its = 5;
	} else {
		its = atoi(argv[1]);
	}
	runModel(its);

	//double ll = 0.0;
	//double gp = 0.0;
	//double test = 0.0;
	//test = 0.005;
	//gp = gsl_ran_gaussian_pdf((411 - 810), 20.0);
	//printf("gp: %f\n", ll);
	//ll = log(gsl_ran_gaussian_pdf((411 - 810), 20.0));
	//printf("likelihood: %f\n", ll+test);

	gsl_rng_free(r);

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time consuming in total is %f\n", time_spent);

	printf("Running for %i", its);
	printf(" iterations, Done.");
	return 0;
}
*/
