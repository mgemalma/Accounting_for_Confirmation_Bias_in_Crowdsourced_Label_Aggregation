#include <math.h>
#include "data.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

//Put the values of the current values into the optimizer for update (Maximization Step)
void pack (gsl_vector *x, Dataset *data)
{
	int i, j;

	//Put the current value of labeler current stance and labeler skill level
	//for the optimizer
	for (i = 0; i < data->numLabelers; i++) {
		gsl_vector_set(x, i, data->C[i]);
	}

  int save_i = i;
	for (i = 0; i < data->numLabelers; i++)
	{
		gsl_vector_set(x, save_i + i, data->P[i]);
	}
	save_i += i;
	//Put the statement stance into the optimizer
	for (j = 0; j < data->numStatements; j++) {
		gsl_vector_set(x, save_i + j, data->S[j]);
	}
	save_i += j;
	//Put the symmetry coefficient to the optimizer
	gsl_vector_set(x, save_i, data->A);
}


//From the optimizer get calculated parameters
void unpack (const gsl_vector *x, Dataset *data)
{
	int i, j;

	//Get the Curren Stance and the skill level
	//calculations for every labeler from the optimizer.
	for (i = 0; i < data->numLabelers; i++) {
		data->C[i] = gsl_vector_get(x, i);


		//Make sure the parameter values are not defective.
		if (isnan(data->C[i])) {
			printf("Library Failed\n");
			abort();
		}
	}

	int save_i = i;
	for (i = 0; i < data->numLabelers; i++) {
		data->P[i] = gsl_vector_get(x, save_i + i);


		//Make sure the parameter values are not defective.
		if (isnan(data->P[i])) {
			printf("Library Failed\n");
			abort();
		}
	}
	save_i += i;
	//Get the statement stance calculations from the optimizer
	for (j = 0; j < data->numStatements; j++) {
		data->S[j] = gsl_vector_get(x, save_i + j);

		//Make sure the parameter values are not defective.
		if (isnan(data->S[j])) {
			printf("Library Failed\n");
			abort();
		}
	}
	save_i += j;
	//Get the calculated symmetry coefficient A from the optimizer
	data->A = gsl_vector_get(x, save_i);
	if (isnan(data->A)) {
		printf("Library Failed\n");
		abort();
	}


	double prob_new = 0;
	for (j = 0; j < data->numStatements; j++)
	{
			prob_new += data->probZ1[j];

	}
	prob_new /= (double) data->numStatements;

	for (j = 0; j < data->numStatements; j++)
	{
		data->priorZ1[j] = prob_new;
	}


}

//Function to clear the dataset
void clearData(Dataset *data)
{
	int i,j;
	for (j = 0; j < data->numStatements; j++) {
		data->probZ1[j] = 0;
		data->probZ0[j] = 0;
		data->priorZ1[j] = 0;
		data->S[j] = 0;
	}

	for (i = 0; i < data->numLabelers; i++)
	{
			data->C[i] = 0;
			data->P[i] = 0;
	}
	data->A = 0;
}


//Function to read the data in, read the labels of the labelers
void readData (char *filename, Dataset *data)
{
	int i, j, label_id;

	//Open the Dataset File
	FILE *fp = fopen(filename, "rt");

	//Scan the number of items
	fscanf(fp, "%d %d %d\n", &data->numLabels, &data->numLabelers, &data->numStatements);

	//Initilize the priors
	data->priorZ1 = (double *) malloc(sizeof(double) * data->numStatements);
	data->probZ1 = (double *) malloc(sizeof(double) * data->numStatements);
	data->probZ0 = (double *) malloc(sizeof(double) * data->numStatements);
	data->S = (double *) malloc(sizeof(double) * data->numStatements);
	data->C = (double *) malloc(sizeof(double) * data->numLabelers);
	data->P = (double *) malloc(sizeof(double) * data->numLabelers);
	data->labels = (Label *) malloc(sizeof(Label) * data->numLabels);


	//Read labels
	label_id = 0;
	while (fscanf(fp, "%d %d %d\n", &(data->labels[label_id].statementId),
		      &(data->labels[label_id].labelerId), &(data->labels[label_id].label)) == 3) {
		label_id++;
	}
	fclose(fp);
}

/*
	Function to save the results
*/
void outputResults (Dataset *data)
{
	int i, j;

	//Save the parameters calculated for the labelers.
	FILE *fptr;
	fptr = fopen("./worker_parameter_results.txt","w");

	for (i = 0; i < data->numLabelers; i++) {
		fprintf(fptr,"%f\n",  fabs(sin(data->C[i])));
		fprintf(fptr,"%f\n",  fabs(sin(data->P[i])));
	}
	for (j = 0; j < data->numStatements; j++) {
		fprintf(fptr,"%f\n", fabs(sin(data->S[j])));
	}

	fprintf(fptr,"%f\n", exp(data->A));
	fclose(fptr);

	//Save the labels calculated for the statements, take a threshold of 0.5 to make a decision.
	FILE *fptr1;
	fptr1 = fopen("./results.txt","w");

	for (j = 0; j < data->numStatements; j++) {
		//printf("%f\n",data->probZ1[j]);
		if( data->probZ1[j] > 0.5)
		{
				fprintf(fptr1,"%d\n",1);
		}else if (data->probZ0[j] > 0.5)
		{
			fprintf(fptr1,"%d\n",0);
		} else
		{
			//Set the Randomization Seed
			srand(time(0));
			//Make a Random Prediction
			fprintf(fptr1,"%d\n",(((double)rand() / (double)RAND_MAX) > 0.5) ? 1 : 0);
		}
	}

	fclose(fptr1);

}

//In case want to use multiple EM Rounds

// void saveResults (Dataset *max_data,Dataset *data)
// {
// 	int i, j;
//
// 	max_data->probZ1 = (double *) malloc(sizeof(double) * data->numStatements);
// 	max_data->probZ0 = (double *) malloc(sizeof(double) * data->numStatements);
// 	max_data->S = (double *) malloc(sizeof(double) * data->numStatements);
// 	max_data->C = (double *) malloc(sizeof(double) * data->numLabelers);
// 	max_data->P = (double *) malloc(sizeof(double) * data->numLabelers);
// 	max_data->priorC = (double *) malloc(sizeof(double) * data->numLabelers);
// 	max_data->priorP = (double *) malloc(sizeof(double) * data->numLabelers);
// 	max_data->priorS = (double *) malloc(sizeof(double) * data->numStatements);
// 	max_data->priorZ1 = (double *) malloc(sizeof(double) * data->numStatements);
// 	max_data->numLabelers = data->numLabelers;
// 	max_data->numStatements = data->numStatements;
// 	max_data->A = data->A;
//
//
// 	for (i = 0; i < data->numLabelers; i++) {
// 		max_data->C[i] = data->C[i];
// 		max_data->P[i] = data->P[i];
// 	}
// 	for (j = 0; j < data->numStatements; j++) {
// 		max_data->S[j] = data->S[j];
// 	}
// 	for (j = 0; j < data->numStatements; j++) {
// 		max_data->probZ1[j] = data->probZ1[j];
// 	}
//
// }
