#ifndef DATA_H
#define DATA_H

#include <gsl/gsl_vector.h>

//Define the struct of a single label
typedef struct {
	int statementId;
	int labelerId;
	int label;
} Label;

//Define the struct of the entire dataset
typedef struct {
	Label *labels;
	int numLabels;
	int numLabelers;
	int numStatements;
	double *probZ1, *probZ0;
	double * S, *C, *P;
	double A;
	double *priorZ1;
} Dataset;

//Functions for managing the Dataset
void pack (gsl_vector *x, Dataset *data);
void unpack (const gsl_vector *x, Dataset *data);
void clearData(Dataset *data);
void readData (char *filename, Dataset *data);
void outputResults (Dataset *data);
//void saveResults (Dataset *max_data,Dataset *data);

#endif
