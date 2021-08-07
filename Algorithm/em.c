#include <gsl/gsl_vector.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "prob_functions.h"
#include <float.h>
double * c_array;
double * s_array;

float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}


double EM (Dataset *data, int index)
{
	int i, j;

	//Threshold to terminate the EM Algorithm
	const double THRESHOLD = 1E-5;
	double Q_value, lastQ;

	//Initialize all the parameters from the same starting point since we don't have any prior information about them.
	for (i = 0; i < data->numLabelers; i++)
  {
		data->P[i] = 0.5;
    if (index == 0)
    {
		    data->C[i] = 0.5;
    }else if (index == 1)
    {
       data->C[i] = 1;
    } else
    {
      data->C[i] = c_array[i];
    }
	}


  //Initialize the Statement Stances
  for (j = 0; j < data->numStatements; j++)
  {
      if (index == 0)
      {
			     data->S[j] = 1;
      } else if (index == 1)
      {
        data->S[j] = 0.5;
      }else
      {
        data->S[j] = s_array[j];
      }
	}

  //Initialize the symmetry coefficient
	data->A = -0.69;

	//Initialize the Hidden Variable Prior
	for (j = 0; j < data->numStatements; j++) {
		data->priorZ1[j] = 0.5;
	}

	//Start the EM Algorithm
	Q_value = 0;

	int counter = 0;
  double saveL = -INFINITY;
	do {

    //Starter for the EM Algorithm
    if (counter == 0)
    {
      Expectation(data);
    	Q_value = Q_computer(data);
    }
		lastQ = Q_value;


		Expectation(data);
		Q_value = Q_computer(data);

		Maximization(data);

		Q_value = Q_computer(data);

		//Make sure the Likelihood doesn't fluctuate, if it does then terminate the algorithm
    double epsilon = 0.000001f;
    if (Likelihood(data) < saveL && fabs(Likelihood(data) - saveL) > epsilon)
    {
      printf("Likelihood Fluct! L %f save %f\n",Likelihood(data), saveL );
			abort();

    }

    saveL = Likelihood(data);
		counter++;

	} while (fabs((Q_value - lastQ)/lastQ) > THRESHOLD && counter < 1000);

  return Likelihood(data);
}

int main (int argc, char *argv[])
{
  //Initialize the dataset storage
	Dataset * data =  (Dataset *) malloc(sizeof(Dataset));

	if (argc < 2) {
		fprintf(stdout, "Usage format is \"./em data_file\" \n");
		exit(1);
	}


  //Number of EM Rounds Selection
	int em_rounds = 3;
  int j;
  //int save_ind = 0;
  //Maximum Likelihood Found in the multiple EM Iteration
	double max_likelihood = -(RAND_MAX - 1);
	int i = 0;

  //For Loop That will manage the number of EM Rounds
	for(i=0; i < em_rounds; i++)
	{
    
  	//Read the Data
  	readData(argv[1], data);

  	//Run the EM Algorithm
  	double likelihood = EM(data,i);

  	//Error Condition, make sure Likelihood is not Infinity or nan
  	if ((likelihood == INFINITY || likelihood == -INFINITY) || isnan(likelihood) ){
  		printf("max_likelihood is inf or nan\n");
  		abort();

  	}

    if (i == 0)
    {
      c_array = (double *) malloc(sizeof(double)*data->numLabelers);
      for (j = 0; j < data->numLabelers; j++)
      {
        c_array[j] = data->C[j];
      }
    }
    if (i == 1)
    {
      s_array = (double *) malloc(sizeof(double)*data->numStatements);
      for (j = 0; j < data->numStatements; j++)
      {
        s_array[j] = data->S[j];
      }
    }
    // printf("index\n");
    // for (j = 0; j < data->numStatements; j++) {
    //
  	// 	printf("%f\n",data->probZ1[j]);}
  	//Output the calculated results
  	if (likelihood > max_likelihood)
    {
      max_likelihood =  likelihood;
      // save_ind = i;
      outputResults (data);

    }
  }
  // printf("%d\n",max_likelihood);
  // FILE *fp;
	// fp = fopen("./ground_truth_parameters.txt","rt");
  // fscanf(fp, "%d %d\n", &data->numLabelers, &data->numStatements);
  //
	// for (i = 0; i < data->numLabelers; i++) {
  //   fscanf(fp, "%lf\n", &data->C[i]);
	// 	fscanf(fp, "%lf\n", &data->P[i]);
	// }
	// for (j = 0; j < data->numStatements; j++) {
	// 	fscanf(fp, "%lf\n", &data->S[j]);
	// }
  //
	// fscanf(fp, "%lf\n", &data->A);
	// fclose(fp);
  // printf("Likelihood 2: %f\n", Likelihood(data));
  // FILE *fptr;
  // fptr = fopen("./freq.txt","a");
  // fprintf(fptr,"%d\n", save_ind);
  // fclose(fptr);

}
