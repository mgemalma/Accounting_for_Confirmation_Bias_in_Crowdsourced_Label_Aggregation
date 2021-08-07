#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_erf.h>
#include "prob_functions.h"
#include "data.h"

//Error rate for the denominator of the gradients in case they become 0 (optional parameter value, if 0 then unused
double error_rate = 0.00000000000000001;

//Define the required functions for the optimizer,
//this one puts back the calculated parameters from optimizer to the dataset
double my_f(const gsl_vector *x, void *params)
{
	//Get the Dataset
	Dataset *data = (Dataset *) params;

	//Put the parameter values into the dataset.
	unpack(x, data);

	//We are using gradient descent, since we want to maximize we return negative of Q function.
	return - Q_computer(data);
}

//Calculate the gradients for this function
void my_df(const gsl_vector *x, void *params, gsl_vector *g)
{
	int i, j;

	//Get the dataset
	Dataset *data = (Dataset *) params;

	//For each parameter, we will have to have gradient direction so that we can update the parameters,
	//thus for each parameter we have an array entry to do so.
	//Initialize the gradient parameter arrays
	double *gradient_c = (double *) malloc(sizeof(double) * data->numLabelers);
	double *gradient_p = (double *) malloc(sizeof(double) * data->numLabelers);
	double *gradient_s = (double *) malloc(sizeof(double) * data->numStatements);
	double * gradient_a = (double *) malloc(sizeof(double));

	//Make sure the memory location is actually zero
	for (i = 0; i < data->numLabelers; i++) {
		gradient_c[i] =  0;
		gradient_p[i] =  0;
	}
	for (j = 0; j < data->numStatements; j++) {
		gradient_s[j] =  0;
	}
	gradient_a[0] =  0;

	//Make sure we have the latest of the parameters
	unpack(x, data);

	//Calculate the gradient for all the parameters
	param_gradients(data, gradient_c, gradient_p, gradient_s, gradient_a );

	//Since we are using gradient descent, flip the sign since we want to minimize
	for (i = 0; i < data->numLabelers; i++) {
		gsl_vector_set(g, i, - gradient_c[i]);

	}
	int save_i = i;

	for (i = 0; i < data->numLabelers; i++)
	{
		gsl_vector_set(g, save_i + i, - gradient_p[i]);
	}
	save_i += i;
	for (j = 0; j < data->numStatements; j++) {
		gsl_vector_set(g, save_i + j, - gradient_s[j]);
	}
  save_i += j;

	gsl_vector_set(g, save_i, -gradient_a[0]);


	free(gradient_c);
	free(gradient_s);
	free(gradient_p);
}

void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *g)
{
	//unpack the params into the dataset and have the negative Q value since we are using gradient descent.
	//my_f returns -Q
	*f = my_f(x, params);
	//calculate the gradients
	my_df(x, params, g);

}

//Gradient calculator function for every parameter
void param_gradients(Dataset *data, double *gradient_c, double *gradient_p, double * gradient_s, double * gradient_a)
{
	int i, j, label_id;

	//Go through each of the labels
	for (label_id = 0; label_id < data->numLabels; label_id++) {
		int i = data->labels[label_id].labelerId;
		int j = data->labels[label_id].statementId;
		int lij = data->labels[label_id].label;

		//if the label is 0 then gradients are calculated accordingly
		if (lij == 0)
		{
			//Get the parameters into the range [0,1] through the function |sin(x)|
			double p = fabs(sin(data->P[i]));
			double c = fabs(sin(data->C[i]));
			double s = fabs(sin(data->S[j]));
			double a = exp(data->A);
			int z = 0;
			//Calculate the gradient direction for statement stance
			double gradS = (-2*(1-p)*(s- c  )/((double) a  ));
			//Calculate the gradient direction for the current stance
			double gradC = (2*(1- p)*(s- c  )/((double)a));
			//Calculate the gradient direction for the skill level with z value 0
			double gradP0 = -(-pow(s, 2) + 2*c*s - pow(c, 2) + z) / ((double)a);
			z = 1;
			//Calculate the gradient direction for the skill level with z value 1
			double gradP1 = -(-pow(s, 2) + 2*c*s - pow(c, 2) + z) / ((double)a);
			z = 0;
			//Calculate the gradient direction for the symmetry coefficient with z value 0
			double gradA0 = ((1-p)*(pow(s-c, 2)) +  p*z) / (pow(a,2));
			z = 1;
			//Calculate the gradient direction for the symmetry coefficient with z value 1
			double gradA1 =((1-p)*(pow(s-c, 2)) +  p*z) / (pow(a,2));

			//add the gradient direction into the array
			gradient_s[j] += data->probZ1[j]*gradS + data->probZ0[j]*gradS;

			//Make sure that it is not nan
			if (isnan(gradient_s[j]))
			{
				printf("Gradient Function Failed\n");
				abort();
			}
			//add the gradient direction into the array
			gradient_c[i] += data->probZ1[j]*gradC + data->probZ0[j]*gradC;
			//Make sure that it is not nan
			if (isnan(gradient_c[i]))
			{
				printf("Gradient Function Failed\n");
				abort();
			}
			//add the gradient direction into the array
			gradient_p[i] += data->probZ1[j]*gradP1 + data->probZ0[j]*gradP0;
			//Make sure that it is not nan
			if (isnan(gradient_p[i]))
			{
				printf("Gradient Function Failed\n");
				abort();
			}
			//add the gradient direction into the array
			gradient_a[0] += data->probZ1[j]*gradA1 + data->probZ0[j]*gradA0;
			//Make sure that it is not nan
			if (isnan(gradient_a[0]))
			{
				printf("Gradient Function Failed\n");
				abort();
			}

		}else
		{
				//Get the parameters into the range [0,1] through the function |sin(x)|
				double p = fabs(sin(data->P[i]));
				double c = fabs(sin(data->C[i]));
				double s = fabs(sin(data->S[j]));
				double a = exp(data->A);

				//Calculate the gradient direction for statement stance with z = 0
				int z = 0;
				double gradS0 = (2*(1-p)*(s-c)) / (a*(exp(((1-p)*(pow(s-c, 2))  + p*z) / ((double)a))-1 + error_rate));
				//Calculate the gradient direction for statement stance with z = 1
				z = 1;
				double gradS1 = (2*(1-p)*(s-c)) / (a*(exp(((1-p)*(pow(s-c, 2))  + p*z) / ((double)a))-1 + error_rate));
				//Calculate the gradient direction for the current stance with z = 0
				z = 0;
				double gradC0  = (-2*(1-p)*(s-c)) / (a*(exp(((1-p)*(pow(s-c, 2))  + p*z) / ((double)a))-1 + error_rate));
				//Calculate the gradient direction for the current stance with z = 1
				z = 1;
				double gradC1  = (-2*(1-p)*(s-c)) / (a*(exp(((1-p)*(pow(s-c, 2))  + p*z) / ((double)a))-1 + error_rate));
				//Calculate the gradient direction for the skill level with z = 0
				z = 0;
				double gradP0 = (-(pow(s,2)) + 2*c*s - (pow(c, 2)) + z) / (((double)a)*(exp(((1-p)*(pow(s-c, 2))  + p*z) / ((double)a) )-1 + error_rate));
				//Calculate the gradient direction for the skill level with z = 1
				z = 1;
				double gradP1 = (-(pow(s,2)) + 2*c*s - (pow(c, 2)) + z) / (((double)a)*(exp(((1-p)*(pow(s-c, 2))  + p*z) / ((double)a) )-1 + error_rate));
				//Calculate the gradient direction for the symmetry coefficient with z value 0
				z = 0;
				double gradA0 = -(p*z + (-p + 1)*(pow(s-c, 2)) ) / ((pow(a,2))*( exp(   ((1-p)*(pow(s-c, 2))  + p*z) / ((double)a)   ) - 1  + error_rate ));
				//Calculate the gradient direction for the symmetry coefficient with z value 1
				z = 1;
				double gradA1 = -(p*z + (-p + 1)*(pow(s-c, 2)) ) / ((pow(a,2))*( exp(   ((1-p)*(pow(s-c, 2))  + p*z) / ((double)a)   ) - 1  + error_rate ));

				//add the gradient direction into the array
				gradient_s[j] += data->probZ1[j]*gradS1 + data->probZ0[j]*gradS0;
				//Make sure that it is not nan
				if (isnan(gradient_s[j]))
				{
					printf("Gradient Function Failed\n");
					abort();
				}

				//add the gradient direction into the array
				gradient_c[i] += data->probZ1[j]*gradC1 + data->probZ0[j]*gradC0;
				//Make sure that it is not nan
				if (isnan(gradient_c[i]))
				{
					printf("Gradient Function Failed\n");
					abort();
				}

				//add the gradient direction into the array
				gradient_p[i] += data->probZ1[j]*gradP1 + data->probZ0[j]*gradP0;
				//Make sure that it is not nan
				if (isnan(gradient_p[i]))
				{
					printf("Gradient Function Failed\n");
					abort();
				}

				//add the gradient direction into the array
				gradient_a[0] += data->probZ1[j]*gradA1 + data->probZ0[j]*gradA0;
				//Make sure that it is not nan
				if (isnan(gradient_a[0]))
				{
					printf("Gradient Function Failed\n");
					abort();
				}

		}
	}
}

//Helper function to compare two floats
int cmpf(float A, float B)
{
	float epsilon = 0.000001f;
    return (fabs(A - B) < epsilon);
}


//Function to compute likelihood
double Likelihood(Dataset *data)
{

	int i, j, label_id;

	//Initialize the parameters
	double Likelihood = 0;

	//Access the parameters
	double *c = data->C, *s = data->S, *p = data->P;
	double a = data->A;

	//For all of the statements
	for (j = 0; j < data->numStatements; j++) {

		//Initialize the multiplications for both Z probabilities
		double P1 = data->priorZ1[j];
		double P0 = (1 - data->priorZ1[j]);

		//Go through all of the labels
		for (label_id = 0; label_id < data->numLabels; label_id++) {

			//Get the statement information
			int i = data->labels[label_id].labelerId;
			int sid = data->labels[label_id].statementId;
			int lij = data->labels[label_id].label;

			//if the statement id is the same as the one for which we are calculating the probabilities then;
			if (j == sid){
				double Z1 = 0;
				double Z0 = 0;
				double A = exp(data->A);
				double P = fabs(sin(data->P[i]));
				double C = fabs(sin(data->C[i]));
				double S = fabs(sin(data->S[j]));

				//With respect to the label calculate the probability
				if (lij == 1)
				{
					Z1 = 1 - (1.0 / exp((1.0 / A) *((1 - P)*pow(S - C,2) + P*1)));
					Z0 = 1 - (1.0 / exp((1.0 / A) *((1 - P)*pow(S - C,2) + P*0)));
				}

				if (lij == 0)
				{
					Z1 = (1.0 / exp((1.0 / A) *((1-P)*pow(S - C,2) + P*1)));
					Z0 = (1.0 / exp((1.0 / A) *((1-P)*pow(S - C,2) + P*0)));
				}

				//Multiply with the total probability
				P1 *= Z1;
				P0 *= Z0;
			}
		}
		//Add the statement j likelihood to the total likelihood
		Likelihood += log(P1 + P0);

	}

	return Likelihood;

}

//Function to compute the Q value
double Q_computer(Dataset *data)
{
	int i, j, label_id;

	//Initialize Q value
	double Q_value = 0;

	double *c = data->C, *s = data->S, *p = data->P;
	double a = data->A;

	//For each of the statements, calculate the first summation of the Q function which includes,
	//posterior z probability and prior z probability.
	for (j = 0; j < data->numStatements; j++) {

		double z1 = data->priorZ1[j];

		//Tweak the prior a little bit if it is 0 since log(0) is -inf
		if (cmpf(z1, 0.0) || cmpf(1-z1, 0.0))
		{

			if (cmpf(z1, 0.0) )
			{
				z1 = 0.000001f;
			}else
			{
				z1 = 1 - 0.000001f;
			}
		}

		Q_value += data->probZ1[j] * log(z1);
		Q_value += data->probZ0[j] * log(1 - z1);

		//Make sure at this stage the Q function value is not nan.
		if (isnan(Q_value))
		{
			printf("Z1 %f Z0 %f Z1 calc %f Z0 calc %f\n", data->priorZ1[j], 1 - data->priorZ1[j],  log(data->priorZ1[j]) ,log(1 - data->priorZ1[j]) );
			abort();
		}
	}

	//Go through every labels, do the second summation part of the Q function
	for (label_id = 0; label_id < data->numLabels; label_id++) {

		//Get the label information
		int i = data->labels[label_id].labelerId;
		int j = data->labels[label_id].statementId;
		int lij = data->labels[label_id].label;

		//According to the label value pick the respective probability formula
		double Z1 = 0;
		double Z0 = 0;

		double A = exp(data->A);
		double P = fabs(sin(data->P[i]));
		double C = fabs(sin(data->C[i]));
		double S = fabs(sin(data->S[j]));

		if (lij == 1)
		{
			Z1 = 1 - (1.0 / exp((1.0 / A ) * ((1-P)*pow(S-C,2) + P*1)));
			Z0 = 1 - (1.0 / exp((1.0 / A ) * ((1-P)*pow(S-C,2) + P*0)));
		}

		if (lij == 0)
		{
			Z1 = (1.0 / exp((1.0 / A ) *((1-P)*pow(S-C,2) + P*1)));
			Z0 = (1.0 / exp((1.0 / A ) *((1-P)*pow(S-C,2) + P*0)));
		}

		double old_Q = Q_value;

		//Add the probabilites to the total Q function value
		//Make sure Z0 or Z1 not zero as log(0) is invalid
		if (cmpf(Z1, 0.0) || cmpf(Z0, 0.0))
		{
			Z1 += 1E-10;
			Z0 += 1E-10;
		}
		Q_value += data->probZ1[j] * log(Z1) + data->probZ0[j] * log(Z0);

		//Check the value of Q and make sure it's not nan
		if (isnan(Q_value)) {
			printf("Q is nan A is %f P is %f C is %f S is %f old_Q %f is inner calc %f Z0 %f Z1 %f label %d Z0 %f Z1 %f\n",fabs(sin(data->A)),fabs(sin(data->P[i])),fabs(sin(data->C[i])),fabs(sin(data->S[j])), old_Q, (1.0/exp((1.0 / fabs(sin(data->A)) ) *( (1-fabs(sin(data->P[i])))*pow(fabs(sin(data->S[j]))-fabs(sin(data->C[i])),2) + fabs(sin(data->P[i]))*0))),Z0,Z1,lij, 1.0/exp((1.0 / fabs(sin(data->A)) ) *( (1-fabs(sin(data->P[i])))*pow(fabs(sin(data->S[j]))-fabs(sin(data->C[i])),2) + fabs(sin(data->P[i]))*0)), 1.0/exp((1.0 / fabs(sin(data->A)) ) *( (1-fabs(sin(data->P[i])))*pow(fabs(sin(data->S[j]))-fabs(sin(data->C[i])),2) + fabs(sin(data->P[i]))*1)));
			abort();
		}
	}

	return Q_value;
}

//Log probably calculator with respect to the label, just taking the log of the probability function
double Prob (int l, int z, double c, double p ,double a, double s, Dataset * data,int label )
{
	double prob;
	float epsilon = 0.0001f;

	//If label is 0 then use the 0 label probability formula otherwise use label 1 probability formula
	if (l == 0) {

		prob = 1.0 / (exp((1.0 / a ) * ((1-p)*pow(s-c,2) + p*z)));

	} else {

		prob = (1 - (1.0 / (exp((1.0 / a ) * ((1-p)*pow(s-c,2) + p*z)))));

	}

	//Return the log probability
	return log(prob);
}

//Calculate the E step
void Expectation (Dataset *data)
{
	int j, label_id;

	//Initialize the starting points for the z probabilites, which are the prior probabilites
	for (j = 0; j < data->numStatements; j++) {
		//Use log function to make computations easier
		data->probZ1[j] = log(data->priorZ1[j]) ;
		data->probZ0[j] = log(1 - data->priorZ1[j]) ;
	}

	//Go through every label and compute the respective probabilites
	for (label_id = 0; label_id < data->numLabels; label_id++) {

		//Get the label information
		int i = data->labels[label_id].labelerId;
		int j = data->labels[label_id].statementId;
		int lij = data->labels[label_id].label;

		//Put the parameters into [0,1] range
		//We are scaling the parameter values using the function |sin(x)| to the range [0,1]
		//Labeler Parameters
		double c = fabs(sin(data->C[i]));
		double p = fabs(sin(data->P[i]));

		//Statement Parameters
		double s = fabs(sin(data->S[j]));

		//Symmetry Coefficient
		double a = exp(data->A);

		//Calculate the respective probabilites and add them to the respective statement j probability
		data->probZ1[j] += Prob(lij, 1, c, p, a, s,data, label_id);
		data->probZ0[j] += Prob(lij, 0, c, p, a, s,data, label_id);
	}

	//Normalize the probabilities
	for (j = 0; j < data->numStatements; j++) {

		//For Debugging purposes save the values
		double save = exp(data->probZ1[j]);
		double save1 = exp(data->probZ0[j]);

		//Exponentiate and renormalize
		data->probZ1[j] = exp(data->probZ1[j]);
		data->probZ0[j] = exp(data->probZ0[j]);
		data->probZ1[j] = data->probZ1[j] / (data->probZ1[j] + data->probZ0[j]);
		data->probZ0[j] = 1 - data->probZ1[j];

		//Make sure the probabilites are not nan
		if (isnan(data->probZ1[j]) || isnan(data->probZ0[j])) {
			printf("data->probZ1 is nan! Z1 %f Z0 %f z1 now %f\n",save,save1,exp(data->probZ1[j]));
			abort();
		}

	}
}


//Maximization step
void Maximization(Dataset *data)
{

	//Dataset information
	int numLabelers = data->numLabelers;
	int numStatements = data->numStatements;

	//Gradient descent optimizer library components initialization
	gsl_vector * parameter_array;
	gsl_multimin_fdfminimizer *s;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_function_fdf function_struct;

	//Create a vector entry for all of the parameters that we will optimize
	//skill level and current stance for every labeler + statement stances for every statements +
	//symmetry coefficient
	parameter_array = gsl_vector_alloc(2*numLabelers + numStatements + 1);

	//Put the parameters into optimizer
	pack(parameter_array, data);

	/* Initialize objective function */
	//optimized parameters total
	function_struct.n = 2*numLabelers + numStatements + 1;

	//required functions for gsl optimizer which is defined above
	function_struct.f = &my_f;
	function_struct.df = &my_df;
	function_struct.fdf = &my_fdf;

	//where the parameters are located
	function_struct.params = data;

	//Minimizer
	T = gsl_multimin_fdfminimizer_conjugate_pr;
	s = gsl_multimin_fdfminimizer_alloc(T, 2*numLabelers + numStatements + 1);

	//Last two parameters from left to right: Learning Rate and orthogonality to the search direction requirement
	gsl_multimin_fdfminimizer_set(s, &function_struct, parameter_array, 0.01, 0.01);


	int counter = 0;
	double old_q_val;
	double threshold_val = 0.01;
	int status;
	//Actual gradient descent
	do {
		old_q_val = s->f;
		counter++;

		//Optimize
		status = gsl_multimin_fdfminimizer_iterate(s);
		if (status) {
			break;
		}
		//test gradient
		status = gsl_multimin_test_gradient(s->gradient, 1e-3);

	} while (old_q_val - s->f > threshold_val && status == GSL_CONTINUE && counter < 25);

	//get the calculated parameters into the data
	unpack(s->x, data);

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(parameter_array);

}
