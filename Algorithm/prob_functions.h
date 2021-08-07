#ifndef PROB_FUNCTIONS_H
#define PROB_FUNCTIONS_H

#include <gsl/gsl_vector.h>
#include "data.h"

//required gsl functions
double my_f (const gsl_vector *x, void *params);
void my_df (const gsl_vector *x, void *params, gsl_vector *g);
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *g);
//function that calculates gradient direction for every parameter
void param_gradients(Dataset *data, double *gradient_c, double *gradient_p, double * gradient_s, double * gradient_a);
//likelihood and Q value calculator
double Likelihood (Dataset *data);
double Q_computer (Dataset *data);
//log probability calculator
double Prob(int l, int z, double c, double p ,double a, double s, Dataset * data, int label);
void Expectation (Dataset *data);
void Maximization (Dataset *data);

#endif
