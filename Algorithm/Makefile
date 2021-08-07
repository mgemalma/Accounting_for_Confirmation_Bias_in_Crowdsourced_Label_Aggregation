default: em

clean:
	rm -f *.o em

# Set these variables to correct locations
CBLAS_LIB=/Users/altuggemalmaz/Desktop/gsl-2.6/cblas/.libs
GSL_LIB=/Users/altuggemalmaz/Desktop/gsl-2.6/.libs
GSL_INCLUDE=/Users/altuggemalmaz/Desktop/gsl-2.6

em: em.c data.c data.h prob_functions.c prob_functions.h
	gcc -m64 -O3 $(GSL_LIB)/libgsl.a $(CBLAS_LIB)/libgslcblas.a -I$(GSL_INCLUDE) -std=c99 -o em em.c data.c prob_functions.c
