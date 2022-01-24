// Last modified on August 16th, 2017 -- F. Gelis


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>
#include <gsl/gsl_errno.h>

#include "fourier.h"


// Example of function to be Fourier transformed. Extra parameters are
// passed as a pointer to void. A proper cast must be performed before
// using the pointer.

double F(double k,void *params){
  double *p=(double *)params;

  return p[0];
}


// The result of the Fourier transform, which is known in closed form
// in this case. For testing purposes...

double G(double x,void *params){
  double *p=(double *)params;

  return p[0]/x;
}


// Main program

int main(int argc, char *argv[]){
  double x,y;
  double params[1];


  // Some initialisation stuff -- in principle nothing should be
  // changed here (1000 is the number of zeros of J0(x) that have been
  // encoded in the table in the file fourier.c)

  //set_fpu_state();
  gsl_set_error_handler_off();
  init_workspace_fourier(1000);
  set_fourier_precision(1.0e-12,1.0e-12);


  params[0]=1.0;
  x=0.02;

  do {
    y=fourier_j0(x,F,(void *)params);
    fprintf(stdout,"%.16e %.16e %.16e\n",x,y,G(x,params));
    x*=1.02;
  } while(x<10.0);

  exit(0);
}
