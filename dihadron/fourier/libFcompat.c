// Last modified on August 16th, 2017 -- F. Gelis


/* 
   This file contains wrapper functions that are callable directly
   from a FORTRAN program -- Mainly, it's a rewrite/translation of the
   data structures so that they comply with the FORTRAN argument
   passing conventions.

   The 'makefile' gives details regarding the compilation. Note that
   it may not work with combinations of compilers other than 'gcc' and
   'g77/gfortran' due to non trivial differences in the binary layout
   of code produced by different compilers. If you have 'gcc' and
   'g77/gfortran', just type 'make' and it should work...

*/


#include <gsl/gsl_errno.h>
#include "fourier.h"


void init_workspace_fourier_(int *N){
  init_workspace_fourier(*N);
}


void set_fourier_precision_(double *e,double *e1){
  set_fourier_precision(*e,*e1);
}


void set_fpu_state_(){
  set_fpu_state();
}

void gsl_set_error_handler_off_(){
  gsl_set_error_handler_off();
}



double fourier_j0_(double *x,double (*func)(double*,void*),double* param){
  return fourier_j0(*x,func,(void *)param);
}

double fourier_j0_i_(double *x,int *Ni,double (*func)(double*,void*),double* param){
  return fourier_j0_i(*x,*Ni,func,(void *)param);
}

double fourier_j0_f_(double *x,int *Nf,double (*func)(double*,void*),double* param){
  return fourier_j0_f(*x,*Nf,func,(void *)param);
}

double fourier_j0_if_(double *x,int *Ni,int *Nf,double (*func)(double*,void*),double* param){
  return fourier_j0_if(*x,*Ni,*Nf,func,(void *)param);
}

double fourier_j1_(double *x,double (*func)(double*,void*),double* param){
  return fourier_j1(*x,func,(void *)param);
}

double fourier_j1_i_(double *x,int *Ni,double (*func)(double*,void*),double* param){
  return fourier_j1_i(*x,*Ni,func,(void *)param);
}

double fourier_j1_f_(double *x,int *Nf,double (*func)(double*,void*),double* param){
  return fourier_j1_f(*x,*Nf,func,(void *)param);
}

double fourier_j1_if_(double *x,int *Ni,int *Nf,double (*func)(double*,void*),double* param){
  return fourier_j1_if(*x,*Ni,*Nf,func,(void *)param);
}
