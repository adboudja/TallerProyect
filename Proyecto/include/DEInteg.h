//
// Created by Adam on 25/05/2024.
//

#ifndef PROYECTO_DEINTEG_H
#define PROYECTO_DEINTEG_H


double* DEInteg(double* (*func)(double,double *),double t,double tout,double relerr,double abserr,int n_eqn,double* y);


#endif //PROYECTO_DEINTEG_H
