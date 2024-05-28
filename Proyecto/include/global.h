//
// Created by adboudja on 24/04/2024.
//

#ifndef UNTITLED_GLOBAL_H
#define UNTITLED_GLOBAL_H


#include "Matrix.h"

class global {
public:
    static Matrix *eopdate;
    static double* *geos3;
    static Matrix *Cnm;
    static Matrix *Snm;
    static Matrix *temp;
    static Matrix *PC;
    static double Mjd_UTC;
    static double Mjd_TT;
    static int n;
    static int m;
    static int sun;
    static int moon;
    static int planets;
    static void eop19620101();
    static void GGM03S();
    static void GEOS3(int nobs);
    static void DE430Coeff();
    static void AuxParam();

};


#endif //UNTITLED_GLOBAL_H
