//
// Created by Adam on 25/05/2024.
//

#include "VarEqn.h"
#include "global.h"
#include "IERS.h"
#include "timediff.h"
#include "SAT_Const.h"
#include "NutMatrix.h"
#include "PrecMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "auxFunc.h"

/*%------------------------------------------------------------------------------
%
% VarEqn.m
%
% Purpose:
%   Computes the variational equations, i.e. the derivative of the state vector
%   and the state transition matrix
%
% Input:
%   x           Time since epoch in [s]
%   yPhi        (6+36)-dim vector comprising the state vector (y) and the
%               state transition matrix (Phi) in column wise storage order
%
% Output:
%   yPhip       Derivative of yPhi
%
% Last modified:   2015/08/12   M. Mahooti
%
%------------------------------------------------------------------------------*/
/**
 * @file VarEqn.cpp
 * @brief Implementación de las ecuaciones variacionales.
 *
 * Este archivo contiene la implementación de las ecuaciones variacionales, es decir,
 * la derivada del vector de estado y la matriz de transición de estado.
 *
 * @param x Tiempo desde el epoch en segundos.
 * @param yPhi Vector de dimensión (6+36) que comprende el vector de estado (y) y la
 * matriz de transición de estado (Phi) en orden de almacenamiento por columnas.
 * @return Puntero a un array de tamaño 42 que representa la derivada de yPhi.
 */
double* VarEqn(double x,double* yPhi){

double x_pole,y_pole,UT1_UTC,dpsi,LOD,deps,dx_pole,dy_pole,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;

IERS(*global::eopdate,global::Mjd_UTC,'l',x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
double Mjd_UT1 = global::Mjd_TT + (UT1_UTC-TT_UTC)/86400;

//% Transformation matrix
Matrix P = PrecMatrix(MJD_J2000,global::Mjd_TT + x/86400);
Matrix N = NutMatrix(global::Mjd_TT + x/86400);
Matrix T = N * P;
Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

//% State vector components

int s;
    double* r = extractSubarray(yPhi,0,2,s);
    double* v = extractSubarray(yPhi,3,5,s);

    Matrix Phi(6,6);

//% State transition matrix

for(int j=0;j<6;j++){
    for(int i=0;i<6;i++){
        Phi(i+1,j+1) = yPhi[(6*(j+1))+i];
    }

}
    Phi.print();


//% Acceleration and gradient
double* a = AccelHarmonic ( r, E, global::n, global::m );
Matrix G = G_AccelHarmonic ( r, E, global::n, global::m );

//% Time derivative of state transition matrix
auto* yPhip = new double[42];
Matrix dfdy(6,6);

for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
        dfdy(i+1,j+1) = 0.0;                 //% dv/dr(i,j)
        dfdy(i+4,j+1) = G(i+1,j+1);            //% da/dr(i,j)
        if ( i==j ){
            dfdy(i+1,j+3+1) = 1;
        }else{
            dfdy(i+1,j+3+1) = 0;
        }
        //% dv/dv(i,j)
        dfdy(i+3+1,j+3+1) = 0.0;
    }
                //% da/dv(i,j)

}



Matrix Phip = dfdy*Phi;

//% Derivative of combined state vector and state transition matrix
for(int i=0;i<3;i++){
    yPhip[i]   = v[i];                 //% dr/dt(i)
    yPhip[i+3] = a[i];
}
                //% dv/dt(i)

for(int i=0;i<5;i++){
    for(int j=0;j<5;j++){
        yPhip[6*j+i] = Phip(i+1,j+1);     //% dPhi/dt(i,j)
    }

}

    return yPhip;
}