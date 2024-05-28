//
// Created by Adam on 28/04/2024.
//

#include "Geodetic.h"
#include "SAT_Const.h"
#include "norm.h"
#include <cmath>
#include <stdexcept>
#include <ostream>
#include <iostream>

/*%--------------------------------------------------------------------------
%
% Geodetic.m
%
% Purpose:
%   geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
%   from given position vector (r [m])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file Geodetic.cpp
 * @brief Calcula las coordenadas geodésicas (Longitud [rad], latitud [rad], altitud [m]) a partir de un vector de posición dado (r [m]).
 */
/**
* @brief Calcula las coordenadas geodésicas (Longitud [rad], latitud [rad], altitud [m]) a partir de un vector de posición dado (r [m]).
*
* Esta función calcula las coordenadas geodésicas (Longitud, Latitud y Altitud) a partir de un vector de posición en coordenadas cartesianas.
*
* @param r Vector de posición en coordenadas cartesianas [m].
* @param lon Referencia a la variable para almacenar la longitud resultante [rad].
* @param lat Referencia a la variable para almacenar la latitud resultante [rad].
* @param h Referencia a la variable para almacenar la altitud resultante [m].
*
* @throw std::invalid_argument Si el vector de posición es nulo.
*
* @version 1.0
* @date Fecha de creación
*/
void Geodetic(double* r, double& lon, double& lat, double& h){
    double N,Nh,dZ_new;
    double ZdZ;
    double R_equ = R_Earth;
    double  f  = f_Earth;
    double epsRequ = 2.220446049250313e-12*R_equ;        //% Convergence criterion
    double e2      = f*(2.0-f);        //% Square of eccentricity
    double    X = r[0];                  // % Cartesian coordinates
    double Y = r[1];
    double Z = r[2];
    double rho2 = X*X + Y*Y;          // % Square of distance from z-axis

//% Check validity of input data
if (norm(r,3)==0.0){
    std::cerr << "Invalid input in Geodetic constructor" << std::endl;
    lon = 0.0;
    lat = 0.0;
    h   = -R_Earth;
    return;
}

//% Iteration
    double  dZ = e2*Z;

while(true){
     ZdZ    =  Z + dZ;
     Nh     =  sqrt ( rho2 + ZdZ*ZdZ );
    double SinPhi =  ZdZ / Nh;                   // % Sine of geodetic latitude
     N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
    dZ_new =  N*e2*SinPhi;
    if ( fabs(dZ-dZ_new) < epsRequ ){
        break;
    }
    dZ = dZ_new;
}




//% Longitude, latitude, altitude
     lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;


}