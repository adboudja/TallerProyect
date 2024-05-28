//
// Created by Adam on 28/04/2024.
//

#include <vector>
#include "Position.h"
#include "SAT_Const.h"
#include <cmath>
/*%--------------------------------------------------------------------------
%
% Position.m
%
% Purpose:
%   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
%   latitude [rad], altitude [m])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file Position.cpp
 * @brief Este archivo contiene la definición de la función para calcular el vector de posición.
 */
/**
 * @brief Calcula el vector de posición desde coordenadas geodésicas.
 * @param lon Longitud en radianes.
 * @param lat Latitud en radianes.
 * @param h Altitud en metros.
 * @return Un puntero a un array de tamaño 3 que representa el vector de posición.
 */
double* Position(double lon,double lat,double h){
    double  R_equ,f,e2,CosLat,SinLat,N;
    auto* r = new double[3];
    R_equ = R_Earth;
    f     = f_Earth;

    e2     = f*(2.0-f);   //% Square of eccentricity
    CosLat = cos(lat);    //% (Co)sine of geodetic latitude
    SinLat = sin(lat);

    //% Position vector
    N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

    r[0] =  (N+h)*CosLat*cos(lon);
    r[1] =  (N+h)*CosLat*sin(lon);
    r[2] =  ((1.0-e2)*N+h)*SinLat;

    return r;
}



