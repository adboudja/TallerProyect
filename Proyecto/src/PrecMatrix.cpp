//
// Created by adboudja on 08/05/2024.
//

#include "PrecMatrix.h"
#include "SAT_Const.h"
#include "Matrix.h"
#include "R_z.h"
#include "R_y_01.h"

/*%--------------------------------------------------------------------------
%
% PrecMatrix.m
%
% Purpose:
%
%   Precession transformation of equatorial coordinates
%
% Input:
%   Mjd_1     Epoch given (Modified Julian Date TT)
%   MjD_2     Epoch to precess to (Modified Julian Date TT)
%
% Output:
%   PrecMat   Precession transformation matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/*
 * @file PrecMatrix.cpp
 * @brief Este archivo contiene la implementación de la función para calcular la matriz de precesión.
 */
/**
 * @brief Calcula la matriz de transformación de precesión de coordenadas ecuatoriales.
 * @param Mjd_1 Fecha de época dada (Modified Julian Date TT).
 * @param Mjd_2 Fecha de época a la que precesar (Modified Julian Date TT).
 * @return La matriz de transformación de precesión.
 */
Matrix PrecMatrix (double Mjd_1,double  Mjd_2){


    Matrix PrecMat(3,3);

    double T  = (Mjd_1-MJD_J2000)/36525;
    double dT = (Mjd_2-Mjd_1)/36525;

//% Precession angles
    double zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
    double z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs;
    double theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;

//% Precession matrix
    PrecMat = R_z(-z) * R_y(theta) * R_z(-zeta);

    return PrecMat;
}