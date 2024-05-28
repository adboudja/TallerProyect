//
// Created by adboudja on 08/05/2024.
//

#include "NutMatrix.h"
#include "MeanObliquity.h"
#include "R_x_01.h"
#include "R_z.h"
#include "NutAngles.h"

/*%--------------------------------------------------------------------------
%
% NutMatrix.m
%
% Purpose:
%   Transformation from mean to true equator and equinox
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   NutMat    Nutation matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file NutMatrix.h
 * @brief Este archivo contiene la declaración de la función para calcular la matriz de nutación.
 */
/**
* @brief Calcula la matriz de nutación.
* @param Mjd_TT El Modified Julian Date (Tiempo Terrestre).
* @return La matriz de nutación.
*/
Matrix NutMatrix (double Mjd_TT){

    double dpsi,deps;
        //% Mean obliquity of the ecliptic
double eps = MeanObliquity (Mjd_TT);

//% Nutation in longitude and obliquity
NutAngles (Mjd_TT,dpsi,deps);

//% Transformation from mean to true equator and equinox
Matrix NutMat = R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);

return NutMat;
}