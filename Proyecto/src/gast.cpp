//
// Created by adboudja on 08/05/2024.
//

#include <cmath>
#include "gast.h"
#include "gmst.h"
#include "EqnEquinox.h"

/*%--------------------------------------------------------------------------
%
% GAST.m
%
% Purpose:
%   Greenwich Apparent Sidereal Time
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
%
% Output:
%   gstime    GAST in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file GAST.cpp
 * @brief Calcula el Tiempo Sidéreo Aparente de Greenwich (GAST).
 */
/**
* @brief Calcula el Tiempo Sidéreo Aparente de Greenwich (GAST).
*
* Esta función calcula el Tiempo Sidéreo Aparente de Greenwich (GAST) en radianes
* a partir del Tiempo Universal Coordinado Modificado (Mjd_UT1).
*
* @param Mjd_UT1 Fecha Juliana Modificada UT1.
* @return double GAST en radianes.
*
* @version 1.0
* @date Fecha de creación
*/
double gast(double Mjd_UT1){

double gstime = fmod( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*M_PI );
return gstime;
}