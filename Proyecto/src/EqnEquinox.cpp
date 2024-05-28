//
// Created by adboudja on 08/05/2024.
//

#include "EqnEquinox.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include <cmath>
/*%--------------------------------------------------------------------------
%
% EqnEquinox.m
%
% Purpose:
%   Computation of the equation of the equinoxes
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%    EqE      Equation of the equinoxes
%
% Notes:
%   The equation of the equinoxes dpsi*cos(eps) is the right ascension of
%   the mean equinox referred to the true equator and equinox and is equal
%   to the difference between apparent and mean sidereal time.
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file EqnEquinox.cpp
 * @brief Cálculo de la ecuación de los equinoccios.
 */

/**
 * @brief Calcula la ecuación de los equinoccios.
 *
 * Esta función calcula la ecuación de los equinoccios, que es la ascensión recta del equinoccio medio referida al verdadero ecuador y equinoccio,
 * y es igual a la diferencia entre el tiempo sidéreo aparente y el tiempo sidéreo medio.
 *
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre).
 * @return double Ecuación de los equinoccios.
 *
 * @note La ecuación de los equinoccios dpsi*cos(eps) es la ascensión recta del equinoccio medio referida al verdadero ecuador y equinoccio, y es igual
 * a la diferencia entre el tiempo sidéreo aparente y el tiempo sidéreo medio.
 *
 * @version 1.0
 * @date Fecha de creación
 */
double EqnEquinox (double Mjd_TT){
    //% Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles (Mjd_TT,dpsi,deps);

    //% Equation of the equinoxes
    double EqE = dpsi * cos ( MeanObliquity(Mjd_TT) );

    return EqE;
}

