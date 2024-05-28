//
// Created by Adam on 27/04/2024.
//

#include <cmath>
#include "Frac.h"
/*%--------------------------------------------------------------------------
%
%  Fractional part of a number (y=x-[x])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file Frac.cpp
 * @brief Parte fraccionaria de un número (y = x - [x]).
 */

/**
 * @brief Calcula la parte fraccionaria de un número.
 *
 * Esta función calcula la parte fraccionaria de un número, es decir, la diferencia entre el número y su parte entera.
 *
 * @param x El número del cual se calculará la parte fraccionaria.
 * @return double Parte fraccionaria del número.
 *
 * @version 1.0
 * @date 27/04/2024
 */
double Frac(double x){
    double res;
    res = x-floor(x);
    return res;
}


