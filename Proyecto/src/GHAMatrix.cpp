//
// Created by adboudja on 15/05/2024.
//

#include "GHAMatrix.h"
#include "Matrix.h"
#include "gast.h"
#include "R_z.h"

/*%--------------------------------------------------------------------------
%
% GHAMatrix.m
%
% Purpose:
%   Transformation from true equator and equinox to Earth equator and
%   Greenwich meridian system
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
%
% Output:
%   GHAmat    Greenwich Hour Angle matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file GHAMatrix.cpp
 * @brief Calcula la matriz de ángulo horario de Greenwich (GHAmat) para la transformación del verdadero ecuador y equinoccio al ecuador terrestre y al sistema del meridiano de Greenwich.
 */
/**
* @brief Calcula la matriz de ángulo horario de Greenwich (GHAmat) para la transformación del verdadero ecuador y equinoccio al ecuador terrestre y al sistema del meridiano de Greenwich.
*
* Esta función calcula la matriz de ángulo horario de Greenwich (GHAmat) para la transformación del verdadero ecuador y equinoccio al ecuador terrestre y al sistema del meridiano de Greenwich, en función del Tiempo Universal (UT1) modificado.
*
* @param Mjd_UT1 Fecha Juliana Modificada UT1.
*
* @return Matrix Matriz de ángulo horario de Greenwich (GHAmat).
*
* @version 1.0
* @date Fecha de creación
*/
Matrix GHAMatrix (double Mjd_UT1){
    Matrix GHAmat = R_z( gast(Mjd_UT1) );
    return GHAmat;
}


