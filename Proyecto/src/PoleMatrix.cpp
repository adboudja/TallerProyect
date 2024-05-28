//
// Created by adboudja on 08/05/2024.
//

#include "PoleMatrix.h"
#include "Matrix.h"
#include "R_y_01.h"
#include "R_x_01.h"

/*%--------------------------------------------------------------------------
%
% PoleMatrix.m
%
% Purpose:
%   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
%   for a given date
%
% Input:
%   Pole coordinte(xp,yp)
%
% Output:
%   PoleMat   Pole matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file PoleMatrix.h
 * @brief Este archivo contiene la declaración de la función para calcular la matriz de polo.
 */
/**
* @brief Calcula la matriz de polo para una coordenada de polo dada.
* @param xp La coordenada x del polo.
* @param yp La coordenada y del polo.
* @return La matriz de polo.
*/
Matrix PoleMatrix (double xp,double yp){
    Matrix PoleMat = R_y(-xp) * R_x(-yp);
    return PoleMat;
}