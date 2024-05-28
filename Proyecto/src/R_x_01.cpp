//
// Created by adboudja on 11/04/2024.
//

#include "R_x_01.h"
#include "Matrix.h"
#include "math.h"

/*%--------------------------------------------------------------------------
%  input:
%    angle       - angle of rotation [rad]
%
%  output:
%    rotmat      - vector result
%--------------------------------------------------------------------------*/
/*
 * @file R_x_01.cpp
 * @brief Este archivo contiene la implementación de la función para la rotación alrededor del eje x.
 */
/**
 * @brief Realiza una rotación alrededor del eje x.
 * @param alpha El ángulo de rotación en radianes.
 * @return La matriz de rotación resultante.
 */
Matrix R_x(double alpha){

double C,S;
Matrix rotmat(3,3);

C = cos(alpha);
S = sin(alpha);


rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;


return rotmat;
}