//
// Created by adboudja on 18/04/2024.
//

#include "R_z.h"
#include "math.h"
/*%--------------------------------------------------------------------------
%  input:
%    angle       - angle of rotation [rad]
%
%  output:
%    rotmat      - vector result
%--------------------------------------------------------------------------*/
/*
 * @file R_z.cpp
 * @brief Este archivo contiene la implementación de la función para la rotación alrededor del eje z.
 */
/**
 * @brief Realiza una rotación alrededor del eje z.
 * @param alpha El ángulo de rotación en radianes.
 * @return La matriz de rotación resultante.
 */
Matrix R_z(double alpha){

double C,S;

C = cos(alpha);
S = sin(alpha);
Matrix rotmat(3,3);

rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;

return rotmat;
}