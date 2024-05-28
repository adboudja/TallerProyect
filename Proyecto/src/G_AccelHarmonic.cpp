//
// Created by Adam on 24/05/2024.
//

#include "G_AccelHarmonic.h"
#include "Matrix.h"
#include "AccelHarmonic.h"

/*%--------------------------------------------------------------------------
%
% G_AccelHarmonic.m
%
% Purpose:
%   Computes the gradient of the Earth's harmonic gravity field
%
% Inputs:
%   r           Satellite position vector in the true-of-date system
%   U           Transformation matrix to body-fixed system
%   n           Gravity model degree
%   m 			Gravity model order
%
% Output:
%   G    		Gradient (G=da/dr) in the true-of-date system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file G_AccelHarmonic.cpp
 * @brief Calcula la matriz de gradiente de la aceleración armónica en coordenadas cartesianas.
 */
/**
* @brief Calcula la matriz de gradiente de la aceleración armónica en coordenadas cartesianas.
*
* Esta función calcula la matriz de gradiente de la aceleración armónica en coordenadas cartesianas
* para una posición dada, utilizando la diferencia finita para aproximar las derivadas parciales.
*
* @param r Puntero al vector de posición.
* @param U Matriz de coeficientes armónicos.
* @param n_max Grado máximo de los términos armónicos.
* @param m_max Orden máximo de los términos armónicos.
* @return Matrix Matriz de gradiente de la aceleración armónica.
*
* @version 1.0
* @date Fecha de creación
*/
Matrix G_AccelHarmonic( double* r,Matrix U,int n_max,int m_max ){

    double d = 1.0;  // % Position increment [m]
    // Inicializar matriz G de ceros
    Matrix G(3, 3);

    // Inicializar vector columna dr de ceros
    double dr[3] = {0.0, 0.0, 0.0};
    double* aux = r;
//% Gradient
    for (int i = 0; i < 3; ++i) {
//% Set offset in i-th component of the position vector
        for (int j = 0; j < 3; ++j) {
            dr[j] = 0.0;
        }
        dr[i] = d;
//% Acceleration difference
        for (int j = 0; j < 3; ++j) {
            r[j] = r[j] + dr[j];
        }
        double* da1 = AccelHarmonic(r, U, n_max, m_max);
        for (int j = 0; j < 3; ++j) {
            dr[j] = -dr[j]; // Invertir el vector de desplazamiento para calcular la otra diferencia de aceleración
        }
        for (int j = 0; j < 3; ++j) {
            aux[j] = aux[j] + dr[j];
        }
        double* da2 = AccelHarmonic(aux, U, n_max, m_max);
//% Derivative with respect to i-th axis
        for (int j = 0; j < 3; ++j) {
            G(j + 1, i + 1) = (da1[j] - da2[j]) / d;
        }
    }
    return G;
}
