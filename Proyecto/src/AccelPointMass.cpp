//
// Created by Adam on 27/04/2024.
//

#include <vector>
#include "AccelPointMass.h"
#include "norm.h"
#include <cmath>
/*%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file AccelPointMass.h
 * @brief Función para calcular la aceleración debida a una masa puntual.
 */

/**
 * @brief Calcula la aceleración perturbacional debido a una masa puntual.
 *
 * Esta función computa la aceleración de un satélite debido a una masa puntual, dados los vectores de posición
 * del satélite y de la masa puntual, así como el coeficiente gravitacional de la masa puntual.
 *
 * @param r Vector de posición del satélite.
 * @param s Vector de posición de la masa puntual.
 * @param GM Coeficiente gravitacional de la masa puntual.
 * @return double* Aceleración (a=d^2r/dt^2).
 *
 * @details
 * La función calcula el vector de aceleración sumando los efectos gravitacionales de la masa puntual sobre el satélite:
 * - `d`: Vector de posición relativa del satélite con respecto a la masa puntual.
 * - `a`: Vector de aceleración debida a la masa puntual.
 *
 * @note Esta función asume que los vectores de posición están correctamente dimensionados y que no hay desbordamientos de memoria.
 *
 * @version 1.0
 * @date 2024-04-27
 *
 * @bug Asegúrese de que los vectores estén correctamente dimensionados y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
double* AccelPointMass(double* r,double* s,float GM){


    auto* a = new double[3];
    auto* d = new double[3];
            //% Relative position vector of satellite w.r.t. point mass
    for (int i = 0; i < sizeof(r); i++) {
        d[i] = r[i] - s[i];
    }


//% Acceleration
    for (int i = 0; i < sizeof(r); i++) {
        a[i] = -GM * ( (d[i]/pow(norm(d,3),3)) + (s[i]/pow(norm(s,3),3)) );
    }
        return a;
}