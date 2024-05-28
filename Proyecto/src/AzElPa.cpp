//
// Created by Adam on 27/04/2024.
//

#include <vector>
#include "AzElPa.h"
#include "norm.h"

#include <cmath>
/*%--------------------------------------------------------------------------
%
% Purpose:
%  Computes azimuth, elevation and partials from local tangent coordinates
%
% Input:
%   s      Topocentric local tangent coordinates (East-North-Zenith frame)
%
% Outputs:
%   A      Azimuth [rad]
%   E      Elevation [rad]
%   dAds   Partials of azimuth w.r.t. s
%   dEds   Partials of elevation w.r.t. s
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/

/**
 * @file AzElPa.h
 * @brief Función para calcular la azimut y la elevación a partir de un vector de posición.
 */

/**
 * @brief Calcula el azimut y la elevación a partir de un vector de posición.
 *
 * Esta función calcula el azimut y la elevación a partir de un vector de posición dado, así como las derivadas
 * parciales de los ángulos respecto a las componentes del vector de posición.
 *
 * @param s Vector de posición.
 * @param[out] dAds Derivadas parciales del azimut respecto a las componentes del vector de posición.
 * @param[out] dEds Derivadas parciales de la elevación respecto a las componentes del vector de posición.
 * @param[out] Azim Azimut resultante en radianes.
 * @param[out] Elev Elevación resultante en radianes.
 *
 * @details
 * La función calcula el azimut y la elevación a partir del vector de posición utilizando las siguientes fórmulas:
 * - Azimut: \f$ \text{Az} = \text{atan2}(s[0], s[1]) \f$
 * - Elevación: \f$ \text{El} = \text{atan} \left( \frac{s[2]}{\sqrt{s[0]^2 + s[1]^2}} \right) \f$
 *
 * Además, calcula las derivadas parciales del azimut y la elevación respecto a las componentes del vector de posición:
 * - \f$ \frac{\partial \text{Az}}{\partial s_i} = \frac{s_1}{\text{rho}^2} \f$
 * - \f$ \frac{\partial \text{El}}{\partial s_i} = \frac{-s_0 \cdot s_2 / \text{rho}}{\text{dot}(s, s)} \f$
 *
 * @note Esta función asume que el vector de posición está correctamente dimensionado y que no hay desbordamientos de memoria.
 *
 * @version 1.0
 * @date Fecha de creación
 *
 * @bug Asegúrese de que el vector de posición esté correctamente dimensionado y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
void  AzElPa(double* s, double*& dAds, double*& dEds,double& Azim, double& Elev){


double pi2 = 2.0*M_PI;

double rho = sqrt(s[0]*s[0]+s[1]*s[1]);

//% Angles
        double Az = atan2(s[0],s[1]);

if (Az<0.0){
    Az = Az+pi2;
}

    double El = atan ( s[2] / rho );

//% Partials
     dAds[0] = s[1]/(rho*rho);
     dAds[1] =-s[0]/(rho*rho);
     dAds[2] = 0.0 ;
     dEds[0] = (-s[0]*s[2]/rho)/ dot(s,3,s,3);
     dEds[1] = (-s[1]*s[2]/rho)/ dot(s,3,s,3);
     dEds[2] = rho / dot(s,3,s,3) ;

     Azim = Az;
     Elev = El;

}