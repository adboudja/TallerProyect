//
// Created by Adam on 27/04/2024.
//

#include "EccAnom.h"
#include <cmath>
#include <stdexcept>
#include <limits>

/*%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the eccentric anomaly for elliptic orbits
%
% Inputs:
%   M         Mean anomaly in [rad]
%   e         Eccentricity of the orbit [0,1]
%
% Output:
%             Eccentric anomaly in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @brief Calcula la anomalía excéntrica para una órbita elíptica.
 *
 * Esta función calcula la anomalía excéntrica para una órbita elíptica utilizando el método de iteración de Kepler.
 *
 * @param M Anomalía media.
 * @param e Excentricidad de la órbita.
 * @return double Anomalía excéntrica.
 *
 * @details
 * La anomalía excéntrica es un parámetro angular que define la posición de un objeto en una órbita elíptica en función del tiempo.
 * Este método utiliza iteraciones de Kepler para encontrar la anomalía excéntrica correspondiente a la anomalía media y la excentricidad dadas.
 *
 * @note Esta función asume que la anomalía media y la excentricidad están correctamente definidas y en el rango adecuado.
 *
 * @version 1.0
 * @date 27/04/2024
 *
 * @throws std::runtime_error Si hay problemas de convergencia durante las iteraciones.
 * @warning Asegúrese de proporcionar valores válidos para la anomalía media y la excentricidad para obtener resultados precisos.
 */