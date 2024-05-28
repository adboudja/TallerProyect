//
// Created by Adam on 27/04/2024.
//

#include "unit.h"
#include "norm.h"
/*%--------------------------------------------------------------------------
%
%  unit.m
%
%  this function calculates a unit vector given the original vector. if a
%  zero vector is input, the vector is set to zero.
%
%  input:
%    vec         - vector
%
%  output:
%    outvec      - unit vector
%
%--------------------------------------------------------------------------
**/
#include <vector>
/**
 * @file unit.cpp
 * @brief Implementación de la función unit.
 *
 * Este archivo contiene la implementación de la función unit, que calcula un vector unitario
 * dado un vector original. Si se proporciona un vector nulo como entrada, el vector resultante
 * se establece en cero.
 *
 * Última modificación: 27/04/2024
 * Autor: Adam
 */
/**
* @brief Calcula un vector unitario dado un vector original.
*
* Esta función calcula un vector unitario dado un vector original. Si se proporciona un vector nulo
* como entrada, el vector resultante se establece en cero.
*
* @param vec Vector original.
* @return Puntero a un array de tamaño 3 que representa el vector unitario calculado.
*/

double * unit(double* vec){



double small = 0.000001;
double magv = norm(vec,3);
double* outvec = new double[3];


if ( magv > small ){
    for (int i=0;i<3;i++){
        outvec[i]= vec[i]/magv;
    }

}else{
    for (int i=0;i<3;i++) {
        outvec[i] = 0.0;
    }
}

    return  outvec;

}