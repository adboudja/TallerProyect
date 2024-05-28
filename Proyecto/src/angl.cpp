//
// Created by adboudja on 08/05/2024.
//

#include <cmath>
#include "angl.h"
#include "norm.h"
#include "sign_.h"
#include "sign.h"

/*%--------------------------------------------------------------------------
%
%  inputs:
%    vec1         - vector 1
%    vec2         - vector 2
%
%  output:
%    theta        - angle between the two vectors  -pi to pi
%
%--------------------------------------------------------------------------*/
/**
 * @file angl.h
 * @brief Función para calcular el ángulo entre dos vectores.
 */

/**
 * @brief Calcula el ángulo entre dos vectores.
 *
 * Esta función calcula el ángulo entre dos vectores dados, utilizando la fórmula del producto escalar.
 *
 * @param vec1 Primer vector.
 * @param vec2 Segundo vector.
 * @return double Ángulo en radianes.
 *
 * @details
 * La función calcula el ángulo entre dos vectores utilizando la fórmula:
 * \f$ \theta = \arccos \left( \frac{\vec{v}_1 \cdot \vec{v}_2}{\left\| \vec{v}_1 \right\| \left\| \vec{v}_2 \right\|} \right) \f$
 *
 * @note Esta función asume que los vectores están correctamente dimensionados y que no hay desbordamientos de memoria.
 *
 * @version 1.0
 * @date Fecha de creación
 *
 * @bug Asegúrese de que los vectores estén correctamente dimensionados y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
double angl (double* vec1,double* vec2 ){


    double theta;
    double small     = 0.00000001;
double undefined = 999999.1;

    double magv1 = norm(vec1,3);
    double magv2 = norm(vec2,3);

if (magv1*magv2 > small^2){
    double temp= dot(vec1,3,vec2,3) / (magv1*magv2);
    if (fabs( temp ) > 1.0){
        temp= sign(temp) * 1.0;
    }
    theta= acos( temp );
}else{
    theta= undefined;
}

return theta;

}