//
// Created by Adam on 27/04/2024.
//

#include <vector>
#include <cstdio>
#include "Cheb3D.h"



/*%--------------------------------------------------------------------------
%
% Chebyshev approximation of 3-dimensional vectors
%
% Inputs:
%     N       Number of coefficients
%     Ta      Begin interval
%     Tb      End interval
%     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
%     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
%     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file Cheb3D.h
 * @brief Función para evaluar un polinomio de Chebyshev tridimensional.
 */

/**
 * @brief Evalúa un polinomio de Chebyshev tridimensional en un punto dado.
 *
 * Esta función evalúa un polinomio de Chebyshev tridimensional en un punto dado utilizando el algoritmo de Clenshaw.
 *
 * @param t Punto en el que se evaluará el polinomio.
 * @param N Orden del polinomio.
 * @param Ta Extremo inferior del intervalo de definición del polinomio.
 * @param Tb Extremo superior del intervalo de definición del polinomio.
 * @param Cx Coeficientes del polinomio en la dimensión x.
 * @param Cy Coeficientes del polinomio en la dimensión y.
 * @param Cz Coeficientes del polinomio en la dimensión z.
 * @return double* Resultado de evaluar el polinomio en el punto dado.
 *
 * @details
 * La función evalúa un polinomio de Chebyshev tridimensional en el punto especificado utilizando el algoritmo de Clenshaw.
 *
 * @note Esta función asume que los coeficientes del polinomio están correctamente dimensionados y que no hay desbordamientos de memoria.
 *
 * @version 1.0
 * @date Fecha de creación
 *
 * @bug Asegúrese de que los coeficientes del polinomio estén correctamente dimensionados y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
double* Cheb3D(double t,double N,double Ta,double Tb,double* Cx,double* Cy,double* Cz){
    double old_f1[3];
    double* ChebApp=new double[3];
                   //% Check validity

if((t<Ta) || (Tb<t) ){
    printf("ERROR: Time out of range in Cheb3D::Value\n");
}

//% Clenshaw algorithm
double tau = (2*t-Ta-Tb)/(Tb-Ta);

    double f1[3]={0.0, 0.0, 0.0};
    double f2[3]={0.0, 0.0, 0.0};

for (int i=N;i>0;i--){
        // Copy f1 to old_f1
        for (int j = 0; j < 3; j++) {
            old_f1[j] = f1[j];
        }
        f1[0] = 2*tau*f1[0]-f2[0]+Cx[i];
        f1[1] = 2*tau*f1[1]-f2[1]+Cy[i];
        f1[2] = 2*tau*f1[2]-f2[2]+Cz[i];
        for (int j = 0; j < 3; j++) {
            f2[j] = old_f1[j];
        }
}



        ChebApp[0] = tau*f1[0]-f2[0]+Cx[0];
        ChebApp[1] = tau*f1[1]-f2[1]+Cy[0];
        ChebApp[2] = tau*f1[2]-f2[2]+Cz[0];
        return ChebApp;
}