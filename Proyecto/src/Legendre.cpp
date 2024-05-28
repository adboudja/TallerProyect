//
// Created by Adam on 28/04/2024.
//

#include "Legendre.h"
#include "Matrix.h"
#include <cmath>
//% fi [rad]
/**
 * @brief Calcula los polinomios de Legendre y sus derivadas asociadas.
 *
 * @param n Orden máximo de los polinomios de Legendre.
 * @param m Grado máximo de los polinomios de Legendre.
 * @param fi Ángulo en radianes.
 * @param pnm Matriz de polinomios de Legendre (salida).
 * @param dpnm Matriz de derivadas de los polinomios de Legendre (salida).
 *
 * Calcula los polinomios de Legendre y sus derivadas asociadas hasta el orden y grado dados, evaluados en un ángulo específico.
 */
void Legendre(int n,int m,double fi, Matrix& pnm,Matrix& dpnm){



pnm(1,1)=1;
dpnm(1,1)=0;
pnm(2,2)=sqrt(3)*cos(fi);
dpnm(2,2)=-sqrt(3)*sin(fi);
//% diagonal coefficients
for(int i = 2; i <= n; ++i){
    pnm(i+1,i+1)= sqrt((2.*i+1.)/(2.*i))*cos(fi)*pnm(i,i);
    dpnm(i + 1, i + 1) = sqrt((2. * i + 1.) / (2. * i)) * ((cos(fi) * dpnm(i, i))-(sin(fi) * pnm(i, i)));
}

//% horizontal first step coefficients
    for (int i = 1; i <= n; ++i) {
        pnm(i + 1, i) = sqrt(2. * i + 1.) * sin(fi) * pnm(i, i);
        dpnm(i + 1, i) = sqrt(2. * i + 1.) * ((cos(fi) * pnm(i, i)) + (sin(fi) * dpnm(i, i)));
    }

//% horizontal second step coefficients
    for (int j = 0, k = 2; j <= m; ++j, ++k) {
        int i = k;
        while (i <= n) {
            pnm(i + 1, j + 1) = sqrt((2. * i + 1.) / ((i - j) * (i + j))) *((sqrt(2. * i - 1.) * sin(fi) * pnm(i, j + 1))-(sqrt(((i + j - 1.) * (i - j - 1.)) / (2. * i - 3.)) * pnm(i - 1, j + 1)));
            dpnm(i+1,j+1)=sqrt((2.*i+1.)/((i-j)*(i+j)))*((sqrt(2.*i-1.)*sin(fi)*dpnm(i,j+1))+(sqrt(2.*i-1.)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1.)*(i-j-1.))/(2.*i-3.))*dpnm(i-1,j+1)));
            ++i;
        }


    }
}