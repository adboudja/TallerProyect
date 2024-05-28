//
// Created by Adam on 24/05/2024.
//
#include "Matrix.h"
#include "norm.h"
#include "Legendre.h"
#include "global.h"
#include <cmath>
/*%--------------------------------------------------------------------------
%
% AccelHarmonic.m
%
% Purpose:
%   Computes the acceleration due to the harmonic gravity field of the
%   central body
%
% Inputs:
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file AccelHarmonic.h
 * @brief Función para calcular la aceleración debido al campo gravitacional armónico de un cuerpo central.
 */

/**
 * @brief Calcula la aceleración debida al campo gravitacional armónico de un cuerpo central.
 *
 * Esta función computa la aceleración de un satélite debido al campo gravitacional armónico de un cuerpo central,
 * considerando un máximo grado y orden específico.
 *
 * @param r Vector de posición del satélite en el sistema inercial.
 * @param E Matriz de transformación al sistema centrado en el cuerpo central.
 * @param n_max Máximo grado del campo armónico.
 * @param m_max Máximo orden del campo armónico (m_max <= n_max; m_max = 0 solo para armónicos zonales).
 * @return double* Aceleración en el sistema inercial.
 *
 * @details
 * La función utiliza las siguientes subrutinas y constantes globales:
 * - `norm`: Para calcular la norma del vector.
 * - `Legendre`: Para calcular los polinomios de Legendre asociados.
 *
 * @note Esta función asume que las matrices y vectores están correctamente dimensionados y que no hay desbordamientos de memoria.
 *
 * @version 1.0
 * @date 2024-05-24
 *
 * @bug Asegúrese de que las matrices y vectores estén correctamente dimensionados y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
double* AccelHarmonic(double* r,Matrix E,int n_max,int m_max){



double r_ref = 6378.1363e3;   //% Earth's radius [m]; GGM03S
double gm    = 398600.4415e9; //% [m^3/s^2]; GGM03S

//% Body-fixed position
double* r_bf =E.multiply(E,r,3);

//% Auxiliary quantities
double d = norm(r_bf,3);                    // % distance
double latgc = asin(r_bf[2]/d);
double lon = atan2(r_bf[1],r_bf[0]);


    Matrix pnm(n_max+1,m_max+1);
    Matrix dpnm(n_max+1,m_max+1);
    Legendre(n_max, m_max, latgc, pnm, dpnm);

double dUdr = 0;
double dUdlatgc = 0;
double dUdlon = 0;
double q3 = 0;double q2 = q3;double q1 = q2;


    double b1,b2,b3;
for(int n=0;n<=n_max;n++){
    b1 = (-gm/pow(d,2))*pow((r_ref/d),n)*(n+1);
    b2 =  (gm/d)*pow((r_ref/d),n);
    b3 =  (gm/d)*pow((r_ref/d),n);
    for (int m=0;m<=m_max;m++){
        q1 = q1 + pnm(n+1,m+1)*((*global::Cnm)(n+1,m+1)*cos(m*lon)+(*global::Snm)(n+1,m+1)*sin(m*lon));
        q2 = q2 + dpnm(n+1,m+1)*((*global::Cnm)(n+1,m+1)*cos(m*lon)+(*global::Snm)(n+1,m+1)*sin(m*lon));
        q3 = q3 + m*pnm(n+1,m+1)*((*global::Snm)(n+1,m+1)*cos(m*lon)-(*global::Cnm)(n+1,m+1)*sin(m*lon));
    }
    dUdr     = dUdr     + q1*b1;
    dUdlatgc = dUdlatgc + q2*b2;
    dUdlon   = dUdlon   + q3*b3;
    q3 = 0; q2 = q3; q1 = q2;
}







//% Body-fixed acceleration
double r2xy = pow(r_bf[0],2)+pow(r_bf[1],2);

double ax = (1/d*dUdr-r_bf[2]/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf[0]-(1/r2xy*dUdlon)*r_bf[1];
double ay = (1/d*dUdr-r_bf[2]/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf[1]+(1/r2xy*dUdlon)*r_bf[0];
double az =  1/d*dUdr*r_bf[2]+sqrt(r2xy)/pow(d,2)*dUdlatgc;

    auto* a_bf = new double[3];
    a_bf[0] = ax;
    a_bf[1] = ay;
    a_bf[2] = az;

//% Inertial acceleration
 double* a = E.multiply(E,a_bf,3);

 return a;

}
