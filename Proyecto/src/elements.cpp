//
// Created by adboudja on 08/05/2024.
//

#include <cmath>
#include "elements.h"
#include "norm.h"
#include "SAT_Const.h"

/*%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the osculating Keplerian elements from the satellite state
%   vector for elliptic orbits
%
% Input:
%    y        State vector (x,y,z,vx,vy,vz)
%
% Outputs:
%    p        semilatus rectum [m]
%    a        Semimajor axis
%    e        Eccentricity
%    i        Inclination [rad]
%    Omega    Longitude of the ascending node [rad]
%    omega    Argument of pericenter [rad]
%    M        Mean anomaly [rad]
%
% Notes:
%   The function cannot be used with state vectors describing a circular
%   or non-inclined orbit.
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file elements.cpp
 * @brief Cálculo de los elementos keplerianos osculantes a partir del vector de estado del satélite para órbitas elípticas.
 */

/**
 * @brief Calcula los elementos keplerianos osculantes a partir del vector de estado del satélite para órbitas elípticas.
 *
 * Esta función toma el vector de estado del satélite (posición y velocidad) y calcula los elementos keplerianos osculantes correspondientes
 * para una órbita elíptica.
 *
 * @param y Vector de estado del satélite (x,y,z,vx,vy,vz).
 * @param p Semilatus rectum [m].
 * @param a Semieje mayor.
 * @param e Excentricidad.
 * @param i Inclinación [rad].
 * @param Omega Longitud del nodo ascendente [rad].
 * @param omega Argumento del pericentro [rad].
 * @param M Anomalía media [rad].
 *
 * @note La función no puede utilizarse con vectores de estado que describan una órbita circular o no inclinada.
 *
 * @version 1.0
 * @date 08/05/2024
 *
 * @warning Los resultados pueden no ser válidos para órbitas circulares o no inclinadas.
 */
void elements (double* y,double& p, double& a,double& e, double& i, double& Omega, double& omega,double& M){



    double pi2 = 2*M_PI;
    auto* h = new double[3];
    auto* r=new double[3];;
    auto* v=new double[3];;
    r[0]= y[0];
    r[1]= y[1];
    r[2]= y[2];             // % Position
    v[0]= y[3];
    v[1]= y[4];
    v[2]= y[5];                                       //% Velocity

    cross(r,3,v,3,h);                                   // % Areal velocity
    double magh = norm(h,3);
    p = magh*magh/GM_Earth;
    double H = norm(h,3);

    Omega = atan2 ( h[0], -h[1] );                     //% Long. ascend. node
    Omega = fmod(Omega,pi2);
    i     = atan2 ( sqrt(h[0]*h[0]+h[1]*h[1]), h[2] ); //% Inclination
    double u     = atan2 ( r[2]*H, -r[0]*h[1]+r[0]*h[0] );    //% Arg. of latitude

    double R  = norm(r,3);                                      //% Distance

    a = 1/(2/R-dot(v,3,v,3)/GM_Earth);               //% Semi-major axis

    double eCosE = 1-R/a;                                     //% e*cos(E)
    double eSinE = dot(r,3,v,3)/sqrt(GM_Earth*a);           //% e*sin(E)

    double e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                    // % Eccentricity
    double E  = atan2(eSinE,eCosE);                           //% Eccentric anomaly

    M  = fmod(E-eSinE,pi2);                             //% Mean anomaly

    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          //% True anomaly

    omega = fmod(u-nu,pi2);                             //% Arg. of perihelion
}