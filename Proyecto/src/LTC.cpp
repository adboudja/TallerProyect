//
// Created by adboudja on 15/05/2024.
//

#include "LTC.h"
#include "Matrix.h"
#include "R_y_01.h"
#include "R_z.h"

/*%--------------------------------------------------------------------------
%
% LTC.m
%
% Purpose:
%   Transformation from Greenwich meridian system to
%   local tangent coordinates
%
% Inputs:
%   lon      -Geodetic East longitude [rad]
%   lat      -Geodetic latitude [rad]
%
% Output:
%   M        -Rotation matrix from the Earth equator and Greenwich meridian
%             to the local tangent (East-North-Zenith) coordinate system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @brief Transformación del sistema de meridiano de Greenwich a coordenadas tangentes locales.
 *
 * @param lon Longitud este geodésica [rad].
 * @param lat Latitud geodésica [rad].
 * @return Matrix Matriz de rotación desde el ecuador terrestre y el meridiano de Greenwich
 *                al sistema de coordenadas tangentes locales (Este-Norte-Zenit).
 *
 * Esta función calcula una matriz de rotación que transforma las coordenadas de un sistema de meridiano de Greenwich
 * a coordenadas tangentes locales basadas en una longitud y latitud geodésicas dadas.
 */
Matrix LTC(double lon,double  lat){

    Matrix M = R_y(-1.0*lat)*R_z(lon);
    for (int j=1;j<=3;j++){
        double Aux = M(1,j); M(1,j)=M(2,j); M(2,j)=M(3,j); M(3,j)= Aux;
    }
    return M;
}