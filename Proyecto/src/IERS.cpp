//
// Created by Adam on 28/04/2024.
//

#include <cmath>
#include "IERS.h"
#include "Matrix.h"
#include "SAT_Const.h"

/**
 * @file IERS.h
 * @brief Función para interpolar parámetros de rotación terrestre a partir de datos IERS.
 */

/**
 * @brief Interpola parámetros de rotación terrestre a partir de datos IERS.
 * @param eop Matriz que contiene los parámetros de orientación terrestre (datos IERS).
 * @param Mjd_UTC Fecha Juliana Modificada (UTC).
 * @param interp Método de interpolación ('l' para lineal, 'n' para vecino más cercano).
 * @param x_pole Coordenada del polo en radianes (salida).
 * @param y_pole Coordenada del polo en radianes (salida).
 * @param UT1_UTC Diferencia de tiempo UT1 - UTC en segundos (salida).
 * @param LOD Longitud del día en segundos (salida).
 * @param dpsi Desplazamiento del polo celestial en segundos de arco (salida).
 * @param deps Desplazamiento del polo celestial en segundos de arco (salida).
 * @param dx_pole Tasa de cambio de la coordenada del polo en segundos de arco por día (salida).
 * @param dy_pole Tasa de cambio de la coordenada del polo en segundos de arco por día (salida).
 * @param TAI_UTC Diferencia de tiempo TAI - UTC en segundos (salida).
 *
 * Esta función interpola los parámetros de rotación terrestre (x_pole, y_pole, UT1_UTC, LOD,
 * dpsi, deps, dx_pole, dy_pole, TAI_UTC) a partir de los datos IERS proporcionados en la matriz `eop`.
 * El método de interpolación puede especificarse como lineal ('l') o vecino más cercano ('n').
 * Los valores interpolados se devuelven en los parámetros de salida.
 */

void IERS(Matrix eop, double Mjd_UTC, char interp,
          double& x_pole, double& y_pole, double& UT1_UTC, double& LOD,
          double& dpsi, double& deps, double& dx_pole, double& dy_pole,
          double& TAI_UTC) {

    Matrix preeop(1,13 ); // Matriz para almacenar preeop
    Matrix nexteop(1,13 ); // Matriz para almacenar nexteop
    int i;
    if (interp == 'l'){


    //% linear interpolation
    double mjd = (floor(Mjd_UTC));
        for (i = 0; i < eop.getCol(); i++) {
            if (eop(4, i + 1) == mjd) {
                mjd = i + 1;
                break;

            }
        }

        for (int row = 1; row <= 13; row++) {
            preeop(1, row) = eop(row, mjd);
        }
        for (int row = 1; row <= 13; row++) {
            nexteop(1, row) = eop(row, mjd+1);
        }

        preeop.print();
    double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
    double fixf = mfme / 1440;
    /*% Setting
    of
    IERS
    Earth
    rotation
    parameters
    % (UT1 - UTC[s], TAI - UTC[s], x["], y ["])*/
    preeop.print();
        x_pole = preeop(1, 5) + (nexteop(1, 5) - preeop(1, 5)) * fixf;
        y_pole = preeop(1, 6) + (nexteop(1, 6) - preeop(1, 6)) * fixf;
        UT1_UTC = preeop(1, 7) + (nexteop(1, 7) - preeop(1, 7)) * fixf;
        LOD = preeop(1, 8) + (nexteop(1, 8) - preeop(1, 8)) * fixf;
        dpsi = preeop(1, 9) + (nexteop(1, 9) - preeop(1, 9)) * fixf;
        deps = preeop(1, 10) + (nexteop(1, 10) - preeop(1, 10)) * fixf;
        dx_pole = preeop(1, 11) + (nexteop(1, 11) - preeop(1, 11)) * fixf;
        dy_pole = preeop(1, 12) + (nexteop(1, 12) - preeop(1, 12)) * fixf;
        TAI_UTC = preeop(1, 13);
    x_pole = x_pole /Arcs;  //% Pole coordinate[rad]
    y_pole = y_pole /Arcs;  //% Pole coordinate[rad]
    dpsi = dpsi /Arcs;
    deps = deps /Arcs;
    dx_pole = dx_pole /Arcs; //% Pole coordinate[rad]
    dy_pole = dy_pole /Arcs; //% Pole coordinate[rad]
    }else{
        if(interp == 'n'){
            double mjd = (floor(Mjd_UTC));
            double i = eop.find(eop,mjd);
            for (int col = 1; col <= eop.getCol(); col++) {
                eop(i, col) = eop(i, col);
            }


            /*% Setting
            of
            IERS
            Earth

            rotation
            parameters
            % (UT1 - UTC[s], TAI - UTC[s], x["], y ["])*/
            x_pole = eop(1, 5) / Arcs;  // Pole coordinate [rad]
            y_pole = eop(1, 6) / Arcs;  // Pole coordinate [rad]
            UT1_UTC = eop(1, 7);        // UT1 - UTC time difference [s]
            LOD = eop(1, 8);            // Length of day [s]
            dpsi = eop(1, 9) / Arcs;
            deps = eop(1, 10) / Arcs;
            dx_pole = eop(1, 11) / Arcs; // Pole coordinate [rad]
            dy_pole = eop(1, 12) / Arcs; // Pole coordinate [rad]
            TAI_UTC = eop(1, 13);           //% TAI - UTC time difference[s]
        }
    }

}