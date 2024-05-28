//
// Created by adboudja on 16/05/2024.
//

#include "Accel.h"
#include "Matrix.h"
#include "global.h"
#include "IERS.h"
#include "timediff.h"
#include "NutMatrix.h"
#include "SAT_Const.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "PrecMatrix.h"
#include "Mjday_TDB.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "auxFunc.h"

/*%--------------------------------------------------------------------------
%
% Accel.m
%
% Purpose:
%   Computes the acceleration of an Earth orbiting satellite due to
%    - the Earth's harmonic gravity field,
%    - the gravitational perturbations of the Sun and Moon
%    - the solar radiation pressure and
%    - the atmospheric drag
%
% Inputs:
%   Mjd_TT      Terrestrial Time (Modified Julian Date)
%   Y           Satellite state vector in the ICRF/EME2000 system
%
% Output:
%   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/

/**
 * @file Accel.h
 * @brief Función para calcular la aceleración de un satélite en órbita terrestre.
 */

/**
 * @brief Calcula la aceleración de un satélite en órbita terrestre debido a varios factores.
 *
 * Esta función computa la aceleración de un satélite en órbita terrestre debido a:
 * - El campo gravitacional armónico de la Tierra,
 * - Las perturbaciones gravitacionales del Sol y la Luna,
 * - La presión de radiación solar,
 * - La resistencia atmosférica.
 *
 * @param x Tiempo en segundos desde la época de referencia.
 * @param Y Vector de estado del satélite en el sistema ICRF/EME2000.
 * @return double* Aceleración (a=d^2r/dt^2) en el sistema ICRF/EME2000.
 *
 * @details
 * La función utiliza varias subrutinas y datos globales para calcular la aceleración:
 * - `IERS`: Para obtener parámetros de rotación de la Tierra.
 * - `timediff`: Para calcular diferencias de tiempo entre varios marcos temporales.
 * - `PrecMatrix`, `NutMatrix`, `PoleMatrix`, `GHAMatrix`: Para calcular matrices de precesión, nutación y otros efectos.
 * - `Mjday_TDB`: Para convertir la fecha a tiempo dinámico baricéntrico.
 * - `JPL_Eph_DE430`: Para obtener posiciones planetarias del efeméride JPL DE430.
 * - `AccelHarmonic`: Para calcular la aceleración debida al campo gravitacional armónico de la Tierra.
 * - `AccelPointMass`: Para calcular las aceleraciones debidas a cuerpos puntuales (Sol, Luna, planetas).
 *
 * @note Esta función asume que los datos globales y funciones auxiliares están correctamente definidas y accesibles.
 *
 * @version 1.0
 * @date 2024-05-16
 *
 * @bug Asegúrese de que las matrices y vectores estén correctamente dimensionados y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
double* Accel(double x, double* Y) {
    double x_pole, y_pole, UT1_UTC, dpsi, LOD, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    double* r_Mercury;
    double *r_Venus;
    double* r_Earth;
    double* r_Mars;
    double* r_Jupiter;
    double* r_Saturn;
    double* r_Uranus;
    double* r_Neptune;
    double* r_Pluto;
    double* r_Moon;
    double* r_Sun;

    IERS(*global::eopdate, global::Mjd_UTC + x/86400, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_UT1 = global::Mjd_UTC + x/86400 + UT1_UTC/86400;
    double Mjd_TT = global::Mjd_UTC + x/86400 + TT_UTC/86400;

    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    JPL_Eph_DE430(MJD_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    // Acceleration due to harmonic gravity field
    int s;
    double* a = AccelHarmonic(extractSubarray(Y, 0, 2, s), E, global::n, global::m);
    double* aux;

    // Luni-solar perturbations
    if (global::sun) {
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Sun, GM_Sun);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
    }

    if (global::moon) {
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Moon);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
    }

    // Planetary perturbations
    if (global::planets) {
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Mars);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Mercury);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Venus);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Jupiter);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Saturn);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Uranus);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Neptune);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
        aux = AccelPointMass(extractSubarray(Y, 0, 2, s), r_Moon, GM_Pluto);
        for (int i = 0; i < 3; i++) {
            a[i] = a[i] + aux[i];
        }
    }
    int f;
    double* dY = concatenateVectors(extractSubarray(Y, 3, 5, f), a, 3, 3);

    return dY;
}