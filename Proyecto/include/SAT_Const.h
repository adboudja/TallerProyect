//
// Created by Adam on 28/04/2024.
//

#ifndef PROYECTO_SAT_CONST_H
#define PROYECTO_SAT_CONST_H

extern const double pi2;
extern const double Rad;
extern const double Deg;
extern const double Arcs;

//% General
extern const double MJD_J2000;
extern const double T_B1950;
extern const double c_light;
extern const double AU;

//% Physical parameters of the Earth, Sun and Moon

//% Equatorial radius and flattening
extern const double R_Earth;
extern const double f_Earth   ;
extern const double R_Sun    ;
extern const double R_Moon  ;

//% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
extern const double omega_Earth ;

//% Gravitational coefficients
extern const double GM_Earth  ;
extern const double GM_Sun      ;
extern const double GM_Moon     ;
extern const double GM_Mercury  ;
extern const double GM_Venus   ;
extern const double GM_Mars   ;
extern const double GM_Jupiter  ;
extern const double GM_Saturn  ;
extern const double GM_Uranus  ;
extern const double GM_Neptune  ;
extern const double GM_Pluto;

//% Solar radiation pressure at 1 AU
extern const double P_Sol;

#endif //PROYECTO_SAT_CONST_H
