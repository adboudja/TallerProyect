//
// Created by adboudja on 15/05/2024.
//

#include <fstream>
#include <iostream>
#include <vector>
#include "EKF_GEOS3.h"
#include "global.h"
#include "Mjday.h"
#include "SAT_Const.h"
#include "Position.h"
#include "Accel.h"
#include "DEInteg.h"
#include "LTC.h"
#include "IERS.h"
#include "timediff.h"
#include "VarEqn.h"
#include "gmst.h"
#include "R_z.h"
#include "TimeUpdate.h"
#include "AzElPa.h"
#include "MeasUpdate.h"
#include "norm.h"
#include "auxFunc.h"

/*%--------------------------------------------------------------------------
%
%  Initial Orbit Determination using Gauss and Extended Kalman Filter methods
%
% References:
%   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
%   Applications", Springer Verlag, Heidelberg, 2000.
%
%   D. Vallado, "Fundamentals of Astrodynamics and Applications",
%   4th Edition, 2013.
%
%   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
%
% Last modified:   2020/03/16   Meysam Mahooti
%--------------------------------------------------------------------------*/
/**
 * @file main.cpp
 * @brief Determinación de órbita inicial utilizando métodos de Gauss y Filtro de Kalman Extendido.
 */

/**
 * @brief Determinación de órbita inicial utilizando métodos de Gauss y Filtro de Kalman Extendido.
 *
 * Este programa realiza la determinación de órbita inicial de un satélite utilizando métodos de Gauss y el Filtro de Kalman Extendido.
 *
 * @details
 * El programa lee observaciones de posición de un archivo de entrada, realiza la propagación de la órbita utilizando modelos dinámicos y de observación,
 * y aplica el Filtro de Kalman Extendido para estimar la posición y velocidad del satélite.
 *
 * @note Este programa asume que se han proporcionado observaciones válidas y que los modelos dinámicos y de observación son adecuados para el problema.
 *
 * @version 1.0
 * @date Fecha de creación
 *
 * @warning El usuario debe asegurarse de proporcionar observaciones precisas y ajustar los parámetros de los modelos de acuerdo con las necesidades del problema.
 */

void EKF_GEOS3(){
    double x_pole,y_pole,UT1_UTC,dpsi,LOD,deps,dx_pole,dy_pole,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    double Mjd_TT;
    global::DE430Coeff();
    Matrix PC = *global::PC;

//% Model parameters

//% read Earth orientation parameters
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s
%  ----------------------------------------------------------------------------------------------------
global::eop19620101(6);

Matrix eopdata = *global::eopdate;
    double s,az,el,Dist;
    int M,D,h,m;
    double* yPhi;
    int nobs = 46;
Matrix obs(nobs,4);

//% read observations

    std::ifstream file("../data/GEOS3.txt");
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo GEOS3.txt" << std::endl;
        return 1;
    }


    std::string line;
    int i=1;
    while (std::getline(file, line)) {

        if (line.empty()) {
            break;
        }

        int Y = std::stoi(line.substr(0, 4));
        M = std::stoi(line.substr(5, 2));
        D = std::stoi(line.substr(8, 2));
        h = std::stoi(line.substr(12, 2));
        m = std::stoi(line.substr(15, 2));
        s = std::stod(line.substr(18, 6));
        az = std::stod(line.substr(25, 8));
        el = std::stod(line.substr(35, 7));
        Dist = std::stod(line.substr(44, 10));

        obs(i,1) = Mjday(Y,M,D,h,m,s);
        obs(i,2) = Rad*az;
        obs(i,3) = Rad*el;
        obs(i,4) = 1e3*Dist;
        i++;
    }

    file.close();





double sigma_range = 92.5;         // % [m]
double sigma_az = 0.0224*Rad; //% [rad]
double sigma_el = 0.0139*Rad; //% [rad]

//% Kaena Point station
     double   lat = Rad*21.5748;     //% [rad]
    double lon = Rad*(-158.2706);// % [rad]
    double alt = 300.20;              //  % [m]

double* Rss = Position(lon, lat, alt);
Matrix Rs(1,3,Rss,3);

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);


double r2[3]={6221397.62857869,2867713.77965738,3006155.98509949};
double v2[3]={4645.04725161806,-2752.21591588204,-7507.99940987031};
//% [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
//%                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

double* Y0_apr = concatenateVectors(r2,v2,3,3);

double Mjd0 = Mjday(1995,1,29,02,38,0);

double Mjd_UTC = obs(9,1);

    global::Mjd_UTC = Mjd_UTC;
    global::n      = 20;
    global::m      = 20;
    global::sun     = 1;
    global::moon    = 1;
    global::planets = 1;

int n_eqn  = 6;


double* Y = DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

Matrix P(6,6);

    for (int i = 1; i <= 3; ++i) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i <= 6; ++i) {
        P(i, i) = 1e3;
    }
    Matrix LT = LTC(lon,lat);

    yPhi = new double[42];
    Matrix Phi(6,6);

//% Measurement loop
double t = 0;

for(i=1;i<=nobs;i++){
    double t_old = t;
    double* Y_old = Y;

    //% Previous step


    //% Time increment and propagation
    Mjd_UTC = obs(i,1);                       //% Modified Julian Date
    t   = (Mjd_UTC-Mjd0)*86400.0;         //% Time since epoch [s]

    IERS(eopdata,Mjd_UTC,'l',x_pole,y_pole,UT1_UTC,LOD,dpsi,deps, dx_pole,dy_pole, TAI_UTC);
    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
    global::Mjd_UTC = Mjd_UTC;
    global::Mjd_TT = Mjd_TT;

    for (int ii = 0; ii < 6; ++ii) {
        yPhi[ii] = Y_old[ii];
        for (int j = 0; j < 6; ++j) {
            if (ii == j) {
                yPhi[6 * j + ii + 6] = 1; // +6 para compensar los primeros 6 elementos
            } else {
                yPhi[6 * j + ii + 6] = 0; // +6 para compensar los primeros 6 elementos
            }
        }
    }

    yPhi = DEInteg(VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);

    //% Extract state transition matrices
    for (int j = 0; j < 6; ++j) {
        for (i = 0; i < 6; ++i) {
            Phi(i+1,j+1) = yPhi[6 * j + i];
        }
    }

    Y = DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old);

    //% Topocentric coordinates
    double theta = gmst(Mjd_UT1);                    //% Earth rotation
    Matrix U = R_z(theta);
    int o;
    double* rz = extractSubarray(Y, 0, 2, o) ;
    Matrix r(1,3,rz,3);
    Matrix s = LT*(U*r-Rs);                          //% Topocentric position [m]

    //% Time update
    TimeUpdate(P, Phi);

    //% Azimuth and partials
    double Az,El;
    double* dAds;
    double* dEds;
    auto *ss = new double[3];
    ss[0]=s(1,1);
    ss[1]=s(1,2);
    ss[2]=s(1,3);

    AzElPa(ss,dAds, dEds,Az, El);     //% Azimuth, Elevation
    int fu;
    Matrix dAds2(3,3,dAds,fu);
    Matrix dAdY = s.concatenateHorizontal(dAds2*LT*U, Matrix(1,3));

    //% Measurement update [K, Y, P]
    Matrix Azim(1,1);
    Matrix K(6,6);
    Azim(1,1)=Az;
    Matrix Y2(42,1);
    for(int y=0;y<42;y++){
        Y2(y+1,1)=Y[y];
    }
    Matrix z(1,1);
    z(1,1)=obs(i,2);
    Matrix sigaz(1,1);
    sigaz(1,1)=sigma_az;
    MeasUpdate(Y2, z, Azim, sigaz, dAdY, P, 6 ,K);

    //% Elevation and partials
    double* rr = extractSubarray(Y, 0, 2, o) ;
    Matrix r2(1,3,rr,3);
    s = LT*(U*r2-Rs);
    ss[0]=s(1,1);
    ss[1]=s(1,2);
    ss[2]=s(1,3);//% Topocentric position [m]
    AzElPa(ss,dAds, dEds,Az, El);     //% Azimuth, Elevation
    Matrix dEds2(3,3,dEds,fu);
    Matrix dEdY = s.concatenateHorizontal(dEds2*LT*U, Matrix(1,3));
    Matrix Elev(1,1);
    Elev(1,1)=El;
    //% Measurement update


    z(1,1)=obs(i,3);
    Matrix sigel(1,1);
    sigel(1,1)=sigma_el;
    MeasUpdate ( Y2, z, Elev, sigel, dEdY, P, 6,K );

    //% Range and partials
    rr = extractSubarray(Y, 0, 2, o) ;
    Matrix r3(1,3,rr,3);
    s = LT*(U*r3-Rs);
    ss[0]=s(1,1);
    ss[1]=s(1,2);
    ss[2]=s(1,3);//% Topocentric position [m]
    Dist = norm(ss,3);
    auto *dDds = new double[3];
    dDds[0] = (ss[0]/Dist);
    dDds[1] = (ss[1]/Dist);
    dDds[2] = (ss[2]/Dist);//% Range
    Matrix dDds2(3,3,dDds,fu);
    Matrix dDdY = s.concatenateHorizontal(dDds2*LT*U, Matrix(1,3));

    //% Measurement update

    z(1,1)=obs(i,4);
    Matrix sigran(1,1);
    sigran(1,1)=sigma_range;
    Matrix Dist2(1,1);
    Dist2(1,1)=El;
    MeasUpdate ( Y2, z, Dist2,sigran , dDdY, P, 6 ,K);
}

 IERS(eopdata,obs(46,1),'l',x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
 timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
global::Mjd_UTC = Mjd_UTC;
    global::Mjd_TT = Mjd_TT;

double* Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

double Y_true[6] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};

    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n", Y0[0] - Y_true[0]);
    printf("dY%10.1f [m]\n", Y0[1] - Y_true[1]);
    printf("dZ%10.1f [m]\n", Y0[2] - Y_true[2]);

    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n", Y0[3] - Y_true[3]);
    printf("dVy%8.1f [m/s]\n", Y0[4] - Y_true[4]);
    printf("dVz%8.1f [m/s]\n", Y0[5] - Y_true[5]);

}