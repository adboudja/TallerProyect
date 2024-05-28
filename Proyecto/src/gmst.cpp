//
// Created by adboudja on 08/05/2024.
//

#include <cmath>
#include "gmst.h"
#include "Frac.h"

/*%--------------------------------------------------------------------------
%
% Purpose:
%   Greenwich Mean Sidereal Time
%
% Input:
%   Mjd_UT1    Modified Julian Date UT1
%
% Output:
%   gmstime	   GMST in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @file gmst.h
 * @brief Function to calculate Greenwich Mean Sidereal Time.
 */

/**
 * @brief Calculate Greenwich Mean Sidereal Time.
 * @param Mjd_UT1 Modified Julian Date UT1.
 * @return GMST in radians.
 *
 * This function calculates the Greenwich Mean Sidereal Time (GMST) for a given
 * Modified Julian Date (Mjd_UT1) in radians.
 *
 * The algorithm used here is based on a series of terms that account for the
 * rotation of the Earth and temporal evolution.
 *
 * The GMST is normalized to ensure it lies within the range of 0 to 2*pi radians.
 */
double gmst(double Mjd_UT1){



double Secs = 86400.0;                     //  % Seconds per day
        double MJD_J2000 = 51544.5;

double Mjd_0 = floor(Mjd_UT1);
double UT1   = Secs*(Mjd_UT1-Mjd_0);        // % [s]
double T_0   = (Mjd_0  -MJD_J2000)/36525.0;
double T     = (Mjd_UT1-MJD_J2000)/36525.0;

double gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1+ (0.093104-6.2e-6*T)*T*T;    //% [s]

double gmstime = 2*M_PI*Frac(gmst/Secs);      // % [rad], 0..2pi

return gmstime;
}