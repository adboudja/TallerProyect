//
// Created by Adam on 28/04/2024.
//

#include "MeanObliquity.h"
#include "SAT_Const.h"

/*%--------------------------------------------------------------------------
%
% MeanObliquity.m
%
% Purpose:
%   Computes the mean obliquity of the ecliptic
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   MOblq     Mean obliquity of the ecliptic [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
/**
 * @brief Computes the mean obliquity of the ecliptic.
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 * @return Mean obliquity of the ecliptic [rad].
 */
double MeanObliquity (double Mjd_TT){
double MOblq;

    double T = (Mjd_TT-MJD_J2000)/36525;
    MOblq = Rad *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );
    return MOblq;
}