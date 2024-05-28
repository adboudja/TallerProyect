//
// Created by adboudja on 24/04/2024.
//

#include "timediff.h"
#include <vector>
/**
 * @file timediff.cpp
 * @brief Implementation of the timediff function.
 *
 * This file contains the implementation of the timediff function, which calculates various time differences
 * based on input parameters related to time standards.

 */
/**
* @brief Calculates various time differences based on input parameters.
*
* This function calculates the following time differences based on the input parameters:
* - UT1-TAI time difference
* - UTC-GPS time difference
* - UT1-GPS time difference
* - TT-UTC time difference
* - GPS-UTC time difference
*
* The function uses predefined constants for TT-TAI and GPS-TAI time differences.
*
* @param UT1_UTC Difference between UT1 and UTC time standards [s].
* @param TAI_UTC Difference between TAI and UTC time standards [s].
* @param UT1_TAI Output parameter for the UT1-TAI time difference [s].
* @param UTC_GPS Output parameter for the UTC-GPS time difference [s].
* @param UT1_GPS Output parameter for the UT1-GPS time difference [s].
* @param TT_UTC Output parameter for the TT-UTC time difference [s].
* @param GPS_UTC Output parameter for the GPS-UTC time difference [s].
*/
void timediff(double UT1_UTC, double TAI_UTC,double& UT1_TAI,double& UTC_GPS,double& UT1_GPS,double& TT_UTC,double& GPS_UTC){


    double TT_TAI  = +32.184;          //% TT-TAI time difference [s]

    double GPS_TAI = -19.0;            //% GPS-TAI time difference [s]

    double TT_GPS  =  TT_TAI-GPS_TAI;  //% TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI;         //% TAI-GPS time difference [s]

     UT1_TAI = UT1_UTC-TAI_UTC;  //% UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC;         //% UTC-TAI time difference [s]

    UTC_GPS = UTC_TAI-GPS_TAI;  //% UTC_GPS time difference [s]

    UT1_GPS = UT1_TAI-GPS_TAI;  //% UT1-GPS time difference [s]

    TT_UTC  = TT_TAI-UTC_TAI;  // %  TT-UTC time difference [s]

    GPS_UTC = GPS_TAI-UTC_TAI;  //% GPS-UTC time difference [s]




}