//
// Created by Adam on 28/04/2024.
//

#include <cmath>
#include "Mjday.h"
/*%--------------------------------------------------------------------------
%  inputs:
%    year        - year
%    mon         - month
%    day         - day
%    hr          - universal time hour
%    min         - universal time min
%    sec         - universal time sec
%
%  output:
%    Mjd         - Modified julian date
%--------------------------------------------------------------------------*/
/**
 * @brief Calculate Modified Julian Date (MJD) from the given date and time.
 * @param yr Year.
 * @param mon Month.
 * @param day Day.
 * @param hr Hour (optional, default value is 0).
 * @param min Minute (optional, default value is 0).
 * @param sec Second (optional, default value is 0).
 * @return Modified Julian Date (MJD).
 */
double Mjday(int yr,int  mon,int day,int hr= 0,int min= 0,int sec= 0){
double Mjd;


        double jd = 367.0 * yr
- floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )
+ floor( 275 * mon / 9.0 )
+ day + 1721013.5
+ ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;

return Mjd = jd-2400000.5;
};