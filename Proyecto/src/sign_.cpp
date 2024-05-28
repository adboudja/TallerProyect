//
// Created by adboudja on 24/04/2024.
//

#include <cmath>
#include "sign_.h"
/**
 * @file sign_.cpp
 * @brief Implementation of the sign_ function.
 *
 * This file contains the implementation of the sign_ function, which returns the signed magnitude of a number.
 */
/**
* @brief Returns the signed magnitude of a number.
*
* This function returns the signed magnitude of the first parameter, using the sign of the second parameter:
* - If the second parameter (b) is greater than 0, the function returns the absolute value of the first parameter (a).
* - If the second parameter (b) is less than or equal to 0, the function returns the negative of the absolute value of the first parameter (a).
*
* @param a The number whose signed magnitude is to be determined.
* @param b The sign parameter. If greater than 0, the sign of a will be preserved; otherwise, the sign will be inverted.
* @return The signed magnitude of the input number.
*/
double sign_(double a, double b){
    if(b > 0.0){
        return fabs(a);
    }else{
        return -1 * fabs(a);
    }
}