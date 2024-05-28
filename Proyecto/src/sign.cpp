//
// Created by adboudja on 08/05/2024.
//

#include "sign.h"
#include <iostream>
#include <cmath>
/**
 * @file sign.cpp
 * @brief Implementation of the sign function.
 *
 * This file contains the implementation of the sign function, which returns the sign of a number.
 *
 */

/**
 * @brief Returns the sign of a number.
 *
 * This function returns the sign of the given number:
 * - 1.0 if the number is positive.
 * - -1.0 if the number is negative.
 * - 0.0 if the number is zero.
 *
 * @param x The input number.
 * @return The sign of the input number.
 */
double sign(double x) {
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}