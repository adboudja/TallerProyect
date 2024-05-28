//
// Created by adboudja on 08/05/2024.
//

#include "TimeUpdate.h"
#include "Matrix.h"
/**
 * @file TimeUpdate.cpp
 * @brief Implementación de la función TimeUpdate.
 *
 * Este archivo contiene la implementación de la función TimeUpdate, que actualiza la matriz de covarianza
 * de error y propagación del estado temporal utilizando la ecuación de actualización del tiempo de Kalman.
 */
/**
* @brief Actualiza la matriz de covarianza de error y propagación del estado temporal.
*
* Esta función actualiza la matriz de covarianza de error y propagación del estado temporal utilizando
* la ecuación de actualización del tiempo de Kalman:
*
* P = Phi * P * Phi^T + Qdt
*
* Donde:
* - P es la matriz de covarianza de error y propagación del estado (entrada y salida).
* - Phi es la matriz de propagación del estado.
* - Qdt es la matriz de covarianza del ruido del proceso multiplicada por el intervalo de tiempo.
*
* @param P Matriz de covarianza de error y propagación del estado [entrada y salida].
* @param Phi Matriz de propagación del estado.
* @param Qdt Matriz de covarianza del ruido del proceso multiplicada por el intervalo de tiempo.
*/
void TimeUpdate(Matrix& P,Matrix Phi, double Qdt) {
    P = Phi * P * Phi.transpose()+ Qdt;
}