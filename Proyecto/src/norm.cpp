//
// Created by Adam on 27/04/2024.
//

#include "norm.h"
#include <cmath>
#include <vector>
/**
 * @file norm.h
 * @brief Este archivo contiene las declaraciones de las funciones para el cálculo de la norma y productos vectoriales.
 */

/**
 * @brief Calcula la norma euclidiana de un vector.
 * @param v El vector.
 * @param n El tamaño del vector.
 * @return La norma euclidiana del vector.
 */
double norm(double v[], int n){

    double sol = 0;

    for(int i=0;i<n;i++){
        sol += v[i]*v[i];
    }
    return sqrt(sol);
}
/**
 * @brief Calcula el producto punto de dos vectores.
 * @param v1 El primer vector.
 * @param n1 El tamaño del primer vector.
 * @param v2 El segundo vector.
 * @param n2 El tamaño del segundo vector.
 * @return El producto punto de los dos vectores.
 * @throw Se lanza una excepción si las dimensiones de los vectores no son iguales.
 */
double dot(double v1[],int n1, double v2[], int n2){

    double sol = 0;
    if(n1!=n2){
        throw "Not equal dimensions";
    }else{
        for(int i=0;i<n1;i++){
            sol += v1[i]*v2[i];
        }
    }

    return sol;
}
/**
 * @brief Calcula el producto cruz de dos vectores tridimensionales.
 * @param v1 El primer vector.
 * @param n1 El tamaño del primer vector.
 * @param v2 El segundo vector.
 * @param n2 El tamaño del segundo vector.
 * @param sol2 El vector resultante del producto cruz.
 * @throw Se lanza una excepción si las dimensiones de los vectores no son iguales.
 */
void cross(double v1[],int n1, double v2[], int n2, double sol2[]){


    if(n1!=n2){
        throw "Not equal dimensions";
    }else{

        sol2[0] = v1[1]*v2[2]-v1[2]*v2[1];
        sol2[1] = v1[2]*v2[0]-v1[0]*v2[2];
        sol2[2] = v1[0]*v2[1]-v1[1]*v2[0];

    }


}