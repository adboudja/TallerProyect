//
// Created by Adam on 25/05/2024.
//

#include "auxFunc.h"
/**
 * @file auxFunc.h
 * @brief Funciones auxiliares para manipulación de vectores.
 */

/**
 * @brief Extrae un subarreglo de un arreglo dado.
 *
 * Esta función extrae un subarreglo de un arreglo dado, comenzando desde el índice de inicio hasta el índice de finalización.
 *
 * @param source Arreglo de entrada del que se extraerá el subarreglo.
 * @param startIndex Índice de inicio del subarreglo (inclusive).
 * @param endIndex Índice de finalización del subarreglo (exclusivo).
 * @param[out] newSize Tamaño del subarreglo extraído.
 * @return double* Subarreglo extraído.
 *
 * @details
 * La función extrae un subarreglo del arreglo de entrada, comenzando desde el índice de inicio (inclusive) hasta el índice de finalización (exclusivo).
 *
 * @note Esta función asume que el arreglo de entrada está correctamente dimensionado y que no hay desbordamientos de memoria.
 *
 * @version 1.0
 * @date 2024-05-25
 *
 * @bug Asegúrese de que el arreglo de entrada esté correctamente dimensionado y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
double* extractSubarray(const double* source, int startIndex, int endIndex, int& newSize) {
    newSize = endIndex - startIndex;
    auto* subarray = new double[newSize+1];
    for (int i = 0; i < newSize+1; ++i) {
        subarray[i] = source[startIndex + i];
    }
    return subarray;
}
/**
 * @brief Concatena dos vectores en uno nuevo.
 *
 * Esta función concatena dos vectores en uno nuevo, manteniendo el orden de los elementos.
 *
 * @param v1 Primer vector.
 * @param v2 Segundo vector.
 * @param size1 Tamaño del primer vector.
 * @param size2 Tamaño del segundo vector.
 * @return double* Vector resultante de la concatenación.
 *
 * @details
 * La función crea un nuevo vector que contiene todos los elementos del primer vector seguidos por todos los elementos del segundo vector.
 *
 * @note Esta función asume que los vectores están correctamente dimensionados y que no hay desbordamientos de memoria.
 *
 * @version 1.0
 * @date 2024-05-25
 *
 * @bug Asegúrese de que los vectores estén correctamente dimensionados y que no haya desbordamientos de memoria.
 * @warning Verificar la precisión de los parámetros de entrada para obtener resultados precisos.
 */
double* concatenateVectors(const double* v1, const double* v2, int size1, int size2) {
    auto* result = new double[size1 + size2];
    for (int i = 0; i < size1; ++i) {
        result[i] = v1[i];
    }
    for (int i = 0; i < size2; ++i) {
        result[size1 + i] = v2[i];
    }
    return result;
}