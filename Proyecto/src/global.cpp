//
// Created by adboudja on 24/04/2024.
//

#include <cstdio>
#include <cstdlib>
#include "global.h"
/**
 * @file global.h
 * @brief Define una clase global y funciones globales para cargar datos y parámetros globales.
 */

/**
 * @brief Clase global para almacenar datos y parámetros globales.
 *
 * Esta clase proporciona un espacio de nombres para almacenar datos y parámetros globales.
 */ /**
     * @brief Matriz para almacenar los datos de la tabla de parámetros EOP.
     */
Matrix *global::eopdate;
/**
     * @brief Vector para almacenar datos del satélite GEOS3.
     */
double* *global::geos3;
Matrix *global::PC;
/**
    * @brief Matriz para almacenar los coeficientes de la expansión de armónicos del campo gravitatorio terrestre.
    */
Matrix *global::Cnm;
/**
    * @brief Matriz para almacenar los coeficientes de la expansión de armónicos del campo gravitatorio terrestre.
    */
Matrix *global::Snm;

/**
     * @brief Fecha Juliana Modificada (UT1) actual.
     */
double global::Mjd_UTC;
/**
    * @brief Fecha Juliana Modificada (TT) actual.
    */
double global::Mjd_TT;
/**
   * @brief Grado máximo para la expansión de armónicos.
   */
int global::n;
/**
     * @brief Orden máximo para la expansión de armónicos.
     */
int global::m;
/**
    * @brief Indicador para incluir el efecto de la influencia del Sol.
    */
int global::sun;
/**
     * @brief Indicador para incluir el efecto de la influencia de la Luna.
     */
int global::moon;
/**
     * @brief Indicador para incluir el efecto de la influencia de los planetas.
     */
int global::planets;
/**
     * @brief Carga los datos de los parámetros EOP desde un archivo.
     */
void global::eop19620101() {
    global::eopdate = new Matrix(13,21413);

    FILE *fid =fopen("../data/eop19620101.txt","r");

    if(fid== nullptr){
        printf("Error");
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= 21413; i++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &(*global::eopdate)(1, i),
               &(*global::eopdate)(2, i), &(*global::eopdate)(3, i),
               &(*global::eopdate)(4, i), &(*global::eopdate)(5, i), &(*global::eopdate)(6, i),
               &(*global::eopdate)(7, i), &(*global::eopdate)(8, i), &(*global::eopdate)(9, i),
               &(*global::eopdate)(10, i), &(*global::eopdate)(11, i), &(*global::eopdate)(12, i),
               &(*global::eopdate)(13, i));
    }
    fclose(fid);
}
/**
     * @brief Carga los coeficientes de la expansión de armónicos GGM03S desde un archivo.
     */
void global::GGM03S(){
        global::Cnm = new Matrix(181,181);
        global::Snm = new Matrix(181,181);
        Matrix *temp = new Matrix(6,1);
    FILE *fid =fopen("../data/GGM03S.txt","r");

    if(fid== nullptr){
        printf("Error");
        exit(EXIT_FAILURE);
    }
    for(int n=0;n<=180;n++){
        for(int m=0;m<=n;m++){
            fscanf(fid,"%lf %lf %lf %lf %lf %lf",&(*temp)(1,1),&(*temp)(2,1),&(*temp)(3,1),&(*temp)(4,1),&(*temp)(5,1),&(*temp)(6,1));
            (*global::Cnm)(n+1,m+1)=(*temp)(3,1);
            (*global::Snm)(n+1,m+1)=(*temp)(4,1);
        }
    }
    fclose(fid);
}
/**
     * @brief Carga los datos del satélite GEOS3 desde un archivo.
     *
     * @param nobs Número de observaciones.
     */
void global::GEOS3(int nobs){
    FILE *fid = fopen("../data/GEOS3.txt","r");

    if(fid== nullptr){
        printf("Error");
        exit(EXIT_FAILURE);
    }

    double* tline;


    for(int i=0;i<=nobs;i++){
        double value;
        fscanf(fid,"%d",&value);
        (*global::geos3)[i]=value;
    }

    fclose(fid);
}
/**
     * @brief Inicializa los parámetros auxiliares globales.
     */
void global::AuxParam(){
    global::Mjd_UTC=0;

    global::Mjd_TT=0;
    global::m=0;
    global::n=0;
    global::n=0;
    global::n=0;
    global::n=0;
}
/**
    * @brief Carga los coeficientes DE430 desde un archivo.
    */
void global::DE430Coeff(){
    global::PC = new Matrix(2285,1020);
    FILE *fid =fopen("../data/DE430Coeff.txt","r");

    if(fid== nullptr){
        printf("Error");
        exit(EXIT_FAILURE);
    }
    double value;
    for(int i = 1;i<global::PC->getRows();i++){
        for(int j = 1;j<global::PC->getCol();j++){
            fscanf(fid,"%d",&value);
            (*global::PC)(i,j)=value;
        }
    }
}