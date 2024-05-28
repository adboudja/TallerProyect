//
// Created by adboudja on 16/05/2024.
//

#include <cstring>
#include "JPL_Eph_DE430.h"
#include "global.h"
#include "Cheb3D.h"

/*%--------------------------------------------------------------------------
%
% JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
%   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Notes: Light-time is already taken into account
%
% Last modified:   2018/01/11   M. Mahooti
%
%--------------------------------------------------------------------------*/

double* extractSubarray(double* source, int startIndex, int endIndex, int& newSize) {
    newSize = endIndex - startIndex;
    auto* subarray = new double[newSize];
    for (int i = 0; i < newSize; ++i) {
        subarray[i] = source[startIndex + i];
    }
    return subarray;
}
void appendArrays(double*& dest, int destSize, const double* src, int srcSize) {
    // Calculate the new size
    int newSize = destSize + srcSize;

    // Allocate new memory
    auto* newArray = new double[newSize];

    // Copy the contents of the old array
    if (destSize > 0) {
        std::memcpy(newArray, dest, destSize * sizeof(double));
    }

    // Copy the contents of the source array
    std::memcpy(newArray + destSize, src, srcSize * sizeof(double));

    // Deallocate the old array
    delete[] dest;

    // Update the pointer and size
    dest = newArray;
    destSize = newSize;
}
/**
 * @brief Calcula la posición ecuatorial del sol, la luna y los nueve planetas principales utilizando las efemérides de JPL.
 *
 * @param Mjd_TDB Fecha Juliana Modificada de TDB.
 * @param r_Mercury Puntero a la posición de Mercurio (salida).
 * @param r_Venus Puntero a la posición de Venus (salida).
 * @param r_Earth Puntero a la posición de la Tierra (salida).
 * @param r_Mars Puntero a la posición de Marte (salida).
 * @param r_Jupiter Puntero a la posición de Júpiter (salida).
 * @param r_Saturn Puntero a la posición de Saturno (salida).
 * @param r_Uranus Puntero a la posición de Urano (salida).
 * @param r_Neptune Puntero a la posición de Neptuno (salida).
 * @param r_Pluto Puntero a la posición de Plutón (salida).
 * @param r_Moon Puntero a la posición de la Luna (salida).
 * @param r_Sun Puntero a la posición del Sol (salida).
 *
 * Calcula la posición ecuatorial del sol, la luna y los nueve planetas principales utilizando las efemérides de JPL para una fecha dada.
 */
void JPL_Eph_DE430(double Mjd_TDB,double*& r_Mercury,double*& r_Venus,double*& r_Earth,double*& r_Mars,double*& r_Jupiter,double*& r_Saturn,double*& r_Uranus, double*&r_Neptune,double*& r_Pluto,double*& r_Moon,double*& r_Sun){


    double Mjd0;
    auto* PCtemp = new double[global::PC->getCol()];
    double JD = Mjd_TDB + 2400000.5;
    int i;
    for (i = 0; i < 2285; i++) {
        if((*global::PC)(i+1,1)<=JD && JD<=(*global::PC)(i+1,2)){
            break;

        }
    }

    for(int in=0;in<global::PC->getCol();in++){

        PCtemp[in] = (*global::PC)(i,in);
    }


    double t1 = PCtemp[1]-2400000.5; //% MJD at start of interval

    double dt = Mjd_TDB - t1;
    auto* temp = new int[4];
    int j=0;
    for (i = 231; i <= 270; i += 13) {
        temp[j]=i;
        j++;
    }
    int k;
    double* Cx_Earth = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
    double* Cy_Earth = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
    double* Cz_Earth = extractSubarray(PCtemp,temp[2],temp[3]-1,k);
    for (i=0; i <= 3; i ++) {
        temp[i]+=39;

    }

    double* Cx = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
    double* Cy = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
    double* Cz = extractSubarray(PCtemp,temp[2],temp[3]-1,k);
    appendArrays(Cx_Earth, sizeof(Cx_Earth), Cx, sizeof(Cx));
    appendArrays(Cy_Earth, sizeof(Cy_Earth), Cx, sizeof(Cy));
    appendArrays(Cz_Earth, sizeof(Cz_Earth), Cz, sizeof(Cz));
    if (0<=dt && dt<=16){
        j=0;
         Mjd0 = t1;
    }else{
        if(16<dt && dt<=32){
            j=1;
            Mjd0 = t1+16*j;
        }
    }
    double* aux = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, extractSubarray(Cx_Earth,13*j+1,13*j+13,k),extractSubarray(Cy_Earth,13*j+1,13*j+13,k), extractSubarray(Cz_Earth,13*j+1,13*j+13,k));
    for(i=0;i<3;i++){
        r_Earth[i] = 1e3*aux[i];
    }


    temp = new int[4];
    j=0;
    for (i = 441; i <= 480; i += 13) {
        temp[j]=i;
        j++;
    }
    double* Cx_Moon = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
    double* Cy_Moon = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
    double* Cz_Moon = extractSubarray(PCtemp,temp[2],temp[3]-1,k);
    for(int i=1;i<=7;i++){
        for (j=0; j <= 3; j ++) {
            temp[j]+=39;

        }
        Cx = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
        Cy = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
        Cz = extractSubarray(PCtemp,temp[2],temp[3]-1,k);
        appendArrays(Cx_Moon, sizeof(Cx_Moon), Cx, sizeof(Cx));
        appendArrays(Cy_Moon, sizeof(Cy_Moon), Cx, sizeof(Cy));
        appendArrays(Cz_Moon, sizeof(Cz_Moon), Cz, sizeof(Cz));
    }
    if (0<=dt && dt<=4){
        j=0;
        Mjd0 = t1;
    }
    else{
        if(4<dt && dt<=8){
            j=1;
            Mjd0 = t1+4*j;
        }else{
            if(8<dt && dt<=12){
                j=2;
                Mjd0 = t1+4*j;
            }else{
                if(12<dt && dt<=16){
                    j=3;
                    Mjd0 = t1+4*j;
                }else{
                    if(16<dt && dt<=20){
                        j=4;
                        Mjd0 = t1+4*j;
                    }else{
                        if(20<dt && dt<=24){
                            j=5;
                            Mjd0 = t1+4*j;
                        }else{
                            if(24<dt && dt<=28){
                                j=6;
                                Mjd0 = t1+4*j;
                            }else{
                                if(28<dt && dt<=32){
                                    j=7;
                                    Mjd0 = t1+4*j;
                                }
                            }
                        }
                    }
                }
            }
        }
    }








    aux = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, extractSubarray(Cx_Moon,13*j+1,13*j+13,k),extractSubarray(Cy_Moon,13*j+1,13*j+13,k), extractSubarray(Cz_Moon,13*j+1,13*j+13,k));
    for(i=0;i<3;i++){
        r_Moon[i] = 1e3*aux[i];
    }
    temp = new int[4];
    j=0;
    for (i = 753; i <= 786; i += 11) {
        temp[j]=i;
        j++;
    }

    double* Cx_Sun = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
    double* Cy_Sun = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
    double* Cz_Sun = extractSubarray(PCtemp,temp[2],temp[3]-1,k);
    for (i=0; i <= 3; i ++) {
        temp[i]+=33;

    }
    Cx = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
    Cy = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
    Cz = extractSubarray(PCtemp,temp[2],temp[3]-1,k);
    appendArrays(Cx_Sun, sizeof(Cx_Sun), Cx, sizeof(Cx));
    appendArrays(Cy_Sun, sizeof(Cy_Sun), Cx, sizeof(Cy) );
    appendArrays(Cz_Sun, sizeof(Cz_Sun), Cz, sizeof(Cz) );

    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }else{
        if(16<dt && dt<=32){
            j=1;
            Mjd0 = t1+16*j;
        }
    }



    aux = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, extractSubarray(Cx_Sun,11*j+1,11*j+11,k),extractSubarray(Cy_Sun,11*j+1,11*j+11,k), extractSubarray(Cz_Sun,11*j+1,11*j+11,k));
    for(i=0;i<3;i++){
        r_Sun[i] = 1e3*aux[i];
    }

    temp = new int[4];
    j=0;
    for (i = 3; i <= 45; i += 14) {
        temp[j]=i;
        j++;
    }

    double* Cx_Mercury = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
    double* Cy_Mercury = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
    double* Cz_Mercury = extractSubarray(PCtemp,temp[2],temp[3]-1,k);

    for (i=1;i<=3;i++){
        for (i=0; i <= 3; i ++) {
            temp[i]+=42;

        }
        Cx = extractSubarray(PCtemp,temp[0],temp[1]-1,k);
        Cy = extractSubarray(PCtemp,temp[1],temp[2]-1,k);
        Cz = extractSubarray(PCtemp,temp[2],temp[3]-1,k);
        appendArrays(Cx_Mercury, sizeof(Cx_Mercury) , Cx, sizeof(Cx) );
        appendArrays(Cy_Mercury, sizeof(Cy_Mercury), Cx, sizeof(Cy) );
        appendArrays(Cz_Mercury, sizeof(Cz_Mercury), Cz, sizeof(Cz));
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }else{
        if(8<dt && dt<=16){
            j=1;
            Mjd0 = t1+8*j;
        }else{
            if(16<dt && dt<=24){
                j=2;
                Mjd0 = t1+8*j;
            }else{
                if(24<dt && dt<=32){
                    j=3;
                    Mjd0 = t1+8*j;
                }
            }
        }
    }







    aux = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, extractSubarray(Cx_Mercury,14*j+1,14*j+14,k),extractSubarray(Cy_Mercury,14*j+1,14*j+14,k), extractSubarray(Cz_Mercury,14*j+1,14*j+14,k));
    for( i=0;i<3;i++){
        r_Mercury[i] = 1e3*aux[i];
    }

    temp = new int[4];
    j = 0;
    for (i = 171; i <= 201; i += 10) {
        temp[j] = i;
        j++;
    }

    double* Cx_Venus = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Venus = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Venus = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

    for (i=0; i <= 3; i ++) {
        temp[i]+=30;

    }
     Cx = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
     Cy = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
     Cz = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

    appendArrays(Cx_Venus, sizeof(Cx_Venus) , Cx, sizeof(Cx) );
    appendArrays(Cy_Venus, sizeof(Cy_Venus), Cy, sizeof(Cy) );
    appendArrays(Cz_Venus, sizeof(Cz_Venus), Cz, sizeof(Cz));

// Calculate Mjd0 and compute r_Venus
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }


    aux = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, extractSubarray(Cx_Venus,10*j+1,10*j+10,k),extractSubarray(Cy_Venus,10*j+1,10*j+10,k), extractSubarray(Cz_Venus,10*j+1,10*j+10,k));
    for(i=0;i<3;i++){
        r_Venus[i] = 1e3*aux[i];
    }


    temp = new int[4];
    j = 0;
    for (i = 309; i <= 342; i += 11) {
        temp[j] = i;
        j++;
    }

    double* Cx_Mars = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Mars = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Mars = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

// Calculate Mjd0 and compute r_Mars
    j = 0;
    Mjd0 = t1;

    aux = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, extractSubarray(Cx_Mars,11*j+1,11*j+11,k),extractSubarray(Cy_Mars,11*j+1,11*j+11,k), extractSubarray(Cz_Mars,11*j+1,11*j+11,k));
    for(i=0;i<3;i++){
        r_Mars[i] = 1e3*aux[i];
    }


    temp = new int[4];
    j = 0;
    for (i = 342; i <= 366; i += 8) {
        temp[j] = i;
        j++;
    }

    double* Cx_Jupiter = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Jupiter = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Jupiter = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

// Calculate Mjd0 and compute r_Jupiter
    j = 0;
    Mjd0 = t1;

    aux = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, extractSubarray(Cx_Jupiter,8*j+1,8*j+8,k),extractSubarray(Cy_Jupiter,8*j+1,8*j+8,k), extractSubarray(Cz_Jupiter,8*j+1,8*j+8,k));
    for(i=0;i<3;i++){
        r_Jupiter[i] = 1e3*aux[i];
    }

    temp = new int[4];
    j = 0;
    for (i = 366; i <= 387; i += 7) {
        temp[j] = i;
        j++;
    }

    double* Cx_Saturn = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Saturn = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Saturn = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

    j = 0;
    Mjd0 = t1;


    aux = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, extractSubarray(Cx_Saturn,7*j+1,7*j+7,k),extractSubarray(Cy_Saturn,7*j+1,7*j+7,k), extractSubarray(Cz_Saturn,7*j+1,7*j+7,k));
    for(i=0;i<3;i++){
        r_Saturn[i] = 1e3*aux[i];
    }



    temp = new int[4];
    for (i = 387; i <= 405; i += 6) {
        temp[j] = i;
        j++;
    }

    double* Cx_Uranus = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Uranus = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Uranus = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

    j = 0;
    Mjd0 = t1;


    aux = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, extractSubarray(Cx_Uranus,6*j+1,6*j+6,k),extractSubarray(Cy_Uranus,6*j+1,6*j+6,k), extractSubarray(Cz_Uranus,6*j+1,6*j+6,k));
    for(i=0;i<3;i++){
        r_Uranus[i] = 1e3*aux[i];
    }


    temp = new int[4];
    for (i = 405; i <= 423; i += 6) {
        temp[j] = i;
        j++;
    }

    double* Cx_Neptune = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Neptune = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Neptune = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

    j = 0;
    Mjd0 = t1;

    aux = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, extractSubarray(Cx_Neptune,6*j+1,6*j+6,k),extractSubarray(Cy_Neptune,6*j+1,6*j+6,k), extractSubarray(Cz_Neptune,6*j+1,6*j+6,k));
    for(i=0;i<3;i++){
        r_Neptune[i] = 1e3*aux[i];
    }

    temp = new int[4];
    for (i = 423; i <= 441; i += 6) {
        temp[j] = i;
        j++;
    }

    double* Cx_Pluto = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Pluto = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Pluto = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

    j = 0;
    Mjd0 = t1;

    aux = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, extractSubarray(Cx_Pluto,6*j+1,6*j+6,k),extractSubarray(Cy_Pluto,6*j+1,6*j+6,k), extractSubarray(Cz_Pluto,6*j+1,6*j+6,k));
    for(i=0;i<3;i++){
        r_Pluto[i] = 1e3*aux[i];
    }



    temp = new int[4];
    for (i = 819; i <= 839; i += 10) {
        temp[j] = i;
        j++;
    }

    double* Cx_Nutations = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Nutations = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);

    for (i = 0; i < 3; i++) {
        for (j=0; j <= 3; j ++) {
            temp[j]+=20;

        }
        Cx = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
        Cy = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
        appendArrays(Cx_Nutations, sizeof(Cx_Nutations) , Cx, sizeof(Cx));
        appendArrays(Cy_Nutations, sizeof(Cy_Nutations) , Cy, sizeof(Cy));
    }

    if (0 <= dt && dt <= 8)
        j = 0;
    else if (8 < dt && dt <= 16)
        j = 1;
    else if (16 < dt && dt <= 24)
        j = 2;
    else if (24 < dt && dt <= 32)
        j = 3;


    double* Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, extractSubarray(Cx_Nutations,10*j+1,10*j+10,k),extractSubarray(Cy_Nutations,10*j+1,10*j+10,k), new double[10]);




    temp = new int[4];
    j=0;
    for (i = 899; i <= 929; i += 10) {
        temp[j] = i;
        j++;
    }

    double* Cx_Librations = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
    double* Cy_Librations = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
    double* Cz_Librations = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);

    for (i = 0; i < 3; i++) {
        for (j=0; j <= 3; j ++) {
            temp[j]+=20;

        }
        double* Cx = extractSubarray(PCtemp, temp[0], temp[1] - 1, k);
        double* Cy = extractSubarray(PCtemp, temp[1], temp[2] - 1, k);
        double* Cz = extractSubarray(PCtemp, temp[2], temp[3] - 1, k);
        appendArrays(Cx_Librations, sizeof(Cx_Librations), Cx, sizeof(Cx));
        appendArrays(Cy_Librations, sizeof(Cy_Librations), Cy, sizeof(Cy));
        appendArrays(Cz_Librations, sizeof(Cz_Librations), Cz, sizeof(Cz));
    }
    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    double* Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, extractSubarray(Cx_Librations,10*j+1,10*j+10,k),extractSubarray(Cy_Librations,10*j+1,10*j+10,k), extractSubarray(Cx_Librations,10*j+1,10*j+10,k));
    double EMRAT = 81.30056907419062; //% DE430
    double  EMRAT1 = 1/(1+EMRAT);

    for(i=0;i<3;i++){
        r_Earth[i] = r_Earth[i]-EMRAT1*r_Moon[i];
        r_Mercury[i] = -r_Earth[i]+r_Mercury[i];
        r_Venus[i] = -r_Earth[i]+r_Venus[i];
        r_Mars[i] = -r_Earth[i]+r_Mars[i];
        r_Jupiter[i] = -r_Earth[i]+r_Jupiter[i];
        r_Saturn[i] = -r_Earth[i]+r_Saturn[i];
        r_Uranus[i] = -r_Earth[i]+r_Uranus[i];
        r_Neptune[i] = -r_Earth[i]+r_Neptune[i];
        r_Pluto[i] = -r_Earth[i]+r_Pluto[i];
        r_Sun[i] = -r_Earth[i]+r_Sun[i];
    }


}


