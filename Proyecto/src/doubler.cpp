//
// Created by adboudja on 09/05/2024.
//

#include "doubler.h"
#include "norm.h"
#include "SAT_Const.h"
#include <cmath>
/**
 * @file doubler.h
 * @brief Función para realizar el proceso de doublers.
 */

/**
 * @brief Realiza el proceso de doublers para determinar la trayectoria de un objeto en un sistema tridimensional.
 *
 * Esta función calcula la trayectoria de un objeto en un sistema tridimensional utilizando el proceso de doublers.
 *
 * @param cc1 Coeficiente de la ecuación del radar en el primer sitio.
 * @param cc2 Coeficiente de la ecuación del radar en el segundo sitio.
 * @param magrsite1 Magnitud del vector de posición del primer sitio.
 * @param magrsite2 Magnitud del vector de posición del segundo sitio.
 * @param magr1in Magnitud del vector de posición inicial.
 * @param magr2in Magnitud del vector de posición final.
 * @param los1 Vector de línea de visión del primer sitio.
 * @param los2 Vector de línea de visión del segundo sitio.
 * @param los3 Vector de línea de visión del tercer sitio.
 * @param rsite1 Vector de posición del primer sitio.
 * @param rsite2 Vector de posición del segundo sitio.
 * @param rsite3 Vector de posición del tercer sitio.
 * @param t1 Tiempo de observación del primer sitio.
 * @param t3 Tiempo de observación del tercer sitio.
 * @param direct Indicador de dirección de la trayectoria.
 * @param r2 Vector de posición del objeto en el segundo sitio.
 * @param r3 Vector de posición del objeto en el tercer sitio.
 * @param f1 Valor de la primera función de corrección.
 * @param f2 Valor de la segunda función de corrección.
 * @param q1 Valor de la distancia de corrección.
 * @param magr1 Magnitud del vector de posición del objeto en el primer sitio.
 * @param magr2 Magnitud del vector de posición del objeto en el segundo sitio.
 * @param a Parámetro semimayor de la órbita.
 * @param deltae32 Ángulo de cambio de la anomalía excéntrica entre los sitios 2 y 3.
 *
 * @details
 * Esta función realiza el proceso de doublers para determinar la trayectoria de un objeto en un sistema tridimensional utilizando la información proporcionada.
 *
 * @note Esta función asume que los vectores de posición y los coeficientes están correctamente definidos y dimensionados.
 *
 * @version 1.0
 * @date 09/05/2024
 *
 * @bug Esta función puede arrojar resultados inesperados si los parámetros de entrada no están correctamente definidos.
 * @warning Asegúrese de proporcionar correctamente los parámetros de entrada para obtener resultados precisos.
 */
void doubler(double cc1,double cc2,double magrsite1,double magrsite2,double magr1in,double magr2in,double* los1,double* los2,double* los3,double* rsite1,double* rsite2,double* rsite3,double t1,double t3,double direct,double* r2,double* r3,double f1,double f2,double q1,double magr1,double magr2,double a,double deltae32){


    double p,esinv2,n,c,s,sinde32,sinde21,cosde32,cosde21,deltam32,deltam12;
    double rho1 = (-cc1 + sqrt(pow(cc1,2)-4*(pow(magrsite1,2)-pow(magr1in,2)))) / 2.0;
    double rho2 = (-cc2 + sqrt(pow(cc2,2)-4*(pow(magrsite2,2)-pow(magr2in,2)))) / 2.0;

    double* r1;
    r1[0] = rho1*los1[0] + rsite1[0];
    r2[0] = rho2*los2[0] + rsite2[0];
    r1[1] = rho1*los1[1] + rsite1[1];
    r2[1] = rho2*los2[1] + rsite2[1];
    r1[2] = rho1*los1[2] + rsite1[2];
    r2[2] = rho2*los2[2] + rsite2[2];

    magr1 = norm(r1,3);
    magr2 = norm(r2,3);
    double* w;
    if (direct == 'y'){
        cross(r1,3,r2,3,w);
        w[0]=w[0]/(magr1*magr2);
        w[1]=w[1]/(magr1*magr2);
        w[2]=w[2]/(magr1*magr2);

    }else{
        cross(r1,3,r2,3,w);
        w[0]=-w[0]/(magr1*magr2);
        w[1]=-w[1]/(magr1*magr2);
        w[2]=-w[2]/(magr1*magr2);

    }

    double rho3 =  -dot(rsite3,3,w,3)/dot(los3,3,w,3);
    r3[0] = rho3*los3[0] + rsite3[0];
    r3[2] = rho3*los3[2] + rsite3[2];
    r3[1] = rho3*los3[1] + rsite3[1];
    double magr3 = norm(r3,3);

    double cosdv21 = dot(r2,3,r1,3)/(magr2*magr1);
    double* sol;
    cross(r2,3,r1,3,sol);
    double sindv21 = norm(sol,3)/(magr2*magr1);
    double dv21 = atan2(sindv21,cosdv21);

    double cosdv31 = dot(r3,3,r1,3)/(magr3*magr1);
    double sindv31 = sqrt(1.0 - pow(cosdv31,2));
    double dv31 = atan2(sindv31,cosdv31);

    double cosdv32 = dot(r3,3,r2,3)/(magr3*magr2);
    double* sol2;
    cross(r3,3,r2,3,sol2);
    double sindv32 = norm(sol2,3)/(magr3*magr2);
    double dv32 = atan2(sindv32,cosdv32);

    if (dv31 > M_PI){
        double c1 = (magr2*sindv32)/(magr1*sindv31);
        double c3 = (magr2*sindv21)/(magr3*sindv31);
        double p = (c1*magr1+c3*magr3-magr2)/(c1+c3-1);
    }else{
        double c1 = (magr1*sindv31)/(magr2*sindv32);
        double c3 = (magr1*sindv21)/(magr3*sindv32);
        p = (c3*magr3-c1*magr2+magr1)/(-c1+c3+1);
    }

    double ecosv1 = p/magr1-1;
    double ecosv2 = p/magr2-1;
    double ecosv3 = p/magr3-1;

    if (dv21 != M_PI){
        esinv2 = (-cosdv21*ecosv2+ecosv1)/sindv21;
    }else{
        esinv2 = (cosdv32*ecosv2-ecosv3)/sindv31;
    }

    double e = sqrt(pow(ecosv2,2)+pow(esinv2,2));
    a = p/(1-pow(e,2));

    if (e*e < 0.99){
        n = sqrt(GM_Earth/pow(a,3));

        s = magr2/p*sqrt(1-pow(e,2))*esinv2;
        c = magr2/p*(pow(e,2)+ecosv2);

         sinde32 = magr3/sqrt(a*p)*sindv32-magr3/p*(1-cosdv32)*s;
         cosde32 = 1-magr2*magr3/(a*p)*(1-cosdv32);
         deltae32 = atan2(sinde32,cosde32);

         sinde21 = magr1/sqrt(a*p)*sindv21+magr1/p*(1-cosdv21)*s;
         cosde21 = 1-magr2*magr1/(a*p)*(1-cosdv21);
         double deltae21 = atan2(sinde21,cosde21);

         deltam32 = deltae32+2*s*pow((sin(deltae32/2)),2)-c*sin(deltae32);
         deltam12 = -deltae21+2*s*pow((sin(deltae21/2)),2)+c*sin(deltae21);
    }else{

        n = sqrt(GM_Earth/pow(-a,3));

        s = magr2/p*sqrt(pow(e,2)-1)*esinv2;
        c = magr2/p*(pow(e,2)+ecosv2);

        double sindh32 = magr3/sqrt(-a*p)*sindv32-magr3/p*(1-cosdv32)*s;
        double sindh21 = magr1/sqrt(-a*p)*sindv21+magr1/p*(1-cosdv21)*s;

        double deltah32 = log( sindh32 + sqrt(pow(sindh32,2+1)) );
        double deltah21 = log( sindh21 + sqrt(pow(sindh21,2+1)) );

        deltam32 = -deltah32+2*s*pow((sinh(deltah32/2)),2)+c*sinh(deltah32);
        deltam12 = deltah21+2*s*pow((sinh(deltah21/2)),2)-c*sinh(deltah21);

        deltae32=deltah32;
    }


    f1 = t1-deltam12/n;
    f2 = t3-deltam32/n;

    q1 = sqrt(pow(f1,2)+pow(f2,2));

}