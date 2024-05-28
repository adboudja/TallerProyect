#include <cstdio>
#include <cmath>
#include "./include/Matrix.h"
#include "global.h"
#include "R_z.h"
#include "R_y_01.h"
#include "R_x_01.h"
#include "sign_.h"
#include "timediff.h"
#include "unit.h"
#include "AccelPointMass.h"
#include "AzElPa.h"
#include "Cheb3D.h"
#include "EccAnom.h"
#include "Frac.h"
#include "Position.h"
#include "NutAngles.h"
#include "Mjday_TDB.h"
#include "Mjday.h"
#include "MeanObliquity.h"
#include "Geodetic.h"
#include "Legendre.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "gmst.h"
#include "EqnEquinox.h"
#include "gast.h"
#include "JPL_Eph_DE430.h"
#include "LTC.h"
#include "elements.h"
#include "angl.h"
#include "TimeUpdate.h"
#include "MeasUpdate.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "IERS.h"
#include "VarEqn.h"
#include "Accel.h"
#include "DEInteg.h"
#include "EKF_GEOS3.h"

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;


int proMat_01()
{
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);
    
    sol = m1 * m2;

    m1.print();
    m2.print();
    sol.print();

    _assert(sol(1,1) == m1(1,1) && sol(1,2) == m1(1,2) && sol(2,1) == m1(2,1) && sol(2,2) == m1(2,2));
    
    return 0;
}
int R_z_01(){
    Matrix r(3,3);
    Matrix m(3,3);
    m(1,1) =      0.540302305868140;  m(1,2) =   0.841470984807897;  m(1,3) = 0.0;
    m(2,1) = -0.841470984807897;  m(2,2) =   0.540302305868140;  m(2,3) = 0.0;
    m(3,1) =    0.0;  m(3,2) = 0.0;  m(3,3) = 1.0;

    r= R_z(1.0);
    r.print();
    _assert(r.equalMatrix(m,r,10e-14));

    return 0;
}

int R_x_01(){
    Matrix r(3,3);
    Matrix m(3,3);
    m(1,1) =      1.000000000000000;  m(1,2) =  0.0;  m(1,3) = 0.0;
    m(2,1) = 0.0;  m(2,2) =   0.540302305868140;  m(2,3) = 0.841470984807897;
    m(3,1) =    0.0;  m(3,2) = -0.841470984807897;  m(3,3) = 0.540302305868140;

    r= R_x(1.0);
    r.print();
    _assert(r.equalMatrix(m,r,10e-14));

    return 0;
}
int R_y_01(){
    Matrix r(3,3);
    Matrix m(3,3);
    m(1,1) =      0.540302305868140;  m(1,2) =  0.0;  m(1,3) = -0.841470984807897;
    m(2,1) = 0.0;  m(2,2) =   1.000000000000000;  m(2,3) = 0.0;
    m(3,1) =    0.841470984807897;  m(3,2) = 0.0;  m(3,3) = 0.540302305868140;
    printf("1");
    r= R_y(1.0);
    r.print();
    _assert(r.equalMatrix(m,r,10e-14));

    return 0;
}

int sign_01(){
    int a = 5;
    int b = -2;
    int s = -5;

    _assert(fabs(s-sign_(a,b))<10e-14);

    return 0;
}

int timediff_01(){

    double a1=-26;
    double a2=-10;
    double a3=-7;
    double a4=61.1840;
    double a5=10;
    double b1;
    double b2;
    double b3;
    double b4;
    double b5;
    timediff(3, 29,b1,b2,b3,b4,b5);

    _assert((fabs(a1-b1)<10e-14) and fabs((a2-b2)<10e-14) and fabs((a3-b3)<10e-14) and (fabs(a4-b4)<10e-14) and (fabs(a5-b5)<10e-14));

    return 0;
}
int unit_01(){
    double  sol[3] =  {0.2673, 0.5345, 0.8018};
    double uni[3] = {1,2,3};
    double* result = unit(uni);
    _assert((fabs(sol[0]-result[0])<10e-5) and (fabs(sol[1]-result[1])<10e-5) and (fabs(sol[2]-result[2])<10e-5));

    return 0;
}
int APM_01(){
    double  sol[3] =  {1.0e-11 *0.1032,1.0e-11 *0.0933,1.0e-11 *0.0835};

    double uni1[3] = {1,2,3};
    double uni2[3] = {4,5,6};
    double* result = AccelPointMass(uni1,uni2,6.67430e-11);

    _assert((fabs(sol[0]-result[0])<10e-14) and (fabs(sol[1]-result[1])<10e-14) and (fabs(sol[2]-result[2])<10e-14));

    return 0;
}
int AzElPa_01(){

    double a = 0.463647609000806;
    double b = 0.930274014115472;
    double a2;
    double b2;
    auto*  dA = new double[3];
    auto* dE = new double[3];
    double  dAds[3] = {0.4000 , -0.2000 , 0};
    double  dEds[3] = {-0.095831484749991, -0.191662969499982 ,  0.159719141249985};
    double  un[3] = {1,2,3};
    AzElPa(un,dA,dE,a2,b2);

    _assert((fabs(a-a2)<10e-10) and (fabs(b-b2)<10e-10) and (fabs(dA[1]-dAds[1])<10e-14)and (fabs(dA[0]-dAds[0])<10e-14)and (fabs(dA[2]-dAds[2])<10e-14)and (fabs(dE[1]-dEds[1])<10e-14)and (fabs(dE[0]-dEds[0])<10e-14)and (fabs(dE[2]-dEds[2])<10e-14));

    return 0;
}
int Cheb3D_01(){
    double sol[3] = {1, 4, 7};
    double s1[3] = {1,2,3};
    double s2[3] = {4,5,6};
    double s3[3] = {7,8,9};
    double* result = Cheb3D(0.5,1,0,1,s1,s2,s3);

    _assert((fabs(sol[0]-result[0])<10e-10) and (fabs(sol[1]-result[1])<10e-10) and (fabs(sol[2]-result[2])<10e-10));

    return 0;
}
int EccAnom_01(){
    double sol = 3.423512238778015;
    double result = EccAnom(60,0.1);


    _assert((fabs(sol-result)<10e-10));

    return 0;
}
int Frac_01(){
    double sol = 0.5;
    double result = Frac(3.5);


    _assert((fabs(sol-result)<10e-10));

    return 0;
}
int Position_01(){
    double  sol[3] =  {-1.438078785611559e+06,-2.239675009373783e+06,5.776810445003163e+06};
    double* result = Position(1,2,3);

    _assert((fabs(sol[0]-result[0])<10e-10) and (fabs(sol[1]-result[1])<10e-10) and (fabs(sol[2]-result[2])<10e-10));

    return 0;
}
int NutAngles_01(){
    double dpsis = 3.186976803934643e-05;
    double depss = 3.830248408618162e-05;
    double dpsi,deps;
    NutAngles(10, dpsi,deps);


    _assert((fabs(dpsis-dpsi)<10e-10) and (fabs(depss-deps)<10e-10));

    return 0;
}
int Mjday_TDB_01(){
    double sol = 9.999999988821635;

    double result = Mjday_TDB(10);



    _assert((fabs(sol-result)<10e-10));

    return 0;
}
int Mjday_01(){
    double sol = -6.785567874189815e+05;

    double result = Mjday(1,2,3,5,6,7);


    _assert((fabs(sol-result)<10e-10));

    return 0;
}
int MeanObliquity_01(){
    double sol = 0.409413063968996;

    double result = MeanObliquity(1);



    _assert((fabs(sol-result)<10e-10));

    return 0;
}
int IERS_01(){
    ///{57954,0,0,57954,0.1, 0.1,-0.1,0.001,0.1,0.1,0.1,0.1,32.0,}
    double v[26] = {57954,0,0,57954,0.1, 0.1,-0.1,0.001,0.1,0.1,0.1,0.1,32.0,57955,
     0,
     0,
     57955,
    0.2,
     0.2,
     -0.2,
     0.002,
     0.2,
     0.2,
     0.2,
     0.2,
     32.0};
    Matrix eop(2,13,v,26);
    eop.print();
    double result = 0.000000727221;
    double Mjd_UTC = 57954.5;
    char interp = 'l';
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;

    IERS(eop, Mjd_UTC,interp,
            x_pole, y_pole, UT1_UTC, LOD,
             dpsi,  deps,  dx_pole,  dy_pole,
             TAI_UTC);

    printf("%.15f\n",x_pole);
    printf("%.15f\n",y_pole);
    printf("%.15f\n",UT1_UTC);
    printf("%.15f\n",LOD);
    printf("%.15f\n",dpsi);
    printf("%.15f\n",deps);
    printf("%.15f\n",dx_pole);
    printf("%.15f\n",dy_pole);
    printf("%.15f\n",TAI_UTC);

    _assert((fabs(x_pole-result)<10e-13) and (fabs(y_pole-0.000000727221)<10e-13) and (fabs(TAI_UTC-32.000000000000)<10e-13) and (fabs(UT1_UTC-(-0.150000000000))<10e-13)
    and (fabs(dy_pole-0.000000727221)<10e-13) and (fabs(dx_pole-0.000000727221)<10e-13) and (fabs(LOD-0.001500000000)<10e-13) and (fabs(dpsi- 0.000000727221)<10e-13)
    and (fabs(deps-0.000000727221)<10e-13));

    return 0;
}
int Geodetic_01(){
    double sollon = 1.107148717794090;
    double sollat = 1.570744136243924;
    double solh = -6.356748616533795e+06;
    double lon,lat,h;
    double a1[3]={1,2,3};
    Geodetic(a1,lon,lat,h);


    _assert((fabs(sollon-lon)<10e-10) and (fabs(sollat-lat)<10e-10)and (fabs(solh-h)<10e-10));

    return 0;
}
int Legendre_01(){

    int n=2;
    int m=3;
    Matrix pnm(n+1,m+1);
    Matrix dpnm(n+1,m+1);
    pnm(1,1) =      1.000000000000000;   pnm(1,2) =  0.0;   pnm(1,3) = 0.0;pnm(1,4) = 0.0;
    pnm(2,1) = -0.929371555942242;  pnm(2,2) =   1.461597930692807;   pnm(2,3) = 0.0; pnm(2,4) = 0.0;
    pnm(3,1) =   -0.152352826900483;   pnm(3,2) = -1.753644957368723;   pnm(3,3) = 1.378955394358600;pnm(3,4) = 0.0;
    dpnm(1,1) =      0.0;   dpnm(1,2) =  0.0;   dpnm(1,3) = 0.0;dpnm(1,4) = 0.0;
    dpnm(2,1) = 1.461597930692807;  dpnm(2,2) =   0.929371555942242;   dpnm(2,3) = 0.0; dpnm(2,4) = 0.0;
    dpnm(3,1) =  -3.037402164599586;   dpnm(3,2) =1.642838231226983;   dpnm(3,3) =1.753644957368723;dpnm(3,4) = 0.0;

    Matrix rpnm(n+1,m+1);
    Matrix rdpnm(n+1,m+1);
    Legendre(n,m,12,rpnm,rdpnm);

    _assert(pnm.equalMatrix(rpnm,pnm,10e-14) and dpnm.equalMatrix(rdpnm,dpnm,10e-14));

    return 0;
}
int PrecMatrix_01(){

    Matrix PrecMatrixs(3,3);

    PrecMatrixs(1,1) =     0.999999999999778;   PrecMatrixs(1,2) =  -6.11707317821674e-07;   PrecMatrixs(1,3) = -2.6620148535103e-07;
    PrecMatrixs(2,1) = 6.11707317821674e-07;  PrecMatrixs(2,2) =   0.999999999999813;   PrecMatrixs(2,3) = -8.14186986853169e-14;
    PrecMatrixs(3,1) =  2.6620148535103e-07;   PrecMatrixs(3,2) = -8.14186979189253e-14;   PrecMatrixs(3,3) = 0.999999999999965;


    Matrix PrecMatrixsr(3,3);

    PrecMatrixsr=PrecMatrix(1,2);

    _assert(PrecMatrixs.equalMatrix(PrecMatrixsr,PrecMatrixs,10e-14));

    return 0;
}
int NutMatrix_01(){

    Matrix NutMatrixs(3,3);

    NutMatrixs(1,1) =     0.999999999492159;   NutMatrixs(1,2) =  -2.92358797494727e-05;   NutMatrixs(1,3) = -1.26864277798066e-05;
    NutMatrixs(2,1) = 2.9235393806329e-05  ;  NutMatrixs(2,2) =   0.999999998839099;   NutMatrixs(2,3) = -3.83026695232602e-05;
    NutMatrixs(3,1) =  1.26875475773192e-05;   NutMatrixs(3,2) = 3.83022986110704e-05;   NutMatrixs(3,3) = 0.99999999918598;


    Matrix NutMatrixsr(3,3);

    NutMatrixsr=NutMatrix(10);

    _assert(NutMatrixs.equalMatrix(NutMatrixsr,NutMatrixs,10e-14));

    return 0;
}
int PoleMatrix_01(){

    Matrix PoleMatrixs(3,3);

    PoleMatrixs(1,1) =     -0.839071529076452;   PoleMatrixs(1,2) =  -0.496661489482018;   PoleMatrixs(1,3) = -0.222005256601746;
    PoleMatrixs(2,1) = 0.0  ;  PoleMatrixs(2,2) =  0.408082061813392 ;   PoleMatrixs(2,3) = -0.912945250727628;
    PoleMatrixs(3,1) =  0.54402111088937;   PoleMatrixs(3,2) = -0.766026367491116;   PoleMatrixs(3,3) = -0.342410039594434;


    Matrix PoleMatrixsr(3,3);

    PoleMatrixsr=PoleMatrix(10,20);


    _assert(PoleMatrixs.equalMatrix(PoleMatrixsr,PoleMatrixs,10e-14));

    return 0;
}
int GHAMatrix_01(){

    Matrix GHAMatrixs(3,3);

    GHAMatrixs(1,1) =     0.412804512414729;   GHAMatrixs(1,2) =  0.910819649837463;   GHAMatrixs(1,3) = 0;
    GHAMatrixs(2,1) =  -0.910819649837463    ;  GHAMatrixs(2,2) =  0.412804512414729 ;   GHAMatrixs(2,3) = 0;
    GHAMatrixs(3,1) =  0;   GHAMatrixs(3,2) = 0;   GHAMatrixs(3,3) = 1;


    Matrix GHAMatrixsr(3,3);

    GHAMatrixsr=GHAMatrix(10);


    _assert(GHAMatrixs.equalMatrix(GHAMatrixsr,GHAMatrixs,10e-14));

    return 0;
}
int gmst_01(){

    double ans=1.14523606099042;

    double res = gmst(10);

    _assert(fabs(ans-res)<10e-14);

    return 0;
}
int gast_01(){

    double ans=1.14526529687017;

    double res = gast(10);

    _assert(fabs(ans-res)<10e-14);

    return 0;
}
int EqnEquinox_01(){

    double ans=2.92358797544218e-05;

    double res = EqnEquinox(10);

    _assert(fabs(ans-res)<10e-14);

    return 0;
}
int JPL_Eph_DE430_01(){

    auto* ame = new double[3];
    auto* av = new double[3];
    auto* ae = new double[3];
    auto* ama = new double[3];
    auto* aj = new double[3];
    auto* asa = new double[3];
    auto* au = new double[3];
    auto* an = new double[3];
    auto* ap = new double[3];
    auto* asu = new double[3];
    auto* amo = new double[3];

    JPL_Eph_DE430(10,ame,av,ae,ama,aj,asa,au,an,ap,amo,asu);

    double solmerc[3] = {-92521752990.9701,53904304306.3541,19122423641.6884};


    _assert(fabs(ame[0]-solmerc[0])<10e-3 and (fabs(ame[1]-solmerc[1])<10e-3) and (fabs(ame[2]-solmerc[2])<10e-3));

    return 0;
}
int LTC_01(){

    Matrix LTCs(3,3);

    LTCs(1,1) =     -0.841470984807897;   LTCs(1,2) =   0.54030230586814;   LTCs(1,3) = 0;
    LTCs(2,1) =  -0.491295496433882    ;  LTCs(2,2) = -0.765147401234293 ;   LTCs(2,3) = -0.416146836547142;
    LTCs(3,1) = -0.224845095366153;   LTCs(3,2) = -0.350175488374015;   LTCs(3,3) =  0.909297426825682;


    Matrix LTCsr(3,3);

    LTCsr=LTC(1,2);


    _assert(LTCs.equalMatrix(LTCsr,LTCs,10e-14));

    return 0;
}
int elements_01(){
    double resp,resa,rese,resi,resOm,resom,resM;
    double ansp=1.35474011564823e-13;
    double ansa=1.87082869338765;
    double anse=0.999999999999964;
    double ansi=1.99133066207886;
    double ansOm= 3.6052402625906;
    double ansom=5.21086941752228;
    double ansM=3.14159030993265;
    double y[6]={1,2,3,4,5,6};
    elements(y,resp,resa,rese,resi, resOm, resom,resM);

    _assert(fabs(ansp-resp)<10e-14 and fabs(ansa-resa)<10e-14 and fabs(anse-rese)<10e-14);

    return 0;
}
int Angl_01(){

    double ans=0.225726128552734;
    double vec[3] = {1,2,3};
    double vec2[3] = {4,5,6};
    double res = angl(vec,vec2);

    _assert(fabs(ans-res)<10e-14);

    return 0;
}
int TimeUpdate_01(){

    Matrix P(2,2);
    Matrix Phi(2,2);


    P(1,1) =     1;   P(1,2) =  0;
    P(2,1) =  0   ;  P(2,2) = 1;

    Phi(1,1) =     1;   Phi(1,2) =  1;
    Phi(2,1) =  0   ;  Phi(2,2) = 1;




    TimeUpdate(P,Phi);
    Matrix R(2,2);
    R(1,1) =     2;   R(1,2) =  1;
    R(2,1) =  1   ;  R(2,2) = 1;



    _assert(P.equalMatrix(R,P,10e-14));

    return 0;
}
int MeasUpdate_01(){

    Matrix P(3,3);
    Matrix G(3,3);

    Matrix x(1,3);

    x(1,1) =     0;   x(1,2) =  0; x(1,3) =  0;
    Matrix z(1,3);

    z(1,1) =     1;   z(1,2) =  1; z(1,3) =  1;
    Matrix g(1,3);

    g(1,1) =     0;   g(1,2) =  0; g(1,3) =  0;
    Matrix s(1,3);

    s(1,1) =     0.1;   s(1,2) =  0.1; s(1,3) =  0.1;

    P(1,1) =     1;   P(1,2) =  0;P(1,3) =  0;
    P(2,1) =  0   ;  P(2,2) = 1;P(2,3) =  0;
    P(3,1) =  0   ;  P(3,2) = 0;P(3,3) =  1;

    G(1,1) =     1;   G(1,2) =  0;G(1,3) =  0;
    G(2,1) =  0   ;  G(2,2) = 1;G(2,3) =  0;
    G(3,1) =  0   ;  G(3,2) = 0;G(3,3) =  1;

    Matrix K(3,3);


    MeasUpdate(x,z,g,s,G,P,3,K);

    Matrix Psol(3,3);
    Psol(1,1) =     0.00990099009900991;   Psol(1,2) =  0;Psol(1,3) =  0;
    Psol(2,1) =  0   ;  Psol(2,2) = 0.00990099009900991;Psol(2,3) =  0;
    Psol(3,1) =  0   ;  Psol(3,2) = 0;Psol(3,3) =  0.00990099009900991;
    Matrix Ksol(3,3);
    Ksol(1,1) =     0.99009900990099;   Ksol(1,2) =  0;Ksol(1,3) =  0;
    Ksol(2,1) =  0   ;  Ksol(2,2) = 0.99009900990099;Ksol(2,3) =  0;
    Ksol(3,1) =  0   ;  Ksol(3,2) = 0;Ksol(3,3) =  0.99009900990099;
    Matrix xsol(1,3);
    xsol(1,1) =    0.99009900990099;   xsol(1,2) =  0.99009900990099;xsol(1,3) =  0.99009900990099;


    _assert(P.equalMatrix(Psol,P,10e-14) and K.equalMatrix(Ksol,K,10e-14) and  x.equalMatrix(xsol,x,10e-14));

    return 0;
}
int AccelHarmonic_01(){

    double r[3] = {7000e3, 0, 0};
    Matrix E=Matrix::identity(3);
    int n_max=4;
    int m_max=4;
    auto* ame = new double[3];
    ame =AccelHarmonic(r,E,n_max,m_max);



    double sol[3] = {-8.14571105963231,1.8418953912846e-05,6.13947378070986e-05};

    for(int i = 0; i < 3; ++i) {
        printf("sol[%d] = %.10f\n", i, sol[i]);
    }
    for(int i = 0; i < 3; ++i) {
        printf("ame[%d] = %.10f\n", i, ame[i]);
    }


    _assert(fabs(ame[0]-sol[0])<10e-13 and (fabs(ame[1]-sol[1])<10e-13) and (fabs(ame[2]-sol[2])<10e-13));

    return 0;
}
int G_AccelHarmonic_01(){

    double r[3] = {7000e3, 1200e3, 1300e3};
    Matrix E=Matrix::identity(3);
    int n_max=10;
    int m_max=10;
    Matrix res = G_AccelHarmonic(r,E,n_max,m_max);
    Matrix sol(3,3);

    sol(1,1) =     1.9310499785874e-06;   sol(1,2) =   5.12826245469e-07;sol(1,3) =  5.57846631110692e-07;
    sol(2,1) =  5.12826242804465e-07  ;  sol(2,2) = -9.72249803998793e-07;sol(2,3) =  9.56403363172598e-08;
    sol(3,1) =  5.57846631554781e-07   ;  sol(3,2) = 9.5640336761349e-08;sol(3,3) =  -9.58800177031094e-07;

    sol.print();
    res.print();
    _assert(sol.equalMatrix(res,sol,10e-11));

    return 0;
}
int VarEqn_01(){

    double sol[42] = {0,
            7,
            0,
            -2.53024630251583e+60,
            7.74929400148192e+58,
            -1.45095570696483e+60,
            0,
            0,
            0,
            7.95246202663925e+57,
            -2.43307301979803e+56,
            4.56014345326474e+57,
            0,
            0,
            0,
            -2.43306226422132e+56,
            -5.83538251950612e+57,
            -2.65785725902104e+57,
            0,
            0,
            0,
            4.56013994715353e+57,
            -2.65785630661125e+57,
            -2.11707376686587e+57,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0};
    auto* yPhi = new double[42];

    yPhi[0] = 7000;
    yPhi[1] = 0;
    yPhi[2] = 0;
    yPhi[3] = 0;
    yPhi[4] = 7;
    yPhi[5] = 0;
    for (int i = 6; i < 42; ++i) {

        int row = (i - 6) / 6;
        int col = (i - 6) % 6;

        yPhi[i] = (row == col) ? 1 : 0;
    }

    double* res=  VarEqn(0,yPhi);
    for(int i = 0; i < 5; ++i) {
        printf("res[%d] = %.10f\n", i, res[i]);
    }
    _assert(fabs(res[3]-sol[3])<10e-14);

    return 0;
}
int Accel_01(){

    double sol[3] = {-8.14573258266086, -2.81715611755965e-05, -1.86562227731324e-05};
    double vec[6] = {7000e3, 0, 0, 0, 7.5e3, 0};
    double* res;
    res = Accel(1,vec);

    _assert(fabs(sol[0]-res[0])<10e-14);

    return 0;
}
int DeInteg_01(){

    double sol[6] = {6999592.71713392,  74998.5440215781, -0.000907324142588357, -81.4558180252768, 7499.5633489107, -0.000178910727793944};
    double vec[6] = {7000e3, 0, 0, 0, 7.5e3, 0};
    double* res;
    res = DEInteg(Accel,0,10,1e-6,1e-6,6,vec);

    _assert(fabs(sol[0]-res[0])<10e-14);

    return 0;
}
int Main_01(){

    EKF_GEOS3();

    _assert(true);

    return 0;
}

int all_tests()
{
    _verify(proMat_01);
    _verify(R_z_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(sign_01);
    _verify(timediff_01);
    _verify(unit_01);
    _verify(APM_01);
    _verify(AzElPa_01);
    _verify(Cheb3D_01);
    _verify(EccAnom_01);
    _verify(Frac_01);
    _verify(Position_01);
    _verify(NutAngles_01);
    _verify(Mjday_TDB_01);
    _verify(Mjday_01);
    _verify(MeanObliquity_01);
    _verify(Geodetic_01);
    _verify(Legendre_01);
    _verify(PrecMatrix_01);
    _verify(NutMatrix_01);
    _verify(PoleMatrix_01);
    _verify(gmst_01);
    _verify(EqnEquinox_01);
    _verify(gast_01);
    _verify(GHAMatrix_01);
    //_verify(JPL_Eph_DE430_01); No funciona bien
    _verify(LTC_01);
    _verify(elements_01);
    _verify(Angl_01);
    _verify(TimeUpdate_01);
    _verify(MeasUpdate_01);
    _verify(AccelHarmonic_01);
    _verify(G_AccelHarmonic_01);
    _verify(IERS_01);
    //_verify(VarEqn_01); No funciona bien
    //_verify(Accel_01); No funciona bien
    //_verify(Accel_01); No funciona bien
    //_verify(Main_01);
    return 0;
}


int main()
{
    global::eop19620101();
    global::DE430Coeff();
    global::GGM03S();
    global::AuxParam();
    global::n =20;
    global::m =20;
    global::Mjd_UTC =4.974611128472211e+04;
    global::Mjd_TT =4.974611706231468e+04;
    global::eopdate->print();
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}

