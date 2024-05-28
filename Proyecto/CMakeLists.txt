cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(Proyecto
        include/global.h
        include/Matrix.h
        include/R_x_01.h
        include/R_y_01.h
        include/R_z.h
        src/global.cpp
        src/Matrix.cpp
        src/R_x_01.cpp
        src/R_y_01.cpp
        src/R_z.cpp
        EKF_Test.cpp src/sign_.cpp include/sign_.h src/timediff.cpp include/timediff.h src/unit.cpp include/unit.h src/norm.cpp include/norm.h src/AccelPointMass.cpp include/AccelPointMass.h src/AzElPa.cpp include/AzElPa.h src/Cheb3D.cpp include/Cheb3D.h src/EccAnom.cpp include/EccAnom.h src/Frac.cpp include/Frac.h src/SAT_Const.cpp include/SAT_Const.h src/Position.cpp include/Position.h src/NutAngles.cpp include/NutAngles.h src/Mjday_TDB.cpp include/Mjday_TDB.h src/Mjday.cpp include/Mjday.h src/MeanObliquity.cpp include/MeanObliquity.h src/IERS.cpp include/IERS.h src/Geodetic.cpp include/Geodetic.h src/Legendre.cpp include/Legendre.h src/TimeUpdate.cpp include/TimeUpdate.h src/NutMatrix.cpp include/NutMatrix.h src/PoleMatrix.cpp include/PoleMatrix.h src/PrecMatrix.cpp include/PrecMatrix.h src/angl.cpp include/angl.h src/sign.cpp include/sign.h src/elements.cpp include/elements.h src/gmst.cpp include/gmst.h src/gast.cpp include/gast.h src/EqnEquinox.cpp include/EqnEquinox.h src/doubler.cpp include/doubler.h src/LTC.cpp include/LTC.h src/GHAMatrix.cpp include/GHAMatrix.h src/MeasUpdate.cpp include/MeasUpdate.h src/EKF_GEOS3.cpp include/EKF_GEOS3.h src/Accel.cpp include/Accel.h src/JPL_Eph_DE430.cpp include/JPL_Eph_DE430.h
        src/AccelHarmonic.cpp
        include/AccelHarmonic.h
        src/G_AccelHarmonic.cpp
        include/G_AccelHarmonic.h
        src/DEInteg.cpp
        include/DEInteg.h
        src/VarEqn.cpp
        include/VarEqn.h
        src/auxFunc.cpp
        include/auxFunc.h)
