//
// Created by adboudja on 15/05/2024.
//

#include "MeasUpdate.h"
#include "Matrix.h"

/**
 * @brief Update step for a measurement in a Kalman filter.
 * @param x State vector to be updated.
 * @param z Measured state vector.
 * @param g Predicted state vector from the system model.
 * @param s Measurement noise.
 * @param G Measurement sensitivity matrix.
 * @param P Covariance matrix of the state estimate.
 * @param n Size of the state vector.
 * @param K Kalman gain (output).
 */
void MeasUpdate(Matrix& x,Matrix z,Matrix g,Matrix s,Matrix G,Matrix& P,int n,Matrix& K){

    int m = z.getCol();

    // Create Inv_W matrix
    Matrix Inv_W(m, m);
    for (int i = 1; i <= m; ++i) {
        Inv_W(i, i) = s(1, i) * s(1, i); // Assuming s is a column vector and s(i, 1) is the ith element
    }

    // Kalman gain
    K = P * G.transpose() * (Inv_W + (G * P * G.transpose())).inverse();

    // State update
    Matrix u = (z - g);
    Matrix v = u*K ;
    x =  x + v ;
    Matrix p = Matrix::identity(n);
    // Covariance update
    P = (p - K * G) * P;
}


