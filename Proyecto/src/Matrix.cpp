#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}
 
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}
 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}
 
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    if (this == &matrix2)
        return *this;



    // Allocate new resources based on matrix2's dimensions
    this->fil = matrix2.fil;
    this->col = matrix2.col;
    this->matrix = new double*[this->fil];
    for (int i = 0; i < this->fil; ++i)
        this->matrix[i] = new double[this->col];

    // Copy the data from matrix2
    for (int i = 0; i < this->fil; ++i)
        for (int j = 0; j < this->col; ++j)
            this->matrix[i][j] = matrix2.matrix[i][j];

    return *this;
}
 
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}

bool Matrix::equalMatrix(const Matrix& m1, const Matrix& m2,double TOL){

    for (int i = 0; i < m1.fil; i++){
        for (int j = 0; j < m1.col; j++){
            if((fabs(m1.matrix[i][j] - m2.matrix[i][j])>TOL)){
                return false;
            }
        }
    }

    return true;

}
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}
 
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
int Matrix::getRows() const
{
    return this->fil;
}
int Matrix::getCol() const
{
    return this->col;
}
int Matrix::find(const Matrix& matrix, double value) {
    int rows = matrix.fil;
    int cols = matrix.col;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (matrix(i+1, j+1) == value) {
                return i * cols + j; // Convertir las coordenadas (i, j) a un índice único
            }
        }
    }
    return -1; // Si el valor no se encuentra en la matriz
}
Matrix Matrix::inverse() const{
    if (fil != col) {
        std::cerr << "La matriz no es cuadrada, no se puede calcular la inversa." << std::endl;
        return Matrix(0, 0); // Devolver una matriz vacía
    }

    // Crear una matriz identidad del mismo tamaño que la matriz original
    Matrix identity(fil, col);
    for (int i = 0; i < this->fil; ++i)
        identity.matrix[i][i] = 1.0;

    // Copiar la matriz original para no modificarla
    Matrix temp(*this);

    // Algoritmo de eliminación gaussiana para obtener la matriz inversa
    for (int i = 0; i < fil; ++i) {
        // Buscar el máximo elemento en la columna actual
        double maxVal = temp.matrix[i][i];
        int maxRow = i;
        for (int j = i + 1; j < fil; ++j) {
            if (std::abs(temp(j, i)) > std::abs(maxVal)) {
                maxVal = temp(j, i);
                maxRow = j;
            }
        }

        // Intercambiar filas para poner el máximo elemento en la diagonal
        if (maxRow != i) {
            temp.swapRows(i, maxRow);
            identity.swapRows(i, maxRow);
        }

        // Dividir la fila por el elemento diagonal para obtener un 1 en la diagonal
        double pivot = temp.matrix[i][i];
        for (int j = 0; j < this->col; ++j) {
            temp.matrix[i][j] /= pivot;
            identity.matrix[i][j] /= pivot;
        }

        // Restar múltiplos de la fila actual para hacer ceros en la columna
        for (int j = 0; j < this->fil; ++j) {
            if (j != i) {
                double factor = temp.matrix[j][i];
                for (int k = 0; k < col; ++k) {
                    temp.matrix[j][k] -= factor * temp.matrix[i][k];
                    identity.matrix[j][k]  -= factor * identity.matrix[i][k] ;
                }
            }
        }
    }

    return identity;
}

Matrix Matrix::identity(int size) {
    Matrix identity(size, size);
    for (int i = 0; i < size; ++i)
        identity(i+1, i+1) = 1.0;
    return identity;
}

void Matrix::swapRows(int row1, int row2) {
    if (row1 == row2)
        return;

    double* temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
}

Matrix Matrix::operator+(double scalar) const {
    Matrix result(fil, col);
    for (int i = 0; i <  this->fil; ++i) {
        for (int j = 0; j <  this->col; ++j) {
            result.matrix[i][j] = this->matrix[i][j] + scalar;
        }
    }
    return result;
}
Matrix Matrix::concatenateHorizontal(Matrix A, Matrix B) {
    if (A.getRows() != B.getRows()) {
        throw std::invalid_argument("Matrices must have the same number of rows for horizontal concatenation.");
    }

    Matrix result(A.getRows(), A.getCol() + B.getCol());

    // Copy A into result
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getCol(); ++j) {
            result(i+1, j+1) = A(i+1, j+1);
        }
    }

    // Copy B into result
    for (int i = 0; i < B.getRows(); ++i) {
        for (int j = 0; j < B.getCol(); ++j) {
            result(i+1, j + A.getCol() + 1) = B(i+1, j+1);
        }
    }

    return result;
}
double* Matrix::multiply(Matrix E, const double* r, int size) {
    auto* result = new double[E.getRows()];
    for (int i = 0; i < E.getRows(); ++i) {
        for (int j = 0; j < E.getCol(); ++j) {
            result[i] += E(i + 1, j + 1) * r[j];
        }
    }
    return result;
}
Matrix Matrix::transpose() const {
    // Create a new matrix with swapped dimensions
    Matrix transposed(col, fil);

    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            // Transpose by swapping rows and columns during copying
            transposed.matrix[j][i] = matrix[i][j];
        }
    }

    return transposed;
}