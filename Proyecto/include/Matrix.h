#ifndef _MATRIX_
#define _MATRIX_
/**
 * @class Matrix
 * @brief Represents a mathematical matrix.
 */
class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);
        Matrix(const Matrix& m);
        ~Matrix();
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator-(const Matrix& matrix2);
        Matrix  operator*(const Matrix& matrix2);
        Matrix operator+(double scalar) const;
        Matrix concatenateHorizontal(Matrix A, Matrix B);
    double& operator()(const int i, const int j) const;
        bool equalMatrix(const Matrix& m1,const Matrix& m2, double TOL);
        int find(const Matrix& matrix, double value);
        double* multiply(Matrix E, const double* r, int size);
        void print();
        int getRows() const;
        Matrix inverse() const;
        Matrix transpose() const;

        static Matrix identity(int size);

        void swapRows(int row1, int row2);
        int getCol() const;
 
    private:
        void initMatrix();
 
    private:
        int fil;
        int col;
        double **matrix;



};

#endif
