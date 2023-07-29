#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

# include <iomanip>
# include <iostream>
# include <cstdlib>
# include <cmath>

using namespace std;

/// I create the Matrix class for square matrixes.
class Matrix {
/// I state the parameters to initialise the matrix: a list of lists and the order of the matrix.
        double **matrix;
        int m;

    public:
/// Constructor.
        Matrix(int m);
/// Copy constructor.
        Matrix(const Matrix&);
/// Destructor.
        ~Matrix();
/// Assignment operator (=).
        Matrix& operator=(const Matrix&);
/// Addition operator.
        Matrix& operator+=(const Matrix&);
/// Multiplication operator.
        Matrix& operator*=(const Matrix&);
/// Multiplication operator for constants.
        Matrix& operator*=(const double);
/// Find the norm of a matrix.
        double norm() const;
/// Print out a matrix.
        void printMatrix() const;
/// Fill up a matrix with random values.
        void fillMatrix();
/// Create an identity matrix.
        void identityMatrix();
/// Create an array from a matrix.
        double* toArray();
};

#endif // MATRIX_H_INCLUDED
