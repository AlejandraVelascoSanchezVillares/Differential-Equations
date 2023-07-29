/// I include the files I will need to execute my program.

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "Matrix.h"

/// I create the constructor with parameters which allocates memory for an mxm matrix.
Matrix::Matrix(int n) : m(n)
{
    m = n;
    if(m >= 0)
        matrix = (double **)calloc(m, sizeof(double **));
        for (int i = 0; i < m; i++) { matrix[i] = (double *) calloc(m, sizeof(double));}
}

/// I create the copy constructor.
Matrix::Matrix (const Matrix& M): m(M.m)
{
    matrix = new double*[m];
    for(int i = 0; i < m; ++i){matrix[i]= new double[m];}
    for(int i=0; i<m; i++) {
        for(int j=0; j<m; j++) {
            matrix[i][j] = M.matrix[i][j];
        }
    }
}

/// I create the destructor.
Matrix::~Matrix()
{
    for (int i = 0; i < m; i++) { free(matrix[i]);}
    free(matrix);
}

/// Assignment operator (=)
Matrix& Matrix::operator=(const Matrix& M)
{
    if (this != &M) {
        for (int i = 0; i < m; ++i) {free(matrix[i]);}
        free(matrix);
        m = M.m;
        matrix = new double*[m];
        for(int i = 0; i < m; ++i){matrix[i]= new double[m];}
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = M.matrix[i][j];
            }
        }
    }
    return *this;
}

/// I create a function to fill up a matrix with random values.
void Matrix::fillMatrix() {
    for(int i=0; i<m; i++) {
        for(int j=0; j<m; j++) {
            matrix[i][j] = rand() % 101 - 50;
        }
    }
}

/// I create a function to print out a matrix.
void Matrix::printMatrix() const {
    std::cout << " New matrix " << std::endl;
    for (int i= 0 ; i<m ; i++) {
        for(int j = 0 ; j<m ; j++) {
            std::cout << matrix[i][j] << ", " << std::endl;
        }
    }
}

/// I create the addition operation between matrixes.
Matrix& Matrix::operator+=(const Matrix& M)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] += M.matrix[i][j];
        }
    }
    return *this;
}

/// I create the multiplication operator between matrixes.
Matrix& Matrix::operator*=(const Matrix& M)
{
    Matrix T(M.m);
    for (int i = 0; i < M.m; i++) {
        for (int j = 0; j < M.m; j++) {
            for (int k = 0; k < M.m; k++) {
                T.matrix[i][j] += (matrix[i][k] * M.matrix[k][j]);
            }
        }
    }
    return (*this = T);
}

/// I create the operator between matrix and constant.
Matrix& Matrix::operator*=(const double k)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] *= k;
        }
    }
    return *this;
}

/// I create the euclidean norm from a square matrix.
double Matrix::norm() const
{
    double norm, s=0.0;
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++)
            s += matrix[i][j]*matrix[i][j];
    }
    norm = sqrt(s);
    return norm;
}


/// I create an identity matrix.
void Matrix::identityMatrix() {
    for(int i=0; i<m; i++) {
        for(int j=0; j<m; j++) {
            if(i==j){
                matrix[i][j] = 1;
            }
            else {
                matrix[i][j] = 0;
            }
        }
    }
}

/// I create an array from a matrix.
double* Matrix::toArray() {
	double* arr = (double*) calloc(m * m, sizeof(double));
	int counter = 0;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			arr[counter] = matrix[j][i];
			counter++;
		};
	};
	return arr;
};
