/** LU Decomposition.
   <P>
   For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
   unit lower triangular matrix L, an n-by-n upper triangular matrix U,
   and a permutation vector piv of length m so that A(piv,:) = L*U.
   If m < n, then L is m-by-m and U is m-by-n.
   <P>
   The LU decompostion with pivoting always exists, even if the matrix is
   singular, so the constructor will never fail.  The primary use of the
   LU decomposition is in the solution of square systems of simultaneous
   linear equations.  This will fail if isNonsingular() returns false.
   */
#pragma once

#include <cmath>
#include <algorithm>
#include "Matrix.h"

namespace jama{
namespace impact {
class LUDecomposition {
public:
    LUDecomposition(Matrix A);
    ~LUDecomposition();
    double det ();
    double* getDoublePivot ();
    int getPivLen();
    Matrix getL ();
    int* getPivot ();
    Matrix getU ();
    bool isNonsingular ();
    Matrix solve (Matrix B);


private:
    /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   double** LU;

   /** Row and column dimensions, and pivot sign.
   @serial column dimension.
   @serial row dimension.
   @serial pivot sign.
   */
   int m, n, pivsign; 

   /** Internal storage of pivot vector.
   @serial pivot vector.
   */
   int* piv;
};
}
}