#pragma once

#include <cmath>
#include <algorithm>

#include "Matrix.h"

/** Cholesky Decomposition.
   <P>
   For a symmetric, positive definite matrix A, the Cholesky decomposition
   is an lower triangular matrix L so that A = L*L'.
   <P>
   If the matrix is not symmetric or positive definite, the constructor
   returns a partial decomposition and sets an internal flag that may
   be queried by the isSPD() method.
   */
namespace jama {
namespace impact {
class CholeskyDecomposition {
public:
    CholeskyDecomposition (Matrix Arg);
    Matrix getL ();
    bool isSPD ();
    Matrix solve (Matrix B);

private:
 /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   double** L;

   /** Row and column dimension (square matrix).
   @serial matrix dimension.
   */
   int n;

   /** Symmetric and positive definite flag.
   @serial is symmetric and positive definite flag.
   */
   bool isspd;

};
}
}