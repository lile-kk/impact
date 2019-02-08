/** QR Decomposition.
<P>
   For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
   orthogonal matrix Q and an n-by-n upper triangular matrix R so that
   A = Q*R.
<P>
   The QR decompostion always exists, even if the matrix does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.  This will fail if isFullRank()
   returns false.
*/
#pragma once

#include "Matrix.h"

namespace jama {
namespace impact {
class QRDecomposition {
private:
    /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   double** QR;

   /** Row and column dimensions.
   @serial column dimension.
   @serial row dimension.
   */
   int m, n;

   /** Array for internal storage of diagonal of R.
   @serial diagonal of R.
   */
   double* Rdiag;

public:
    QRDecomposition(Matrix A);
    ~QRDecomposition();
    Matrix getH ();
    Matrix getQ ();
    Matrix getR ();
    bool isFullRank ();
    Matrix solve (Matrix B);

};
}
}