#pragma once

#include <cmath>

#include "Matrix.h"
#include "algorithm"

namespace jama {
namespace impact {
class SingularValueDecomposition {
public:
    SingularValueDecomposition(Matrix Arg);
    ~SingularValueDecomposition();
    double cond ();
    Matrix getS ();
    double* getSingularValues ();
    Matrix getU ();
    Matrix getV ();
    double norm2 ();
    int rank ();
    int get_s_size(); // s数组大小


private:
    /** Arrays for internal storage of U and V.
   @serial internal storage of U.
   @serial internal storage of V.
   */
    double** U;
    double** V;

   /** Array for internal storage of singular values.
   @serial internal storage of singular values.
   */
    double* s;

   /** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
   */
    int m, n;
};
}
}