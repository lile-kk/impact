#include "CholeskyDecomposition.h"

namespace jama {
namespace impact {
/** Cholesky algorithm for symmetric and positive definite matrix.
   @param  A   Square, symmetric matrix.
   @return     Structure to access L and isspd flag.
   */

CholeskyDecomposition::CholeskyDecomposition (Matrix Arg) {
	// Initialize.
	double** A = Arg.getArray();
	n = Arg.getRowDimension();
	//动态开辟空间
    L = new double*[n]; //开辟行
	for(int i = 0; i < n; i++) {
		L[i] = new double[n]; //开辟列
    }
	isspd = (Arg.getColumnDimension() == n);
	// Main loop.
	for (int j = 0; j < n; j++) {
	    double* Lrowj = L[j];
		double d = 0.0;
		for (int k = 0; k < j; k++) {
		    double* Lrowk = L[k];
			double s = 0.0;
			for (int i = 0; i < k; i++) {
			   s += Lrowk[i]*Lrowj[i];
			}
			Lrowj[k] = s = (A[j][k] - s)/L[k][k];
			d = d + s*s;
			isspd = isspd & (A[k][j] == A[j][k]); 
		 }
		 d = A[j][j] - d;
		 isspd = isspd & (d > 0.0);
		 L[j][j] = std::sqrt(std::max(d,0.0));
		 for (int k = j+1; k < n; k++) {
			L[j][k] = 0.0;
		 }
	  }
   }  

/** Return triangular factor.
   @return     L
   */
Matrix CholeskyDecomposition::getL () {
    return *(new Matrix(L,n,n));
}  

/** Is the matrix symmetric and positive definite?
   @return     true if A is symmetric and positive definite.
   */
bool CholeskyDecomposition::isSPD () {
	return isspd;
}   

/** Solve A*X = B
   @param  B   A Matrix with as many rows as A and any number of columns.
   @return     X so that L*L'*X = B
   @exception  IllegalArgumentException  Matrix row dimensions must agree.
   @exception  RuntimeException  Matrix is not symmetric positive definite.
   */
  Matrix CholeskyDecomposition::solve (Matrix B) {
	if (B.getRowDimension() != n) {
		throw ("Matrix row dimensions must agree.");
	}
	if (!isspd) {
		throw ("Matrix is not symmetric positive definite.");
	}

	// Copy right hand side.
	double** X = B.getArrayCopy();
	int nx = B.getColumnDimension();

	// Solve L*Y = B;
	for (int k = 0; k < n; k++) {
		for (int i = k+1; i < n; i++) {
		    for (int j = 0; j < nx; j++) {
			   X[i][j] -= X[k][j]*L[i][k];
			}
		 }
		for (int j = 0; j < nx; j++) {
			X[k][j] /= L[k][k];
		}
	}

	// Solve L'*X = Y;
	for (int k = n-1; k >= 0; k--) {
		for (int j = 0; j < nx; j++) {
			X[k][j] /= L[k][k];
		}
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < nx; j++) {
			   X[i][j] -= X[k][j]*L[k][i];
			}
		}
	}
	return *(new Matrix(X,n,nx));
   }   

}
}