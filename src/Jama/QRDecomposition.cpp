#include "QRDecomposition.h"

namespace jama {
namespace impact {
/** QR Decomposition, computed by Householder reflections.
   @param A    Rectangular matrix
   @return     Structure to access R and the Householder vectors and compute Q.
   */
QRDecomposition::QRDecomposition (Matrix A) {
	  // Initialize.
	  QR = A.getArrayCopy();
	  m = A.getRowDimension();
	  n = A.getColumnDimension();
	  Rdiag = new double[n];

	  // Main loop.
	  for (int k = 0; k < n; k++) {
		 // Compute 2-norm of k-th column without under/overflow.
		 double nrm = 0;
		 for (int i = k; i < m; i++) {
			nrm = jama::util::Maths::hypot(nrm,QR[i][k]);
		 }

		 if (nrm != 0.0) {
			// Form k-th Householder vector.
			if (QR[k][k] < 0) {
			   nrm = -nrm;
			}
			for (int i = k; i < m; i++) {
			   QR[i][k] /= nrm;
			}
			QR[k][k] += 1.0;

			// Apply transformation to remaining columns.
			for (int j = k+1; j < n; j++) {
			   double s = 0.0; 
			   for (int i = k; i < m; i++) {
				  s += QR[i][k]*QR[i][j];
			   }
			   s = -s/QR[k][k];
			   for (int i = k; i < m; i++) {
				  QR[i][j] += s*QR[i][k];
			   }
			}
		 }
		 Rdiag[k] = -nrm;
	  }
}

/** Return the Householder vectors
   @return     Lower trapezoidal matrix whose columns define the reflections
   */
Matrix QRDecomposition::getH () {
	  Matrix* X = new Matrix(m,n);
	  double** H = X->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			if (i >= j) {
			   H[i][j] = QR[i][j];
			} else {
			   H[i][j] = 0.0;
			}
		 }
	  }
	  return *X;
}   
   /** Generate and return the (economy-sized) orthogonal factor
   @return     Q
   */
Matrix QRDecomposition::getQ () {
	  Matrix* X = new Matrix(m,n);
	  double** Q = X->getArray();
	  for (int k = n-1; k >= 0; k--) {
		 for (int i = 0; i < m; i++) {
			Q[i][k] = 0.0;
		 }
		 Q[k][k] = 1.0;
		 for (int j = k; j < n; j++) {
			if (QR[k][k] != 0) {
			   double s = 0.0;
			   for (int i = k; i < m; i++) {
				  s += QR[i][k]*Q[i][j];
			   }
			   s = -s/QR[k][k];
			   for (int i = k; i < m; i++) {
				  Q[i][j] += s*QR[i][k];
			   }
			}
		 }
	  }
	  return *X;
}   
   /** Return the upper triangular factor
   @return     R
   */
Matrix QRDecomposition::getR () {
	  Matrix* X = new Matrix(n,n);
	  double** R = X->getArray();
	  for (int i = 0; i < n; i++) {
		 for (int j = 0; j < n; j++) {
			if (i < j) {
			   R[i][j] = QR[i][j];
			} else if (i == j) {
			   R[i][j] = Rdiag[i];
			} else {
			   R[i][j] = 0.0;
			}
		 }
	  }
	  return *X;
} 
/** Is the matrix full rank?
   @return     true if R, and hence A, has full rank.
   */
bool QRDecomposition::isFullRank () {
	  for (int j = 0; j < n; j++) {
		 if (Rdiag[j] == 0)
			return false;
	  }
	  return true;
   }   
   /** Least squares solution of A*X = B
   @param B    A Matrix with as many rows as A and any number of columns.
   @return     X that minimizes the two norm of Q*R*X-B.
   @exception  IllegalArgumentException  Matrix row dimensions must agree.
   @exception  RuntimeException  Matrix is rank deficient.
   */
Matrix QRDecomposition::solve (Matrix B) {
	  if (B.getRowDimension() != m) {
		 throw ("Matrix row dimensions must agree.");
	  }
	  if (!this->isFullRank()) {
		 throw ("Matrix is rank deficient.");
	  }
	  
	  // Copy right hand side
	  int nx = B.getColumnDimension();
	  double** X = B.getArrayCopy();

	  // Compute Y = transpose(Q)*B
	  for (int k = 0; k < n; k++) {
		 for (int j = 0; j < nx; j++) {
			double s = 0.0; 
			for (int i = k; i < m; i++) {
			   s += QR[i][k]*X[i][j];
			}
			s = -s/QR[k][k];
			for (int i = k; i < m; i++) {
			   X[i][j] += s*QR[i][k];
			}
		 }
	  }
	  // Solve R*X = Y;
	  for (int k = n-1; k >= 0; k--) {
		 for (int j = 0; j < nx; j++) {
			X[k][j] /= Rdiag[k];
		 }
		 for (int i = 0; i < k; i++) {
			for (int j = 0; j < nx; j++) {
			   X[i][j] -= X[k][j]*QR[i][k];
			}
		 }
	  }
	  return ((new Matrix(X,n,nx))->getMatrix(0,n-1,0,nx-1));
}     

QRDecomposition::~QRDecomposition() {
    //释放开辟的资源
    for(int i = 0; i < m; i++){
        delete[] QR[i];
    }
    delete[] QR;
    delete[] Rdiag;
}
}
}