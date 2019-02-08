#include <exception>
#include <cstdlib>

#include "util/common_utils.cpp"
#include "CholeskyDecomposition.h"
#include "EigenvalueDecomposition.h"
#include "SingularValueDecomposition.h"
#include "LUDecompositio.h"
#include "QRDecomposition.h"
#include "util/Maths.h"
/**
   Jama = Java Matrix class.
<P>
   The Java Matrix Class provides the fundamental operations of numerical
   linear algebra.  Various constructors create Matrices from two dimensional
   arrays of double precision floating point numbers.  Various "gets" and
   "sets" provide access to submatrices and matrix elements.  Several methods 
   implement basic matrix arithmetic, including matrix addition and
   multiplication, matrix norms, and element-by-element array operations.
   Methods for reading and printing matrices are also included.  All the
   operations in this version of the Matrix Class involve real matrices.
   Complex matrices may be handled in a future version.
<P>
   Five fundamental matrix decompositions, which consist of pairs or triples
   of matrices, permutation vectors, and the like, produce results in five
   decomposition classes.  These decompositions are accessed by the Matrix
   class to compute solutions of simultaneous linear equations, determinants,
   inverses and other matrix functions.  The five decompositions are:
<P><UL>
   <LI>Cholesky Decomposition of symmetric, positive definite matrices.
   <LI>LU Decomposition of rectangular matrices.
   <LI>QR Decomposition of rectangular matrices.
   <LI>Singular Value Decomposition of rectangular matrices.
   <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
</UL>
<DL>
<DT><B>Example of use:</B></DT>
<P>
<DD>Solve a linear system A x = b and compute the residual norm, ||b - A x||.
<P><PRE>
	  double[][] vals = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
	  Matrix A = new Matrix(vals);
	  Matrix b = Matrix.random(3,1);
	  Matrix x = A.solve(b);
	  Matrix r = A.times(x).minus(b);
	  double rnorm = r.normInf();
</PRE></DD>
</DL>

@author The MathWorks, Inc. and the National Institute of Standards and Technology.
@version 5 August 1998
*/
namespace jama {
namespace impact {
class Matrix {
   /* ------------------------
   Class variables
 * ------------------------ */
private:
/** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
   */
   int m, n;

   /** Array for internal storage of elements.
   @serial internal array storage.
   */
   double** A;

public:
   Matrix(int m, int n);
   Matrix(int m, int n, double s);
   Matrix (double vals[], int m, int size); // 需要船数组的大小
   //Matrix (double** A); // 这个不应该支持
   Matrix (double** A, int m, int n);
   ~Matrix();

public:
   Matrix arrayLeftDivide (Matrix B);
   double** getArray ();
   Matrix arrayLeftDivideEquals (Matrix B);
   Matrix arrayRightDivide (Matrix B);
   Matrix arrayRightDivideEquals (Matrix B);
   Matrix arrayTimes (Matrix B);
   Matrix arrayTimesEquals (Matrix B);
   int getRowDimension();
   int getColumnDimension ();
   double** getArrayCopy ();
   Matrix clone ();
   double cond ();
   Matrix copy ();
   static Matrix constructWithCopy(double** A, int m, int n); // m,n
   double det ();
   EigenvalueDecomposition eig ();
   double get (int i, int j);
   Matrix getMatrix (int* r, int r_size, int j0, int j1); // r_size
   double** getArray ();
   double* getColumnPackedCopy ();
   Matrix getMatrix (int i0, int i1, int j0, int j1);
   Matrix getMatrix (int i0, int i1, int* c, int c_size); //c_size
   Matrix getMatrix (int* r, int* c, int r_size, int c_size); //r_size, c_size
   double* getRowPackedCopy ();
   Matrix inverse ();
   double length();
   LUDecomposition lu ();
   Matrix minus (Matrix B);
   Matrix minusEquals (Matrix B);
   double norm1 ();
   double norm2 ();
   double normF ();
   double normInf ();
   Matrix plus (Matrix B);
   Matrix plusEquals (Matrix B);
   QRDecomposition qr ();
   int rank ();
   void set (int i, int j, double s);
   void setMatrix (int i0, int i1, int j0, int j1, Matrix X);
   void setMatrix (int i0, int i1, int* c, int c_size, Matrix X); // c_size
   void setMatrix (int* r, int r_size, int j0, int j1, Matrix X); // r_size
   void setMatrix (int* r, int* c, int r_size, int c_size, Matrix X); // r_size, c_size
   Matrix solve (Matrix B);
   Matrix solveTranspose (Matrix B);
   SingularValueDecomposition svd ();
   Matrix times (double s);
   Matrix times (Matrix B);
   Matrix timesEquals (double s);
   double trace ();
   Matrix transpose ();
   Matrix uminus ();
   Matrix vectorProduct(Matrix B);

// static func
public:
   // m, n
   static Matrix constructWithCopy(double** A, int m, int n) {
      Matrix* X = new Matrix(m,n);
	   double** C = X->getArray();
	   for (int i = 0; i < m; i++) {
		// if (A[i].length != n) {
		// 	throw new IllegalArgumentException
		// 	   ("All rows must have the same length.");
		// }
		   for (int j = 0; j < n; j++) {
			   C[i][j] = A[i][j];
		   }
	   }
	  return *X;
   }
   /** Generate identity matrix
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
   */
   static Matrix identity (int m, int n) {
	  Matrix* A = new Matrix(m,n);
	  double** X = A->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			X[i][j] = (i == j ? 1.0 : 0.0);
		 }
	  }
	  return *A;
   }  
/** Generate matrix with random elements
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with uniformly distributed random elements.
   */
static Matrix random (int m, int n) {
	  Matrix* A = new Matrix(m,n);
	  double** X = A->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			X[i][j] = std::rand();
		 }
	  }
	  return *A;
   } 
//TODO
static Matrix read(std::string input_path) {

}   

private:
   void checkMatrixDimensions (Matrix const B);
   CholeskyDecomposition chol ();
   

};
}// namespace impact
}// namespace jama
class Matrix {}