
#pragma once 
#include "Matrix.h"
#include <algorithm>
#include <cmath>
/** Eigenvalues and eigenvectors of a real matrix. 
<P>
	If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
	diagonal and the eigenvector matrix V is orthogonal.
	I.e. A = V.times(D.times(V.transpose())) and 
	V.times(V.transpose()) equals the identity matrix.
<P>
	If A is not symmetric, then the eigenvalue matrix D is block diagonal
	with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
	lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
	columns of V represent the eigenvectors in the sense that A*V = V*D,
	i.e. A.times(V) equals V.times(D).  The matrix V may be badly
	conditioned, or even singular, so the validity of the equation
	A = V*D*inverse(V) depends upon V.cond().
**/

namespace jama {
namespace impact {
class EigenvalueDecomposition {
private:
    /** Row and column dimension (square matrix).
   @serial matrix dimension.
   */
   int n;
   /** Symmetry flag.
   @serial internal symmetry flag.
   */
   bool issymmetric;
   /** Arrays for internal storage of eigenvalues.
   @serial internal storage of eigenvalues.
   */
   double* d;
   double* e;

   /** Array for internal storage of eigenvectors.
   @serial internal storage of eigenvectors.
   */
   double** V;

   /** Array for internal storage of nonsymmetric Hessenberg form.
   @serial internal storage of nonsymmetric Hessenberg form.
   */
   double** H;

   /** Working storage for nonsymmetric algorithm.
   @serial working storage for nonsymmetric algorithm.
   */
   double* ort;

   // Complex scalar division.
    double cdivr, cdivi;

public:
    EigenvalueDecomposition(Matrix Arg);
    ~EigenvalueDecomposition();
    Matrix getD ();
    double* getImagEigenvalues ();
    double* getRealEigenvalues ();
    Matrix getV ();

private:
    void orthes ();
    void tred2 ();
    void tql2 ();
    void hqr2 ();
    void cdiv(double xr, double xi, double yr, double yi);

};
}
}