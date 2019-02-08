#include "Matrix.h"

namespace jama {
namespace impact {
Matrix::Matrix(int m_t, int n_t) {
    m = m_t;
    n = n_t;
    //动态开辟空间
    A = new double*[m]; //开辟行
	for(int i = 0; i < m; i++){
		A[i] = new double[n]; //开辟列
    }
}

/** Construct an m-by-n constant matrix.
   @param m    Number of rows.
   @param n    Number of colums.
   @param s    Fill the matrix with this scalar value.
   */
Matrix::Matrix(int m_t, int n_t, double s_t) {
    m = m_t;
    n = n_t;
    //动态开辟空间
    A = new double*[m]; //开辟行
	for(int i = 0; i < m; i++) {
		A[i] = new double[n]; //开辟列
    }
    for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			A[i][j] = s_t;
		 }
	  }
}

/** Construct a matrix from a one-dimensional packed array
   @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
   @param m    Number of rows.
   @exception  IllegalArgumentException Array length must be a multiple of m.
   */
Matrix::Matrix (double vals_t[], int m_t, int size) {
    m = m_t;
    n = (m != 0 ? size/m : 0);
    if (m*n != size) {
        throw ("Array length must be a multiple of m.");
    }
    //动态开辟空间
    A = new double*[m]; //开辟行
	for(int i = 0; i < m; i++) {
		A[i] = new double[n]; //开辟列
    }
	for (int i = 0; i < m; i++) {
	    for (int j = 0; j < n; j++) {
			A[i][j] = vals_t[i+j*m];
		 }
	  }
}

/** Construct a matrix from a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
   @see        #constructWithCopy
   */
// Matrix::Matrix (double** A_t) {
//     m = jama::util::get_array_len(A_t);
//     n = jama::util::get_array_len(A_t[0]);
//     for (int i = 0; i < m; i++) {
//         if (jama::util::get_array_len(A_t[i]) != n) {
//             throw("All rows must have the same length.");
//         }
//     }
//     A = A_t;
// }

/** Construct a matrix quickly without checking arguments.
   @param A    Two-dimensional array of doubles.
   @param m    Number of rows.
   @param n    Number of colums.
   */
Matrix::Matrix (double** A_t, int m_t, int n_t) {
    A = A_t;
    m = m_t;
    n = n_t;
}

Matrix::~Matrix() {
    //释放开辟的资源
    for(int i = 0; i < m; i++)
		delete[] A[i];
	delete[] A;
}

/** Element-by-element left division, C = A.\B
	 @param B    another matrix
	 @return     A.\B
	 */
Matrix Matrix::arrayLeftDivide (Matrix B) {
	checkMatrixDimensions(B);
	Matrix* X = new Matrix(m,n);
	double** C = X->getArray();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		    C[i][j] = B.A[i][j] / A[i][j];
		}
	}
	return *X;
} 

 /** Element-by-element left division in place, A = A.\B
   @param B    another matrix
   @return     A.\B
   */
Matrix Matrix::arrayLeftDivideEquals (Matrix B) {
	checkMatrixDimensions(B);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		A[i][j] = B.A[i][j] / A[i][j];
		}
	}
	return *this;
}   

/** Element-by-element right division, C = A./B
   @param B    another matrix
   @return     A./B
   */
Matrix Matrix::arrayRightDivide (Matrix B) {
	checkMatrixDimensions(B);
	Matrix* X = new Matrix(m,n);
	double** C = X->getArray();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		    C[i][j] = A[i][j] / B.A[i][j];
		}
	}
	return *X;
}   

/** Element-by-element right division in place, A = A./B
   @param B    another matrix
   @return     A./B
   */
Matrix Matrix::arrayRightDivideEquals (Matrix B) {
	checkMatrixDimensions(B);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		    A[i][j] = A[i][j] / B.A[i][j];
		}
	}
	return *this;
}   

/** Element-by-element multiplication, C = A.*B
   @param B    another matrix
   @return     A.*B
   */

Matrix Matrix::arrayTimes (Matrix B) {
	checkMatrixDimensions(B);
	Matrix* X = new Matrix(m,n);
	double** C = X->getArray();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		    C[i][j] = A[i][j] * B.A[i][j];
		}
	}
	return *X;
}   

  /** Element-by-element multiplication in place, A = A.*B
   @param B    another matrix
   @return     A.*B
   */
  Matrix Matrix::arrayTimesEquals (Matrix B) {
	checkMatrixDimensions(B);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		    A[i][j] = A[i][j] * B.A[i][j];
		}
	}
	return *this;
}   

/** Cholesky Decomposition
   @return     CholeskyDecomposition
   @see CholeskyDecomposition
   */
CholeskyDecomposition Matrix::chol () {
	return *(new CholeskyDecomposition(*this));
}  

/** Clone the Matrix object.
   */
Matrix Matrix::clone () {
	return this->copy();
}   
   
/** Matrix condition (2 norm)
@return     ratio of largest to smallest singular value.
*/
double Matrix::cond () {
   SingularValueDecomposition* singular = new SingularValueDecomposition(*this);
   return singular->cond();
}  

/** Make a deep copy of a matrix
   */
Matrix Matrix::copy () {
	Matrix* X = new Matrix(m,n);
	double** C = X->getArray();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		    C[i][j] = A[i][j];
		}
	}
	return *X;
}   

/** Matrix determinant
   @return     determinant
   */
double Matrix::det () {
	LUDecomposition* luDecomp = new LUDecomposition(*this);
	return luDecomp->det();
}   
   /** Eigenvalue Decomposition
   @return     EigenvalueDecomposition
   @see EigenvalueDecomposition
   */
 EigenvalueDecomposition Matrix::eig () {
	 EigenvalueDecomposition* eiDecomp = new EigenvalueDecomposition(*this);
	return *eiDecomp;
}   
   /** Get a single element.
   @param i    Row index.
   @param j    Column index.
   @return     A(i,j)
   @exception  ArrayIndexOutOfBoundsException
   */
double Matrix::get (int i, int j) {
	return A[i][j];
}   
   /** Access the internal two-dimensional array.
   @return     Pointer to the two-dimensional array of matrix elements.
   */

/** Get row dimension.
   @return     m, the number of rows.
   */
int Matrix::getRowDimension () {
	return m;
}  

/** Get column dimension.
   @return     n, the number of columns.
   */
int Matrix::getColumnDimension () {
	return n;
} 

/** Copy the internal two-dimensional array.
   @return     Two-dimensional array copy of matrix elements.
   */
double** Matrix::getArrayCopy () {
    //动态开辟空间
    double** C = new double*[m]; //开辟行
	for(int i = 0; i < m; i++) {
		C[i] = new double[n]; //开辟列
    }
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
		    C[i][j] = A[i][j];
		}
	}
	return C;
}  

/** Access the internal two-dimensional array.
   @return     Pointer to the two-dimensional array of matrix elements.
   */
double** Matrix::getArray () {
	  return A;
}  

/** Make a one-dimensional column packed copy of the internal array.
   @return     Matrix elements packed in a one-dimensional array by columns.
   */
double* Matrix::getColumnPackedCopy () {
	  double* vals = new double[m*n];
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			vals[i+j*m] = A[i][j];
		 }
	  }
	  return vals;
} 

/** Get a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @return     A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
 Matrix Matrix::getMatrix (int i0, int i1, int j0, int j1) {
	  Matrix* X = new Matrix(i1-i0+1,j1-j0+1);
	  double** B = X->getArray();
	  try {
		 for (int i = i0; i <= i1; i++) {
			for (int j = j0; j <= j1; j++) {
			   B[i-i0][j-j0] = A[i][j];
			}
		 }
	  } catch(std::exception& e) {
		std::cerr << "exception caught: " << e.what() << '\n'; 
	  }
	  return *X;
}

/** Get a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @return     A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
Matrix Matrix::getMatrix (int i0, int i1, int* c, int c_size) {
	  Matrix* X = new Matrix(i1-i0+1,c_size);
	  double** B = X->getArray();
	  try {
		 for (int i = i0; i <= i1; i++) {
			for (int j = 0; j < c_size; j++) {
			   B[i-i0][j] = A[i][c[j]];
			}
		 }
	  } catch(std::exception& e) {
		  std::cerr << "exception caught: " << e.what() << '\n';
	  }
	  return *X;
   }   

/** Get a submatrix.
   @param r    Array of row indices.
   @param i0   Initial column index
   @param i1   Final column index
   @return     A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
Matrix Matrix::getMatrix (int* r, int r_size, int j0, int j1) {
	  Matrix* X = new Matrix(r_size,j1-j0+1);
	  double** B = X->getArray();
	  try {
		 for (int i = 0; i < r_size; i++) {
			for (int j = j0; j <= j1; j++) {
			   B[i][j-j0] = A[r[i]][j];
			}
		 }
	  } catch(std::exception& e) {
		 std::cerr << "exception caught: " << e.what() << '\n';
	  }
	  return *X;
}

/** Get a submatrix.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @return     A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
Matrix Matrix::getMatrix (int* r, int* c, int r_size, int c_size) {
	  Matrix* X = new Matrix(r_size, c_size);
	  double** B = X->getArray();
	  try {
		 for (int i = 0; i < r_size; i++) {
			for (int j = 0; j < c_size; j++) {
			   B[i][j] = A[r[i]][c[j]];
			}
		 }
	  } catch(std::exception& e) {
		 std::cerr << "exception caught: " << e.what() << '\n';
	  }
	  return *X;
}

/** Make a one-dimensional row packed copy of the internal array.
   @return     Matrix elements packed in a one-dimensional array by rows.
   */
double* Matrix::getRowPackedCopy () {
	  double* vals = new double[m*n];
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			vals[i*n+j] = A[i][j];
		 }
	  }
	return vals;
} 


   /** Matrix inverse or pseudoinverse
   @return     inverse(A) if A is square, pseudoinverse otherwise.
   */
Matrix Matrix::inverse () {
	  return solve(identity(m,m));
}   
 /**
 * Returns the length of the vector
 * Creation date: (2001-10-21 20.57.43)
 * @return double
 * @exception java.lang.IllegalArgumentException The exception description.
 */
double Matrix::length() {
	double length;
	if (this->getColumnDimension() != 1)
		throw ("Vector Lengths can only be applied to vectors. There can only be one column.");
	if (this->getRowDimension() != 3)
		throw ("The length method currently applies only to three dimensional vectors");
//
	length = std::sqrt(this->get(0,0)*this->get(0,0) + this->get(1,0)*this->get(1,0)+this->get(2,0)*this->get(2,0));
		
	return length;
}
   /** LU Decomposition
   @return     LUDecomposition
   @see LUDecomposition
   */
LUDecomposition Matrix::lu () {
	  return *(new LUDecomposition(*this));
   }   
   /** C = A - B
   @param B    another matrix
   @return     A - B
   */
Matrix Matrix::minus (Matrix B) {
	  checkMatrixDimensions(B);
	  Matrix* X = new Matrix(m,n);
	  double** C = X->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			C[i][j] = A[i][j] - B.A[i][j];
		 }
	  }
	  return *X;
   }   
   /** A = A - B
   @param B    another matrix
   @return     A - B
   */
Matrix Matrix::minusEquals (Matrix B) {
	  checkMatrixDimensions(B);
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			A[i][j] = A[i][j] - B.A[i][j];
		 }
	  }
	  return *this;
   }   
   /** One norm
   @return    maximum column sum.
   */
double Matrix::norm1 () {
	  double f = 0;
	  for (int j = 0; j < n; j++) {
		 double s = 0;
		 for (int i = 0; i < m; i++) {
			s += std::abs(A[i][j]);
		 }
		 f = std::max(f,s);
	  }
	  return f;
   }   
   /** Two norm
   @return    maximum singular value.
   */
double Matrix::norm2 () {
	  return ((new SingularValueDecomposition(*this))->norm2());
}   
   /** Frobenius norm
   @return    sqrt of sum of squares of all elements.
   */
double Matrix::normF () {
	  double f = 0;
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			f = jama::util::Maths::hypot(f,A[i][j]);
		 }
	  }
	  return f;
   }   
   /** Infinity norm
   @return    maximum row sum.
   */
double Matrix::normInf () {
	  double f = 0;
	  for (int i = 0; i < m; i++) {
		 double s = 0;
		 for (int j = 0; j < n; j++) {
			s += std::abs(A[i][j]);
		 }
		 f = std::max(f,s);
	  }
	  return f;
   }   
   /** C = A + B
   @param B    another matrix
   @return     A + B
   */
Matrix Matrix::plus (Matrix B) {
	  checkMatrixDimensions(B);
	  Matrix* X = new Matrix(m,n);
	  double** C = X->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			C[i][j] = A[i][j] + B.A[i][j];
		 }
	  }
	  return *X;
   }   
   /** A = A + B
   @param B    another matrix
   @return     A + B
   */
Matrix Matrix::plusEquals (Matrix B) {
	  checkMatrixDimensions(B);
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			A[i][j] = A[i][j] + B.A[i][j];
		 }
	  }
	  return *this;
}
/** QR Decomposition
   @return     QRDecomposition
   @see QRDecomposition
   */
QRDecomposition Matrix::qr () {
	  return *(new QRDecomposition(*this));
   }   

/** Matrix rank
   @return     effective numerical rank, obtained from SVD.
   */
int Matrix::rank () {
	  return (new SingularValueDecomposition(*this))->rank();
}   

void Matrix::set (int i, int j, double s) {
	  A[i][j] = s;
   }   
   /** Set a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
void Matrix::setMatrix (int i0, int i1, int j0, int j1, Matrix X) {
	  try {
		 for (int i = i0; i <= i1; i++) {
			for (int j = j0; j <= j1; j++) {
			   A[i][j] = X.get(i-i0,j-j0);
			}
		 }
	  } catch(std::exception& e) {
		 std::cerr << "exception caught: " << e.what() << '\n';
	  }
}

 /** Set a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @param X    A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
void Matrix::setMatrix (int i0, int i1, int* c, int c_size, Matrix X) {
	  try {
		 for (int i = i0; i <= i1; i++) {
			for (int j = 0; j < c_size; j++) {
			   A[i][c[j]] = X.get(i-i0,j);
			}
		 }
	  } catch(std::exception& e) {
		 std::cerr << "exception caught: " << e.what() << '\n';
	  }
   }   
   /** Set a submatrix.
   @param r    Array of row indices.
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
void Matrix::setMatrix (int* r, int r_size, int j0, int j1, Matrix X) {
	  try {
		 for (int i = 0; i < r_size; i++) {
			for (int j = j0; j <= j1; j++) {
			   A[r[i]][j] = X.get(i,j-j0);
			}
		 }
	  } catch(std::exception& e) {
		std::cerr << "exception caught: " << e.what() << '\n';	  
		}
   }   
   /** Set a submatrix.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @param X    A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */
void Matrix::setMatrix (int* r, int* c, int r_size, int c_size, Matrix X) {
	  try {
		 for (int i = 0; i < r_size; i++) {
			for (int j = 0; j < c_size; j++) {
			   A[r[i]][c[j]] = X.get(i,j);
			}
		 }
	  } catch(std::exception& e) {
		 std::cerr << "exception caught: " << e.what() << '\n';
	  }
   }   
   /** Solve A*X = B
   @param B    right hand side
   @return     solution if A is square, least squares solution otherwise
   */
Matrix Matrix::solve (Matrix B) {
	  return (m == n ? (new LUDecomposition(*this))->solve(B) :
					   (new QRDecomposition(*this))->solve(B));
   }   
   /** Solve X*A = B, which is also A'*X' = B'
   @param B    right hand side
   @return     solution if A is square, least squares solution otherwise.
   */
Matrix Matrix::solveTranspose (Matrix B) {
	  return transpose().solve(B.transpose());
   }

/** Singular Value Decomposition
   @return     SingularValueDecomposition
   @see SingularValueDecomposition
   */
SingularValueDecomposition Matrix::svd () {
	  return *(new SingularValueDecomposition(*this));
   }   
   /** Multiply a matrix by a scalar, C = s*A
   @param s    scalar
   @return     s*A
   */
Matrix Matrix::times (double s) {
	  Matrix* X = new Matrix(m,n);
	  double** C = X->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			C[i][j] = s*A[i][j];
		 }
	  }
	  return *X;
}    

  /** Linear algebraic matrix multiplication, A * B
   @param B    another matrix
   @return     Matrix product, A * B
   @exception  IllegalArgumentException Matrix inner dimensions must agree.
   */
Matrix Matrix::times (Matrix B) {
	  if (B.m != n) {
		 throw ("Matrix inner dimensions must agree.");
	  }
	  Matrix* X = new Matrix(m,B.n);
	  double** C = X->getArray();
	  double* Bcolj = new double[n];
	  for (int j = 0; j < B.n; j++) {
		 for (int k = 0; k < n; k++) {
			Bcolj[k] = B.A[k][j];
		 }
		 for (int i = 0; i < m; i++) {
			double* Arowi = A[i];
			double s = 0;
			for (int k = 0; k < n; k++) {
			   s += Arowi[k]*Bcolj[k];
			}
			C[i][j] = s;
		 }
	  }
	  return *X;
   }   
   /** Multiply a matrix by a scalar in place, A = s*A
   @param s    scalar
   @return     replace A by s*A
   */
Matrix Matrix::timesEquals (double s) {
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			A[i][j] = s*A[i][j];
		 }
	  }
	  return *this;
   }   
   /** Matrix trace.
   @return     sum of the diagonal elements.
   */
double Matrix::trace () {
	  double t = 0;
	  for (int i = 0; i < std::min(m,n); i++) {
		 t += A[i][i];
	  }
	  return t;
   }   
   /** Matrix transpose.
   @return    A'
   */
Matrix Matrix::transpose () {
	  Matrix* X = new Matrix(n,m);
	  double** C = X->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			C[j][i] = A[i][j];
		 }
	  }
	  return *X;
   }   
   /**  Unary minus
   @return    -A
   */
Matrix Matrix::uminus () {
	  Matrix* X = new Matrix(m,n);
	  double** C = X->getArray();
	  for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++) {
			C[i][j] = -A[i][j];
		 }
	  }
	  return *X;
}

/**
 * This method calculates the vectorProduct between two matrices. The matrices must be three dimensional.
 * Creation date: (2001-10-21 19.46.07) Jonas Forssell
 * @return result_vector Jama.Matrix
 * @param a Jama.Matrix
 * @param b Jama.Matrix
 * @exception java.lang.IllegalArgumentException The exception description.
 */
Matrix Matrix::vectorProduct(Matrix B) {
	//
	Matrix* result_vector = new Matrix(3,1);
	//
	if ((B.getColumnDimension() != 1) || (this->getColumnDimension() != 1))
		throw ("Vector Products can only be applied to vectors. There can only be one column.");
	if ((B.getRowDimension() != 3) || (this->getRowDimension() != 3))
		throw ("Matrix inner dimensions must agree.");
	//
	// Calculate the vector product term by term.
	result_vector->set(0, 0, this->get(1, 0)*B.get(2, 0)-B.get(1, 0)*this->get(2, 0));
	result_vector->set(1, 0, this->get(2, 0)*B.get(0, 0)-B.get(2, 0)*this->get(0, 0));
	result_vector->set(2, 0, this->get(0, 0)*B.get(1, 0)-B.get(0, 0)*this->get(1, 0));
	// Now, return the completed vector.
	return *result_vector;
}



/******** private method *************/

/** Access the internal two-dimensional array.
   @return     Pointer to the two-dimensional array of matrix elements.
   */
double** Matrix::getArray () {
    return A;
}     

/** Check if size(A) == size(B) **/
void Matrix::checkMatrixDimensions (Matrix B) {
	if (B.m != m || B.n != n) {
		throw("Matrix dimensions must agree.");
	}
}   
}// namespace impact
}// namespace jama