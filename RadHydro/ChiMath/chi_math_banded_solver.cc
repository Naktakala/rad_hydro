#include "chi_math_banded_solvers.h"

#include <stddef.h>

namespace chi_math
{
  /**Wrapper function for TDMA that takes a dense matrix and
   * a right-hand side.
   *
   * \param A A matrix with square size N.
   * \param b A rhs-vector of length N.
   *
   * \return X A vector of length N, the solution to the system.*/
  std::vector<double> TDMA(MatDbl& dense_A,
                           VecDbl& b)
  {
    const size_t N = dense_A.size();

    VecDbl  A(N,0.0);
    VecDbl  B(N,0.0);
    VecDbl  C(N,0.0);
    VecDbl& D = b;

    for (size_t i=0; i<N; ++i)
    {
      if (i>0)     A[i] = dense_A[i][i-1];
                   B[i] = dense_A[i][i];
      if (i<(N-1)) C[i] = dense_A[i][i+1];
    }

    return TDMA(A,B,C,D);
  }

  /**Thomas algorithm for Gaussian elimination of a
   * tridiagonal system. Assume the system size is N;
   *
   * \param A A vector length N for the band below the diagonal.
   * \param B A vector length N for the diagonal.
   * \param C A vector length N for the band above the diagonal.
   * \param D A vector length N for the right-hand side.
   *
   * \return X A vector of length N, the solution to the system.*/
  std::vector<double> TDMA(std::vector<double>& A,
                           std::vector<double>& B,
                           std::vector<double>& C,
                           std::vector<double>& D)
  {
    const int N = static_cast<int>(D.size());
    std::vector<double> X(N,0.0);

    for (int i=1; i<N; ++i)
    {
      double W = A[i] / B[i-1];
      B[i] = B[i] - W*C[i-1];
      D[i] = D[i] - W*D[i-1];
    }

    X[N-1] = D[N-1] / B[N-1];
    for (int i=(N-2); i>=0; --i)
    {
      X[i] = (D[i] - C[i]*X[i+1])/B[i];
    }

    return X;
  }
}//namespace chi_math
