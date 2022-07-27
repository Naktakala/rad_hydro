#include "chi_math_banded_solvers.h"

#include <stddef.h>
#include <iostream>
#include <cmath>

#include <map> //TODO: Remove

namespace chi_math
{
//###################################################################
/**Wrapper function for TDMA that takes a dense matrix and
 * a right-hand side.
 *
 * \param A A matrix with square size N.
 * \param b A rhs-vector of length N.
 *
 * \return X A vector of length N, the solution to the system.*/
std::vector<double> TDMA(const MatDbl& dense_A,
                         const VecDbl& b)
{
  const size_t N = dense_A.size();

  VecDbl  A(N,0.0);
  VecDbl  B(N,0.0);
  VecDbl  C(N,0.0);
  VecDbl  D = b;

  for (size_t i=0; i<N; ++i)
  {
    if (i>0)     A[i] = dense_A[i][i-1];
                 B[i] = dense_A[i][i];
    if (i<(N-1)) C[i] = dense_A[i][i+1];
  }

  return TDMA(A,B,C,D);
}

//###################################################################
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

//###################################################################
/**Wrapper function for the banded solver given a square matrix
 * A of square-size N, and a right hand side vector of size N.
 * All the bands in this matrix must be closely packed around the
 * diagonal (no gaps).
 *
 * \param dense_A A dense square matrix which size matches the rhs.
 * \param b       The right hand side vector.
 * \param num_bands_above_diag The number of bands above the diagonal.
 * \param num_bands_below_diag The number of bands below the diagonal.
 *
 * \return A solution vector matching the size of b.*/
std::vector<double> BandedSolver(const MatDbl& dense_A, const VecDbl& b,
                                 int num_bands_above_diag,
                                 int num_bands_below_diag)
{
  const size_t N  = b.size();
  const int& NbBD = num_bands_below_diag;
  const int& NbAD = num_bands_above_diag;
  const int  Nb   = NbAD + 1 + NbBD;

  std::vector<VecDbl> bands(Nb+1, VecDbl(N,0.0));
  bands[Nb] = b;

  for (int j=0; j<N; ++j)
  {
    int offset = -NbAD;
    for (int k=0; k<Nb; ++k)
    {
      int i = j + offset;
      if (i >= 0 and i < N)
        bands[k][j] = dense_A[i][j];
      ++offset;
    }//for k
  }//for i

  BandedSolver(bands, NbAD, NbBD);

  return bands.back();
}

//###################################################################
/**Banded solver given a square matrix of square-size N, and a right
 * hand side vector of size N.
 * All the bands in this matrix must be closely packed around the
 * diagonal (no gaps).
 *
 * \param bands Starting with the band furthest above the diagonal downward.
 *              The last band contains the rhs.
 * \param num_bands_above_diag The number of bands above the diagonal.
 * \param num_bands_below_diag The number of bands below the diagonal.
 *
 * \return The solution will be contained within the last band.*/
void BandedSolver(std::vector<VecDbl>& bands,
                  int num_bands_above_diag,
                  int num_bands_below_diag)
{
  const int& NbBD = num_bands_below_diag;
  const int& NbAD = num_bands_above_diag;
  const int  Nb   = NbAD + 1 + NbBD;

  auto&      rhs  = bands[Nb];
  const int  N    = static_cast<int>(rhs.size());

  const int diag_b_id = 0 + NbAD;
  const auto& diag_band = bands[diag_b_id];

  //============================================= Forward elimination
  for (int i=0; i<(N-1); ++i)
  {
    for (int jb=(diag_b_id+1); jb<Nb; ++jb)
    {
      double fac = bands[jb][i]/diag_band[i];
      for (int kb=jb; (kb-(jb-diag_b_id))>=0; --kb)
      {
        int j = i + jb - kb;
        bands[kb][j] -= fac*bands[kb-(jb-diag_b_id)][j];
      }
      rhs[i+jb-diag_b_id] -= fac*rhs[i];
    }//for jb
  }//for j

  //============================================= Backward substitution
  for (int i=N-1; i>=0; --i)
  {
    for (int jb=diag_b_id-1; jb>=0; --jb)
    {
      int j = i+diag_b_id-jb;
      if (j < N)
        rhs[i] -= bands[jb][j]*rhs[j];
    }//for jb
    rhs[i] /= diag_band[i];
  }//for i

}

}//namespace chi_math
