#ifndef CHI_MATH_BANDED_SOLVERS_H
#define CHI_MATH_BANDED_SOLVERS_H

#include <vector>

namespace chi_math
{
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;
  std::vector<double> TDMA(const MatDbl& A, const VecDbl& b);
  std::vector<double> TDMA(std::vector<double>& A,
                           std::vector<double>& B,
                           std::vector<double>& C,
                           std::vector<double>& D);

  std::vector<double> BandedSolver(const MatDbl& A, const VecDbl& b,
                                   int num_bands_above_diag,
                                   int num_bands_below_diag);

  void BandedSolver(std::vector<VecDbl>& bands,
                    int num_bands_above_diag,
                    int num_bands_below_diag);
}


#endif //CHI_MATH_BANDED_SOLVERS_H
