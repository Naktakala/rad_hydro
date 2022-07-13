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
}


#endif //CHI_MATH_BANDED_SOLVERS_H
