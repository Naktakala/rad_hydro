#ifndef chi_math_MatrixNXxNX_h
#define chi_math_MatrixNXxNX_h

#include <iostream>

#include <vector>
#include <cmath>
#include <sstream>
#include <array>
#include <type_traits>
#include "ChiMath/chi_math_vectorNX.h"

#include "chi_log.h"

namespace chi_math{
  //#################################################################
  /**Generalized fixed size square matrix using chi_math::VectorNX for rows.
   * \author Jerry, Jan.*/

  template <int N, class NumberFormat>
  struct MatrixNXxNX
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for VectorNX." );
    std::array<chi_math::VectorNX<N,NumberFormat>,N> elements;

    /**Constructor for a matrix filled with 0s.*/
    MatrixNXxNX()
    {
      chi_math::VectorNX<N,NumberFormat> newVector;
      for(int i = 0; i<N;++i)
        elements[i] = newVector;
    }

    /**Constructor for a matrix filled with value.*/
    explicit MatrixNXxNX(NumberFormat value)
    {
      chi_math::VectorNX<N,NumberFormat> newVector(value);
      for(int i = 0; i<N;++i)
        elements[i] = newVector;
    }

    /**Constructor for a matrix for std::vector of std::vector.*/
    explicit MatrixNXxNX(const std::vector<std::vector<NumberFormat>>& A)
    {
      bool checks_passed = true;
      if (A.size() != N) checks_passed = false;
      for (const auto& row : A)
        if (row.size() != N) checks_passed = false;

      if (not checks_passed)
        throw std::logic_error("MatrixNXxNX constructor from vecofvec "
                               "encountered incompatible dimensions.");

      chi_math::VectorNX<N,NumberFormat> newVector;
      for(int i = 0; i<N;++i)
        for (int j=0; j<N; ++j)
        elements[i][j] = A[i][j];
    }

    /**Component-wise copy.*/
    MatrixNXxNX& operator=(const MatrixNXxNX& rhs)
    {
      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          elements[i](j) = rhs.elements[i][j];

      return *this;
    }

    /**Component-wise assignment from initializer list.*/
    MatrixNXxNX& operator=(const std::initializer_list<std::initializer_list<NumberFormat>>& lst)
    {
      if (lst.size() != N)
        throw std::logic_error("MatrixNXxNX constructor from initializer list "
                               "encountered incompatible dimensions.");

      int i = 0;
      for (auto &row: lst)
      {
        if (row.size() != N)
          throw std::logic_error("MatrixNXxNX constructor from initializer list "
                                 "encountered incompatible dimensions.");

        int j = 0;
        for (NumberFormat value: row)
          elements[i](j++) = value;

        ++i;
      }

      return *this;
    }

    //============================================ Addition
    /** Component-wise addition by a scalar.*/
    MatrixNXxNX& Shift( const NumberFormat value )
    {
      for(int i = 0;i<N;++i)
        for(int j = 0;  j<N; ++j)
          elements[i](j) += value;

      return *this;
    }

    /** Component-wise addition of two matrices.*/
    MatrixNXxNX operator+(const MatrixNXxNX& rhs) const
    {
      MatrixNXxNX<N,NumberFormat> newMatrix;
      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          newMatrix.elements[i](j) = elements[i][j] + rhs.elements[i][j];

      return newMatrix;
    }

    /** In-place component-wise addition of two matrices.*/
    MatrixNXxNX& operator +=(const MatrixNXxNX& rhs)
    {
      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          elements[i](j) += rhs.elements[i][j];

      return *this;
    }
    //============================================ Subtraction
    /** Component-wise subtraction of two matrices.*/
    MatrixNXxNX operator-(const MatrixNXxNX& rhs)
    {
      MatrixNXxNX<N,NumberFormat> newMatrix;
      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          newMatrix.elements[i](j) = elements[i][j]- rhs.elements[i][j];


      return newMatrix;
    }
    /** In-place componet wise subtraction of two matrices.*/
    MatrixNXxNX& operator -=(const MatrixNXxNX& rhs)
    {
      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          elements[i](j) -= rhs.elements[i][j];

      return *this;
    }
    //============================================ Multiplication
    /** Component-wise multiplication by a scalar.*/
    MatrixNXxNX operator*(NumberFormat value) const
    {
      MatrixNXxNX<N,NumberFormat> newMatrix;

      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          newMatrix.elements[i](j) = elements[i][j] * value;

      return newMatrix;
    }


    /** Pure multiplication of two matricies.*/
    MatrixNXxNX operator*(const MatrixNXxNX& rhs) const
    {
      MatrixNXxNX<N,NumberFormat> result;
      for(int i = 0; i < elements.size(); ++i)
        for(int j = 0; j < N; ++j)
          for(int k = 0; k < N; ++k)
            result.elements[i](j) += elements[i][k] * rhs.elements[k][j];

      return result;
    }

    /** Multiplication of matrix and vector.*/
    chi_math::VectorNX<N,NumberFormat>
        operator*(const chi_math::VectorNX<N,NumberFormat>& rhs) const
    {
      chi_math::VectorNX<N,NumberFormat> result;
      for (int i=0; i<N; ++i)
        for(int j=0; j<N; ++j)
          result(i) += elements[i][j]*rhs.elements[j];

      return result;
    }

    /** Component-wise multiplication of two matrices.*/
    MatrixNXxNX MultComponentWise(const MatrixNXxNX& rhs) const
    {
      MatrixNXxNX<N,NumberFormat> newMatrix;

      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          newMatrix.elements[i](j) = elements[i][j] * rhs.elements[i][j];

      return newMatrix;
    }

    //=============================================Division
    /** Component-wise division by a scalar.*/
    MatrixNXxNX operator/(NumberFormat value) const
    {
      MatrixNXxNX<N,NumberFormat> newMatrix;

      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          newMatrix.elements[i](j) = elements[i][j] / value;

      return newMatrix;
    }

    /** Component-wise division of two matrices.*/
    MatrixNXxNX DivComponentWise(const MatrixNXxNX& rhs) const
    {
      MatrixNXxNX<N,NumberFormat> newMatrix;

      for(int i = 0;i<N;++i)
        for(int j = 0;j<N;++j)
          newMatrix.elements[i](j) = elements[i][j] / rhs.elements[i][j];

      return newMatrix;
    }

    //=============================================Element access
    /**Returns the value stored at (i,j).*/
    NumberFormat GetIJ(int i, int j)
    {
      return elements[i][j];
    }

    /**Sets the position (i,j) equal to value.*/
    void SetIJ(int i, int j, NumberFormat value)
    {
      elements[i](j) = value;
    }

    /**Adds value to the exisiting value at (i,j).*/
    void AddIJ(int i, int j, NumberFormat value)
    {
      elements[i](j) += value;
    }

    //=============================================Utilities
    /**Sets the diagonal of the matrix.*/
    void SetDiagonalVec(std::vector<NumberFormat> vec)
    {
      if (vec.size()!= N){
        auto& log = chi_objects::ChiLog::GetInstance();
        log.LogAllError()<<"MatrixNXxNX::SetDiagonalVec, incompatible vector size";
        exit(EXIT_FAILURE);
      }
      for(int i = 0;i<N;++i)
        elements[i](i) = vec[i];
    }
    /**Set row i using a vector.*/
    void SetRowIVec(int i, std::vector<NumberFormat> vec){
      if (N!=vec.size()){
        auto& log = chi_objects::ChiLog::GetInstance();
        log.LogAllError()<<"error"<<std::endl;
      }

      for (int j =0; j<N;++j)
        elements[i](j) = vec[j];

    }

    /** Sets column j using a vector.*/
    void SetColJVec(int j, std::vector<NumberFormat> vec){
      if (N!=vec.size()){
        auto& log = chi_objects::ChiLog::GetInstance();
        log.LogAllError()<<"error"<<std::endl;
      }

      for (int i =0; i<N;++i)
        elements[i](j) = vec[i];

    }

    /* Recursive function for finding determinant of matrix.
    n is current dimension of mat[][]. */
    NumberFormat Det() const
    {
      NumberFormat determinant = 0.0;

      if (N <= 0)
        return 0.0;
      else if (N == 1)
        return elements[0][0];
      else if (N == 2)
        return elements[0][0]*elements[1][1] - elements[0][1]*elements[1][0];
      else if (N>=3)
      {
        MatrixNXxNX<N,NumberFormat> newMatrix;
        newMatrix.elements = elements;

        int sign = -1;
        for (int i=0; i<N; i++)
        {
          for (int j=0; j<N; j++)
          {
            sign*=-1;
            newMatrix.SetIJ(i,j,newMatrix.GetIJ(i,j)*sign);
          }
        }
        for(int j = 0; j<N; ++j)
          determinant = determinant + newMatrix.elements[0][j]*MinorIJ(0,j);
      }

      return determinant;

    }

    /**Gets the minor value associated with row ir and column jr
    this is being used Recursively to find the determiant.*/
    NumberFormat MinorIJ(int ir, int jr) const
    {
      MatrixNXxNX<N-1,NumberFormat> minor;
      int a = 0;
      int b = 0;
      for (int i = 0; i<N;++i){
        if (i == ir) continue;
        for (int j = 0; j<N; ++j){
          if (j == jr) continue;

          minor.elements[a](b) = elements[i][j];
          ++b;
        }
        b = 0;
        ++a;
      }
      return minor.Det();
    }

    /**Compute the matrix transpose.*/
    MatrixNXxNX Transpose() const
    {
      MatrixNXxNX<N,NumberFormat> newMatrix;
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
          newMatrix.elements[j](i) = elements[i][j];

      return newMatrix;
    }
    /**Computes the matrix inverse.*/
    MatrixNXxNX Inverse() const
    {
      MatrixNXxNX<N,NumberFormat> invertedMatrix;
      // ======================= compute the matrix of minors
      for (int i = 0; i<N; ++i)
        for (int j = 0; j<N; ++j)
          invertedMatrix.SetIJ(i,j,MinorIJ(i,j));

      // ======================= compute the matrix of cofactors
      int sign = -1;
      for (int i=0; i<N; i++)
      {
        for (int j=0; j<N; j++)
        {
          sign*=-1;
          invertedMatrix.SetIJ(i,j,invertedMatrix.GetIJ(i,j)*sign);
        }
      }

      // ======================= compute the transpose
      invertedMatrix = invertedMatrix.Transpose();

      //======================== compute the determinant
      double determiant = Det();

      return invertedMatrix*(1/determiant);
    }

    /** Computes the L2-norm of the matrix.
     * \sqrt{\sum a_ij^2} */
    NumberFormat Norm()
    {
      NumberFormat adderupper = 0.0;
      for (int i = 0; i<N; ++i)
        for (int j = 0; j<N;++j)
          adderupper += (elements[i][j]*elements[i][j]);
      return sqrt(adderupper);
    }

    void Print() const
    {
      for(int i = 0;i<N;++i){
        for(int j = 0;j<N-1;++j)
          std::cout<<elements[i][j]<<" ";

        std::cout<<elements[i][N-1]<<std::endl;
      }


    }

    /**Outputs matrix to stringstream.*/
    std::string PrintStr() const
    {
      std::stringstream out;

      for (int i=0; i<N; ++i){
        out << "[";
        for(int j=0; j<N; ++j){
          out<<elements[i][j];

          if(j!=N-1)
            out << " ";
        }

        if (i <= (N-1))
          out << "]\n";
        else
          out << "]" << std::endl;
      }


      return out.str();
    }



  };

  /** This prevents the det() from creating a matrix of size -1.*/
  template <class NumberFormat>
  struct MatrixNXxNX<0,NumberFormat>
  {
    std::array<VectorNX<1,NumberFormat>,1> elements;

    double Minor(int i, int j)
    {
      return 0.0;
    }

    double Det()
    {
      return 0.0;
    }
  };
}

/**Multiplication by a scalar from the left.*/
template<int N,class NumberFormat>
chi_math::MatrixNXxNX<N,NumberFormat>
operator*(const double value, const chi_math::MatrixNXxNX<N,NumberFormat>& rhs)
{
  chi_math::MatrixNXxNX<N,NumberFormat> newMatrix;
  for (int i = 0; i<N;++i)
    newMatrix.elements[i] = rhs.elements[i]*value;
  return newMatrix;
}
#endif