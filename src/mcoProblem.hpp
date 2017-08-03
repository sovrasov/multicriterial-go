#pragma once

#include <vector>
#include <functional>
#include <initializer_list>

class MCOProblem
{
protected:
  using functionType = std::function<double(const double*)>;

  int mDimension;
  int mCriterionsNumber;
  std::vector<functionType> mCriterions;
  std::vector<functionType> mConstraints;
  std::vector<double> mLeftBound;
  std::vector<double> mRightBound;

public:

  void AddCriterion(const functionType& criterion);
  void AddConstraint(const functionType& criterion);
  void SetDomain(int dimension, const double* lb, const double* ub);
  void SetDomain(int dimension, std::initializer_list<double> lb,
      std::initializer_list<double> ub);

  int GetDimension() const;
  int GetCriterionsNumber() const;
  int GetConstraintsNumber() const;
  void GetBounds(double* lb, double* ub) const;
  double CalculateCriterion(int fNumber, const double* y) const;
  double CalculateConstraint(int fNumber, const double* y) const;
};
