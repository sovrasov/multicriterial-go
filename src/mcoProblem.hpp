#pragma once

#include <vector>
#include <functional>
#include <initializer_list>

class MCOProblem
{
protected:
  int mDimension;
  int mCriterionsNumber;
  std::vector<double> mLeftBound;
  std::vector<double> mRightBound;
public:

  void AddCriterion(std::function<double(const double*)>);
  void SetDomain(int dimension, const double* lb, const double* ub);
  void SetDomain(int dimension, std::initializer_list<double> lb,
      std::initializer_list<double> ub);

  int GetDimansion() const;
  int GetCriterionsNumber() const;
  void GetBounds(double* lb, double* ub) const;
};
