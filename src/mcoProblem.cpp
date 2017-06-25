#include "mcoProblem.hpp"

void MCOProblem::AddCriterion(const functionType& criterion)
{
  mCriterions.push_back(criterion);
}

void MCOProblem::SetDomain(int dimension, const double* lb, const double* ub)
{
  mDimension = dimension;
  mLeftBound.resize(mDimension);
  mRightBound.resize(mDimension);

  for(int i = 0; i < mDimension; i++)
  {
    mLeftBound[i] = lb[i];
    mRightBound[i] = ub[i];
  }
}
void MCOProblem::SetDomain(int dimension, std::initializer_list<double> lb,
    std::initializer_list<double> ub)
{
  mDimension = dimension;
  mLeftBound.reserve(mDimension);
  mRightBound.reserve(mDimension);

  for(auto it = lb.begin(); it != lb.end(); ++it)
    mLeftBound.push_back(*it);

  for(auto it = ub.begin(); it != ub.end(); ++it)
    mRightBound.push_back(*it);
}

int MCOProblem::GetDimension() const
{
  return mDimension;
}

int MCOProblem::GetCriterionsNumber() const
{
  return mCriterions.size();
}

void MCOProblem::GetBounds(double* lb, double* ub) const
{
  for(int i = 0; i < mDimension; i++)
  {
    lb[i] = mLeftBound[i];
    ub[i] = mRightBound[i];
  }
}

double MCOProblem::CalculateFunction(int fNumber, const double* y) const
{
  return mCriterions[fNumber](y);
}
