#include "multicriterialSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

const double zeroHLevel = 1e-12;

template< typename T >
typename std::vector<T>::iterator
   insert_sorted( std::vector<T> & vec, T const& item )
{
    return vec.insert
        (
            std::upper_bound( vec.begin(), vec.end(), item ),
            item
        );
}

bool isVectorLess(const double* v1, const double* v2, int dim)
{
  for(int i = 0; i < dim; i++)
    if(v1[i] >= v2[i])
      return false;

  return true;
}

void MCOSolver::SetParameters(const SolverParameters& params)
{
  mParameters = params;
}

void MCOSolver::SetProblem(const MCOProblem& problem)
{
  mProblem = problem;
}

void MCOSolver::Solve()
{
  bool needStop = false;
  InitDataStructures();
  FirstIteration();
  RecalcZ();

  do {
    CalculateNextPoints();
    InsertNextPoints();
    RecalcZ();
    needStop = CheckStopCondition();
    mIterationsCounter++;
  } while(mNumberOfTrials < mParameters.trialsLimit && !needStop);

  ClearDataStructures();
}

void MCOSolver::InitDataStructures()
{
  double leftDomainBound[solverMaxDim], rightDomainBound[solverMaxDim];
  mProblem.GetBounds(leftDomainBound, rightDomainBound);
  mEvolvent = Evolvent(mProblem.GetDimension(), mParameters.evloventTightness, leftDomainBound, rightDomainBound);
  mSearchData.reserve(mParameters.trialsLimit);
  mNextPoints.resize(mParameters.numThreads);
  mNextIntervals.resize(mParameters.numThreads);

  mHEstimations.resize(mProblem.GetCriterionsNumber());
  std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);
}

void MCOSolver::FirstIteration()
{
  for(size_t i = 0; i <= mParameters.numThreads; i++)
  {
    mSearchData.emplace_back((double)i / mParameters.numThreads);
    mEvolvent.GetImage(mSearchData.back().x, mSearchData.back().y);
    for(int j = 0; j < mProblem.GetCriterionsNumber(); j++)
      mSearchData.back().z[j] = mProblem.CalculateFunction(j, mSearchData.back().y);

    if(i > 0)
      UpdateH(mSearchData[i - 1], mSearchData[i]);
  }

  mNeedFullRecalc = true;
  mIterationsCounter = 1;
  mNumberOfTrials = mParameters.numThreads + 1;
}

void MCOSolver::UpdateH(const Trial& left, const Trial& right)
{
  for(int j = 0; j < mProblem.GetCriterionsNumber(); j++)
  {
    double oldH = mHEstimations[j];
    double newH = fabs(right.z[j] - left.z[j]) /
      pow(right.x - left.x, 1. / mProblem.GetDimension());
    if (newH > oldH || (oldH == 1.0 && newH > zeroHLevel))
    {
      mHEstimations[j] = newH;
      mNeedFullRecalc = true;
    }
  }
}

void MCOSolver::RecalcZ()
{
  for(size_t i = 0; i < mSearchData.size(); i++)
  {
    if(!mNeedFullRecalc)
    {
      for(size_t j = 0; j < mNextPoints.size(); j++)
        mSearchData[i].h = fmax(mSearchData[i].h, ComputeH(mSearchData[i], mNextPoints[j]));
    }
    else
    {
      mSearchData[i].h = std::numeric_limits<double>::min();
      for(size_t j = 0; j < mSearchData.size(); j++)
        mSearchData[i].h = fmax(mSearchData[i].h, ComputeH(mSearchData[i], mSearchData[j]));
    }

    mNextIntervals.clear();
    if(i > 0)
    {
      Interval currentInt(mSearchData[i-1], mSearchData[i]);
      currentInt.delta = pow(currentInt.pr.x - currentInt.pl.x, 1. / mProblem.GetDimension());
      currentInt.R = currentInt.delta + pow(currentInt.pr.h - currentInt.pl.h, 2) / currentInt.delta -
          2.*(currentInt.pr.h + currentInt.pl.h)/mParameters.r;
      mNextIntervals.pushWithPriority(currentInt);
    }
  }
}

double MCOSolver::ComputeH(const Trial& x1, const Trial& x2)
{
  double value = std::numeric_limits<double>::max();

  for(int i = 0; i < mProblem.GetCriterionsNumber(); i++)
  {
    value = fmin(value, (x1.z[i] - x2.z[i]) / mHEstimations[i]);
  }

  return value;
}

void MCOSolver::InsertNextPoints()
{
  for(size_t i = 0; i < mNextPoints.size(); i++)
  {
    insert_sorted(mSearchData, mNextPoints[i]);
  }
}

void MCOSolver::CalculateNextPoints()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    Interval nextInterval = mNextIntervals.pop();
    mNextPoints[i] = 0.5 * (nextInterval.pr.x + nextInterval.pl.x) -
    0.5*(((nextInterval.pr.h - nextInterval.pl.h) > 0.) ? 1. : -1.) *
      pow(fabs(nextInterval.pr.h - nextInterval.pl.h), mProblem.GetDimension()) / mParameters.r;

    mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
      for(int j = 0; j < mProblem.GetCriterionsNumber(); j++)
        mNextPoints[i].z[j] = mProblem.CalculateFunction(j, mNextPoints[i].y);
    UpdateH(nextInterval.pl, nextInterval.pr);
  }
}

void MCOSolver::ClearDataStructures()
{
  mNextIntervals.clear();
}

bool MCOSolver::CheckStopCondition()
{
  auto nextIntervals = mNextIntervals.getElements();
  for(size_t i = 0; i < nextIntervals.size(); i++)
    if(nextIntervals[i].delta < mParameters.eps)
      return true;

  return false;
}

std::vector<Trial> MCOSolver::GetWeakOptimalPoints()
{
  std::vector<Trial> optTrials;

  for(size_t i = 0; i < mSearchData.size(); i++)
  {
    bool isWeakOptimal = true;
    for(size_t j = 0; j < mSearchData.size(); j++)
    {
      if(i != j)
      {
        if(isVectorLess(mSearchData[j].z, mSearchData[i].z, mProblem.GetCriterionsNumber()))
          isWeakOptimal = false;
      }
    }
    if(isWeakOptimal)
      optTrials.push_back(mSearchData[i]);
  }

  return mSearchData;
}
int MCOSolver::GetIterationsNumber() const
{
  return mIterationsCounter;
}
