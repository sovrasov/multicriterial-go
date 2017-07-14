#include "multicriterialSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  const double zeroHLevel = 1e-12;
  bool isVectorLess(const double* v1, const double* v2, int dim)
  {
    for (int i = 0; i < dim; i++)
      if (v1[i] >= v2[i])
        return false;

    return true;
  }
}

MCOSolver::MCOSolver() : mLocalOffset(pow(1.5, -15)) {}

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
  RecalcZandR();

  do {
    CalculateNextPoints();
    InsertNextPoints();
    RecalcZandR();
    needStop = CheckStopCondition();
    mIterationsCounter++;
  } while(mIterationsCounter < mParameters.iterationsLimit && !needStop);

  ClearDataStructures();
}

void MCOSolver::InitDataStructures()
{
  double leftDomainBound[solverMaxDim], rightDomainBound[solverMaxDim];
  mProblem.GetBounds(leftDomainBound, rightDomainBound);
  mEvolvent = Evolvent(mProblem.GetDimension(), mParameters.evloventTightness, leftDomainBound, rightDomainBound);
  mNextPoints.resize(mParameters.numThreads);
  mNextIntervals.resize(mParameters.numThreads);

  mHEstimations.resize(mProblem.GetCriterionsNumber());
  std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);
}

void MCOSolver::FirstIteration()
{
  for(size_t i = 0; i <= mParameters.numThreads; i++)
  {
    auto insertedIt = mSearchData.emplace((double)i / mParameters.numThreads).first;
    mEvolvent.GetImage(insertedIt->x, insertedIt->y);
    for(int j = 0; j < mProblem.GetCriterionsNumber(); j++)
      insertedIt->z[j] = mProblem.CalculateFunction(j, insertedIt->y);

    if(i > 0)
    {
      auto rightPoint = insertedIt--;
      UpdateH(*insertedIt, *rightPoint);
    }
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

void MCOSolver::RecalcZandR()
{
  mNextIntervals.clear();
  bool isLocal = IsLocalIteration();
  for(auto trial = mSearchData.begin(); trial != mSearchData.end(); ++trial)
  {
    if(!mNeedFullRecalc)
    {
      for (size_t j = 0; j < mNextPoints.size(); j++)
        trial->h = fmax(trial->h, ComputeH(*trial, mNextPoints[j]));
    }
    else
    {
      trial->h = std::numeric_limits<double>::min();
      for(auto otherTrial = mSearchData.begin(); otherTrial != mSearchData.end(); ++otherTrial)
        trial->h = fmax(trial->h, ComputeH(*trial, *otherTrial));
    }

    if(trial != mSearchData.begin() && !isLocal)
    {
      auto leftTrial = trial;
      auto rightTrial = leftTrial--;
      Interval currentInt(*leftTrial, *rightTrial);
      currentInt.delta = pow(currentInt.pr.x - currentInt.pl.x, 1. / mProblem.GetDimension());
      currentInt.R = CalculateR(currentInt);
      mNextIntervals.pushWithPriority(currentInt);
    }
  }

  if(isLocal)
  {
    double zOpt = std::numeric_limits<double>::max();
    for(const auto& point : mSearchData)
      zOpt = fmin(zOpt, point.h);

    auto trial = mSearchData.begin();
    for(++trial; trial != mSearchData.end(); ++trial)
    {
      auto leftTrial = trial;
      auto rightTrial = leftTrial--;
      Interval currentInt(*leftTrial, *rightTrial);
      currentInt.delta = pow(currentInt.pr.x - currentInt.pl.x, 1. / mProblem.GetDimension());
      currentInt.R = CalculateLocalR(currentInt, zOpt);
      mNextIntervals.pushWithPriority(currentInt);
    }
  }

  mNeedFullRecalc = false;
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
    mNextPoints[i].h = ComputeH(mNextPoints[i], mNextPoints[i]);
    for(const auto& point : mSearchData)
      mNextPoints[i].h = fmax(mNextPoints[i].h, ComputeH(mNextPoints[i], point));
    mSearchData.insert(mNextPoints[i]);
  }
}

void MCOSolver::CalculateNextPoints()
{
#pragma omp parallel for num_threads(mParameters.numThreads)
  for(int i = 0; i < (int)mParameters.numThreads; i++)
  {
    Interval nextInterval;
#pragma omp critical
    {
      nextInterval = mNextIntervals.pop();
    }

    double dh = nextInterval.pr.h - nextInterval.pl.h;
    mNextPoints[i].x = 0.5 * (nextInterval.pr.x + nextInterval.pl.x) -
      0.5*((dh > 0.) ? 1. : -1.) * pow(fabs(dh), mProblem.GetDimension()) / mParameters.r;

    if (mNextPoints[i].x >= nextInterval.pr.x || mNextPoints[i].x <= nextInterval.pl.x)
      throw std::runtime_error("Point is outside of interval\n");

    mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
    for(int j = 0; j < mProblem.GetCriterionsNumber(); j++)
      mNextPoints[i].z[j] = mProblem.CalculateFunction(j, mNextPoints[i].y);

#pragma omp critical
    {
      UpdateH(mNextPoints[i], nextInterval.pr);
      UpdateH(nextInterval.pl, mNextPoints[i]);
    }
  }
}

void MCOSolver::ClearDataStructures()
{
  mNextIntervals.clear();
}

bool MCOSolver::CheckStopCondition() const
{
  auto nextIntervals = mNextIntervals.getElements();
  for(size_t i = 0; i < nextIntervals.size(); i++)
    if(nextIntervals[i].delta < mParameters.eps)
      return true;

  return false;
}

std::vector<Trial> MCOSolver::GetWeakOptimalPoints() const
{
  std::vector<Trial> optTrials;

  for(const auto& point1 : mSearchData)
  {
    bool isWeakOptimal = true;
    for(const auto& point2 : mSearchData)
    {
      if(point1.x != point2.x)
      {
        if(isVectorLess(point2.z, point1.z, mProblem.GetCriterionsNumber()))
          isWeakOptimal = false;
      }
    }
    if(isWeakOptimal)
      optTrials.push_back(point1);
  }

  return optTrials;
}

int MCOSolver::GetIterationsNumber() const
{
  return mIterationsCounter;
}

bool MCOSolver::IsLocalIteration() const
{
  if(mIterationsCounter < 100)
    return false;

  bool isLocal = false;

  if (mParameters.localMix > 0) {
    int localMixParameter = mParameters.localMix + 1;

    if (mIterationsCounter % localMixParameter != 0)
      isLocal = false;
    else
      isLocal = true;
  }
  else if (mParameters.localMix < 0) {
    int localMixParameter = -mParameters.localMix;
    localMixParameter++;

    if (mIterationsCounter % localMixParameter != 0)
      isLocal = true;
    else
      isLocal = false;
  }
  else //mParameters.localMix == 0
    isLocal = false;

  return isLocal;
}

double MCOSolver::CalculateR(const Interval& i) const
{
  return i.delta + pow((i.pr.h - i.pl.h) / mParameters.r, 2) / i.delta -
      2.*(i.pr.h + i.pl.h) / mParameters.r;
}

double MCOSolver::CalculateLocalR(const Interval& i, double minH) const
{
  return CalculateR(i) / (sqrt((i.pr.h - minH)*(i.pl.h - minH)) + mLocalOffset);
}
