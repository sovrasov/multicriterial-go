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

  template<typename T>
  size_t insert_sorted(std::vector<T> & vec, T const& item)
  {
    return std::distance(vec.begin(), vec.insert
    (std::upper_bound(vec.begin(), vec.end(), item), item));
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
  mSearchData.reserve(mParameters.iterationsLimit);
  mNextPoints.resize(mParameters.numThreads);
  mNextIntervals.resize(mParameters.numThreads);

  mZEstimations.resize(mProblem.GetConstraintsNumber() + 1);
  std::fill(mZEstimations.begin(), mZEstimations.end(), std::numeric_limits<double>::max());

  mHEstimations.resize(mProblem.GetCriterionsNumber());
  std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);

  mMuEstimations.resize(mProblem.GetConstraintsNumber());
  std::fill(mMuEstimations.begin(), mMuEstimations.end(), 1.0);
}

void MCOSolver::FirstIteration()
{
  for(size_t i = 0; i <= mParameters.numThreads; i++)
  {
    mSearchData.emplace_back((double)i / mParameters.numThreads);
    MakeTrial(mSearchData.back());
    if(i > 0)
      CalculateConstsEstimationsAfterInsert(i, false);
      //UpdateH(mSearchData[i - 1], mSearchData[i]);
  }

  mNeedFullRecalc = true;
  mIterationsCounter = 1;
  mNumberOfTrials = mParameters.numThreads + 1;
}

void MCOSolver::MakeTrial(Trial& trial)
{
  mEvolvent.GetImage(trial.x, trial.y);

  trial.v = 0;
  while(trial.v < mProblem.GetConstraintsNumber())
  {
    trial.g[trial.v] = mProblem.CalculateConstraint(trial.v, trial.y);
    if(trial.g[trial.v] > 0)
      break;
    else
      trial.v++;
  }

  if(trial.v < mProblem.GetConstraintsNumber())
  {
    trial.h = trial.g[trial.v] / mMuEstimations[trial.v];
    mZEstimations[trial.v] = fmin(mZEstimations[trial.v], trial.g[trial.v]);
  }

  if(trial.v == mProblem.GetConstraintsNumber())
    for(int j = 0; j < mProblem.GetCriterionsNumber(); j++)
      trial.z[j] = mProblem.CalculateCriterion(j, trial.y);
}

void MCOSolver::CalculateConstsEstimationsAfterInsert(size_t idx, bool searchRight)
{
  Trial& currentPoint = mSearchData[idx];
  if(currentPoint.v < mProblem.GetConstraintsNumber()) //update Mu for constraints
  {
    int left_idx = idx - 1;
    while(left_idx > 0 && mSearchData[left_idx].v != currentPoint.v)
      left_idx--;
    if(left_idx != (int)idx && mSearchData[left_idx].v == mSearchData[idx].v)
      UpdateZ(mSearchData[left_idx], mSearchData[idx]);

    if(searchRight)
    {
      size_t right_idx = idx + 1;
      while(right_idx < mSearchData.size() - 1 && mSearchData[right_idx].v != currentPoint.v)
        right_idx++;
      if(right_idx != idx && mSearchData[right_idx].v == mSearchData[idx].v)
        UpdateZ(mSearchData[idx], mSearchData[right_idx]);
    }
  }
  else //update H for criterions
  {
    int left_idx = idx - 1;
    while(left_idx > 0 && mSearchData[left_idx].v != currentPoint.v)
      left_idx--;
    if(left_idx != (int)idx && mSearchData[left_idx].v == mSearchData[idx].v)
      UpdateH(mSearchData[left_idx], mSearchData[idx]);

    if(searchRight)
    {
      size_t right_idx = idx + 1;
      while(right_idx < mSearchData.size() - 1 && mSearchData[right_idx].v != currentPoint.v)
        right_idx++;
      if(right_idx != idx && mSearchData[right_idx].v == mSearchData[idx].v)
        UpdateH(mSearchData[idx], mSearchData[right_idx]);
    }
  }
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

void MCOSolver::UpdateZ(const Trial& left, const Trial& right)
{
  //assert(left.v == right.v);
  double oldZ = mZEstimations[left.v];
  double newZ = fabs(right.g[right.v] - left.g[left.v]) /
    pow(right.x - left.x, 1. / mProblem.GetDimension());
  if (newZ > oldZ || (oldZ == 1.0 && newZ > zeroHLevel))
  {
    mZEstimations[left.v] = newZ;
    //mNeedFullRecalc = true;
  }
}

void MCOSolver::RecalcZandR()
{
  mNextIntervals.clear();
  bool isLocal = IsLocalIteration();
  const int constraintsNumber = mProblem.GetConstraintsNumber();
  for(size_t i = 0; i < mSearchData.size(); i++)
  {
    if(mSearchData[i].v == constraintsNumber)
    {
      if(!mNeedFullRecalc)
      {
        for (size_t j = 0; j < mNextPoints.size(); j++)
          if(mNextPoints[j].v == constraintsNumber)
            mSearchData[i].h = fmax(mSearchData[i].h, ComputeH(mSearchData[i], mNextPoints[j]));
      }
      else
      {
        mSearchData[i].h = std::numeric_limits<double>::min();
        for(size_t j = 0; j < mSearchData.size(); j++)
          if(mSearchData[j].v == constraintsNumber)
            mSearchData[i].h = fmax(mSearchData[i].h, ComputeH(mSearchData[i], mSearchData[j]));
      }
      mZEstimations[constraintsNumber] = fmin(mZEstimations[constraintsNumber], mSearchData[i].h);
    }

    if(i > 0 && !isLocal)
    {
      Interval currentInt(mSearchData[i - 1], mSearchData[i]);
      currentInt.delta = pow(currentInt.pr.x - currentInt.pl.x, 1. / mProblem.GetDimension());
      currentInt.R = CalculateR(currentInt);
      mNextIntervals.pushWithPriority(currentInt);
    }
  }

  if(isLocal)
  {
    for(size_t i = 1; i < mSearchData.size(); i++)
    {
      Interval currentInt(mSearchData[i - 1], mSearchData[i]);
      currentInt.delta = pow(currentInt.pr.x - currentInt.pl.x, 1. / mProblem.GetDimension());
      currentInt.R = CalculateLocalR(currentInt);
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
  const int constraintsNumber = mProblem.GetConstraintsNumber();
  for(size_t i = 0; i < mNextPoints.size(); i++)
  {
    if(mNextPoints[i].v == constraintsNumber)
    {
      mNextPoints[i].h = ComputeH(mNextPoints[i], mNextPoints[i]);
      for (size_t j = 0; j < mSearchData.size(); j++)
        if(mSearchData[j].v == constraintsNumber)
          mNextPoints[i].h = fmax(mNextPoints[i].h, ComputeH(mNextPoints[i], mSearchData[j]));
    }
    size_t insert_idx = insert_sorted(mSearchData, mNextPoints[i]);
    CalculateConstsEstimationsAfterInsert(insert_idx);
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

    if(nextInterval.pr.v == nextInterval.pl.v)
    {
      double dh = nextInterval.pr.h - nextInterval.pl.h;
      mNextPoints[i].x = 0.5 * (nextInterval.pr.x + nextInterval.pl.x) -
        0.5*((dh > 0.) ? 1. : -1.) * pow(fabs(dh), mProblem.GetDimension()) / mParameters.r;
    }
    else
      mNextPoints[i].x = 0.5 * (nextInterval.pr.x + nextInterval.pl.x);

    if (mNextPoints[i].x >= nextInterval.pr.x || mNextPoints[i].x <= nextInterval.pl.x)
      throw std::runtime_error("Point is outside of interval\n");

    MakeTrial(mNextPoints[i]);
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

  for(size_t i = 0; i < mSearchData.size(); i++)
  {
    if(mSearchData[i].v == mProblem.GetConstraintsNumber())
    {
      bool isWeakOptimal = true;
      for(size_t j = 0; j < mSearchData.size(); j++)
      {
        if(i != j && mSearchData[j].v == mProblem.GetConstraintsNumber())
        {
          if(isVectorLess(mSearchData[j].z, mSearchData[i].z, mProblem.GetCriterionsNumber()))
            isWeakOptimal = false;
        }
      }
      if(isWeakOptimal)
        optTrials.push_back(mSearchData[i]);
    }
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
  if(i.pl.v == i.pr.v)
    return i.delta + pow((i.pr.h - i.pl.h) / mParameters.r, 2) / i.delta -
      2.*(i.pr.h + i.pl.h - 2*mZEstimations[i.pr.v]) / mParameters.r;
  else if(i.pl.v < i.pr.v)
    return 2*i.delta - 4*(i.pr.h - mZEstimations[i.pr.v]);
  else
    return 2*i.delta - 4*(i.pl.h - mZEstimations[i.pl.v]);
}

double MCOSolver::CalculateLocalR(const Interval& i) const
{
  double value;
  if(i.pl.v == i.pr.v)
    value = CalculateR(i) / (sqrt((i.pr.h - mZEstimations[i.pr.v])*
        (i.pl.h - mZEstimations[i.pr.v])) + mLocalOffset);
  else if(i.pl.v < i.pr.v)
    value = CalculateR(i) / (i.pr.h - mZEstimations[i.pr.v] + mLocalOffset);
  else
    value = CalculateR(i) / (i.pl.h - mZEstimations[i.pl.v] + mLocalOffset);

#ifdef WIN32
  if (!_finite(value))
#else
  if (!std::isfinite(value))
#endif
  {
    throw std::runtime_error("Infinite R!");
  }
  return value;
}
