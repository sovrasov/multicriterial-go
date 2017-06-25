#include "multicriterialSolver.hpp"

#include <algorithm>
#include <cmath>

const double zeroHLevel = 1e-12;

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
    double y[solverMaxDim];
    mEvolvent.GetImage(mSearchData.back().x, y);
    for(int j = 0; j < mProblem.GetCriterionsNumber(); j++)
      mSearchData.back().z[j] = mProblem.CalculateFunction(j, y);

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

}

void MCOSolver::InsertNextPoints()
{

}

void MCOSolver::CalculateNextPoints()
{

}

void MCOSolver::ClearDataStructures()
{

}

bool MCOSolver::CheckStopCondition()
{
  return true;
}

std::vector<Trial> MCOSolver::GetWeakOptimalPoints()
{
  return std::vector<Trial>();
}
int MCOSolver::GetIterationsNumber() const
{
  return 0;
}
