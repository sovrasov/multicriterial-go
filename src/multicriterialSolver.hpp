#pragma once

#include "dataTypes.hpp"
#include "mcoProblem.hpp"
#include "evolvent.hpp"
#include "queue.hpp"

struct SolverParameters
{
  double eps;
  double r;
  unsigned numThreads;
  unsigned trialsLimit;
  unsigned evloventTightness = 12;
  bool verbose = false;

  SolverParameters() {}
  SolverParameters(double _eps, double _r,
      unsigned _numThreads, unsigned _trialsLimit) :
        eps(_eps), r(_r), numThreads(_numThreads), trialsLimit(_trialsLimit)
  {}
};

class MCOSolver
{
protected:

  Evolvent mEvolvent;
  std::vector<Trial> mSearchData;
  std::vector<Trial> mNextPoints;
  std::vector<double> mHEstimations;
  IntervalsQueue mNextIntervals;
  SolverParameters mParameters;
  MCOProblem mProblem;
  bool mNeedFullRecalc;
  unsigned mIterationsCounter;
  unsigned mNumberOfTrials;

  void InitDataStructures();
  void FirstIteration();
  void UpdateH(const Trial& left, const Trial& right);
  double ComputeH(const Trial& x1, const Trial& x2);
  void RecalcZ();
  void InsertNextPoints();
  void ClearDataStructures();
  void CalculateNextPoints();
  bool CheckStopCondition();

public:

  void SetParameters(const SolverParameters& params);
  void SetProblem(const MCOProblem& problem);

  void Solve();

  std::vector<Trial> GetWeakOptimalPoints();
  int GetIterationsNumber() const;
};
