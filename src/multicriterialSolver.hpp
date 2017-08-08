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
  unsigned iterationsLimit;
  unsigned evloventTightness = 12;
  int localMix;
  bool verbose = false;

  SolverParameters() {}
  SolverParameters(double _eps, double _r,
      unsigned _numThreads, unsigned _itersLimit, int _localMix = 0) :
        eps(_eps), r(_r), numThreads(_numThreads), iterationsLimit(_itersLimit),
        localMix(_localMix)
  {}
};

class MCOSolver
{
protected:

  Evolvent mEvolvent;
  std::vector<Trial> mSearchData;
  std::vector<Trial> mNextPoints;
  std::vector<double> mHEstimations;
  std::vector<double> mMuEstimations;
  std::vector<double> mZEstimations;
  IntervalsQueue mNextIntervals;
  SolverParameters mParameters;
  MCOProblem mProblem;
  bool mNeedFullRecalc;
  unsigned mIterationsCounter;
  unsigned mNumberOfTrials;
  const double mLocalOffset;

  void InitDataStructures();
  void MakeTrial(Trial& trial);
  void FirstIteration();
  void UpdateH(const Trial& left, const Trial& right);
  void UpdateMu(const Trial& left, const Trial& right);
  void CalculateConstsEstimationsAfterInsert(size_t idx, bool searchRight = true);
  double ComputeH(const Trial& x1, const Trial& x2);
  void RecalcZandR();
  void InsertNextPoints();
  void ClearDataStructures();
  void CalculateNextPoints();
  bool CheckStopCondition() const;
  bool IsLocalIteration() const;
  double CalculateR(const Interval&) const;
  double CalculateLocalR(const Interval&) const;

public:
  MCOSolver();

  void SetParameters(const SolverParameters& params);
  void SetProblem(const MCOProblem& problem);

  void Solve();

  std::vector<Trial> GetWeakOptimalPoints() const;
  int GetIterationsNumber() const;
};
