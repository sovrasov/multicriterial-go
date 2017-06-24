#pragma once

#include "dataTypes.hpp"
#include "mcoProblem.hpp"
#include "evolvent.hpp"

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

  SolverParameters mParameters;

public:

  void SetParameters(const SolverParameters& params);
  void SetProblem(const MCOProblem& problem);

  void Solve();

  std::vector<Trial> GetWeakOptimalPoints();
  int GetIterationsNumber() const;
};
