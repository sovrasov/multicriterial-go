#pragma once

#include "dataTypes.hpp"
#include "mcoProblem.hpp"

class MCOSolver
{
protected:
public:
  void SetParameters();
  void SetProblem(const MCOProblem& problem);

  void Solve();

  std::vector<Trial> GetWeakOptimalPoints();
  int GetIterationsNumber() const;
};
