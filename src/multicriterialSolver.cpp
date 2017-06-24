#include "multicriterialSolver.hpp"

void MCOSolver::SetParameters(const SolverParameters& params)
{
  mParameters = params;
}
void MCOSolver::SetProblem(const MCOProblem& problem)
{

}

void MCOSolver::Solve()
{

}

std::vector<Trial> MCOSolver::GetWeakOptimalPoints()
{
  return std::vector<Trial>();
}
int MCOSolver::GetIterationsNumber() const
{
  return 0;
}
