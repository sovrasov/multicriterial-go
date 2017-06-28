#include "multicriterialSolver.hpp"
#include <iostream>

int main(int argc, const char** argv)
{
  MCOProblem problem;
  problem.AddCriterion(
    [](const double* y) -> double { return y[0]*y[0];  }
  );
  problem.AddCriterion(
    [](const double* y) -> double { return (y[0] - 2)*(y[0] - 2);  }
  );
  problem.SetDomain(1, {-100.}, {100.});

  MCOSolver solver;
  solver.SetParameters(SolverParameters(0.0001, 3, 1, 10000));
  solver.SetProblem(problem);
  solver.Solve();
  auto solution = solver.GetWeakOptimalPoints();

  std::cout << "Iterations performed: " << solver.GetIterationsNumber() << std::endl;

  for(auto x : solution)
  {
    std::cout << x.h << "\n";
  }

  return 0;
}
