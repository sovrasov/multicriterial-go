#include "multicriterialSolver.hpp"
#include <iostream>
#include <cmdline.h>

int main(int argc, const char** argv)
{
  cmdline::parser parser;
  parser.add<int>("threadsNum", 'p', "number of threads in parallel method", false, 1,
     cmdline::range(1, 32));
  parser.add<int>("evolventTightness", 't', "", false, 12,
        cmdline::range(9, 16));
  parser.add<double>("reliability", 'r', "reliability parameter for the method",
    false, 4.5, cmdline::range(1., 1000.));
  parser.add<double>("accuracy", 'e', "accuracy of the method", false, 0.01);
  parser.add<int>("itersLimit", 'l', "limit of iterations for the method", false, 2000);
  parser.parse_check(argc, argv);

  MCOProblem problem;
  problem.AddCriterion(
    [](const double* y) -> double { return y[0]*y[0];  }
  );
  problem.AddCriterion(
    [](const double* y) -> double { return (y[0] - 2)*(y[0] - 2);  }
  );
  problem.SetDomain(1, {-10.}, {10.});

  MCOSolver solver;
  solver.SetParameters(SolverParameters(parser.get<double>("accuracy"),
    parser.get<double>("reliability"),
    parser.get<int>("threadsNum"),
    parser.get<int>("itersLimit")));
  solver.SetProblem(problem);
  solver.Solve();
  auto solution = solver.GetWeakOptimalPoints();

  std::cout << "Iterations performed: " << solver.GetIterationsNumber() << std::endl;

  for(auto& x : solution)
  {
    std::cout << x.y[0] << "\n";
  }

  return 0;
}
