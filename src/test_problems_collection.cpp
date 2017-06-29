#include "test_problems_collection.hpp"
#include <cmath>

MCOProblem TestMCOProblems::create(const std::string& name, int dimension)
{
  MCOProblem problem;
  if(name == "schaffer")
  {
    problem.AddCriterion(
      [](const double* y) -> double { return y[0]*y[0];  }
    );
    problem.AddCriterion(
      [](const double* y) -> double { return (y[0] - 2)*(y[0] - 2);  }
    );
    problem.SetDomain(1, {-10.}, {10.});
  }
  else if(name == "fonseca")
  {
    if(dimension == -1) dimension = 2;

    problem.AddCriterion(
      [dimension](const double* y) -> double {
        double shift = 1. / sqrt(dimension);
        double eArg = 0;
        for(int i = 0; i < dimension; i++)
          eArg -= pow(y[i] - shift, 2);
        return 1 - exp(eArg);
      }
    );
    problem.AddCriterion(
      [dimension](const double* y) -> double {
        double shift = 1. / sqrt(dimension);
        double eArg = 0;
        for(int i = 0; i < dimension; i++)
          eArg -= pow(y[i] + shift, 2);
        return 1 - exp(eArg);
      }
    );
    problem.SetDomain(dimension, {-4., -4.}, {4., 4.});
  }

  return problem;
}
