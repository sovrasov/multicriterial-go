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
    std::vector<double> lb(dimension, -4.0);
    std::vector<double> ub(dimension, 4.0);
    problem.SetDomain(dimension, lb.data(), ub.data());
  }
  else if(name == "kursawe")
  {
    dimension = 3;
    problem.SetDomain(dimension, {-5, -5, -5}, {5, 5, 5});

    problem.AddCriterion(
      [dimension](const double* y) -> double {
        double val = 0.;
        for(int i = 0; i < dimension - 1; i++)
          val -= 10.*exp(-0.2*sqrt(y[i]*y[i] + y[i + 1]*y[i + 1]));
        return val;
      }
    );
    problem.AddCriterion(
      [dimension](const double* y) -> double {
        double val = 0.;
        for(int i = 0; i < dimension; i++)
          val += pow(fabs(y[i]), 0.8) + 5.*sin(pow(y[i], 3));
        return val;
      }
    );
  }

  return problem;
}
