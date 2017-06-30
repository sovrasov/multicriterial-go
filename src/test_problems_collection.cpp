#include "test_problems_collection.hpp"
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
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
          val -= 10.*exp(-0.2*hypot(y[i], y[i + 1]));
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
  if(name == "schaffer2")
  {
    problem.AddCriterion(
      [](const double* y) -> double {
        double x = y[0];
        if(x <= 1)
          return -x;
        else if(x <= 3)
          return x - 2;
        else if(x <= 4)
          return 4 - x;

        return x - 4;
       }
    );
    problem.AddCriterion(
      [](const double* y) -> double { return (y[0] - 5)*(y[0] - 5);  }
    );
    problem.SetDomain(1, {-5.}, {10.});
  }
  if(name == "poloni")
  {
    double A1 = 0.5*sin(1.) - 2.*cos(1) + sin(2.) - 1.5*cos(2.);
    double A2 = 1.5*sin(1.) - cos(1.) + 2.*sin(2.) - 0.5*cos(2.);
    problem.AddCriterion(
      [A1, A2](const double* y) -> double {
        double B1 = 0.5*sin(y[0]) - 2.*cos(y[0]) + sin(y[1]) - 1.5*cos(y[1]);
        double B2 = 1.5*sin(y[0]) - cos(y[0]) + 2.*sin(y[1]) - 0.5*cos(y[1]);
        return 1 + (A1 - B1)*(A1 - B1) + (A2 - B2)*(A2 - B2);
       }
    );
    problem.AddCriterion(
      [](const double* y) -> double { return (y[0] + 3)*(y[0] + 3) + (y[1] + 1)*(y[1] + 1);  }
    );
    problem.SetDomain(2, {-M_PI, -M_PI}, {M_PI, M_PI});
  }

  return problem;
}
