#include "test_problems_collection.hpp"
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

namespace
{
  double compute_load(int complexity = 100)
  {
    double value = 0.;

    for (int i = 0; i < complexity * 1000; i++)
    {
      double a1 = fma(value, 2., value + 1.);
      double a2 = fma(value, 1., a1);
      value = a2 - a1;
    }

    return value + 1.;
  }
}

MCOProblem TestMCOProblems::create(const std::string& name, int dimension, bool useComputeLoad)
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
      [dimension, useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;
        double shift = 1. / sqrt(dimension);
        double eArg = 0;
        for(int i = 0; i < dimension; i++)
          eArg -= pow(y[i] - shift, 2);
        return dummy*(1 - exp(eArg)) / dummy2;
      }
    );
    problem.AddCriterion(
      [dimension, useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;
        double shift = 1. / sqrt(dimension);
        double eArg = 0;
        for(int i = 0; i < dimension; i++)
          eArg -= pow(y[i] + shift, 2);
        return dummy*(1 - exp(eArg)) / dummy2;
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
      [useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;

        double x = dummy*y[0] / dummy2;
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
      [useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;
        return dummy*(y[0] - 5)*(y[0] - 5) / dummy2;
       }
    );
    problem.SetDomain(1, {-5.}, {10.});
  }
  if(name == "poloni")
  {
    double A1 = 0.5*sin(1.) - 2.*cos(1) + sin(2.) - 1.5*cos(2.);
    double A2 = 1.5*sin(1.) - cos(1.) + 2.*sin(2.) - 0.5*cos(2.);
    problem.AddCriterion(
      [A1, A2, useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;
        double B1 = 0.5*sin(y[0]) - 2.*cos(y[0]) + sin(y[1]) - 1.5*cos(y[1]);
        double B2 = 1.5*sin(y[0]) - cos(y[0]) + 2.*sin(y[1]) - 0.5*cos(y[1]);
        return dummy*(1 + (A1 - B1)*(A1 - B1) + (A2 - B2)*(A2 - B2)) / dummy2;
       }
    );
    problem.AddCriterion(
      [useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;
        return dummy*(y[0] + 3)*(y[0] + 3) + (y[1] + 1)*(y[1] + 1) / dummy2;
      }
    );
    problem.SetDomain(2, {-M_PI, -M_PI}, {M_PI, M_PI});
  }
  if(name == "strongin")
  {
    problem.AddCriterion(
      [useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;
        return dummy*fmin(hypot(y[0], y[1]) + .5, hypot(y[0] - 1.5, y[1] + 1.5)) / dummy2;
      }
    );
    problem.AddCriterion(
      [useComputeLoad](const double* y) -> double {
        double dummy = useComputeLoad ? compute_load() : 1.;
        double dummy2 = useComputeLoad ? compute_load() : 1.;
        return dummy*hypot(y[0] + .5, y[1] - .5) / dummy2;
      }
    );
    problem.SetDomain(2, {-1., -2.}, {2., 1.});
  }
  if(name == "oka2")
  {
    problem.AddCriterion(
      [](const double* y) -> double {
        return y[0];
      }
    );
    problem.AddCriterion(
      [](const double* y) -> double {
        return 1 - pow((y[0] + M_PI)/ (2*M_PI), 2) +
          cbrt(fabs(y[1] - 5*cos(y[0]))) + cbrt(fabs(y[2] - 5*sin(y[1])));
      }
    );
    problem.SetDomain(3, {-M_PI, -5, -5}, {M_PI, 5., 5.});
  }
  if(name == "oka1")
  {
    const double cos_coeff = cos(M_PI / 12);
    const double sin_coeff = sin(M_PI / 12);

    problem.AddCriterion(
      [cos_coeff, sin_coeff](const double* y) -> double {
        return cos_coeff*y[0] - sin_coeff*y[1];
      }
    );
    problem.AddCriterion(
      [cos_coeff, sin_coeff](const double* y) -> double {
        const double x1 = cos_coeff*y[0] - sin_coeff*y[1];
        const double x2 = sin_coeff*y[0] + cos_coeff*y[1];

        return sqrt(2*M_PI) - sqrt(fabs(x1)) + 2*cbrt(fabs(x2 - 3*cos(x1) -3.));
      }
    );
    problem.SetDomain(2, {6*sin_coeff, -2*M_PI*sin_coeff},
       {6*sin_coeff + 2*M_PI*cos_coeff, 6*cos_coeff});
  }
  if(name == "evtushenko")
  {
    problem.AddCriterion(
      [](const double* y) -> double {
        return (y[0] - 1)*y[1]*y[1] + 1;
      }
    );
    problem.AddCriterion(
      [](const double* y) -> double {
        return y[1];
      }
    );
    problem.SetDomain(2, {0, 0}, {1, 1});
  }
  if(name == "mop6")
  {
    problem.AddCriterion(
      [](const double* y) -> double {
        return y[0];
      }
    );
    problem.AddCriterion(
      [](const double* y) -> double {
        const double t = 1 + 10*y[1];
        return t*(1 - pow(y[1] / t, 2) - y[0]* sin(8*M_PI*y[0]) / t);
      }
    );
    problem.SetDomain(2, {0, 0}, {1, 1});
  }
  if(name == "viennet")
  {
    problem.AddCriterion(
      [](const double* y) -> double {
        const double sqr = y[0]*y[0] + y[1]*y[1];
        return 0.5*sqr + sin(sqr);
      }
    );
    problem.AddCriterion(
      [](const double* y) -> double {
        return pow( 3*y[0] - 2*y[1] + 4, 2) / 8 + pow(y[0] - y[1] + 1, 2) / 27 + 15;
      }
    );
    problem.AddCriterion(
      [](const double* y) -> double {
        const double sqr = y[0]*y[0] + y[1]*y[1];
        return 1. / (sqr + 1) - 1.1*exp(-sqr);
      }
    );
    problem.SetDomain(2, {-3, -3}, {3, 3});
  }

  return problem;
}
