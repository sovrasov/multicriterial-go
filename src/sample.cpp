#include "multicriterialSolver.hpp"
#include "test_problems_collection.hpp"

#include <iostream>
#include <fstream>
#include <cmdline.h>
#include <chrono>
#include <string>

void saveOptimalPoints(const std::string& filename, std::vector<Trial> points, int dim, int m)
{
  std::ofstream fout;
  fout.open(filename, std::ios_base::out);

  fout << dim << "\n" << m << "\n";

  for(auto& x : points)
  {
    for(int i = 0; i < dim; i++)
      fout << x.y[i] << ", ";
    for(int i = 0; i < m - 1; i++)
      fout << x.z[i] << ", ";
    fout << x.z[m - 1] << ";\n";
  }

  fout.close();
}

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
  parser.add<double>("reserves", 'E', "eps-reserves for all constraints", false, 0);
  parser.add<int>("itersLimit", 'l', "limit of iterations for the method", false, 2000);
  parser.add<int>("dim", 'd', "test problem dimension (will be set if supported)", false, 2);
  parser.add<int>("localMix", 'm', "local mix parameter", false, 0);
  parser.add<std::string>("outFile", 'f', "name of the file to write solution",
    false, "solution.csv");
  parser.add<std::string>("probName", 'n', "name of the test problem",
    false, "fonseca");
  parser.add("saveSolution", 's', "determines whether the method will "
    "save solution into a .csv file");
  parser.add("computeLoad", 'c', "make test functions hard to calculate (for experiments)");
  parser.parse_check(argc, argv);

  MCOProblem problem = TestMCOProblems::create(
    parser.get<std::string>("probName"), parser.get<int>("dim"), parser.exist("computeLoad"));
  if(problem.GetCriterionsNumber() == 0)
  {
    std::cerr << "Wrong test problem name\n";
    exit(0);
  }

  MCOSolver solver;
  solver.SetParameters(SolverParameters(parser.get<double>("accuracy"),
    parser.get<double>("reserves"),
    parser.get<double>("reliability"),
    parser.get<int>("threadsNum"),
    parser.get<int>("itersLimit"),
    parser.get<int>("localMix")));
  solver.SetProblem(problem);

  auto start = std::chrono::system_clock::now();
  solver.Solve();
  auto solution = solver.GetWeakOptimalPoints();
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Time elapsed: " << elapsed_seconds.count() << "s\n";

  std::cout << "Iterations performed: " << solver.GetIterationsNumber() << std::endl;
  std::cout << "Number of weak optimal points: " << solution.size() << std::endl;

  if(parser.exist("saveSolution"))
    saveOptimalPoints(parser.get<std::string>("outFile"), solution,
      problem.GetDimension(), problem.GetCriterionsNumber());

  return 0;
}
