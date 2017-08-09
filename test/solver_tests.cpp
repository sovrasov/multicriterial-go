#include "gtest/gtest.h"

#include <multicriterialSolver.hpp>
#include <test_problems_collection.hpp>

#define RUN_SOLVER(problem) \
MCOProblem problem = TestMCOProblems::create(#problem); \
MCOSolver solver; \
solver.SetParameters(parameters); \
solver.SetProblem(problem); \
solver.Solve(); \
auto solution = solver.GetWeakOptimalPoints();

TEST(Solver, unconstrained_run_smoke)
{
	auto parameters = SolverParameters(0.01, 0,	4, 1,	2000,	0); \
	RUN_SOLVER(strongin);

	ASSERT_LE(solver.GetIterationsNumber(), 2000);
	ASSERT_GE(solution.size(), 100u);
}

TEST(Solver, unconstrained_parallel_run_smoke)
{
	auto parameters = SolverParameters(0.01, 0,	4, 4,	2000,	0); \
	RUN_SOLVER(strongin);

	ASSERT_LE(solver.GetIterationsNumber(), 500);
	ASSERT_GE(solution.size(), 100u);
}

TEST(Solver, constrained_run_smoke)
{
	auto parameters = SolverParameters(0.01, 0,	4, 1,	2000,	0); \
	RUN_SOLVER(chakong);

	ASSERT_LE(solver.GetIterationsNumber(), 1000);
	ASSERT_GE(solution.size(), 80u);
}

TEST(Solver, constrained_parallel_run_smoke)
{
	auto parameters = SolverParameters(0.01, 0,	4, 4,	2000,	0); \
	RUN_SOLVER(chakong);

	ASSERT_LE(solver.GetIterationsNumber(), 250);
	ASSERT_GE(solution.size(), 80u);
}
