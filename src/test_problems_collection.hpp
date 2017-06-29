#pragma once

#include "mcoProblem.hpp"
#include <string>

class TestMCOProblems
{
public:
  static MCOProblem create(const std::string& name, int dimension = -1);
};
