#pragma once

const unsigned solverMaxDim = 5;
const unsigned solverMaxCriterions = 5;

struct Trial
{
  double x;
  double y[solverMaxDim];
  double h;
  double z[solverMaxCriterions];
  Trial() {}
  Trial(double _x) : x(_x){}
};

inline bool operator<(const Trial& t1, const Trial& t2)
{
  return t1.x < t2.x;
}

struct Interval
{
  Trial pl;
  Trial pr;
  double R;
  double delta;
  Interval() {}
  Interval(const Trial& _xl, const Trial& _xr) : pl(_xl), pr(_xr) {}
};

inline bool operator<(const Interval& i1, const Interval& i2)
{
  return i1.R < i2.R;
}
