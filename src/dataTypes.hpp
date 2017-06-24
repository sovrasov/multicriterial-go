#ifndef DATA_TYPES_HPP
#define DATA_TYPES_HPP

const unsigned solverMaxDim = 5;

struct Trial
{
  double x;
  double y[solverMaxDim];
  double z;
  Trial() {}
  Trial(double _x, double _z) : x(_x), z(_z) {}
};

struct Interval
{
  double xl;
  double xr;
  double zl;
  double zr;
  double R;
  double delta;
  Interval() {}
  Interval(double _xl, double _xr) : xl(_xl), xr(_xr) {}
};

inline bool operator<(const Interval& i1, const Interval& i2)
{
  return i1.xl < i2.xl;
}

class CompareByR
{
public:
  bool operator() (const Interval* i1, const Interval* i2)
  {
    return i1->R < i2->R;
  }
};

#endif
