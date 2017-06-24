#ifndef EVOLVENT_HPP
#define EVOLVENT_HPP

#define MAX_PREIMAGES 32

enum MapType {
  Simple = 1, Linear = 2, Noninjective = 3
};

class Evolvent
{
protected:
  int mDimension;
  int mTightness;
  bool mIsInitialized;

private:
  MapType mMapType;
  int mMapKey;

public:
  Evolvent();
  Evolvent(int dimension, int tightness, MapType type = Simple);
  ~Evolvent();

  void GetImage(double x, double y[]);
  void GetImage(double x, double y[], double lb[], double ub[]);
  int GetAllPreimages(double* p, double xp[]);
};


#endif
