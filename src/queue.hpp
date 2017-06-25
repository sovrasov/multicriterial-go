#pragma once

#include "minmaxheap.hpp"
#include "dataTypes.hpp"
#include <memory>

class IntervalsQueue
{
protected:
  std::shared_ptr<MinMaxHeap<Interval>> mHeapPtr;

public:
  void resize(size_t size);
  size_t size() const;
  void clear();

  void pushWithPriority(const Interval& i);
  Interval pop();
};
