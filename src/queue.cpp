#include "queue.hpp"

void IntervalsQueue::resize(size_t size)
{
  mHeapPtr = std::make_shared<MinMaxHeap<Interval>>(size);
}

void IntervalsQueue::pushWithPriority(const Interval& i)
{
  if(!mHeapPtr->empty())
  {
    if(i.R > mHeapPtr->findMin().R)
    {
      if(mHeapPtr->is_full())
        mHeapPtr->popMin();
      mHeapPtr->push(i);
    }
  }
  else
    mHeapPtr->push(i);
}

void IntervalsQueue::clear()
{
  mHeapPtr->clear();
}

size_t IntervalsQueue::size() const
{
  return mHeapPtr->size();
}

Interval IntervalsQueue::pop()
{
  return mHeapPtr->pop();
}
