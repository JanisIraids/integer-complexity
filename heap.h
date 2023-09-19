#pragma once

#include <utility>

// A simple min heap with dummy elements for min and max.
template <class KEY, class VALUE, bool (*COMP)(const KEY &, const KEY &)>
class MinHeap
{
  struct mp
  {
    KEY k;
    VALUE v;
  };
  mp* data;
  KEY minguard, maxguard;
  unsigned int currsize;

  static unsigned int parent(const unsigned int index) __attribute__((const))
  {
    return (index >> 1);
  }

  static unsigned int leftch(const unsigned int index) __attribute__((const))
  {
    return (index<<1);
  }

  static unsigned int rightch(const unsigned int index) __attribute__((const))
  {
    return (index<<1) + 1;
  }

  unsigned int moveup(unsigned int index)
  {
    unsigned int pindex = parent(index);
    mp vatindex = data[index];
    while (COMP(vatindex.k, data[pindex].k))
    {
      data[index] = data[pindex];
      index = pindex;
      pindex = parent(index);
    }
    data[index] = vatindex;
    return index;
  }

  unsigned int movedown(unsigned int index)
  {
    mp vatindex = data[index];
    unsigned int mcindex;
    if (COMP(data[leftch(index)].k, data[rightch(index)].k))
      mcindex = leftch(index);
    else
      mcindex = rightch(index);
    while (COMP(data[mcindex].k, vatindex.k))
    {
      data[index] = data[mcindex];
      index = mcindex;
      if (COMP(data[leftch(index)].k, data[rightch(index)].k))
	mcindex = leftch(index);
      else
	mcindex = rightch(index);
    }
    data[index] = vatindex;
    return index;
  }

public:
  MinHeap(const KEY max, const KEY min, unsigned int capacity = 1000000)
  {
    minguard = min;
    maxguard = max;
    unsigned int realcap = 2*capacity+2;
    data = new mp[realcap];
    currsize = 0;
    data[0].k = minguard;
    for (unsigned int i = 1; i < realcap; i++)
      data[i].k = maxguard;
  }

  ~MinHeap()
  {
    delete[] data;
  }

  KEY getK(const unsigned int index) const
  {
    return data[index].k;
  }

  VALUE getV(const unsigned int index) const
  {
    return data[index].v;
  }

  unsigned int top() const
  {
    return 1;
  }

  unsigned int size() const
  {
    return currsize;
  }

  unsigned int insert(const KEY k, const VALUE v)
  {
    ++currsize;
    data[currsize].k = k;
    data[currsize].v = v;
    return moveup(currsize);
  }

  void removeTop()
  {
    data[1] = data[currsize];
    data[currsize].k = maxguard;
    currsize--;
    movedown(1);
  }

  void remove(const unsigned int index)
  {
    data[index] = data[currsize];
    data[currsize].k = maxguard;
    currsize--;
    movedown(moveup(index));
  }
  
  unsigned int decrease(const unsigned int index, KEY newvalue)
  {
    data[index].k = newvalue;
    return moveup(index);
  }

  unsigned int increase(const unsigned int index, KEY newvalue)
  {
    data[index].k = newvalue;
    return movedown(index);
  }
};
