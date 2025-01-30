/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2022 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_COUNTED_PTR_H_
#define _BRESEQ_COUNTED_PTR_H_

#include "libbreseq/common.h"

using namespace std;

namespace breseq
{
// counted_ptr keeps track of number of references 
  
  template <class X> class counted_ptr
  {
  public:
    
    //We have to keep a global here, otherwise we deallocate the mutex!
    
    typedef X element_type;
    
    explicit counted_ptr(X* p = 0) // allocate a new counter
    : itsCounter(0) {if (p) itsCounter = new counter(p);}
    ~counted_ptr()
    {release();}
    counted_ptr(const counted_ptr& r) throw()
    {
      acquire(r.itsCounter);
    }
    
    counted_ptr& operator=(const counted_ptr& r)
    {
      if (this != &r) {
        release();
        acquire(r.itsCounter);
      }
      return *this;
    }
    
    // Full suite of comparison operators defined for what we point to
    friend bool operator== (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
      {return *(lhs.itsCounter->ptr) == *(rhs.itsCounter->ptr);}
    friend bool operator!= (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
      {return !(lhs == rhs);}
    friend bool operator< (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
      {return *(lhs.itsCounter->ptr) < *(rhs.itsCounter->ptr);}
    friend bool operator> (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
      {return (rhs < lhs);}
    friend bool operator<= (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
      {return !(lhs > rhs);}
    friend bool operator>= (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
      {return !(lhs < rhs);}
    
    X& operator*()  const throw()   {return *itsCounter->ptr;}
    X* operator->() const throw()   {return itsCounter->ptr;}
    X* get()        const throw()   {return itsCounter ? itsCounter->ptr : 0;}
    bool unique()   const throw()
      {return (itsCounter ? itsCounter->count == 1 : true);}
    
  private:
    
    struct counter {
      counter(X* p = 0, unsigned c = 1) : ptr(p), count(c) {}
      X*          ptr;
      unsigned    count;
      mutex       mtx;
    }* itsCounter;
    
    void acquire(counter* c) throw()
    { // increment the count

      if (c) {
        c->mtx.lock();
        c->count++;
        c->mtx.unlock();
      }
      
      itsCounter = c;
    }
    
    void release()
    { // decrement the count, delete if it is 0
      if (itsCounter) {
        itsCounter->mtx.lock();
        itsCounter->count--;
        bool im_the_deleter = itsCounter->count == 0;
        itsCounter->mtx.unlock();
        if (im_the_deleter) {
          delete itsCounter->ptr;
          delete itsCounter;
          itsCounter = 0;
        }
      }
    }
  };
  
} // breseq

#endif
