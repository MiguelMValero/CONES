#ifndef __SINGLETON_H__
#define __SINGLETON_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>

namespace cwipi {
  template <typename T>
  class Singleton
  {
  protected:
    Singleton () { }
    virtual ~Singleton () {
#if defined(DEBUG) && 0
      std::cout << "destroying singleton." << std::endl;
#endif
 }
    Singleton (const Singleton & other);
    Singleton & operator=(const Singleton & other);

  public:
    static T &getInstance ()
    {
      if (NULL == _singleton)
        _singleton = new T;

      return *(static_cast<T*> (_singleton));
    }

    static void kill ()
    {
      if (NULL != _singleton) {
        delete _singleton;
        _singleton = NULL;
      }
    }

  private:
    static T *_singleton;
  };

  template <typename T>
  T *Singleton<T>::_singleton = NULL;
} // namespace cwipi

#endif
