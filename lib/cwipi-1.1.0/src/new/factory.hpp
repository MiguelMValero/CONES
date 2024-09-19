#ifndef __CWIPI_FACTORY_
#define __CWIPI_FACTORY_

/*
  This file is part of the CWIPI library. 

  Copyright (C) 2021-2023  ONERA

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

#include <map>
#include "singleton.hpp"

namespace cwipi {

/**
 * \defgroup	FactoryErrorPoliciesGroup Factory Error Policies
 * \ingroup		FactoryGroup
 * \brief		Manages the "Unknown Type" error in an object factory
 *
 
 *
 * \class DefaultFactoryError
 * \ingroup		FactoryErrorPoliciesGroup
 * \brief		Default policy that throws an exception		
 * 
 */

  template <typename IdentifierType, class AbstractProduct>
  struct DefaultFactoryError {
  
 /**
 * \struct Exception
 * \ingroup		FactoryErrorPoliciesGroup
 * \brief		Throws an exception		
 * 
 */ 
    struct Exception : public std::exception
    {
    
     /**
      *  \brief	what function		
      * 
      *  \return     "Unknown Type"
      */ 
      const char* what() const throw() { return "Unknown Type"; }
    };

 /**
 *  \brief	Throws an exception		
 * 
 *  \return     AbstractProduct*
 */     
    static AbstractProduct* OnUnknownType(IdentifierType)
    {
      throw Exception();
    }
  };

/**
 *  \class Factory
 *
 *  \ingroup FactoryGroup
 *  Implements a generic object factory.
 *
 */

  template <
    class                           AbstractProduct,
    typename                        IdentifierType,
    typename                        ProductCreator = AbstractProduct* (*)(),
    template<typename, class> class FactoryErrorPolicy = DefaultFactoryError
    >

  class Factory : 
    public Singleton < Factory <
                         AbstractProduct,
                         IdentifierType,
                         ProductCreator,
                         FactoryErrorPolicy > >,
    public FactoryErrorPolicy < IdentifierType, AbstractProduct >  
  {
    friend class Singleton  < Factory <
                                AbstractProduct,
                                IdentifierType,
                                ProductCreator,
                                FactoryErrorPolicy > >;
 
  public:

    /**
     *  \brief Register a new couple id / concrete product
     *
     *  \tparam     ConcreteProduct   Concrete product            
     *  \param [in] id                Identifier
     *
     */

    template < class ConcreteProduct >
    bool Register(const IdentifierType& id) 
    {
      return associations_.insert(
        typename IdToProductMap::value_type(id, 
                               newCreator<ConcreteProduct>)).second != 0;
    }
    
    /**
     *  \brief Unregister a concrete product
     *
     *  \param [in] id                Identifier
     *
     */

    bool Unregister(const IdentifierType& id)
    {
      return associations_.erase(id) != 0;
    }
    
    /**
     *  \brief Create an object of Concrete product from it identifier
     *
     *  \param [in] id                Identifier
     *
     */

    AbstractProduct* CreateObject(const IdentifierType& id)
    {
      typename IdToProductMap::iterator i = associations_.find(id);
      if (i != associations_.end())
        {
          return (i->second)();
        }
      return this->OnUnknownType(id);
    }

  private:
    
    typedef std::map<IdentifierType, ProductCreator> IdToProductMap;   /*!< container type (map) */
    IdToProductMap associations_;                                      /*!< container */
  
    /**
     *  \brief Create an object of Concrete product from it identifier
     *
     *  \tparam     ConcreteProduct   Concrete product            
     *  \return     A new concrete product with the abstract product interface
     */

    template <class ConcreteProduct>
    static AbstractProduct * newCreator() {
      return new ConcreteProduct;
  }

  };  
  
}

#endif
