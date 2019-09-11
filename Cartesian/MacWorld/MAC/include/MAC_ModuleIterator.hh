#ifndef MAC_MODULE_ITERATOR_HH
#define MAC_MODULE_ITERATOR_HH

#include <MAC_ListIterator.hh>
#include <MAC_Module.hh>

/*
   Iterators on MAC_Module objects.

   Return the `MAC_Module::' objects stored in a `MAC_List::'. 
*/

class MAC_ModuleIterator : public MAC_ListIterator
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Create and return an instance for traversing `a_list' whose
      // items are `MAC_Module::' objects.
      static MAC_ModuleIterator* create( MAC_Object* a_owner,
                                         MAC_List const* a_list ) ;
      
   //-- Access

      virtual MAC_Module* item( void ) const ;
      
   protected: //------------------------------------------------------------
      
   private: //--------------------------------------------------------------

      MAC_ModuleIterator( void ) ;
     ~MAC_ModuleIterator( void ) ;
      MAC_ModuleIterator( MAC_ModuleIterator const& other ) ;
      MAC_ModuleIterator const& operator=( MAC_ModuleIterator const& other ) ;
      MAC_ModuleIterator( MAC_Object* a_owner,
                          MAC_List const* a_list ) ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
} ;

#endif
