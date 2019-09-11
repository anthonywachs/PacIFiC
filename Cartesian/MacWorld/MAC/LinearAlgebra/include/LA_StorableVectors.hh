#ifndef LA_STORABLEVECTORS_HH
#define LA_STORABLEVECTORS_HH

#include <MAC_Object.hh>

#include <string>
#include <list>
using std::list;

class LA_Vector;

/*
Object to store (distributed) vectors for reload

PUBLISHED
*/

class LA_StorableVectors : public MAC_Object
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      static LA_StorableVectors* create( MAC_Object * a_owner ) ; 
      
   //-- Persistence

      // Use `writer' to store `self' so that ist can be retrieved
      // with `::restore_state'.
      virtual void save_state( MAC_ObjectWriter* writer ) const ;
      
      // Retrieve `self' from `reader'.
      virtual void restore_state( MAC_ObjectReader* reader ) ;
      
      // Add a vector to be stored
      void add_vector_to_store( LA_Vector* vec );
                       
   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      LA_StorableVectors( void ) ;
     ~LA_StorableVectors( void ) ;
      LA_StorableVectors( LA_StorableVectors const& other ) ;
      LA_StorableVectors& operator=( LA_StorableVectors const& other ) ;

      LA_StorableVectors( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      list<LA_Vector*> STORED_VECTORS;
      
} ;


#endif

