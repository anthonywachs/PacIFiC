#ifndef MAC_VECTOR_HH
#define MAC_VECTOR_HH

#include <MAC_Sequence.hh>

#include <MAC_VectorIterator.hh>

/*
Sequences whose items access from their index if efficient
*/

class MAC_Vector : public MAC_Sequence
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create an return an instance.
      static MAC_Vector* create( MAC_Object* a_owner, size_t size ) ;

      // Reinitialize the internal state, as if the `::create' method
      // was just completed.
      void re_initialize( size_t size ) ;

      virtual MAC_Vector* create_clone( MAC_Object* a_owner ) const ;

      // Reinitialize by copying all the items of `other'.
      void copy( MAC_Vector const* other ) ;
      
   //-- Resizing

      // Change the exclusive upper limit for indices, without losing
      // previously entered items whose index was within the new bound,
      // nor changing their coordinate.
      void resize( size_t size ) ;

  //-- Measurement

      virtual size_t index_limit( void ) const ;

      virtual size_t count( void ) const ;
      
   //-- Access
      
      virtual MAC_Object* item( MAC_Object const* object ) const ;

      virtual MAC_Object* at( size_t i ) const ;

      virtual size_t index_of( MAC_Object const* object ) const ;

      virtual MAC_VectorIterator* create_iterator( MAC_Object* a_owner ) const ;

   //-- Element change

      // Make all items be the 0 pointer.
      void nullify( void ) ;
      
      virtual void append( MAC_Object* object ) ;

      virtual void prepend( MAC_Object* object ) ;

      virtual void set_at( size_t i, MAC_Object* object ) ;

      virtual void insert_at( size_t i, MAC_Object* object ) ; 

   //-- Removal
      
      virtual void remove_at( size_t i ) ;
      
      virtual void remove_section( size_t iFirst, size_t length ) ;

      virtual void destroy_items_and_remove_section( size_t iFirst, 
                                                 size_t length ) ;
      virtual void clear( void ) ;

         
   protected: //--------------------------------------------------------

      virtual ~MAC_Vector( void );

      MAC_Vector( MAC_Object* a_owner, size_t size ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

      virtual bool remove_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t old_state_id ) const ;

      virtual bool remove_section_POST( size_t old_index_limit,
                                        size_t old_count,
                                        size_t length,
                                        size_t old_state_id  ) const ;

  private: //----------------------------------------------------------

      MAC_Vector( void ) ;
      MAC_Vector( MAC_Vector const& other ) ;
      MAC_Vector& operator=( MAC_Vector const& other ) ;
   
   //-- Attributes

      MAC_Object** VECTOR ;
      size_t LENGTH ;
      size_t NB_ENTRIES ;
      size_t CAPACITY ;
} ;

#endif
