#ifndef MAC_INDEX_SET_HH
#define MAC_INDEX_SET_HH

#include <MAC_Object.hh>

#include <size_t_vector.hh>

class MAC_IndexSet : public MAC_Object
{
      
   public : //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_IndexSet* create( MAC_Object* a_owner ) ;
      
      static MAC_IndexSet* create( MAC_Object* a_owner,
                                   size_t_vector const& vec,
                                   size_t a_id ) ;
      
      void re_initialize( size_t_vector const& vec, 
                          size_t a_id ) ;
      
   //-- Identifier

      size_t id( void ) const ;
      
   //-- Element access

      // elements in increasingly order
      size_t_vector const& elements( void ) const ;

   //-- Comparison

      virtual bool is_equal( MAC_Object const* other ) const ;
      
      virtual int three_way_comparison( MAC_Object const* other ) const ;
      
      virtual size_t hash_code( void ) const ;    

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected ://------------------------------------------------------------

   private :  //------------------------------------------------------------

      MAC_IndexSet( void ) ;
     ~MAC_IndexSet( void ) ;
      MAC_IndexSet( MAC_IndexSet const& other ) ;
      MAC_IndexSet operator=( MAC_IndexSet const& other ) ;

      MAC_IndexSet( MAC_Object* a_owner ) ;
      
      MAC_IndexSet( MAC_Object* a_owner,
                    size_t_vector const& vec,
                    size_t a_id ) ;

   //-- Attributes

      size_t ID ;
      size_t_vector SET ;
} ;

#endif
