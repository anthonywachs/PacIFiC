#ifndef LA_SEQ_IMPLEMENTATION_HH
#define LA_SEQ_IMPLEMENTATION_HH

#include <LA_Implementation.hh>

/*
Implementation of MAC sequential matrices.

PUBLISHED
*/

class LA_SeqImplementation : public LA_Implementation
{
   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization
      
      static LA_SeqImplementation const* object( void ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

     ~LA_SeqImplementation( void ) ;
      LA_SeqImplementation( void ) ;
      
      LA_SeqImplementation( LA_SeqImplementation const& other ) ; 
      LA_SeqImplementation& operator=( LA_SeqImplementation const& other ) ;
     
} ;

#endif

