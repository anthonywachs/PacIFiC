#ifndef MAC_EXCEPTIONS_HH
#define MAC_EXCEPTIONS_HH

/*
PUBLISHED
*/

class MAC_Exceptions
{
   public: //----------------------------------------------------------------
   
      class Error { } ;
      
      class InternalError : public Error { } ;
      
      class UserError : public Error { } ;
      
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

      MAC_Exceptions( void ) ;
     ~MAC_Exceptions( void ) ;
      MAC_Exceptions( MAC_Exceptions const& other ) ;
      MAC_Exceptions& operator=( MAC_Exceptions const& other ) ;
} ;

#endif
