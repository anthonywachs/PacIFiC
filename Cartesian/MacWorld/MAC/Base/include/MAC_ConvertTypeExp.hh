#ifndef MAC_CONVERT_TYPE_EXP_HH
#define MAC_CONVERT_TYPE_EXP_HH

#include <MAC_Expression.hh>

/* 
Expressions for type conversion

--
name     : double
argument : Int
type     : Double

double( i ) converts the integer i into a double using the standard
C conversion

example :
   double( 1 ) : value is 1.0

--
name     : int
argument : Double
type     : Int

int( x ) converts the double x into an int

example :
   int( 1.0 )   : value is 1
   int( 1.01 )  : value is 1
   int( 1.99 )  : value is 1
   int( -2.01 ) : value is -2
   int( -2.9 )  : value is -2
--

PUBLISHED
*/

class MAC_ConvertTypeExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( MAC_Context const* ct = 0 ) const ;

      virtual int to_int( MAC_Context const* ct = 0 ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_ConvertTypeExp( void ) ;
     ~MAC_ConvertTypeExp( void ) ;
      MAC_ConvertTypeExp( MAC_ConvertTypeExp const& other ) ;
      MAC_ConvertTypeExp& operator=( MAC_ConvertTypeExp const& other ) ;

      MAC_ConvertTypeExp( MAC_Object* a_owner,
                          std::string const& a_name,
                          MAC_Sequence const* argument_list ) ;

   //-- Plug in

      MAC_ConvertTypeExp( std::string const& a_name ) ;

      virtual MAC_ConvertTypeExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      
      static MAC_ConvertTypeExp const* PROTOTYPE_DOUBLE ;
      static MAC_ConvertTypeExp const* PROTOTYPE_INT ;
      
   //-- Attributes

      MAC_Data const* const P0 ;
} ;

#endif
