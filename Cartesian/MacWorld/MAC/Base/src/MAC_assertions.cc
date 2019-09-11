#include <MAC_assertions.hh>

#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC.hh>

#include <iostream>

#ifdef OUTLINE
#define inline
#include <MAC_assertions.icc>
#undef inline
#endif

//no doc---------------------------------------------------------------------
MAC_Assertion::CheckType MAC_Assertion::current_check = None ;
MAC_Assertion MAC_Assertion::unique_instance ;
bool MAC_Assertion::negation ;
bool MAC_Assertion::result ;
bool MAC_Assertion::eval ;
bool MAC_Assertion::short_cut ;
bool MAC_Assertion::bool_table[ MAXBOOL ] ;
int MAC_Assertion::checking_level = 0 ;
size_t MAC_Assertion::nb_bool = 0 ;
//no doc---------------------------------------------------------------------

//no doc---------------------------------------------------------------------
bool
MAC_Assertion:: test_implement_handling( char const* file,
                                         int line,
                                         std::string const& text )
//no doc---------------------------------------------------------------------
{
   MAC_Error::object()->raise_not_tested( file, line, text ) ;
   return false ;
}

//no doc---------------------------------------------------------------------
bool
MAC_Assertion:: action( const char* file, int line, const char* text )
//no doc---------------------------------------------------------------------
{
   CheckType check = current_check ;
   current_check = None ;
   switch( check )
   {
      case Check :
         MAC_Error::object()->raise_assertion_violation( file, line, text ) ;
         break ;
      case Precondition :
	 MAC_Error::object()->raise_precondition_violation( text ) ;
         break ;
      case Postcondition :
	 MAC_Error::object()->raise_postcondition_violation( file, line, text ) ;
         break ;
      case Invariant :
	 MAC_Error::object()->raise_invariant_violation( file, line, text ) ;
         break ;
      default :
	 MAC_Error::object()->raise_assertion_violation( file, line, text ) ;
         break ;
   }
   return false ;
}

//no doc---------------------------------------------------------------------
const char* MAC_Marker::ring[ 256 ] ;
size_t MAC_Marker::ring_pos = 0 ;
//no doc---------------------------------------------------------------------

//no doc---------------------------------------------------------------------
void
MAC_Assertion:: add_handled_check( CheckType a_check )
//no doc---------------------------------------------------------------------
{
   checking_level |= a_check ;
   MAC_ASSERT( true && ( checking_level & a_check ) ) ;
}

//no doc---------------------------------------------------------------------
bool
MAC_Marker:: is_collective( int line )
//no doc---------------------------------------------------------------------
{   
   static int order = 1 ;
   
   bool result = MAC_Exec::communicator()->same_value_everywhere((int)line*order) ;
   order++ ;

   return result ;
}
