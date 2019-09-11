#include <MAC_Object.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Iterator.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_Module.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>

#include <typeinfo>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <map>

using std::map ;
using std::string ;

static map<MAC_Object const*,size_t>* _MAC_OBJECT_LIST = 0 ;
static size_t _RANK = 0 ;
static MAC_Object const* catched_object = 0 ;
static size_t catched_object_rank = (size_t) ~0 ;

static map<MAC_Object const*,size_t>* _MAC_OBJECT_LIST2 = 0 ;
static bool trace_allocation = false ;

//----------------------------------------------------------------------
size_t MAC_Object:: ALLOCATED = 0 ;
//----------------------------------------------------------------------




//----------------------------------------------------------------------
MAC_Object:: MAC_Object( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MY_OWNER( 0 ),
     POSSESSIONS( 0 )   
{
   
   // Construct the object
   ALLOCATED++ ;
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Objects ) )
   {
      static bool prem = true ;
      if( prem )
      {
         _MAC_OBJECT_LIST2 = new map<MAC_Object const*,size_t> ;
         _MAC_OBJECT_LIST = new map<MAC_Object const*,size_t> ;
         _RANK = 0 ;
         prem = false ;
      }
      if( a_owner==0 )
      {
         _MAC_OBJECT_LIST->operator[]( this ) = _RANK ;
      }
      if( catched_object==this )
      {
         MAC_Error::object()->trace( "Catched object creation" ) ;
      }
      if( _RANK==catched_object_rank )
      {
         MAC_Error::object()->trace( "Catched object creation" ) ;
         catched_object=this ;
      }
      if( trace_allocation )
      {
         _MAC_OBJECT_LIST2->operator[]( this ) = _RANK ;
      }
      _RANK++ ;
   }
   if( a_owner != 0 )
   {
      a_owner->insert_possession( this ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_Object:: ~MAC_Object( void )
//----------------------------------------------------------------------
{
   if( MY_OWNER!=0 &&  owner()->POSSESSIONS!=0 )
   {
      MY_OWNER->POSSESSIONS->remove( this ) ;
      MY_OWNER=0 ;
   }
   
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Objects ) )
   {
      if( owner() == 0 )
      {
         _MAC_OBJECT_LIST->erase( this ) ;
      }
      if( catched_object==this )
      {
         MAC_Error::object()->trace( "Catched object deletion" ) ;
      }
      if( trace_allocation )
      {
         _MAC_OBJECT_LIST2->erase( this ) ;
      }
   }
   if( POSSESSIONS != 0 )
   {
      MAC_ListIterator *it = POSSESSIONS->create_iterator( 0 ) ;
      for( ; it->is_valid() ; it->go_next() )
      {
         MAC_Object* child = it->item() ;
         child->MY_OWNER = 0 ;
         delete child ;
      }
      it->destroy() ;
      POSSESSIONS->destroy() ; POSSESSIONS=0 ;
   }
   ALLOCATED-- ;
}




//----------------------------------------------------------------------
MAC_Object*
MAC_Object:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: create_clone" ) ;
   MAC_Error::object()->raise_not_implemented( this, "create_clone" ) ;

   MAC_Object* result = 0 ;
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: update( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: update" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Error::object()->raise_not_implemented( this, "update" ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: destroy( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: destroy" ) ;
   MAC_CHECK_PRE( owner()==0 ) ;
   
   delete this ;
}




//----------------------------------------------------------------------
void
MAC_Object:: destroy_possession( MAC_Object const* a_possession )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: destroy_possession" ) ;
   MAC_CHECK_PRE( a_possession->owner()==this ) ;

   delete a_possession ;
}




//----------------------------------------------------------------------
size_t
MAC_Object:: address( void ) const
//----------------------------------------------------------------------
{
   return( (size_t)(long) this  ) ; //????????????????????????????????
}

bool valid_first( char c )
{
   return isalpha(c)!=0 ;
}




//----------------------------------------------------------------------
std::string const&
MAC_Object:: type_name( void ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: type_name" ) ;
   MAC_CHECK_INV( invariant() ) ;

   static string result ;

   string nn = typeid(*this).name() ;

   // on supprime quelques chiffres qui sont accolés au début
   // du nom de la classe par le compilateur gcc
   string::iterator it = std::find_if( nn.begin(), nn.end(), valid_first ) ;
   if( it == nn.end() )
   {
      result = nn ;
   }
   else
   {
      result.assign( it, nn.end() ) ;
   }
   // on supprime egalement le "class" ajoute par le compilateur VC++
   size_t idx ;
   if( ( idx = result.find( "class " ) ) < result.length() )
   {
      result = result.substr( idx+6 ) ;
   }
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Object const*
MAC_Object:: owner( void ) const
//----------------------------------------------------------------------
{
   return( MY_OWNER ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: is_under_ownership_of( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: is_under_ownership_of" ) ;
   MAC_CHECK_PRE( other!=this ) ;
   MAC_CHECK_INV( invariant() ) ;
   bool result = false ;
   if( MY_OWNER==0 )
   {
      result = ( other==0 ) ;
   }
   else if( MY_OWNER==other )
   {
      result = true ;
   }
   else
   {
      result = MY_OWNER->is_under_ownership_of( other ) ;
   }
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: same_type( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: same_type" ) ;
   MAC_CHECK_PRE( other != 0 ) ;

   bool result = true ;
// ( typeid(*this) == typeid(*other) ) ;

   MAC_CHECK_POST( 
       FORMAL( EQUIVALENT( result, typeid(*this) == typeid(*other ) ) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: set_owner( MAC_Object* a_owner )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: set_owner" ) ;
   MAC_CHECK_PRE( a_owner != 0 ) ;
   MAC_CHECK_PRE( owner() == 0 ) ;

   a_owner->insert_possession( this ) ;

   if( MAC_Assertion::is_handling_check( MAC_Assertion::Objects ) )
   {
      _MAC_OBJECT_LIST->erase( this ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( owner() == a_owner ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: change_owner( MAC_Object* a_owner,
                           MAC_Object* a_possession )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: change_owner" ) ;
   MAC_CHECK_PRE( a_possession->is_under_ownership_of(this) ) ;

   MAC_Object* father = a_possession->MY_OWNER ;
   
   father->POSSESSIONS->remove( a_possession ) ;
   a_possession->MY_OWNER = 0 ;
   if( a_owner!=0 )
   {
      a_possession->set_owner( a_owner ) ;
   }

   MAC_CHECK_POST( a_possession->owner()==a_owner ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: comparable( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   return( same_type( other ) ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: is_equal( MAC_Object const* other ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: is_equal" ) ;
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;

   bool result = has_same_address( other ) ;

   MAC_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
int
MAC_Object:: three_way_comparison( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: three_way_comparison" ) ;
   MAC_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   MAC_Error::object()->raise_not_implemented( this, "three_way_comparison" ) ;
   int result = -1 ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t
MAC_Object:: hash_code( void ) const 
//----------------------------------------------------------------------
{
   size_t result = address() ;
   return( result  ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: has_same_address( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   return( other == this ) ;
}




//----------------------------------------------------------------------
int
MAC_Object:: address_comparison( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   int result = 0 ;
   if( this > other )
   {
      result = 1 ;
   }
   else if( this < other )
   {
      result = -1 ;
   }
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Error::object()->raise_not_implemented( this, "save_state" ) ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Error::object()->raise_not_implemented( this, "restore_state" ) ;
   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: update_for_restore_state( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: update_for_restore_state" ) ;
   MAC_CHECK_INV( invariant() ) ;
   update() ;
}




//----------------------------------------------------------------------
void
MAC_Object:: display_info( std::ostream& os, size_t indent_width ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: display_info" ) ;
   std::string const s( indent_width, ' ' ) ;
   os << s << "type    : " << type_name() << std::endl ;
   if( owner()==0 )
   {
      os << s << "owner   : nil" << std::endl ;
   }
   else
   {
      os << s << "owner   : " << std::endl ;
      os << s << "   type    : " << owner()->type_name() << std::endl ;
      os << s << "   address : "
         << std::setiosflags( std::ios::hex )
         << std::setiosflags( std::ios::showbase )
         << owner()
         << std::resetiosflags( std::ios::showbase )
         << std::resetiosflags( std::ios::hex )
         << std::endl ;
   }
   if( POSSESSIONS!=0 && POSSESSIONS->count()!=0 )
   {
      os << s << "owned objects : " << std::endl ;
      for( size_t i=0 ; i<POSSESSIONS->count() ; i++ )
      {
         POSSESSIONS->at(i)->display_info(os,indent_width+3) ;
      }
   }
   os << s << "address : "
      << std::setiosflags( std::ios::hex )
      << std::setiosflags( std::ios::showbase )
      << this
      << std::resetiosflags( std::ios::showbase )
      << std::resetiosflags( std::ios::hex )
      << std::endl ;
}




//----------------------------------------------------------------------
void
MAC_Object:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: print" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
int
MAC_Object:: GetNumberOf_MAC_objects( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: GetNumberOf_MAC_objects" ) ;
   return( (int)ALLOCATED ) ;
}




//----------------------------------------------------------------------
std::ostream&  
MAC_Object:: TraceRemainingObjects( std::ostream& out  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Object:: TraceRemainingObjects" ) ;
   out << std::endl ;
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Objects ) )
   {
      for( map<MAC_Object const*,size_t>::const_iterator it=
              _MAC_OBJECT_LIST->begin() ;
           it!=_MAC_OBJECT_LIST->end() ;
           ++it ) 
      {
         out << "[" << (*it).second << "] :" << std::endl ;
         (*it).first->display_info( out, 3 ) ;
         (*it).first->print( out, 3 ) ;
         out << std::endl ;
      }
   }
   else
   {
      out << "To enable traceRemainingObjects facilities do run \n"
          << " application with -Cobjects command line option.\n" ;
   }
   return out ;
}




//----------------------------------------------------------------------
void
MAC_Object:: catch_object( MAC_Object const* obj )
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( MAC_Assertion::is_handling_check( MAC_Assertion::Objects )) ;
   catched_object = obj ;
}




//----------------------------------------------------------------------
void
MAC_Object:: catch_object_by_rank( size_t rank )
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( MAC_Assertion::is_handling_check( MAC_Assertion::Objects )) ;
   catched_object_rank = rank ;
}




//----------------------------------------------------------------------
void
MAC_Object:: start_trace_allocating( void )
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( !trace_allocating() ) ;
   trace_allocation = true ;
   MAC_CHECK_POST( trace_allocating() ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: trace_allocating( void )
//----------------------------------------------------------------------
{
   return( trace_allocation ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: stop_trace_allocating( void )
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( trace_allocating() ) ;
   trace_allocation = false ;
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Objects ) )
   {
      _MAC_OBJECT_LIST2->clear() ;
   }
   MAC_CHECK_POST( !trace_allocating() ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: trace_not_destroyed_object( std::ostream& out )
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( trace_allocating() ) ;
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Objects ) )
   {
      for( map<MAC_Object const*,size_t>::const_iterator it=
                       _MAC_OBJECT_LIST2->begin() ;
           it!=_MAC_OBJECT_LIST2->end() ;
           ++it )
      {
         out << "[" << it->second << "] :" << std::endl ;
         it->first->display_info( out, 3 ) ;
         it->first->print( out, 3 ) ;
         out << std::endl ;
      }
   }
}




//----------------------------------------------------------------------
bool
MAC_Object:: invariant( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: create_clone_POST( MAC_Object const* result,
                                MAC_Object const* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: save_state_PRE( MAC_ObjectWriter const* writer ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( writer != 0 ) ;
   MAC_ASSERT( writer->has_an_opened_cycle() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: save_state_POST( MAC_ObjectWriter const* writer ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( writer->has_an_opened_cycle() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: restore_state_PRE( MAC_ObjectReader const* reader ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( reader != 0 ) ;
   MAC_ASSERT( reader->positioned_in_a_valid_cycle() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: restore_state_POST( MAC_ObjectReader const* reader ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( reader->positioned_in_a_valid_cycle() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: is_equal_PRE( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( comparable( other ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: is_equal_POST( bool result, MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( IMPLIES( result, hash_code()==other->hash_code() ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: three_way_comparison_PRE( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( comparable( other ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Object:: three_way_comparison_POST( int result, 
                                        MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( EQUIVALENT( result==0, is_equal(other) ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
void
MAC_Object:: insert_possession( MAC_Object* obj )
//----------------------------------------------------------------------
{
   MAC_CHECK( obj!=0 ) ;

   if( POSSESSIONS == 0 )
   {
       POSSESSIONS = MAC_ListIdentity:: create( 0 ) ;  
   }
   POSSESSIONS->prepend( obj ) ;
   obj->MY_OWNER = this ;

   MAC_CHECK_POST( obj->has_same_address( POSSESSIONS->item( obj ) ) ) ;
}

