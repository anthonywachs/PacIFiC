#include <MAC_Context.hh>

#include <MAC_assertions.hh>
#include <MAC_List.hh>
#include <MAC_Variable.hh>

#include <iostream>
#include <string>

//----------------------------------------------------------------------
MAC_Context:: MAC_Context( MAC_Object* a_owner ) 
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , OBSERVERS( MAC_List::create( this ) )
{
}

//----------------------------------------------------------------------
MAC_Context:: ~MAC_Context( void ) 
//----------------------------------------------------------------------
{
   OBSERVERS = 0 ;
}

//----------------------------------------------------------------------
void
MAC_Context:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Context:: print" ) ;
   print( os, indent_width, (MAC_Context const*) 0 ) ;
}

//----------------------------------------------------------------------
void
MAC_Context:: print( std::ostream& os, size_t indent_width,
                     MAC_Context const* ctx ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Context:: print" ) ;
   
   std::string bl( indent_width, ' ' ) ;
   for( size_t i=0 ; i<nb_variables() ; ++i  )
   {
      MAC_Variable const* var = variable(i) ;
      MAC_Data const* dat = value( var ) ;
      os << bl <<  "$" << var->name() ;
      if( ctx == 0 || !ctx->has_variable( var ) )
      {
         os << " = " ;
      }
      else
      {
         os << " == " ;
      }
      dat->print( os, 0 ) ;
      os << std::endl ;
   }   
}

//----------------------------------------------------------------------
void
MAC_Context:: attach_observer( MAC_Context* observer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Context:: attach_observer" ) ;
   MAC_ASSERT( observer != this ) ;
   OBSERVERS->extend( observer ) ;
}

//----------------------------------------------------------------------
void
MAC_Context:: detach_observer( MAC_Context* observer ) const
//----------------------------------------------------------------------
{
   OBSERVERS->remove( observer ) ;
}


//----------------------------------------------------------------------
void
MAC_Context:: update_observers( void ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<OBSERVERS->index_limit() ; i++ )
   {
      MAC_Context* ct = static_cast<MAC_Context*>( OBSERVERS->at(i) ) ;
      ct->update() ;
   }
}

//----------------------------------------------------------------------
void
MAC_Context:: notify_observers_of_my_destruction( void ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<OBSERVERS->index_limit() ; i++ )
   {
      MAC_Context* ct = static_cast<MAC_Context* >(OBSERVERS->at(i)) ;
      ct->update_for_destruction_of( this ) ;
   }
}

//----------------------------------------------------------------------
bool
MAC_Context:: variable_PRE( size_t i ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( i < nb_variables() ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
MAC_Context:: variable_POST( MAC_Variable const* result ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( has_variable( result ) ) ;
   MAC_ASSERT( result!=0 ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
MAC_Context:: has_variable_PRE( MAC_Variable const* var ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( var!=0 ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
MAC_Context:: value_PRE( MAC_Variable const* var ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( var!=0 ) ;
   MAC_ASSERT( has_variable( var ) ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
MAC_Context:: value_POST( MAC_Data* result,
                          MAC_Variable const* var ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result!=0 ) ;
   MAC_ASSERT( var->data_type()==result->data_type() ) ;
   return true ;
}

