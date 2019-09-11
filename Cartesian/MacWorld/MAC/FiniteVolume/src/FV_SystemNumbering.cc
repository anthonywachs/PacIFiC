#include <FV_SystemNumbering.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ModuleExplorer.hh>
#include <LA_Scatter.hh>
#include <LA_Vector.hh>
#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
using std::endl ;
using std::string ;
using std::ostringstream ;


//----------------------------------------------------------------------
FV_SystemNumbering*
FV_SystemNumbering::create( MAC_Object* a_owner,
	FV_DiscreteField const* the_field,
	size_t a_verbose_level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_SystemNumbering::create" ) ;
   MAC_CHECK_PRE( the_field != 0 ) ;

   FV_SystemNumbering* result = new FV_SystemNumbering( a_owner,
   	the_field, a_verbose_level ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_SystemNumbering:: FV_SystemNumbering( MAC_Object* a_owner,
	FV_DiscreteField const* the_field,
	size_t a_verbose_level )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , FIELD( the_field )
   , IDX_LOCS( 0 )
   , IDX_GLOBS( 0 )
   , OK_SCATTER( false )
   , SCATTER( 0 )
   , NB_HANDLED_UNK( 0 )
   , MAT_SIZE( MAC::bad_index() )
   , VERB( a_verbose_level )
{
   MAC_LABEL( "FV_SystemNumbering:: FV_SystemNumbering" ) ;

   MAT_SIZE = FIELD->nb_global_unknowns() ;
   NB_HANDLED_UNK = FIELD->nb_local_unknowns_handled_by_proc() ;
   
   size_t nlocunk = FIELD->nb_local_unknowns() ;
   IDX_LOCS = new size_t_vector( nlocunk ) ;
   IDX_GLOBS = new size_t_vector( nlocunk ) ;   

   FIELD->build_system_numbering( IDX_LOCS, IDX_GLOBS ) ;     
}




//----------------------------------------------------------------------
FV_SystemNumbering:: ~FV_SystemNumbering( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_SystemNumbering:: ~FV_SystemNumbering" ) ;

   delete IDX_LOCS ;
   delete IDX_GLOBS ;   

   MAC_CHECK_INV( invariant() ) ;
}





//----------------------------------------------------------------------
FV_DiscreteField const*
FV_SystemNumbering:: field( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_SystemNumbering:: field" ) ;

   return ( FIELD ) ;
}




//----------------------------------------------------------------------
void
FV_SystemNumbering:: define_scatter( LA_Vector const* vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_SystemNumbering:: define_scatter" ) ;
   MAC_CHECK_PRE( vec->nb_rows() == nb_global_unknowns() ) ;

   SCATTER = vec->create_scatter( this, *IDX_GLOBS, *IDX_LOCS );
   OK_SCATTER = true ;

   MAC_CHECK_POST( scatter_is_defined() ) ;
   MAC_CHECK_POST( SCATTER->implementation() == vec->implementation() ) ;
}




//----------------------------------------------------------------------
bool
FV_SystemNumbering:: scatter_is_defined( void ) const
//----------------------------------------------------------------------
{
   return( OK_SCATTER ) ;
}




//----------------------------------------------------------------------
LA_Scatter const*
FV_SystemNumbering:: scatter( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_SystemNumbering:: scatter" ) ;
   MAC_CHECK_PRE( scatter_is_defined() ) ;

   LA_Scatter const* result = SCATTER ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_SystemNumbering:: nb_global_unknowns( void ) const
//----------------------------------------------------------------------
{
   return( MAT_SIZE );
}




//----------------------------------------------------------------------
size_t
FV_SystemNumbering:: nb_unknowns_on_current_process( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_SystemNumbering:: nb_unknowns_on_current_process" ) ;

   return( NB_HANDLED_UNK ) ;
}




//-----------------------------------------------------------------------------
void
FV_SystemNumbering:: print( std::ostream& os, size_t indent_width ) const
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SystemNumbering:: print" ) ;

   std::string space( indent_width, ' ' ) ;

   os << space << "nb global unknowns = " << MAT_SIZE << endl ;
   size_t nlocunk = IDX_LOCS->size();
   os << space << "nb local unknowns = " << nlocunk << endl ;
   os << space << "# unkNumLoc unkNumGlob" << endl;
   for (size_t i=0;i<nlocunk;++i)
     os << space << (*IDX_LOCS)(i) << " " << (*IDX_GLOBS)(i) << endl;  
}
