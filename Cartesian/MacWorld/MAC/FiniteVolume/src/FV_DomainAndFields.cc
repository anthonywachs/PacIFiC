#include <FV_DomainAndFields.hh>
#include <FV_PostProcessingWriter.hh>
#include <FV_DomainBuilder.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
#include <MAC.hh>
#include <MAC_Bool.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_Int.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_Root.hh>
#include <MAC_String.hh>
#include <MAC_assertions.hh>
#include <fstream>
#include <iostream>
#include <sstream>
using std::cout ;
using std::endl ;
using std::ostringstream ; 

struct FV_DomainAndFields_ERROR
{
   static void n1( std::string const& model_name,
	std::string const& name_of_new_field ) ;   
} ;


//-------------------------------------------------------------------------
FV_DomainAndFields*
FV_DomainAndFields:: create( MAC_Object* a_owner,
                              MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   FV_DomainAndFields* result = new FV_DomainAndFields( a_owner, exp, 0 ) ;
  
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( !result->is_distributed() ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_DomainAndFields*
FV_DomainAndFields:: create( MAC_Object* a_owner,
                              MAC_ModuleExplorer const* exp,
                              MAC_Communicator const* com )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( com != 0 ) ;

   FV_DomainAndFields* result = new FV_DomainAndFields( a_owner, exp, com ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->is_distributed() ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_DomainAndFields:: FV_DomainAndFields( MAC_Object* a_owner,
                                           MAC_ModuleExplorer const* exp,
                                           MAC_Communicator const* com )
//-------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , EXP( exp->create_clone( this ) )
   , BUILDER_FV( 0 )
   , COM( com )
   , BUILDER_WRITER( 0 )
{
   MAC_LABEL( "FV_DomainAndFields:: FV_DomainAndFields" ) ;
   
   std::string my_name = "unnamed" ;
   if( exp->has_entry( "name" ) ) my_name = exp->string_data( "name" ) ;

   BUILDER_FV = FV_DomainBuilder::create( this, exp, my_name ) ;

   FVFIELDS = *(BUILDER_FV->list_primary_fields()) ;

   MAC_CHECK_INV( invariant() ) ;
}




//-------------------------------------------------------------------------
FV_DomainAndFields:: ~FV_DomainAndFields( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: ~FV_DomainAndFields" ) ;
}




//-------------------------------------------------------------------------
std::string const&
FV_DomainAndFields:: name( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: name" ) ;

   std::string const& result = BUILDER_FV->name() ;
   return( result ) ;
}




//-------------------------------------------------------------------------
void
FV_DomainAndFields:: duplicate_field(
		std::string const& model_name,
		std::string const& name_of_new_field ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: duplicate_field(with copy)" ) ;

   list< FV_DiscreteField* >::const_iterator il ;
   FV_DiscreteField const* modelField = NULL ;
   bool found = false ;

   for (il=FVFIELDS.begin(); il!=FVFIELDS.end() && !found; il++)
     if ( (*il)->name() == model_name )
     {
       modelField = *il ;
       found = true ; 
     } 
     
   if ( !found ) 
     FV_DomainAndFields_ERROR:: n1( model_name, name_of_new_field ) ;
   else
   {
     FV_DiscreteField* newField = modelField->create_clone( BUILDER_FV,
      		name_of_new_field ) ; 
     FVFIELDS.push_back( newField ) ;
   }           
}




//----------------------------------------------------------------------
void
FV_DomainAndFields:: append_field( FV_DiscreteField* newField )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: append_field" ) ;
   MAC_CHECK_PRE( newField != 0 ) ;

   FVFIELDS.push_back( newField ) ;
}




//-------------------------------------------------------------------------
bool
FV_DomainAndFields:: is_distributed( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: is_distributed" ) ;

   return( COM != 0 ) ;
}




//-------------------------------------------------------------------------
size_t
FV_DomainAndFields:: nb_space_dimensions( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: nb_space_dimensions" ) ;
   
   size_t result = BUILDER_FV->nb_space_dimensions() ;
   return( result ) ;
}



//-------------------------------------------------------------------------
FV_PostProcessingWriter*
FV_DomainAndFields:: post_processing_writer( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "PDE_DomainAndFields:: post_processing_writer" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   list< FV_DiscreteField*>::const_iterator il;
   bool isBinary = false;
   std::string postprocessingwriter = "" ;
   std::string writing_mode = "" ;

   if ( BUILDER_WRITER == 0 )
   {
     MAC_ModuleExplorer* sexp =
                           EXP->create_subexplorer( 0, "FV_ResultSaver" ) ;
   
     if ( sexp->has_entry( "postprocessingwriter" ) )
 	postprocessingwriter = sexp->string_data( "postprocessingwriter" );
     else postprocessingwriter = "paraview"; //by default 
     if ( postprocessingwriter != "paraview" 
     		&& postprocessingwriter != "matlab" )
     {
       string error_message = "  - paraview\n  - matlab\n ";
       MAC_Error::object()->raise_bad_data_value(sexp, 
	"postprocessingwriter", 
	error_message );
     }
     
     if ( sexp->has_entry( "writing_mode" ) )
       writing_mode = sexp->string_data( "writing_mode" );
     if ( writing_mode == "binary" ) isBinary = true ;
     if ((writing_mode != "text") && (writing_mode != "binary") )
     {
       string error_message = "  - text\n  - binary\n ";
       MAC_Error::object()->raise_bad_data_value(sexp, 
	"writing_mode", 
	error_message );
     }

 			
     sexp->start_module_iterator() ;
     for( ; sexp->is_valid_module() ; sexp->go_next_module() )
     {
       MAC_ModuleExplorer* ssexp = sexp->create_subexplorer( 0 ) ;
       std::string const& save_name = ssexp->string_data( "entry_name" ) ;
       std::string const& location = ssexp->string_data( "where_to_save" ) ;
       ssexp->test_data_in( "where_to_save", "at_vertices," "at_cell_centers" );
       if( ssexp->has_entry( "field" ) )
       {
         std::string const& field_name = ssexp->string_data( "field" ) ;
	 FV_DiscreteField* field = discrete_field( field_name ) ;
	 if (field)
	 {
	   BUILDER_FVFIELDS.push_back( field ) ;
	   field->set_postprocessing_options( location, save_name ) ;
	 }
	 else
	   MAC_Error::object()->raise_plain( 
	 		"Bad syntax in module " + ssexp->name() ) ;
       }
       else
       {
         MAC_Error::object()->raise_plain( 
	 		"Bad syntax in module " + ssexp->name() ) ;
       }
       ssexp->destroy() ;

     }

     list< FV_DiscreteField const* > constFields ;
     for ( il=BUILDER_FVFIELDS.begin(); il!=BUILDER_FVFIELDS.end(); il++)
       constFields.push_back(*il);
      
     BUILDER_WRITER =
       FV_PostProcessingWriter::make(
     			const_cast<FV_DomainAndFields*>( this ), 
     			postprocessingwriter,
     			sexp,
     			COM,
     			constFields,
     			primary_grid(),
     			isBinary ) ;


     sexp->destroy() ;
   }
   FV_PostProcessingWriter* result = BUILDER_WRITER ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_DomainAndFields:: append_post_processing_field( 
	FV_DiscreteField* newField, std::string const& location, 
	std::string const& save_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: append_post_processing_field" ) ;
   MAC_CHECK_PRE( newField != 0 ) ;
   MAC_CHECK_PRE( location == "at_vertices" || location == "at_cell_centers" ) ;
   MAC_CHECK_PRE( !save_name.empty() ) ;  

   newField->set_postprocessing_options( location, save_name ) ;
   BUILDER_FVFIELDS.push_back( newField ) ;
}




//-------------------------------------------------------------------------
void
FV_DomainAndFields:: save_state( MAC_ObjectWriter* writer ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "FV_DomainAndFields" ) ;

   // Discretization status :
   writer->add_entry( "nb_space_dimensions",
                      MAC_Int::create( 0, nb_space_dimensions() ) ) ;

   // Mesh translation features storing
   BUILDER_FV->save_state( writer ) ;

   // Fields storing :
   for ( list< FV_DiscreteField*>::const_iterator
   		il=FVFIELDS.begin(); il!=FVFIELDS.end(); il++)
     (*il)->save_state( writer ) ;

   writer->finalize_object() ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}



//-------------------------------------------------------------------------
void
FV_DomainAndFields:: restore_state( MAC_ObjectReader* reader )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "FV_DomainAndFields" ) ;

   // Does some checks
   MAC_ASSERT( reader->data_of_entry( "nb_space_dimensions" )->to_int()==
                                            (int) nb_space_dimensions() ) ;

   // Restore mesh translation features
   BUILDER_FV->restore_state( reader ) ;

   // Restore fields :
   std::string stored_name; 
   bool found = false;   
   // If a field is defined in the input data file but was not stored in the 
   // restart files (because of a change of application), it is initialized
   // to the values in the input data file
   while( reader->next_object_class_in_current_module() == "FV_DiscreteField" )
   {
     stored_name = reader->next_object_name_in_current_module();
     found = false;
     for ( list< FV_DiscreteField*>::const_iterator
   		il=FVFIELDS.begin(); il!=FVFIELDS.end() && !found; il++)
     {
       if ( (*il)->name() == stored_name )
       { 
	 found = true;	 
	 (*il)->restore_state( reader ) ;
         (*il)->restore_translated_field_mesh() ;
         
	 // Boundary conditions: we reset them in case they have been modified
	 // by hand in the restart header file
	 (*il)->set_imposed_DOF_values();
	 (*il)->set_neumann_DOF_values();	 
       }
     }
     
     // If the field was stored but is not used in the new application
     // We still need to do a "fake" read for the reader to proceed properly
     // i.e. essentially to be able to move to the next object
     if ( !found )
     {
       FV_DiscreteField:: read_state_nonrestored( reader );
       if ( MAC_Exec::communicator()->rank() == 0 )
         MAC::out() << "      Field \"" << stored_name << "\" stored but not"
	 	<< " reloaded" << endl;
     }  
   }
   reader->end_object_retrieval() ;

   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//-------------------------------------------------------------------------
bool
FV_DomainAndFields:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: invariant" ) ;
   MAC_ASSERT( MAC_Object::invariant() ) ;
   MAC_ASSERT( nb_space_dimensions()==2 ||
               nb_space_dimensions()==3 ) ;
   return( true ) ;
}




//-------------------------------------------------------------------------
void
FV_DomainAndFields:: print_grid( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: print_grid" ) ;
      
   BUILDER_FV->primary_grid()->print( os, indent_width ) ;
}




//-------------------------------------------------------------------------
FV_DiscreteField*
FV_DomainAndFields:: discrete_field( std::string const& field_name ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: discrete_field" ) ;
   
   FV_DiscreteField* result = 0;

   list<FV_DiscreteField*>::const_iterator il;
   bool found = false;
   
   for (il=FVFIELDS.begin(); il!=FVFIELDS.end() && !found; il++)
     if ( (*il)->name() == field_name ) 
     {
       result = *il;
       found = true;
     }

   MAC_CHECK_POST( result != 0 ) ;   
   return( result ) ;
}




//-------------------------------------------------------------------------
bool
FV_DomainAndFields:: has_discrete_field( std::string const& field_name ) 
	const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: has_discrete_field" ) ;
   
   bool result = false;
   list<FV_DiscreteField*>::const_iterator il;
   
   for (il=FVFIELDS.begin(); il!=FVFIELDS.end() && !result; il++)
     if ( (*il)->name() == field_name ) result = true;

   return( result ) ;
}




//----------------------------------------------------------------------
FV_Mesh const*
FV_DomainAndFields:: primary_grid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainAndFields:: primary_grid" ) ;
   
   return( BUILDER_FV->primary_grid() ) ;
}




//internal--------------------------------------------------------------
void
FV_DomainAndFields_ERROR:: n1( std::string const& model_name,
	std::string const& name_of_new_field )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Cannot duplicate field \"" << model_name << "\" in "
   	<< "field \"" << name_of_new_field << "\" because field \""
	<< model_name << "\" does not exist in FV_DomainAndFields !! " 
	<< endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
