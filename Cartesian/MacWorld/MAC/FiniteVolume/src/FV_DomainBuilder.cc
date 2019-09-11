#include <FV_DomainBuilder.hh>
#include <FV_Mesh.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DiscreteField.hh>
#include <FV.hh>
#include <MAC.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_IndexSet.hh>
#include <MAC_List.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Root.hh>
#include <MAC_Vector.hh>
#include <MAC_VectorIterator.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_ObjectRegister.hh>
#include <doubleArray2D.hh>
#include <stringVector.hh>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
using std::cout ; 
using std::endl ;
using std::ostringstream ; 
using std::setw ;
using std::map ;

struct FV_DomainBuilder_ERROR
{
   static void n1( std::string const& name ) ;   
   static void n2( std::string const& name ) ;
   static void n3( std::string const& name ) ;   
   static void n4( size_t const& colorID ) ;
} ;

stringVector* FV_DomainBuilder::allowedDOFcolors = NULL;
stringVector* FV_DomainBuilder::allowedDOFstatus = NULL;
size_t FV_DomainBuilder::DIM = 2;
list< pair< size_t, list<size_t> > >* FV_DomainBuilder::subMainColorIds
	= NULL ;
list< pair< size_t, FV_SHIFT_TRIPLET > >* 
	FV_DomainBuilder::mainColorShiftMacTriplets = NULL ;	


//----------------------------------------------------------------------
std::map< std::string, FV_DomainBuilder* >&
FV_DomainBuilder::INSTANCES( void )
//----------------------------------------------------------------------
{
   static std::map< std::string, FV_DomainBuilder* > result ;
   return result ;
}




//----------------------------------------------------------------------
FV_DomainBuilder*
FV_DomainBuilder:: object( std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: object" ) ;

   FV_DomainBuilder* result = 0 ;

   map<string,FV_DomainBuilder*>::const_iterator it =
                                                  INSTANCES().find( a_name ) ;
   if( it == INSTANCES().end() )
   {
      std::ostringstream ss ;
      ss << "domain \"" << a_name << "\" does not exists" ;
      MAC_Error::object()->raise_plain( ss.str() ) ;
   }
   else
   {
      result = (*it).second ;
   }
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_DomainBuilder*
FV_DomainBuilder:: create( MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp,
                            std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: create" ) ;

   if( INSTANCES().count( a_name ) != 0 )
   {
      std::ostringstream ss ;
      ss << "domain \"" << a_name << "\" already exists" ;
      MAC_Error::object()->raise_plain( ss.str() ) ;
   }

   FV_DomainBuilder* result = new FV_DomainBuilder( a_owner, exp, a_name ) ;

   INSTANCES()[ a_name ] = result ;

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_DomainBuilder:: FV_DomainBuilder( MAC_Object* a_owner,
                                       MAC_ModuleExplorer const* exp,
                                       std::string const& a_name )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , EXP( exp->create_clone( this ) )
   , NAME( a_name )
   , VERB( exp->int_data( "verbose_level" ) )
   , GRID( 0 )
{
   MAC_LABEL( "FV_DomainBuilder:: FV_DomainBuilder" ) ;

   exp->test_data_in( "nb_space_dimensions", "2,3" ) ;
   DIM = exp->int_data( "nb_space_dimensions" ) ;

   build_allowed_DOFstatus();

   build_allowed_DOFcolors();
   
   build_special_colors();

   build_all();   
}




//----------------------------------------------------------------------
FV_DomainBuilder:: ~FV_DomainBuilder( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: ~FV_DomainBuilder" ) ;

   if ( allowedDOFcolors ) 
   {
     delete allowedDOFcolors;
     allowedDOFcolors = NULL; 
   }
   
   if ( subMainColorIds ) 
   {
     delete subMainColorIds;
     subMainColorIds = NULL; 
   } 
   
   if ( mainColorShiftMacTriplets ) 
   {
     delete mainColorShiftMacTriplets;
     mainColorShiftMacTriplets = NULL; 
   }      

   map<string,FV_DomainBuilder*>::iterator it = INSTANCES().find( NAME ) ;
   MAC_CHECK( it != INSTANCES().end() ) ;
   INSTANCES().erase( it ) ;
}




//----------------------------------------------------------------------
std::string const&
FV_DomainBuilder:: name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: name" ) ;

   return( NAME ) ;
}




//----------------------------------------------------------------------
bool
FV_DomainBuilder:: has_explorer( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: has_explorer" ) ;

   return( EXP->has_module( path_and_name ) ) ;
}




//----------------------------------------------------------------------
MAC_ModuleExplorer*
FV_DomainBuilder:: create_explorer( MAC_Object* a_owner,
                                     std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: create_explorer" ) ;

   return( EXP->create_subexplorer( a_owner, path_and_name ) ) ;
}




//----------------------------------------------------------------------
size_t
FV_DomainBuilder:: nb_space_dimensions( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: nb_space_dimensions" ) ;

   return( DIM ) ;
}




//----------------------------------------------------------------------
FV_Mesh const*
FV_DomainBuilder:: primary_grid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: primary_grid" ) ;
   
   if( GRID == 0 )
   {
      MAC_Error::object()->raise_plain( "Finite Volume Grid not available" ) ;
   }
   FV_Mesh* result = GRID ;
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner()==this ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_DomainBuilder:: build_all( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: build_all" ) ;

   if( VERB!=0 && name()!="unnamed" )
      FV::out() << "*** Domain : \"" << name() << "\"" << endl << endl ;

   MAC_Exec::communicator()->barrier();
   if ( MAC_Exec::communicator()->rank() == 0 ) FV::out() << endl;
   
   build_grid();   
   build_discrete_fields();   
}




//----------------------------------------------------------------------
void
FV_DomainBuilder:: build_grid( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: build_grid" ) ;

   if( VERB!=0 ) FV::out() << "*** Building grid" << endl ;
   GRID = FV_Mesh::create( this, EXP, DIM ) ;
}




//----------------------------------------------------------------------
void
FV_DomainBuilder:: build_discrete_fields( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "PDE_DomainBuilder:: build_discrete_fields" ) ;
   
   if( !EXP->has_module( "interior_fields" ) )
   {
      FV_DomainBuilder_ERROR::n2( EXP->name() ) ;
   }

   if ( MAC_Exec::communicator()->rank() == 0 ) 
     FV::out() << "*** Building FV discrete fields" << endl;
   MAC_Exec::communicator()->barrier();
   
   // fields defined on the interior of the domain
   // --------------------------------------------
   if( EXP->has_module( "interior_fields" ) )
   {  
      bool b_dualBC = false;
      MAC_ModuleExplorer* ee = EXP->create_subexplorer( 0, "interior_fields" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {  
         MAC_ModuleExplorer const* se = ee->create_subexplorer( 0 ) ;
	 if ( se->has_entry( "dualBC" ) ) b_dualBC = se->bool_data( "dualBC" );
	 size_t storage_depth = b_dualBC == true ? 
	   ( se->int_data( "storage_depth" ) + 1 ) :
	   se->int_data( "storage_depth");
         FV_DiscreteField* ff = 
             FV_DiscreteField::create( this, GRID,
		se->string_data( "name" ),
	     	se->string_data( "discretization" ),
		se->int_data( "nb_components" ),
		storage_depth) ;

	 ff->build_BCs( se, this ) ;
	 ff->build_field_numbering() ;	 	 
	 ff->set_imposed_DOF_values() ;	 
	 ff->initialize_DOFs( se );
	 ff->set_BC_values_modif_status( b_dualBC );
	 b_dualBC = false;
	 	 
         PRIMARY_FVFIELDS.push_back( ff );
	 
	 ff->out_endOfBuilding( FV::out(), 6, 
	 	MAC_Exec::communicator()->rank() ) ;
         
         se->destroy() ; se = 0 ;
      }
      ee->destroy() ; ee=0 ;
   }  
}




//----------------------------------------------------------------------
void
FV_DomainBuilder:: build_special_colors( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: build_special_colors" ) ;

   // macro-colors
   // ------------
   if( EXP->has_module( "macro_colors" ) )
   {
      MAC_ModuleExplorer* se = EXP->create_subexplorer( 0, "macro_colors" ) ;
      se->start_entry_iterator() ;
      for( ; se->is_valid_entry() ; se->go_next_entry() )
      {
         string const& nn = se->keyword() ;
         stringVector const& name_list = se->data(this)->to_string_vector(0) ;
         for (size_t i=0;i<name_list.size();++i)
	   if ( !does_primary_color_exist( name_list(i) ) )
             MAC_Error::object()->raise_data_error( se, nn, 
		"contains "+name_list(i)+", a primary color that does not exist"
		) ; 
	 pair< string, stringVector > pp( nn, name_list ) ;
	 FVRO_COLORS.push_back( pp );	   	   
      }
      se->destroy() ;
   }
}




//----------------------------------------------------------------------
list< FV_DiscreteField* > const*
FV_DomainBuilder:: list_primary_fields( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: list_primary_fields" ) ;
   
   MAC_CHECK_PRE( !PRIMARY_FVFIELDS.empty() ) ;

   return( &PRIMARY_FVFIELDS );
}




//----------------------------------------------------------------------
void
FV_DomainBuilder:: build_allowed_DOFstatus( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: build_allowed_DOFstatus" ) ;
   
   if ( !allowedDOFstatus )
   {
     allowedDOFstatus = new stringVector( 6 ) ;
     (*allowedDOFstatus)(FV_DOF_ONPROC) = "onproc";
     (*allowedDOFstatus)(FV_DOF_BUFFERZONE) = "bufferzone";
     (*allowedDOFstatus)(FV_DOF_HALOZONE) = "halozone";
     (*allowedDOFstatus)(FV_DOF_PERIODIC_BUFFERZONE) = "periodic_bufferzone";
     (*allowedDOFstatus)(FV_DOF_PERIODIC_HALOZONE) = "periodic_halozone";
     (*allowedDOFstatus)(FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE) = 
     	"bufferzone_periodic_bufferzone";
   }
}




//----------------------------------------------------------------------
void
FV_DomainBuilder:: build_allowed_DOFcolors( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: build_allowed_DOFcolors" ) ;
   
   if ( !allowedDOFcolors )
   {
     allowedDOFcolors = new stringVector( DIM == 2 ? 10 : 28 ) ;
     (*allowedDOFcolors)(FV_BC_INTERIOR) = "interior";
     (*allowedDOFcolors)(FV_BC_PERIODIC) = "periodic";
     (*allowedDOFcolors)(FV_BC_LEFT) = "left";      
     (*allowedDOFcolors)(FV_BC_RIGHT) = "right";
     (*allowedDOFcolors)(FV_BC_BOTTOM) = "bottom";      
     (*allowedDOFcolors)(FV_BC_TOP) = "top";
     (*allowedDOFcolors)(FV_BC_BOTTOM_LEFT) = "bottom_left";      
     (*allowedDOFcolors)(FV_BC_BOTTOM_RIGHT) = "bottom_right";
     (*allowedDOFcolors)(FV_BC_TOP_LEFT) = "top_left";      
     (*allowedDOFcolors)(FV_BC_TOP_RIGHT) = "top_right";
     if ( DIM == 3 )
     {
       (*allowedDOFcolors)(FV_BC_BEHIND) = "behind"; 
       (*allowedDOFcolors)(FV_BC_FRONT) = "front";                   
       (*allowedDOFcolors)(FV_BC_BEHIND_LEFT) = "behind_left";
       (*allowedDOFcolors)(FV_BC_BEHIND_RIGHT) = "behind_right";
       (*allowedDOFcolors)(FV_BC_FRONT_LEFT) = "front_left";
       (*allowedDOFcolors)(FV_BC_FRONT_RIGHT) = "front_right"; 
       (*allowedDOFcolors)(FV_BC_BEHIND_BOTTOM) = "behind_bottom";
       (*allowedDOFcolors)(FV_BC_BEHIND_TOP) = "behind_top";
       (*allowedDOFcolors)(FV_BC_FRONT_BOTTOM) = "front_bottom";
       (*allowedDOFcolors)(FV_BC_FRONT_TOP) = "front_top";        
       (*allowedDOFcolors)(FV_BC_BEHIND_BOTTOM_LEFT) = "behind_bottom_left";
       (*allowedDOFcolors)(FV_BC_BEHIND_BOTTOM_RIGHT) = "behind_bottom_right";
       (*allowedDOFcolors)(FV_BC_BEHIND_TOP_LEFT) = "behind_top_left";
       (*allowedDOFcolors)(FV_BC_BEHIND_TOP_RIGHT) = "behind_top_right";
       (*allowedDOFcolors)(FV_BC_FRONT_BOTTOM_LEFT) = "front_bottom_left";
       (*allowedDOFcolors)(FV_BC_FRONT_BOTTOM_RIGHT) = "front_bottom_right";
       (*allowedDOFcolors)(FV_BC_FRONT_TOP_LEFT) = "front_top_left";
       (*allowedDOFcolors)(FV_BC_FRONT_TOP_RIGHT) = "front_top_right";
     }  
   }
   
   if ( !subMainColorIds )
   {
     subMainColorIds = new list< pair< size_t, list<size_t> > > ;
     pair< size_t, list<size_t> > pp;

     pp.first = get_color_number("bottom_left");
     pp.second.push_back(get_color_number("bottom"));
     pp.second.push_back(get_color_number("left"));
     subMainColorIds->push_back(pp);
     pp.second.clear();
     
     pp.first = get_color_number("bottom_right");
     pp.second.push_back(get_color_number("bottom"));
     pp.second.push_back(get_color_number("right"));
     subMainColorIds->push_back(pp);
     pp.second.clear();

     pp.first = get_color_number("top_left");
     pp.second.push_back(get_color_number("top"));
     pp.second.push_back(get_color_number("left"));
     subMainColorIds->push_back(pp);
     pp.second.clear();     
     
     pp.first = get_color_number("top_right");
     pp.second.push_back(get_color_number("top"));
     pp.second.push_back(get_color_number("right"));
     subMainColorIds->push_back(pp);
     pp.second.clear();     

     if ( DIM == 3 )
     {     
       pp.first = get_color_number("behind_left");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("left"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
     
       pp.first = get_color_number("behind_right");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("right"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
      
       pp.first = get_color_number("front_left");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("left"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
     
       pp.first = get_color_number("front_right");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("right"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
     
       pp.first = get_color_number("behind_bottom");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("bottom"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
     
       pp.first = get_color_number("behind_top");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("top"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
     
       pp.first = get_color_number("front_bottom");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("bottom"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
     
       pp.first = get_color_number("front_top");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("top"));
       subMainColorIds->push_back(pp);
       pp.second.clear(); 
     
       pp.first = get_color_number("behind_bottom_left");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("bottom"));
       pp.second.push_back(get_color_number("left"));     
       subMainColorIds->push_back(pp);
       pp.second.clear();
     
       pp.first = get_color_number("behind_bottom_right");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("bottom"));
       pp.second.push_back(get_color_number("right"));     
       subMainColorIds->push_back(pp);
       pp.second.clear();
     
       pp.first = get_color_number("behind_top_left");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("top"));
       pp.second.push_back(get_color_number("left"));     
       subMainColorIds->push_back(pp);
       pp.second.clear();
     
       pp.first = get_color_number("behind_top_right");
       pp.second.push_back(get_color_number("behind"));
       pp.second.push_back(get_color_number("top"));
       pp.second.push_back(get_color_number("right"));     
       subMainColorIds->push_back(pp);
       pp.second.clear();
     
       pp.first = get_color_number("front_bottom_left");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("bottom"));
       pp.second.push_back(get_color_number("left"));     
       subMainColorIds->push_back(pp);
       pp.second.clear();
     
       pp.first = get_color_number("front_bottom_right");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("bottom"));
       pp.second.push_back(get_color_number("right"));     
       subMainColorIds->push_back(pp);
       pp.second.clear();
     
       pp.first = get_color_number("front_top_left");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("top"));
       pp.second.push_back(get_color_number("left"));     
       subMainColorIds->push_back(pp);
       pp.second.clear();
     
       pp.first = get_color_number("front_top_right");
       pp.second.push_back(get_color_number("front"));
       pp.second.push_back(get_color_number("top"));
       pp.second.push_back(get_color_number("right")); 
       subMainColorIds->push_back(pp);
       pp.second.clear();
     }
   } 
   
   if ( !mainColorShiftMacTriplets )
   {
     mainColorShiftMacTriplets = new list< pair< size_t, FV_SHIFT_TRIPLET > > ;
     FV_SHIFT_TRIPLET mst ;
     pair< size_t, FV_SHIFT_TRIPLET > pp(0,mst) ;
     
     pp.first = get_color_number("left");
     pp.second.i = 1 ;     
     pp.second.j = 0 ;     
     pp.second.k = 0 ;
     mainColorShiftMacTriplets->push_back(pp);
     
     pp.first = get_color_number("right");
     pp.second.i = -1 ;     
     pp.second.j = 0 ;     
     pp.second.k = 0 ;
     mainColorShiftMacTriplets->push_back(pp);
     
     pp.first = get_color_number("bottom");
     pp.second.i = 0 ;     
     pp.second.j = 1 ;     
     pp.second.k = 0 ;
     mainColorShiftMacTriplets->push_back(pp);
     
     pp.first = get_color_number("top");
     pp.second.i = 0 ;     
     pp.second.j = -1 ;     
     pp.second.k = 0 ;
     mainColorShiftMacTriplets->push_back(pp);
     
     if ( DIM == 3 )
     {
       pp.first = get_color_number("front");
       pp.second.i = 0 ;     
       pp.second.j = 0 ;     
       pp.second.k = -1 ;
       mainColorShiftMacTriplets->push_back(pp);
     
       pp.first = get_color_number("behind");
       pp.second.i = 0 ;     
       pp.second.j = 0 ;     
       pp.second.k = 1 ;
       mainColorShiftMacTriplets->push_back(pp); 
     }                   
   }  
}




//----------------------------------------------------------------------
size_t
FV_DomainBuilder:: get_color_number( string const& color_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: get_color_number" ) ;
   MAC_CHECK_PRE( allowedDOFcolors != 0 ) ;   

   bool found = false;
   size_t ncolors = allowedDOFcolors->size(),result=0;
   
   for (size_t i=0;i<ncolors && !found;++i)
     if ( (*allowedDOFcolors)(i) == color_name ) 
     {
       result = i;
       found = true;
     }
   
   if (!found)
     FV_DomainBuilder_ERROR:: n1( color_name );
     
   return ( result );     
}




//----------------------------------------------------------------------
string
FV_DomainBuilder:: get_color_name( size_t const& color_id )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: get_color_name" ) ;
   MAC_CHECK_PRE( allowedDOFcolors != 0 ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 2, color_id < 10 ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 3, color_id < 28 ) ) ; 
   
   string result = (*allowedDOFcolors)(color_id) ;
   
   return ( result ) ;
}     




//----------------------------------------------------------------------
string
FV_DomainBuilder:: get_status_name( size_t const& status_id )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: get_status_name" ) ;
   MAC_CHECK_PRE( allowedDOFstatus != 0 ) ;
   MAC_CHECK_PRE( status_id < 6 ) ; 
   
   string result = (*allowedDOFstatus)(status_id) ;
   
   return ( result ) ;
} 




//----------------------------------------------------------------------
bool
FV_DomainBuilder:: does_primary_color_exist( string const& color_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: does_primary_color_exist" ) ;
   MAC_CHECK_PRE( allowedDOFcolors != 0 ) ;   

   bool found = false;
   size_t ncolors = allowedDOFcolors->size();
   
   for (size_t i=0;i<ncolors && !found;++i)
     if ( (*allowedDOFcolors)(i) == color_name ) found = true;
     
   return ( found );     
}




//----------------------------------------------------------------------
bool
FV_DomainBuilder:: is_main_color( size_t const& colorID )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: is_main_color" ) ;
   MAC_CHECK_PRE( allowedDOFcolors != 0 ) ;   

   bool result = ( ( colorID >= FV_BC_LEFT && colorID <= FV_BC_TOP )
   	|| colorID == FV_BC_BEHIND || colorID == FV_BC_FRONT ) ;
     
   return ( result );     
}




//----------------------------------------------------------------------
list<size_t> const* 
FV_DomainBuilder:: get_sub_main_color_ids( size_t const& color_id )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: get_sub_main_color_ids" ) ;
   MAC_CHECK_PRE( subMainColorIds != 0 ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 2, color_id >= FV_BC_BOTTOM_LEFT 
   		&& color_id <= FV_BC_TOP_RIGHT ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	( color_id >= FV_BC_BOTTOM_LEFT && color_id <= FV_BC_TOP_RIGHT ) ||
	( color_id >= FV_BC_BEHIND_LEFT 
		&& color_id <= FV_BC_FRONT_TOP_RIGHT ) ) ) ; 

   list<size_t> const* result = NULL;
   list< pair< size_t, list<size_t> > >::iterator il ;
   bool found = false;   
   
   for (il=subMainColorIds->begin();il!=subMainColorIds->end() && !found;il++)
     if ( il->first == color_id )
     {
       found = true ;
       result = &il->second ;
     }   
   
   return ( result ) ;
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET const*
FV_DomainBuilder:: get_shift_MacTriplet( size_t const& color_id )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: get_shift_MacTriplet" ) ;
   MAC_CHECK_PRE( mainColorShiftMacTriplets != 0 ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 2, color_id >= FV_BC_LEFT 
   		&& color_id <= FV_BC_TOP ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	( color_id >= FV_BC_LEFT && color_id <= FV_BC_TOP ) ||
	( color_id >= FV_BC_BEHIND && color_id <= FV_BC_FRONT ) ) ) ; 

   FV_SHIFT_TRIPLET const* result = NULL;
   list< pair< size_t, FV_SHIFT_TRIPLET > >::iterator il ;
   bool found = false;   
   
   for (il=mainColorShiftMacTriplets->begin();
   	il!=mainColorShiftMacTriplets->end() && !found;il++)
     if ( il->first == color_id )
     {
       found = true ;
       result = &il->second ;
     }   
   
   return ( result ) ;
}




//----------------------------------------------------------------------
stringVector const*
FV_DomainBuilder:: get_colors_from_macro_color( string const& color_name )
	const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: get_colors_from_macro_color" ) ;
   MAC_CHECK_PRE( allowedDOFcolors != 0 ) ;   

   stringVector const* result = NULL ;
   list< pair< string, stringVector > >::const_iterator il;
   bool found = false;   
   
   for (il=FVRO_COLORS.begin();il!=FVRO_COLORS.end() && !found;il++)
     if (il->first == color_name )
     {
       found = true ;
       result = &il->second ;
     }
   
   if ( !found )
     FV_DomainBuilder_ERROR::n3( color_name ) ;
   
   return ( result ) ;  
}




//----------------------------------------------------------------------
pair<size_t,string>
FV_DomainBuilder:: normal_direction_to_main_color( size_t const& colorID )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: normal_direction_to_main_color" ) ;
   MAC_CHECK_PRE( allowedDOFcolors != 0 ) ;   

   if ( !FV_DomainBuilder:: is_main_color( colorID ) )
     FV_DomainBuilder_ERROR:: n4( colorID ) ;
     
   pair<size_t,string> result(0,"min");
   switch(colorID)
   {
     case FV_BC_LEFT:
       result.first = 0;
       result.second = "min";
       break;
       
     case FV_BC_RIGHT:
       result.first = 0;
       result.second = "max";
       break;
       
     case FV_BC_BOTTOM:
       result.first = 1;
       result.second="min";
       break;
       
     case FV_BC_TOP:
       result.first = 1;
       result.second = "max";
       break;
       
     case FV_BC_BEHIND:
       result.first = 2;
       result.second = "min";
       break;
       
     case FV_BC_FRONT:
       result.first = 2;
       result.second = "max";
       break;                          
   }
     
   return ( result );     
}



//----------------------------------------------------------------------
void
FV_DomainBuilder:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "FV_DomainBuilder" ) ;
   
   GRID->save_state( writer ) ;
   
   writer->finalize_object() ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
FV_DomainBuilder:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DomainBuilder:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "FV_DomainBuilder" ) ;

   GRID->restore_state( reader ) ;   

   reader->end_object_retrieval() ;

   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//internal--------------------------------------------------------------
void
FV_DomainBuilder_ERROR:: n1( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Wrong color name \"" << name << "\" in "
   	<< "FV_DomainBuilder:: get_color_number" << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_DomainBuilder_ERROR:: n2( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "module \"" << name << "\" :" << endl
        << "   there should be at least a submodule called" << endl
        << "   \"interior_fields\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_DomainBuilder_ERROR:: n3( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "macro color \"" << name << "\" :" << endl
        << "   does not exist in "
	<< "FV_DomainBuilder:: get_colors_from_macro_color" << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_DomainBuilder_ERROR:: n4( size_t const& colorID )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "macro color \"" << colorID << "\" :" << endl
        << "   is not a macro color in "
	<< "FV_DomainBuilder:: normal_direction_to_main_color" << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
