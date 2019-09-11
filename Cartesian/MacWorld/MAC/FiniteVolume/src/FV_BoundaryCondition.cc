#include <FV_BoundaryCondition.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <FV_DomainBuilder.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Variable.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Communicator.hh>
#include <stringVector.hh>
#include <boolVector.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <math.h>
using std::cout ; 
using std::endl ;
using std::ostringstream ; 
using std::setw ;


struct FV_BoundaryCondition_ERROR
{
   static void n1( string const& field_name, size_t component, 
   	string const& color_name ) ;   
} ;


//----------------------------------------------------------------------
FV_BoundaryCondition*
FV_BoundaryCondition:: create( MAC_Object* a_owner, 
 		size_t number_of_components, size_t color )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: create" ) ;
   MAC_CHECK_PRE( number_of_components > 0 );

   FV_BoundaryCondition* result = new FV_BoundaryCondition( 
   	a_owner, number_of_components, color ) ;

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_BoundaryCondition:: FV_BoundaryCondition( MAC_Object* a_owner, 
 		size_t number_of_components, size_t color )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , BC_TYPE( 0 )
   , UNKNOWN_ON_BC( 0 )   
   , BC_VALUES_ALL_COMPS( 0 )
   , NB_COMPS( number_of_components )
   , BC_COLOR( color )
   , CTX( 0 )
   , COORDS( 0 )
   , READ( 0 )
{
   MAC_LABEL( "FV_BoundaryCondition:: FV_BoundaryCondition" ) ;
 
   list<FV_TRIPLET> empty_list;
   nodes_localStructuredNumbering.reserve( NB_COMPS );
   for (size_t i=0;i<NB_COMPS;++i)
     nodes_localStructuredNumbering.push_back( empty_list );
     
   string Neumann = "neumann";
   BC_TYPE = new stringVector( NB_COMPS );
   for (size_t i=0;i<NB_COMPS;++i)
     (*BC_TYPE)(i) = Neumann ;
     
   READ = new boolVector( NB_COMPS, false );
   UNKNOWN_ON_BC = new boolVector( NB_COMPS, false ); 
   
   SHIFT_UPDATE_BC_VALUES.reserve( NB_COMPS );
   list< pair< double,FV_SHIFT_TRIPLET> > mst ;    
   for (size_t i=0;i<NB_COMPS;++i)
     SHIFT_UPDATE_BC_VALUES.push_back( mst ) ;
               
   MAC_DataWithContext const* pc = 0;  
   BC_VALUES_PER_COMP.reserve( NB_COMPS );
   for (size_t i=0;i<NB_COMPS;++i)
     BC_VALUES_PER_COMP.push_back( pc );
     
   MAC_ContextSimple* c = MAC_ContextSimple::create( this ) ;
   COORDS = MAC_DoubleVector::create( c, doubleVector(0) ) ;
   c->extend( MAC_Variable::object( "DV_X" ), COORDS ) ;
   CTX = c ;               
}




//----------------------------------------------------------------------
FV_BoundaryCondition:: ~FV_BoundaryCondition( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: ~FV_BoundaryCondition" ) ;
   
   vector< list<FV_TRIPLET> >::iterator iv ;
   for (iv=nodes_localStructuredNumbering.begin();
     iv!=nodes_localStructuredNumbering.end();iv++) iv->clear();
   vector< list< pair< double,FV_SHIFT_TRIPLET> > >::iterator ivv ;
   for (ivv=SHIFT_UPDATE_BC_VALUES.begin();ivv!=SHIFT_UPDATE_BC_VALUES.end();
     ivv++) ivv->clear();
   if ( BC_TYPE ) delete BC_TYPE;
   if ( READ ) delete READ ;
   if ( UNKNOWN_ON_BC ) delete UNKNOWN_ON_BC ;   
}




//----------------------------------------------------------------------
FV_BoundaryCondition*
FV_BoundaryCondition:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: create_clone" ) ;

   FV_BoundaryCondition* result = new FV_BoundaryCondition( 
   	a_owner, NB_COMPS, BC_COLOR ) ;

   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     for (list<FV_TRIPLET>::const_iterator il=
     	nodes_localStructuredNumbering[comp].begin();il!=
     	nodes_localStructuredNumbering[comp].end();il++)
       result->nodes_localStructuredNumbering[comp].push_back( *il );
       
     (*result->BC_TYPE)(comp) = (*BC_TYPE)(comp) ;
     (*result->UNKNOWN_ON_BC)(comp) = (*UNKNOWN_ON_BC)(comp) ;

     for (list< pair<double,FV_SHIFT_TRIPLET> >::const_iterator il=
     	SHIFT_UPDATE_BC_VALUES[comp].begin();il!=
     	SHIFT_UPDATE_BC_VALUES[comp].end();il++)
       result->SHIFT_UPDATE_BC_VALUES[comp].push_back( *il );
       
     if ( BC_VALUES_PER_COMP[comp] != 0 )
       result->BC_VALUES_PER_COMP[comp] = 
       	BC_VALUES_PER_COMP[comp]->create_clone( result->CTX ) ;

     (*result->READ)(comp) = (*READ)(comp) ;     
   } 
   
   if ( BC_VALUES_ALL_COMPS != 0 )
     result->BC_VALUES_ALL_COMPS = 
       BC_VALUES_ALL_COMPS->create_clone( result->CTX ) ;          

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: add_MacTriplet( size_t component, 
	FV_TRIPLET const& mt )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: add_MacTriplet" ) ;

   if ( component < NB_COMPS )
     nodes_localStructuredNumbering[component].push_back( mt );
   else
     for (size_t i=0;i<NB_COMPS;++i)
       nodes_localStructuredNumbering[i].push_back( mt );    
    
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: add_MacTriplet( size_t component, 
	size_t i, size_t j, size_t k )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: add_MacTriplet" ) ;
   
   FV_TRIPLET mt;
   mt.i = i;
   mt.j = j;
   mt.k = k;   

   if ( component < NB_COMPS )
     nodes_localStructuredNumbering[component].push_back( mt );
   else
     for (size_t l=0;l<NB_COMPS;++l)
       nodes_localStructuredNumbering[l].push_back( mt );
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_BC_type( string const& type_of_boundary_condition,
	size_t component )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_BC_type" ) ;
   MAC_CHECK_PRE( type_of_boundary_condition == "dirichlet"
   	|| type_of_boundary_condition == "neumann" ) ;
   
   if ( component < NB_COMPS )
     (*BC_TYPE)(component) = type_of_boundary_condition ;
   else
     for (size_t i=0;i<NB_COMPS;++i)
       (*BC_TYPE)(i) = type_of_boundary_condition ;   
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: print( std::ostream& os, size_t indent_width,
      	size_t nb_space_dimensions ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: print" ) ;
   
   std::string space( indent_width, ' ' ) ;
   list<FV_TRIPLET>::const_iterator il ; 
   list< pair< double,FV_SHIFT_TRIPLET> >::const_iterator im ; 
    
   os << space << "Color = " << FV_DomainBuilder::get_color_name(BC_COLOR)
   	<< endl;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     os << space << "Component " << comp << endl;
     os << space << "Type of BC = " << (*BC_TYPE)(comp) << endl;
     os << space << "Unknown on BC = " << 
     	((*UNKNOWN_ON_BC)(comp) == 0 ? "no" : "yes") << endl;
     os << space << "Shift triplets for neumann BC values =";	
     if ( SHIFT_UPDATE_BC_VALUES[comp].empty() )
       os << space << " none";
     else
       for (im=SHIFT_UPDATE_BC_VALUES[comp].begin();
     		im!=SHIFT_UPDATE_BC_VALUES[comp].end();im++) 
       {
         os << space << " (w=" << im->first << "," << im->second.i << ","
       	  	<< im->second.j;
         if ( nb_space_dimensions == 3 ) os << "," << im->second.k;
         os << ")";
       }
     os << endl;              
     os << space << "Local node numbering : ";
     if ( nodes_localStructuredNumbering[comp].empty() )
         os << space << "none";
     else
       for (il=nodes_localStructuredNumbering[comp].begin();
       		il!=nodes_localStructuredNumbering[comp].end();il++)
       {
	 os << endl << space << "  " << il->i << " " << il->j;
	 if ( nb_space_dimensions == 3 ) os << " " << il->k;
       }
     os << endl;        
   }
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: read_dirichlet_BC( MAC_ModuleExplorer* sse,
	size_t nb_space_dimensions, string const& field_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: read_dirichlet_BC" ) ;
   MAC_CHECK_PRE( sse != 0 ) ;

   size_t comp = sse->int_data( "component" ) ;

   if ( comp < NB_COMPS ) 
   {  
     if ( (*READ)(comp) )
       FV_BoundaryCondition_ERROR:: n1( field_name, comp,  
     	FV_DomainBuilder::get_color_name(BC_COLOR) ) ;
   }
   else
     for (size_t j=0;j<NB_COMPS;++j)
       if ( (*READ)(j) )
         FV_BoundaryCondition_ERROR:: n1( field_name, j, 
     		FV_DomainBuilder::get_color_name(BC_COLOR) ) ;
   
   set_BC_type( "dirichlet", comp ); 

   if ( !sse->has_entry( "value" ) )
   {
     MAC_Error::object()->raise_missing_keyword( sse, "value" ) ; 
   }
   MAC_DataWithContext const* DEFAULT_FORMULA = sse->abstract_data( 
   	CTX, "value", CTX ) ;
   if( DEFAULT_FORMULA->data_type() != MAC_Data::DoubleVector )
   {
      MAC_Error::object()->raise_bad_data_type( sse, "value", 
		MAC_Data::DoubleVector ) ;
   }
   doubleVector dof_coordinates( nb_space_dimensions, 0. );
   COORDS->set( dof_coordinates ) ;
   doubleVector vv = DEFAULT_FORMULA->to_double_vector( CTX );
   if ( ( comp < NB_COMPS && vv.size() != 1 ) 
   	|| ( comp >= NB_COMPS && vv.size() != NB_COMPS ) )
   {
      MAC_Error::object()->raise_data_error( sse, "value", 
		"number of values does not match components number" ) ;
   }
   
   if ( vv.size() == 1 )
   {
     BC_VALUES_PER_COMP[comp] = DEFAULT_FORMULA;
     (*READ)(comp) = true ;
   }
   else
   {
     BC_VALUES_ALL_COMPS = DEFAULT_FORMULA;
     for (size_t j=0;j<NB_COMPS;++j) (*READ)(j) = true ; 
   }      
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_imposed_DOF_values( FV_DiscreteField* ff )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_imposed_DOF_values" ) ;
   
   size_t field_storage_depth = ff->storage_depth(), formula_index ;
   size_t dim = ff->primary_grid()->nb_space_dimensions();
   list<FV_TRIPLET>::iterator il;
   doubleVector dof_coordinates( dim, 0. );
   MAC_DataWithContext const* CURRENT_FORMULA = NULL ;

   // Up to now only dirichlet and homogeneous Neumann are handled
   // BC_TYPE(i) == "dirichlet" <=> BC_VALUES_PER_COMP(i) != 0
   // BC_VALUES_ALL_COMPS != 0 <=> BC_TYPE(i) == "dirichlet" for all i
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     if ( BC_VALUES_ALL_COMPS ) 
     {
       CURRENT_FORMULA = BC_VALUES_ALL_COMPS;
       formula_index = comp ;
     }
     else if ( BC_VALUES_PER_COMP[comp] )
     { 
       CURRENT_FORMULA = BC_VALUES_PER_COMP[comp] ;
       formula_index = 0 ;
     }       
     else CURRENT_FORMULA = NULL ;
     
     if ( CURRENT_FORMULA )
       for (il=nodes_localStructuredNumbering[comp].begin();
       	il!=nodes_localStructuredNumbering[comp].end();il++)
       {
         for (size_t j=0;j<dim;++j) 
           dof_coordinates(j) = ff->get_DOF_coordinate( il->get_index(j), 
	   	comp, j ); 
	 COORDS->set( dof_coordinates ) ;
	 doubleVector const& val = CURRENT_FORMULA->to_double_vector() ;
         for (size_t lev=0;lev<field_storage_depth;++lev)
	   ff->set_DOF_value( il->i, il->j, il->k, comp, 
	   	lev, val( formula_index ) ) ;
       }
   }   
}




//----------------------------------------------------------------------
double 
FV_BoundaryCondition:: get_imposed_DOF_values( size_t & component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: get_imposed_DOF_values" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   MAC_CHECK_PRE( BC_TYPE != 0 ) ;
   
   MAC_DataWithContext const* CURRENT_FORMULA = NULL ;
   
   if ( BC_VALUES_ALL_COMPS ) 
     CURRENT_FORMULA = BC_VALUES_ALL_COMPS;
   else if ( BC_VALUES_PER_COMP[component] )
     CURRENT_FORMULA = BC_VALUES_PER_COMP[component] ;
   else CURRENT_FORMULA = NULL ;
   
   doubleVector val = CURRENT_FORMULA->to_double_vector() ;
   
     
   return( val(component) );
   
}





//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_free_DOF_values( 
	FV_DiscreteField* ff, size_t level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_free_DOF_values" ) ;

   list<FV_TRIPLET>::iterator il;
   double neighboring_DOF_value = 0. ; 
   list< pair< double,FV_SHIFT_TRIPLET> >::iterator is;
     
   for (size_t comp=0;comp<NB_COMPS;++comp)
     if ( (*BC_TYPE)(comp) == "neumann" && !(*UNKNOWN_ON_BC)(comp) )     
     {
       for (il=nodes_localStructuredNumbering[comp].begin();
       	il!=nodes_localStructuredNumbering[comp].end();il++)
       {
         neighboring_DOF_value = 0. ;
	 for (is=SHIFT_UPDATE_BC_VALUES[comp].begin();
	 	is!=SHIFT_UPDATE_BC_VALUES[comp].end();is++)
	   neighboring_DOF_value += is->first * ff->DOF_value( 
	 	il->i+is->second.i, 
		il->j+is->second.j, 
		il->k+is->second.k, 
	 	comp, level ) ;
		
	 ff->set_DOF_value( il->i, il->j, il->k, comp, level, 
	 	neighboring_DOF_value ) ;
       }     
     }      
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_free_DOF_values( 
	FV_DiscreteField* ff, size_t component, size_t level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_free_DOF_values" ) ;

   list<FV_TRIPLET>::iterator il;
   double neighboring_DOF_value = 0. ; 
   list< pair< double,FV_SHIFT_TRIPLET> >::iterator is;
     
   if ( (*BC_TYPE)(component) == "neumann" && !(*UNKNOWN_ON_BC)(component) )     
   {
     for (il=nodes_localStructuredNumbering[component].begin();
       	il!=nodes_localStructuredNumbering[component].end();il++)
     {
       neighboring_DOF_value = 0. ;
       for (is=SHIFT_UPDATE_BC_VALUES[component].begin();
	 	is!=SHIFT_UPDATE_BC_VALUES[component].end();is++)
	 neighboring_DOF_value += is->first * ff->DOF_value( 
	 	il->i+is->second.i, 
		il->j+is->second.j, 
		il->k+is->second.k, 
	 	component, level ) ;
		
        ff->set_DOF_value( il->i, il->j, il->k, component, level, 
	 	neighboring_DOF_value ) ;
     }     
   }      
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_periodic( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_periodic" ) ;
   MAC_CHECK_PRE( BC_TYPE != 0 ) ;      

   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     (*BC_TYPE)(comp) = "periodic" ;
     (*READ)(comp) = true ;
   }
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_none( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_none" ) ;
   MAC_CHECK_PRE( BC_TYPE != 0 ) ;      

   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     (*BC_TYPE)(comp) = "none" ;
     (*READ)(comp) = true ;
   }
}




//----------------------------------------------------------------------
bool
FV_BoundaryCondition:: is_dirichlet( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: is_dirichlet" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   MAC_CHECK_PRE( BC_TYPE != 0 ) ;      

   return ( (*BC_TYPE)(component) == "dirichlet" ) ;
}




//----------------------------------------------------------------------
bool
FV_BoundaryCondition:: is_neumann( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: is_neumann" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   MAC_CHECK_PRE( BC_TYPE != 0 ) ;      

   return ( (*BC_TYPE)(component) == "neumann" ) ;
}




//----------------------------------------------------------------------
bool
FV_BoundaryCondition:: is_periodic( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: is_periodic" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   MAC_CHECK_PRE( BC_TYPE != 0 ) ;      

   return ( (*BC_TYPE)(component) == "periodic" ) ;
}




//----------------------------------------------------------------------
bool
FV_BoundaryCondition:: is_none( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: is_none" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   MAC_CHECK_PRE( BC_TYPE != 0 ) ;      

   return ( (*BC_TYPE)(component) == "none" ) ;
}




//----------------------------------------------------------------------
size_t
FV_BoundaryCondition:: get_color_ID( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: get_color_ID" ) ;     

   return ( BC_COLOR ) ;
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_unknown_on_BC( size_t component )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_unknown_on_BC" ) ;     
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   MAC_CHECK_PRE( UNKNOWN_ON_BC != 0 ) ;   

   (*UNKNOWN_ON_BC)(component) = true ;      
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_shift_MacTriplet( size_t component,
	int i, int j, int k, double coefficient )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_shift_MacTriplet" ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   
   if ( SHIFT_UPDATE_BC_VALUES[component].empty() )
   {
     FV_SHIFT_TRIPLET mst;
     mst.i = 0 ;
     mst.j = 0 ;
     mst.k = 0 ;
     pair< double,FV_SHIFT_TRIPLET> ppp( 0., mst ) ;
     SHIFT_UPDATE_BC_VALUES[component].push_back( ppp ) ;
   }  
   
   SHIFT_UPDATE_BC_VALUES[component].front().first = coefficient ;  
   SHIFT_UPDATE_BC_VALUES[component].front().second.i = i ;
   SHIFT_UPDATE_BC_VALUES[component].front().second.j = j ;
   SHIFT_UPDATE_BC_VALUES[component].front().second.k = k ;        
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: set_shift_MacTriplet( size_t component,
	FV_SHIFT_TRIPLET const* mst, double coefficient )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: set_shift_MacTriplet" ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   
   if ( SHIFT_UPDATE_BC_VALUES[component].empty() )
   {
     pair< double,FV_SHIFT_TRIPLET> ppp( 0., *mst ) ;
     SHIFT_UPDATE_BC_VALUES[component].push_back( ppp ) ;
   }  
   
   SHIFT_UPDATE_BC_VALUES[component].front().first = coefficient ;  
   SHIFT_UPDATE_BC_VALUES[component].front().second.i = mst->i ;
   SHIFT_UPDATE_BC_VALUES[component].front().second.j = mst->j ;
   SHIFT_UPDATE_BC_VALUES[component].front().second.k = mst->k ;        
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: add_shift_MacTriplet( size_t component,
	int i, int j, int k, double coefficient )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: add_shift_MacTriplet" ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   
   FV_SHIFT_TRIPLET mst;
   mst.i = i ;
   mst.j = j ;
   mst.k = k ;
   pair< double,FV_SHIFT_TRIPLET> ppp( coefficient, mst ) ;
   SHIFT_UPDATE_BC_VALUES[component].push_back( ppp ) ;   
}




//----------------------------------------------------------------------
void
FV_BoundaryCondition:: add_shift_MacTriplet( size_t component,
	FV_SHIFT_TRIPLET const* mst, double coefficient )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: add_shift_MacTriplet" ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ) ; 
   
   pair< double,FV_SHIFT_TRIPLET> ppp( coefficient, *mst ) ;
   SHIFT_UPDATE_BC_VALUES[component].push_back( ppp ) ;   
}




//----------------------------------------------------------------------
bool
FV_BoundaryCondition:: has_unknown( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: has_unknown" ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ) ;
   MAC_CHECK_PRE( UNKNOWN_ON_BC != 0 ) ;    

   return ( (*UNKNOWN_ON_BC)(component) ) ;
}




//----------------------------------------------------------------------
bool
FV_BoundaryCondition:: has_DOF_on_proc( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_BoundaryCondition:: has_DOF_on_proc" ) ; 

   bool result = false ;
   vector< list<FV_TRIPLET> >::const_iterator iv ;
   
   for (iv=nodes_localStructuredNumbering.begin();
	iv!=nodes_localStructuredNumbering.end() && !result;iv++)
     if ( !iv->empty() ) result = true ;	   

   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_BoundaryCondition:: compute_boundary_cell_centered_DOF_integral( 
	FV_DiscreteField const* ff, size_t component, size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_BoundaryCondition:: compute_boundary_cell_centered_DOF_integral" ) ;
   MAC_CHECK( ff->discretization_type() == "centered"
   	|| ff->discretization_type() == "staggered" ); 
   MAC_CHECK( FV_DomainBuilder::is_main_color( BC_COLOR ) );
   MAC_CHECK( IMPLIES( ff->discretization_type() == "staggered", 
	IMPLIES( BC_COLOR == FV_BC_LEFT || BC_COLOR == FV_BC_RIGHT, 
   		component == 0 ) ) ) ;
   MAC_CHECK( IMPLIES( ff->discretization_type() == "staggered", 
	IMPLIES( BC_COLOR == FV_BC_BOTTOM || BC_COLOR == FV_BC_TOP, 
   		component == 1 ) ) ) ;
   MAC_CHECK( IMPLIES( ff->discretization_type() == "staggered", 
	IMPLIES( BC_COLOR == FV_BC_BEHIND || BC_COLOR == FV_BC_FRONT, 
   		component == 2 ) ) ) ;				   	

   // DOFs centered in the boundary cells of a face always have their
   // coordinate in the normal direction to the face equal to that of the face
   // Hence the integral is simply the sum of the products of the boundary cell
   // measure and the field value at the boundary cell center, restricted to DOF
   // that are located in the domain on the current processor  
   
   list<FV_TRIPLET>::const_iterator il;
   double result = 0. ;
   size_t normal_dir = FV_DomainBuilder:: normal_direction_to_main_color( 
   	BC_COLOR ).first ;
   bool in = true ;
   FV_Mesh const* primary_mesh = ff->primary_grid() ;
   size_t dimens = primary_mesh->nb_space_dimensions() ;   

   for (il=nodes_localStructuredNumbering[component].begin();
   	il!=nodes_localStructuredNumbering[component].end();il++)
   {
     in = true ;
     for (size_t i=0;i<dimens && in;++i)
       if ( i != normal_dir )
         in = primary_mesh->is_in_domain_on_current_processor( 
	 	ff->get_DOF_coordinate( il->get_index(i), component, i ), 
		i, 0. ) ;
		
     if ( in )
       result += ff->DOF_value( il->i, il->j, il->k, component, level )
       		* ff->get_face_perp_to_direction_measure( il->i, il->j, il->k, 
		component, normal_dir ) ;       
   }
   
   double collective_result = MAC_Exec::communicator()->sum( result ) ;
   
   return ( collective_result ) ;
}




//----------------------------------------------------------------------
double
FV_BoundaryCondition:: compute_boundary_mean_normal_derivative( 
	FV_DiscreteField const* ff, size_t component, size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_BoundaryCondition:: compute_boundary_cell_centered_DOF_integral" ) ;
   MAC_CHECK( FV_DomainBuilder::is_main_color( BC_COLOR ) );
   MAC_CHECK( (*BC_TYPE)(component) == "dirichlet" );
   
   list<FV_TRIPLET>::const_iterator il = 
   	nodes_localStructuredNumbering[component].begin();
   double result = 0., surface = 0., face_area = 0. ;
   size_t normal_dir = FV_DomainBuilder:: normal_direction_to_main_color( 
   	BC_COLOR ).first ;
   bool in = true ;
   FV_Mesh const* primary_mesh = ff->primary_grid() ;
   size_t dimens = primary_mesh->nb_space_dimensions() ;
   FV_SHIFT_TRIPLET normalshift ;
   normalshift.i = 0 ;
   normalshift.j = 0 ;   
   normalshift.k = 0 ;   
   
   switch( BC_COLOR )
   {
     case FV_BC_LEFT : normalshift.i = 1 ; break ; 
     case FV_BC_RIGHT : normalshift.i = -1 ; break ;      
     case FV_BC_BOTTOM : normalshift.j = 1 ; break ;      
     case FV_BC_TOP : normalshift.j = -1 ; break ;      
     case FV_BC_BEHIND : normalshift.k = 1 ; break ;      
     case FV_BC_FRONT : normalshift.k = -1 ; break ;           
   }
      
   // Determine signed local normal grid spacing
   double hnormal = 0.;
   if ( !nodes_localStructuredNumbering[component].empty() ) 
     hnormal = ff->get_DOF_coordinate( il->get_index(normal_dir), 
   	component, normal_dir ) - ff->get_DOF_coordinate( 
	il->get_index(normal_dir) + normalshift.get_index(normal_dir), 
	component, normal_dir ) ;

   // Compute mean normal derivative
   for ( ;il!=nodes_localStructuredNumbering[component].end();il++)
   {
     in = true ;
     for (size_t i=0;i<dimens && in;++i)
       if ( i != normal_dir )
         in = primary_mesh->is_in_domain_on_current_processor( 
	 	ff->get_DOF_coordinate( il->get_index(i), component, i ), 
		i, 0. ) ;
		
     if ( in )
     {
       face_area = ff->get_face_perp_to_direction_measure( il->i + normalshift.i,
			 il->j + normalshift.j, il->k + normalshift.k, 
			component, normal_dir ) ;
       result += ( ff->DOF_value( il->i, il->j, il->k, component, level )
       		-  ff->DOF_value( il->i + normalshift.i, il->j + normalshift.j, 
			il->k + normalshift.k, component, level ) )
       		* face_area / hnormal ; 
       surface +=  face_area ;			
     }     
   }
   
   double collective_result = MAC_Exec::communicator()->sum( result ) ;
   double collective_surface = MAC_Exec::communicator()->sum( surface ) ; 
   
   collective_result /= collective_surface ;  
   
   return ( collective_result ) ;
}




//internal--------------------------------------------------------------
void
FV_BoundaryCondition_ERROR:: n1( string const& field_name, 
	size_t component, string const& color_name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Boundary condition of color \"" << color_name 
   	<< "\" for component " << component << " of field \"" << field_name
    	<< "\" is specified twice or is both Dirichlet and periodic; "
	<< "check the boundary conditions" << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
