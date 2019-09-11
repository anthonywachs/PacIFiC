#include <FV_DiscreteField_Centered.hh>
#include <FV_Mesh.hh>
#include <FV_BoundaryCondition.hh>
#include <FV_DomainBuilder.hh>
#include <FV.hh>
#include <MAC_Exec.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleArray2D.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Error.hh>
#include <MAC_String.hh>
#include <MAC_Variable.hh>
#include <LA_Vector.hh>
#include <LA_Matrix.hh>
#include <doubleArray2D.hh>
#include <intArray3D.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>

using std::cout ; 
using std::endl ;
using std::string ; 
using std::ostringstream ;


FV_DiscreteField_Centered const* 
FV_DiscreteField_Centered::PROTOTYPE = new FV_DiscreteField_Centered( 
	"centered" ) ;


//----------------------------------------------------------------------
FV_DiscreteField_Centered:: FV_DiscreteField_Centered( 
	std::string const& a_type )
//----------------------------------------------------------------------	
   : FV_DiscreteField( a_type ) 
{
   MAC_LABEL( "FV_DiscreteField_Centered:: FV_DiscreteField_Centered") ;
}




//----------------------------------------------------------------------
FV_DiscreteField_Centered:: FV_DiscreteField_Centered( 
	MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type )
//----------------------------------------------------------------------	
   : FV_DiscreteField( a_owner, a_primary_mesh, a_name, a_type ) 
{
   MAC_LABEL( "FV_DiscreteField_Centered:: FV_DiscreteField_Centered") ;
}




//----------------------------------------------------------------------
FV_DiscreteField*
FV_DiscreteField_Centered:: create_replica( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type,
	size_t a_nb_components,
	size_t a_depth ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: create_replica" ) ;

   FV_DiscreteField* result = new FV_DiscreteField_Centered( a_owner,
	a_primary_mesh, a_name, a_type, a_nb_components, a_depth ) ;
    
   return( result ) ;
}




//----------------------------------------------------------------------
FV_DiscreteField*
FV_DiscreteField_Centered:: create_clone_replica( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: create_clone_replica" ) ;

   FV_DiscreteField* result = new FV_DiscreteField_Centered( a_owner, 
   	a_primary_mesh, a_name, a_type ) ;
    
   return( result ) ;
}




//----------------------------------------------------------------------
FV_DiscreteField_Centered:: FV_DiscreteField_Centered( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type,
	size_t a_nb_components,
	size_t a_depth )
//----------------------------------------------------------------------
    : FV_DiscreteField( a_owner, a_primary_mesh, a_name, a_type, 
    	a_nb_components, a_depth )
{
   MAC_LABEL( "FV_DiscreteField_Centered:: FV_DiscreteField_Centered") ;

   vector< boolVector > const* primary_mesh_handles_node = 
   	PRIMARY_GRID->get_nodes_owner();
   size_t security_bandwidth = PRIMARY_GRID->get_security_bandwidth();

   // Allocate index vectors
   size_t_vector work_size_t( DIM, 0 ) ;
   global_max_index = new vector<size_t_vector>( 1, work_size_t ) ;
   local_min_index_in_global = new vector<size_t_vector>( 1, work_size_t ) ;
   local_max_index_in_global = new vector<size_t_vector>( 1, work_size_t ) ;
   local_dof_number = new vector<size_t_vector>(1, work_size_t) ;
   min_index_unknown_handled_by_proc = new vector<size_t_vector>( 1, 
   	work_size_t ) ;
   for (size_t i=0;i<DIM;++i) 
     (*min_index_unknown_handled_by_proc)[0](i) = 1000000;
   max_index_unknown_handled_by_proc = new vector<size_t_vector>( 1, 
   	work_size_t ) ;
   min_index_unknown_on_proc = new vector<size_t_vector>( 1, 
   	work_size_t ) ;
   for (size_t i=0;i<DIM;++i) 
     (*min_index_unknown_on_proc)[0](i) = 1000000;
   max_index_unknown_on_proc = new vector<size_t_vector>( 1, 
   	work_size_t ) ;	 
   
   // Build global coordinates
   vector< doubleVector > const* primary_mesh_global_main_coordinates = 
	PRIMARY_GRID->get_global_main_coordinates();
   size_t_vector np(DIM,0);
   for (size_t i=0;i<DIM;++i) 
   {
     np(i) = (*primary_mesh_global_main_coordinates)[i].size(); 
     if ( (*PERIODIC)(i) ) (*global_max_index)[0](i) = np(i) - 2 ;
     else (*global_max_index)[0](i) = np(i);
   }  
   	
   doubleVector work_double(1,0.);
   vector< doubleVector > vwork_double( DIM, work_double ) ;   
   global_main_coordinates = new vector < vector< doubleVector > >;
   global_main_coordinates->reserve(1);
   global_main_coordinates->push_back( vwork_double ) ;
   for (size_t i=0;i<DIM;++i)
     if ( (*PERIODIC)(i) ) 
     {
       (*global_main_coordinates)[0][i].re_initialize(np(i)-1,0.);
       for( size_t j=0;j<np(i)-1;++j)
        (*global_main_coordinates)[0][i](j) = 0.5 *
       		( (*primary_mesh_global_main_coordinates)[i](j)
		+ (*primary_mesh_global_main_coordinates)[i](j+1) ); 
     }
     else
     {
       (*global_main_coordinates)[0][i].re_initialize(np(i)+1,0.);
       (*global_main_coordinates)[0][i](0) = 
       	(*primary_mesh_global_main_coordinates)[i](0);
       for( size_t j=1;j<np(i);++j)
         (*global_main_coordinates)[0][i](j) = 0.5 *
       		( (*primary_mesh_global_main_coordinates)[i](j)
		+ (*primary_mesh_global_main_coordinates)[i](j-1) );
       (*global_main_coordinates)[0][i](np(i)) = 
       	(*primary_mesh_global_main_coordinates)[i](np(i)-1);
     }	   
   
   // Build local coordinates
   local_main_coordinates = new vector < vector< doubleVector > >;
   local_main_coordinates->reserve(1);
   local_main_coordinates->push_back( vwork_double ) ;      
   size_t_vector primary_mesh_local_max_index_in_global = 
   	*PRIMARY_GRID->get_local_max_index_in_global();
   size_t_vector primary_mesh_local_min_index_in_global = 
   	*PRIMARY_GRID->get_local_min_index_in_global();
   for (size_t i=0;i<DIM;++i)
   {
     if ( (*PERIODIC)(i) ) 
     {
       (*local_min_index_in_global)[0](i) =
       		primary_mesh_local_min_index_in_global(i) ;
       (*local_max_index_in_global)[0](i) = 
       	        primary_mesh_local_max_index_in_global(i) - 1 ;
     }
     else
     {     
       (*local_min_index_in_global)[0](i) = 
     		primary_mesh_local_min_index_in_global(i) == 0 ?
     		0 : primary_mesh_local_min_index_in_global(i)+1 ;
       (*local_max_index_in_global)[0](i) = 
     		primary_mesh_local_max_index_in_global(i) == np(i)-1 ?
     		np(i) : primary_mesh_local_max_index_in_global(i) ;
     }
     (*local_dof_number)[0](i) = (*local_max_index_in_global)[0](i) -
     		(*local_min_index_in_global)[0](i) + 1 ;      
     (*local_main_coordinates)[0][i].re_initialize( (*local_dof_number)[0](i), 
     	0.);
     for (size_t j=0;j<(*local_dof_number)[0](i);++j)
       (*local_main_coordinates)[0][i](j) = 
       	(*global_main_coordinates)[0][i]( (*local_min_index_in_global)[0](i) 
		+ j );
   }
      		
   // Nodes handled by current process
   vector< size_t_vector >* vwork_size_t = new vector< size_t_vector >( DIM, 
   	work_size_t ) ;   
   on_current_processor = new vector< vector< size_t_vector > >( 1 ,
   	*vwork_size_t );
   delete vwork_size_t;  	
   for (size_t i=0;i<DIM;++i)
   {
     (*on_current_processor)[0][i].re_initialize( 
     	(*local_dof_number)[0](i), 0 );
     for (size_t j=0;j<security_bandwidth;++j)
     {
       if ( !(*primary_mesh_handles_node)[i](j) )
         if ( (*local_min_index_in_global)[0](i) == 0 )
	   (*on_current_processor)[0][i](j) = 4 ;
	 else
           (*on_current_processor)[0][i](j) = 2 ;
	 
       if ( !(*primary_mesh_handles_node)[i](
		(*primary_mesh_handles_node)[i].size()-1-j) )
	 if ( (*local_max_index_in_global)[0](i) == (*global_max_index)[0](i) )
           (*on_current_processor)[0][i]((*local_dof_number)[0](i)-1-j) = 4 ;
	 else
	   (*on_current_processor)[0][i]((*local_dof_number)[0](i)-1-j) = 2 ;  
     }
   }
     
   // Bufferzone dof
   for (size_t i=0;i<DIM;++i)
   {  
     if ( (*on_current_processor)[0][i](0) == 2 && 
     	(*local_min_index_in_global)[0](i) != 0 )
       for (size_t j=0;j<security_bandwidth;++j)
         (*on_current_processor)[0][i](security_bandwidth+j) = 1 ;
     
     size_t index_max = (*on_current_processor)[0][i].size()-1;
     if ( (*on_current_processor)[0][i](index_max) == 2 && 
     	(*local_max_index_in_global)[0](i) != (*global_max_index)[0](i) )
       for (size_t j=0;j<security_bandwidth;++j)
         (*on_current_processor)[0][i](index_max-security_bandwidth-j) = 1 ;
   }
   
   // Periodic bufferzone
   for (size_t i=0;i<DIM;++i)
   {  
     if ( (*local_min_index_in_global)[0](i) == 0 && (*PERIODIC)(i) )
       for (size_t j=0;j<security_bandwidth;++j)
         if ( (*on_current_processor)[0][i](security_bandwidth+j) == 1 )
           (*on_current_processor)[0][i](security_bandwidth+j) = 5 ;
	 else
	   (*on_current_processor)[0][i](security_bandwidth+j) = 3 ;
     
     size_t index_max = (*on_current_processor)[0][i].size()-1;
     if ( (*local_max_index_in_global)[0](i) == (*global_max_index)[0](i) 
     	&& (*PERIODIC)(i) )
       for (size_t j=0;j<security_bandwidth;++j)
         if ( (*on_current_processor)[0][i](index_max-security_bandwidth-j) 
	 	== 1 )
           (*on_current_processor)[0][i](index_max-security_bandwidth-j) = 5 ;
	 else
	   (*on_current_processor)[0][i](index_max-security_bandwidth-j) = 3 ;
   }       

   // Allocate field values, colors and status
   doubleArray3D* work_doubleArray3D = new doubleArray3D(
   	(*local_main_coordinates)[0][0].size(),
   	(*local_main_coordinates)[0][1].size(),
	DIM == 2 ? 1 : (*local_main_coordinates)[0][2].size());
   vector< doubleArray3D >* work_values = new vector< doubleArray3D >
   	(NB_COMPS,*work_doubleArray3D);
   delete work_doubleArray3D;
   VALUES = new vector< vector< doubleArray3D > >(STO_DEPTH,*work_values);
   work_values->clear();
   delete work_values;
   intArray3D* work_intArray3D = new intArray3D(
   	(*local_main_coordinates)[0][0].size(),
   	(*local_main_coordinates)[0][1].size(),
	DIM == 2 ? 1 : (*local_main_coordinates)[0][2].size(), 
	FV_BC_INTERIOR);
   DOFcolors = new vector< intArray3D >( 1, *work_intArray3D );
   work_intArray3D->set( FV_DOF_ONPROC );
   DOFstatus = new vector< intArray3D >( 1, *work_intArray3D );      
   delete work_intArray3D;    
	
   // Set dof colors    
   set_DOF_colors( 0 );

   // Set periodic features
   size_t_vector* periodic_depth = NULL ;
   for (size_t i=0;i<DIM && !periodic_depth;++i) 
     if ( (*PERIODIC)(i) )  
     {
       periodic_depth = new size_t_vector( DIM, 0 );
       FV_SHIFT_TRIPLET mst;
       mst.i = 0 ;
       mst.j = 0 ;
       mst.k = 0 ;
       PERIODIC_SHIFT = new vector< FV_SHIFT_TRIPLET >( 1, mst ) ;
     }
   for (size_t i=0;i<DIM;++i) 
     if ( (*PERIODIC)(i) ) 
     {
       (*periodic_depth)(i) = security_bandwidth ;
       (*PERIODIC_SHIFT)[0].index(i) = (*global_main_coordinates)[0][i].size()
       		- 2 * security_bandwidth ;
     }     
   set_PERIODIC_default( 0, periodic_depth );
   if ( periodic_depth ) delete periodic_depth ;
   
   // Set dof status
   set_DOF_status(0);   
   
   // Compute cell size
   local_cell_size = new vector < vector< doubleVector > >;
   local_cell_size->reserve( 1 );
   local_cell_size->push_back( vwork_double ) ;
   vector< doubleVector > const* primary_coordinates = PRIMARY_GRID->
   	get_local_main_coordinates() ; 
   size_t shift = 0, start_index = 0, end_index = 0 ;        
   for (size_t i=0;i<DIM;++i)
   {     
     (*local_cell_size)[0][i].re_initialize( (*local_dof_number)[0](i), 0. );
     if ( (*local_min_index_in_global)[0](i) == 0 
     	&& !(*PERIODIC)(i) )
     {
       shift = 1 ;
       start_index = 1 ;       
       (*local_cell_size)[0][i](0) = (*local_main_coordinates)[0][i](1)
       	- (*local_main_coordinates)[0][i](0) ;
     }
     else
     { 
       shift = 0 ;
       start_index = 0 ;
     }
     
     if ( (*local_max_index_in_global)[0](i) == (*global_max_index)[0](i) 
     	&& !(*PERIODIC)(i) )
     {
       end_index = (*local_dof_number)[0](i)-2 ;
       (*local_cell_size)[0][i](end_index+1) = 
       	(*local_main_coordinates)[0][i](end_index+1)
       	- (*local_main_coordinates)[0][i](end_index) ;     
     }
     else end_index = (*local_dof_number)[0](i)-1 ;       
     
     for (size_t j=start_index;j<=end_index;++j)
       (*local_cell_size)[0][i](j) = 
       		(*primary_coordinates)[i](j+1-shift)
     		- (*primary_coordinates)[i](j-shift) ;
   }               
}




//----------------------------------------------------------------------
FV_DiscreteField_Centered:: ~FV_DiscreteField_Centered( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: ~FV_DiscreteField_Centered") ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: print( std::ostream& os, size_t indent_width,
      	bool b_values ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: print" ) ;

   std::string space( indent_width, ' ' ) ;   
   
   os << space << "FV_DiscreteField_Centered n°" << ID << endl;
   os << space << "   Name = " << FNAME << endl;   
   os << space << "   Discretization type = " << FDISCRETIZATION << endl;
   os << space << "   Number of components = " << NB_COMPS << endl;
   os << space << "   Level of storage = " << STO_DEPTH << endl;
   os << space << "Field mesh (on all processors)" << endl;
   os << space << "-------------------------------" << endl;
   for (size_t i=0;i<DIM;++i)
   {
     os << space << "Global number of points in direction " << 
     	(*FV_Mesh::directionName)(i) <<
     	" = " << (*global_main_coordinates)[0][i].size() << endl;
     os << space << "Min index = " << 0 << "  Max index = " 
     	<< (*global_max_index)[0](i) << endl;
     os << space << "Coordinate values = ";
     for (size_t j=0;j<(*global_main_coordinates)[0][i].size();++j)
       os << space << " " << (*global_main_coordinates)[0][i](j);
     os << endl;
   }
   os << space << "Local field on processor " << 
   	MAC_Exec::communicator()->rank() << endl;
   os << space << "---------------------------------" << endl;
   os << space << "Mesh" << endl;
   os << space << "----" << endl;   
   for (size_t i=0;i<DIM;++i)
   {
     os << space << "Local number of points in direction " << 
     	(*FV_Mesh::directionName)(i) <<
     	" = " << (*local_main_coordinates)[0][i].size() << endl;
     os << space << "Min index = " << (*local_min_index_in_global)[0](i) 
     	<< "  Max index = " <<
     	(*local_max_index_in_global)[0](i) << endl;	
     os << space << "Coordinate values = ";
     for (size_t j=0;j<(*local_main_coordinates)[0][i].size();++j)
       os << space << " " << (*local_main_coordinates)[0][i](j);
     os << endl;
     os << space << "Status = ";
     for (size_t j=0;j<(*on_current_processor)[0][i].size();++j)
       os << space << " " << (*on_current_processor)[0][i](j);
     os << endl;     
     os << space << "Coordinates in halo zone = ";
     for (size_t j=0;j<(*on_current_processor)[0][i].size();++j)
       if ( (*on_current_processor)[0][i](j) == 2 
       	|| (*on_current_processor)[0][i](j) == 4 ) 
         os << space << " " << (*local_main_coordinates)[0][i](j);	 
     os << endl;
     os << space << "Coordinates in buffer zone = ";
     for (size_t j=0;j<(*on_current_processor)[0][i].size();++j)
       if ( (*on_current_processor)[0][i](j) == 1 
       	|| (*on_current_processor)[0][i](j) == 3 
       	|| (*on_current_processor)[0][i](j) == 5 ) 
         os << space << " " << (*local_main_coordinates)[0][i](j);	 
     os << endl;
     if ( PERIODIC_SHIFT )
       os << "Periodic shift = " << (*PERIODIC_SHIFT)[0].get_index(i) << endl; 
     os << space << "Cell size = ";
     for (size_t j=0;j<(*local_cell_size)[0][i].size();++j)
       os << space << " " << (*local_cell_size)[0][i](j);
     os << endl;         	     	 
   } 
   
   print_BC_FieldValues( os, indent_width, b_values ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: build_BCs( MAC_ModuleExplorer const* exp, 
      	FV_DomainBuilder const* DB ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: build_BCs" ) ;

   // Add triplets to BCs   
   size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[0](2);
   for (size_t i=0;i<(*local_dof_number)[0](0);++i)
     for (size_t j=0;j<(*local_dof_number)[0](1);++j) 
       for (size_t k=0;k<kmax;++k) 
	 if ( (*DOFcolors)[0](i,j,k) > 1 )
	   (*V_BCS)[(*DOFcolors)[0](i,j,k)]->add_MacTriplet( NB_COMPS, i, 
	   	j, k );
   
   // Read BCs
   read_BCs( exp, DB ) ;    
   
   // Set free_DOF_on_boundary update features 
   set_freeDOFonBC_update_features() ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Centered:: get_DOF_coordinate( size_t i, size_t component, 
      	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: get_DOF_coordinate" ) ;
   MAC_CHECK_PRE( i < (*local_dof_number)[0](direction) );
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( direction < DIM ) ;   
   
   return( (*local_main_coordinates)[0][direction](i) );
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Centered:: get_DOF_coordinate_Assembling( int i, 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: get_DOF_coordinate_Assembling" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( direction < DIM ) ;   

   double result = i < 0 ? 0. : i >= int((*local_dof_number)[0](direction)) ?
   	0. : (*local_main_coordinates)[0][direction](i) ;
   
   return( result );
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Centered:: is_global_triplet( 
	int i, int j, int k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: is_global_triplet" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );  
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   

   bool result = 
   	( i >= 0 && i <= int((*global_max_index)[0](0))
   	&& j >= 0 && j <= int((*global_max_index)[0](1)) );
   if ( DIM == 3 && result )
     result = ( k >= 0 && k <= int((*global_max_index)[0](2)) );
   
   return ( result ) ;  
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Centered:: is_global_triplet_local_DOF( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: is_global_triplet_local_DOF" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );  
   MAC_CHECK_PRE( i <= (*global_max_index)[0](0) );
   MAC_CHECK_PRE( j <= (*global_max_index)[0](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, k <= (*global_max_index)[0](2) ) ) ;  
   
   bool result = 
   	( i >= (*local_min_index_in_global)[0](0) 
	&& i <= (*local_max_index_in_global)[0](0)
   	&& j >= (*local_min_index_in_global)[0](1) 
	&& j <= (*local_max_index_in_global)[0](1) );
   if ( DIM == 3 && result )
     result = ( k >= (*local_min_index_in_global)[0](2) 
     			&& k <= (*local_max_index_in_global)[0](2) );
   
   return ( result ) ;  	
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: set_postprocessing_options( 
		std::string const& a_location,
		std::string const& a_paraview_fname ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_post_processing_options" ) ;
   MAC_CHECK_PRE( a_location == "at_cell_centers"
   		|| a_location == "at_vertices" ) ;   
   
   LOCATION = a_location ;
   PARAVIEW_FNAME = a_paraview_fname ;

   // Set the number of dof written in post-processing output files
   set_nb_dof_post_processing() ;

}   




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: copy_DOFs_value(
      size_t source_level, size_t target_level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: copy_DOFs_value" ) ;
   MAC_CHECK_PRE( source_level < STO_DEPTH ); 
   MAC_CHECK_PRE( target_level < STO_DEPTH );
   
   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[0](2) ;
     for (size_t i=0;i<(*local_dof_number)[0](0);++i)
       for (size_t j=0;j<(*local_dof_number)[0](1);++j)
         for (size_t k=0;k<kmax;++k)
     (*VALUES)[target_level][comp]( i, j, k )
      = (*VALUES)[source_level][comp]( i, j, k ) ;
   }          
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: add_to_DOFs_value( size_t component, size_t level,
          double const& val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: add_to_DOFs_value" ) ;
   MAC_CHECK_PRE( level < STO_DEPTH ); 
   MAC_CHECK_PRE( component < NB_COMPS );
   
   size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[0](2) ;
   for (size_t i=0;i<(*local_dof_number)[0](0);++i)
     for (size_t j=0;j<(*local_dof_number)[0](1);++j)
       for (size_t k=0;k<kmax;++k)
         if ( (*UNK_LOCAL_NUMBERING)[component](i,j,k) != -1 )
     (*VALUES)[level][component]( i, j, k )
      += val ;
    
   bool set_bc_values = SET_BC_VALUES_ALLOWED ;
   SET_BC_VALUES_ALLOWED = true ;        
   for ( list< FV_BoundaryCondition* >::iterator il=SET_OF_BCS.begin();
    il!=SET_OF_BCS.end();il++)
     (*il)->set_free_DOF_values( this, component, level ) ; 
   SET_BC_VALUES_ALLOWED =  set_bc_values ;          
}   




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: write_field( MAC_Module* point_data,
                               MAC_Module* cell_data ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: write_field" ) ;

   size_t iter = 0;
   bool scalar = ( NB_COMPS == 1 ) ;

   MAC_Module* target = cell_data ;

   // Surement mieux à faire, mais je ne sais pas faire simple avec un pointeur
   // pour éviter la recopie (GVi)
   doubleArray2D X( NB_COMPS, NB_DOF_POSTPROCESSING_PER_COMP );
   size_t nelem0, nelem1, nelem2 ;
   vector< boolVector > const* primary_mesh_handles_node = 
   	PRIMARY_GRID->get_nodes_owner();

   // In case of periodic flow, periodic part of the pressure is added
   // at the post processing stage, the reconstruction implies that pressure is
   // 0 at the upper face of the domain in the periodic direction
   size_t dir = 0 ;
   doubleVector dpdl( DIM, 0. ), cmax( DIM, 0. ) ;
   bool b_periodic_pressure = false ;
   if ( FNAME == "pressure" && PRIMARY_GRID->is_periodic_flow() 
   	&& NB_COMPS == 1 )
   {
     dir = PRIMARY_GRID->get_periodic_flow_direction() ;
     dpdl(dir) = PRIMARY_GRID->get_periodic_pressure_drop() /
   	( PRIMARY_GRID->get_main_domain_max_coordinate( dir )
	- PRIMARY_GRID->get_main_domain_min_coordinate( dir ) ) ;
     cmax(dir) = PRIMARY_GRID->get_main_domain_max_coordinate( dir ) ;
     b_periodic_pressure = true ;
   }
	
   if ( LOCATION == "at_cell_centers" )
   {
     nelem0 = (*DOFstatus)[0].index_bound(0) ;
     nelem1 = (*DOFstatus)[0].index_bound(1) ;
     nelem2 = (*DOFstatus)[0].index_bound(2) ;

     for (size_t comp=0;comp<NB_COMPS;++comp)
     {
       iter = 0 ;
       for (size_t k=0;k<nelem2;++k)
         for (size_t j=0;j<nelem1;++j) 
           for (size_t i=0;i<nelem0;++i)
             if ( (*DOFcolors)[0](i,j,k) < 2
	   	&& (*DOFstatus)[0](i,j,k) != FV_DOF_HALOZONE
	   	&& (*DOFstatus)[0](i,j,k) != FV_DOF_PERIODIC_HALOZONE )
	     {
	       X( comp, iter ) = (*VALUES)[0][comp](i,j,k);
	       iter++;
	     }
     }
     
     if ( b_periodic_pressure )
     {
       if ( DIM == 2 )
       {
         iter = 0 ;
         for (size_t j=0;j<nelem1;++j) 
           for (size_t i=0;i<nelem0;++i)
             if ( (*DOFcolors)[0](i,j,0) < 2
	   	&& (*DOFstatus)[0](i,j,0) != FV_DOF_HALOZONE
	   	&& (*DOFstatus)[0](i,j,0) != FV_DOF_PERIODIC_HALOZONE )
	       {
	         X( 0, iter ) += 
	dpdl(0) * ( (*local_main_coordinates)[0][0](i) - cmax(0) )
	+ dpdl(1) * ( (*local_main_coordinates)[0][1](j) - cmax(1) ) ;
	         iter++;
	       }
       }
       else
       {
         iter = 0 ;
         for (size_t k=0;k<nelem2;++k)
           for (size_t j=0;j<nelem1;++j) 
             for (size_t i=0;i<nelem0;++i)
               if ( (*DOFcolors)[0](i,j,k) < 2
	   	&& (*DOFstatus)[0](i,j,k) != FV_DOF_HALOZONE
	   	&& (*DOFstatus)[0](i,j,k) != FV_DOF_PERIODIC_HALOZONE )
	       {
	         X( 0, iter ) += 
	dpdl(0) * ( (*local_main_coordinates)[0][0](i) - cmax(0) )
	+ dpdl(1) * ( (*local_main_coordinates)[0][1](j) - cmax(1) )
	+ dpdl(2) * ( (*local_main_coordinates)[0][2](k) - cmax(2) ) ;
	         iter++;
	       }       
       }       
     }     
   }
   else if ( LOCATION == "at_vertices" ) 
   {
     target = point_data ;

     nelem0 = PRIMARY_GRID->get_local_nb_points(0) ;
     nelem1 = PRIMARY_GRID->get_local_nb_points(1) ;
     nelem2 =  DIM == 2 ?
   		1 : PRIMARY_GRID->get_local_nb_points(2) ;

     size_t shift_i = (*local_min_index_in_global)[0](0) == 0 
     	&& !(*PERIODIC)(0) ? 1 : 0;
     size_t shift_j = (*local_min_index_in_global)[0](1) == 0 
     	&& !(*PERIODIC)(1) ? 1 : 0;
     size_t shift_k = DIM == 2 ?
   		0 : (*local_min_index_in_global)[0](2) == 0 
     	&& !(*PERIODIC)(2) ? 1 : 0;

     if (DIM == 2)
     {
       for (size_t comp=0;comp<NB_COMPS;++comp)
       {
         iter = 0 ;
         for (size_t k=0;k<nelem2;++k)
           for (size_t j=0;j<nelem1;++j) 
             for (size_t i=0;i<nelem0;++i)
	       if ( (*primary_mesh_handles_node)[0](i)
	       		&& (*primary_mesh_handles_node)[1](j) )
	       {
	         X( comp, iter ) = interpolate_values_at_nodes( 
				i, shift_i, j, shift_j, k, shift_k, comp, 0 );
	         iter++;
	       }
       }
       
       if ( b_periodic_pressure )
       {
         vector< doubleVector > const* primary_mesh_local_coordinates = 
	 	PRIMARY_GRID->get_local_main_coordinates() ;
	 iter = 0 ;
         for (size_t j=0;j<nelem1;++j) 
           for (size_t i=0;i<nelem0;++i)
	     if ( (*primary_mesh_handles_node)[0](i)
		&& (*primary_mesh_handles_node)[1](j) )
	     {
	       X( 0, iter ) += 
	dpdl(0) * ( (*primary_mesh_local_coordinates)[0](i) - cmax(0) )
	+ dpdl(1) * ( (*primary_mesh_local_coordinates)[1](j) - cmax(1) ) ;
	       iter++;
	     }
       }            
     }     
     else
     {
       for (size_t comp=0;comp<NB_COMPS;++comp)
       {
         iter = 0 ;
         for (size_t k=0;k<nelem2;++k)
           for (size_t j=0;j<nelem1;++j) 
             for (size_t i=0;i<nelem0;++i)
	       if ( (*primary_mesh_handles_node)[0](i)
	     		&& (*primary_mesh_handles_node)[1](j)
	     		&& (*primary_mesh_handles_node)[2](k) )
	       {
	         X( comp, iter ) = interpolate_values_at_nodes( 
				i, shift_i, j, shift_j, k, shift_k, comp, 0 );
	         iter++;
	       }
       } 
       
       if ( b_periodic_pressure )
       {
         vector< doubleVector > const* primary_mesh_local_coordinates = 
	 	PRIMARY_GRID->get_local_main_coordinates() ;
         iter = 0 ;
         for (size_t k=0;k<nelem2;++k)
           for (size_t j=0;j<nelem1;++j) 
             for (size_t i=0;i<nelem0;++i)
	       if ( (*primary_mesh_handles_node)[0](i)
	     		&& (*primary_mesh_handles_node)[1](j)
	     		&& (*primary_mesh_handles_node)[2](k) )
	       {
	         X( 0, iter ) += 
	dpdl(0) * ( (*primary_mesh_local_coordinates)[0](i) - cmax(0) )
	+ dpdl(1) * ( (*primary_mesh_local_coordinates)[1](j) - cmax(1) )
	+ dpdl(2) * ( (*primary_mesh_local_coordinates)[2](k) - cmax(2) ) ;
	         iter++;
	       }
       }         
     }     
   }
   else
   {
     ostringstream mesg ;
     mesg << "*** Unable to save field of name \""
        << FNAME << "\" at location \"" << LOCATION << "\"." << endl << endl ;
     mesg << "    allowed locations are:" << endl;
     mesg << "       - \"at_vertices\"" ;
     mesg << "       - \"at_cell_centers\"" ;
     MAC_Error::object()->raise_plain( mesg.str() ) ;
   }


   if ( scalar && !target->has_entry("Scalars") ) 
      target->add_entry( "Scalars",
      				MAC_String::create( target, PARAVIEW_FNAME ) ) ;
   else if ( !scalar && !target->has_entry("Vectors") ) 
      target->add_entry( "Vectors",
      				MAC_String::create( target, PARAVIEW_FNAME ) ) ;
 
   target->add_entry( PARAVIEW_FNAME, MAC_DoubleArray2D::create( target, X ) ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: set_nb_dof_post_processing( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: set_nb_dof_post_processing") ;

   if( LOCATION == "at_cell_centers" )
   {
     NB_DOF_POSTPROCESSING_PER_COMP = 0 ;
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[0](2);
     for (size_t i=0;i<(*local_dof_number)[0](0);++i)
       for (size_t j=0;j<(*local_dof_number)[0](1);++j) 
         for (size_t k=0;k<kmax;++k)
           if ( (*DOFcolors)[0](i,j,k) < 2
	   	&& (*DOFstatus)[0](i,j,k) != FV_DOF_HALOZONE
	   	&& (*DOFstatus)[0](i,j,k) != FV_DOF_PERIODIC_HALOZONE ) 
	     NB_DOF_POSTPROCESSING_PER_COMP++ ;
   }
   else if( LOCATION == "at_vertices" )
   {
     NB_DOF_POSTPROCESSING_PER_COMP = 1 ;
     for (size_t i=0;i<DIM;++i)
       NB_DOF_POSTPROCESSING_PER_COMP
       		*= PRIMARY_GRID->get_local_nb_points_on_current_proc( i );
   }
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Centered:: interpolate_values_at_nodes( 
      		size_t i, size_t shift_i,
		size_t j, size_t shift_j,
		size_t k, size_t shift_k,
		size_t component,
               	size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: interpolate_values_at_nodes") ;
//   MAC_CHECK_PRE( LOCATION == "at_vertices" )
   
   double result = 0.;
   double coeff_bottom_left, coeff_top_left, 
   		coeff_bottom_right, coeff_top_right;
   double coeff_behind_bottom_left, coeff_behind_top_left, 
   		coeff_behind_bottom_right, coeff_behind_top_right,
		coeff_front_bottom_left, coeff_front_top_left, 
   		coeff_front_bottom_right, coeff_front_top_right;
   
   if ( DIM == 2 )
   {
     // bottom_left
     if ( (*DOFcolors)[0](i,j,k) == FV_BC_BOTTOM_LEFT )
       result = (*VALUES)[level][component](i,j,k) ;
     // bottom_right
     else if ( (*DOFcolors)[0](i+shift_i,j,k) == FV_BC_BOTTOM_RIGHT )
       result = (*VALUES)[level][component](i+shift_i,j,k) ;
     // top_left
     else if ( (*DOFcolors)[0](i,j+shift_j,k) == FV_BC_TOP_LEFT )
       result = (*VALUES)[level][component](i,j+shift_j,k) ;
     // top_right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j,k) == FV_BC_TOP_RIGHT )
       result = (*VALUES)[level][component](i+shift_i,j+shift_j,k) ;
     // left
     else if ( (*DOFcolors)[0](i,j+shift_j-1,k) == FV_BC_LEFT )
     {
       coeff_top_left = (*local_cell_size)[0][1](j+shift_j);
       coeff_bottom_left = (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_top_left
       			* (*VALUES)[level][component](i,j+shift_j,k)
       		+ coeff_bottom_left
			* (*VALUES)[level][component](i,j+shift_j-1,k) )
		 / ( coeff_top_left + coeff_bottom_left ) ;
     }
     // right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j-1,k) == FV_BC_RIGHT )
     {
       coeff_top_right = (*local_cell_size)[0][1](j+shift_j);
       coeff_bottom_right = (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_top_right 
			* (*VALUES)[level][component](i+shift_i,j+shift_j,k)
       		+ coeff_bottom_right
			* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k) )
		/ ( coeff_top_right + coeff_bottom_right ) ;
     }
     // bottom
     else if ( (*DOFcolors)[0](i+shift_i-1,j,k) == FV_BC_BOTTOM )
     {
       coeff_bottom_right = (*local_cell_size)[0][0](i+shift_i);
       coeff_bottom_left = (*local_cell_size)[0][0](i+shift_i-1);
       result = ( coeff_bottom_right 
       			* (*VALUES)[level][component](i+shift_i,j,k)
  		+ coeff_bottom_left 
			* (*VALUES)[level][component](i+shift_i-1,j,k) )
		/ ( coeff_bottom_right + coeff_bottom_left ) ;
     }
     // top
     else if ( (*DOFcolors)[0](i+shift_i-1,j+shift_j,k) == FV_BC_TOP )
     {
       coeff_top_right = (*local_cell_size)[0][0](i+shift_i);
       coeff_top_left = (*local_cell_size)[0][0](i+shift_i-1);
       result = ( coeff_top_right 
       			* (*VALUES)[level][component](i+shift_i,j+shift_j,k)
  		+ coeff_top_left 
			* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k) )
		/ ( coeff_top_right + coeff_top_left ) ;
     }
     // interior
     else
     {
       coeff_top_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j);
       coeff_top_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j);
       coeff_bottom_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j-1);
       coeff_bottom_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_top_right
       			* (*VALUES)[level][component](i+shift_i,j+shift_j,k)
		+ coeff_top_left
			* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k)
		+ coeff_bottom_right
			* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k)
		+ coeff_bottom_left
			* (*VALUES)[level][component](i+shift_i-1,j+shift_j-1,k) )
		/ ( coeff_top_right + coeff_top_left
			+ coeff_bottom_right + coeff_bottom_left );
     }
   }
   else
   {
     // behind_bottom_left
     if ( (*DOFcolors)[0](i,j,k) == FV_BC_BEHIND_BOTTOM_LEFT )
       result = (*VALUES)[level][component](i,j,k) ;
     // behind_bottom_right
     else if ( (*DOFcolors)[0](i+shift_i,j,k) == FV_BC_BEHIND_BOTTOM_RIGHT )
       result = (*VALUES)[level][component](i+shift_i,j,k) ;
     // behind_top_left
     else if ( (*DOFcolors)[0](i,j+shift_j,k) == FV_BC_BEHIND_TOP_LEFT )
       result = (*VALUES)[level][component](i,j+shift_j,k) ;
     // behind_top_right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j,k) == 
     	FV_BC_BEHIND_TOP_RIGHT )
       result = (*VALUES)[level][component](i+shift_i,j+shift_j,k) ;
     // front_bottom_left
     else if ( (*DOFcolors)[0](i,j,k+shift_k) == FV_BC_FRONT_BOTTOM_LEFT )
       result = (*VALUES)[level][component](i,j,k+shift_k) ;
     // front_bottom_right
     else if ( (*DOFcolors)[0](i+shift_i,j,k+shift_k) == 
     	FV_BC_FRONT_BOTTOM_RIGHT )
       result = (*VALUES)[level][component](i+shift_i,j,k+shift_k) ;
     // front_top_left
     else if ( (*DOFcolors)[0](i,j+shift_j,k+shift_k) == FV_BC_FRONT_TOP_LEFT )
       result = (*VALUES)[level][component](i,j+shift_j,k+shift_k) ;
     // front_top_right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j,k+shift_k) == 
     	FV_BC_FRONT_TOP_RIGHT )
       result = (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k) ;
     // bottom_left
     else if ( (*DOFcolors)[0](i,j,k+shift_k-1) == FV_BC_BOTTOM_LEFT )
     {
       coeff_front_bottom_left = (*local_cell_size)[0][2](k+shift_k);
       coeff_behind_bottom_left = (*local_cell_size)[0][2](k+shift_k-1);
       result = ( coeff_front_bottom_left
       			*(*VALUES)[level][component](i,j,k+shift_k)
        	+ coeff_behind_bottom_left 
			*(*VALUES)[level][component](i,j,k+shift_k-1) )
		/ ( coeff_front_bottom_left + coeff_behind_bottom_left ) ;
     }
     // bottom_right
     else if ( (*DOFcolors)[0](i+shift_i,j,k+shift_k-1) == FV_BC_BOTTOM_RIGHT )
     {
       coeff_front_bottom_right = (*local_cell_size)[0][2](k+shift_k);
       coeff_behind_bottom_right = (*local_cell_size)[0][2](k+shift_k-1);
       result = ( coeff_front_bottom_right
       			* (*VALUES)[level][component](i+shift_i,j,k+shift_k)
       		+ coeff_behind_bottom_right
			* (*VALUES)[level][component](i+shift_i,j,k+shift_k-1) )
		/ ( coeff_front_bottom_right + coeff_behind_bottom_right ) ;
     }
     // top_left
     else if ( (*DOFcolors)[0](i,j+shift_j,k+shift_k-1) == FV_BC_TOP_LEFT )
     {
       coeff_front_top_left = (*local_cell_size)[0][2](k+shift_k);
       coeff_behind_top_left = (*local_cell_size)[0][2](k+shift_k-1);
       result = ( coeff_front_top_left
       			* (*VALUES)[level][component](i,j+shift_j,k+shift_k)
       		+ coeff_behind_top_left
			* (*VALUES)[level][component](i,j+shift_j,k+shift_k-1) )
		/ ( coeff_front_top_left + coeff_behind_top_left ) ;
     }
     // top_right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j,k+shift_k-1) == 
     	FV_BC_TOP_RIGHT )
     {
       coeff_front_top_right = (*local_cell_size)[0][2](k+shift_k);
       coeff_behind_top_right = (*local_cell_size)[0][2](k+shift_k-1);
       result = ( coeff_front_top_right
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k)
       		+ coeff_behind_top_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k-1) )
		/ ( coeff_front_top_right + coeff_behind_top_right ) ;
     }
     // behind_left
     else if ( (*DOFcolors)[0](i,j+shift_j-1,k) == FV_BC_BEHIND_LEFT )
     {
       coeff_behind_top_left = (*local_cell_size)[0][1](j+shift_j);
       coeff_behind_bottom_left = (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_behind_top_left
       			* (*VALUES)[level][component](i,j+shift_j,k)
       		+ coeff_behind_bottom_left
			* (*VALUES)[level][component](i,j+shift_j-1,k) )
		/ ( coeff_behind_top_left + coeff_behind_bottom_left ) ;
     }
     // behind_right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j-1,k) == FV_BC_BEHIND_RIGHT )
     {
       coeff_behind_top_right = (*local_cell_size)[0][1](j+shift_j);
       coeff_behind_bottom_right = (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_behind_top_right
       			* (*VALUES)[level][component](i+shift_i,j+shift_j,k)
       		+ coeff_behind_bottom_right
			* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k) )
		/ ( coeff_behind_top_right + coeff_behind_bottom_right ) ;
     }
     // front_left
     else if ( (*DOFcolors)[0](i,j+shift_j-1,k+shift_k) == FV_BC_FRONT_LEFT )
     {
       coeff_front_top_left = (*local_cell_size)[0][1](j+shift_j);
       coeff_front_bottom_left = (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_front_top_left
       			* (*VALUES)[level][component](i,j+shift_j,k+shift_k)
       		+ coeff_front_bottom_left
			* (*VALUES)[level][component](i,j+shift_j-1,k+shift_k) )
		/ ( coeff_front_top_left + coeff_front_bottom_left ) ;
     }
     // front_right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j-1,k+shift_k) == 
     	FV_BC_FRONT_RIGHT )
     {
       coeff_front_top_right = (*local_cell_size)[0][1](j+shift_j);
       coeff_front_bottom_right = (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_front_top_right
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k)
       		+ coeff_front_bottom_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k+shift_k) )
		/ ( coeff_front_top_right + coeff_front_bottom_right ) ;
     }
     // behind_bottom
     else if ( (*DOFcolors)[0](i+shift_i-1,j,k) == FV_BC_BEHIND_BOTTOM )
     {
       coeff_behind_bottom_right = (*local_cell_size)[0][0](i+shift_i);
       coeff_behind_bottom_left = (*local_cell_size)[0][0](i+shift_i-1);
       result = ( coeff_behind_bottom_right
       			* (*VALUES)[level][component](i+shift_i,j,k)
       		+ coeff_behind_bottom_left
			* (*VALUES)[level][component](i+shift_i-1,j,k) )
		/ ( coeff_behind_bottom_right + coeff_behind_bottom_left ) ;
     }
     // behind_top
     else if ( (*DOFcolors)[0](i+shift_i-1,j+shift_j,k) == FV_BC_BEHIND_TOP )
     {
       coeff_behind_top_right = (*local_cell_size)[0][0](i+shift_i);
       coeff_behind_top_left = (*local_cell_size)[0][0](i+shift_i-1);
       result = ( coeff_behind_top_right
       			* (*VALUES)[level][component](i+shift_i,j+shift_j,k)
       		+ coeff_behind_top_left
			* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k) )
		/ ( coeff_behind_top_right + coeff_behind_top_left ) ;
     }
     // front_bottom
     else if ( (*DOFcolors)[0](i+shift_i-1,j,k+shift_k) == 
     	FV_BC_BEHIND_BOTTOM )
     {
       coeff_front_bottom_right = (*local_cell_size)[0][0](i+shift_i);
       coeff_front_bottom_left = (*local_cell_size)[0][0](i+shift_i-1);
       result = ( coeff_front_bottom_right
       			* (*VALUES)[level][component](i+shift_i,j,k+shift_k)
       		+ coeff_front_bottom_left
			* (*VALUES)[level][component](i+shift_i-1,j,k+shift_k) )
		/ ( coeff_front_bottom_right + coeff_front_bottom_left ) ;
     }
     // front_top
     else if ( (*DOFcolors)[0](i+shift_i-1,j+shift_j,k+shift_k) == 
     	FV_BC_FRONT_TOP )
     {
       coeff_front_top_right = (*local_cell_size)[0][0](i+shift_i);
       coeff_front_top_left = (*local_cell_size)[0][0](i+shift_i-1);
       result = ( coeff_front_top_right
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k)
       		+ coeff_front_top_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k+shift_k) )
		/ ( coeff_front_top_right + coeff_front_top_left ) ;
     }
     // left
     else if ( (*DOFcolors)[0](i,j+shift_j-1,k+shift_k-1) == FV_BC_LEFT )
     {
       coeff_behind_bottom_left = (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_behind_top_left = (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_front_bottom_left = (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_front_top_left = (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k);
       result = ( coeff_behind_bottom_left
       		* (*VALUES)[level][component](i,j+shift_j,k+shift_k)
 		+ coeff_behind_top_left
		* (*VALUES)[level][component](i,j+shift_j-1,k+shift_k)
 		+ coeff_front_bottom_left
		* (*VALUES)[level][component](i,j+shift_j,k+shift_k-1)
 		+ coeff_front_top_left
		* (*VALUES)[level][component](i,j+shift_j-1,k+shift_k-1) )  
		/ ( coeff_behind_bottom_left + coeff_behind_top_left
		+ coeff_front_bottom_left + coeff_front_top_left ) ;
     }
     // right
     else if ( (*DOFcolors)[0](i+shift_i,j+shift_j-1,k+shift_k-1) == 
     	FV_BC_RIGHT )
     {
       coeff_behind_bottom_right = (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_behind_top_right = (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_front_bottom_right = (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_front_top_right = (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k);
       result = ( coeff_behind_bottom_right
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k)
 		+ coeff_behind_top_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k+shift_k)
 		+ coeff_front_bottom_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k-1)
 		+ coeff_front_top_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k+shift_k-1) )
		/ ( coeff_behind_bottom_right + coeff_behind_top_right
		+ coeff_front_bottom_right + coeff_front_top_right ) ;
     }
     // bottom
     else if ( (*DOFcolors)[0](i+shift_i-1,j,k+shift_k-1) == FV_BC_BOTTOM )
     {
       coeff_behind_bottom_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_behind_bottom_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_front_bottom_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_front_bottom_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][2](k+shift_k);
       result = ( coeff_behind_bottom_left
       		* (*VALUES)[level][component](i+shift_i,j,k+shift_k)
 		+ coeff_behind_bottom_right
		* (*VALUES)[level][component](i+shift_i-1,j,k+shift_k)
 		+ coeff_front_bottom_left
		* (*VALUES)[level][component](i+shift_i,j,k+shift_k-1)
 		+ coeff_front_bottom_right
		* (*VALUES)[level][component](i+shift_i-1,j,k+shift_k-1) )
		/ ( coeff_behind_bottom_left + coeff_behind_bottom_right
		+ coeff_front_bottom_left + coeff_front_bottom_right ) ;
     }
     // top
     else if ( (*DOFcolors)[0](i+shift_i-1,j+shift_j,k+shift_k-1) == 
     	FV_BC_TOP )
     {
       coeff_behind_top_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_behind_top_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_front_top_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_front_top_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][2](k+shift_k);
       result = ( coeff_behind_top_left
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k)
 		+ coeff_behind_top_right
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k+shift_k)
 		+ coeff_front_top_left
		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k-1)
 		+ coeff_front_top_right
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k+shift_k-1) )
		/ ( coeff_behind_top_left + coeff_behind_top_right
		+ coeff_front_top_left + coeff_front_top_right ) ;
     }
     // behind
     else if ( (*DOFcolors)[0](i+shift_i-1,j+shift_j-1,k) == FV_BC_BEHIND )
     {
       coeff_behind_top_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j);
       coeff_behind_top_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j);
       coeff_behind_bottom_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j-1);
       coeff_behind_bottom_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_behind_top_right
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k)
 		+ coeff_behind_top_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k)
 		+ coeff_behind_bottom_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k)
 		+ coeff_behind_bottom_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j-1,k) )
		/ ( coeff_behind_top_right + coeff_behind_top_left
		+ coeff_behind_bottom_right + coeff_behind_bottom_left ) ;
     }
     // front
     else if ( (*DOFcolors)[0](i+shift_i-1,j+shift_j-1,k+shift_k) == 
     	FV_BC_FRONT )
     {
       coeff_front_top_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j);
       coeff_front_top_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j);
       coeff_front_bottom_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j-1);
       coeff_front_bottom_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j-1);
       result = ( coeff_front_top_right
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k)
 		+ coeff_front_top_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k+shift_k)
 		+ coeff_front_bottom_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k+shift_k)
 		+ coeff_front_bottom_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j-1,k+shift_k) )
		/ ( coeff_front_top_right + coeff_front_top_left
		+ coeff_front_bottom_right + coeff_front_bottom_left ) ;
     }
     // interior
     else
     {
       coeff_front_top_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_front_bottom_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_behind_top_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_behind_bottom_right = (*local_cell_size)[0][0](i+shift_i)
       			* (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_front_top_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_front_bottom_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k);
       coeff_behind_top_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       coeff_behind_bottom_left = (*local_cell_size)[0][0](i+shift_i-1)
       			* (*local_cell_size)[0][1](j+shift_j-1)
       			* (*local_cell_size)[0][2](k+shift_k-1);
       result = ( coeff_front_top_right
       		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k)
 		+ coeff_front_bottom_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k+shift_k)
 		+ coeff_behind_top_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j,k+shift_k-1)
 		+ coeff_behind_bottom_right
		* (*VALUES)[level][component](i+shift_i,j+shift_j-1,k+shift_k-1)
		+ coeff_front_top_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k+shift_k)
 		+ coeff_front_bottom_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j-1,k+shift_k)
 		+ coeff_behind_top_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j,k+shift_k-1)
 		+ coeff_behind_bottom_left
		* (*VALUES)[level][component](i+shift_i-1,j+shift_j-1,k+shift_k-1) )
		/ ( coeff_front_top_right + coeff_front_bottom_right
		+ coeff_behind_top_right + coeff_behind_bottom_right
		+ coeff_front_top_left + coeff_front_bottom_left
		+ coeff_behind_top_left + coeff_behind_bottom_left ) ;
     }
   }
   
   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Centered:: get_cell_size( size_t i, size_t component, 
      	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: get_cell_size") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS );       
   MAC_CHECK_PRE( check_invariant_cell_features( i, direction ) ) ;   

   double result = (*local_cell_size)[0][direction](i) ;
	         
   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Centered:: get_cell_measure( size_t i, size_t j, size_t k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: get_cell_measure") ;   
   MAC_CHECK_PRE( component < NB_COMPS );    
   MAC_CHECK_PRE( check_invariant_cell_features( i, 0 ) ) ;
   MAC_CHECK_PRE( check_invariant_cell_features( j, 1 ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, check_invariant_cell_features( k, 2 ) ) ) ;
   
   double result = (*local_cell_size)[0][0](i) * (*local_cell_size)[0][1](j) ;
   if ( DIM == 3 ) 
     result *= (*local_cell_size)[0][2](k) ;
	
   return ( result ) ;
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Centered:: DOF_offset( int &i, int &j, int &k,
        size_t_vector center, size_t_vector stencil,
	vector<double> &offset, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: DOF_offset" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 

   bool outdomain_or_BC = false ;
   
   if ( i + int((*local_min_index_in_global)[0](0))< 0 )
   {
     i += - 2 * i ;
     offset[0] = double(  2 * i  - 1)
     	 * (*local_cell_size)[0][0](i);
     outdomain_or_BC = true;
   }
   else if ( i+ int((*local_min_index_in_global)[0](0)) == 0 )
   {
    i = center(0) - stencil(0) - 1;
    i += - 2 * i;
    offset[0] = double(  2 * i - 1)
     	* (*local_cell_size)[0][0](i);
    outdomain_or_BC = true;
   } 
   else if ( i + int((*local_min_index_in_global)[0](0)) ==
   	int((*global_max_index)[0](0)) ) 
   {
    i = center(0) + stencil(0) + 1;
    i += - 2 * ( i + int((*local_min_index_in_global)[0](0))
     - int((*global_max_index)[0](0)) );
    offset[0] = double( 1 + 2 * ( i + int((*local_min_index_in_global)[0](0))
       - int((*global_max_index)[0](0)) ) ) * (*local_cell_size)[0][0](i); 
    outdomain_or_BC = true;
   }	
   else if ( i + int((*local_min_index_in_global)[0](0)) >
   	int((*global_max_index)[0](0)) )
   {
     i += - 2 * ( i + int((*local_min_index_in_global)[0](0))
      - int((*global_max_index)[0](0)) );
     offset[0] = double( 1 + 2 * ( i + int((*local_min_index_in_global)[0](0))
     	- int((*global_max_index)[0](0)) ) ) * (*local_cell_size)[0][0](i); 
     outdomain_or_BC = true;
   }
       
   if ( j + int((*local_min_index_in_global)[0](1)) < 0 )
   {
     j += - 2 * j;
     offset[1] = double( 2 * j - 1 )
     	 * (*local_cell_size)[0][1](j);
     outdomain_or_BC = true;
   }
   else if ( j + int((*local_min_index_in_global)[0](1)) == 0 )
   {
     j = center(1) - stencil(1) - 1;
     j += - 2 * j;
     offset[1] = double( 2 * j - 1 )
     	 * (*local_cell_size)[0][1](j);
     outdomain_or_BC = true;
   } 
   else if ( j + int((*local_min_index_in_global)[0](1)) ==
   	int((*global_max_index)[0](1)) ) 
   {
     j = center(1) + stencil(1) + 1;
     j += - 2 * ( j + int((*local_min_index_in_global)[0](1))
      - int((*global_max_index)[0](1)) );
     offset[1] = double(1 + 2 * ( j + int((*local_min_index_in_global)[0](1))
     	- int((*global_max_index)[0](1)) ) ) * (*local_cell_size)[0][1](j); 
     outdomain_or_BC = true;
   }	
   else if ( j + int((*local_min_index_in_global)[0](1)) >
   	int((*global_max_index)[0](1)) )
   {
     j += - 2 * ( j + int((*local_min_index_in_global)[0](1))
      - int((*global_max_index)[0](1)) );
     offset[1] = double(1 + 2 * ( j + int((*local_min_index_in_global)[0](1))
     	- int((*global_max_index)[0](1)) ) ) * (*local_cell_size)[0][1](j); 
     outdomain_or_BC = true;
   } 
   
   if ( DIM == 3 )
   {
    if ( k + int((*local_min_index_in_global)[0](2)) < 0 )
    {
      k += - 2 * k;
      offset[2] = double(2 * k -1)
     	 * (*local_cell_size)[0][2](k);
      outdomain_or_BC = true;
    }
    else if ( k + int((*local_min_index_in_global)[0](2)) == 0 )
    {
     k = center(2) - stencil(2) - 1;
     k += - 2 * k;
      offset[2] = double(2 * k -1)
     	 * (*local_cell_size)[0][2](k);
     outdomain_or_BC = true;
    }
    else if ( k + int((*local_min_index_in_global)[0](2)) == 
   	int((*global_max_index)[0](2)) )
    {
     k = center(2) + stencil(2) + 1;
     k += - 2 * ( k + int((*local_min_index_in_global)[0](2))
       - int((*global_max_index)[0](2)) );
      offset[2] = double( 1 + 2 * ( k + int((*local_min_index_in_global)[0](2))
     	- int((*global_max_index)[0](2)) ) ) * (*local_cell_size)[0][2](k); 
     outdomain_or_BC = true;
    }	
    else if ( k + int((*local_min_index_in_global)[0](2)) >
   	int((*global_max_index)[0](2)) )
    {
      k += - 2 * ( k + int((*local_min_index_in_global)[0](2))
       - int((*global_max_index)[0](2)) );
      offset[2] = double( 1 + 2 * ( k + int((*local_min_index_in_global)[0](2))
     	- int((*global_max_index)[0](2)) ) ) * (*local_cell_size)[0][2](k);  
      outdomain_or_BC = true;
    }
   }
   
   return ( outdomain_or_BC ) ;

}




//----------------------------------------------------------------------
double
FV_DiscreteField_Centered:: get_face_perp_to_direction_measure( 
	size_t i, size_t j, size_t k, size_t component, size_t direction ) 
	const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Centered:: get_face_perp_to_direction_measure") ;

   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( IMPLIES( direction != 0, 
   	check_invariant_cell_features( i, 0 ) ) ) ;
   MAC_CHECK_PRE( IMPLIES( direction != 1, 
   	check_invariant_cell_features( j, 1 ) ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, IMPLIES( direction != 2, 
   	check_invariant_cell_features( k, 2 ) ) ) ) ;
   
   double result = 0. ;
 
   if ( DIM == 2 )
   {
     switch( direction )
     {
       case 0: 
         result = (*local_cell_size)[0][1](j) ;
         break;
       case 1: 
         result = (*local_cell_size)[0][0](i) ;
         break;       
     }
   }
   else
   {
     switch( direction )
     {
       case 0: 
         result = (*local_cell_size)[0][1](j) * (*local_cell_size)[0][2](k) ;
         break;
       case 1: 
         result = (*local_cell_size)[0][0](i) * (*local_cell_size)[0][2](k) ;
         break;       
       case 2: 
         result = (*local_cell_size)[0][0](i) * (*local_cell_size)[0][1](j) ; 
         break;
     }   
   }
                  
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Centered:: get_min_index_unknown_handled_by_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Centered:: get_min_index_unknown_handled_by_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*min_index_unknown_handled_by_proc)[0](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Centered:: get_max_index_unknown_handled_by_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Centered:: get_max_index_unknown_handled_by_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*max_index_unknown_handled_by_proc)[0](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Centered:: get_min_index_unknown_on_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Centered:: get_min_index_unknown_on_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*min_index_unknown_on_proc)[0](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Centered:: get_max_index_unknown_on_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Centered:: get_max_index_unknown_on_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*max_index_unknown_on_proc)[0](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Centered:: get_local_nb_dof( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Centered:: get_local_nb_dof") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*local_dof_number)[0](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Centered:: DOF_has_imposed_Dirichlet_value( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: DOF_has_imposed_Dirichlet_value" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( i <= (*local_dof_number)[0](0) );
   MAC_CHECK_PRE( j <= (*local_dof_number)[0](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, k <= (*local_dof_number)[0](2) ) ) ;      

   bool result = false ;
   if ( (*DOFcolors)[0](i,j,k) > 1 )
     result = get_BC( (*DOFcolors)[0](i,j,k) )->is_dirichlet( component ) ;
   
   return ( result ) ; 
}	




//----------------------------------------------------------------------
bool
FV_DiscreteField_Centered:: DOF_on_BC( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: DOF_on_BC" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( i <= (*local_dof_number)[0](0) );
   MAC_CHECK_PRE( j <= (*local_dof_number)[0](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, k <= (*local_dof_number)[0](2) ) ) ;      

   bool result = (*DOFcolors)[0](i,j,k) > 1 ;
   
   return ( result ) ; 
}	




//----------------------------------------------------------------------
bool
FV_DiscreteField_Centered:: DOF_is_unknown_handled_by_proc( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: DOF_is_unknown_handled_by_proc" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );
   MAC_CHECK_PRE( i <= (*local_dof_number)[0](0) );
   MAC_CHECK_PRE( j <= (*local_dof_number)[0](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, k <= (*local_dof_number)[0](2) ) ) ;

   bool result = (*UNK_LOCAL_NUMBERING)[component]( i, j, k ) != -1 
   	&& (*DOFstatus)[0](i,j,k) != FV_DOF_HALOZONE
	&& (*DOFstatus)[0](i,j,k) != FV_DOF_PERIODIC_HALOZONE ; 
   
   return ( result ) ; 
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET
FV_DiscreteField_Centered:: shift_staggeredToCentered( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: shift_staggeredToCentered" ) ;

   // If (i,j,k) is the location of the centered dof, 
   // (i+result.i,j+result.j,k+result.k) is the location of the staggered
   // dof on the right, top, front of the centered dof
   
   FV_SHIFT_TRIPLET result ;
   result.i = 1 ;
   result.j = 1 ;   
   result.k = 1 ;
   
   for (size_t i=0;i<DIM;++i)
     if ( (*local_min_index_in_global)[0](i) == 0  && !(*PERIODIC)(i) ) 
       result.index( i ) = 0 ;
       
   return ( result ) ;
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Centered:: check_invariant_cell_features( size_t i,
      	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: check_invariant_cell_features") ;
   MAC_ASSERT( IMPLIES( !(*PERIODIC)(direction),
   	IMPLIES( (*local_min_index_in_global)[0](direction) == 0, 
   	i > 0 ) ) ) ;
   MAC_ASSERT( i < (*local_dof_number)[0](direction) ) ; 
   MAC_ASSERT( IMPLIES( !(*PERIODIC)(direction),
   	IMPLIES( (*local_max_index_in_global)[0](direction) == 
		(*global_max_index)[0](direction) , 
   	i < (*local_dof_number)[0](direction) - 1 ) ) ) ;
	
	
   return ( true ) ;
}





//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: compute_normLinf( double time ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: compute_normLinf" ) ;

}




//----------------------------------------------------------------------
doubleVector const*
FV_DiscreteField_Centered:: get_DOF_coordinates_vector( size_t component, 
      	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: get_DOF_coordinates_vector") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS );
   MAC_CHECK_PRE( local_main_coordinates != NULL ) ;    
   
   doubleVector const* result = &((*local_main_coordinates)[0][direction]) ;

   return ( result ) ;      
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: create_transproj_interpolation( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: create_transproj_interpolation" ) ;

   if ( STO_DEPTH < 2 ) FV_DiscreteField_ERROR::n1( FNAME ) ;

   size_t trans_dir = PRIMARY_GRID->get_translation_direction() ;
   double trans_mag = PRIMARY_GRID->get_translation_magnitude() ;
   size_t ncoor = (*local_main_coordinates)[0][trans_dir].size(), imin ;

   size_t_vector work( ncoor, 0 ) ;
   transproj_interpolation = new vector< size_t_vector >( 1, work ) ;
   
   for (size_t i=0;i<ncoor;++i)
     if ( FV_Mesh::between(&(*local_main_coordinates)[0][trans_dir], 
     	(*local_main_coordinates)[0][trans_dir](i) + trans_mag, imin ) )
       (*transproj_interpolation)[0](i) = imin ;
     else (*transproj_interpolation)[0](i) = OUT_OF_DOMAIN_PROJTRANS ;          
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: translation_projection( size_t const& level, 
      	size_t const& temporary_level, bool translate_mesh,
	doubleVector const* outOfDomain_values )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: translation_projection" ) ;

   if ( !transproj_interpolation ) FV_DiscreteField_ERROR::n2( FNAME ) ;

   size_t m = 0, mmin, i0, i1, j0, j1, k0, k1 ;
   size_t trans_dir = PRIMARY_GRID->get_translation_direction() ;
   double trans_mag = PRIMARY_GRID->get_translation_magnitude(),
   	translated_coordinate = 0., c0 = 0., c1 = 0. ;
   doubleVector odv( NB_COMPS, 0. );
   if ( outOfDomain_values )
     for (size_t comp=0;comp<NB_COMPS;++comp) 
       odv( comp ) = (*outOfDomain_values)( comp ) ;
     
   // Interpolate field at level "level" and store at level "temporary_level"
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t nelem0 = (*VALUES)[0][comp].index_bound(0) ;
     size_t nelem1 = (*VALUES)[0][comp].index_bound(1) ;
     size_t nelem2 = (*VALUES)[0][comp].index_bound(2) ; 
     for (size_t i=0;i<nelem0;++i)
       for (size_t j=0;j<nelem1;++j) 
         for (size_t k=0;k<nelem2;++k)
	   if ( !DOF_has_imposed_Dirichlet_value( i, j, k, comp ) )
	   {
             switch( trans_dir )
	     {
	       case 0: m = i ; break ;
	       case 1: m = j ; break ;
	       case 2: m = k ; break ;
	     }
	     	     	     
             mmin = (*transproj_interpolation)[0](m) ;
	     if ( mmin == OUT_OF_DOMAIN_PROJTRANS )
	       (*VALUES)[temporary_level][comp](i,j,k) = odv( comp ) ;
	     else
	     {
	       translated_coordinate = 
	       	(*local_main_coordinates)[0][trans_dir](m) + trans_mag ;
	       
	       c1 = ( translated_coordinate - 
	       	(*local_main_coordinates)[0][trans_dir](mmin) ) /
		( (*local_main_coordinates)[0][trans_dir](mmin+1)
		- (*local_main_coordinates)[0][trans_dir](mmin) ) ;
	       c0 = 1. - c1 ;
	       
	       i0 = i1 = i ;
	       j0 = j1 = j ;
	       k0 = k1 = k ;
	       switch( trans_dir )
	       {
	         case 0: i0 = mmin ; i1 = mmin+1 ; break ;
	         case 1: j0 = mmin ; j1 = mmin+1 ; break ;
	         case 2: k0 = mmin ; k1 = mmin+1 ; break ;
	       }
	       
	       (*VALUES)[temporary_level][comp](i,j,k) = 
	       		c0 * (*VALUES)[level][comp](i0,j0,k0)
			+ c1 * (*VALUES)[level][comp](i1,j1,k1) ;
	     }
	   }
   }

   // Copy from level "temporary_level" to level "level"
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t nelem0 = (*VALUES)[0][comp].index_bound(0) ;
     size_t nelem1 = (*VALUES)[0][comp].index_bound(1) ;
     size_t nelem2 = (*VALUES)[0][comp].index_bound(2) ; 
     for (size_t i=0;i<nelem0;++i)
       for (size_t j=0;j<nelem1;++j) 
         for (size_t k=0;k<nelem2;++k)
	   (*VALUES)[level][comp](i,j,k) = 
	   	(*VALUES)[temporary_level][comp](i,j,k) ;
   }  
   
   // Translate field mesh
   if ( translate_mesh )
   {
     size_t ncoor = (*global_main_coordinates)[0][trans_dir].size() ;
     for (size_t i=0;i<ncoor;++i)
       (*global_main_coordinates)[0][trans_dir](i) += trans_mag ;
     ncoor = (*local_main_coordinates)[0][trans_dir].size() ;
     for (size_t i=0;i<ncoor;++i)
       (*local_main_coordinates)[0][trans_dir](i) += trans_mag ;
   }  
       
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: restore_translated_field_mesh( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: restore_translated_field_mesh" ) ;

   size_t trans_dir = PRIMARY_GRID->get_translation_direction() ;
   double trans_dist = PRIMARY_GRID->get_translation_distance() ;

   if ( fabs( trans_dist ) > 1e-12 )
   {      
     if ( MAC_Exec::communicator()->rank() == 0 ) 
        FV::out() << "      Restore translated mesh of field " << FNAME
		<< endl ;
     size_t ncoor = (*global_main_coordinates)[0][trans_dir].size() ;
     for (size_t i=0;i<ncoor;++i)
       (*global_main_coordinates)[0][trans_dir](i) += trans_dist ;
     ncoor = (*local_main_coordinates)[0][trans_dir].size() ;
     for (size_t i=0;i<ncoor;++i)
       (*local_main_coordinates)[0][trans_dir](i) += trans_dist ; 
   }    
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: translate_field_mesh( const size_t& trans_dir, 
      	const double &trans_dist )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: translate_field_mesh" ) ;

   size_t ncoor = (*global_main_coordinates)[0][trans_dir].size() ;
   for (size_t i=0;i<ncoor;++i)
     (*global_main_coordinates)[0][trans_dir](i) += trans_dist ;
   ncoor = (*local_main_coordinates)[0][trans_dir].size() ;
   for (size_t i=0;i<ncoor;++i)
     (*local_main_coordinates)[0][trans_dir](i) += trans_dist ; 
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Centered:: compute_boundary_cell_centered_DOF_integral( 
	size_t component, size_t level, std::string const& boundary_name ) 
	const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "FV_DiscreteField_Centered:: compute_boundary_cell_centered_DOF_integral" );
   
   size_t bc_color = FV_DomainBuilder::get_color_number( boundary_name ) ;   
   MAC_CHECK( FV_DomainBuilder::is_main_color( bc_color ) );   
   
   double result = (*V_BCS)[bc_color]
   	->compute_boundary_cell_centered_DOF_integral( this, component, 
	level ) ;
            
   return ( result ) ;
}
