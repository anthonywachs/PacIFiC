#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <FV_BoundaryCondition.hh>
#include <FV_DomainBuilder.hh>
#include <FV_TimeIterator.hh>
#include <FV.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Data.hh>
#include <MAC_DoubleArray3D.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Int.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_String.hh>
#include <MAC_Vector.hh>
#include <MAC_Root.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Variable.hh>
#include <stringVector.hh>
#include <LA_SeqVector.hh>
#include <LA_Vector.hh>
#include <LA_Matrix.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using std::cout ; 
using std::endl ;
using std::string ; 
using std::ostringstream ;


size_t FV_DiscreteField:: NB_INSTANCES = 0 ;


//----------------------------------------------------------------------
FV_DiscreteField*
FV_DiscreteField:: create( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type,
	size_t a_nb_components,
	size_t a_depth )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: create" ) ;
   MAC_CHECK_PRE( a_primary_mesh != 0 ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( a_type == "centered" || a_type == "staggered"
   	|| a_type == "vertex" || a_type == "tensor" || a_type == "vorticity" ) ;   
   MAC_CHECK_PRE( a_nb_components>0 ) ;

   FV_DiscreteField const* proto =
      static_cast<FV_DiscreteField const*>( 
      	plugins_map()->item( a_type ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;

   FV_DiscreteField* result = proto->create_replica( a_owner,
	a_primary_mesh, a_name, a_type, a_nb_components, a_depth ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   MAC_CHECK_POST( result->storage_depth() == a_depth ) ;
   MAC_CHECK_POST( result->nb_components() == a_nb_components ) ;
   MAC_CHECK_POST( result->primary_grid() == a_primary_mesh ) ;
   MAC_CHECK_POST( result->discretization_type() == a_type ) ;      
   return( result ) ;
}




//----------------------------------------------------------------------
FV_DiscreteField:: FV_DiscreteField( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type,
	size_t a_nb_components,
	size_t a_depth )	
//----------------------------------------------------------------------
    : MAC_Object( a_owner ) 
   , PRIMARY_GRID( a_primary_mesh )
   , FNAME( a_name )
   , FDISCRETIZATION( a_type )
   , ID( NB_INSTANCES++ )
   , DIM( a_primary_mesh->nb_space_dimensions() )
   , NB_COMPS( a_nb_components )
   , NB_CELLS( 0 )   
   , STO_DEPTH( a_depth ) 
   , VALUES( 0 ) 
   , V_BCS( 0 )
   , UNK_LOCAL_NUMBERING( 0 )
   , UNK_GLOBAL_NUMBERING( 0 )
   , NB_LOCAL_UNKNOWNS( 0 )
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC( 0 )
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE( 0 )
   , NB_LOCAL_DOF( 0 )
   , NB_GLOBAL_UNKNOWNS( 0 )
   , INITIALIZER( 0 )
   , CTX( 0 )
   , COORDS( 0 )
   , global_max_index( 0 )
   , local_max_index_in_global( 0 )    
   , local_min_index_in_global( 0 )
   , local_dof_number( 0 )
   , max_index_unknown_handled_by_proc( 0 )
   , min_index_unknown_handled_by_proc( 0 )
   , max_index_unknown_on_proc( 0 )
   , min_index_unknown_on_proc( 0 )   
   , global_main_coordinates( 0 )
   , local_main_coordinates( 0 )
   , local_cell_size( 0 )   
   , DOFcolors( 0 )
   , DOFstatus( 0 )   
   , PERIODIC( 0 )
   , PERIODIC_SHIFT( 0 )
   , PARAVIEW_FNAME( "" )
   , LOCATION( "" )
   , NB_DOF_POSTPROCESSING_PER_COMP( 0 )  
   , ALL_COMPS_SAME_LOCATION( true ) 
   , SET_BC_VALUES_ALLOWED( false )      	   
   , transproj_interpolation( 0 )
   , synchronization_ready( false )
   , IS_PROTO( false )           
{
   MAC_LABEL( "FV_DiscreteField:: FV_DiscreteField" ) ;

   // Create set of boundary conditions
   FV_BoundaryCondition* pbc = NULL ;
   size_t nbcs = DIM == 2 ? 8 : 26 ;
   V_BCS = new vector< FV_BoundaryCondition* >( nbcs+2, pbc );
   for( size_t i=0;i<nbcs;++i) 
   {
     SET_OF_BCS.push_back( pbc );
     SET_OF_BCS.back() = FV_BoundaryCondition:: create( this, NB_COMPS, i+2 ) ;
     (*V_BCS)[i+2] = SET_OF_BCS.back();
   }

   // Set periodicity
   PERIODIC = new boolVector( DIM, false ) ; 	
   *PERIODIC = *(PRIMARY_GRID->get_periodic_directions()) ;
   for (size_t i=0;i<DIM && !PERIODIC_SHIFT;++i)
     if ( (*PERIODIC)(i) ) 
     {
       FV_SHIFT_TRIPLET mst;
       mst.i = 0 ;
       mst.j = 0 ;
       mst.k = 0 ;
       PERIODIC_SHIFT = new vector< FV_SHIFT_TRIPLET >( NB_COMPS, mst ) ;     
     }
   
   MAC_CHECK_POST( owner() == a_owner ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;
}




//----------------------------------------------------------------------
FV_DiscreteField:: FV_DiscreteField( 
	std::string const& a_type )
//----------------------------------------------------------------------   
   : MAC_Object( plugins_map() )
   , PRIMARY_GRID( 0 )
   , FNAME( "" )
   , FDISCRETIZATION( "" )
   , ID( 0 )
   , NB_COMPS( 0 )
   , NB_CELLS( 0 )  
   , STO_DEPTH( 0 )
   , VALUES( 0 ) 
   , V_BCS( 0 )   
   , UNK_LOCAL_NUMBERING( 0 )
   , UNK_GLOBAL_NUMBERING( 0 )
   , NB_LOCAL_UNKNOWNS( 0 )
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC( 0 )
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE( 0 ) 
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_PERIODIC_BUFFERZONE( 0 )         
   , NB_LOCAL_DOF( 0 )
   , NB_GLOBAL_UNKNOWNS( 0 ) 
   , INITIALIZER( 0 )
   , CTX( 0 )
   , COORDS( 0 ) 
   , global_max_index( 0 )
   , local_max_index_in_global( 0 )    
   , local_min_index_in_global( 0 )
   , local_dof_number( 0 )
   , max_index_unknown_handled_by_proc( 0 )
   , min_index_unknown_handled_by_proc( 0 )
   , max_index_unknown_on_proc( 0 )
   , min_index_unknown_on_proc( 0 )       
   , global_main_coordinates( 0 )
   , local_main_coordinates( 0 )
   , local_cell_size( 0 )    
   , DOFcolors( 0 )
   , DOFstatus( 0 ) 
   , PERIODIC( 0 )
   , PERIODIC_SHIFT( 0 )         
   , PARAVIEW_FNAME( "" )
   , LOCATION( "" )
   , NB_DOF_POSTPROCESSING_PER_COMP( 0 ) 
   , ALL_COMPS_SAME_LOCATION( true )
   , SET_BC_VALUES_ALLOWED( false )
   , transproj_interpolation( 0 ) 
   , synchronization_ready( false )                      	   
   , IS_PROTO( true )    
{
   MAC_LABEL( "FV_DiscreteField:: FV_DiscreteField" ) ;

   plugins_map()->register_item( a_type, this ) ;

   MAC_CHECK_POST( owner() == plugins_map() ) ;
   MAC_CHECK_POST( is_a_prototype() ) ;   
}




//----------------------------------------------------------------------
FV_DiscreteField:: FV_DiscreteField( MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type )
//----------------------------------------------------------------------   
   : MAC_Object( a_owner )
   , PRIMARY_GRID( a_primary_mesh )
   , FNAME( a_name )
   , FDISCRETIZATION( a_type )
   , ID( NB_INSTANCES++ )
   , NB_COMPS( 0 )
   , NB_CELLS( 0 )  
   , STO_DEPTH( 0 )
   , VALUES( 0 ) 
   , V_BCS( 0 )   
   , UNK_LOCAL_NUMBERING( 0 )
   , UNK_GLOBAL_NUMBERING( 0 )
   , NB_LOCAL_UNKNOWNS( 0 )
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC( 0 )
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE( 0 )
   , NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_PERIODIC_BUFFERZONE( 0 )          
   , NB_LOCAL_DOF( 0 )
   , NB_GLOBAL_UNKNOWNS( 0 ) 
   , INITIALIZER( 0 )
   , CTX( 0 )
   , COORDS( 0 ) 
   , global_max_index( 0 )
   , local_max_index_in_global( 0 )    
   , local_min_index_in_global( 0 )
   , local_dof_number( 0 )
   , max_index_unknown_handled_by_proc( 0 )
   , min_index_unknown_handled_by_proc( 0 ) 
   , max_index_unknown_on_proc( 0 )
   , min_index_unknown_on_proc( 0 )      
   , global_main_coordinates( 0 )
   , local_main_coordinates( 0 )
   , local_cell_size( 0 )    
   , DOFcolors( 0 )
   , DOFstatus( 0 ) 
   , PERIODIC( 0 )
   , PERIODIC_SHIFT( 0 )      
   , PARAVIEW_FNAME( "" )
   , LOCATION( "" )
   , NB_DOF_POSTPROCESSING_PER_COMP( 0 ) 
   , ALL_COMPS_SAME_LOCATION( true ) 
   , SET_BC_VALUES_ALLOWED( false ) 
   , transproj_interpolation( 0 )
   , synchronization_ready( false )                  	   
   , IS_PROTO( false )    
{
   MAC_LABEL( "FV_DiscreteField:: FV_DiscreteField" ) ;

   MAC_CHECK_POST( owner() == a_owner ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;   
}





//----------------------------------------------------------------------
FV_DiscreteField:: ~FV_DiscreteField( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: ~FV_DiscreteField" ) ;

   if ( V_BCS )
   {
     V_BCS->clear() ;
     delete V_BCS ;
   }
   
   if ( VALUES )
   { 
     vector< vector< doubleArray3D > >::iterator il;
     for (il=VALUES->begin();il!=VALUES->end();il++)
       il->clear();
     VALUES->clear();
     delete VALUES;
   }

   if ( UNK_LOCAL_NUMBERING )
   {
     UNK_LOCAL_NUMBERING->clear() ;
     delete UNK_LOCAL_NUMBERING ;
   }

   if ( UNK_GLOBAL_NUMBERING )
   {
     UNK_GLOBAL_NUMBERING->clear() ;
     delete UNK_GLOBAL_NUMBERING ;
   }
   
   if ( global_main_coordinates )
   {
     global_max_index->clear();
     delete global_max_index;
     
     local_min_index_in_global->clear();
     delete local_min_index_in_global;
     
     local_max_index_in_global->clear();
     delete local_max_index_in_global;
     
     local_dof_number->clear();
     delete local_dof_number; 
       
     min_index_unknown_handled_by_proc->clear();
     delete min_index_unknown_handled_by_proc;
            
     max_index_unknown_handled_by_proc->clear();
     delete max_index_unknown_handled_by_proc; 
     
     min_index_unknown_on_proc->clear();
     delete min_index_unknown_on_proc;
            
     max_index_unknown_on_proc->clear();
     delete max_index_unknown_on_proc;      
               
     vector < vector< doubleVector > >::iterator iv ;
     for (iv=global_main_coordinates->begin();
     	iv!=global_main_coordinates->end();iv++)
       iv->clear();
     global_main_coordinates->clear();
     delete global_main_coordinates;
     
     for (iv=local_main_coordinates->begin();
     	iv!=local_main_coordinates->end();iv++)
       iv->clear();
     local_main_coordinates->clear();
     delete local_main_coordinates;
     
     for (iv=local_cell_size->begin();iv!=local_cell_size->end();iv++)
       iv->clear();
     local_cell_size->clear();
     delete local_cell_size;     
     
     vector < vector< size_t_vector > >::iterator ivv ;
     for (ivv=on_current_processor->begin();
     	ivv!=on_current_processor->end();ivv++)
       ivv->clear();
     on_current_processor->clear();
     delete on_current_processor;
     
     DOFcolors->clear();
     delete DOFcolors;
     
     DOFstatus->clear();
     delete DOFstatus; 
     
     delete PERIODIC;
     
     if ( PERIODIC_SHIFT )
     {
       PERIODIC_SHIFT->clear();
       delete PERIODIC_SHIFT;
     }                    
   } 
   
   if ( transproj_interpolation )
   {
     transproj_interpolation->clear();
     delete transproj_interpolation ; 
   } 
   
   if ( synchronization_ready )
   {
     list< vector< vector< FV_TRIPLET > > >::iterator ilv;
     vector< vector< FV_TRIPLET > >::iterator iv;
     for (ilv=halozone_received.begin();ilv!=halozone_received.end();ilv++)
     {
       for (iv=ilv->begin();iv!=ilv->end();iv++)       
         iv->clear();
       ilv->clear();
     }
     halozone_received.clear();
     
     for (ilv=bufferzone_sent.begin();ilv!=bufferzone_sent.end();ilv++)
     {
       for (iv=ilv->begin();iv!=ilv->end();iv++)       
         iv->clear();
       ilv->clear();
     }
     bufferzone_sent.clear();
     
     synchronization_MPI_rank_neighbors.clear(); 
     
     list< double* >::iterator ild;
     for (ild=halozone_received_data.begin();ild!=halozone_received_data.end();
     	ild++)
       delete [] (*ild);
     halozone_received_data.clear();
     for (ild=bufferzone_sent_data.begin();ild!=bufferzone_sent_data.end();
     	ild++)
       delete [] (*ild);  
     bufferzone_sent_data.clear();
     
     halozone_received_data_size.clear();  
     bufferzone_sent_data_size.clear();         
   }    
}




//----------------------------------------------------------------------
FV_DiscreteField*
FV_DiscreteField:: create_clone( MAC_Object* a_owner,
      		std::string const& name_of_new_field ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: create_clone" ) ;

   FV_DiscreteField const* proto =
      static_cast<FV_DiscreteField const*>( 
      	plugins_map()->item( FDISCRETIZATION ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;

   FV_DiscreteField* result = proto->create_clone_replica( a_owner,
   	PRIMARY_GRID, name_of_new_field, FDISCRETIZATION ) ;
   
   result->DIM = DIM ;
   result->NB_COMPS = NB_COMPS ;
   result->NB_CELLS = NB_CELLS ;
   result->STO_DEPTH = STO_DEPTH ;

   doubleArray3D work_double( 1, 1, 1, 0. ) ;
   vector< doubleArray3D > vwork_double( (*VALUES)[0].size(), work_double ) ;
   result->VALUES = new vector< vector< doubleArray3D > >( VALUES->size(),
   	vwork_double ) ;
   for (size_t lev=0;lev<STO_DEPTH;++lev)
     for (size_t comp=0;comp<NB_COMPS;++comp)
       (*result->VALUES)[lev][comp] = (*VALUES)[lev][comp] ;
 
   
   FV_BoundaryCondition* pbc = NULL ;
   result->V_BCS = new vector< FV_BoundaryCondition* >( V_BCS->size(), 
   	pbc );
   size_t bcnum = 0 ;	
   for (list< FV_BoundaryCondition* >::const_iterator il=SET_OF_BCS.begin();
   	il!=SET_OF_BCS.end();il++,++bcnum)
   {
     pbc = (*il)->create_clone( result );
     result->SET_OF_BCS.push_back( pbc ) ;
     (*result->V_BCS)[bcnum+2] = result->SET_OF_BCS.back() ;
   }

   intArray3D work_int( 1, 1, 1, 0 ) ; 
   result->UNK_LOCAL_NUMBERING = new vector< intArray3D >(
   	UNK_LOCAL_NUMBERING->size(), work_int ) ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
     (*result->UNK_LOCAL_NUMBERING)[comp] = (*UNK_LOCAL_NUMBERING)[comp] ;
  
   longLongIntArray3D work_llint( 1, 1, 1, 0 ) ;
   result->UNK_GLOBAL_NUMBERING = new vector< longLongIntArray3D >(
   	UNK_GLOBAL_NUMBERING->size(), work_llint ) ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
     (*result->UNK_GLOBAL_NUMBERING)[comp] = (*UNK_GLOBAL_NUMBERING)[comp] ;

   result->NB_LOCAL_UNKNOWNS = NB_LOCAL_UNKNOWNS ;
   result->NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC = 
   	NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC ;
   result->NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE =
   	NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE ;
   result->NB_LOCAL_DOF = NB_LOCAL_DOF ;
   result->NB_GLOBAL_UNKNOWNS = NB_GLOBAL_UNKNOWNS ;

   if ( INITIALIZER != 0 )
   {
     MAC_ContextSimple* c = MAC_ContextSimple::create( result ) ; 
     result->COORDS = MAC_DoubleVector::create( c, doubleVector(0) ) ;    
     c->extend( MAC_Variable::object( "DV_X" ), result->COORDS ) ;    
     result->CTX = c ;
     result->INITIALIZER = INITIALIZER->create_clone( result->CTX ) ;
   }
       
   size_t_vector work_sizet( 1, 0 ) ; 
   result->global_max_index = new vector< size_t_vector >( 
   	global_max_index->size(), work_sizet ) ;
   for (size_t i=0;i<global_max_index->size();++i)
     (*result->global_max_index)[i] = (*global_max_index)[i] ;
     
   result->local_max_index_in_global = new vector< size_t_vector >( 
   	local_max_index_in_global->size(), work_sizet ) ;
   for (size_t i=0;i<local_max_index_in_global->size();++i)
     (*result->local_max_index_in_global)[i] = (*local_max_index_in_global)[i] ;
     
   result->local_min_index_in_global = new vector< size_t_vector >( 
   	local_min_index_in_global->size(), work_sizet ) ;
   for (size_t i=0;i<local_min_index_in_global->size();++i)
     (*result->local_min_index_in_global)[i] = (*local_min_index_in_global)[i] ;
     
   result->local_dof_number = new vector< size_t_vector >( 
   	local_dof_number->size(), work_sizet ) ;
   for (size_t i=0;i<local_dof_number->size();++i)
     (*result->local_dof_number)[i] = (*local_dof_number)[i] ;     
     
   result->max_index_unknown_handled_by_proc = new vector< size_t_vector >( 
   	max_index_unknown_handled_by_proc->size(), work_sizet ) ;
   for (size_t i=0;i<max_index_unknown_handled_by_proc->size();++i)
     (*result->max_index_unknown_handled_by_proc)[i] = 
     	(*max_index_unknown_handled_by_proc)[i] ;     
     
   result->min_index_unknown_handled_by_proc = new vector< size_t_vector >( 
   	min_index_unknown_handled_by_proc->size(), work_sizet ) ;
   for (size_t i=0;i<min_index_unknown_handled_by_proc->size();++i)
     (*result->min_index_unknown_handled_by_proc)[i] = 
     	(*min_index_unknown_handled_by_proc)[i] ;
	
   result->max_index_unknown_on_proc = new vector< size_t_vector >( 
   	max_index_unknown_on_proc->size(), work_sizet ) ;
   for (size_t i=0;i<max_index_unknown_on_proc->size();++i)
     (*result->max_index_unknown_on_proc)[i] = 
     	(*max_index_unknown_on_proc)[i] ;     
     
   result->min_index_unknown_on_proc = new vector< size_t_vector >( 
   	min_index_unknown_on_proc->size(), work_sizet ) ;
   for (size_t i=0;i<min_index_unknown_on_proc->size();++i)
     (*result->min_index_unknown_on_proc)[i] = 
     	(*min_index_unknown_on_proc)[i] ;	     
     
   doubleVector work_vd( 1 , 0. ) ; 
   vector< doubleVector > work_vvd( (*global_main_coordinates)[0].size(), 
   	work_vd ) ;     
   result->global_main_coordinates = new vector < vector< doubleVector > >
   	( global_main_coordinates->size(), work_vvd ) ;
   for (size_t i=0;i<global_main_coordinates->size();++i)
     for (size_t j=0;j<(*global_main_coordinates)[i].size();++j)
       (*result->global_main_coordinates)[i][j] = 
       			(*global_main_coordinates)[i][j] ;
			
   result->local_main_coordinates = new vector < vector< doubleVector > >
   	( local_main_coordinates->size(), work_vvd ) ;
   for (size_t i=0;i<local_main_coordinates->size();++i)
     for (size_t j=0;j<(*local_main_coordinates)[i].size();++j)
       (*result->local_main_coordinates)[i][j] = 
       			(*local_main_coordinates)[i][j] ;			
     
   result->local_cell_size = new vector < vector< doubleVector > >
   	( local_cell_size->size(), work_vvd ) ;
   for (size_t i=0;i<local_cell_size->size();++i)
     for (size_t j=0;j<(*local_cell_size)[i].size();++j)
       (*result->local_cell_size)[i][j] = 
       			(*local_cell_size)[i][j] ;

   vector< size_t_vector > work_vsizet( (*on_current_processor)[0].size(), 
   	work_sizet ) ;
   result->on_current_processor = new vector < vector< size_t_vector > >
   	( on_current_processor->size(), work_vsizet ) ;
   for (size_t i=0;i<on_current_processor->size();++i)
     for (size_t j=0;j<(*on_current_processor)[i].size();++j)
       (*result->on_current_processor)[i][j] = 
       			(*on_current_processor)[i][j] ;	
			
   result->DOFcolors = new vector< intArray3D >(
   	DOFcolors->size(), work_int ) ;
   for (size_t i=0;i<DOFcolors->size();++i)
     (*result->DOFcolors)[i] = (*DOFcolors)[i] ;

   result->DOFstatus = new vector< intArray3D >(
   	DOFstatus->size(), work_int ) ;
   for (size_t i=0;i<DOFstatus->size();++i)
     (*result->DOFstatus)[i] = (*DOFstatus)[i] ;  

   result->PERIODIC = new boolVector( PERIODIC->size(), false ) ;
   *result->PERIODIC = *PERIODIC ;
   
   if ( PERIODIC_SHIFT )
   {
     FV_SHIFT_TRIPLET mst;
     mst.i = 0 ;
     mst.j = 0 ;
     mst.k = 0 ;
     result->PERIODIC_SHIFT = new vector< FV_SHIFT_TRIPLET >( 
     	PERIODIC_SHIFT->size(), mst ) ;
     for (size_t i=0;i<PERIODIC_SHIFT->size();++i)
     {
       (*result->PERIODIC_SHIFT)[i].i = (*PERIODIC_SHIFT)[i].i ; 
       (*result->PERIODIC_SHIFT)[i].j = (*PERIODIC_SHIFT)[i].j ; 
       (*result->PERIODIC_SHIFT)[i].k = (*PERIODIC_SHIFT)[i].k ;
     }   	
   }
        
   result->NB_DOF_POSTPROCESSING_PER_COMP = NB_DOF_POSTPROCESSING_PER_COMP ; 
   
   result->ALL_COMPS_SAME_LOCATION = ALL_COMPS_SAME_LOCATION ; 
   
   result->SET_BC_VALUES_ALLOWED = SET_BC_VALUES_ALLOWED ;       

   if ( transproj_interpolation )
   {
     size_t_vector work( 1, 0 ) ;
     result->transproj_interpolation = new vector< size_t_vector >( NB_COMPS,
     	work ) ;
     for (size_t comp=0;comp<NB_COMPS;++comp) 
       (*result->transproj_interpolation)[comp] = 
       		(*transproj_interpolation)[comp] ;
   }   
   
   // Copy synchronization features
   if ( synchronization_ready )
   {
     list< vector< vector< FV_TRIPLET > > >::const_iterator ihr, ibs;
     for (ihr=halozone_received.begin();ihr!=halozone_received.end();ihr++)
       result->halozone_received.push_back( *ihr );
       
     for (ibs=bufferzone_sent.begin();ibs!=bufferzone_sent.end();ibs++)
       result->halozone_received.push_back( *ibs ); 
       
     list< size_t >::const_iterator il;
     for (il=synchronization_MPI_rank_neighbors.begin();
     	il!=synchronization_MPI_rank_neighbors.end();il++)
       result->synchronization_MPI_rank_neighbors.push_back( *il );
     for (il=halozone_received_data_size.begin();
     	il!=halozone_received_data_size.end();il++)
       result->halozone_received_data_size.push_back( *il );
     for (il=bufferzone_sent_data_size.begin();
     	il!=bufferzone_sent_data_size.end();il++)
       result->bufferzone_sent_data_size.push_back( *il ); 
       
     double* pnull = NULL;
     list< double* >::const_iterator ild;
     list< double* >::iterator ildr;
     size_t nhal = halozone_received_data.size();
     for (size_t i=0;i<nhal;++i) 
       result->halozone_received_data.push_back( pnull );
     for (ild=halozone_received_data.begin(),
     	ildr=result->halozone_received_data.begin(),
	il=halozone_received_data_size.begin();	
	ild!=halozone_received_data.end();
     	ild++,il++,ildr++)
     {
       size_t vecsize = *il;
       (*ildr) = new double[vecsize];
       for (size_t i=0;i<vecsize;++i) (*ildr)[i] = (*ild)[i];  
     }
     
     size_t nbuf = bufferzone_sent.size();
     for (size_t i=0;i<nbuf;++i) 
       result->bufferzone_sent_data.push_back( pnull );
     for (ild=bufferzone_sent_data.begin(),
     	ildr=result->bufferzone_sent_data.begin(),
	il=bufferzone_sent_data_size.begin();	
	ild!=bufferzone_sent_data.end();
     	ild++,il++,ildr++)
     {
       size_t vecsize = *il;
       (*ildr) = new double[vecsize];
       for (size_t i=0;i<vecsize;++i) (*ildr)[i] = (*ild)[i];  
     }     
   }   
			
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->storage_depth() == STO_DEPTH ) ;
   MAC_CHECK_POST( result->nb_components() == NB_COMPS ) ;
   MAC_CHECK_POST( result->primary_grid() == PRIMARY_GRID ) ;
   MAC_CHECK_POST( result->discretization_type() == FDISCRETIZATION ) ;      
   MAC_CHECK_POST( result->paraview_location() == LOCATION ) ;      
   return( result ) ;
}




//----------------------------------------------------------------------
bool 
FV_DiscreteField:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: is_a_prototype" ) ;
   return( IS_PROTO ) ;   
}




//----------------------------------------------------------------------
MAC_ObjectRegister* 
FV_DiscreteField:: plugins_map( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: plugins_map" ) ;
   static MAC_ObjectRegister* result =
      MAC_ObjectRegister::create( MAC_Root::object(),
                                  "FV_DiscreteField descendant" ) ;
   return( result ) ;
   
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: nb_objects( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: nb_objects" ) ;
   return( NB_INSTANCES ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: id_number( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: id_number" ) ;
   
   size_t result = ID ;

   MAC_CHECK_POST( result < nb_objects() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
string const& 
FV_DiscreteField:: name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: name" ) ;

   return( FNAME ) ;
}




//----------------------------------------------------------------------
string const& 
FV_DiscreteField:: discretization_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: discretization_type" ) ;

   return( FDISCRETIZATION ) ;
}




//----------------------------------------------------------------------
string const& 
FV_DiscreteField:: paraview_location( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: paraview_location" ) ;

   return( LOCATION ) ;
}




//----------------------------------------------------------------------
FV_Mesh const*
FV_DiscreteField:: primary_grid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: primary_grid" ) ;

   return( PRIMARY_GRID ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: nb_cells( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: nb_cells" ) ;

   return( NB_CELLS ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: nb_components( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: nb_components" ) ;

   return( NB_COMPS ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: storage_depth( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: storage_depth" ) ;

   return( STO_DEPTH ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: read_BCs( MAC_ModuleExplorer const* exp, 
      	FV_DomainBuilder const* DB ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: read_BCs" ) ;
   
   // Read BCs from data file
   if ( exp->has_module( "DOFs_imposed_value" ) )
   {
      MAC_ModuleExplorer* sexp =
                exp->create_subexplorer( 0, "DOFs_imposed_value" ) ;
      sexp->start_module_iterator() ;
      for( ; sexp->is_valid_module() ; sexp->go_next_module() )
      {
	 MAC_ModuleExplorer* sse = sexp->create_subexplorer( 0 ) ;
         string bc_color = sse->string_data( "color" ) ;
	 if ( FV_DomainBuilder::does_primary_color_exist( bc_color ) )
	 {	 
	   size_t color_number = FV_DomainBuilder:: get_color_number( 
	   	bc_color ) ;
	   (*V_BCS)[color_number]->read_dirichlet_BC( sse, DIM, FNAME ) ;
	 }
	 else
	 {
	   stringVector const* primary_colors = 
	   	DB->get_colors_from_macro_color( bc_color ) ;
		
	   for (size_t i=0;i<primary_colors->size();++i)
	   {
	     size_t color_number = FV_DomainBuilder:: get_color_number( 
	   	(*primary_colors)(i) ) ;
	     (*V_BCS)[color_number]->read_dirichlet_BC( sse, DIM, FNAME ) ;	   
	   }  	
	 }	
         sse->destroy() ;
      }
      sexp->destroy() ; sexp = 0 ;
   }     
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_freeDOFonBC_update_features( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_freeDOFonBC_update_features" ) ;
   
   // Set free_DOF_on_boundary update features
   // Rule for main colors
   //   If main color is Neumann without unknowns, use neighboring values in
   //   domain
   // Rules for blended colors
   //   1. if one sub main color is Dirichlet, use this Dirichlet value
   //   2. if more than one sub main color are Dirichlet, use a weighted value 
   //   of all Dirichlet values
   //   3. if all sub main colors are Neumann and none have unknowns, sum the 
   //   shift Mac triplets
   //   4. if at least one of the sub main colors is Neumann with unknowns, 
   //   apply rules 1 and 2 (with Dirichlet value replaced by Neumann with 
   //   unknowns value) 
   list<size_t>::const_iterator ic; 
   for (list< FV_BoundaryCondition* >::iterator il=SET_OF_BCS.begin();
   	il!=SET_OF_BCS.end();il++)
   {
     size_t bc_color = (*il)->get_color_ID() ;
     for (size_t comp=0;comp<NB_COMPS;++comp)
     {
       if ( (*il)->is_neumann( comp ) && !(*il)->has_unknown( comp ) )
       {
         if ( FV_DomainBuilder::is_main_color( bc_color ) )
	 {
	   FV_SHIFT_TRIPLET const* mst = 
	     FV_DomainBuilder::get_shift_MacTriplet( bc_color ) ;
	   (*il)->set_shift_MacTriplet( comp, mst, 1. ) ;
	 }
         else
	 {
	   list<size_t> const* submaincolors = 
	     FV_DomainBuilder::get_sub_main_color_ids( bc_color ) ;
	     
           FV_SHIFT_TRIPLET summst;
	   summst.i = 0 ;
	   summst.j = 0 ;
	   summst.k = 0 ;	     	     
	   for (ic=submaincolors->begin();ic!=submaincolors->end();ic++)
	   {
	     FV_SHIFT_TRIPLET const* mst = 
	     		FV_DomainBuilder::get_shift_MacTriplet( *ic ) ;
	     summst.i += mst->i ;
	     summst.j += mst->j ;		 
	     summst.k += mst->k ;		  
	   }
	   	     
	   size_t nDirichlet = 0 ;
	   for (ic=submaincolors->begin();ic!=submaincolors->end();ic++)
	     nDirichlet += (*V_BCS)[*ic]->is_dirichlet( comp ) ;
	     
	   if ( nDirichlet )
	   {
	     double weight = 1. / double(nDirichlet) ;     	     
	     for (ic=submaincolors->begin();ic!=submaincolors->end();ic++)
	       if ( (*V_BCS)[*ic]->is_dirichlet( comp ) )
	       {
	         FV_SHIFT_TRIPLET const* mst = 
	     		FV_DomainBuilder::get_shift_MacTriplet( *ic ) ;
	         (*il)->add_shift_MacTriplet( comp, 
		 	summst.i - mst->i,
			summst.j - mst->j,
			summst.k - mst->k,
			weight ) ;
	       }	       
	   }
	   else
	   {
             size_t nNeumannWithUnknown = 0 ;
	     for (ic=submaincolors->begin();ic!=submaincolors->end();ic++)
	       nNeumannWithUnknown += (*V_BCS)[*ic]->has_unknown( comp ) ;
	       
	     if ( nNeumannWithUnknown )
	     {
	       double weight = 1. / double(nNeumannWithUnknown) ;     	     
	       for (ic=submaincolors->begin();ic!=submaincolors->end();ic++)
	         if ( (*V_BCS)[*ic]->has_unknown( comp ) )
		 {
	           FV_SHIFT_TRIPLET const* mst = 
	     		FV_DomainBuilder::get_shift_MacTriplet( *ic ) ;
		   (*il)->add_shift_MacTriplet( comp, 
		 	summst.i - mst->i,
			summst.j - mst->j,
			summst.k - mst->k,
			weight ) ;
		 }	     
	     }
	     else
	       (*il)->set_shift_MacTriplet( comp, &summst, 1. ) ;
	   }  
	 }
       }
     }
   }         
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_imposed_DOF_values( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_imposed_DOF_values" ) ;
   MAC_CHECK_PRE( VALUES != 0 ) ;
   MAC_CHECK_PRE( SET_OF_BCS.size() != 0 ) ;   
   
   bool set_bc_values = SET_BC_VALUES_ALLOWED ;
   SET_BC_VALUES_ALLOWED = true ;
   for (list< FV_BoundaryCondition* >::iterator il=SET_OF_BCS.begin();
   	il!=SET_OF_BCS.end();il++)
     (*il)->set_imposed_DOF_values( this );
   SET_BC_VALUES_ALLOWED = set_bc_values ; 
}      




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_neumann_DOF_values( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_neumann_DOF_values" ) ;
   MAC_CHECK_PRE( VALUES != 0 ) ;
   MAC_CHECK_PRE( SET_OF_BCS.size() != 0 ) ;   

   bool set_bc_values = SET_BC_VALUES_ALLOWED ;
   SET_BC_VALUES_ALLOWED = true ;
   for (size_t level=0;level<STO_DEPTH;++level)
     for ( list< FV_BoundaryCondition* >::iterator il=SET_OF_BCS.begin();
   	il!=SET_OF_BCS.end();il++)
       (*il)->set_free_DOF_values( this, level ) ;       
   SET_BC_VALUES_ALLOWED = set_bc_values ; 
}      




//----------------------------------------------------------------------
void 
FV_DiscreteField:: extract_unknown_DOFs_value( size_t level, 
	LA_SeqVector* vec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: extract_unknown_DOFs_value" ) ;
   MAC_CHECK_PRE( level < STO_DEPTH ) ;
   MAC_CHECK_PRE( vec != 0 ) ;
   MAC_CHECK_PRE( UNK_LOCAL_NUMBERING->size() == NB_COMPS );
   MAC_CHECK_PRE( UNK_GLOBAL_NUMBERING->size() == NB_COMPS );   
   MAC_CHECK_PRE( vec->nb_rows() >= NB_LOCAL_UNKNOWNS ) ;

   vec->nullify() ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t nelem0 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(0) ;
     size_t nelem1 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(1) ;
     size_t nelem2 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(2) ; 
     for (size_t i=0;i<nelem0;++i)
       for (size_t j=0;j<nelem1;++j) 
         for (size_t k=0;k<nelem2;++k)
	 {
	   int locnum = (*UNK_LOCAL_NUMBERING)[comp](i,j,k) ;
	   if ( locnum != -1 ) 
	     vec->set_item( locnum, (*VALUES)[level][comp](i,j,k) );
	 }
   }
   vec->synchronize() ;

   MAC_CHECK_POST( vec->is_synchronized() ) ;
}




//----------------------------------------------------------------------
void 
FV_DiscreteField:: update_free_DOFs_value( size_t level, 
	LA_SeqVector const* vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: update_free_DOFs_value" ) ;
   MAC_CHECK_PRE( level < STO_DEPTH ) ;
   MAC_CHECK_PRE( vec != 0 ) ;
   MAC_CHECK_PRE( UNK_LOCAL_NUMBERING->size() == NB_COMPS );
   MAC_CHECK_PRE( UNK_GLOBAL_NUMBERING->size() == NB_COMPS );   
   MAC_CHECK_PRE( vec->nb_rows() >= NB_LOCAL_UNKNOWNS ) ;  
   
   for( size_t comp=0; comp<NB_COMPS; ++comp )
   {
     size_t nelem0 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(0) ;
     size_t nelem1 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(1) ;
     size_t nelem2 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(2) ; 
     for( size_t i=0; i<nelem0; ++i )
       for( size_t j=0; j<nelem1; ++j ) 
         for( size_t k=0; k<nelem2; ++k )
	 {
	   int locnum = (*UNK_LOCAL_NUMBERING)[comp](i,j,k) ;
	   if( locnum != -1 )
	     (*VALUES)[level][comp](i,j,k) = vec->item( locnum ) ;
	 }
   }

   bool set_bc_values = SET_BC_VALUES_ALLOWED;
   SET_BC_VALUES_ALLOWED = true;
   for( list< FV_BoundaryCondition* >::iterator il=SET_OF_BCS.begin();
   	il!=SET_OF_BCS.end(); il++ )
     (*il)->set_free_DOF_values( this, level );
   SET_BC_VALUES_ALLOWED = set_bc_values;

}




//----------------------------------------------------------------------
void 
FV_DiscreteField:: add_to_free_DOFs_value( size_t level, 
	LA_SeqVector const* vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: add_to_free_DOFs_value" ) ;
   MAC_CHECK_PRE( level < STO_DEPTH ) ;
   MAC_CHECK_PRE( vec != 0 ) ;
   MAC_CHECK_PRE( UNK_LOCAL_NUMBERING->size() == NB_COMPS );
   MAC_CHECK_PRE( UNK_GLOBAL_NUMBERING->size() == NB_COMPS );   
   MAC_CHECK_PRE( vec->nb_rows() >= NB_LOCAL_UNKNOWNS ) ;  
   
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t nelem0 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(0) ;
     size_t nelem1 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(1) ;
     size_t nelem2 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(2) ; 
     for (size_t i=0;i<nelem0;++i)
       for (size_t j=0;j<nelem1;++j) 
         for (size_t k=0;k<nelem2;++k)
	 {
	   int locnum = (*UNK_LOCAL_NUMBERING)[comp](i,j,k) ;
	   if ( locnum != -1 )
	     (*VALUES)[level][comp](i,j,k) += vec->item( locnum ) ; 
	 }
   }

   bool set_bc_values = SET_BC_VALUES_ALLOWED ;
   SET_BC_VALUES_ALLOWED = true ;   
   for ( list< FV_BoundaryCondition* >::iterator il=SET_OF_BCS.begin();
   	il!=SET_OF_BCS.end();il++)
     (*il)->set_free_DOF_values( this, level ) ;
   SET_BC_VALUES_ALLOWED =  set_bc_values ;            
}
	 



//----------------------------------------------------------------------
size_t 
FV_DiscreteField:: nb_global_unknowns( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: nb_global_unknowns" ) ;
 
   return ( NB_GLOBAL_UNKNOWNS ) ;
}




//----------------------------------------------------------------------
size_t 
FV_DiscreteField:: nb_local_unknowns_handled_by_proc( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: nb_local_unknowns_handled_by_proc" ) ;
 
   return ( NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC ) ;
}




//----------------------------------------------------------------------
size_t 
FV_DiscreteField:: nb_local_unknowns( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: nb_local_unknowns" ) ;
 
   return ( NB_LOCAL_UNKNOWNS ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_DOF_colors( size_t component ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_DOF_colors" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ;
   
   // Boundary conditions
   if ( DIM == 2 )
   {  
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
       {
         if ( i + (*local_min_index_in_global)[component](0) == 0 )
	 {
	   if ( j + (*local_min_index_in_global)[component](1) == 0 )
	     // Bottom left
	     (*DOFcolors)[component](i,j,0) = FV_BC_BOTTOM_LEFT;
	   else if ( j + (*local_min_index_in_global)[component](1) == 
	   	(*global_max_index)[component](1) )
	     // Top left
	     (*DOFcolors)[component](i,j,0) = FV_BC_TOP_LEFT;
	   else
	     // Left
	     (*DOFcolors)[component](i,j,0) = FV_BC_LEFT;
	 }
	 else if ( i + (*local_min_index_in_global)[component](0) == 
	 		(*global_max_index)[component](0) )
	 {
	   if ( j + (*local_min_index_in_global)[component](1) == 0 )
	     // Bottom right
	     (*DOFcolors)[component](i,j,0) = FV_BC_BOTTOM_RIGHT;
	   else if ( j + (*local_min_index_in_global)[component](1) == 
	   		(*global_max_index)[component](1) )
	     // Top right
	     (*DOFcolors)[component](i,j,0) = FV_BC_TOP_RIGHT;
	   else
	     // Right
	     (*DOFcolors)[component](i,j,0) = FV_BC_RIGHT; 
	 }
	 else
	 {
	   if ( j + (*local_min_index_in_global)[component](1) == 0 )
	     // Bottom
	     (*DOFcolors)[component](i,j,0) = FV_BC_BOTTOM;
	   else if ( j + (*local_min_index_in_global)[component](1) == 
	   		(*global_max_index)[component](1) )
	     // Top
	     (*DOFcolors)[component](i,j,0) = FV_BC_TOP;	 
	 }
       }	         
   }
   else
   {
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j)
	 for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
         {
           if ( i + (*local_min_index_in_global)[component](0) == 0 )
	   {
	     if ( j + (*local_min_index_in_global)[component](1) == 0 )
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind bottom left
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_BOTTOM_LEFT;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front bottom left
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_BOTTOM_LEFT;
	       else
	         // Bottom left
	         (*DOFcolors)[component](i,j,k) = FV_BC_BOTTOM_LEFT;
	     }
	     else if ( j + (*local_min_index_in_global)[component](1) == 
	     			(*global_max_index)[component](1) )
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind top left
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_TOP_LEFT;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front top left
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_TOP_LEFT;
	       else
	         // Top left
	         (*DOFcolors)[component](i,j,k) = FV_BC_TOP_LEFT;
	     }
	     else
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind left
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_LEFT;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front left
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_LEFT;
	       else
	         // Left
		 (*DOFcolors)[component](i,j,k) = FV_BC_LEFT;	 
	     }
	   }
	   else if ( i + (*local_min_index_in_global)[component](0) == 
	   			(*global_max_index)[component](0) )
	   {
	     if ( j + (*local_min_index_in_global)[component](1) == 0 )
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind bottom right
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_BOTTOM_RIGHT;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front bottom right
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_BOTTOM_RIGHT;
	       else
	         // Bottom right
	         (*DOFcolors)[component](i,j,k) = FV_BC_BOTTOM_RIGHT;
	     }
	     else if ( j + (*local_min_index_in_global)[component](1) == 
	     			(*global_max_index)[component](1) )
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind top right
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_TOP_RIGHT;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front top right
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_TOP_RIGHT;
	       else
	         // Top right
	         (*DOFcolors)[component](i,j,k) = FV_BC_TOP_RIGHT;
	     }
	     else
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind right
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_RIGHT;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front right
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_RIGHT;
	       else
	         // right
		 (*DOFcolors)[component](i,j,k) = FV_BC_RIGHT;	 
	     }
	   }
	   else
	   {
	     if ( j + (*local_min_index_in_global)[component](1) == 0 )
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind bottom
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_BOTTOM;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front bottom
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_BOTTOM;
	       else
	         // Bottom
	         (*DOFcolors)[component](i,j,k) = FV_BC_BOTTOM;
	     }
	     else if ( j + (*local_min_index_in_global)[component](1) == 
	     			(*global_max_index)[component](1) )
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind top
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND_TOP;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front top
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT_TOP;
	       else
	         // Top
	         (*DOFcolors)[component](i,j,k) = FV_BC_TOP;
	     }
	     else
	     {
	       if ( k + (*local_min_index_in_global)[component](2) == 0 )
	         // Behind
	         (*DOFcolors)[component](i,j,k) = FV_BC_BEHIND;
	       else if ( k + (*local_min_index_in_global)[component](2) == 
	       			(*global_max_index)[component](2) )
	         // Front
	         (*DOFcolors)[component](i,j,k) = FV_BC_FRONT;	 
	     }
	   }
	 }	   	   	     	     
   }  
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_PERIODIC_default( size_t component, 
      	size_t_vector const* periodic_depth ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_PERIODIC_default" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ;
   MAC_CHECK_PRE( IMPLIES( ALL_COMPS_SAME_LOCATION, component == 0 ) ) ;
   
   size_t compBC = ALL_COMPS_SAME_LOCATION == false ? component : NB_COMPS ;
   size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;    
   double tol = PRIMARY_GRID->get_smallest_grid_size() / 10. ;    
      
   if ( periodic_depth )
   {
     if ( DIM == 2 )
     {
       if ( (*PERIODIC)(0) )
       {
         // Left plane
	 for (size_t j=0;j<(*local_dof_number)[component](1);++j)
           if ( (*DOFcolors)[component](0,j,0) == FV_BC_LEFT )
	   {
	     for (size_t i=0;i<(*periodic_depth)(0);++i)
	       (*DOFcolors)[component](i,j,0) = FV_BC_PERIODIC ;
	     if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][1](j), 1, tol ) )
	       (*V_BCS)[FV_BC_LEFT]->add_MacTriplet( compBC, 
	     	(*periodic_depth)(0) - 1, j, 0 ) ;  
	   }
	   else if ( (*DOFcolors)[component](0,j,0) == FV_BC_BOTTOM_LEFT )
	     (*DOFcolors)[component](0,j,0) = FV_BC_BOTTOM ;
	   else if ( (*DOFcolors)[component](0,j,0) == FV_BC_TOP_LEFT )
	     (*DOFcolors)[component](0,j,0) = FV_BC_TOP ;
	     
         (*V_BCS)[FV_BC_LEFT]->set_periodic();
	 (*V_BCS)[FV_BC_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_TOP_LEFT]->set_none();
        
         // Right plane
	 for (size_t j=0;j<(*local_dof_number)[component](1);++j)
           if ( (*DOFcolors)[component]((*local_dof_number)[component](0)-1,
	   	j,0) == FV_BC_RIGHT )
	   {
	     for (size_t i=0;i<(*periodic_depth)(0);++i)
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1-i,j,0) = FV_BC_PERIODIC ;
             if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][1](j), 1, tol ) )
	       (*V_BCS)[FV_BC_RIGHT]->add_MacTriplet( compBC, 
	     	(*local_dof_number)[component](0) - (*periodic_depth)(0), j, 
		0 ) ; 
	   }
	   else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,0) == FV_BC_BOTTOM_RIGHT )
	     (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,0) = FV_BC_BOTTOM ;
	   else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,0) == FV_BC_TOP_RIGHT )
	     (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,0) = FV_BC_TOP ;	
	     
         (*V_BCS)[FV_BC_RIGHT]->set_periodic();
	 (*V_BCS)[FV_BC_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_TOP_RIGHT]->set_none();         	     
       }
       		
       if ( (*PERIODIC)(1) )
       {
         // Bottom plane
	 for (size_t i=0;i<(*local_dof_number)[component](0);++i)
           if ( (*DOFcolors)[component](i,0,0) == FV_BC_BOTTOM )
	   {
	     for (size_t j=0;j<(*periodic_depth)(1);++j)
	       (*DOFcolors)[component](i,j,0) = FV_BC_PERIODIC ;
	     if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][0](i), 0, tol ) )	     
	       (*V_BCS)[FV_BC_BOTTOM]->add_MacTriplet( compBC, 
	     	i, (*periodic_depth)(1) - 1, 0 ) ; 	       
	   }
	   else if ( (*DOFcolors)[component](i,0,0) == FV_BC_BOTTOM_LEFT )
	     (*DOFcolors)[component](i,0,0) = FV_BC_LEFT ;
	   else if ( (*DOFcolors)[component](i,0,0) == FV_BC_BOTTOM_RIGHT )
	     (*DOFcolors)[component](i,0,0) = FV_BC_RIGHT ;	   
	     
         (*V_BCS)[FV_BC_BOTTOM]->set_periodic();
	 (*V_BCS)[FV_BC_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BOTTOM_RIGHT]->set_none();    

         // Top plane
	 for (size_t i=0;i<(*local_dof_number)[component](0);++i)
           if ( (*DOFcolors)[component](i,(*local_dof_number)[component](1)-1,0)
	   	 == FV_BC_TOP )
	   {
	     for (size_t j=0;j<(*periodic_depth)(1);++j)
	       (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1-j,0) = FV_BC_PERIODIC ;
	     if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][0](i), 0, tol ) )
	       (*V_BCS)[FV_BC_TOP]->add_MacTriplet( compBC, 
	     	i, (*local_dof_number)[component](1) - (*periodic_depth)(1), 
		0 ) ; 	     
	   }
	   else if ( (*DOFcolors)[component](i,
	   	(*local_dof_number)[component](1)-1,0) == FV_BC_TOP_LEFT )
	    (*DOFcolors)[component](i,(*local_dof_number)[component](1)-1,0) =
			FV_BC_LEFT ;
	   else if ( (*DOFcolors)[component](i,
	   	(*local_dof_number)[component](1)-1,0) == FV_BC_TOP_RIGHT )
	    (*DOFcolors)[component](i,(*local_dof_number)[component](1)-1,0) =
			FV_BC_RIGHT ;

         (*V_BCS)[FV_BC_TOP]->set_periodic();
	 (*V_BCS)[FV_BC_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_TOP_RIGHT]->set_none(); 	     
       }
     }
     else
     {
       if ( (*PERIODIC)(0) )
       {
         // Left plane
	 for (size_t j=0;j<(*local_dof_number)[component](1);++j)
           for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
	     if ( (*DOFcolors)[component](0,j,k) == FV_BC_LEFT )
	     {
	       for (size_t i=0;i<(*periodic_depth)(0);++i)
	         (*DOFcolors)[component](i,j,k) = FV_BC_PERIODIC ;
	       if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][1](j), 1, tol ) 
		&& PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][2](k), 2, tol ) )
	         (*V_BCS)[FV_BC_LEFT]->add_MacTriplet( compBC, 
	     	(*periodic_depth)(0) - 1, j, k ) ;		 
	     }
	     else if ( (*DOFcolors)[component](0,j,k) == FV_BC_BOTTOM_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_BOTTOM ;
	     else if ( (*DOFcolors)[component](0,j,k) == FV_BC_TOP_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_TOP ;
	     else if ( (*DOFcolors)[component](0,j,k) == FV_BC_BEHIND_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_BEHIND ;
	     else if ( (*DOFcolors)[component](0,j,k) == FV_BC_FRONT_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_FRONT ;
	     else if ( (*DOFcolors)[component](0,j,k) 
	     	== FV_BC_BEHIND_BOTTOM_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_BEHIND_BOTTOM ;
	     else if ( (*DOFcolors)[component](0,j,k) 
	     	== FV_BC_BEHIND_TOP_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_BEHIND_TOP ;
	     else if ( (*DOFcolors)[component](0,j,k) 
	     	== FV_BC_FRONT_BOTTOM_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_FRONT_BOTTOM ;
	     else if ( (*DOFcolors)[component](0,j,k) == FV_BC_FRONT_TOP_LEFT )
	     	(*DOFcolors)[component](0,j,k) = FV_BC_FRONT_TOP ;
		
         (*V_BCS)[FV_BC_LEFT]->set_periodic();
	 (*V_BCS)[FV_BC_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP_LEFT]->set_none();

         // Right plane
         for (size_t j=0;j<(*local_dof_number)[component](1);++j)
	   for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
             if ( (*DOFcolors)[component]((*local_dof_number)[component](0)-1,
	   	j,k) == FV_BC_RIGHT )
	     {
	       for (size_t i=0;i<(*periodic_depth)(0);++i)
	         (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1-i,j,k) = FV_BC_PERIODIC ;
	       if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][1](j), 1, tol ) 
		&& PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][2](k), 2, tol ) )
	       (*V_BCS)[FV_BC_RIGHT]->add_MacTriplet( compBC, 
	     	(*local_dof_number)[component](0) - (*periodic_depth)(0), j, 
		k ) ; 		
	     }
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_BOTTOM_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_BOTTOM ;
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_TOP_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_TOP ;
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_BEHIND_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_BEHIND ;
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_FRONT_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_FRONT ;
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_BEHIND_BOTTOM_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_BEHIND_BOTTOM;
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_BEHIND_TOP_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_BEHIND_TOP ;
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_FRONT_BOTTOM_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_FRONT_BOTTOM ;
	     else if ( (*DOFcolors)[component]((*local_dof_number)[component](0)
	   	-1,j,k) == FV_BC_FRONT_TOP_RIGHT )
	       (*DOFcolors)[component](
	       	(*local_dof_number)[component](0)-1,j,k) = FV_BC_FRONT_TOP ;
		
         (*V_BCS)[FV_BC_RIGHT]->set_periodic();
	 (*V_BCS)[FV_BC_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_TOP_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP_RIGHT]->set_none();
       }
		
       if ( (*PERIODIC)(1) )
       {
         // Bottom plane
	 for (size_t i=0;i<(*local_dof_number)[component](0);++i)
	   for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
             if ( (*DOFcolors)[component](i,0,k) == FV_BC_BOTTOM )
	     {
	       for (size_t j=0;j<(*periodic_depth)(1);++j)
	         (*DOFcolors)[component](i,j,k) = FV_BC_PERIODIC ;
	       if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][0](i), 0, tol ) 
		&& PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][2](k), 2, tol ) )	       
	         (*V_BCS)[FV_BC_BOTTOM]->add_MacTriplet( compBC, 
	     	i, (*periodic_depth)(1) - 1, k ) ; 		 
	     }
	     else if ( (*DOFcolors)[component](i,0,k) == FV_BC_BOTTOM_LEFT )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_LEFT ;
	     else if ( (*DOFcolors)[component](i,0,k) == FV_BC_BOTTOM_RIGHT )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_RIGHT ;
	     else if ( (*DOFcolors)[component](i,0,k) == FV_BC_BEHIND_BOTTOM )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_BEHIND ;
	     else if ( (*DOFcolors)[component](i,0,k) == FV_BC_FRONT_BOTTOM )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_FRONT ;
	     else if ( (*DOFcolors)[component](i,0,k) 
	     	== FV_BC_BEHIND_BOTTOM_LEFT )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_BEHIND_LEFT ;
	     else if ( (*DOFcolors)[component](i,0,k) 
	     	== FV_BC_BEHIND_BOTTOM_RIGHT )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_BEHIND_RIGHT ;
	     else if ( (*DOFcolors)[component](i,0,k) 
	     	== FV_BC_FRONT_BOTTOM_LEFT )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_FRONT_LEFT ;
	     else if ( (*DOFcolors)[component](i,0,k) 
	     	== FV_BC_FRONT_BOTTOM_RIGHT )
	     	(*DOFcolors)[component](i,0,k) = FV_BC_FRONT_RIGHT ;

         (*V_BCS)[FV_BC_BOTTOM]->set_periodic();
	 (*V_BCS)[FV_BC_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM_RIGHT]->set_none();		

         // Top plane
	 for (size_t i=0;i<(*local_dof_number)[component](0);++i)
           for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
	     if ( (*DOFcolors)[component](i,
	     	(*local_dof_number)[component](1)-1,k) == FV_BC_TOP )
	     {
	       for (size_t j=0;j<(*periodic_depth)(1);++j)
	         (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1-j,k) = FV_BC_PERIODIC ;
	       if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][0](i), 0, tol ) 
		&& PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][2](k), 2, tol ) )
	         (*V_BCS)[FV_BC_TOP]->add_MacTriplet( compBC, 
	     	i, (*local_dof_number)[component](1) - (*periodic_depth)(1), 
		k ) ; 		
	     }
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_TOP_LEFT )
	         (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_LEFT ;
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_TOP_RIGHT )
	     (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_RIGHT ;
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_BEHIND_TOP )
	     (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_BEHIND ;
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_FRONT_TOP )
	     (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_FRONT ;
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_BEHIND_TOP_LEFT )
	     (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_BEHIND_LEFT;
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_BEHIND_TOP_RIGHT )
	     (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_BEHIND_RIGHT ;
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_FRONT_TOP_LEFT )
	     (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_FRONT_LEFT ;
	     else if ( (*DOFcolors)[component](i,(*local_dof_number)[component]
	     	(1)-1,k) == FV_BC_FRONT_TOP_RIGHT )
	     (*DOFcolors)[component](
	       	i,(*local_dof_number)[component](1)-1,k) = FV_BC_FRONT_RIGHT ;
		
         (*V_BCS)[FV_BC_TOP]->set_periodic();
	 (*V_BCS)[FV_BC_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_TOP_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP_RIGHT]->set_none();	
       }
		
       if ( (*PERIODIC)(2) )
       {
         // Behind plane
	 for (size_t i=0;i<(*local_dof_number)[component](0);++i)
	   for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
             if ( (*DOFcolors)[component](i,j,0) == FV_BC_BEHIND )
	     {
	       for (size_t k=0;k<(*periodic_depth)(2);++k)
	         (*DOFcolors)[component](i,j,k) = FV_BC_PERIODIC ;
	       if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][0](i), 0, tol ) 
		&& PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][1](j), 1, tol ) )	       
	         (*V_BCS)[FV_BC_BEHIND]->add_MacTriplet( compBC, 
	     	i, j, (*periodic_depth)(2) - 1 ) ;	        	 
	     }
	     else if ( (*DOFcolors)[component](i,j,0) == FV_BC_BEHIND_LEFT )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_LEFT ;
	     else if ( (*DOFcolors)[component](i,j,0) == FV_BC_BEHIND_RIGHT )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_RIGHT ;
	     else if ( (*DOFcolors)[component](i,j,0) == FV_BC_BEHIND_BOTTOM )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_BOTTOM ;
	     else if ( (*DOFcolors)[component](i,j,0) == FV_BC_BEHIND_TOP )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_TOP ;
	     else if ( (*DOFcolors)[component](i,j,0) 
	     	== FV_BC_BEHIND_BOTTOM_LEFT )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_BOTTOM_LEFT ;
	     else if ( (*DOFcolors)[component](i,j,0) 
	     	== FV_BC_BEHIND_BOTTOM_RIGHT )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_BOTTOM_RIGHT ;
	     else if ( (*DOFcolors)[component](i,j,0) 
	     	== FV_BC_BEHIND_TOP_LEFT )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_TOP_LEFT ;
	     else if ( (*DOFcolors)[component](i,j,0) 
	     	== FV_BC_BEHIND_TOP_RIGHT )
	     	(*DOFcolors)[component](i,j,0) = FV_BC_TOP_RIGHT ;
		
         (*V_BCS)[FV_BC_BEHIND]->set_periodic();
	 (*V_BCS)[FV_BC_BEHIND_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_BEHIND_TOP_RIGHT]->set_none();

         // Front plane
	 for (size_t i=0;i<(*local_dof_number)[component](0);++i)
	   for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
             if ( (*DOFcolors)[component](i,j,
	     	(*local_dof_number)[component](2)-1) == FV_BC_FRONT )
	     {
	       for (size_t k=0;k<(*periodic_depth)(2);++k)
	         (*DOFcolors)[component](i,j,
		(*local_dof_number)[component](2)-1-k) = FV_BC_PERIODIC ;
	       if ( PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][0](i), 0, tol ) 
		&& PRIMARY_GRID->is_in_main_domain( 
	     	(*local_main_coordinates)[locNum][1](j), 1, tol ) )
	         (*V_BCS)[FV_BC_FRONT]->add_MacTriplet( compBC, 
	     	i, j, (*local_dof_number)[component](2) 
			- (*periodic_depth)(2) ) ; 		
	     }
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) == FV_BC_FRONT_LEFT )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_LEFT ;
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) == FV_BC_FRONT_RIGHT )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_RIGHT ;
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) == FV_BC_FRONT_BOTTOM )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_BOTTOM ;
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) == FV_BC_FRONT_TOP )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_TOP ;
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) 
		== FV_BC_FRONT_BOTTOM_LEFT )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_BOTTOM_LEFT;
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) 
		== FV_BC_FRONT_BOTTOM_RIGHT )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_BOTTOM_RIGHT ;
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) == FV_BC_FRONT_TOP_LEFT )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_TOP_LEFT ;
	     else if ( (*DOFcolors)[component](i,j,
	   	(*local_dof_number)[component](2)-1) == FV_BC_FRONT_TOP_RIGHT )
	       (*DOFcolors)[component](
	       	i,j,(*local_dof_number)[component](2)-1) = FV_BC_TOP_RIGHT ;

         (*V_BCS)[FV_BC_FRONT]->set_periodic();
	 (*V_BCS)[FV_BC_FRONT_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_BOTTOM_RIGHT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP_LEFT]->set_none();
	 (*V_BCS)[FV_BC_FRONT_TOP_RIGHT]->set_none();		
       }	
     }       
   }        
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_DOF_status( size_t component ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_DOF_status" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ) ;
    
   // Halozone dof
   if ( DIM == 2 )
   {  
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         if ( (*on_current_processor)[component][0](i) == 2 
	 	|| (*on_current_processor)[component][1](j) == 2 )
	   (*DOFstatus)[component](i,j,0) = FV_DOF_HALOZONE ;
   }
   else
   {
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
         if ( (*on_current_processor)[component][0](i) == 2
	 	|| (*on_current_processor)[component][1](j) == 2 
		|| (*on_current_processor)[component][2](k) == 2 )
	   (*DOFstatus)[component](i,j,k) = FV_DOF_HALOZONE ; 
   }      

   // Buffer zone
   if ( DIM == 2 )
   {  
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         if ( (*on_current_processor)[component][0](i) == 1 
	 	|| (*on_current_processor)[component][1](j) == 1 )
	   if ( (*DOFstatus)[component](i,j,0) != FV_DOF_HALOZONE ) 
	     (*DOFstatus)[component](i,j,0) = FV_DOF_BUFFERZONE ;
   }
   else
   {
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
           if ( (*on_current_processor)[component][0](i) == 1 
	 	|| (*on_current_processor)[component][1](j) == 1
	 	|| (*on_current_processor)[component][2](k) == 1 )
	     if ( (*DOFstatus)[component](i,j,k) != FV_DOF_HALOZONE ) 
	       (*DOFstatus)[component](i,j,k) = FV_DOF_BUFFERZONE ;
   }
   
   // Halozone periodic dof
   if ( DIM == 2 )
   {  
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         if ( (*on_current_processor)[component][0](i) == 4 
	 	|| (*on_current_processor)[component][1](j) == 4 )
	   (*DOFstatus)[component](i,j,0) = FV_DOF_PERIODIC_HALOZONE ;
   }
   else
   {
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
           if ( (*on_current_processor)[component][0](i) == 4
	 	|| (*on_current_processor)[component][1](j) == 4 
		|| (*on_current_processor)[component][2](k) == 4 )
	     (*DOFstatus)[component](i,j,k) = FV_DOF_PERIODIC_HALOZONE ; 
   }
   
   // Buffer periodic dof
   if ( DIM == 2 )
   {  
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         if ( (*on_current_processor)[component][0](i) == 3 
	 	|| (*on_current_processor)[component][1](j) == 3 )
	   if ( (*DOFstatus)[component](i,j,0) != FV_DOF_HALOZONE )
	   { 
	     if ( (*DOFstatus)[component](i,j,0) == FV_DOF_ONPROC ) 
	       (*DOFstatus)[component](i,j,0) = FV_DOF_PERIODIC_BUFFERZONE ;
	     else if ( (*DOFstatus)[component](i,j,0) == FV_DOF_BUFFERZONE )
	       (*DOFstatus)[component](i,j,0) = 
	       		FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE ; 
	   }
   }
   else
   {
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
           if ( (*on_current_processor)[component][0](i) == 3 
	 	|| (*on_current_processor)[component][1](j) == 3
	 	|| (*on_current_processor)[component][2](k) == 3 )
	     if ( (*DOFstatus)[component](i,j,k) != FV_DOF_HALOZONE ) 
	     {
	       if ( (*DOFstatus)[component](i,j,k) == FV_DOF_ONPROC ) 
	         (*DOFstatus)[component](i,j,k) = FV_DOF_PERIODIC_BUFFERZONE ;
	       else if ( (*DOFstatus)[component](i,j,k) == FV_DOF_BUFFERZONE )
	         (*DOFstatus)[component](i,j,k) = 
	       		FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE ;
	     }
   } 
   
   // Buffer - Buffer periodic dof
   if ( DIM == 2 )
   {  
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         if ( (*on_current_processor)[component][0](i) == 5 
	 	|| (*on_current_processor)[component][1](j) == 5 )
	   if ( (*DOFstatus)[component](i,j,0) != FV_DOF_HALOZONE 
	   	&& (*DOFstatus)[component](i,j,0) != FV_DOF_PERIODIC_HALOZONE )
	     (*DOFstatus)[component](i,j,0) = 
	       		FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE ; 
   }
   else
   {
     for (size_t i=0;i<(*local_dof_number)[component](0);++i)
       for (size_t j=0;j<(*local_dof_number)[component](1);++j) 
         for (size_t k=0;k<(*local_dof_number)[component](2);++k) 
           if ( (*on_current_processor)[component][0](i) == 5 
	 	|| (*on_current_processor)[component][1](j) == 5
	 	|| (*on_current_processor)[component][2](k) == 5 )
	     if ( (*DOFstatus)[component](i,j,k) != FV_DOF_HALOZONE 
	   	&& (*DOFstatus)[component](i,j,k) != FV_DOF_PERIODIC_HALOZONE )
	       (*DOFstatus)[component](i,j,k) = 
	       		FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE ;
   }                      
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: build_field_numbering( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: build_field_numbering" ) ;
   MAC_CHECK_PRE( DOFcolors != 0 );
   MAC_CHECK_PRE( DOFstatus != 0 );   

   size_t locNum = 0 ;
        
   // Allocate numbering arrays
   intArray3D work(1,1,1);
   UNK_LOCAL_NUMBERING = new vector< intArray3D >( NB_COMPS, work ) ;
   longLongIntArray3D llwork(1,1,1);   
   UNK_GLOBAL_NUMBERING = new vector< longLongIntArray3D >( NB_COMPS, llwork ) ;   
   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     (*UNK_LOCAL_NUMBERING)[comp].re_initialize( 
     	(*local_dof_number)[locNum](0),
   	(*local_dof_number)[locNum](1),
	DIM == 2 ? 1 : (*local_dof_number)[locNum](2), 
	-1 );
     (*UNK_GLOBAL_NUMBERING)[comp].re_initialize( 
     	(*local_dof_number)[locNum](0),
   	(*local_dof_number)[locNum](1),
	DIM == 2 ? 1 : (*local_dof_number)[locNum](2), 
	-1 );	
   }
   
   // Set the number of local dofs
   NB_LOCAL_DOF = 0 ;
   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     NB_LOCAL_DOF += (*UNK_LOCAL_NUMBERING)[locNum].index_bound(0)
   	* (*UNK_LOCAL_NUMBERING)[locNum].index_bound(1)
   	* (*UNK_LOCAL_NUMBERING)[locNum].index_bound(2) ;
   }
	
   // Set the number of local unknowns   
   NB_LOCAL_UNKNOWNS = 0 ; 
   NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC = 0 ;
   NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE = 0 ; 
   NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_PERIODIC_BUFFERZONE = 0; 	

   bool numbering = false ;
   int status = 0 ;
   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[locNum](2);
     for (size_t i=0;i<(*local_dof_number)[locNum](0);++i)
       for (size_t j=0;j<(*local_dof_number)[locNum](1);++j) 
         for (size_t k=0;k<kmax;++k)
	 {
           numbering = false ;	   
	   if ( (*DOFcolors)[locNum](i,j,k) < 2 ) numbering = true ;
	   else if ( (*V_BCS)[(*DOFcolors)[locNum](i,j,k)]
	   	->has_unknown( comp ) )
	     numbering = true ;
	   
	   if ( numbering )
	   {
	     status = (*DOFstatus)[locNum](i,j,k) ;
	     (*UNK_LOCAL_NUMBERING)[comp](i,j,k) = NB_LOCAL_UNKNOWNS++;
	     set_min_max_indices_unknown_on_proc( i, j, k, locNum ) ;
	     if ( status != FV_DOF_HALOZONE 
	     	&& status != FV_DOF_PERIODIC_HALOZONE )
	     {
	       NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC++ ;
	       set_min_max_indices_unknown_handled_by_proc( i, j, k, locNum ) ;
	       if ( status == FV_DOF_BUFFERZONE 
	       	|| status == FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE ) 
	         NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE++ ;
	       if ( status == FV_DOF_PERIODIC_BUFFERZONE 
	       	|| status == FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE ) 
	         NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_PERIODIC_BUFFERZONE++ ;
	     }
	   }
	 }
   }

   // Share the local number of unknowns on processor between processes
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   size_t nb_ranks = macCOMM->nb_ranks();
   size_t rank = macCOMM->rank();
   intVector nb_local_unk_per_rank( nb_ranks ) ;
   macCOMM->all_gather( NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC, 
   	nb_local_unk_per_rank ) ;

   // Number of global unknowns
   NB_GLOBAL_UNKNOWNS = 0 ;
   for (size_t i=0;i<nb_ranks;++i) 
     NB_GLOBAL_UNKNOWNS += nb_local_unk_per_rank(i) ;

   // Global numbering for unknowns on processor
   size_t starting_number = 0 ;
   if ( rank )
     for (size_t i=0;i<rank;++i) starting_number += nb_local_unk_per_rank(i) ;
   
   size_t local_numbering_on_proc = 0 ;
   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[locNum](2);
     for (size_t i=0;i<(*local_dof_number)[locNum](0);++i)
       for (size_t j=0;j<(*local_dof_number)[locNum](1);++j) 
         for (size_t k=0;k<kmax;++k)
	 {
	   numbering = false ;
	   if ( (*DOFcolors)[locNum](i,j,k) < 2 ) numbering = true ;
	   else if ( (*V_BCS)[(*DOFcolors)[locNum](i,j,k)]
	   	->has_unknown( comp ) )
	     numbering = true ;	   
	   
	   if ( numbering )
	     if ( (*DOFstatus)[locNum](i,j,k) != FV_DOF_HALOZONE 
	     	&& (*DOFstatus)[locNum](i,j,k) != FV_DOF_PERIODIC_HALOZONE )
	     {
               (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) = starting_number
	     	+ local_numbering_on_proc;
               local_numbering_on_proc++; 
	     }
	 }   
   }
   
   // Global numbering for unknowns in halozone
   size_t width = 5 ;
   bool send_numbering = false ;
   longLongIntVector global_num_comm( 
   	NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE*width, 0 );
   size_t positionIndex = 0 ;

   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[locNum](2);
     for (size_t i=0;i<(*local_dof_number)[locNum](0);++i)
       for (size_t j=0;j<(*local_dof_number)[locNum](1);++j) 
         for (size_t k=0;k<kmax;++k)
	 { 
           send_numbering = false ;
	   if ( (*DOFcolors)[locNum](i,j,k) < 2 ) send_numbering = true ;
	   else if ( (*V_BCS)[(*DOFcolors)[locNum](i,j,k)]
	   	->has_unknown( comp ) )
	     send_numbering = true ;	
	   
	   if ( send_numbering )
	     if ( (*DOFstatus)[locNum](i,j,k) == FV_DOF_BUFFERZONE
	     	|| (*DOFstatus)[locNum](i,j,k) 
			== FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE )
	     {
               global_num_comm(positionIndex) = i + 
	     	(*local_min_index_in_global)[locNum](0) ;
               global_num_comm(positionIndex+1) = j + 
	     	(*local_min_index_in_global)[locNum](1) ;
               global_num_comm(positionIndex+2) = DIM == 2 ? 0 : 
	     	k + (*local_min_index_in_global)[locNum](2) ;
	       global_num_comm(positionIndex+3) = comp ;	     
               global_num_comm(positionIndex+4) = 
	     	(*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ;	     
	       positionIndex += width ;
	     }
	 }
   }         
   
   // Send global numbering in bufferzone to other processors
   // for update if they own the unknown in their own halozone 
   list<size_t> MPI_neighbors_world =
   	*PRIMARY_GRID->get_MPI_neighbors_ranks() ;
   MPI_neighbors_world.push_back( rank ) ;
   MPI_neighbors_world.sort();
   list<size_t>::iterator r,k;
   for (r=MPI_neighbors_world.begin();r!=MPI_neighbors_world.end();r++)
   {
     if ( *r == rank )
     {
       for (k=MPI_neighbors_world.begin();k!=MPI_neighbors_world.end();k++)
         if ( *k != rank )
	   macCOMM->send( *k, global_num_comm ) ; 
     }
     else
     {
       longLongIntVector global_num_from_r(1) ;
       macCOMM->receive( *r, global_num_from_r ) ;
       size_t vecsize = global_num_from_r.size(), ig, jg, kg = 0, comp, gnum; 
       
       for (positionIndex=0;positionIndex<vecsize; )
       {
         ig = global_num_from_r(positionIndex);
         jg = global_num_from_r(positionIndex+1);
         kg = global_num_from_r(positionIndex+2);
	 comp = global_num_from_r(positionIndex+3);
	 gnum = global_num_from_r(positionIndex+4);
	 
         locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
	 if ( is_global_triplet_local_DOF( ig, jg, kg, locNum ) ) 
	   (*UNK_GLOBAL_NUMBERING)[comp](
	   	ig-(*local_min_index_in_global)[locNum](0),
	   	jg-(*local_min_index_in_global)[locNum](1),
		DIM == 2 ? 0 : kg - (*local_min_index_in_global)[locNum](2) )
		= gnum ;
	 	 	 
         positionIndex += width ;
       }
     }
   }   
   
   // Global numbering in periodic halozone
   if ( PERIODIC_SHIFT ) build_periodic_numbering() ;   
} 




//----------------------------------------------------------------------
void
FV_DiscreteField:: build_periodic_numbering( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: build_periodic_numbering" ) ;
   MAC_CHECK_PRE( DOFcolors != 0 );
   MAC_CHECK_PRE( DOFstatus != 0 ); 

   list<FV_TRIPLET> ijk ;
   FV_TRIPLET mt ;
   mt.i = 0 ;
   mt.j = 0 ;   
   mt.k = 0 ; 
   list<size_t> comps, globalnums ;
   size_t locNum ; 

   // Store periodic buffer zone unknown numbering in lists
   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[locNum](2);
     for (size_t i=0;i<(*local_dof_number)[locNum](0);++i)
       for (size_t j=0;j<(*local_dof_number)[locNum](1);++j) 
         for (size_t k=0;k<kmax;++k)
	 {            
	   if ( ( (*DOFstatus)[locNum](i,j,k) == FV_DOF_PERIODIC_BUFFERZONE
		|| (*DOFstatus)[locNum](i,j,k) 
			== FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE )
		&& (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) != -1 ) 
	   {
             // X periodicity
	     if ( (*PERIODIC)(0) )
	     {
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			+ (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			- (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }	        	       
	     }
	     
             // Y periodicity
	     if ( (*PERIODIC)(1) )
	     {
               mt.i = i + (*local_min_index_in_global)[locNum](0) ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }	        	       
	     }
	     
             // X-Y biperiodicity
	     if ( (*PERIODIC)(0) && (*PERIODIC)(1) )
	     {
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			+ (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			+ (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			- (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			- (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ) ;
	       }	       
	     }
	     
	     // 3D cases
	     if ( DIM == 3 )
	     {
               // Z periodicity
	       if ( (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }	        	       
	       }
	       
               // X-Z biperiodicity
	       if ( (*PERIODIC)(0) && (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }		 	        	       
	       }
	       
               // Y-Z biperiodicity
	       if ( (*PERIODIC)(1) && (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }		 	        	       
	       }
	       
               // X-Y-Z triperiodicity
	       if ( (*PERIODIC)(0) && (*PERIODIC)(1) && (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   globalnums.push_back( (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) );
	         }		 		 	        	       
	       }	       	       	       	     
	     }	     	     
	   }
	 }
   }   

   // Transfer the lists into a single vector
   size_t width = 5 ;
   longLongIntVector global_num_comm( ijk.size()*width, 0 );
   size_t positionIndex = 0 ;
   list<FV_TRIPLET>::iterator imac = ijk.begin();
   list<size_t>::iterator icomp = comps.begin(),iglobnum;
   
   for (iglobnum=globalnums.begin();iglobnum!=globalnums.end();iglobnum++,
   	icomp++,imac++)
   {
      global_num_comm(positionIndex) = imac->i ;
      global_num_comm(positionIndex+1) = imac->j ;
      global_num_comm(positionIndex+2) = imac->k ;
      global_num_comm(positionIndex+3) = *icomp ;	     
      global_num_comm(positionIndex+4) = *iglobnum ;	     
      positionIndex += width ;     
   }

   // Send global numbering in periodic bufferzone to other processors
   // for update if they own the unknown in their own periodic halozone
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   size_t rank = macCOMM->rank();
   size_t ig = 0, jg = 0, kg = 0, comp, gnum ;   
   list<size_t> MPI_periodic_neighbors_world =
   	*PRIMARY_GRID->get_MPI_periodic_neighbors_ranks() ;
   MPI_periodic_neighbors_world.push_back( rank ) ;
   MPI_periodic_neighbors_world.sort();
   list<size_t>::iterator r,k;
   for (r=MPI_periodic_neighbors_world.begin();
   	r!=MPI_periodic_neighbors_world.end();r++)
   {
     if ( *r == rank )
     {
       for (k=MPI_periodic_neighbors_world.begin();
       	k!=MPI_periodic_neighbors_world.end();k++)
         if ( *k != rank )
	   macCOMM->send( *k, global_num_comm ) ; 
     }
     else
     {
       longLongIntVector global_num_from_r(1) ;
       macCOMM->receive( *r, global_num_from_r ) ;
       size_t vecsize = global_num_from_r.size() ; 
       
       for (positionIndex=0;positionIndex<vecsize; )
       {
         ig = global_num_from_r(positionIndex);
         jg = global_num_from_r(positionIndex+1);
         kg = global_num_from_r(positionIndex+2);
	 comp = global_num_from_r(positionIndex+3);
	 gnum = global_num_from_r(positionIndex+4);
	 
         locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
	 if ( is_global_triplet_local_DOF( ig, jg, kg, locNum ) ) 
	   (*UNK_GLOBAL_NUMBERING)[comp](
	   	ig-(*local_min_index_in_global)[locNum](0),
	   	jg-(*local_min_index_in_global)[locNum](1),
		DIM == 2 ? 0 : kg - (*local_min_index_in_global)[locNum](2) )
		= gnum ;
	 	 	 
         positionIndex += width ;
       }
     }
   }   
   
   // Check for periodic unknown numbering on local processor
   for (positionIndex=0;positionIndex<global_num_comm.size(); )
   {
     ig = global_num_comm(positionIndex);
     jg = global_num_comm(positionIndex+1);
     kg = global_num_comm(positionIndex+2);
     comp = global_num_comm(positionIndex+3);
     gnum = global_num_comm(positionIndex+4);
	 
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     if ( is_global_triplet_local_DOF( ig, jg, kg, locNum ) ) 
       (*UNK_GLOBAL_NUMBERING)[comp](
	   	ig-(*local_min_index_in_global)[locNum](0),
	   	jg-(*local_min_index_in_global)[locNum](1),
		DIM == 2 ? 0 : kg - (*local_min_index_in_global)[locNum](2) )
		= gnum ;
	 	 	 
      positionIndex += width ;
   }             
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_synchronization_features( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_synchronization_features" ) ;

   synchronization_ready = true;
   
   size_t locNum = 0 ;
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   size_t rank = macCOMM->rank();

   // Fill buffer_features with the (i,j,k,comp) of dof handled by this proc
   // i.e. dof in the buffer zone
   size_t width = 4, positionIndex = 0, i, comp, ig, jg, kg, nijk, m ;
   bool send_features = false ;
   intVector buffer_features( 
   	NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE*width, 0 );
   for (comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[locNum](2);
     for (i=0;i<(*local_dof_number)[locNum](0);++i)
       for (size_t j=0;j<(*local_dof_number)[locNum](1);++j) 
         for (size_t k=0;k<kmax;++k)
	 { 
           send_features = false ;
	   if ( (*DOFcolors)[locNum](i,j,k) < 2 )
	     send_features = true ;
	   else if ( (*V_BCS)[(*DOFcolors)[locNum](i,j,k)]
	   	->has_unknown( comp ) )
	     send_features = true ;	
	   
	   if ( send_features )
	     if ( (*DOFstatus)[locNum](i,j,k) == FV_DOF_BUFFERZONE
	     	|| (*DOFstatus)[locNum](i,j,k) 
			== FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE )
	     {
               buffer_features(positionIndex) = i + 
	     	(*local_min_index_in_global)[locNum](0) ;
               buffer_features(positionIndex+1) = j + 
	     	(*local_min_index_in_global)[locNum](1) ;
               buffer_features(positionIndex+2) = DIM == 2 ? 0 : k + 
	     	(*local_min_index_in_global)[locNum](2) ;
	       buffer_features(positionIndex+3) = comp ;
	       positionIndex += width ;
	     }
	 }
   }

   // Send the bufferzone (i,j,k,comp) features to the neighboring procs
   // and process the received features. We store the (i,j,k,comp) that belong
   // to this proc in halozone_received
   list<size_t> MPI_neighbors_world =
   	*PRIMARY_GRID->get_MPI_neighbors_ranks() ;
   MPI_neighbors_world.push_back( rank ) ;
   MPI_neighbors_world.sort();   
   list<size_t>::iterator r, k, ils;
   FV_TRIPLET mt;
   mt.i = 0 ;
   mt.j = 0 ;
   mt.k = 0 ;
   list<FV_TRIPLET> ijk ;
   vector< list<FV_TRIPLET> > kept_ijk_per_comp( NB_COMPS, ijk );
   list<FV_TRIPLET>::iterator imt;
   
   for (r=MPI_neighbors_world.begin();r!=MPI_neighbors_world.end();r++)
   {
     if ( *r == rank )
     {
       for (k=MPI_neighbors_world.begin();k!=MPI_neighbors_world.end();k++)
         if ( *k != rank )
           macCOMM->send( *k, buffer_features ) ;
     }
     else
     {
       intVector received_buffer_features_from_r(1) ;
       macCOMM->receive( *r, received_buffer_features_from_r ) ;
       size_t vecsize = received_buffer_features_from_r.size() ; 

       for ( positionIndex=0;positionIndex<vecsize; )
       {
         ig = received_buffer_features_from_r(positionIndex);
         jg = received_buffer_features_from_r(positionIndex+1);
         kg = received_buffer_features_from_r(positionIndex+2);
	 comp = received_buffer_features_from_r(positionIndex+3);
	 
         locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;

	 if ( is_global_triplet_local_DOF( ig, jg, kg, locNum ) )
	 {
	   mt.i = ig - (*local_min_index_in_global)[locNum](0);
	   mt.j = jg - (*local_min_index_in_global)[locNum](1);
	   mt.k = DIM == 2 ? 0 : kg - 
	     	(*local_min_index_in_global)[locNum](2) ;
	   kept_ijk_per_comp[comp].push_back( mt );
	 }
	   
         positionIndex += width ;
       }
       
       // Transfer the vector of list into a vector of vector
       vector< FV_TRIPLET > vvv;
       vector< vector< FV_TRIPLET > > vec_kept_ijk_per_comp( NB_COMPS, vvv );
       for (comp=0;comp<NB_COMPS;comp++)
       {
         size_t vsize = kept_ijk_per_comp[comp].size();
	 vec_kept_ijk_per_comp[comp].reserve( vsize );
	 for (i=0;i<vsize;++i) vec_kept_ijk_per_comp[comp].push_back( mt ); 
       }    
       for (comp=0;comp<NB_COMPS;comp++)
       {
         size_t m = 0;
         for (imt=kept_ijk_per_comp[comp].begin();
	 	imt!=kept_ijk_per_comp[comp].end();imt++,++m)
	   vec_kept_ijk_per_comp[comp][m] = *imt;
       }	
	            
       // Store the vector of vector into halozone_received
       halozone_received.push_back( vec_kept_ijk_per_comp );
       
       // Store the MPI rank of the neighbor
       synchronization_MPI_rank_neighbors.push_back( *r );
       
       // Clear the content of the lists per comp
       for (comp=0;comp<NB_COMPS;comp++)
         kept_ijk_per_comp[comp].clear();
     }
   }  

   // Send back the halozone_received (i,j,k,comp) features to neighbors
   // such that they can create the corresponding bufferzone_sent (i,j,k,comp)
   // features
   intVector hz_features( 0, 0 );
   list< vector< vector< FV_TRIPLET > > >::iterator ihr, ibs;
   for (r=MPI_neighbors_world.begin();r!=MPI_neighbors_world.end();r++)
   {
     if ( *r == rank )
     {
       for ( k=MPI_neighbors_world.begin(),ihr=halozone_received.begin();
       		k!=MPI_neighbors_world.end(); k++ )
         if ( *k != rank )
	 {
           // Number of ijk to be sent
	   nijk = 0;
	   for (comp=0;comp<NB_COMPS;comp++)
	     nijk += (*ihr)[comp].size();
	   hz_features.resize( nijk * width ) ;

           // Fill the buffer
	   positionIndex = 0;
	   for (comp=0;comp<NB_COMPS;comp++)
	   {
	     nijk = (*ihr)[comp].size();
	     
	     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
	     
	     for (size_t m=0;m<nijk;++m)
	     {
	       hz_features(positionIndex) = (*ihr)[comp][m].i + 
	     	(*local_min_index_in_global)[locNum](0) ;			
               hz_features(positionIndex+1) = (*ihr)[comp][m].j + 
	     	(*local_min_index_in_global)[locNum](1) ;
               hz_features(positionIndex+2) = DIM == 2 ? 0 :
	       	(*ihr)[comp][m].k + (*local_min_index_in_global)[locNum](2) ;
	       hz_features(positionIndex+3) = comp ;
	       positionIndex += width ;	       
	     }
	   }

	   // Send message
	   macCOMM->send( *k, hz_features ) ;
	   ihr++;
	 }
     }
     else
     {
       intVector received_bufferzone_features_from_r(1) ;
       macCOMM->receive( *r, received_bufferzone_features_from_r ) ;
       size_t vecsize = received_bufferzone_features_from_r.size();
       vector<size_t> nijk_per_comp( NB_COMPS, 0 );
       
       // Get the number of ijk per comp
       for ( positionIndex=0; positionIndex<vecsize; )
       {
         comp = received_bufferzone_features_from_r(positionIndex+3);
	 nijk_per_comp[comp]++;
	 positionIndex += width;
       } 
       
       // Allocate bufferzone_sent for this neighboring proc
       vector< FV_TRIPLET > vvv;
       vector< vector< FV_TRIPLET > > vec_ijk_per_comp( NB_COMPS, vvv );
       for (comp=0;comp<NB_COMPS;comp++)
       {
	 vec_ijk_per_comp[comp].reserve( nijk_per_comp[comp] );
	 for (i=0;i<nijk_per_comp[comp];++i) 
	   vec_ijk_per_comp[comp].push_back( mt ); 
       }        
       
       positionIndex = 0;
       for (comp=0;comp<NB_COMPS;comp++)
       {
         for (m=0;m<nijk_per_comp[comp];++m)
	 {
           ig = received_bufferzone_features_from_r(positionIndex);
           jg = received_bufferzone_features_from_r(positionIndex+1);
           kg = received_bufferzone_features_from_r(positionIndex+2);
	 
           locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;

           vec_ijk_per_comp[comp][m].i = ig 
	 	- (*local_min_index_in_global)[locNum](0);
           vec_ijk_per_comp[comp][m].j = jg 
	 	- (*local_min_index_in_global)[locNum](1);		
           vec_ijk_per_comp[comp][m].k = DIM == 2 ? 0 : kg - 
	     	(*local_min_index_in_global)[locNum](2); 
	   
	   positionIndex += width ;	 
	 }
       }
                    
       // Store the vector of vector into bufferzone_sent
       bufferzone_sent.push_back( vec_ijk_per_comp );
     }
   }

   // If periodic, set the periodic synchronization features
   if ( PERIODIC_SHIFT ) set_periodic_synchronization_features() ; 
   
   // Allocate 1D arrays to send and receive data 
   size_t nneighbors = halozone_received.size();
   double* pdouble = NULL;
   for ( i=0;i<nneighbors;++i) 
   {
     halozone_received_data.push_back( pdouble );
     bufferzone_sent_data.push_back( pdouble );
   }
   
   list< double* >::iterator ild = halozone_received_data.begin();
   for ( ihr=halozone_received.begin();ihr!=halozone_received.end();
   	ihr++,ild++ )
   {
     nijk = 0;
     for (comp=0;comp<NB_COMPS;comp++) nijk += (*ihr)[comp].size();
     (*ild) = new double[nijk];
     halozone_received_data_size.push_back(nijk); 
   } 
   ild = bufferzone_sent_data.begin();
   for ( ibs=bufferzone_sent.begin();ibs!=bufferzone_sent.end();
   	ibs++,ild++ )
   {
     nijk = 0;
     for (comp=0;comp<NB_COMPS;comp++) nijk += (*ibs)[comp].size();
     (*ild) = new double[nijk];
     bufferzone_sent_data_size.push_back(nijk); 
   }
   
   
//    // Debug
//    size_t nb_ranks = macCOMM->nb_ranks() ;
//    size_t RANK = macCOMM->rank() ;
//    for (size_t i=0;i<nb_ranks;++i)
//    {
//      if ( i == RANK )
//      {       
//        cout << "Rank " << RANK << endl ;
//        cout << "   Buffer data" << endl;
//        for ( r=synchronization_MPI_rank_neighbors.begin(),
//   	ibs=bufferzone_sent.begin(),
// 	ils=bufferzone_sent_data_size.begin();
//   	r!=synchronization_MPI_rank_neighbors.end();
// 	r++,ibs++,ils++ )
//        {
//          cout << "   To rank " << *r << endl;
// 	 cout << "   Buffer data size = " << *ils << endl;
// 	 cout << "   Number of triplets per comp = " << *ils << endl;
// 	 for (comp=0;comp<NB_COMPS;comp++)
// 	 {
// 	   cout << "      Comp = " << comp << " n = " << (*ibs)[comp].size() 
// 	   	<< endl;
// 	   nijk = (*ibs)[comp].size();
// 	   for (m=0;m<nijk;++m)
// 	     cout << "      " << (*ibs)[comp][m].i << " " << (*ibs)[comp][m].j 
// 	     	<< " " << (*ibs)[comp][m].k << endl;
// 	 }	 
//        }
//        cout << endl; 
//        cout << "   Halozone data" << endl;
//        for ( r=synchronization_MPI_rank_neighbors.begin(),
//   	ihr=halozone_received.begin(),
// 	ils=halozone_received_data_size.begin();
//   	r!=synchronization_MPI_rank_neighbors.end();
// 	r++,ihr++,ils++ )
//        {
//          cout << "   From rank " << *r << endl;
// 	 cout << "   Received data size = " << *ils << endl;
// 	 cout << "   Number of triplets per comp = " << *ils << endl;
// 	 for (comp=0;comp<NB_COMPS;comp++)
// 	 {
// 	   cout << "      Comp = " << comp << " n = " << (*ihr)[comp].size() 
// 	   	<< endl;
// 	   nijk = (*ihr)[comp].size();
// 	   for (m=0;m<nijk;++m)
// 	     cout << "      " << (*ihr)[comp][m].i << " " << (*ihr)[comp][m].j 
// 	     	<< " " << (*ihr)[comp][m].k << endl;
// 	 }
//        }
//        cout << endl;                     
//      }
//      macCOMM->barrier();
//    }
//    if ( RANK == 0 ) cout << endl;       
} 




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_periodic_synchronization_features( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_periodic_synchronization_features" ) ;
   
   size_t locNum = 0, comp, i, j, m ;
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   size_t rank = macCOMM->rank();
   list<FV_TRIPLET> ijk, signed_periodic_shift ;
   FV_TRIPLET mt ;
   mt.i = 0 ;
   mt.j = 0 ;   
   mt.k = 0 ; 
   list<size_t> comps ;
   
   // Fill buffer_features with the (i,j,k,comp) of dof handled by this proc
   // i.e. dof in the buffer zone
   for (comp=0;comp<NB_COMPS;comp++)
   {
     locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[locNum](2);
     for (i=0;i<(*local_dof_number)[locNum](0);++i)
       for (size_t j=0;j<(*local_dof_number)[locNum](1);++j) 
         for (size_t k=0;k<kmax;++k)
	 {            
	   if ( ( (*DOFstatus)[locNum](i,j,k) == FV_DOF_PERIODIC_BUFFERZONE
		|| (*DOFstatus)[locNum](i,j,k) 
			== FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE )
		&& (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) != -1 ) 
	   {
             // X periodicity
	     if ( (*PERIODIC)(0) )
	     {
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			+ (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ; 
		 mt.i = (*PERIODIC_SHIFT)[locNum].i;
		 mt.j = 0 ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ;
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			- (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ;
		 mt.i = - (*PERIODIC_SHIFT)[locNum].i;
		 mt.j = 0 ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ;		  
	       }	        	       
	     }
	     
             // Y periodicity
	     if ( (*PERIODIC)(1) )
	     {
               mt.i = i + (*local_min_index_in_global)[locNum](0) ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ;
		 mt.i = 0;
		 mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ; 
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ;
		 mt.i = 0;
		 mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ; 
	       }	        	       
	     }
	     
             // X-Y biperiodicity
	     if ( (*PERIODIC)(0) && (*PERIODIC)(1) )
	     {
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			+ (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ;
		 mt.i = (*PERIODIC_SHIFT)[locNum].i;
		 mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ; 
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			+ (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ;
		 mt.i = (*PERIODIC_SHIFT)[locNum].i;
		 mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ; 
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			- (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ;
		 mt.i = - (*PERIODIC_SHIFT)[locNum].i;
		 mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ; 
	       }
	       
               mt.i = i + (*local_min_index_in_global)[locNum](0) 
			- (*PERIODIC_SHIFT)[locNum].i ;
               mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
               mt.k = DIM == 2 ? 0 : 
	       		k + (*local_min_index_in_global)[locNum](2) ;	       
	       if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	       {
	         ijk.push_back( mt ) ;
		 comps.push_back( comp ) ;
		 mt.i = - (*PERIODIC_SHIFT)[locNum].i;
		 mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		 mt.k = 0 ;
		 signed_periodic_shift.push_back( mt ) ; 
	       }	       
	     }
	     
	     // 3D cases
	     if ( DIM == 3 )
	     {
               // Z periodicity
	       if ( (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
                   mt.i = 0 ;
		   mt.j = 0 ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;		    
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = 0 ;
		   mt.j = 0 ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;
	         }	        	       
	       }
	       
               // X-Z biperiodicity
	       if ( (*PERIODIC)(0) && (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = 0 ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = 0 ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;  
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = - (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = 0 ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = - (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = 0 ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }		 	        	       
	       }
	       
               // Y-Z biperiodicity
	       if ( (*PERIODIC)(1) && (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = 0 ;
		   mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = 0 ;
		   mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = 0 ;
		   mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;  
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = 0 ;
		   mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }		 	        	       
	       }
	       
               // X-Y-Z triperiodicity
	       if ( (*PERIODIC)(0) && (*PERIODIC)(1) && (*PERIODIC)(2) )
	       {
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		+ (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = - (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		+ (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ;
		   mt.i = - (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ; 
	         }
		 
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		+ (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = - (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;
	         }
	       
                 mt.i = i + (*local_min_index_in_global)[locNum](0) 
	       		- (*PERIODIC_SHIFT)[locNum].i ;
                 mt.j = j + (*local_min_index_in_global)[locNum](1) 
	       		- (*PERIODIC_SHIFT)[locNum].j ;   
                 mt.k = k + (*local_min_index_in_global)[locNum](2) 
	       		- (*PERIODIC_SHIFT)[locNum].k ;	       
	         if ( is_global_triplet( mt.i, mt.j, mt.k, comp ) )	       
	         {
	           ijk.push_back( mt ) ;
		   comps.push_back( comp ) ; 
		   mt.i = - (*PERIODIC_SHIFT)[locNum].i ;
		   mt.j = - (*PERIODIC_SHIFT)[locNum].j ;
		   mt.k = - (*PERIODIC_SHIFT)[locNum].k ;
		   signed_periodic_shift.push_back( mt ) ;
	         }		 		 	        	       
	       }	       	       	       	     
	     }	     	     
	   }
	 }
   }   

   // Transfer the lists into a single vector
   size_t width = 7 ;
   intVector buffer_features( ijk.size()*width, 0 );
   size_t positionIndex = 0 ;
   list<FV_TRIPLET>::iterator imac, ishift;
   list<size_t>::iterator icomp;
   
   for (imac=ijk.begin(),icomp=comps.begin(),
   	ishift=signed_periodic_shift.begin();
   	imac!=ijk.end();icomp++,imac++,ishift++)
   {
      buffer_features(positionIndex) = imac->i ;
      buffer_features(positionIndex+1) = imac->j ;
      buffer_features(positionIndex+2) = imac->k ;
      buffer_features(positionIndex+3) = *icomp ;
      buffer_features(positionIndex+4) = ishift->i ;
      buffer_features(positionIndex+5) = ishift->j ;
      buffer_features(positionIndex+6) = ishift->k ;      	          
      positionIndex += width ;     
   }      

   // Send the bufferzone (i,j,k,comp) features to the neighboring procs
   // and process the received features. We store the (i,j,k,comp) that belong
   // to this proc in halozone_received
   size_t ig = 0, jg = 0, kg = 0 ;   
   list<size_t> MPI_periodic_neighbors_world =
   	*PRIMARY_GRID->get_MPI_periodic_neighbors_ranks() ;
   MPI_periodic_neighbors_world.push_back( rank ) ;
   MPI_periodic_neighbors_world.sort();
   list<size_t>::iterator r, k, ils;
   list<FV_TRIPLET> trip ;
   vector< list<FV_TRIPLET> > kept_ijk_per_comp( NB_COMPS, trip );
   list<FV_TRIPLET>::iterator imt;
   list<size_t> stl;
   vector< list<size_t> > position_bufferzone_list(
   	MPI_periodic_neighbors_world.size(), stl );

   for (r=MPI_periodic_neighbors_world.begin(),i=0;
   	r!=MPI_periodic_neighbors_world.end();r++,++i)
   {
     if ( *r == rank )
     {
       // Send messages to the neighboring procs
       for (k=MPI_periodic_neighbors_world.begin();
       	k!=MPI_periodic_neighbors_world.end();k++)
         if ( *k != rank )
	   macCOMM->send( *k, buffer_features ) ;
	   
       // Check on the proc itself (in case of a single proc in the periodic
       // direction)
       size_t vecsize = buffer_features.size() ; 
       for (positionIndex=0,m=0;positionIndex<vecsize;++m )
       {
         ig = buffer_features(positionIndex);
         jg = buffer_features(positionIndex+1);
         kg = buffer_features(positionIndex+2);
	 comp = buffer_features(positionIndex+3);
	 
         locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;

	 if ( is_global_triplet_local_DOF( ig, jg, kg, locNum ) ) 
	 {
	   mt.i = ig - (*local_min_index_in_global)[locNum](0);
	   mt.j = jg - (*local_min_index_in_global)[locNum](1);
	   mt.k = DIM == 2 ? 0 : kg - 
	     	(*local_min_index_in_global)[locNum](2) ;
	   kept_ijk_per_comp[comp].push_back( mt );
	   position_bufferzone_list[i].push_back( m );
	 }
	 	 	 
         positionIndex += width ;
       }        
     }
     else
     {
       intVector received_buffer_features_from_r(1) ;
       macCOMM->receive( *r, received_buffer_features_from_r ) ;
       size_t vecsize = received_buffer_features_from_r.size() ; 

       for (positionIndex=0,m=0;positionIndex<vecsize;++m )
       {
         ig = received_buffer_features_from_r(positionIndex);
         jg = received_buffer_features_from_r(positionIndex+1);
         kg = received_buffer_features_from_r(positionIndex+2);
	 comp = received_buffer_features_from_r(positionIndex+3);
	 
         locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;

	 if ( is_global_triplet_local_DOF( ig, jg, kg, locNum ) ) 
	 {
	   mt.i = ig - (*local_min_index_in_global)[locNum](0);
	   mt.j = jg - (*local_min_index_in_global)[locNum](1);
	   mt.k = DIM == 2 ? 0 : kg - 
	     	(*local_min_index_in_global)[locNum](2) ;
	   kept_ijk_per_comp[comp].push_back( mt );
	   position_bufferzone_list[i].push_back( m );
	 }
	 	 	 
         positionIndex += width ;
       }
     }
     
     if ( !position_bufferzone_list[i].empty() )
     {
       // Transfer the vector of list into a vector of vector
       vector< FV_TRIPLET > vvv;
       vector< vector< FV_TRIPLET > > vec_kept_ijk_per_comp( NB_COMPS, vvv );
       for (comp=0;comp<NB_COMPS;comp++)
       {
         size_t vsize = kept_ijk_per_comp[comp].size();
	 vec_kept_ijk_per_comp[comp].reserve( vsize );
	 for (m=0;m<vsize;++m) vec_kept_ijk_per_comp[comp].push_back( mt ); 
       }    
       for (comp=0;comp<NB_COMPS;comp++)
       {
         m = 0;
         for (imt=kept_ijk_per_comp[comp].begin();
	 	imt!=kept_ijk_per_comp[comp].end();imt++,++m)
	   vec_kept_ijk_per_comp[comp][m] = *imt;
       }	
	            
       // Important remark: 
       // In the case of 2 procs only in the periodic direction,
       // these 2 procs already exchange data through their interior boundary
       // Now they also exchange data through their periodic boundary
       // Consequently, a MPI rank of a neighboring proc can appear multiple
       // times (once for interior boundary and once for periodic boundary)
       // We therefore send-receive a few more messages as 2 neighboring procs
       // may exchange 2 (or more ) messages instead of 1, but this simplifies 
       // the programming
       
       // Store the vector of vector into halozone_received
       halozone_received.push_back( vec_kept_ijk_per_comp );
       
       // Store the MPI rank of the neighbor
       synchronization_MPI_rank_neighbors.push_back( *r );
       
       // Clear the content of the lists per comp
       for (comp=0;comp<NB_COMPS;comp++)
         kept_ijk_per_comp[comp].clear();       
     }
   }    

   // Send back the halozone_received (i,j,k,comp) features to neighbors
   // as a list of positions in the initial buffer sent by the neighbor
   // such that they can create the corresponding bufferzone_sent (i,j,k,comp)
   // features 
   size_t_vector position_bufferzone( 0, 0 );
   size_t is, js, ks;
   for (r=MPI_periodic_neighbors_world.begin(),i=0;
   	r!=MPI_periodic_neighbors_world.end();r++,++i)
   {
     if ( *r == rank )
     {
       for (k=MPI_periodic_neighbors_world.begin(),j=0;
       	k!=MPI_periodic_neighbors_world.end();k++,++j)
         if ( *k != rank )
	 {
           // Number of positions in the bufferzone list
	   size_t npos = position_bufferzone_list[j].size();
	   
	   // Transfer list into vector
	   position_bufferzone.resize( npos ) ;
	   for (ils=position_bufferzone_list[j].begin(),m=0;
	   	ils!=position_bufferzone_list[j].end();ils++,++m)
	     position_bufferzone(m) = *ils;
	    
	   // Send message
	   macCOMM->send( *k, position_bufferzone ) ; 	   
	 }
	 
       // Check on the proc itself (in case of a single proc in the periodic
       // direction) 
       if ( !position_bufferzone_list[i].empty() )
       { 
         vector<size_t> nijk_per_comp( NB_COMPS, 0 );

	 // Get the number of ijk per comp
         for (ils=position_bufferzone_list[i].begin();
	 	ils!=position_bufferzone_list[i].end();ils++)
         {
           comp = buffer_features((*ils)*width+3);
	   nijk_per_comp[comp]++;
         } 
 
         // Allocate bufferzone_sent for this neighboring proc
         vector< FV_TRIPLET > vvv;
         vector< vector< FV_TRIPLET > > vec_ijk_per_comp( NB_COMPS, vvv );
         for (comp=0;comp<NB_COMPS;comp++)
         {
	   vec_ijk_per_comp[comp].reserve( nijk_per_comp[comp] );
	   for (m=0;m<nijk_per_comp[comp];++m) 
	     vec_ijk_per_comp[comp].push_back( mt ); 
         }

         m = 0;
         size_t previous_comp = 0;
         for (ils=position_bufferzone_list[i].begin();
	 	ils!=position_bufferzone_list[i].end();ils++)
         {
	   ig = buffer_features((*ils)*width);
           jg = buffer_features((*ils)*width+1);
           kg = buffer_features((*ils)*width+2);
           comp = buffer_features((*ils)*width+3);
	   is = buffer_features((*ils)*width+4);
           js = buffer_features((*ils)*width+5);
           ks = buffer_features((*ils)*width+6);		
			 
           locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
	 
	   if ( comp != previous_comp ) m = 0; // data are ordered per comp

           vec_ijk_per_comp[comp][m].i = ig - is  
		- (*local_min_index_in_global)[locNum](0);
           vec_ijk_per_comp[comp][m].j = jg - js 
		- (*local_min_index_in_global)[locNum](1);		
           vec_ijk_per_comp[comp][m].k = DIM == 2 ? 0 : kg - ks 
	     	- (*local_min_index_in_global)[locNum](2); 
	   ++m;
	   previous_comp = comp;	 
         }

         // Store the vector of vector into bufferzone_sent
         bufferzone_sent.push_back( vec_ijk_per_comp );  	 	       
       }       
     }
     else
     {
       size_t_vector received_buffer_features_from_r(1) ;
       macCOMM->receive( *r, received_buffer_features_from_r ) ;
       size_t vecsize = received_buffer_features_from_r.size() ; 
       vector<size_t> nijk_per_comp( NB_COMPS, 0 );

       // Get the number of ijk per comp
       for (positionIndex=0;positionIndex<vecsize;++positionIndex)
       {
         comp = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width+3);
	 nijk_per_comp[comp]++;
       }       

       // Allocate bufferzone_sent for this neighboring proc
       vector< FV_TRIPLET > vvv;
       vector< vector< FV_TRIPLET > > vec_ijk_per_comp( NB_COMPS, vvv );
       for (comp=0;comp<NB_COMPS;comp++)
       {
	 vec_ijk_per_comp[comp].reserve( nijk_per_comp[comp] );
	 for (m=0;m<nijk_per_comp[comp];++m) 
	   vec_ijk_per_comp[comp].push_back( mt ); 
       }
              
       m = 0;
       size_t previous_comp = 0;
       for (positionIndex=0;positionIndex<vecsize;++positionIndex)
       {
	 ig = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width);
         jg = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width+1);
         kg = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width+2);
         comp = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width+3);
	 is = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width+4);
         js = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width+5);
         ks = buffer_features(received_buffer_features_from_r(positionIndex)
	 	*width+6);		
			 
         locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
	 
	 if ( comp != previous_comp ) m = 0; // data are ordered per comp

         vec_ijk_per_comp[comp][m].i = ig - is  
		- (*local_min_index_in_global)[locNum](0);
         vec_ijk_per_comp[comp][m].j = jg - js 
		- (*local_min_index_in_global)[locNum](1);		
         vec_ijk_per_comp[comp][m].k = DIM == 2 ? 0 : kg - ks 
	     	- (*local_min_index_in_global)[locNum](2); 
	 ++m;
	 previous_comp = comp;	 
       }
       
       // Store the vector of vector into bufferzone_sent
       bufferzone_sent.push_back( vec_ijk_per_comp );       
     }
   }     
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: synchronize( size_t level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: synchronize_field" ) ;
   MAC_ASSERT( level < STO_DEPTH );

  if ( !synchronization_ready ) set_synchronization_features();

  MAC_Communicator const* macCOMM = MAC_Exec::communicator();
  list<size_t>::iterator r, k;
  list< double* >::iterator ild;
  list< size_t >::iterator ils;
  list< vector< vector< FV_TRIPLET > > >::iterator ihr, ibs;
  size_t comp, ntriplets, m, positionIndex;
  void* sreq = NULL; 
  vector<void*> idreq( synchronization_MPI_rank_neighbors.size(), sreq );
  size_t ireq; 

  // Important remark:
  // When a proc sends 2 separate messages to another proc with the same tag and
  // matching receives, the MPI standard guarantees that messages are received 
  // in the order that they were sent (Messages are non-overtaking) provided 
  // (i) multi-threading is no performed 
  // (ii) MPI_ANY_SOURCE is 
  // Here, the send-receive method implemented in MAC_Communicator uses a tag 
  // per data type, hence below double vectors are all sent with the same tag.
  // We use neither multi-threading nor MPI_ANY_SOURCE anywhere, therefore the 
  // programming model below is robust while same tag + (i) + (ii) is fulfilled.
   
  // Send bufferzone data to neighboring procs
  for ( r=synchronization_MPI_rank_neighbors.begin(),
  	ibs=bufferzone_sent.begin(),ild=bufferzone_sent_data.begin(),
	ils=bufferzone_sent_data_size.begin(),ireq=0;
  	r!=synchronization_MPI_rank_neighbors.end();
	r++,ibs++,ild++,ils++,ireq++ )
  {
    // Fill the bufferzone_sent_data array 
    positionIndex = 0;
    for (comp=0;comp<NB_COMPS;comp++)
    {
      ntriplets = (*ibs)[comp].size();
      for (m=0;m<ntriplets;++m)
      {
        (*ild)[positionIndex] = (*VALUES)[level][comp]( 
		(*ibs)[comp][m].i, (*ibs)[comp][m].j, (*ibs)[comp][m].k );
	++positionIndex;   
      }  
    }
    
    // Send message
    idreq[ireq] = macCOMM->Isend( *r, *ild, *ils );      	
  }

  // Receive data in halozone from neighboring procs 
  for ( r=synchronization_MPI_rank_neighbors.begin(),
  	ihr=halozone_received.begin(),ild=halozone_received_data.begin(),
	ils=halozone_received_data_size.begin(),ireq=0;
  	r!=synchronization_MPI_rank_neighbors.end();
	r++,ihr++,ild++,ils++ )
  {
    // Receive the halozone_received_data array
    macCOMM->receive( *r, *ild, *ils );
    
    // Update the halozone values    
    positionIndex = 0;
    for (comp=0;comp<NB_COMPS;comp++)
    {
      ntriplets = (*ihr)[comp].size();
      for (m=0;m<ntriplets;++m)
      {
        (*VALUES)[level][comp]( (*ihr)[comp][m].i, (*ihr)[comp][m].j, 
		(*ihr)[comp][m].k ) = (*ild)[positionIndex];
	++positionIndex;   
      }  
    }     	
  }  
  
  // Check that all non-blocking messages are complete
  for (ireq=0;ireq<idreq.size();++ireq)
    macCOMM->wait( idreq[ireq] ); 
    
    
//    // Debug
//    size_t nb_ranks = macCOMM->nb_ranks() ;
//    size_t RANK = macCOMM->rank() ;
//    size_t nijk = 0;
//    for (size_t i=0;i<nb_ranks;++i)
//    {
//      if ( i == RANK )
//      {       
//        cout << "Rank " << RANK << endl ;
//        cout << "   Buffer data" << endl;
//        for ( r=synchronization_MPI_rank_neighbors.begin(),
//   	ibs=bufferzone_sent.begin(),
// 	ils=bufferzone_sent_data_size.begin();
//   	r!=synchronization_MPI_rank_neighbors.end();
// 	r++,ibs++,ils++ )
//        {
//          cout << "   To rank " << *r << endl;
// 	 cout << "   Buffer data size = " << *ils << endl;
// 	 cout << "   Number of triplets per comp = " << *ils << endl;
// 	 for (comp=0;comp<NB_COMPS;comp++)
// 	 {
// 	   cout << "      Comp = " << comp << " n = " << (*ibs)[comp].size() 
// 	   	<< endl;
// 	   nijk = (*ibs)[comp].size();
// 	   for (m=0;m<nijk;++m)
// 	     cout << "      " << (*ibs)[comp][m].i << " " << (*ibs)[comp][m].j 
// 	     	<< " " << (*ibs)[comp][m].k << " " << 
// 		MAC::doubleToString( std::ios::scientific, 8, 
// 		(*VALUES)[level][comp]( 
// 		(*ibs)[comp][m].i, (*ibs)[comp][m].j, (*ibs)[comp][m].k ) )
// 		<< endl;
// 	 }	 
//        }
//        cout << endl; 
//        cout << "   Halozone data" << endl;
//        for ( r=synchronization_MPI_rank_neighbors.begin(),
//   	ihr=halozone_received.begin(),
// 	ils=halozone_received_data_size.begin();
//   	r!=synchronization_MPI_rank_neighbors.end();
// 	r++,ihr++,ils++ )
//        {
//          cout << "   From rank " << *r << endl;
// 	 cout << "   Received data size = " << *ils << endl;
// 	 cout << "   Number of triplets per comp = " << *ils << endl;
// 	 for (comp=0;comp<NB_COMPS;comp++)
// 	 {
// 	   cout << "      Comp = " << comp << " n = " << (*ihr)[comp].size() 
// 	   	<< endl;
// 	   nijk = (*ihr)[comp].size();
// 	   for (m=0;m<nijk;++m)
// 	     cout << "      " << (*ihr)[comp][m].i << " " << (*ihr)[comp][m].j 
// 	     	<< " " << (*ihr)[comp][m].k << " " << 
// 		MAC::doubleToString( std::ios::scientific, 8, 
// 		(*VALUES)[level][comp]( (*ihr)[comp][m].i, (*ihr)[comp][m].j, 
// 		(*ihr)[comp][m].k ) ) << endl;
// 	 }
//        }
//        cout << endl;                     
//      }
//      macCOMM->barrier();
//    }
//    if ( RANK == 0 ) cout << endl;        
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: print_BC_FieldValues( std::ostream& os, 
	size_t indent_width, bool b_values ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: print_BC_FieldValues" ) ;

   std::string space( indent_width, ' ' ) ;      
   size_t locNum = 0 ;
         
   os << space << "Boundary conditions" << endl;
   os << space << "-------------------" << endl;   
   list< FV_BoundaryCondition* >::const_iterator ibc;
   size_t color=0;
   for (ibc=SET_OF_BCS.begin();ibc!=SET_OF_BCS.end();ibc++,++color)
   {
     (*ibc)->print( os, 0, DIM );
     os << endl;
   }
   
   if ( b_values )
   {
     os << endl;     
     os << space << "Field values" << endl;     
     os << space << "------------" << endl; 
     os << space << "Nb local dofs = " << NB_LOCAL_DOF << endl;
     os << space << "Nb local unknowns = " << NB_LOCAL_UNKNOWNS << endl;     
     os << space << "Nb local unknowns handled by processor = " 
     	<< NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC << endl;
     os << space << "Nb local unknowns handled by processor in buffer zone = " 
     	<< NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE << endl;	      
     os << space << "Nb global unknowns = " << NB_GLOBAL_UNKNOWNS << endl;
     for (size_t comp=0;comp<NB_COMPS;++comp)
     {
       locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;       
       os << space << "Component " << comp << endl;
       os << space << "----------- " << endl;
       os << space << "Min index of unknown handled by processor:";
       for (size_t i=0;i<DIM;++i) 
         os << " " << (*FV_Mesh::directionName)(i) << "=" << 
	 	(*min_index_unknown_handled_by_proc)[locNum](i);
       os << endl;
       os << space << "Max index of unknown handled by processor:";
       for (size_t i=0;i<DIM;++i) 
         os << " " << (*FV_Mesh::directionName)(i) << "=" << 
	 	(*max_index_unknown_handled_by_proc)[locNum](i);
       os << endl;       
       if ( DIM == 2 )
         os << space << "# i_loc j_loc i_glob j_glob color status unkNumLoc "
	  	<< "unkNumGlob Value" << endl;
       else
         os << space << "# i_loc j_loc k_loc i_glob j_glob k_glob color status "
	 	<< "unkNumLoc unkNumGlob Value" << endl;  
       
       size_t nelem0 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(0) ;
       size_t nelem1 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(1) ;
       size_t nelem2 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(2) ;        
             
       for (size_t lev=0;lev<STO_DEPTH;++lev)
       {
         os << space << "# Level " << lev << endl;	          
         if ( DIM == 2 )
         { 
           for (size_t i=0;i<nelem0;++i)
             for (size_t j=0;j<nelem1;++j) 
               os << space << i << " " << j << " " 
	       	<< i+(*local_min_index_in_global)[locNum](0) << " " 
	       	<< j+(*local_min_index_in_global)[locNum](1) << " " 
		<< FV_DomainBuilder::get_color_name( 
			(*DOFcolors)[locNum](i,j,0) ) << " " 
		<< FV_DomainBuilder::get_status_name( 
			(*DOFstatus)[locNum](i,j,0) ) << " "		
		<< (*UNK_LOCAL_NUMBERING)[comp](i,j,0) << " "
		<< (*UNK_GLOBAL_NUMBERING)[comp](i,j,0) << " "
		<< (*VALUES)[lev][comp](i,j,0) << endl;	         
         }
         else
         {
           for (size_t i=0;i<nelem0;++i)
             for (size_t j=0;j<nelem1;++j)
	       for (size_t k=0;k<nelem2;++k) 
                 os << space << i << " " << j << " " << k << " " 
	   	<< i+(*local_min_index_in_global)[locNum](0) << " "
	   	<< j+(*local_min_index_in_global)[locNum](1) << " " 
		<< k+(*local_min_index_in_global)[locNum](2) << " " 
		<< FV_DomainBuilder::get_color_name( 
			(*DOFcolors)[locNum](i,j,k) ) << " " 
		<< FV_DomainBuilder::get_status_name( 
			(*DOFstatus)[locNum](i,j,k) ) << " "		
		<< (*UNK_LOCAL_NUMBERING)[comp](i,j,k) << " "
		<< (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) << " "
		<< (*VALUES)[lev][comp](i,j,k) << endl; 
	 }
       }     
     }       
   }     
}




//----------------------------------------------------------------------
FV_BoundaryCondition const*
FV_DiscreteField:: get_BC( size_t color ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: get_BC" ) ;
   MAC_CHECK_PRE( color > 1 );     
   MAC_CHECK_PRE( IMPLIES( DIM == 2, color < 10 ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 3, color < 28 ) ) ; 
   
   list< FV_BoundaryCondition* >::const_iterator il;
   bool found = false ;
   FV_BoundaryCondition const* result = NULL ;
   
   for (il=SET_OF_BCS.begin(); il!=SET_OF_BCS.end() && !found; il++)
     if ( (*il)->get_color_ID() == color )
     {
       result = *il ;
       found = true ;
     }

   return ( result ) ;
}   




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_DOF_value( size_t i, size_t j, size_t k,
      	size_t component, size_t level, double val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_DOF_value" ) ;
   MAC_CHECK_PRE( i < (*VALUES)[level][component].index_bound(0) );
   MAC_CHECK_PRE( j < (*VALUES)[level][component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*VALUES)[level][component].index_bound(2) );   
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( level < STO_DEPTH );
   
   // Non-unknown DOF cannot be modified, except if forced ( 
   // SET_BC_VALUES_ALLOWED = true )
   MAC_ASSERT( IMPLIES( (*UNK_LOCAL_NUMBERING)[component](i,j,k) == -1, 
   	SET_BC_VALUES_ALLOWED == true ) ) ;           
     
   (*VALUES)[level][component]( i, j, k ) = val ;  
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_DOFs_value( size_t component, size_t level, 
	double value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_DOFs_value" ) ;
   
   size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[component](2) ;
     
   for (size_t i=0;i<(*local_dof_number)[component](0);++i)
     for (size_t j=0;j<(*local_dof_number)[component](1);++j)
       for (size_t k=0;k<kmax;++k)
	 (*VALUES)[level][component]( i, j, k ) = value ;         

}




//----------------------------------------------------------------------
void
FV_DiscreteField:: add_value_to_DOF( size_t i, size_t j, size_t k,
      	size_t component, size_t level, double value )
//----------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: add_value_to_DOF" ) ;
   
  (*VALUES)[level][component]( i, j, k ) += value ;      

}




//----------------------------------------------------------------------
double
FV_DiscreteField:: DOF_value( size_t i, size_t j, size_t k,
      	size_t component, size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: DOF_value" ) ;
   MAC_CHECK_PRE( i < (*VALUES)[level][component].index_bound(0) );
   MAC_CHECK_PRE( j < (*VALUES)[level][component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*VALUES)[level][component].index_bound(2) );  
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( level < STO_DEPTH );         

   double result = (*VALUES)[level][component]( i, j, k ) ; 
   
   return ( result ) ; 
}




//----------------------------------------------------------------------
int
FV_DiscreteField:: DOF_color( size_t i, size_t j, size_t k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: DOF_color" ) ;
   MAC_CHECK_PRE( i < (*DOFcolors)[component].index_bound(0) );
   MAC_CHECK_PRE( j < (*DOFcolors)[component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*DOFcolors)[component].index_bound(2) );  

   size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
   int result = (*DOFcolors)[locNum]( i, j, k ) ; 
   
   return ( result ) ; 
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: build_system_numbering( size_t_vector* idx_locs, 
      	size_t_vector* idx_globs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: build_system_numbering" ) ;
   MAC_CHECK_PRE( idx_locs != 0 );
   MAC_CHECK_PRE( idx_globs != 0 );   
   MAC_CHECK_PRE( idx_locs->size() == NB_LOCAL_UNKNOWNS );
   MAC_CHECK_PRE( idx_globs->size() == NB_LOCAL_UNKNOWNS );

   size_t idx = 0;   
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t nelem0 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(0) ;
     size_t nelem1 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(1) ;
     size_t nelem2 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(2) ; 
     for (size_t i=0;i<nelem0;++i)
       for (size_t j=0;j<nelem1;++j) 
         for (size_t k=0;k<nelem2;++k)
	 {
	   int locidx = (*UNK_LOCAL_NUMBERING)[comp](i,j,k) ;
	   if ( locidx != -1 )
	   {
	     (*idx_locs)(idx) = locidx ;
	     (*idx_globs)(idx) = (*UNK_GLOBAL_NUMBERING)[comp](i,j,k) ;
	     ++idx ;	     
	   }
	 }
   } 
}   
   
   
   

//----------------------------------------------------------------------
void
FV_DiscreteField:: initialize_DOFs( 
	MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: initialize_DOFs" ) ;

   if ( exp->has_module( "DOFs_values" ) )
   {
     // Read formula
     MAC_ContextSimple* c = MAC_ContextSimple::create( this ) ;
     COORDS = MAC_DoubleVector::create( c, doubleVector(0) ) ;
     c->extend( MAC_Variable::object( "DV_X" ), COORDS ) ;
     CTX = c ; 

     MAC_ModuleExplorer* sse = exp->create_subexplorer( 0, "DOFs_values" ) ;   
     
     if ( !sse->has_entry( "value" ) )
     {
       MAC_Error::object()->raise_missing_keyword( sse, "value" ) ; 
     }
     
     INITIALIZER = sse->abstract_data( CTX, "value", CTX ) ;
     if( INITIALIZER->data_type() != MAC_Data::DoubleVector )
     {
       MAC_Error::object()->raise_bad_data_type( sse, "value", 
		MAC_Data::DoubleVector ) ;
     }
     doubleVector dof_coordinates( DIM, 0. );
     COORDS->set( dof_coordinates ) ;
     doubleVector vv = INITIALIZER->to_double_vector( CTX );
     if ( vv.size() != NB_COMPS )
     {
        MAC_Error::object()->raise_data_error( sse, "value", 
		"number of values does not match components number" ) ;
     }
     
     sse->destroy() ; 
     
     // Initialize DOFs
     bool b_setDOF = false ;
     for (size_t comp=0;comp<NB_COMPS;comp++)
     {     
       size_t locNum = ALL_COMPS_SAME_LOCATION == false ? comp : 0 ;
       size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[locNum](2) ;
       for (size_t i=0;i<(*local_dof_number)[locNum](0);++i)
       {
         dof_coordinates(0) = get_DOF_coordinate( i, locNum, 0 ); 
         for (size_t j=0;j<(*local_dof_number)[locNum](1);++j)
         {
           dof_coordinates(1) = get_DOF_coordinate( j, locNum, 1 );  
           for (size_t k=0;k<kmax;++k) 
	   {
	     if ( DIM == 3 ) 
	       dof_coordinates(2) = get_DOF_coordinate( k, locNum, 2 );
	     COORDS->set( dof_coordinates ) ;
	     doubleVector const& val = INITIALIZER->to_double_vector() ;

	     if ( (*DOFcolors)[locNum](i,j,k) < 2 ) b_setDOF = true ;
	     else if ( !(*V_BCS)[(*DOFcolors)[locNum](i,j,k)]
	     	->is_dirichlet(comp) )
	       b_setDOF = true ;
	     else b_setDOF = false ;
	     
	     if ( b_setDOF )
	       for (size_t lev=0;lev<STO_DEPTH;++lev)
	         (*VALUES)[lev][comp]( i, j, k ) = val( comp ) ;
	   }
	 }
       }
     }
   }          
}




//----------------------------------------------------------------------
bool
FV_DiscreteField:: DOF_is_unknown( size_t i, size_t j, size_t k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: DOF_is_unknown" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );
   MAC_CHECK_PRE( i < (*UNK_LOCAL_NUMBERING)[component].index_bound(0) );
   MAC_CHECK_PRE( j < (*UNK_LOCAL_NUMBERING)[component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*UNK_LOCAL_NUMBERING)[component].index_bound(2) );

   bool result = (*UNK_LOCAL_NUMBERING)[component]( i, j, k ) != -1 ; 
   
   return ( result ) ; 
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: DOF_local_number( size_t i, size_t j, size_t k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: DOF_local_number" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );
   MAC_CHECK_PRE( i < (*UNK_LOCAL_NUMBERING)[component].index_bound(0) );
   MAC_CHECK_PRE( j < (*UNK_LOCAL_NUMBERING)[component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*UNK_LOCAL_NUMBERING)[component].index_bound(2) );

   size_t result = (*UNK_LOCAL_NUMBERING)[component]( i, j, k ) ; 
   
   return ( result ) ; 
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: DOF_global_number( size_t i, size_t j, size_t k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: DOF_global_number" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( i < (*UNK_GLOBAL_NUMBERING)[component].index_bound(0) );
   MAC_CHECK_PRE( j < (*UNK_GLOBAL_NUMBERING)[component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*UNK_GLOBAL_NUMBERING)[component].index_bound(2) );

   size_t result = (*UNK_GLOBAL_NUMBERING)[component]( i, j, k ) ; 
   
   return ( result ) ; 
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_min_max_indices_unknown_handled_by_proc( 
	size_t i, size_t j, size_t k, size_t component )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField:: set_min_max_indices_unknown_handled_by_proc" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );     
   MAC_CHECK_PRE( i < (*UNK_GLOBAL_NUMBERING)[component].index_bound(0) );
   MAC_CHECK_PRE( j < (*UNK_GLOBAL_NUMBERING)[component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*UNK_GLOBAL_NUMBERING)[component].index_bound(2) );   

   size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;

   (*min_index_unknown_handled_by_proc)[locNum](0) = 
   	i < (*min_index_unknown_handled_by_proc)[locNum](0) ? 
	i : (*min_index_unknown_handled_by_proc)[locNum](0) ;
   (*min_index_unknown_handled_by_proc)[locNum](1) = 
   	j < (*min_index_unknown_handled_by_proc)[locNum](1) ? 
	j : (*min_index_unknown_handled_by_proc)[locNum](1) ;
   (*max_index_unknown_handled_by_proc)[locNum](0) = 
   	i > (*max_index_unknown_handled_by_proc)[locNum](0) ? 
	i : (*max_index_unknown_handled_by_proc)[locNum](0) ;
   (*max_index_unknown_handled_by_proc)[locNum](1) = 
   	j > (*max_index_unknown_handled_by_proc)[locNum](1) ? 
	j : (*max_index_unknown_handled_by_proc)[locNum](1) ;		
   if ( DIM == 3 )
   {
     (*min_index_unknown_handled_by_proc)[locNum](2) = 
   	k < (*min_index_unknown_handled_by_proc)[locNum](2) ? 
	k : (*min_index_unknown_handled_by_proc)[locNum](2) ;
     (*max_index_unknown_handled_by_proc)[locNum](2) = 
   	k > (*max_index_unknown_handled_by_proc)[locNum](2) ? 
	k : (*max_index_unknown_handled_by_proc)[locNum](2) ;	
   }
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_min_max_indices_unknown_on_proc( 
	size_t i, size_t j, size_t k, size_t component )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField:: set_min_max_indices_unknown_on_proc" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );     
   MAC_CHECK_PRE( i < (*UNK_GLOBAL_NUMBERING)[component].index_bound(0) );
   MAC_CHECK_PRE( j < (*UNK_GLOBAL_NUMBERING)[component].index_bound(1) );   
   MAC_CHECK_PRE( k < (*UNK_GLOBAL_NUMBERING)[component].index_bound(2) );   

   size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;

   (*min_index_unknown_on_proc)[locNum](0) = 
   	i < (*min_index_unknown_on_proc)[locNum](0) ? 
	i : (*min_index_unknown_on_proc)[locNum](0) ;
   (*min_index_unknown_on_proc)[locNum](1) = 
   	j < (*min_index_unknown_on_proc)[locNum](1) ? 
	j : (*min_index_unknown_on_proc)[locNum](1) ;
   (*max_index_unknown_on_proc)[locNum](0) = 
   	i > (*max_index_unknown_on_proc)[locNum](0) ? 
	i : (*max_index_unknown_on_proc)[locNum](0) ;
   (*max_index_unknown_on_proc)[locNum](1) = 
   	j > (*max_index_unknown_on_proc)[locNum](1) ? 
	j : (*max_index_unknown_on_proc)[locNum](1) ;		
   if ( DIM == 3 )
   {
     (*min_index_unknown_on_proc)[locNum](2) = 
   	k < (*min_index_unknown_on_proc)[locNum](2) ? 
	k : (*min_index_unknown_on_proc)[locNum](2) ;
     (*max_index_unknown_on_proc)[locNum](2) = 
   	k > (*max_index_unknown_on_proc)[locNum](2) ? 
	k : (*max_index_unknown_on_proc)[locNum](2) ;	
   }
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: out_endOfBuilding( std::ostream& os, 
	size_t indent_width, size_t rank ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: out_endOfBuilding" ) ;

   std::string space( indent_width, ' ' ) ;
   std::string plus( 3, ' ' ) ;
   size_t nb_ranks = MAC_Exec::communicator()->nb_ranks() ; 
        
   if ( rank == 0 )
   {
     os << endl ;
     os << space << FNAME << endl ;
     os << space << plus << "Field ID = " << ID << endl ;
     os << space << plus << "Discretization = " << FDISCRETIZATION 
   	<< endl ;
     os << space << plus << "Number of components = " << NB_COMPS 
     	<< endl ;
     os << space << plus << "Storage depth = " << STO_DEPTH 
     	<< endl ;
     os << space << plus << "Global number of unknowns = " << 
     	NB_GLOBAL_UNKNOWNS << endl << endl;
   }
				
   for (size_t i=0;i<nb_ranks;++i)
   {
     if ( i == rank ) 
        os << space << plus << "Local number of"
		<< " unknowns handled by rank " << i << " = " 
		<< NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC << endl;  
     MAC_Exec::communicator()->barrier();
   }
   
   if ( rank == 0 ) FV::out() << endl ;       
}




//----------------------------------------------------------------------
bool
FV_DiscreteField:: DOF_in_domain( int i, int j, int k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: DOF_in_domain" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 

   bool result = true ;
   size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;

   if ( i + int((*local_min_index_in_global)[locNum](0)) < 0
   	|| j + int((*local_min_index_in_global)[locNum](1)) < 0 )
     result = false ;
   else if ( i + int((*local_min_index_in_global)[locNum](0)) >
   	int((*global_max_index)[locNum](0))
   	|| j + int((*local_min_index_in_global)[locNum](1)) >
	int((*global_max_index)[locNum](1)) ) result = false ;
	
   if ( result && DIM == 3 )
   {
     if ( k + int((*local_min_index_in_global)[locNum](2)) < 0 )
       result = false ;
     else if ( k + int((*local_min_index_in_global)[locNum](2)) >
   	int((*global_max_index)[locNum](2)) ) result = false ;
   }

   return ( result ) ; 
}




//----------------------------------------------------------------------
bool
FV_DiscreteField:: DOF_on_proc( int i, int j, int k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: DOF_on_proc" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 

   bool result = true ;

   if ( i < 0 || j < 0 || k < 0 ) result = false ;
   else 
   {
     size_t ii = abs(i), jj = abs(j), kk = abs(k) ; 
     if (  ii >= (*UNK_GLOBAL_NUMBERING)[component].index_bound(0) 
     	|| jj >= (*UNK_GLOBAL_NUMBERING)[component].index_bound(1) ) 
       result = false ; 
     else if ( DIM == 3 )
       if ( kk >= (*UNK_GLOBAL_NUMBERING)[component].index_bound(2) ) 
         result = false ;
   } 
        
   return ( result ) ; 
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: copy_DOFs_value(
			size_t source_level, size_t target_level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: copy_DOFs_value" ) ;
   MAC_CHECK_PRE( source_level < STO_DEPTH ); 
   MAC_CHECK_PRE( target_level < STO_DEPTH );
   
   for (size_t comp=0;comp<NB_COMPS;comp++)
   {
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[comp](2) ;
     for (size_t i=0;i<(*local_dof_number)[comp](0);++i)
       for (size_t j=0;j<(*local_dof_number)[comp](1);++j)
         for (size_t k=0;k<kmax;++k)
	   (*VALUES)[target_level][comp]( i, j, k )
	   	= (*VALUES)[source_level][comp]( i, j, k ) ;
   }          
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: add_to_DOFs_value( size_t component, size_t level,
      		double const& val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: add_to_DOFs_value" ) ;
   MAC_CHECK_PRE( level < STO_DEPTH ); 
   MAC_CHECK_PRE( component < NB_COMPS );
   
   size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[component](2) ;
   for (size_t i=0;i<(*local_dof_number)[component](0);++i)
     for (size_t j=0;j<(*local_dof_number)[component](1);++j)
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
FV_SHIFT_TRIPLET
FV_DiscreteField:: shift_staggeredToCentered( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: shift_staggeredToCentered" ) ;

   FV_SHIFT_TRIPLET result ;
   result.i = 0 ;
   result.j = 0 ;   
   result.k = 0 ;   
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the FV_DiscreteField_Centered type"
   	<<"; method \"shift_staggeredToCentered\" can be called by a field of "
	<< "type FV_DiscreteField_Centered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET
FV_DiscreteField:: shift_vertexToStaggered( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: shift_vertexToStaggered" ) ;

   FV_SHIFT_TRIPLET result ;
   result.i = 0 ;
   result.j = 0 ;   
   result.k = 0 ;   
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the FV_DiscreteField_Vertex type"
   	<<"; method \"shift_vertexToStaggered\" can be called by a field of "
	<< "type FV_DiscreteField_Vertex only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET
FV_DiscreteField:: shift_staggeredToStaggered( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: shift_staggeredToStaggered" ) ;

   FV_SHIFT_TRIPLET result ;
   result.i = 0 ;
   result.j = 0 ;   
   result.k = 0 ;   
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the FV_DiscreteField_Staggered "
   	<<"type; method \"shift_staggeredToStaggered\" can be called by "
	<< "a field of type FV_DiscreteField_Staggered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET
FV_DiscreteField:: shift_vorticityToStaggered( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: shift_vorticityToStaggered" ) ;

   FV_SHIFT_TRIPLET result ;
   result.i = 0 ;
   result.j = 0 ;   
   result.k = 0 ;   
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the FV_DiscreteField_Vorticity "
   	<<"type; method \"shift_vorticityToStaggered\" can be called by "
	<< "a field of type FV_DiscreteField_Vorticity only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET
FV_DiscreteField:: shift_staggeredToTensor( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: shift_staggeredToTensor" ) ;

   FV_SHIFT_TRIPLET result ;
   result.i = 0 ;
   result.j = 0 ;   
   result.k = 0 ;   
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the FV_DiscreteField_tensor "
   	<<"type; method \"shift_staggeredToTensor\" can be called by "
	<< "a field of type FV_DiscreteField_Tensor only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET
FV_DiscreteField:: shift_tensorToTensor( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: shift_tensorToTensor" ) ;

   FV_SHIFT_TRIPLET result ;
   result.i = 0 ;
   result.j = 0 ;   
   result.k = 0 ;   
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the FV_DiscreteField_tensor "
   	<<"type; method \"shift_tensorToTensor\" can be called by "
	<< "a field of type FV_DiscreteField_Tensor only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//----------------------------------------------------------------------
double 
FV_DiscreteField:: interpolateOneCompOnAnotherComp(
		size_t i, size_t j,
		size_t comp1, size_t comp2, size_t level,
		FV_SHIFT_TRIPLET shift ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: interpolateOneCompOnAnotherComp" ) ;

   ostringstream mesg ;
   mesg << "This method is not implemented for Field " << FNAME << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( 0. ) ;
}




//----------------------------------------------------------------------
double 
FV_DiscreteField:: interpolateOneCompOnAnotherComp(
		size_t i, size_t j, size_t k,
		size_t comp1, size_t comp2, size_t level,
		FV_SHIFT_TRIPLET shift ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: interpolateOneCompOnAnotherComp" ) ;

   ostringstream mesg ;
   mesg << "This method is not implemented for Field " << FNAME << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( 0. ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: set_BC_values_modif_status( bool const& allowed ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_BC_values_modif_status" ) ;

   SET_BC_VALUES_ALLOWED = allowed ;
}       
 



//----------------------------------------------------------------------
void
FV_DiscreteField:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "FV_DiscreteField" ) ;
   
   writer->add_entry( "name", MAC_String::create( 0, FNAME ) ) ;
   writer->add_entry( "nb_DOF", MAC_Int::create( 0, NB_LOCAL_DOF ) ) ;
   writer->add_entry( "nb_components", MAC_Int::create( 0, NB_COMPS ) ) ;
   writer->add_entry( "nb_levels", MAC_Int::create( 0, STO_DEPTH ) ) ;
   
   // Saving values
   doubleArray3D OutputValues( NB_COMPS, STO_DEPTH, NB_LOCAL_DOF, 0. ) ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t nelem0 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(0) ;
     size_t nelem1 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(1) ;
     size_t nelem2 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(2) ;        
             
     for (size_t lev=0;lev<STO_DEPTH;++lev)
     {
       size_t idx = 0 ;
       if ( DIM == 2 )
       {
         for (size_t i=0;i<nelem0;++i)
           for (size_t j=0;j<nelem1;++j) 
             OutputValues( comp, lev, idx++ ) = (*VALUES)[lev][comp](i,j,0) ;
       }
       else
       {
         for (size_t i=0;i<nelem0;++i)
           for (size_t j=0;j<nelem1;++j)
	     for (size_t k=0;k<nelem2;++k) 
               OutputValues( comp, lev, idx++ ) = (*VALUES)[lev][comp](i,j,k) ;
       }
     }
   }
   
   writer->add_entry( "values", MAC_DoubleArray3D::create( 0, OutputValues ) ) ;
   
   writer->finalize_object() ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "FV_DiscreteField" ) ;

   // Retrieving stored datas :
   string const& read_name = reader->data_of_entry( "name" )->to_string() ;
   int const read_nb_DOF = reader->data_of_entry( "nb_DOF" )->to_int() ;
   int const read_nb_comps = reader->data_of_entry( "nb_components" )->to_int();
   int const read_nb_levels = reader->data_of_entry( "nb_levels" )->to_int() ;

   // Does some checks
   MAC_ASSERT( read_name==FNAME ) ;
   MAC_ASSERT( read_nb_DOF==(int) NB_LOCAL_DOF ) ;
   MAC_ASSERT( read_nb_comps==(int) NB_COMPS ) ;
   MAC_ASSERT( read_nb_levels==(int) STO_DEPTH ) ;

   // Retrieving values
   doubleArray3D OutputValues
   		= reader->data_of_entry( "values" )->to_double_array3D() ;
   MAC_ASSERT( OutputValues.index_bound(0) == NB_COMPS ) ;
   MAC_ASSERT( OutputValues.index_bound(1) == STO_DEPTH ) ;
   MAC_ASSERT( OutputValues.index_bound(2) == NB_LOCAL_DOF ) ;    
   
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t nelem0 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(0) ;
     size_t nelem1 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(1) ;
     size_t nelem2 = (*UNK_LOCAL_NUMBERING)[comp].index_bound(2) ;        

     for (size_t lev=0;lev<STO_DEPTH;++lev)
     {
       size_t idx = 0 ;
       if ( DIM == 2 )
       {
         for (size_t i=0;i<nelem0;++i)
           for (size_t j=0;j<nelem1;++j)
             (*VALUES)[lev][comp](i,j,0) = OutputValues( comp, lev, idx++ ) ;
       }
       else
       {
         for (size_t i=0;i<nelem0;++i)
           for (size_t j=0;j<nelem1;++j)
	     for (size_t k=0;k<nelem2;++k) 
               (*VALUES)[lev][comp](i,j,k) = OutputValues( comp, lev, idx++ ) ;
       }
     }
   }   

   reader->end_object_retrieval() ;

   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: read_state_nonrestored( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: read_state_nonrestored" ) ;

   reader->start_object_retrieval( "FV_DiscreteField" ) ;
   reader->end_object_retrieval() ;
   
}




//----------------------------------------------------------------------
bool
FV_DiscreteField:: all_BCs_nonDirichlet( size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: all_BCs_nonDirichlet" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
      
   bool result = true ;
   
   for ( list< FV_BoundaryCondition* >::const_iterator il=SET_OF_BCS.begin();
   	il!=SET_OF_BCS.end() && result == true ;il++) 
     result = !(*il)->is_dirichlet( component ) ;  	         
   
   return ( result ) ;
}




//----------------------------------------------------------------------
list< pair<size_t,double> >
FV_DiscreteField:: main_geometric_boundaries_on_proc( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: main_geometric_boundaries_on_proc" ) ;
   
   pair<size_t,double> geomB(0,0.);
   list< pair<size_t,double> > result ;
   list< FV_BoundaryCondition* >::const_iterator il;
   
   for ( il=SET_OF_BCS.begin();il!=SET_OF_BCS.end();il++)
     if ( FV_DomainBuilder::is_main_color( (*il)->get_color_ID() ) )
       if ( (*il)->has_DOF_on_proc() && !(*il)->is_periodic(0) )
       {
         pair<size_t,string> geoFeatures = FV_DomainBuilder:: 
	 	normal_direction_to_main_color( (*il)->get_color_ID() ) ;
	 geomB.first = geoFeatures.first ;
	 if ( geoFeatures.second == "min" )
	   geomB.second = PRIMARY_GRID
	   	->get_min_coordinate_on_current_processor( geomB.first ) ;
	 else
	   geomB.second = PRIMARY_GRID
	   	->get_max_coordinate_on_current_processor( geomB.first ) ;
	 result.push_back( geomB ) ;
       }
           
   return ( result ) ;    
}




//----------------------------------------------------------------------
void
FV_DiscreteField:: check_field_primary_meshes_coincide( 
	bool const &force, std::ostream& os, size_t indent_width )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: check_field_primary_meshes_coincide" ) ;

   size_t trans_dir = PRIMARY_GRID->get_translation_direction() ;
   double trans_dist = 
   	(*PRIMARY_GRID->get_global_main_coordinates())[trans_dir](0)
   	- (*global_main_coordinates)[0][trans_dir](0) ;

   if ( fabs( trans_dist ) > 1e-12 )
   {       
     if ( force )
     {
       std::string space( indent_width, ' ' ) ; 
       if ( MAC_Exec::communicator()->rank() == 0 ) 
       {
          os << space << "Translate field mesh of field " << FNAME
		<< " to coincide with primary grid" << endl ;
	  os << space << "   Translation direction = " << trans_dir <<
	  	"   Translation distance = " << trans_dist << endl;
       }	
       translate_field_mesh( trans_dir, trans_dist ) ;
     }
     else
       	FV_DiscreteField_ERROR:: n3( FNAME ) ;	
   }       
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField:: global_index_from_local( size_t i, 
      	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: global_index_from_local" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );     
   MAC_CHECK_PRE( i < (*UNK_LOCAL_NUMBERING)[ALL_COMPS_SAME_LOCATION == false ? 
   	component : 0].index_bound(direction) );

   size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;

   return ( i + (*local_min_index_in_global)[locNum](direction) );
}







//----------------------------------------------------------------------
double
FV_DiscreteField:: compute_boundary_cell_centered_DOF_integral( 
	size_t component, size_t level, std::string const& boundary_name ) 
	const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField:: compute_boundary_cell_centered_DOF_integral" ) ;

   double result = 0.;
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the FV_DiscreteField_Centered "
   	<< "or FV_DiscreteField_Staggered type; "
	<< "method \"compute_boundary_cell_centered_DOF_integral\" is "
	<< "implemented for a field of type FV_DiscreteField_Centered "
	<< "or FV_DiscreteField_Staggered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField:: compute_boundary_mean_normal_derivative( 
	size_t component, size_t level, std::string const& boundary_name ) 
	const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "FV_DiscreteField:: compute_boundary_mean_normal_derivative");

   size_t bc_color = FV_DomainBuilder::get_color_number( boundary_name ) ;   
   MAC_CHECK( FV_DomainBuilder::is_main_color( bc_color ) );   	   

   double result = (*V_BCS)[bc_color]
   	->compute_boundary_mean_normal_derivative( this, component, 
	level ) ;
            
   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField:: compute_grad_at_cell_center(
      		size_t i, size_t j,
		FV_SHIFT_TRIPLET shift,
		size_t component, size_t direction,
		double dxC, double dyC )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField:: compute_grad_at_cell_center" ) ;

   double result = 0.;
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the "
   	<< "FV_DiscreteField_Staggered type; "
	<< "method \"compute_grad_at_cell_center\" is "
	<< "implemented for a field of type "
	<< "FV_DiscreteField_Staggered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}



//----------------------------------------------------------------------
double
FV_DiscreteField:: compute_grad_at_cell_center(
      		size_t i, size_t j, size_t k,
		FV_SHIFT_TRIPLET shift,
		size_t component, size_t direction,
		double dxC, double dyC, double dzC )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField:: compute_grad_at_cell_center" ) ;

   double result = 0.;
   
   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the "
   	<< "FV_DiscreteField_Staggered type; "
	<< "method \"compute_grad_at_cell_center\" is "
	<< "implemented for a field of type "
	<< "FV_DiscreteField_Staggered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;   
       
   return ( result ) ;
}




//internal--------------------------------------------------------------
void
FV_DiscreteField_ERROR:: n1( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Translation-projection requires at least 2 levels of storage" 
   	<< " for field \"" << name << "\"";
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_DiscreteField_ERROR:: n2( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Translation-projection requires to create the interpolation"
   	<< " structure first; this has not been done" 
   	<< " for field \"" << name << "\"" << endl 
	<< "Check your implementation !!" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_DiscreteField_ERROR:: n3( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Primary and field meshes do not coincide"
   	<< " for field \"" << name << "\"" << endl 
	<< "Check your implementation !!" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_DiscreteField_ERROR:: n4( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Global numbering problem"
   	<< " for field \"" << name << "\"" << endl 
	<< "Check your implementation !!" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//----------------------------------------------------------------------------
// interpolate field values in 3D
double FV_DiscreteField:: interpolateFieldValues(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: interpolateFieldValues" ) ;	

  size_t dimens = 3, i, j, k ; 
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  	
  size_t_vector indices( dimens ) ;
  double result = 0. ;
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
            
  // Localize in structured mesh
  for (i=0;i<dimens;++i)
  {    
    FV_Mesh:: between( &(*local_main_coordinates)[locNum][i], point(i), 
    	indices(i) ) ;
    cellsize(i) = (*local_main_coordinates)[locNum][i](indices(i)+1) 
    	- (*local_main_coordinates)[locNum][i](indices(i)) ;	    
  }
  
  // Loop on surrounding nodes
  // Weight of nodes that contribute to the field constraint are computed
  // using multi-linear functions
  for (i=indices(0); i<indices(0)+2; ++i)
    for (j=indices(1); j<indices(1)+2; ++j)
      for (k=indices(2); k<indices(2)+2; ++k)
      {                    
	result += ( 1. - fabs( point(0) 
	- (*local_main_coordinates)[locNum][0](i) ) / cellsize(0) )
 	    	* ( 1. - fabs( point(1)
	- (*local_main_coordinates)[locNum][1](j) ) / cellsize(1) )
		* ( 1. - fabs( point(2)
	- (*local_main_coordinates)[locNum][2](k) ) / cellsize(2) )
		* (*VALUES)[level][component](i,j,k) ;
      }
  return result ;
  
}





//----------------------------------------------------------------------------
// reconstruct field values in 3D
double FV_DiscreteField:: GaussianReconstructionFieldValues(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level, 
	const double kernelWidth, const double dp ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: GaussianReconstructionFieldValues" ) ;	

  size_t dimens = 3; 
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. );
  
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  	
  size_t_vector indices( dimens ), stencil( dimens ) ;
  double result = 0.;
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
            
  // Localize in structured mesh (the grid closest to the particle center is
  // returned)
  for (size_t ii=0;ii<dimens;++ii)
  {    
    // By default "between" returns left bottom behind, but we want the closest point:
    if (point(ii) - (*local_main_coordinates)[locNum][ii](indices(ii)) > cellsize(ii)/2.)
      indices(ii) = indices(ii) + 1;
   // We re-estimate the delta x from the middle of the domain, to make sure 
   // we are not in the case of near boundary where grid is half size, because the 
   // kernel requires the full grid size 
    cellsize(ii) = (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2 + 1) 
    	- (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2) ;
    stencil(ii) = int(kernelWidth*dp/cellsize(ii) + 0.5);
  }

  
  // Loop on surrounding nodes
  // Weight of nodes that contribute to the field constraint are computed
  // using Gaussian function
  double GM=0., distance=0., cellMeasure=0.;
  int ip, jp, kp;
  double sigma = kernelWidth*dp/(2*sqrt(2*log(2)));
  for( int i=int(indices(0)-stencil(0)); i<=int(indices(0)+stencil(0)); i++ )
   for( int j=int(indices(1)-stencil(1)); j<=int(indices(1)+stencil(1)); j++ )
     for( int k=int(indices(2)-stencil(2)); k<=int(indices(2)+stencil(2)); k++ )
       if ( DOF_in_domain( i, j, k, locNum ) )
	if( DOF_is_unknown( i, j, k, locNum) )
	{
        cellMeasure = get_cell_measure(
       	    i, j, k, locNum ) ;
 	distance = sqrt(
	    pow(point(0)-(*local_main_coordinates)[locNum][0](i),2) +
	    pow(point(1)-(*local_main_coordinates)[locNum][1](j),2) +
	    pow(point(2)-(*local_main_coordinates)[locNum][2](k),2) );
	  GM = (exp(-pow(distance,2)/(2*pow(sigma,2))))/
	      	    pow(sigma*sqrt(2*acos(-1)),3);
	  result +=	GM*cellMeasure*(*VALUES)[level][locNum](i,j,k);
//DEBUG              SUM += cellMeasure*GM;
//DEBUG              no += 1;
	}
	else
	{
         // the contribution of grids on boundary are given to
	 // grids inside the domain in a symetric fashion 
	 ip = i;
	 jp = j;
	 kp = k;
	 vector<double>  offset( 3 );
	 if ( DOF_offset ( ip, jp, kp, indices, stencil, offset, 0 ) )
	 {
	 // GC - Point_out = GC - Point_in + (Point_in - Point_out) =
	 // GC - Point_in + offset
	 distance = sqrt(
	    pow(point(0)-(*local_main_coordinates)[locNum][0](ip)+offset[0],2) +
	    pow(point(1)-(*local_main_coordinates)[locNum][1](jp)+offset[1],2) +
	    pow(point(2)-(*local_main_coordinates)[locNum][2](kp)+offset[2],2) );
	 GM = (exp(-pow(distance,2)/(2*pow(sigma,2))))/
	      	    pow(sigma*sqrt(2*acos(-1)),3);
	 cellMeasure = get_cell_measure(ip, jp, kp, 0 ) ;
	 result +=	GM*cellMeasure*(*VALUES)[level][locNum](ip,jp,kp);
//DEBUG	     SUM += cellMeasure*GM;
//DEBUG	     no += 1;

	 }
	}  
       else
       {
         // the contribution of grids outside the domain are given to
	 // grids inside the domain in a symetric fashion 
	 ip = i;
	 jp = j;
	 kp = k;
	 vector<double>  offset( 3 );
	 if ( DOF_offset ( ip, jp, kp, indices, stencil, offset, 0 ) )
	 {
	 // GC - Point_out = GC - Point_in + (Point_in - Point_out) =
	 // GC - Point_in + offset
	 distance = sqrt(
            pow(point(0)-(*local_main_coordinates)[locNum][0](ip)
		    +offset[0],2) +
	    pow(point(1)-(*local_main_coordinates)[locNum][1](jp)
		    +offset[1],2) +
	    pow(point(2)-(*local_main_coordinates)[locNum][2](kp)
		    +offset[2],2) );
	 GM = (exp(-pow(distance,2)/(2*pow(sigma,2))))/
	      	    pow(sigma*sqrt(2*acos(-1)),3);
	 cellMeasure = get_cell_measure(ip, jp, kp, 0 ) ;
	 result +=	GM*cellMeasure*(*VALUES)[level][locNum](ip,jp,kp);
//	     SUM += cellMeasure*GM;
//	     no += 1;	
	 }
	 else
	 {
	  cout << "!!!WARNING!!! There is a bug in , \n"<<
	  "FV_DiscreteField:: GaussianReconstructionFieldValues" <<endl;
	  exit(0);
	 }

	}  
  return result ;
}





//----------------------------------------------------------------------------
// reconstruct field values in 3D
double FV_DiscreteField:: CorrectiveKernelAverageFieldValues(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level, 
	const double kernelWidth, const double dp ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: CorrectiveKernelAverageFieldValues" ) ;	

  size_t dimens = 3; 
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  	
  size_t_vector indices( dimens ), stencil( dimens ) ;
  double sum_nom = 0., sum_denom =0.;
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
  double sigma = kernelWidth*dp/(2*sqrt(2*log(2)));

  // Localize in structured mesh (the grid closest to the particle center is
  // returned)
  for (size_t ii=0;ii<dimens;++ii)
  {    
 // By default between returns left bottom behind, but we want the closest point:
 //   if (point(ii) - (*local_main_coordinates)[locNum][ii](indices(ii)) > cellsize(ii)/2.)
 //     indices(ii) = indices(ii) + 1;
 // We re-estimate the delta x from the middle of the domain, to make sure 
 // we are not in the case of near boundary where grid is half size, because the 
 // kernel requires the full grid size 
    cellsize(ii) = (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2 + 1) 
    	- (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2) ;
	
    if(  point(ii)+cellsize(ii)/2.<(*local_main_coordinates)[locNum][ii](
       (*max_index_unknown_on_proc)[locNum](ii)) )
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][ii], point(ii)+
      cellsize(ii)/2., indices(ii) ) ;
    else
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][ii], point(ii),
      	 indices(ii) ) ;	    
    stencil(ii) = int(kernelWidth*dp/cellsize(ii) + 0.5);
  }
  // Loop on surrounding nodes
  // Weight of nodes that contribute to the field constraint are computed
  // using CorrectiveKernelAverage
  double ker=0., r=0., cellMeasure=0.;
  for( int i=int(indices(0)-stencil(0)); i<=int(indices(0)+stencil(0)); i++ )
   for( int j=int(indices(1)-stencil(1)); j<=int(indices(1)+stencil(1)); j++ )
     for( int k=int(indices(2)-stencil(2)); k<=int(indices(2)+stencil(2)); k++ )
       if ( DOF_in_domain( i, j, k, locNum ) )
	if( DOF_is_unknown( i, j, k, locNum) )
	{
	  cellMeasure = get_cell_measure(
       	      i, j, k, locNum ) ;
	  r = sqrt(
	      pow(point(0)-(*local_main_coordinates)[locNum][0](i),2) +
	      pow(point(1)-(*local_main_coordinates)[locNum][1](j),2) +
	      pow(point(2)-(*local_main_coordinates)[locNum][2](k),2) );
	  ker = (exp(-pow(r,2)/(2*pow(sigma,2))))/
	      	    pow(sigma*sqrt(2*acos(-1)),3);
	  sum_nom += ker*cellMeasure*(*VALUES)[level][locNum](i,j,k);
	  sum_denom += ker*cellMeasure;
	}	
  return sum_nom/sum_denom ;
}




//----------------------------------------------------------------------------
// interpolate gradient field values in 3D by using two-cubes averaged on 
// both sides of the particle. This method is not used by default, since the
// corrective kernel method (method below) is a better choice. however, 
// owing to its simplicity I (Amir) do not delete it, since it might be useful 
// for future works (non-spherical particles)
 geomVector FV_DiscreteField:: interpolateGradientWithKernel(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level, 
	const double kernelWidth, const double dp ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: interpolateGradientWithKernel" ) ;
  size_t dimens = 3; 
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  	
  size_t_vector indices( dimens ) ;
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
  geomVector result( dimens ), PosLBB(dimens), PosRTF(dimens),
           cellNbLBB(dimens), cellNbRTF(dimens),ResLBB(dimens),ResRTF(dimens);
  //double SUM = 0.;
  //size_t no = 0;
  
            
  // Localize in structured mesh
  for (size_t ii=0;ii<dimens;++ii)
  { 
    cellsize(ii) = (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2 + 1) 
    	- (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2) ;
    if(  point(ii)+cellsize(ii)/2.<(*local_main_coordinates)[locNum][ii](
       (*max_index_unknown_on_proc)[locNum](ii)) )
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][ii], point(ii)+
      cellsize(ii)/2., indices(ii) ) ;
    else
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][ii], point(ii),
      	 indices(ii) ) ;    	    
  }
  
  // Loop on surrounding nodes
  // Weight of nodes that contribute to the field constraint are computed
  // using Gaussian function
      double ratio = dp/cellsize(0);
      int stencil = int(kernelWidth*ratio + 0.5);
      for( int i=int(indices(0))-stencil; i<=int(indices(0))+stencil; i++ )
       for( int j=int(indices(1))-stencil; j<=int(indices(1))+stencil; j++ )
         for( int k=int(indices(2))-stencil; k<=int(indices(2))+stencil; k++ )
	   if ( DOF_in_domain( i, j, k, locNum ) )
	    if( DOF_is_unknown( i, j, k, locNum) )
	    { 
	      if ( i <= indices(0) )
	      {
                cellNbLBB(0) += 1.;
	        PosLBB(0) += (*local_main_coordinates)[locNum][0](i);
		ResLBB(0) += (*VALUES)[level][locNum](i,j,k);
	      }	
	      if ( i >= indices(0) )
	      {
	        cellNbRTF(0) += 1.;
		PosRTF(0) += (*local_main_coordinates)[locNum][0](i);
		ResRTF(0) += (*VALUES)[level][locNum](i,j,k);
	      }
	      if (j <= indices(1) )
	      {
	        cellNbLBB(1) += 1.;
	        PosLBB(1) += (*local_main_coordinates)[locNum][1](j);
		ResLBB(1) += (*VALUES)[level][locNum](i,j,k);
	      }
	      if (j >= indices(1) )
	      {
	        cellNbRTF(1) += 1.;
		PosRTF(1) += (*local_main_coordinates)[locNum][1](j);
		ResRTF(1) += (*VALUES)[level][locNum](i,j,k);
	      }
	      if (k <= indices(2) )
	      {
	        cellNbLBB(2) += 1.;
		PosLBB(2) += (*local_main_coordinates)[locNum][2](k);
		ResLBB(2) += (*VALUES)[level][locNum](i,j,k);
	      }
	      if (k >= indices(2) ) 
	      {
	        cellNbRTF(2) += 1.;
		PosRTF(2) += (*local_main_coordinates)[locNum][2](k);
		ResRTF(2) += (*VALUES)[level][locNum](i,j,k);
	      }
	      
	    }
      PosLBB(0) = PosLBB(0)/cellNbLBB(0);
      PosLBB(1) = PosLBB(1)/cellNbLBB(1);
      PosLBB(2) = PosLBB(2)/cellNbLBB(2);
      PosRTF(0) = PosRTF(0)/cellNbRTF(0);
      PosRTF(1) = PosRTF(1)/cellNbRTF(1);
      PosRTF(2) = PosRTF(2)/cellNbRTF(2);
      
      ResLBB(0) = ResLBB(0)/cellNbLBB(0);
      ResLBB(1) = ResLBB(1)/cellNbLBB(1);
      ResLBB(2) = ResLBB(2)/cellNbLBB(2);
      ResRTF(0) = ResRTF(0)/cellNbRTF(0);
      ResRTF(1) = ResRTF(1)/cellNbRTF(1);
      ResRTF(2) = ResRTF(2)/cellNbRTF(2);
      
      result(0) = (ResRTF(0) - ResLBB(0))/(PosRTF(0) - PosLBB(0));
      result(1) = (ResRTF(1) - ResLBB(1))/(PosRTF(1) - PosLBB(1));
      result(2) = (ResRTF(2) - ResLBB(2))/(PosRTF(2) - PosLBB(2));
      
//      if ( result(0) > 100000. || result(1) > 100000. || result(2) > 100000.)
//      {
//        cout <<"TOTO x  result(0) = "<<result(0)<<"  result(1) = "<<result(1)
//	<<"  result(2) = "<<result(2)<< endl;
//	cout <<"point 1 2 3 = "<<point(0)<<" "<<point(1)<<" "<<point(2)<<endl;
//	cout <<"indice 1 2 3 = "<<indices(0)<<" "<<indices(1)<<" "<<indices(2)<<endl;
//	cout <<"stencil = "<<stencil<<endl;
//     }	

//     if ( cellNbLBB(1) != cellNbRTF(1) )
//     {
//       cout <<"TOTO y  cellNbLBB(1) = "<<cellNbLBB(1)<<"  cellNbRTF(1) = "<<cellNbRTF(1) << endl;
//	cout <<"point 1 2 3 = "<<point(0)<<" "<<point(1)<<" "<<point(2)<<endl;
//	cout <<"indice 1 2 3 = "<<indices(0)<<" "<<indices(1)<<" "<<indices(2)<<endl;
      
      
	      
  return result ;
} 




//----------------------------------------------------------------------------
// interpolate gradient field values in 3D with a corrective kernel
// c.f. Amir's thesis for detail (in annex)
geomVector FV_DiscreteField:: CorrectiveKernelGradientFieldValues(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level, 
	const double kernelWidth, const double dp ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: CorrectiveKernelGradientFieldValues" ) ;
  size_t dimens = 3; 
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ), dist( dimens, 0. ),
  	derker( dimens, 0. ) ;
  geomVector result( dimens );
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  	
  size_t_vector indices( dimens ), stencil( dimens ) ;
  vector< vector<double> > a( dimens,vector<double>(dimens) );
  vector<double> b(dimens);
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
            
  double sigma = kernelWidth*dp/(2*sqrt(2*log(2)));  

   // Localize in structured mesh (the grid closest to the particle center is
  // returned)
  for (size_t ii=0;ii<dimens;++ii)
  {    
    cellsize(ii) = (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2 + 1) 
    	- (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2) ;
	
    if(  point(ii)+cellsize(ii)/2.<(*local_main_coordinates)[locNum][ii](
       (*max_index_unknown_on_proc)[locNum](ii)) )
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][ii], point(ii)+
      cellsize(ii)/2., indices(ii) ) ;
    else
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][ii], point(ii),
      	 indices(ii) ) ;	    
    stencil(ii) = int(kernelWidth*dp/cellsize(ii) + 0.5);
  }

  // Loop on surrounding nodes
  // Weight of nodes that contribute to the field constraint are computed
  // using CorrectiveKernelAverage
  double r=0., cellMeasure=0., det = 0., GM = 0.;
  for( int i=int(indices(0)-stencil(0)); i<=int(indices(0)+stencil(0)); i++ )
   for( int j=int(indices(1)-stencil(1)); j<=int(indices(1)+stencil(1)); j++ )
     for( int k=int(indices(2)-stencil(2)); k<=int(indices(2)+stencil(2)); k++ )
       if ( DOF_in_domain( i, j, k, locNum ) )
	if( DOF_is_unknown( i, j, k, locNum) )
	{
	  cellMeasure = get_cell_measure(i, j, k, locNum ) ;
	  //The kernel is centered at the position of the closest grid and not
	  //the exact particle position. The reason is that the value of the
	  //kernel at the center is required and we want to avoid interpolation
	  r = sqrt(
	      pow((*local_main_coordinates)[locNum][0](indices(0))-
	      (*local_main_coordinates)[locNum][0](i),2) +
	      pow((*local_main_coordinates)[locNum][1](indices(1))-
	      (*local_main_coordinates)[locNum][1](j),2) +
	      pow((*local_main_coordinates)[locNum][2](indices(2))-
	      (*local_main_coordinates)[locNum][2](k),2) );
	  dist(0) = (*local_main_coordinates)[locNum][0](i) - 
	  	(*local_main_coordinates)[locNum][0](indices(0));
	  dist(1) = (*local_main_coordinates)[locNum][1](j) - 
	  	(*local_main_coordinates)[locNum][1](indices(1));
	  dist(2) = (*local_main_coordinates)[locNum][2](k) - 
	  	(*local_main_coordinates)[locNum][2](indices(2));  
	  GM = (exp(-pow(r,2)/(2*pow(sigma,2))))/
	      	    pow(sigma*sqrt(2*acos(-1)),3);

          derker(0) = (-dist(0)/pow(sigma,2))*GM; 
          derker(1) = (-dist(1)/pow(sigma,2))*GM; 
          derker(2) = (-dist(2)/pow(sigma,2))*GM; 

	  a[0][0] += dist(0) * derker(0) * cellMeasure;
	  a[0][1] += dist(1) * derker(0) * cellMeasure;
	  a[0][2] += dist(2) * derker(0) * cellMeasure;
	  
	  a[1][0] += dist(0) * derker(1) * cellMeasure;
	  a[1][1] += dist(1) * derker(1) * cellMeasure;
	  a[1][2] += dist(2) * derker(1) * cellMeasure;
	  
	  a[2][0] += dist(0) * derker(2) * cellMeasure;
	  a[2][1] += dist(1) * derker(2) * cellMeasure;
	  a[2][2] += dist(2) * derker(2) * cellMeasure;
	  
	  b[0] += ( (*VALUES)[level][locNum](i,j,k) - 
	  	(*VALUES)[level][locNum](indices(0),indices(1),indices(2)) )
		 * derker(0) * cellMeasure;
	  b[1] += ( (*VALUES)[level][locNum](i,j,k) - 
	  	(*VALUES)[level][locNum](indices(0),indices(1),indices(2)) )
		 * derker(1) * cellMeasure;
	  b[2] += ( (*VALUES)[level][locNum](i,j,k) - 
	  	(*VALUES)[level][locNum](indices(0),indices(1),indices(2)) )
		 * derker(2) * cellMeasure;
	}
  det = a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]
    -a[2][0]*a[1][1]*a[0][2]-a[2][1]*a[1][2]*a[0][0]-a[2][2]*a[1][0]*a[0][1];
    
  result(0) = ( (a[2][2]*a[1][1]-a[1][2]*a[2][1])*b[0] 
          +(a[0][2]*a[2][1]-a[2][2]*a[0][1])*b[1] 
	  +(a[0][1]*a[1][2]-a[0][2]*a[1][1])*b[2] )/det;
  result(1) = ( (a[1][2]*a[2][0]-a[1][0]*a[2][2])*b[0] 
	  +(a[0][0]*a[2][2]-a[0][2]*a[2][0])*b[1] 
	  +(a[0][2]*a[1][0]-a[0][0]*a[1][2])*b[2] )/det;
  result(2) = ( (a[1][0]*a[2][1]-a[1][1]*a[2][0])*b[0] 
	  +(a[0][1]*a[2][0]-a[0][0]*a[2][1])*b[1] 
	  +(a[0][0]*a[1][1]-a[0][1]*a[1][0])*b[2] )/det;
  	    
  return result ;
} 





//----------------------------------------------------------------------------
// interpolate gradient field values in 3D 
// A simple method in estimating Laplacian of a value in DEMCFD GK method
// (Epsilon computing method = 3). It is originally implemented only for 
// post-processing & tentative improvements of the method (c.f. IJMF2016 Amir) 
// It might be used in future. 
 double FV_DiscreteField:: interpolateLaplacianWithKernel(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level, 
	const double kernelWidth, const double dp ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: interpolateLaplacianWithKernel" ) ;
  size_t dimens = 3; 
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  	
  size_t_vector indices( dimens ) ;
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
  geomVector result(dimens), PosLBB(dimens), PosRTF(dimens), PosCenter(dimens),
           cellNbLBB(dimens), cellNbRTF(dimens),cellNbCenter(dimens),
	   ResLBB(dimens),ResRTF(dimens),ResCenter(dimens);
  double result_sum;	   
  //double SUM = 0.;
  //size_t no = 0;
  
  double ratio;      
  int stencil;     
  // Localize in structured mesh
  for (size_t ii=0;ii<dimens;++ii)
  { 
    cellsize(ii) = (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2 + 1) 
    	- (*local_main_coordinates)[locNum][ii](
       ((*max_index_unknown_handled_by_proc)[locNum](ii)+
       (*min_index_unknown_handled_by_proc)[locNum](ii))/2) ;
       
    FV_Mesh:: between( &(*local_main_coordinates)[locNum][ii], point(ii),
      	 indices(ii) ) ;    
    ratio = dp/cellsize(ii);
    stencil = int(kernelWidth*ratio + 0.5) + 2;	 
    if (int(indices(ii) -stencil+ 
    	(*local_min_index_in_global)[locNum](ii)) <= 0 )
     indices(ii) = indices(ii) + (*local_min_index_in_global)[locNum](ii)-
      		  (indices(ii) - stencil) + 1 ; 
    else if (int(indices(ii)+stencil +
	        (*local_min_index_in_global)[locNum](ii))>= 
    		(*global_max_index)[locNum](ii))
     indices(ii) = indices(ii) -((indices(ii) + stencil) + 1+
      		(*local_min_index_in_global)[locNum](ii)-
		(*global_max_index)[locNum](ii));
  }
  
  
  // Loop on surrounding nodes
  // Weight of nodes that contribute to the field constraint are computed
  // using Gaussian function
  ratio = dp/cellsize(0);
  stencil = int(kernelWidth*ratio + 0.5) + 2;
  size_t counter = 0;


  for( int i=int(indices(0))-stencil; i<=int(indices(0))+stencil; i++ )
   for( int j=int(indices(1))-stencil; j<=int(indices(1))+stencil; j++ )
     for( int k=int(indices(2))-stencil; k<=int(indices(2))+stencil; k++ )
       if ( DOF_in_domain( i, j, k, locNum ) )
	if( DOF_is_unknown( i, j, k, locNum) )
	{  
	  counter++;
	  if ( i < indices(0) - 1 )
	  {
            cellNbLBB(0) += 1.;
	    PosLBB(0) += (*local_main_coordinates)[locNum][0](i);
	    ResLBB(0) += (*VALUES)[level][locNum](i,j,k);
	  }	
	  else if ( i >= indices(0) -1 && i <= indices(0) + 1 )
	  {
	    cellNbCenter(0) += 1.;			     
	    PosCenter(0) += (*local_main_coordinates)[locNum][0](i);
	    ResCenter(0) += (*VALUES)[level][locNum](i,j,k); 
	  }
	  else if ( i > indices(0) + 1 )
	  {
	    cellNbRTF(0) += 1.;
	    PosRTF(0) += (*local_main_coordinates)[locNum][0](i);
	    ResRTF(0) += (*VALUES)[level][locNum](i,j,k);
	  }
	  
	  if (j < indices(1) - 1 )
	  {
	    cellNbLBB(1) += 1.;
	    PosLBB(1) += (*local_main_coordinates)[locNum][1](j);
	    ResLBB(1) += (*VALUES)[level][locNum](i,j,k);
	  }
	  else if (j >= indices(1) - 1 && j <= indices(1) + 1)
	  {
	    cellNbCenter(1) += 1.;	 
	    PosCenter(1) += (*local_main_coordinates)[locNum][1](j);
 	    ResCenter(1) += (*VALUES)[level][locNum](i,j,k); 
	  }
	  else if (j > indices(1) + 1)
	  {
	    cellNbRTF(1) += 1.;
	    PosRTF(1) += (*local_main_coordinates)[locNum][1](j);
	    ResRTF(1) += (*VALUES)[level][locNum](i,j,k);
	  }
	  
	  if (k < indices(2) - 1 )
	  {
	    cellNbLBB(2) += 1.;
	    PosLBB(2) += (*local_main_coordinates)[locNum][2](k);
	    ResLBB(2) += (*VALUES)[level][locNum](i,j,k);
	  }
	  else if (j >= indices(1) - 1 && j <= indices(1) + 1)
	  {
	    cellNbCenter(2) += 1.;	 
	    PosCenter(2) += (*local_main_coordinates)[locNum][2](k);
 	    ResCenter(2) += (*VALUES)[level][locNum](i,j,k); 
	  }
	  else if (k > indices(2) + 1) 
	  {
	    cellNbRTF(2) += 1.;
	    PosRTF(2) += (*local_main_coordinates)[locNum][2](k);
	    ResRTF(2) += (*VALUES)[level][locNum](i,j,k);
	  }
	}		      

  PosLBB(0) = PosLBB(0)/cellNbLBB(0);
  PosLBB(1) = PosLBB(1)/cellNbLBB(1);
  PosLBB(2) = PosLBB(2)/cellNbLBB(2);
  
  PosCenter(0) = PosCenter(0)/cellNbCenter(0);
  PosCenter(1) = PosCenter(1)/cellNbCenter(1);
  PosCenter(2) = PosCenter(2)/cellNbCenter(2);
  
  PosRTF(0) = PosRTF(0)/cellNbRTF(0);
  PosRTF(1) = PosRTF(1)/cellNbRTF(1);
  PosRTF(2) = PosRTF(2)/cellNbRTF(2);

  ResLBB(0) = ResLBB(0)/cellNbLBB(0);
  ResLBB(1) = ResLBB(1)/cellNbLBB(1);
  ResLBB(2) = ResLBB(2)/cellNbLBB(2);
  
  ResCenter(0) = ResCenter(0)/cellNbCenter(0);
  ResCenter(1) = ResCenter(1)/cellNbCenter(1);
  ResCenter(2) = ResCenter(2)/cellNbCenter(2);
  
  ResRTF(0) = ResRTF(0)/cellNbRTF(0);
  ResRTF(1) = ResRTF(1)/cellNbRTF(1);
  ResRTF(2) = ResRTF(2)/cellNbRTF(2);
  
  result(0) = ((ResRTF(0) - ResCenter(0))/(PosRTF(0) - PosCenter(0))
  	- (ResCenter(0) - ResLBB(0))/(PosCenter(0) - PosLBB(0)) )/
	((PosRTF(0) - PosLBB(0))/2.);
  result(1) = ((ResRTF(1) - ResCenter(1))/(PosRTF(1) - PosCenter(1))
  	- (ResCenter(1) - ResLBB(1))/(PosCenter(1) - PosLBB(1)) )/
	((PosRTF(1) - PosLBB(1))/2.);
  result(2) = ((ResRTF(2) - ResCenter(2))/(PosRTF(2) - PosCenter(2))
  	- (ResCenter(2) - ResLBB(2))/(PosCenter(2) - PosLBB(2)) )/
	((PosRTF(2) - PosLBB(2))/2.);	

  result_sum = result(0) + result(1) + result(2) ;

  return result_sum ;
} 





//----------------------------------------------------------------------------
// interpolate field values in 2D for 1 component
double FV_DiscreteField:: interpolateFieldValues(
	const double &X_coordinate,
	const double &Y_coordinate,
	size_t component, size_t level ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: interpolateFieldValues" ) ;	

  size_t dimens = 2, i, j, k = 0 ;
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  size_t_vector indices( dimens ) ;
  double result = 0. ;
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
            
  // Localize in structured mesh
  for (i=0;i<dimens;++i)
  {    
    FV_Mesh:: between( &(*local_main_coordinates)[locNum][i], point(i), 
    	indices(i) ) ;    
    cellsize(i) = (*local_main_coordinates)[locNum][i](indices(i)+1) 
    	- (*local_main_coordinates)[locNum][i](indices(i)) ;	    
  }
  
  // Loop on surrounding nodes
  // Weight of nodes that contribute to the field constraint are computed
  // using multi-linear functions
  for (i=indices(0); i<indices(0)+2; ++i)
    for (j=indices(1); j<indices(1)+2; ++j)
      result += ( 1. - fabs( point(0) 
      	      - (*local_main_coordinates)[locNum][0](i) ) / cellsize(0) )
      		* ( 1. - fabs( point(1)
	      - (*local_main_coordinates)[locNum][1](j) ) / cellsize(1) )
		* (*VALUES)[level][component](i,j,k) ;
		
  return result ;
}




//----------------------------------------------------------------------------
// compute gradient in a cell in 3D
geomVector FV_DiscreteField:: interpolateGradient(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: interpolateGradient" ) ;	

  size_t dimens = 3, i, j, k, ii=0, jj=0, kk=0 ;
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  
  size_t_vector indices( dimens ) ;
  geomVector result( dimens );
  geomVector X_interpolation( dimens-1 );
  geomVector Y_interpolation( dimens-1 );
  geomVector Z_interpolation( dimens-1 );
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;

  // Localize in structured mesh
  for (i=0;i<dimens;++i)
  {    
    FV_Mesh:: between( &(*local_main_coordinates)[locNum][i], point(i), 
    	indices(i) ) ;    
    cellsize(i) = (*local_main_coordinates)[locNum][i](indices(i)+1) 
    	- (*local_main_coordinates)[locNum][i](indices(i)) ;	
  }

  
  // Interpolate field values on both X faces
  for (i=indices(0); i<indices(0)+2; ++i,++ii)
    for (j=indices(1); j<indices(1)+2; ++j)
      for (k=indices(2); k<indices(2)+2; ++k)
        X_interpolation(ii) += ( 1. - fabs( point(1) 
        	- (*local_main_coordinates)[locNum][1](j) ) / cellsize(1) )
        	* ( 1. - fabs( point(2)
        	- (*local_main_coordinates)[locNum][2](k) ) / cellsize(2) )
        	* (*VALUES)[level][component](i,j,k) ;
    
  // Interpolate field values on both Y faces
  for (j=indices(1); j<indices(1)+2; ++j,++jj)
    for (i=indices(0); i<indices(0)+2; ++i)
      for (k=indices(2); k<indices(2)+2; ++k)
        Y_interpolation(jj) += ( 1. - fabs( point(0) 
        	- (*local_main_coordinates)[locNum][0](i) ) / cellsize(0) )
        	* ( 1. - fabs( point(2)
        	- (*local_main_coordinates)[locNum][2](k) ) / cellsize(2) )
        	* (*VALUES)[level][component](i,j,k) ;


  // Interpolate field values on both Z faces
  for (k=indices(2); k<indices(2)+2; ++k,++kk)
    for (i=indices(0); i<indices(0)+2; ++i)
      for (j=indices(1); j<indices(1)+2; ++j)
        Z_interpolation(kk) += ( 1. - fabs( point(0) 
        	- (*local_main_coordinates)[locNum][0](i) ) / cellsize(0) )
        	* ( 1. - fabs( point(1)
        	- (*local_main_coordinates)[locNum][1](j) ) / cellsize(1) )
        	* (*VALUES)[level][component](i,j,k) ;


              
  result(0)=(X_interpolation(1)-X_interpolation(0)) / cellsize(0);
//        cout <<"result(0)=(X_interpolation(1)-X_interpolation(0)) / cellsize(0);
  result(1)=(Y_interpolation(1)-Y_interpolation(0)) / cellsize(1);
  result(2)=(Z_interpolation(1)-Z_interpolation(0)) / cellsize(2);
  
  return result ;
}




//----------------------------------------------------------------------------
// compute gradient in a cell in 2D
geomVector FV_DiscreteField:: interpolateGradient(
	const double &X_coordinate,
	const double &Y_coordinate,
	size_t component, size_t level ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: interpolateGradient" ) ;	

  size_t dimens = 2, i, j, k = 0, ii = 0, jj = 0 ; 
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  size_t_vector indices( dimens ) ;
  geomVector result(dimens);
  geomVector X_interpolation( dimens );
  geomVector Y_interpolation( dimens );
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
            
  // Localize in structured mesh
  for (i=0;i<dimens;++i)
  {
    FV_Mesh:: between( &(*local_main_coordinates)[locNum][i], point(i), 
    	indices(i) ) ;    
    cellsize(i) = (*local_main_coordinates)[locNum][i](indices(i)+1) 
    	- (*local_main_coordinates)[locNum][i](indices(i)) ;	    
  }
  
  for (i=indices(0); i<indices(0)+2; ++i,++ii)
  {
    jj = 0 ;
    for (j=indices(1); j<indices(1)+2; ++j,++jj)
    {
      // Interpolation on both Y faces
      Y_interpolation(ii) += ( 1. - fabs( point(1) 
	- (*local_main_coordinates)[locNum][1](j) ) / cellsize(1) )
	* (*VALUES)[level][0](i,j,k) ;
		
      // Interpolation on both X faces
      X_interpolation(jj) += ( 1. - fabs( point(0) 
	- (*local_main_coordinates)[locNum][0](i) ) / cellsize(0) )
	* (*VALUES)[level][0](i,j,k) ;
    }
  }  
  
  result(0)=(Y_interpolation(1)-Y_interpolation(0)) / cellsize(0);
  result(1)=(X_interpolation(1)-X_interpolation(0)) / cellsize(1);
  
  return result ;
}




//----------------------------------------------------------------------------
// interpolate gradient field values in 3D 
// A simple method in estimating Laplacian of a value in DEMCFD BC method
// (Epsilon computing method = 1). It is originally implemented only for 
// post-processing & tentative improvements of the method (c.f. IJMF2016 Amir) 
// It might be used in future. 
double FV_DiscreteField:: interpolateLap(
	const double &X_coordinate,
	const double &Y_coordinate,
	const double &Z_coordinate,
	size_t component, size_t level ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: interpolateLap" ) ;	

  size_t dimens = 3, i, j, k, ii=0, jj=0, kk=0 ;
  doubleVector point( dimens, 0. ), cellsize( dimens, 0. ) ;
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  point(2) = Z_coordinate ;  
  size_t_vector indices( dimens ) ;
  geomVector result( dimens );
  geomVector X_interpolation( dimens );
  geomVector Y_interpolation( dimens );
  geomVector Z_interpolation( dimens );
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;
  size_t stencil = 1;
  double result_sum;
  // Localize in structured mesh
  for (i=0;i<dimens;++i)
  {    
    cellsize(i) = (*local_main_coordinates)[locNum][i](
       ((*max_index_unknown_handled_by_proc)[locNum](i)+
       (*min_index_unknown_handled_by_proc)[locNum](i))/2 + 1) 
    	- (*local_main_coordinates)[locNum][i](
       ((*max_index_unknown_handled_by_proc)[locNum](i)+
       (*min_index_unknown_handled_by_proc)[locNum](i))/2) ;
    if(  point(i)+cellsize(i)/2.<(*local_main_coordinates)[locNum][i](
       (*max_index_unknown_on_proc)[locNum](i)) )
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][i], point(i)+
      cellsize(i)/2., indices(i) ) ;
    else
      FV_Mesh:: between( &(*local_main_coordinates)[locNum][i], point(i),
      	 indices(i) ) ; 	
   
	 
    if (int(indices(ii) -stencil+ 
    	(*local_min_index_in_global)[locNum](ii)) <= 0 )
     indices(ii) = indices(ii) + (*local_min_index_in_global)[locNum](ii)-
      		  (indices(ii) - stencil) + 1 ; 
    else if (int(indices(ii)+stencil +
	        (*local_min_index_in_global)[locNum](ii))>= 
    		(*global_max_index)[locNum](ii))
     indices(ii) = indices(ii) -((indices(ii) + stencil) + 1+
      		(*local_min_index_in_global)[locNum](ii)-
		(*global_max_index)[locNum](ii));
  }

  
  // Interpolate field values on both X faces
  for (i=indices(0)-stencil; i<indices(0)+2; ++i,++ii)
    for (j=indices(1)-stencil; j<indices(1)+2; ++j)
      for (k=indices(2)-stencil; k<indices(2)+2; ++k)
        X_interpolation(ii) += ( 1. - fabs( point(1) 
        	- (*local_main_coordinates)[locNum][1](j) ) / cellsize(1) )
        	* ( 1. - fabs( point(2)
        	- (*local_main_coordinates)[locNum][2](k) ) / cellsize(2) )
        	* (*VALUES)[level][component](i,j,k) ;
    
  // Interpolate field values on both Y faces
  for (j=indices(1)-stencil; j<indices(1)+2; ++j,++jj)
    for (i=indices(0)-stencil; i<indices(0)+2; ++i)
      for (k=indices(2)-stencil; k<indices(2)+2; ++k)
        Y_interpolation(jj) += ( 1. - fabs( point(0) 
        	- (*local_main_coordinates)[locNum][0](i) ) / cellsize(0) )
        	* ( 1. - fabs( point(2)
        	- (*local_main_coordinates)[locNum][2](k) ) / cellsize(2) )
        	* (*VALUES)[level][component](i,j,k) ;


  // Interpolate field values on both Z faces
  for (k=indices(2)-stencil; k<indices(2)+2; ++k,++kk)
    for (i=indices(0)-stencil; i<indices(0)+2; ++i)
      for (j=indices(1)-stencil; j<indices(1)+2; ++j)
        Z_interpolation(kk) += ( 1. - fabs( point(0) 
        	- (*local_main_coordinates)[locNum][0](i) ) / cellsize(0) )
        	* ( 1. - fabs( point(1)
        	- (*local_main_coordinates)[locNum][1](j) ) / cellsize(1) )
        	* (*VALUES)[level][component](i,j,k) ;


  result(0)=(X_interpolation(2)-X_interpolation(0)) / pow(cellsize(0),2);
  result(1)=(Y_interpolation(2)-Y_interpolation(0)) / pow(cellsize(1),2);
  result(2)=(Z_interpolation(2)-Z_interpolation(0)) / pow(cellsize(2),2);
  
  result_sum = result(0) + result(1) + result(2);
  return result_sum ;
}




//----------------------------------------------------------------------------
// compute cell size 2D in a cell in 2D
geomVector FV_DiscreteField:: cellsize2D(
	const double &X_coordinate,
	const double &Y_coordinate,
	size_t component, size_t level ) const
//----------------------------------------------------------------------------
{
  MAC_LABEL( "FV_DiscreteField:: cellsize2D" ) ;	

  size_t dimens = 2, i;
  doubleVector point( dimens, 0. );
  point(0) = X_coordinate ;
  point(1) = Y_coordinate ;  
  //point(2) = Z_coordinate ;  
  size_t_vector indices( dimens ) ;
  geomVector result(dimens);
  size_t locNum = ALL_COMPS_SAME_LOCATION == false ? component : 0 ;

  // Localize in structured mesh
  for (i=0;i<dimens;++i)
  {    
    FV_Mesh:: between( &(*local_main_coordinates)[locNum][i], point(i), 
    	indices(i) ) ;    
    result(i) = (*local_main_coordinates)[locNum][i](indices(i)+1) 
    	- (*local_main_coordinates)[locNum][i](indices(i)) ;	
  }

  return result;

}
