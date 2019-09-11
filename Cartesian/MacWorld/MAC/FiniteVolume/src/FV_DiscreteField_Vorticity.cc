#include <FV_DiscreteField_Vorticity.hh>
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


FV_DiscreteField_Vorticity const* 
FV_DiscreteField_Vorticity::PROTOTYPE = new FV_DiscreteField_Vorticity( 
	"vorticity" ) ;


//----------------------------------------------------------------------
FV_DiscreteField_Vorticity:: FV_DiscreteField_Vorticity( 
	std::string const& a_type )
//----------------------------------------------------------------------	
   : FV_DiscreteField( a_type ) 
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: FV_DiscreteField_Vorticity") ;
   
   ALL_COMPS_SAME_LOCATION = false ;   
}




//----------------------------------------------------------------------
FV_DiscreteField_Vorticity:: FV_DiscreteField_Vorticity( 
	MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type )
//----------------------------------------------------------------------	
   : FV_DiscreteField( a_owner, a_primary_mesh, a_name, a_type ) 
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: FV_DiscreteField_Vorticity") ;

   ALL_COMPS_SAME_LOCATION = false ;
}




//----------------------------------------------------------------------
FV_DiscreteField*
FV_DiscreteField_Vorticity:: create_replica( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type,
	size_t a_nb_components,
	size_t a_depth ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: create_replica" ) ;

   FV_DiscreteField* result = new FV_DiscreteField_Vorticity( a_owner,
	a_primary_mesh, a_name, a_type, a_nb_components, a_depth ) ;
    
   return( result ) ;
}




//----------------------------------------------------------------------
FV_DiscreteField*
FV_DiscreteField_Vorticity:: create_clone_replica( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: create_clone_replica" ) ;

   FV_DiscreteField* result = new FV_DiscreteField_Vorticity( a_owner, 
   	a_primary_mesh, a_name, a_type ) ;
    
   return( result ) ;
}




//----------------------------------------------------------------------
FV_DiscreteField_Vorticity:: FV_DiscreteField_Vorticity( MAC_Object* a_owner,
	FV_Mesh const* a_primary_mesh,
	std::string const& a_name,
	std::string const& a_type,
	size_t a_nb_components,
	size_t a_depth )
//----------------------------------------------------------------------
    : FV_DiscreteField( a_owner, a_primary_mesh, a_name, a_type, 
    	a_nb_components, a_depth )
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: FV_DiscreteField_Vorticity") ;
   MAC_ASSERT( IMPLIES( DIM == 2, NB_COMPS == 1 ) ) ; 
   MAC_ASSERT( IMPLIES( DIM == 3, NB_COMPS == 3 ) ) ;

   ALL_COMPS_SAME_LOCATION = false ;
   SET_BC_VALUES_ALLOWED = true ; 
   vector< boolVector > const* primary_mesh_handles_node = 
   	PRIMARY_GRID->get_nodes_owner();
   size_t security_bandwidth = PRIMARY_GRID->get_security_bandwidth();
      
   // Allocate index vectors
   size_t_vector work_size_t( DIM, 0 ) ;
   global_max_index = new vector<size_t_vector>( NB_COMPS, work_size_t ) ;
   local_min_index_in_global = new vector<size_t_vector>( NB_COMPS, 
   	work_size_t ) ;
   local_max_index_in_global = new vector<size_t_vector>( NB_COMPS, 
   	work_size_t ) ;
   local_dof_number = new vector<size_t_vector>( NB_COMPS, work_size_t ) ;
   min_index_unknown_handled_by_proc = 
   	new vector<size_t_vector>( NB_COMPS, work_size_t ) ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
     for (size_t i=0;i<DIM;++i) 
       (*min_index_unknown_handled_by_proc)[comp](i) = 1000000;	
   max_index_unknown_handled_by_proc = 
   	new vector<size_t_vector>( NB_COMPS, work_size_t ) ;
   min_index_unknown_on_proc = 
   	new vector<size_t_vector>( NB_COMPS, work_size_t ) ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
     for (size_t i=0;i<DIM;++i) 
       (*min_index_unknown_on_proc)[comp](i) = 1000000;	
   max_index_unknown_on_proc = 
   	new vector<size_t_vector>( NB_COMPS, work_size_t ) ;	
	
   // Build global coordinates
   doubleVector work_double( 1, 0. ) ;
   vector< doubleVector > const* primary_mesh_global_main_coordinates = 
	PRIMARY_GRID->get_global_main_coordinates();
   size_t_vector const* primary_mesh_global_max_index = 
	PRIMARY_GRID->get_global_max_index() ;		
   size_t_vector np(DIM,0);
   for (size_t i=0;i<DIM;++i) 
     np(i) = (*primary_mesh_global_main_coordinates)[i].size(); 
   
   global_main_coordinates = new vector < vector< doubleVector > >;
   global_main_coordinates->reserve(NB_COMPS);
   vector< doubleVector > vwork_double( DIM, work_double ) ;
   for (size_t comp=0;comp<NB_COMPS;++comp) 
     global_main_coordinates->push_back( vwork_double ) ;   

   size_t dir_ref = DIM == 2 ? 2 : 0 ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {     
     for (size_t i=0;i<DIM;++i)
     {
       if ( i != dir_ref )
       {
         (*global_max_index)[comp](i) = (*primary_mesh_global_max_index)(i) ;
	 (*global_main_coordinates)[comp][i] = 
     	  	(*primary_mesh_global_main_coordinates)[i] ;
       }
       else
       {
         if ( (*PERIODIC)(i) ) 
         {
           (*global_max_index)[comp](i) = np(i) - 2 ;
	   (*global_main_coordinates)[comp][i].re_initialize(np(i)-1,0.);
           for( size_t j=0;j<np(i)-1;++j)
             (*global_main_coordinates)[comp][i](j) = 0.5 *
       		( (*primary_mesh_global_main_coordinates)[i](j)
		+ (*primary_mesh_global_main_coordinates)[i](j+1) ); 
         }
         else
         {
           (*global_max_index)[comp](i) = np(i) ;
	   (*global_main_coordinates)[comp][i].re_initialize(np(i)+1,0.);
           (*global_main_coordinates)[comp][i](0) = 
       		(*primary_mesh_global_main_coordinates)[i](0);
           for( size_t j=1;j<np(i);++j)
             (*global_main_coordinates)[comp][i](j) = 0.5 *
       		( (*primary_mesh_global_main_coordinates)[i](j)
		+ (*primary_mesh_global_main_coordinates)[i](j-1) );
           (*global_main_coordinates)[comp][i](np(i)) = 
           	(*primary_mesh_global_main_coordinates)[i](np(i)-1);   
         }
       }
     }
     dir_ref++; 
   }          

   // Build local coordinates
   vector< doubleVector > const* primary_mesh_local_main_coordinates = 
	PRIMARY_GRID->get_local_main_coordinates();
   local_main_coordinates = new vector < vector< doubleVector > >;
   local_main_coordinates->reserve(NB_COMPS); 
   size_t_vector const* primary_mesh_local_max_index_in_global = 
   	PRIMARY_GRID->get_local_max_index_in_global() ;
   size_t_vector const* primary_mesh_local_min_index_in_global = 
   	PRIMARY_GRID->get_local_min_index_in_global() ;
   for (size_t comp=0;comp<NB_COMPS;++comp) 
     local_main_coordinates->push_back( vwork_double ) ; 

   dir_ref = DIM == 2 ? 2 : 0 ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   { 
     for (size_t i=0;i<DIM;++i)
     {
       if ( i != dir_ref )
       {
         (*local_main_coordinates)[comp][i] = 
     		(*primary_mesh_local_main_coordinates)[i] ;
         (*local_min_index_in_global)[comp](i) = 
     		(*primary_mesh_local_min_index_in_global)(i) ;
         (*local_max_index_in_global)[comp](i) = 
     		(*primary_mesh_local_max_index_in_global)(i) ;	
         (*local_dof_number)[comp](i) = (*local_max_index_in_global)[comp](i) -
     		(*local_min_index_in_global)[comp](i) + 1 ; 
       }
       else
       {
         if ( (*PERIODIC)(i) ) 
         {
           (*local_min_index_in_global)[comp](i) =
       		(*primary_mesh_local_min_index_in_global)(i) ;
           (*local_max_index_in_global)[comp](i) = 
       	        (*primary_mesh_local_max_index_in_global)(i) - 1 ;
         }
         else
         {     
           (*local_min_index_in_global)[comp](i) = 
     		(*primary_mesh_local_min_index_in_global)(i) == 0 ?
     		0 : (*primary_mesh_local_min_index_in_global)(i)+1 ;
           (*local_max_index_in_global)[comp](i) = 
     		(*primary_mesh_local_max_index_in_global)(i) == np(i)-1 ?
     		np(i) : (*primary_mesh_local_max_index_in_global)(i) ;
         }         
         (*local_dof_number)[comp](i) = (*local_max_index_in_global)[comp](i) -
     		(*local_min_index_in_global)[comp](i) + 1 ; 
         (*local_main_coordinates)[comp][i].re_initialize( 
       		(*local_dof_number)[comp](i), 0.);
         for (size_t j=0;j<(*local_dof_number)[comp](i);++j)
           (*local_main_coordinates)[comp][i](j) = 
       		(*global_main_coordinates)[comp][i]( 
		(*local_min_index_in_global)[comp](i) 
		+ j );       
       }
     }
     dir_ref++;    
   }  

   // Nodes handled by current process
   vector< size_t_vector >* vwork_size_t = new vector< size_t_vector >( DIM, 
   	work_size_t ) ; 
   on_current_processor = new vector< vector< size_t_vector > >( NB_COMPS ,
   	*vwork_size_t );
   delete vwork_size_t;  		       
   dir_ref = DIM == 2 ? 2 : 0 ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   { 
     for (size_t i=0;i<DIM;++i)
     {
       if ( i != dir_ref )
       {
         (*on_current_processor)[comp][i].re_initialize( 
     		(*local_dof_number)[comp](i), 0 );
	
         if ( (*local_min_index_in_global)[comp](i) != 0 )
           for (size_t j=0;j<security_bandwidth+1;++j)
             (*on_current_processor)[comp][i]( j ) = 2 ;
         else if ( !(*primary_mesh_handles_node)[i](0) )
           for (size_t j=0;j<security_bandwidth+1;++j)
             (*on_current_processor)[comp][i]( j ) = 4 ;     
	 
         if ( (*local_max_index_in_global)[comp](i) 
	 	!= (*global_max_index)[comp](i) )
           for (size_t j=0;j<security_bandwidth;++j)
             (*on_current_processor)[comp][i]( 
	     	(*on_current_processor)[comp][i].size() - j-1 ) = 2 ;
         else if ( !(*primary_mesh_handles_node)[i](
		(*primary_mesh_handles_node)[i].size()-1) )
           for (size_t j=0;j<security_bandwidth;++j)
             (*on_current_processor)[comp][i]( 
	     	(*on_current_processor)[comp][i].size() - j-1 ) = 4 ;	
       }
       else
       {
         (*on_current_processor)[comp][i].re_initialize( 
     		(*local_dof_number)[comp](i), 0 );

         for (size_t j=0;j<security_bandwidth;++j)
         {
           if ( !(*primary_mesh_handles_node)[i](j) )
             if ( (*local_min_index_in_global)[comp](i) == 0 )
	       (*on_current_processor)[comp][i](j) = 4 ;
	     else
               (*on_current_processor)[comp][i](j) = 2 ;
	 
           if ( !(*primary_mesh_handles_node)[i](
		(*primary_mesh_handles_node)[i].size()-1-j) )
	     if ( (*local_max_index_in_global)[comp](i) 
	     	== (*global_max_index)[comp](i) )
               (*on_current_processor)[comp][i](
	       	(*local_dof_number)[comp](i)-1-j) = 4 ;
	     else
	       (*on_current_processor)[comp][i](
	       	(*local_dof_number)[comp](i)-1-j) = 2 ;  
         }
       }
     }
     dir_ref++;    
   }     
   
   // Bufferzone dof 
   dir_ref = DIM == 2 ? 2 : 0 ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   { 
     for (size_t i=0;i<DIM;++i)
     {
       if ( i != dir_ref )
       {
         if ( (*on_current_processor)[comp][i](0) == 2 && 
     	(*local_min_index_in_global)[comp](i) != 0 )
           for (size_t j=0;j<security_bandwidth;++j)
             (*on_current_processor)[comp][i](security_bandwidth+1+j) = 1 ;
     
         size_t index_max = (*on_current_processor)[comp][i].size()-1;
         if ( (*on_current_processor)[comp][i](index_max) == 2 && 
     	(*local_max_index_in_global)[comp](i) != (*global_max_index)[comp](i) )
           for (size_t j=0;j<security_bandwidth+1;++j)
             (*on_current_processor)[comp][i](index_max-security_bandwidth-j) 
	     	= 1 ;	
       }
       else
       {
         if ( (*on_current_processor)[comp][i](0) == 2 && 
     	(*local_min_index_in_global)[comp](i) != 0 )
           for (size_t j=0;j<security_bandwidth;++j)
             (*on_current_processor)[comp][i](security_bandwidth+j) = 1 ;
     
         size_t index_max = (*on_current_processor)[comp][i].size()-1;
         if ( (*on_current_processor)[comp][i](index_max) == 2 && 
     	(*local_max_index_in_global)[comp](i) != (*global_max_index)[comp](i) )
           for (size_t j=0;j<security_bandwidth;++j)
             (*on_current_processor)[comp][i](index_max-security_bandwidth-j) 
	   	= 1 ;
       }
     }
     dir_ref++;    
   } 
   
   // Periodic bufferzone
   dir_ref = DIM == 2 ? 2 : 0 ;   
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     for (size_t i=0;i<DIM;++i)
     {  
       if ( i != dir_ref ) 
       {  
         if ( (*local_min_index_in_global)[comp](i) == 0 && (*PERIODIC)(i) )
           for (size_t j=0;j<security_bandwidth;++j)
             if ( (*on_current_processor)[comp][i](security_bandwidth+j+1) == 1 )
               (*on_current_processor)[comp][i](security_bandwidth+j+1) = 5 ;
	     else
	       (*on_current_processor)[comp][i](security_bandwidth+j+1) = 3 ;
     
         size_t index_max = (*on_current_processor)[comp][i].size()-1;
         if ( (*local_max_index_in_global)[comp](i) 
	 	== (*global_max_index)[comp](i) && (*PERIODIC)(i) )
           for (size_t j=0;j<security_bandwidth+1;++j)
             if ( (*on_current_processor)[comp][i]
	     	(index_max-security_bandwidth-j) == 1 )
               (*on_current_processor)[comp][i](index_max-security_bandwidth-j)
	       		= 5 ;
	     else
	       (*on_current_processor)[comp][i](index_max-security_bandwidth-j)
	       		= 3 ;     
       }
       else
       {   
         if ( (*local_min_index_in_global)[comp](i) == 0 && (*PERIODIC)(i) )
           for (size_t j=0;j<security_bandwidth;++j)
             if ( (*on_current_processor)[comp][i](security_bandwidth+j) == 1 )
               (*on_current_processor)[comp][i](security_bandwidth+j) = 5 ;
	     else
	       (*on_current_processor)[comp][i](security_bandwidth+j) = 3 ;
     
         size_t index_max = (*on_current_processor)[comp][i].size()-1;
         if ( (*local_max_index_in_global)[comp](i) 
	 	== (*global_max_index)[comp](i) && (*PERIODIC)(i) )
           for (size_t j=0;j<security_bandwidth;++j)
             if ( (*on_current_processor)[comp][i]
	     	(index_max-security_bandwidth-j) == 1 )
               (*on_current_processor)[comp][i](index_max-security_bandwidth-j)
	       		= 5 ;
	     else
	       (*on_current_processor)[comp][i](index_max-security_bandwidth-j)
	       		= 3 ;
       }
     }
     dir_ref++;       
   }      

   // Allocate field values and colors
   doubleArray3D* work_doubleArray3D = new doubleArray3D(1,1,1,0.); 
   vector< doubleArray3D >* work_values = new vector< doubleArray3D >
   	( NB_COMPS , *work_doubleArray3D );   
   delete work_doubleArray3D;   	
   VALUES = new vector< vector< doubleArray3D > >( STO_DEPTH, *work_values );
   work_values->clear();
   delete work_values;
   intArray3D* work_intArray3D = new intArray3D(1,1,1,0); 
   DOFcolors = new vector< intArray3D >( NB_COMPS, *work_intArray3D );
   DOFstatus = new vector< intArray3D >( NB_COMPS, *work_intArray3D );   
   delete work_intArray3D; 
   for (size_t comp=0;comp<NB_COMPS;++comp) 
   {
     (*DOFcolors)[comp].re_initialize( 
     	(*local_main_coordinates)[comp][0].size(),
   	(*local_main_coordinates)[comp][1].size(),
	DIM == 2 ? 1 : (*local_main_coordinates)[comp][2].size() );
     (*DOFstatus)[comp].re_initialize( 
     	(*local_main_coordinates)[comp][0].size(),
   	(*local_main_coordinates)[comp][1].size(),
	DIM == 2 ? 1 : (*local_main_coordinates)[comp][2].size() );
     (*DOFstatus)[comp].set( FV_DOF_ONPROC );	      
     for (size_t lev=0;lev<STO_DEPTH;++lev)   
       (*VALUES)[lev][comp].re_initialize( 
     	(*local_main_coordinates)[comp][0].size(),
   	(*local_main_coordinates)[comp][1].size(),
	DIM == 2 ? 1 : (*local_main_coordinates)[comp][2].size() ); 
   } 
   
   // Set dof colors
   for (size_t comp=0;comp<NB_COMPS;++comp) set_DOF_colors(comp);

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
       PERIODIC_SHIFT = new vector< FV_SHIFT_TRIPLET >( NB_COMPS, mst ) ;
     }
   dir_ref = DIM == 2 ? 2 : 0 ;      
   for (size_t comp=0;comp<NB_COMPS;++comp) 
   {
     for (size_t i=0;i<DIM;++i) 
       if ( (*PERIODIC)(i) ) 
       {
         if ( i != dir_ref )
	 {
	   (*periodic_depth)(i) = security_bandwidth + 1 ;
           (*PERIODIC_SHIFT)[comp].index(i) = 
	   	(*global_main_coordinates)[comp][i].size()
       		- 2 * security_bandwidth - 1 ;	 
	 }
	 else
	 {
           (*periodic_depth)(i) = security_bandwidth ;
           (*PERIODIC_SHIFT)[comp].index(i) = 
	   	(*global_main_coordinates)[comp][i].size()
       		- 2 * security_bandwidth ;
	 }
       }     
     set_PERIODIC_default( comp, periodic_depth );
     dir_ref++;         
   }
   if ( periodic_depth ) delete periodic_depth ; 

   // Set dof status
   for (size_t comp=0;comp<NB_COMPS;++comp) set_DOF_status(comp); 

   // Compute cell size
   local_cell_size = new vector < vector< doubleVector > >;
   local_cell_size->reserve( NB_COMPS ); 
   for (size_t comp=0;comp<NB_COMPS;++comp) 
     local_cell_size->push_back( vwork_double ) ;
   vector< doubleVector > const* primary_coordinates = PRIMARY_GRID->
   	get_local_main_coordinates() ; 
   size_t shift = 0, start_index = 0, end_index = 0 ; 
  
   dir_ref = DIM == 2 ? 2 : 0 ;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   { 
     for (size_t i=0;i<DIM;++i)
     {
       if ( i != dir_ref )
       {
         (*local_cell_size)[comp][i].re_initialize( 
	 	(*local_dof_number)[comp](i), 0. );
         if ( (*local_min_index_in_global)[comp](i) == 0 )
         {
           start_index = 1 ;       
           if ( (*PERIODIC)(i) )
             (*local_cell_size)[comp][i](0) = 
	     	(*local_main_coordinates)[comp][i](1)
       		- (*local_main_coordinates)[comp][i](0) ; 
           else     
             (*local_cell_size)[comp][i](0) = 0.5 * ( 
	 	(*local_main_coordinates)[comp][i](1)
       		- (*local_main_coordinates)[comp][i](0) ) ;
         }
         else start_index = 0 ;
     
         if ( (*local_max_index_in_global)[comp](i) == 
	 	(*global_max_index)[comp](i) )
         {
           end_index = (*local_dof_number)[comp](i)-2 ;
           if ( (*PERIODIC)(i) )
             (*local_cell_size)[comp][i](end_index+1) = 
	 	(*local_main_coordinates)[comp][i](end_index+1)
       		- (*local_main_coordinates)[comp][i](end_index) ; 
           else       
             (*local_cell_size)[comp][i](end_index+1) = 0.5 * (
       		(*local_main_coordinates)[comp][i](end_index+1)
       		- (*local_main_coordinates)[comp][i](end_index) ) ;    
         }
         else end_index = (*local_dof_number)[comp](i)-1 ;       
     
         for (size_t j=start_index;j<=end_index;++j)
           (*local_cell_size)[comp][i](j) = 0.5 * ( 
       	    	(*primary_mesh_global_main_coordinates)[i]( 
		(*local_min_index_in_global)[comp](i) + j + 1 )
     		- (*primary_mesh_global_main_coordinates)[i](
		(*local_min_index_in_global)[comp](i) + j - 1 ) ) ;	
       }
       else
       {
         (*local_cell_size)[comp][i].re_initialize( 
       	  	(*local_dof_number)[comp](i), 0. );
         if ( (*local_min_index_in_global)[comp](i) == 0 
     		&& !(*PERIODIC)(i) )
         {
           shift = 1 ;
           start_index = 1 ;       
           (*local_cell_size)[comp][i](0) = 
	   	(*local_main_coordinates)[comp][i](1)
       		- (*local_main_coordinates)[comp][i](0) ;
         }
         else
         { 
           shift = 0 ;
           start_index = 0 ;
         }
     
         if ( (*local_max_index_in_global)[comp](i) == 
       		(*global_max_index)[comp](i) && !(*PERIODIC)(i) )
         {
           end_index = (*local_dof_number)[comp](i)-2 ;
           (*local_cell_size)[comp][i](end_index+1) = 
       		(*local_main_coordinates)[comp][i](end_index+1)
       		- (*local_main_coordinates)[comp][i](end_index) ;     
         }
         else end_index = (*local_dof_number)[comp](i)-1 ;       
     
         for (size_t j=start_index;j<=end_index;++j)
           (*local_cell_size)[comp][i](j) = 
       		(*primary_coordinates)[i](j+1-shift)
     		- (*primary_coordinates)[i](j-shift) ;
       }
     }
     dir_ref++;    
   }                  
}




//----------------------------------------------------------------------
FV_DiscreteField_Vorticity:: ~FV_DiscreteField_Vorticity( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: ~FV_DiscreteField_Vorticity") ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: print( std::ostream& os, size_t indent_width,
      	bool b_values ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: print" ) ;

   std::string space( indent_width, ' ' ) ;   
   
   os << space << "FV_DiscreteField_Staggered n°" << ID << endl;
   os << space << "   Name = " << FNAME << endl;   
   os << space << "   Discretization type = " << FDISCRETIZATION << endl;
   os << space << "   Number of components = " << NB_COMPS << endl;
   os << space << "   Level of storage = " << STO_DEPTH << endl;
   os << space << "Field mesh (on all processors)" << endl;
   os << space << "-------------------------------" << endl;
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     os << space << "Component " << comp << endl;
     os << space << "----------- " << endl;
     for (size_t i=0;i<DIM;++i)
     {
       os << space << "Global number of points in direction " << 
     	(*FV_Mesh::directionName)(i) <<
     	" = " << (*global_main_coordinates)[comp][i].size() << endl;
       os << space << "Min index = " << 0 << "  Max index = " 
     	<< (*global_max_index)[comp](i) << endl;
       os << space << "Coordinate values = ";
       for (size_t j=0;j<(*global_main_coordinates)[comp][i].size();++j)
         os << space << " " << (*global_main_coordinates)[comp][i](j);
       os << endl;
     }
   } 
   os << space << "Local field on processor " << 
   	MAC_Exec::communicator()->rank() << endl;
   os << space << "---------------------------------" << endl;
   os << space << "Mesh" << endl;
   os << space << "----" << endl;   
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {     
     os << space << "Component " << comp << endl;
     os << space << "----------- " << endl;
     for (size_t i=0;i<DIM;++i)
     {
       os << space << "Local number of points in direction " << 
     	(*FV_Mesh::directionName)(i) <<
     	" = " << (*local_main_coordinates)[comp][i].size() << endl;
       os << space << "Min index = " << (*local_min_index_in_global)[comp](i) 
       	<< "  Max index = " << (*local_max_index_in_global)[comp](i) << endl;	
       os << space << "Coordinate values = ";
       for (size_t j=0;j<(*local_main_coordinates)[comp][i].size();++j)
         os << space << " " << (*local_main_coordinates)[comp][i](j);
       os << endl;
       os << space << "Status = ";
       for (size_t j=0;j<(*on_current_processor)[comp][i].size();++j)
         os << space << " " << (*on_current_processor)[comp][i](j);
       os << endl;
       os << space << "Coordinates in halo zone = ";
       for (size_t j=0;j<(*on_current_processor)[comp][i].size();++j)
         if ( (*on_current_processor)[comp][i](j) == 2 
       	|| (*on_current_processor)[comp][i](j) == 4 ) 
           os << space << " " << (*local_main_coordinates)[comp][i](j);
       os << endl;
       os << space << "Coordinates in buffer zone = ";
       for (size_t j=0;j<(*on_current_processor)[comp][i].size();++j)
         if ( (*on_current_processor)[comp][i](j) == 1 
       	|| (*on_current_processor)[comp][i](j) == 3 
       	|| (*on_current_processor)[comp][i](j) == 5 ) 
           os << space << " " << (*local_main_coordinates)[comp][i](j);	 
       os << endl;
       if ( PERIODIC_SHIFT )
         os << "Periodic shift = " << (*PERIODIC_SHIFT)[comp].get_index(i) 
	 	<< endl;   
       os << space << "Cell size = ";
       for (size_t j=0;j<(*local_cell_size)[comp][i].size();++j)
         os << space << " " << (*local_cell_size)[comp][i](j);
       os << endl;       	     	 
     }
   }

   print_BC_FieldValues( os, indent_width, b_values ) ;      
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: build_BCs( MAC_ModuleExplorer const* exp, 
      	FV_DomainBuilder const* DB ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: build_BCs" ) ;

   // Add triplets to BCs
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     size_t kmax = DIM == 2 ? 1 : (*local_dof_number)[comp](2);   
     for (size_t i=0;i<(*local_dof_number)[comp](0);++i)
       for (size_t j=0;j<(*local_dof_number)[comp](1);++j) 
         for (size_t k=0;k<kmax;++k) 
	   if ( (*DOFcolors)[comp](i,j,k) > 1 )
	     (*V_BCS)[(*DOFcolors)[comp](i,j,k)]->add_MacTriplet( comp, i, 
	     	j, k );
   }
   
   // Set unknown on BC
   for (size_t comp=0;comp<NB_COMPS;++comp) 
     for (list< FV_BoundaryCondition* >::iterator il=SET_OF_BCS.begin();
     		il!=SET_OF_BCS.end();il++)
       if ( (*il)->is_neumann( comp ) ) (*il)->set_unknown_on_BC( comp ) ;
   
   // Set free_DOF_on_boundary update features 
   set_freeDOFonBC_update_features() ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Vorticity:: get_DOF_coordinate( size_t i, size_t component, 
      	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: get_DOF_coordinate" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( direction < DIM ) ;     
   MAC_CHECK_PRE( i < (*local_dof_number)[component](direction) );

   return( (*local_main_coordinates)[component][direction](i) ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Vorticity:: get_DOF_coordinate_Assembling( int i, 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: get_DOF_coordinate_Assembling" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( direction < DIM ) ;   

   double result = i < 0 ? 0. : 
   	i >= int((*local_dof_number)[component](direction)) ?
   	0. : (*local_main_coordinates)[component][direction](i) ;
   
   return( result );
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Vorticity:: is_global_triplet( 
	int i, int j, int k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: is_global_triplet" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );  
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   

   bool result = 
   	( i >= 0 && i <= int((*global_max_index)[component](0))
   	&& j >= 0 && j <= int((*global_max_index)[component](1)) );
   if ( DIM == 3 && result )
     result = ( k >= 0 && k <= int((*global_max_index)[component](2)) );
   
   return ( result ) ;  
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Vorticity:: is_global_triplet_local_DOF( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: is_global_triplet_local_DOF" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );  
   MAC_CHECK_PRE( i <= (*global_max_index)[component](0) );
   MAC_CHECK_PRE( j <= (*global_max_index)[component](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	k <= (*global_max_index)[component](2) ) ) ;

   bool result = 
   	( i >= (*local_min_index_in_global)[component](0) 
	&& i <= (*local_max_index_in_global)[component](0)
   	&& j >= (*local_min_index_in_global)[component](1)
	&& j <= (*local_max_index_in_global)[component](1) );
   if ( DIM == 3 && result )
     result = ( k >= (*local_min_index_in_global)[component](2) 
     	&& k <= (*local_max_index_in_global)[component](2) );
   
   return ( result ) ;  
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: set_postprocessing_options( 
		std::string const& a_location,
		std::string const& a_paraview_fname ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: set_post_processing_options" ) ;
   MAC_CHECK_PRE( a_location == "at_cell_centers"
   		|| a_location == "at_vertices" ) ;   
   
   LOCATION = a_location ;
   PARAVIEW_FNAME = a_paraview_fname ;

   if ( LOCATION != "at_vertices" )
   {
     ostringstream mesg ;
     mesg << "*** Unable to save field of name \""
        << FNAME << "\" at location \"" << LOCATION << "\"." << endl << endl ;
     mesg << "    allowed location is \"at_vertices\" " << endl;
     MAC_Error::object()->raise_plain( mesg.str() ) ;
   }

   // Set the number of dof written in post-processing output files
   set_nb_dof_post_processing() ;

}      




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: write_field(MAC_Module* point_data,
                               MAC_Module* cell_data) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: write_field" ) ;

   size_t iter = 0;
   bool scalar = ( NB_COMPS == 1 ) ;
   MAC_Module* target = point_data ;

   // Surement mieux à faire, mais je ne sais pas faire simple avec un pointeur
   // pour éviter la recopie (GVi)
   doubleArray2D X( NB_COMPS, NB_DOF_POSTPROCESSING_PER_COMP );
   vector< boolVector > const* primary_mesh_handles_node = 
   	PRIMARY_GRID->get_nodes_owner();

   size_t nelem0 = PRIMARY_GRID->get_local_nb_points(0) ;
   size_t nelem1 = PRIMARY_GRID->get_local_nb_points(1) ;
   size_t nelem2 =  DIM == 2 ?
   		1 : PRIMARY_GRID->get_local_nb_points(2) ;

   size_t shift_i = (*local_min_index_in_global)[0](0) == 0 ? 1 : 0;
   size_t shift_j = (*local_min_index_in_global)[0](1) == 0 ? 1 : 0;
   size_t shift_k = DIM == 2 ?
   		0 : (*local_min_index_in_global)[0](2) == 0 ? 1 : 0;

   if ( DIM == 2 )
     for (size_t comp=0;comp<NB_COMPS;++comp)
     {
       iter = 0 ;
       for (size_t k=0;k<nelem2;++k)
         for (size_t j=0;j<nelem1;++j) 
           for (size_t i=0;i<nelem0;++i)
	   {
	     if ( (*primary_mesh_handles_node)[0](i)
	     		&& (*primary_mesh_handles_node)[1](j) )
	     {
	       X( comp, iter ) = interpolate_values_at_nodes(
			i, shift_i, j, shift_j, k, shift_k, comp, 0 ) ;
	       iter++;
	     }
	   }
     }      
   else
     for (size_t comp=0;comp<NB_COMPS;++comp)
     {
       iter = 0 ;
       for (size_t k=0;k<nelem2;++k)
         for (size_t j=0;j<nelem1;++j) 
           for (size_t i=0;i<nelem0;++i)
	   {
	     if ( (*primary_mesh_handles_node)[0](i)
	     		&& (*primary_mesh_handles_node)[1](j)
	     		&& (*primary_mesh_handles_node)[2](k) )
	     {
	       X( comp, iter ) = interpolate_values_at_nodes( 
			i, shift_i, j, shift_j, k, shift_k, comp, 0 ) ;
	       iter++;
	     }
	   }
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
FV_DiscreteField_Vorticity:: set_nb_dof_post_processing( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: set_nb_dof_post_processing") ;

   NB_DOF_POSTPROCESSING_PER_COMP = 1 ;
   for (size_t i=0;i<DIM;++i)
     NB_DOF_POSTPROCESSING_PER_COMP
       	*= PRIMARY_GRID->get_local_nb_points_on_current_proc( i ) ;     
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Vorticity:: interpolate_values_at_nodes( 
      		size_t i, size_t shift_i,
		size_t j, size_t shift_j,
		size_t k, size_t shift_k,
		size_t component,
               	size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: interpolate_values_at_nodes") ;
   MAC_CHECK_PRE( LOCATION == "at_vertices" )

   double result = 0.;
   double coeff_bottom, coeff_top, coeff_left, coeff_right,
   	coeff_behind, coeff_front;

   if ( DIM == 2 ) result = (*VALUES)[level][component](i,j,k) ;
   else
   {
     if ( component == 0 )
     {
       // Left plane: left, bottom_left, top_left, behind_left, front_left, 
       // behind_bottom_left, behind_top_left, front_bottom_left, 
       // front_top_left
       if ( (*DOFcolors)[component](i,j,k) == FV_BC_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_TOP_LEFT 
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_LEFT 		
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP_LEFT )
         result = (*VALUES)[level][component](i,j,k) ;
       // Right plane: right, bottom_right, top_right, behind_right, 
       // front_right, behind_bottom_right, behind_top_right, 
       // front_bottom_right, front_top_right
       if ( (*DOFcolors)[component](i,j,k) == FV_BC_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BOTTOM_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_TOP_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP_RIGHT )
         result = (*VALUES)[level][component](i+shift_i,j,k) ;
       // interior, bottom, top, behind, front, behind_bottom, behind_top, 
       // front_bottom, front_top
       else
       {
         coeff_right = (*local_cell_size)[component][0](i+shift_i);
         coeff_left = (*local_cell_size)[component][0](i+shift_i-1);
         result = ( coeff_right
       			* (*VALUES)[level][component](i+shift_i,j,k)
       		+ coeff_left
			* (*VALUES)[level][component](i+shift_i-1,j,k) )
		 / ( coeff_right + coeff_left ) ;
       }
     }
     else if ( component == 1 )
     {
       // Bottom plane: bottom, bottom_left, bottom_right, behind_bottom,
       // front_bottom, behind_bottom_left, behind_bottom_right,
       // front_bottom_left, front_bottom_right
       if ( (*DOFcolors)[component](i,j,k) == FV_BC_BOTTOM
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BOTTOM_RIGHT 
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM_RIGHT )
         result = (*VALUES)[level][component](i,j,k) ;
       // Top plane: top, top_left, top_right, behind_top,
       // front_top, behind_top_left, behind_top_right,
       // front_top_left, front_top_right
       if ( (*DOFcolors)[component](i,j,k) == FV_BC_TOP
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_TOP_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_TOP_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP 		
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP_RIGHT )
         result = (*VALUES)[level][component](i,j+shift_j,k) ;
       // interior, left, right, behind, front, behind_left, behind_right,
       // front_left, front_right
       else
       {
         coeff_top = (*local_cell_size)[component][1](j+shift_j);
         coeff_bottom = (*local_cell_size)[component][1](j+shift_j-1);
         result = ( coeff_top
       			* (*VALUES)[level][component](i,j+shift_j,k)
       		+ coeff_bottom
			* (*VALUES)[level][component](i,j+shift_j-1,k) )
		 / ( coeff_top + coeff_bottom ) ;
       }
     }
     else if ( component == 2 )
     {
       // Behind plane: behind, behind_left, behind_right, behind_bottom,
       // behind_top, behind_bottom_left, behind_bottom_right, behind_top_left
       // behind_top_right
       if ( (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_RIGHT 
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP 		
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_BOTTOM_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_BEHIND_TOP_RIGHT )
         result = (*VALUES)[level][component](i,j,k) ;
       // Front plane: front, front_left, front_right, front_bottom,
       // front_top, front_bottom_left, front_bottom_right, front_top_left
       // front_top_right
       if ( (*DOFcolors)[component](i,j,k) == FV_BC_FRONT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP		
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_BOTTOM_RIGHT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP_LEFT
       		|| (*DOFcolors)[component](i,j,k) == FV_BC_FRONT_TOP_RIGHT )
         result = (*VALUES)[level][component](i,j,k+shift_k) ;
       // interior, left, right, bottom, top, bottom_left, bottom_right
       // top_left, top_right  
       else
       {
         coeff_front = (*local_cell_size)[component][2](k+shift_k);
         coeff_behind = (*local_cell_size)[component][2](k+shift_k-1);
         result = ( coeff_front
       			* (*VALUES)[level][component](i,j,k+shift_k)
       		+ coeff_behind
			* (*VALUES)[level][component](i,j,k+shift_k-1) )
		 / ( coeff_front + coeff_behind ) ;
       }
     }
   }   
      
   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Vorticity:: get_cell_size( size_t i, size_t component, 
      	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: get_cell_size") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS );       
   MAC_CHECK_PRE( check_invariant_cell_features( i, component, direction ) ) ;

   double result = (*local_cell_size)[component][direction](i) ;
	         
   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Vorticity:: get_cell_measure( size_t i, size_t j, size_t k,
      	size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: get_cell_measure") ;   
   MAC_CHECK_PRE( component < NB_COMPS );    
   MAC_CHECK_PRE( check_invariant_cell_features( i, component, 0 ) ) ;
   MAC_CHECK_PRE( check_invariant_cell_features( j, component, 1 ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	check_invariant_cell_features( k, component, 2 ) ) ) ;
   
   double result = (*local_cell_size)[component][0](i) 
   	* (*local_cell_size)[component][1](j) ;
   if ( DIM == 3 ) 
     result *= (*local_cell_size)[component][2](k) ;
	
   return ( result ) ;
}




//----------------------------------------------------------------------
double
FV_DiscreteField_Vorticity:: get_face_perp_to_direction_measure( 
	size_t i, size_t j, size_t k, size_t component, size_t direction ) 
	const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Vorticity:: get_face_perp_to_direction_measure") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS );    
   MAC_CHECK_PRE( check_invariant_cell_features( i, component, 0 ) ) ;
   MAC_CHECK_PRE( check_invariant_cell_features( j, component, 1 ) ) ;
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	check_invariant_cell_features( k, component, 2 ) ) ) ;
	
   double result = 0. ;
 
   if ( DIM == 2 )
   {
     switch( direction )
     {
       case 0: 
         result = (*local_cell_size)[component][1](j) ;
         break;
       case 1: 
         result = (*local_cell_size)[component][0](i) ;
         break;       
     }
   }
   else
   {
     switch( direction )
     {
       case 0: 
         result = (*local_cell_size)[component][1](j) 
	 	* (*local_cell_size)[component][2](k) ;
         break;
       case 1: 
         result = (*local_cell_size)[component][0](i) 
	 	* (*local_cell_size)[component][2](k) ;
         break;       
       case 2: 
         result = (*local_cell_size)[component][0](i) 
	 	* (*local_cell_size)[component][1](j) ; 
         break;
     }   
   }
                  
   return ( result ) ;	
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Vorticity:: get_min_index_unknown_handled_by_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Vorticity:: get_min_index_unknown_handled_by_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*min_index_unknown_handled_by_proc)[component](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Vorticity:: get_max_index_unknown_handled_by_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Vorticity:: get_max_index_unknown_handled_by_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*max_index_unknown_handled_by_proc)[component](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Vorticity:: get_min_index_unknown_on_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Vorticity:: get_min_index_unknown_on_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*min_index_unknown_on_proc)[component](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Vorticity:: get_max_index_unknown_on_proc( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Vorticity:: get_max_index_unknown_on_proc") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*max_index_unknown_on_proc)[component](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_DiscreteField_Vorticity:: get_local_nb_dof( 
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_DiscreteField_Vorticity:: get_local_nb_dof") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS ); 
     
   size_t result = (*local_dof_number)[component](direction) ;
   
   return ( result ) ;
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Vorticity:: DOF_has_imposed_Dirichlet_value( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: DOF_has_imposed_Dirichlet_value" );
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( i <= (*local_dof_number)[component](0) );
   MAC_CHECK_PRE( j <= (*local_dof_number)[component](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	k <= (*local_dof_number)[component](2) ) ) ;      

   bool result = false ;
   if ( (*DOFcolors)[component](i,j,k) > 1 )
     result = 
       get_BC( (*DOFcolors)[component](i,j,k) )->is_dirichlet( component ) ;
   
   return ( result ) ; 
}	




//----------------------------------------------------------------------
bool
FV_DiscreteField_Vorticity:: DOF_on_BC( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: DOF_on_BC" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 
   MAC_CHECK_PRE( i <= (*local_dof_number)[component](0) );
   MAC_CHECK_PRE( j <= (*local_dof_number)[component](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	k <= (*local_dof_number)[component](2) ) ) ;      

   bool result = (*DOFcolors)[component](i,j,k) > 1 ;
   
   return ( result ) ; 
}	




//----------------------------------------------------------------------
bool
FV_DiscreteField_Vorticity:: DOF_is_unknown_handled_by_proc( 
	size_t i, size_t j, size_t k, size_t component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: DOF_is_unknown_handled_by_proc" ) ;
   MAC_CHECK_PRE( component < NB_COMPS );
   MAC_CHECK_PRE( i <= (*local_dof_number)[component](0) );
   MAC_CHECK_PRE( j <= (*local_dof_number)[component](1) );
   MAC_CHECK_PRE( IMPLIES( DIM == 2, k == 0 ) ) ;   
   MAC_CHECK_PRE( IMPLIES( DIM == 3, 
   	k <= (*local_dof_number)[component](2) ) ) ;

   bool result = (*UNK_LOCAL_NUMBERING)[component]( i, j, k ) != -1 
   	&& (*DOFstatus)[component](i,j,k) != FV_DOF_HALOZONE
	&& (*DOFstatus)[component](i,j,k) != FV_DOF_PERIODIC_HALOZONE ;
   
   return ( result ) ; 
}




//----------------------------------------------------------------------
bool
FV_DiscreteField_Vorticity:: check_invariant_cell_features( size_t i,
	size_t component, size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Staggered:: check_invariant_cell_features") ;
   MAC_ASSERT( i < (*local_dof_number)[component](direction) ) ; 
	
   return ( true ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: compute_normLinf( double time ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: compute_normLinf" ) ;

}




//----------------------------------------------------------------------
doubleVector const*
FV_DiscreteField_Vorticity:: get_DOF_coordinates_vector( size_t component, 
      	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: get_DOF_coordinates_vector") ;
   MAC_CHECK_PRE( direction < DIM ) ; 
   MAC_CHECK_PRE( component < NB_COMPS );
   MAC_CHECK_PRE( local_main_coordinates != NULL ) ;    
   
   doubleVector const* result = 
   	&((*local_main_coordinates)[component][direction]) ;

   return ( result ) ;      
}




//----------------------------------------------------------------------
FV_SHIFT_TRIPLET
FV_DiscreteField_Vorticity:: shift_vorticityToStaggered( size_t component ) 
	const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: shift_vorticityToStaggered" ) ;

   MAC_CHECK_PRE( component < NB_COMPS );
   size_t i ;
   
   // If (i,j,k) is the location of the vorticity dof, 
   // (i+result.i,j+result.j,k+result.k) is the location of the staggered
   // dof on the right, top, front of the vorticity dof
   
   FV_SHIFT_TRIPLET result ;
   result.i = 0 ;
   result.j = 0 ;   
   result.k = 0 ;
   
   if ( DIM == 2 )
   {
     for (i=0;i<DIM;++i)
       if ( (*local_min_index_in_global)[0](i) == 0 && !(*PERIODIC)(i) ) 
         result.index( i ) = 1 ;
   }
   else
   {
     for (i=0;i<DIM;++i)
       if ( (*local_min_index_in_global)[component](i) == 0 && i != component
       		&& !(*PERIODIC)(i) )
         result.index( i ) = 1 ; 
   }        
     
   return ( result ) ;
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: create_transproj_interpolation( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: create_transproj_interpolation" ) ;

   ostringstream mesg ;
   mesg << "Field " << FNAME << " is of the FV_DiscreteField_Vorticity "
   	<<"type; method \"create_transproj_interpolation\" cannot be called by "
	<< "a field of type FV_DiscreteField_Vorticity !!" << endl
	<< " Check your implementation !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ; 
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: translation_projection( size_t const& level, 
      	size_t const& temporary_level, bool translate_mesh,
	doubleVector const* outOfDomain_values )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: translation_projection" ) ;

   if ( MAC_Exec::communicator()->rank() == 0 )
     FV::out() << "Field " << FNAME << " is of the FV_DiscreteField_Vorticity"
   	<< " type; method \"translation_projection\" translates field mesh"
	<< " only" << endl;

   // Translate field mesh
   size_t trans_dir = PRIMARY_GRID->get_translation_direction() ;
   double trans_mag = PRIMARY_GRID->get_translation_magnitude() ;
   if ( translate_mesh )
   {
     size_t ncoor = 0 ;      
     for (size_t comp=0;comp<NB_COMPS;++comp)
     {
       ncoor = (*global_main_coordinates)[comp][trans_dir].size() ;
       for (size_t i=0;i<ncoor;++i)
         (*global_main_coordinates)[comp][trans_dir](i) += trans_mag ;
       ncoor = (*local_main_coordinates)[comp][trans_dir].size() ;
       for (size_t i=0;i<ncoor;++i)
         (*local_main_coordinates)[comp][trans_dir](i) += trans_mag ;
     }
   }   
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: restore_translated_field_mesh( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: restore_translated_field_mesh" ) ;

   size_t trans_dir = PRIMARY_GRID->get_translation_direction() ;
   double trans_dist = PRIMARY_GRID->get_translation_distance() ;
   size_t ncoor = 0 ;

   if ( fabs( trans_dist ) > 1e-12 )
   {       
     if ( MAC_Exec::communicator()->rank() == 0 ) 
        FV::out() << "      Restore translated mesh of field " << FNAME
		<< endl ;
     for (size_t comp=0;comp<NB_COMPS;++comp)
     {
       ncoor = (*global_main_coordinates)[comp][trans_dir].size() ;
       for (size_t i=0;i<ncoor;++i)
         (*global_main_coordinates)[comp][trans_dir](i) += trans_dist ;
       ncoor = (*local_main_coordinates)[comp][trans_dir].size() ;
       for (size_t i=0;i<ncoor;++i)
         (*local_main_coordinates)[comp][trans_dir](i) += trans_dist ;
     } 
   }       
}




//----------------------------------------------------------------------
void
FV_DiscreteField_Vorticity:: translate_field_mesh( const size_t& trans_dir, 
      	const double &trans_dist )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: translate_field_mesh" ) ;

   size_t ncoor = 0 ;
   
   for (size_t comp=0;comp<NB_COMPS;++comp)
   {
     ncoor = (*global_main_coordinates)[comp][trans_dir].size() ;
     for (size_t i=0;i<ncoor;++i)
       (*global_main_coordinates)[comp][trans_dir](i) += trans_dist ;
     ncoor = (*local_main_coordinates)[comp][trans_dir].size() ;
     for (size_t i=0;i<ncoor;++i)
       (*local_main_coordinates)[comp][trans_dir](i) += trans_dist ;
   }
} 




//---------------------------------------------------------------------- 

bool FV_DiscreteField_Vorticity::DOF_offset( int &i, int &j, int &k,
        size_t_vector center, size_t_vector stencil,
	vector<double> &offset, size_t component ) const 
	 
//----------------------------------------------------------------------	 
{
   MAC_LABEL( "FV_DiscreteField_Vorticity:: DOF_offset" ) ;
   MAC_CHECK_PRE( component < NB_COMPS ); 

   bool outdomain_or_BC = false ;
   
   if ( i + int((*local_min_index_in_global)[component](0)) < 0 )
     {
      if ( component == 0 )  i += - 2 * i + 1 ;
      else i += -2 * i ;
      
      offset[0] = double( 2 * i - 1 ) * (*local_cell_size)[component][0](i);
      outdomain_or_BC = true;    
     }
   else if (i + int((*local_min_index_in_global)[component](0)) == 0 )
     {
      if ( component == 0 )     
      {
       i += 1;
       offset[0] = (*local_cell_size)[component][0](i);
      }
      else
      {
       i = center(0) - stencil(0) - 1;
       i += -2 * i;
       offset[0] = double( 2 * i - 1 ) * (*local_cell_size)[component][0](i); 
      } 
      outdomain_or_BC = true; 
     } 
   else if ( i + int((*local_min_index_in_global)[component](0)) ==
   	int((*global_max_index)[component](0)) ) 
     {
      if ( component == 0 )
      {
       i -= 1;
       offset[0] = -(*local_cell_size)[component][0](i);
      }
      else
      {
       i = center(0) + stencil(0) + 1;
       i +=  - 2 * ( i + int((*local_min_index_in_global)[component](0)) 
      		- int((*global_max_index)[component](0)) );
       offset[0] = ( 1 + 2 * ( i + int((*local_min_index_in_global)[component]
     	(0)) - int((*global_max_index)[component](0)) ) ) * 
     	(*local_cell_size)[component][0](i);		
      }
      outdomain_or_BC = true; 
     }	
   else if ( i + int((*local_min_index_in_global)[component](0)) >
   	int((*global_max_index)[component](0)) )
     {
      if ( component == 0 ) 
       i +=  - 2 * ( i + int((*local_min_index_in_global)[component](0)) 
      		- int((*global_max_index)[component](0)) ) - 1;
      else 
       i +=  - 2 * ( i + int((*local_min_index_in_global)[component](0)) 
      		- int((*global_max_index)[component](0)) );
		
      offset[0] = ( 1 + 2 * ( i + int((*local_min_index_in_global)[component]
     	(0)) - int((*global_max_index)[component](0)) ) ) * 
     	(*local_cell_size)[component][0](i); 
	
      outdomain_or_BC = true;
     }
       
   if ( j + int((*local_min_index_in_global)[component](1)) < 0 )
     {
      if ( component == 1 )  j += - 2 * j + 1 ;
      else j += -2 * j ;
      
      offset[1] = double( 2 * j - 1 ) * (*local_cell_size)[component][1](j);    
      outdomain_or_BC = true;    
     }
   else if ( j + int((*local_min_index_in_global)[component](1)) == 0 )
     {
      if ( component == 1 )     
      {
       j += 1;
       offset[1] = (*local_cell_size)[component][1](j);
      }
      else
      {
       j = center(1) - stencil(1) - 1;
       j += -2 * j;
       offset[1] = double( 2 * j - 1 ) * (*local_cell_size)[component][1](j); 
      } 
      outdomain_or_BC = true; 
     }
   else if ( j + int((*local_min_index_in_global)[component](1)) ==
   	int((*global_max_index)[component](1)) )
     {
      if ( component == 1 )
      {
       j -= 1;
       offset[1] = -(*local_cell_size)[component][1](j);
      }
      else
      {
       j = center(1) + stencil(1) + 1;
       j +=  - 2 * ( j + int((*local_min_index_in_global)[component](1)) 
      		- int((*global_max_index)[component](1)) );
       offset[1] = ( 1 + 2 * ( j + int((*local_min_index_in_global)[component]
     	(1)) - int((*global_max_index)[component](1)) ) ) * 
     	(*local_cell_size)[component][1](j);		
      }
      outdomain_or_BC = true; 
     }	
   
   else if ( j + int((*local_min_index_in_global)[component](1)) >
   	int((*global_max_index)[component](1)) )
     {
      if ( component == 1 ) 
       j +=  - 2 * ( j + int((*local_min_index_in_global)[component](1)) 
      		- int((*global_max_index)[component](1)) ) - 1;
      else 
       j +=  - 2 * ( j + int((*local_min_index_in_global)[component](1)) 
      		- int((*global_max_index)[component](1)) );
		
      offset[1] = ( 1 + 2 * ( j + int((*local_min_index_in_global)[component]
     	(1)) - int((*global_max_index)[component](1)) ) ) * 
     	(*local_cell_size)[component][1](j); 
	
      outdomain_or_BC = true;
     } 
   
   if ( DIM == 3 )
   {
    if ( k + int((*local_min_index_in_global)[component](2)) < 0 )
    {
      if ( component == 2 )  k += - 2 * k + 1 ;
      else k += -2 * k ;
      
      offset[2] = double( 2 * k - 1 ) * (*local_cell_size)[component][2](k);    
      outdomain_or_BC = true;    
    }
    else if ( k + int((*local_min_index_in_global)[component](2)) == 0 )
    {
      if ( component == 2 )     
      {
       k += 1;
       offset[2] = (*local_cell_size)[component][2](k);
      }
      else
      {
       k = center(2) - stencil(2) - 1;
       k += -2 * k;
       offset[2] = double( 2 * k - 1 ) * (*local_cell_size)[component][2](k); 
      } 
      outdomain_or_BC = true; 
    }
    else if ( k + int((*local_min_index_in_global)[component](2)) ==
   	int((*global_max_index)[component](2)) )
    {
      if ( component == 2 )
      {
       k -= 1;
       offset[2] = -(*local_cell_size)[component][2](k);
      }
      else
      {
       k = center(2) + stencil(2) + 1;
       k +=  - 2 * ( k + int((*local_min_index_in_global)[component](2)) 
      		- int((*global_max_index)[component](2)) );
       offset[2] = ( 1 + 2 * ( k + int((*local_min_index_in_global)[component]
     	(2)) - int((*global_max_index)[component](2)) ) ) * 
     	(*local_cell_size)[component][2](k);		
      }
      outdomain_or_BC = true; 
     }
    else if ( k + int((*local_min_index_in_global)[component](2)) >
   	int((*global_max_index)[component](2)) )
    {
      if ( component == 2 ) 
       k +=  - 2 * ( k + int((*local_min_index_in_global)[component](2)) 
      		- int((*global_max_index)[component](2)) ) - 1;
      else 
       k +=  - 2 * ( k + int((*local_min_index_in_global)[component](2)) 
      		- int((*global_max_index)[component](2)) );
		
      offset[2] = ( 1 + 2 * ( k + int((*local_min_index_in_global)[component]
     	(2)) - int((*global_max_index)[component](2)) ) ) * 
     	(*local_cell_size)[component][2](k); 
	
      outdomain_or_BC = true;
     }
   }
   
   return ( outdomain_or_BC ) ;

}  
