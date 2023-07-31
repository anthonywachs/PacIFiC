#include <FV_Mesh.hh>
#include <FV.hh>
#include <FV_DomainBuilder.hh>
#include <MAC_assertions.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Module.hh>
#include <MAC_Error.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_GroupExp.hh>
#include <MAC_Variable.hh>
#include <MAC.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Int.hh>
#include <MAC_Bool.hh>
#include <MAC_Double.hh>
#include <iostream>
#include <sstream>
#include <list>
#include <set>
#include <cmath>

using std::cout ; 
using std::endl ;
using std::ostringstream ;
using std::list ;
using std::set ;

struct FV_Mesh_ERROR
{
   static void n1( MAC_ModuleExplorer const* exp ) ;
   static void n2( MAC_ModuleExplorer const* exp ) ; 
   static void n3( MAC_ModuleExplorer const* exp ) ;      
} ;

stringVector* FV_Mesh::directionName = NULL;

//----------------------------------------------------------------------
FV_Mesh*
FV_Mesh:: create( MAC_Object* a_owner,
       	MAC_ModuleExplorer const* exp,
	size_t dim )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: create" ) ;

   if ( !directionName )
   {
     directionName = new stringVector(3);
     (*directionName)(0)="X";
     (*directionName)(1)="Y";   		    
     (*directionName)(2)="Z";
   }

   FV_Mesh* result =
                    new FV_Mesh( a_owner, exp, dim ) ;
      
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_Mesh:: FV_Mesh( MAC_Object* a_owner,
       	MAC_ModuleExplorer const* exp,
	size_t dim )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , global_main_coordinates( NULL )
   , global_number_of_cells( 1 )
   , global_max_index( 0 )
   , global_min_index_in_domain( 0 )
   , global_max_index_in_domain( 0 )
   , local_main_coordinates( NULL )
   , local_number_of_cells( 1 )
   , local_max_index_in_global( 0 )
   , local_min_index_in_global( 0 )
   , local_max_index_in_global_on_current_proc( 0 )    
   , local_min_index_in_global_on_current_proc( 0 )    
   , security_bandwidth( 0 )
   , min_coordinate_on_current_processor( 0 )
   , max_coordinate_on_current_processor( 0 )
   , min_coordinate_with_halozone( 0 )
   , max_coordinate_with_halozone( 0 )
   , main_domain_min_coordinates( 0 )
   , main_domain_max_coordinates( 0 )
   , smallest_grid_size( 1e20 )
   , translation_direction( 0 )
   , translation_magnitude( 0. )
   , translation_distance( 0. )
   , translation_projection( false )
   , periodic( NULL )
   , periodic_flow_direction( 0 )
   , periodic_pressure_drop( NULL )
   , periodic_flow_rate( NULL ) 
   , structured_cartesian_Comm( NULL ) 
   , MPI_coordinates( NULL ) 
   , domain_decomposition ( NULL )     
{
   MAC_LABEL( "FV_Mesh:: FV_Mesh" ) ;

   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   size_t RANK = macCOMM->rank(); 
   size_t NB_RANKS = macCOMM->nb_ranks(); 
   std::string space( 6, ' ' ) ;
   
   if( !exp->has_module( "FV_Mesh" )  )
      FV_Mesh_ERROR::n1( exp ) ;

   if ( RANK == 0 )
     FV::out() << endl << "*** building FV_Mesh" << endl << endl;
   
   // Create FV structured mesh
   // --------------------------
   MAC_ModuleExplorer* ee = exp->create_subexplorer( 0, "FV_Mesh" ) ;
   doubleVector work_double(1,0.);
   global_main_coordinates = new vector< doubleVector >(dim,work_double); 
   main_domain_min_coordinates.re_initialize( dim, 0. );
   main_domain_max_coordinates.re_initialize( dim, 0. );   

   doubleVector const& xcoords__ =
      ee->doubleVector_data( "vertices_coordinate_0" ) ;
   (*global_main_coordinates)[0].re_initialize(xcoords__.size(),0.);
   for (size_t i=0;i<xcoords__.size();++i)
     (*global_main_coordinates)[0](i) = xcoords__(i);
   main_domain_min_coordinates(0) = (*global_main_coordinates)[0](0) ;
   main_domain_max_coordinates(0) = 
   	(*global_main_coordinates)[0]((*global_main_coordinates)[0].size()-1) ;
   global_number_of_cells *= (*global_main_coordinates)[0].size()-1;    

   doubleVector const& ycoords__ =
      ee->doubleVector_data( "vertices_coordinate_1" ) ;   
   (*global_main_coordinates)[1].re_initialize(ycoords__.size(),0.);
   for (size_t i=0;i<ycoords__.size();++i)
     (*global_main_coordinates)[1](i) = ycoords__(i); 
   main_domain_min_coordinates(1) = (*global_main_coordinates)[1](0) ;
   main_domain_max_coordinates(1) = 
   	(*global_main_coordinates)[1]((*global_main_coordinates)[1].size()-1) ; 
   global_number_of_cells *= (*global_main_coordinates)[1].size()-1;

   if ( dim == 3 )
   {
     doubleVector const& zcoords__ =
       ee->doubleVector_data( "vertices_coordinate_2" ) ;
     (*global_main_coordinates)[2].re_initialize(zcoords__.size(),0.);
     for (size_t i=0;i<zcoords__.size();++i)
       (*global_main_coordinates)[2](i) = zcoords__(i); 
     main_domain_min_coordinates(2) = (*global_main_coordinates)[2](0) ;
     main_domain_max_coordinates(2) = 
   	(*global_main_coordinates)[2]((*global_main_coordinates)[2].size()-1) ; 
     global_number_of_cells *= (*global_main_coordinates)[2].size()-1;
   }
   ee->destroy() ; ee=0 ;

   global_max_index.re_initialize( dim, 0 );
   global_min_index_in_domain.re_initialize( dim, 0 );
   global_max_index_in_domain.re_initialize( dim, 0 );
   local_min_index_in_global.re_initialize( dim, 0 );
   local_max_index_in_global.re_initialize( dim, 0 );   
   local_min_index_in_global_on_current_proc.re_initialize( dim, 0 );
   local_max_index_in_global_on_current_proc.re_initialize( dim, 0 );   
   for (size_t i=0;i<dim;++i)
   { 
     global_max_index(i) = (*global_main_coordinates)[i].size()-1;
     global_max_index_in_domain(i) = global_max_index(i);
   }

   // Periodicity
   // -----------
   periodic = new boolVector( dim, false ) ;
   bool build_periodicity = false ;
   size_t_vector periodic_directions( 1 ) ;
   if ( exp->has_module( "periodicity" ) )
   {   
     build_periodicity = true ;
     MAC_ModuleExplorer* eee = exp->create_subexplorer( 0, 
     	"periodicity" ) ;

     periodic_directions = eee->intVector_data( "directions" ) ;
     MAC_ASSERT( periodic_directions.size() <= dim ) ;     
     for (size_t p=0;p<periodic_directions.size();p++)
     {
       MAC_ASSERT( periodic_directions(p) < dim ) ;
       for (size_t j=p+1;j<periodic_directions.size();++j)
         MAC_ASSERT( periodic_directions(p) != periodic_directions(j) ) ;
     }
     
     for (size_t p=0;p<periodic_directions.size();p++)
       (*periodic)(periodic_directions(p)) = true ;
       
     if ( eee->has_entry( "periodic_flow_direction" ) )
     {
       periodic_flow_direction = eee->int_data( "periodic_flow_direction" ) ;
       bool flow_direction_ok = false ;
       for (size_t p=0;p<periodic_directions.size() && !flow_direction_ok;p++)
         if ( periodic_directions(p) == periodic_flow_direction ) 
           flow_direction_ok = true ;
       if ( !flow_direction_ok ) FV_Mesh_ERROR:: n3( eee ); 
     } 
      
     if ( eee->has_entry( "periodic_pressure_drop" ) )
     {
       periodic_pressure_drop = new double ;
       *periodic_pressure_drop = eee->double_data( "periodic_pressure_drop" ) ;
     }               

     if ( eee->has_entry( "periodic_flow_rate" ) )
     {
       periodic_flow_rate = new double ;
       *periodic_flow_rate = eee->double_data( "periodic_flow_rate" ) ;
     }
     
     if ( periodic_pressure_drop && periodic_flow_rate )
       FV_Mesh_ERROR:: n3( eee );
       
     if ( periodic_flow_rate ) 
     {
       periodic_pressure_drop = new double ;
       *periodic_pressure_drop = 0. ;
       // Added by Aashish Goyal on 29Mar2023 for DS solver
       if (eee->has_entry( "initial_pressure_drop" ))
          *periodic_pressure_drop = eee->double_data( "initial_pressure_drop" );
     }
                   
     eee->destroy() ; eee=0 ;          
   }

   // Splitting strategy
   // ------------------
   local_main_coordinates = new vector< doubleVector >(dim,work_double);
   boolVector work_bool(1,false);
   on_current_processor = new vector< boolVector >(dim,work_bool);    
   if ( !exp->has_module( "splitting_strategy" ) )
   {
     if ( NB_RANKS != 1 )
       FV_Mesh_ERROR::n2( exp ) ; 
     for (size_t i=0;i<dim;++i)
     {
       (*local_main_coordinates)[i].re_initialize(
       	(*global_main_coordinates)[i].size(),0.);
       (*on_current_processor)[i].re_initialize(
       	(*global_main_coordinates)[i].size(),true);	
       for (size_t j=0;j<(*global_main_coordinates)[i].size();++j)	
         (*local_main_coordinates)[i](j) = (*global_main_coordinates)[i](j);
     }
     local_number_of_cells = global_number_of_cells;
     for (size_t i=0;i<dim;++i)
     {
       local_max_index_in_global(i) = global_max_index(i);
       local_max_index_in_global_on_current_proc(i) = global_max_index(i);
     }
   }
   else
   {
     MAC_ModuleExplorer* eee = exp->create_subexplorer( 0, 
     	"splitting_strategy" ) ;

     // Halozone width
     if ( eee->has_entry( "security_bandwidth" ) )
       security_bandwidth = eee->int_data( "security_bandwidth" );
     if ( NB_RANKS > 1 )
       eee->test_data( "security_bandwidth", "security_bandwidth>0" ) ;  

     set<size_t> empty_set;
     vector< set<size_t> > local_indices(dim,empty_set);
     
     if ( eee->has_entry( "coordinate_splitting_formula" ) )
     {
//        // Context :
//        MAC_ContextSimple* ctx = MAC_ContextSimple::create( 0 ) ;
//        MAC_DoubleVector* coords = MAC_DoubleVector::create( ctx, dim ) ;
//        ctx->extend( MAC_Variable::object( "DV_X" ), coords ) ;
// 
//        // Check formula expression:
//        MAC_Data* formula =
//           eee->abstract_data( 0, "coordinate_splitting_formula", ctx ) ;
//        if( !formula->value_can_be_evaluated(0) )
//        {
//           MAC_Error::object()->raise_not_evaluable(
//          eee, "coordinate_splitting_formula",
//          formula->undefined_variables(0) ) ;
//        }
//        if( formula->data_type() != MAC_Data::Int )
//        {
//           MAC_Error::object()->raise_bad_data_type(
//          eee, "coordinate_splitting_formula", MAC_Data::Int ) ;
//        }
// 
//        // Optimized evaluation of MAC_Group expressions:
//        MAC_GroupExp::set_optimized_evaluation() ; 
// 
//        // Loop on cells
//        size_t rank_formula; 
//        doubleVector center_coord(dim) ;
//        for (size_t i=0;i<(*global_main_coordinates)[0].size()-1;++i)
//        {
//          center_coord(0) = 0.5 * ((*global_main_coordinates)[0](i)
//        		+ (*global_main_coordinates)[0](i+1));
//          for (size_t j=0;j<(*global_main_coordinates)[1].size()-1;++j)
//          {
//            center_coord(1) = 0.5 * ((*global_main_coordinates)[1](j)
//        		+ (*global_main_coordinates)[1](j+1));
// 	   if ( dim == 3 )
// 	   {
// 	     for (size_t k=0;k<(*global_main_coordinates)[2].size()-1;++k)
//              {
//                center_coord(2) = 0.5 * ((*global_main_coordinates)[2](k)
//        		+ (*global_main_coordinates)[2](k+1));
// 	       coords->set( center_coord ) ;
// 	       rank_formula = point_owner( formula, NB_RANKS ) ;
// 	       if ( rank_formula == RANK )
// 	       {
// 	         local_indices[0].insert(i);
// 	         local_indices[0].insert(i+1);
// 	         local_indices[1].insert(j);
// 	         local_indices[1].insert(j+1);	
// 	         local_indices[2].insert(k);
// 	         local_indices[2].insert(k+1);
// 	       }
// 	     }	              
// 	   }
// 	   else
// 	   {
//              coords->set( center_coord ) ;
// 	     rank_formula = point_owner( formula, NB_RANKS ) ;
// 	     if ( rank_formula == RANK )
// 	     {
// 	       local_indices[0].insert(i);
// 	       local_indices[0].insert(i+1);
// 	       local_indices[1].insert(j);
// 	       local_indices[1].insert(j+1);	
// 	     }	 
// 	   }
//          }
//        }	                
// 
//        MAC_GroupExp::unset_optimized_evaluation() ;
//    
//        ctx->destroy() ; ctx = 0 ;
//        formula->destroy() ; formula = 0 ;

       MAC_Error::object()->raise_module_error(	eee,
	"Splitting strategy \"coordinate_splitting_formula\" is not supported"
	" anymore; use \"domain_decomposition\" instead" ); 
       
     }
     else if ( eee->has_entry( "domain_decomposition" ) )
     {
       size_t_vector dd = eee->intVector_data( "domain_decomposition" );
       MAC_ASSERT( dd.size() == dim ) ;
       size_t total_nb_subdomains = 1 ;
       for (size_t i=0;i<dim;++i) total_nb_subdomains *= dd(i) ;
       MAC_ASSERT( total_nb_subdomains == NB_RANKS ) ;       

       // Create structured cartesian communicator
       structured_cartesian_Comm = new MPI_Comm;
       int reorganisation = 0 ;
       int* period = new int[dim];
       for (size_t i=0;i<dim;++i) period[i] = (*periodic)(i) ;
       domain_decomposition = new int[dim];
       for (size_t i=0;i<dim;++i) domain_decomposition[i] = dd(i) ;       
       MPI_Cart_create( MPI_COMM_WORLD, dim, domain_decomposition, period, 
                reorganisation, structured_cartesian_Comm );
       delete [] period ;
       //delete [] domain_decomposition ;
     
       // Get MPI coordinates & rank in structured cartesian communicator (SSC)
       MPI_coordinates = new int[dim] ;
       MPI_Cart_coords( *structured_cartesian_Comm, RANK, dim, 
       		MPI_coordinates );
		
       // Determine MPI neighbors in SSC
       int* MPI_neighbor_coordinates = new int[dim] ;
       int MPI_neighbor_rank = 0 ;
       boolVector periodic_neighbor( dim, false );
       for (int i=-1;i<2;++i)
       {
         periodic_neighbor(0) = false ;
	 MPI_neighbor_coordinates[0] = MPI_coordinates[0] + i ;
	 if ( (*periodic)(0) )
	 { 
	   if ( MPI_neighbor_coordinates[0] == -1 )
	   {
	     MPI_neighbor_coordinates[0] = int(dd(0)) - 1 ;
	     periodic_neighbor(0) = true ;
	   }
	   else if ( MPI_neighbor_coordinates[0] == int(dd(0)) )
	   {
	     MPI_neighbor_coordinates[0] = 0 ;
	     periodic_neighbor(0) = true ;
	   }
	 }	   

	 if ( MPI_neighbor_coordinates[0] >= 0
	   	 && MPI_neighbor_coordinates[0] < int(dd(0)) )
         {
           for (int j=-1;j<2;++j)
	   {
	     periodic_neighbor(1) = false ;
	     MPI_neighbor_coordinates[1] = MPI_coordinates[1] + j ;
	     if ( (*periodic)(1) ) 
	     {
	       if ( MPI_neighbor_coordinates[1] == -1 )
	       {
	         MPI_neighbor_coordinates[1] = int(dd(1)) - 1 ;
	         periodic_neighbor(1) = true ;
	       }
	       else if ( MPI_neighbor_coordinates[1] == int(dd(1)) )
	       {
	         MPI_neighbor_coordinates[1] = 0 ;
	         periodic_neighbor(1) = true ;
	       }
	     }

	     if ( MPI_neighbor_coordinates[1] >= 0
	   	 && MPI_neighbor_coordinates[1] < int(dd(1)) )
             {
	       if ( dim == 2 )
	       {
	         if ( i || j )
	         {
	           MPI_Cart_rank( *structured_cartesian_Comm, 
	       		MPI_neighbor_coordinates, &MPI_neighbor_rank );
		   if ( MPI_neighbor_rank != int(RANK) )
		     if ( periodic_neighbor(0) || periodic_neighbor(1) )
		       MPI_periodic_neighbors_World.push_back( 
		       	MPI_neighbor_rank ) ;
	             else
		       MPI_neighbors_World.push_back( MPI_neighbor_rank ) ;
	         }
	       }
	       else
	         for (int k=-1;k<2;++k)
	         {
	           periodic_neighbor(2) = false ;
		   MPI_neighbor_coordinates[2] = MPI_coordinates[2] + k ;
		   if ( (*periodic)(2) ) 
	             if ( MPI_neighbor_coordinates[2] == -1 )
	             {
	               MPI_neighbor_coordinates[2] = int(dd(2)) - 1 ;
	               periodic_neighbor(2) = true ;
	             }
	             else if ( MPI_neighbor_coordinates[2] == int(dd(2)) )
	             {
	               MPI_neighbor_coordinates[2] = 0 ;
	               periodic_neighbor(2) = true ;
	             }
		   		   
	           if ( MPI_neighbor_coordinates[2] >= 0
	   	 	&& MPI_neighbor_coordinates[2] < int(dd(2)) )
	             if ( i || j || k )
	             {
	               MPI_Cart_rank( *structured_cartesian_Comm, 
	       		MPI_neighbor_coordinates, &MPI_neighbor_rank );
	               if ( MPI_neighbor_rank != int(RANK) )
		         if ( periodic_neighbor(0) || periodic_neighbor(1) 
			 	|| periodic_neighbor(2) )
		           MPI_periodic_neighbors_World.push_back( 
		       		MPI_neighbor_rank ) ;
	                 else
		           MPI_neighbors_World.push_back( MPI_neighbor_rank ) ;
	             }
	         }
	     }
           }
	 }
	 MPI_neighbors_World.sort();
	 MPI_neighbors_World.unique();
	 if ( !MPI_periodic_neighbors_World.empty() )
	 {
	   MPI_periodic_neighbors_World.sort();
	   MPI_periodic_neighbors_World.unique();
	 }
       }       
	   	       
       // Indices on this process per direction
       for (size_t i=0;i<dim;++i)
       {
         size_t_vector sizes_sequence = subvector_sizes_sequence( 
	 	(*global_main_coordinates)[i].size() - 1 , dd(i) ) ;
	 size_t index_min = 0 ;
	 for (size_t j=0;j<size_t(MPI_coordinates[i]);++j) 
	   index_min += sizes_sequence(j) ;
	 size_t index_max = index_min + sizes_sequence( MPI_coordinates[i] ) ;	
       	 for (size_t j=index_min;j<=index_max;++j)
	   local_indices[i].insert(j);	
       }
     }         
     eee->destroy() ; eee=0 ;

     set<size_t>::iterator is;
     local_number_of_cells = 1;
     for (size_t i=0;i<dim;++i)
     {
       size_t down = 0, up = 0;
       local_min_index_in_global(i) = *local_indices[i].begin();
       local_min_index_in_global_on_current_proc(i) = *local_indices[i].begin();
       for (size_t j=0;j<security_bandwidth;++j)
         if (local_min_index_in_global(i) > 0) 
	 {
	   local_min_index_in_global(i)--;
	   down++;
	 }
	 
       local_max_index_in_global(i) = *local_indices[i].rbegin();	
       local_max_index_in_global_on_current_proc(i) = 
       		*local_indices[i].rbegin();
       for (size_t j=0;j<security_bandwidth;++j)
         if (local_max_index_in_global(i) < global_max_index(i)) 
	 {
	   local_max_index_in_global(i)++;       
           up++;
	 }
	 
       size_t nb_local_coor = local_max_index_in_global(i)
       		- local_min_index_in_global(i)+1;
       (*local_main_coordinates)[i].re_initialize(nb_local_coor,0.);
       (*on_current_processor)[i].re_initialize(nb_local_coor,true);
       local_number_of_cells *= nb_local_coor - 1;		
       for (size_t j=0;j<nb_local_coor;++j)
         (*local_main_coordinates)[i](j) =
	 	(*global_main_coordinates)[i](local_min_index_in_global(i)+j);
       for (size_t j=0;j<down;++j)
         (*on_current_processor)[i](j)=false;
       for (size_t j=0;j<up;++j)
         (*on_current_processor)[i](nb_local_coor-1-j)=false;       		
     }              
   }

   // Periodicity building
   // --------------------
   if( build_periodicity )
   {        
     MAC_ASSERT( security_bandwidth > 0 );
     for (size_t p=0;p<periodic_directions.size();p++)
     {
       double periodic_cell_size = 0.5 * ( 
     	(*global_main_coordinates)[periodic_directions(p)](1) -  
   	(*global_main_coordinates)[periodic_directions(p)](0) +
	(*global_main_coordinates)[periodic_directions(p)](
		(*global_main_coordinates)[periodic_directions(p)].size()-1) -
	(*global_main_coordinates)[periodic_directions(p)](
		(*global_main_coordinates)[periodic_directions(p)].size()-2) ) ;

       // Modify global coordinates vector
       doubleVector global_temp = 
       	(*global_main_coordinates)[periodic_directions(p)] ;
       (*global_main_coordinates)[periodic_directions(p)].resize( 
       		global_temp.size()+ 2 * security_bandwidth, 0. ) ;
       for (size_t j=0;j<security_bandwidth;++j)
         (*global_main_coordinates)[periodic_directions(p)](j) = 
	global_temp(0) - ( security_bandwidth - j ) * periodic_cell_size ;
       for (size_t j=0;j<global_temp.size();++j) 
         (*global_main_coordinates)[periodic_directions(p)]
	 	(j+security_bandwidth) = global_temp(j) ;
       for (size_t j=0;j<security_bandwidth;++j)
         (*global_main_coordinates)[periodic_directions(p)](
	 	j+global_temp.size()+security_bandwidth) = 
			global_temp(global_temp.size()-1) 
			+ ( j + 1 ) * periodic_cell_size ;
	
       global_max_index(periodic_directions(p)) += 2 * security_bandwidth ;
       global_min_index_in_domain(periodic_directions(p)) = security_bandwidth;
       global_max_index_in_domain(periodic_directions(p)) += security_bandwidth;
       global_number_of_cells = 1 ;
       for (size_t i=0;i<dim;++i)
         global_number_of_cells *= (*global_main_coordinates)[i].size()-1;
       
       // Modify local coordinates vector
       if ( local_min_index_in_global(periodic_directions(p)) == 0
     	&& local_max_index_in_global(periodic_directions(p)) 
		== global_temp.size()-1 )
       {
         local_min_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ; 
         local_max_index_in_global(periodic_directions(p)) 
       		+= 2 * security_bandwidth ;
         local_max_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ; 
		
         (*local_main_coordinates)[periodic_directions(p)].resize( 
       	(*global_main_coordinates)[periodic_directions(p)].size(), 0. ) ;
         for (size_t j=0;
	 	j<(*local_main_coordinates)[periodic_directions(p)].size();++j)	
           (*local_main_coordinates)[periodic_directions(p)](j) = 
	 	(*global_main_coordinates)[periodic_directions(p)](
		local_min_index_in_global(periodic_directions(p))+j);
		
         (*on_current_processor)[periodic_directions(p)].resize( 
       	(*global_main_coordinates)[periodic_directions(p)].size(), true );
         for (size_t j=0;j<security_bandwidth;++j)
         {
           (*on_current_processor)[periodic_directions(p)](j) = false ;
           (*on_current_processor)[periodic_directions(p)](
	 	(*on_current_processor)[periodic_directions(p)].size()-j-1) 
		= false ;       
         }
       
         local_number_of_cells = 1 ;
         for (size_t i=0;i<dim;++i)
           local_number_of_cells *= (*local_main_coordinates)[i].size()-1;
	 
       }
       else if ( local_min_index_in_global(periodic_directions(p)) == 0 )
       {
         local_min_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ; 
         local_max_index_in_global(periodic_directions(p)) 
       		+= security_bandwidth ;
         local_max_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ;

         (*local_main_coordinates)[periodic_directions(p)].resize( 
       		(*local_main_coordinates)[periodic_directions(p)].size()
		+ security_bandwidth, 0. ) ;
         for (size_t j=0;
	 	j<(*local_main_coordinates)[periodic_directions(p)].size();++j)	
           (*local_main_coordinates)[periodic_directions(p)](j) = 
	 	(*global_main_coordinates)[periodic_directions(p)](
		local_min_index_in_global(periodic_directions(p))+j);

         boolVector current_temp = 
	 	(*on_current_processor)[periodic_directions(p)] ;
         (*on_current_processor)[periodic_directions(p)].resize( 
       		(*local_main_coordinates)[periodic_directions(p)].size(), 
		true );
         for (size_t j=0;j<security_bandwidth;++j)
           (*on_current_processor)[periodic_directions(p)](j) = false ;
         for (size_t j=0;j<current_temp.size();++j)
           (*on_current_processor)[periodic_directions(p)](j+security_bandwidth)
	   	= current_temp(j) ;

         local_number_of_cells = 1 ;
         for (size_t i=0;i<dim;++i)
           local_number_of_cells *= (*local_main_coordinates)[i].size()-1; 	
       }
       else if ( local_max_index_in_global(periodic_directions(p)) == 
     		global_temp.size()-1 )
       {
         local_min_index_in_global(periodic_directions(p)) 
       		+= security_bandwidth ;
         local_min_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ;       
         local_max_index_in_global(periodic_directions(p)) 
       		+= 2 * security_bandwidth ;
         local_max_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ; 
		
         (*local_main_coordinates)[periodic_directions(p)].resize( 
       		(*local_main_coordinates)[periodic_directions(p)].size()
		+ security_bandwidth, 0. ) ;
         for (size_t j=0;
	 	j<(*local_main_coordinates)[periodic_directions(p)].size();++j)	
           (*local_main_coordinates)[periodic_directions(p)](j) = 
	 	(*global_main_coordinates)[periodic_directions(p)](
		local_min_index_in_global(periodic_directions(p))+j);
		
         boolVector current_temp = 
	 	(*on_current_processor)[periodic_directions(p)] ;
       	 (*on_current_processor)[periodic_directions(p)].resize( 
       		(*local_main_coordinates)[periodic_directions(p)].size(), 
		true );
         for (size_t j=0;j<current_temp.size();++j)
           (*on_current_processor)[periodic_directions(p)](j) 
	   	= current_temp(j) ;
         for (size_t j=0;j<security_bandwidth;++j)
           (*on_current_processor)[periodic_directions(p)](
	 	(*on_current_processor)[periodic_directions(p)].size()-j-1) 
		= false ;
	 
         local_number_of_cells = 1 ;
         for (size_t i=0;i<dim;++i)
           local_number_of_cells *= (*local_main_coordinates)[i].size()-1;
       }
       else
       {
         local_min_index_in_global(periodic_directions(p)) 
       		+= security_bandwidth ;
         local_min_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ;       
         local_max_index_in_global(periodic_directions(p)) 
       		+= security_bandwidth ;
         local_max_index_in_global_on_current_proc(periodic_directions(p)) 
       		+= security_bandwidth ;      
       }       
     }	         				
   }
   
   // Domain size on current processor
   // --------------------------------
   min_coordinate_on_current_processor.re_initialize( dim, 0. );
   max_coordinate_on_current_processor.re_initialize( dim, 0. ); 
   for (size_t i=0;i<dim;++i) 
   {
     min_coordinate_on_current_processor(i) = 
     	(*global_main_coordinates)[i](
	local_min_index_in_global_on_current_proc(i) );
     max_coordinate_on_current_processor(i) = 
     	(*global_main_coordinates)[i](
	local_max_index_in_global_on_current_proc(i) );	
   } 
   
   min_coordinate_with_halozone.re_initialize( dim, 0. );
   max_coordinate_with_halozone.re_initialize( dim, 0. ); 
   for (size_t i=0;i<dim;++i) 
   {
     min_coordinate_with_halozone(i) = 
     	(*global_main_coordinates)[i]( local_min_index_in_global(i) );
     max_coordinate_with_halozone(i) = 
     	(*global_main_coordinates)[i]( local_max_index_in_global(i) );	
   }
   
   // Smallest grid size in all directions
   // ------------------------------------
   double cell_size ;      
   for (size_t i=0;i<dim;++i)
     for (size_t j=1;j<global_max_index(i)+1;++j) 
     {
       cell_size = (*global_main_coordinates)[i](j)-
       	(*global_main_coordinates)[i](j-1) ;
       smallest_grid_size = cell_size < smallest_grid_size ? cell_size : 
       		smallest_grid_size ;
     }   

   // Screen outputs
   // --------------
   if ( RANK == 0 )
   {   
     FV::out() << space << "Global number of cells = " << 
     	global_number_of_cells << endl; 
     for (size_t i=0;i<dim;++i)
       FV::out() << space << "Global number of points in direction " << 
     	(*directionName)(i) << " = " << (*global_main_coordinates)[i].size() 
	<< endl;
     FV::out() << endl; 
   } 

   MAC_Exec::communicator()->barrier();
   for (size_t i=0;i<NB_RANKS;++i)
   {
     if ( i == RANK ) 
       FV::out() << space << "Rank " << RANK << ": local number of cells = "
	 	<< local_number_of_cells << endl;     
     MAC_Exec::communicator()->barrier();
   } 
   if ( RANK == 0 ) FV::out() << endl << endl ;
   MAC_Exec::communicator()->barrier();
   
   // Translation-projection
   // ----------------------
   if( exp->has_module( "FV_TranslationProjection" ) )
   {
     MAC_ModuleExplorer* tp = exp->create_subexplorer( 0, 
     	"FV_TranslationProjection" ) ;
     translation_projection = true ;
     translation_direction = tp->int_data( "direction" ) ;
     MAC_ASSERT( IMPLIES( dim == 2, translation_direction < 2 ) ) ;
     translation_magnitude = tp->double_data( "magnitude" ) ;
     tp->destroy() ; tp=0 ;
     if ( RANK == 0 )
     {
       FV::out() << endl << "*** FV_Mesh translation-projection" << endl 
       	<< endl;
       FV::out() << space << "Direction = " << translation_direction << endl;
       FV::out() << space << "Magnitude = " << translation_magnitude << endl
       	<< endl;       
     }          
   } 
}




//----------------------------------------------------------------------
FV_Mesh:: ~FV_Mesh( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: ~FV_Mesh" ) ;

   global_main_coordinates->clear();
   delete global_main_coordinates;
   local_main_coordinates->clear();
   delete local_main_coordinates;
   on_current_processor->clear();
   delete on_current_processor;    
   if ( directionName ) 
   {
     delete directionName;
     directionName = NULL; 
   }
   delete periodic;
   if ( periodic_pressure_drop ) delete periodic_pressure_drop ;
   if ( periodic_flow_rate ) delete periodic_flow_rate ;
   if ( structured_cartesian_Comm ) 
   {
     delete structured_cartesian_Comm ;
     delete [] MPI_coordinates ;
     delete [] domain_decomposition ;
     MPI_neighbors_World.clear();
     MPI_periodic_neighbors_World.clear();
   }   
}




//----------------------------------------------------------------------
void FV_Mesh:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: print" ) ;

   std::string space( indent_width, ' ' ) ;   
   size_t dim = global_main_coordinates->size();
   size_t rank_world = MAC_Exec::communicator()->rank();

   os << space << "Global mesh (on all processors)" << endl;
   os << space << "-------------------------------" << endl;
   for (size_t i=0;i<dim;++i)
   {
     os << space << "Global number of points in direction " << 
     	(*directionName)(i) <<
     	" = " << (*global_main_coordinates)[i].size() << endl;
     os << space << "Min index = " << 0 << "  Max index = " 
     	<< global_max_index(i) << endl;
     os << space << "Coordinate values = ";
     for (size_t j=0;j<(*global_main_coordinates)[i].size();++j)
       os << space << " " << (*global_main_coordinates)[i](j);
     os << endl;
   }
   os << space << "Global number of cells = " << global_number_of_cells << endl;
   os << endl;
   os << space << "MPI Structured layout" << endl;
   os << space << "---------------------" << endl; 
   os << space << "MPI coordinates = " << MPI_coordinates[0] << " " << 
   	MPI_coordinates[1];
   if ( dim == 3 ) os << " " << MPI_coordinates[2];
   os << endl;   
   os << space << "Rank in MPI_COMM_WORLD = " << rank_world << endl;   	 
   if ( dim == 3 ) 
     os << space << "MPI neighbors: i j k rank_world" << endl;
   else os << space << "MPI neighbors: i j rank_world" << endl;
   int* MPI_coords = new int[dim]; 
   list<size_t>::const_iterator ilw ;
   for (ilw=MPI_neighbors_World.begin();ilw!=MPI_neighbors_World.end();ilw++)
   {
     MPI_Cart_coords( *structured_cartesian_Comm, *ilw, dim, MPI_coords );
     os << space << MPI_coords[0] << " " << MPI_coords[1] << " ";
     if ( dim == 3 ) os << MPI_coords[2] << " "; 
     os << *ilw << endl;
   }
   if ( !MPI_periodic_neighbors_World.empty() )
   {
     if ( dim == 3 ) 
       os << space << "MPI periodic neighbors: i j k rank_world" << endl;
     else os << space << "MPI periodic neighbors: i j rank_world" << endl;
     for (ilw=MPI_periodic_neighbors_World.begin();
     	ilw!=MPI_periodic_neighbors_World.end();ilw++)
     {
       MPI_Cart_coords( *structured_cartesian_Comm, *ilw, dim, MPI_coords );
       os << space << MPI_coords[0] << " " << MPI_coords[1] << " ";
       if ( dim == 3 ) os << MPI_coords[2] << " "; 
       os << *ilw << endl;
     }   
   }
   delete [] MPI_coords;
   os << endl;      
   os << space << "Local mesh on processor " << MAC_Exec::communicator()->rank()
   	<< endl;
   os << space << "------------------------------" << endl;
   for (size_t i=0;i<dim;++i)
   {
     os << space << "Local number of points in direction " << 
     	(*directionName)(i) <<
     	" = " << (*local_main_coordinates)[i].size() << endl;
     os << space << "Min index = " << local_min_index_in_global(i) 
     	<< "  Max index = " << local_max_index_in_global(i) << endl;
     os << space << "On current proc: Min coordinate = " 
     	<< min_coordinate_on_current_processor(i) 
     	<< "  Max coordinate = " << max_coordinate_on_current_processor(i) 
	<< endl; 
     os << space << "With halozone: Min coordinate = " 
     	<< min_coordinate_with_halozone(i) 
     	<< "  Max coordinate = " << max_coordinate_with_halozone(i) 
	<< endl; 	    		
     os << space << "Coordinate values = ";
     for (size_t j=0;j<(*local_main_coordinates)[i].size();++j)
       os << space << " " << (*local_main_coordinates)[i](j);
     os << endl;
     os << space << "Coordinates in halo zone = ";
     for (size_t j=0;j<(*on_current_processor)[i].size();++j)
       if (!(*on_current_processor)[i](j)) 
         os << space << " " << (*local_main_coordinates)[i](j);
     os << endl;	 
   }
   os << space << "Local number of cells = " << local_number_of_cells << endl;
   if ( translation_projection )
   {
     os << space << "Translation-Projection" << endl;
     os << space << "----------------------" << endl;
     os << space << "Direction = " << translation_direction << endl;
     os << space << "Magnitude = " << translation_magnitude << endl;
     os << space << "Distance = " << translation_distance << endl;     
   }
}




//----------------------------------------------------------------------
void
FV_Mesh:: write_grid( MAC_Module* base, bool parallel ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: write_grid" ) ;
   std::string pre = ( parallel ? "P" : "" ) ;

   size_t dim = global_main_coordinates->size();
   size_t iter = 0;
   MAC_Module* Coordinates = MAC_Module::create( base, pre+"Coordinates" ) ;
		
   // Vector is always 3D in Paraview whatever the dimension of the grid
   size_t_vector local_nbPoints(3);
   vector< doubleVector > localCoordinates(3,0);
   doubleVector back_translation(3,0.);
   back_translation(translation_direction) = - translation_distance ;
   
   for (size_t i=0;i<dim;++i)
     local_nbPoints(i) = local_max_index_in_global_on_current_proc(i)
   		- local_min_index_in_global_on_current_proc(i) + 1;

   for (size_t i=0;i<dim;++i)
     localCoordinates[i].re_initialize(local_nbPoints(i),0.);
   
   for (size_t i=0;i<dim;++i)
   {
     iter = 0;
     for (size_t j=0;j<(*on_current_processor)[i].size();++j)
     {
       if ((*on_current_processor)[i](j))
       {
         localCoordinates[i](iter) = (*local_main_coordinates)[i](j)
	 	+ back_translation(i) ;
	 iter++;
       }
     }   
   }
   
   if ( dim == 2 )
     localCoordinates[2].re_initialize(1,0.);

   Coordinates->add_entry("X",
   		MAC_DoubleVector::create(Coordinates,localCoordinates[0]) ) ;
   Coordinates->add_entry("Y",
   		MAC_DoubleVector::create(Coordinates,localCoordinates[1]) ) ;
   Coordinates->add_entry("Z",
   		MAC_DoubleVector::create(Coordinates,localCoordinates[2]) ) ;
        
   base->add_module( Coordinates ) ;
   
}




//----------------------------------------------------------------------
size_t
FV_Mesh:: point_owner( MAC_Data const* formula, size_t const& nb_ranks ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: point_owner" ) ;
   MAC_CHECK( formula != 0 ) ;
   MAC_CHECK( formula->data_type() == MAC_Data::Int ) ;
   MAC_CHECK( formula->value_can_be_evaluated(0) ) ;
   
   int const idx = formula->to_int() ;
   if( idx < 0 )
   {
      MAC_Error::object()->raise_plain(
         "*** GE_CoordinateSplitting error:\n"
         "    the expression of keyword \"coordinate_splitting_formula\"\n"
         "    has negative values" ) ;
   }
   if( idx >= (int)nb_ranks )
   {
      MAC_Error::object()->raise_plain(
         "*** GE_CoordinateSplitting error:\n"
         "    the expression of keyword \"coordinate_splitting_formula\"\n"
         "    has some values equal or greater than the number of ranks" ) ;
   }
   size_t const result = idx ;

   MAC_CHECK( result<nb_ranks ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_Mesh:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: nb_space_dimensions" ) ;
   MAC_CHECK( global_main_coordinates != 0 ) ;

   size_t dim = global_main_coordinates->size();

   MAC_CHECK( dim == 2 || dim == 3 ) ;
   return( dim ) ;
}




//----------------------------------------------------------------------
vector< doubleVector > const* 
FV_Mesh:: get_global_main_coordinates( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_global_main_coordinates" ) ;
   MAC_CHECK( global_main_coordinates != 0 ) ;

   vector< doubleVector > const* result = global_main_coordinates;

   return( result ) ;
}  




//----------------------------------------------------------------------
size_t_vector const*
FV_Mesh::get_global_max_index( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_global_max_index" ) ;
   MAC_CHECK( global_max_index.size() != 0 ) ;
      
   size_t_vector const* result = &global_max_index;
   
   return( result ) ;

}




//----------------------------------------------------------------------
size_t_vector const*
FV_Mesh::get_global_min_index_in_domain( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_global_min_index_in_domain" ) ;
   MAC_CHECK( global_min_index_in_domain.size() != 0 ) ;
      
   size_t_vector const* result = &global_min_index_in_domain;
   
   return( result ) ;

}




//----------------------------------------------------------------------
size_t_vector const*
FV_Mesh::get_global_max_index_in_domain( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_global_max_index_in_domain" ) ;
   MAC_CHECK( global_max_index_in_domain.size() != 0 ) ;
      
   size_t_vector const* result = &global_max_index_in_domain;
   
   return( result ) ;

}




//----------------------------------------------------------------------
size_t_vector const* 
FV_Mesh:: get_local_max_index_in_global( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_local_max_index_in_global" ) ;
   MAC_CHECK( local_max_index_in_global.size() != 0 ) ;

   size_t_vector const* result = &local_max_index_in_global;

   return( result ) ;
}




//----------------------------------------------------------------------
size_t_vector const* 
FV_Mesh:: get_local_min_index_in_global( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_local_min_index_in_global" ) ;
   MAC_CHECK( local_min_index_in_global.size() != 0 ) ;

   size_t_vector const* result = &local_min_index_in_global;

   return( result ) ;
}




//----------------------------------------------------------------------
size_t_vector const* 
FV_Mesh:: get_local_max_index_in_global_on_current_proc( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_local_max_index_in_global_on_current_proc" ) ;
   MAC_CHECK( local_max_index_in_global_on_current_proc.size() != 0 ) ;

   size_t_vector const* result = &local_max_index_in_global_on_current_proc;

   return( result ) ;
}




//----------------------------------------------------------------------
size_t_vector const* 
FV_Mesh:: get_local_min_index_in_global_on_current_proc( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_local_min_index_in_global_on_current_proc" ) ;
   MAC_CHECK( local_min_index_in_global.size() != 0 ) ;

   size_t_vector const* result = &local_min_index_in_global_on_current_proc;

   return( result ) ;
}




//----------------------------------------------------------------------
vector< boolVector > const*
FV_Mesh:: get_nodes_owner( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_nodes_owner" ) ;
   MAC_CHECK( on_current_processor != 0 ) ;

   vector< boolVector > const* result = on_current_processor;

   return( result ) ;
}




//----------------------------------------------------------------------
size_t
FV_Mesh:: get_local_nb_points( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_local_nb_points" ) ;

   return( local_max_index_in_global(direction)
   		- local_min_index_in_global(direction) + 1 ) ;
}


  

//----------------------------------------------------------------------
size_t
FV_Mesh:: get_local_nb_points_on_current_proc( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_local_nb_points_on_current_proc" ) ;

   return( local_max_index_in_global_on_current_proc(direction)
   		- local_min_index_in_global_on_current_proc(direction) + 1 ) ;
}


  

//----------------------------------------------------------------------
size_t
FV_Mesh:: get_security_bandwidth( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_security_bandwidth" ) ;

   return( security_bandwidth ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_smallest_grid_size( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_smallest_grid_size" ) ;

   return( smallest_grid_size ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_smallest_constant_grid_size( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_smallest_constant_grid_size" ) ;

   double cell_size, cell_size_im1, cell_size_ip1, result = 1.e20 ;
   size_t dim = global_main_coordinates->size();
      
   for (size_t i=0;i<dim;++i)
     for (size_t j=2;j<global_max_index(i);++j) 
     {
       cell_size = (*global_main_coordinates)[i](j)-
       	(*global_main_coordinates)[i](j-1) ;
       cell_size_im1 = (*global_main_coordinates)[i](j-1)-
       	(*global_main_coordinates)[i](j-2) ;	
       cell_size_ip1 = (*global_main_coordinates)[i](j+1)-
       	(*global_main_coordinates)[i](j) ;	
       
       // Relative variation less than 10-4 between 3 consecutive
       // cells is considered as a constant grid size zone
       if ( fabs( cell_size - cell_size_im1 ) / fabs( cell_size ) < 1e-4
       	&& fabs( cell_size - cell_size_ip1 ) / fabs( cell_size ) < 1e-4
       	&& cell_size < result )
         result = cell_size ;
     }

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_domain_on_current_processor( double const& x, 
	double const& y ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_domain_on_current_processor" ) ;

   bool result = x <= max_coordinate_on_current_processor(0)
   	&& x >= min_coordinate_on_current_processor(0)
	&& y <= max_coordinate_on_current_processor(1)
   	&& y >= min_coordinate_on_current_processor(1) ;

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_domain_on_current_processor( double const& x, 
	double const& y, double const& z ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_domain_on_current_processor" ) ;

   bool result = x <= max_coordinate_on_current_processor(0)
   	&& x >= min_coordinate_on_current_processor(0)
	&& y <= max_coordinate_on_current_processor(1)
   	&& y >= min_coordinate_on_current_processor(1) 
	&& z <= max_coordinate_on_current_processor(2)
   	&& z >= min_coordinate_on_current_processor(2) ;

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_domain_on_current_processor( double const& coor, 
      	size_t direction, double const& tol ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_domain_on_current_processor" ) ;

   bool result = coor <= max_coordinate_on_current_processor(direction) + tol
   	&& coor >= min_coordinate_on_current_processor(direction) - tol ;

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_domain_with_halozone( double const& x, 
	double const& y ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_domain_with_halozone" ) ;

   bool result = x <= max_coordinate_with_halozone(0)
   	&& x >= min_coordinate_with_halozone(0)
	&& y <= max_coordinate_with_halozone(1)
   	&& y >= min_coordinate_with_halozone(1) ;

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_domain_with_halozone( double const& x, 
	double const& y, double const& z ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_domain_with_halozone" ) ;

   bool result = x <= max_coordinate_with_halozone(0)
   	&& x >= min_coordinate_with_halozone(0)
	&& y <= max_coordinate_with_halozone(1)
   	&& y >= min_coordinate_with_halozone(1) 
	&& z <= max_coordinate_with_halozone(2)
   	&& z >= min_coordinate_with_halozone(2) ;

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_domain_with_halozone_plus_ext( double const& x, 
	double const& y, double const& ext ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_domain_with_halozone_plus_ext" ) ;

   bool result = x <= max_coordinate_with_halozone(0) + ext
   	&& x >= min_coordinate_with_halozone(0) - ext
	&& y <= max_coordinate_with_halozone(1) + ext
   	&& y >= min_coordinate_with_halozone(1) - ext ;

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_domain_with_halozone_plus_ext( double const& x, 
	double const& y, double const& z, double const& ext ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_domain_with_halozone_plus_ext" ) ;

   bool result = x <= max_coordinate_with_halozone(0) + ext
   	&& x >= min_coordinate_with_halozone(0) - ext
	&& y <= max_coordinate_with_halozone(1) + ext
   	&& y >= min_coordinate_with_halozone(1) - ext 
	&& z <= max_coordinate_with_halozone(2) + ext
   	&& z >= min_coordinate_with_halozone(2) - ext ;

   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_in_main_domain( double const& coor, size_t direction,
      	double const& tol ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_in_main_domain" ) ;

   bool result = coor <= main_domain_max_coordinates(direction) + tol
   	&& coor >= main_domain_min_coordinates(direction) - tol ;   

   return( result ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_min_coordinate_on_current_processor( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_min_coordinate_on_current_processor" ) ;

   return( min_coordinate_on_current_processor( direction ) ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_min_coordinate_with_halozone( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_min_coordinate_with_halozone" ) ;

   return( min_coordinate_with_halozone( direction ) ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_max_coordinate_on_current_processor( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_max_coordinate_on_current_processor" ) ;

   return( max_coordinate_on_current_processor( direction ) ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_max_coordinate_with_halozone( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_max_coordinate_with_halozone" ) ;

   return( max_coordinate_with_halozone( direction ) ) ;
}




//----------------------------------------------------------------------
void
FV_Mesh:: translation( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: translation" ) ;

   size_t ng = (*global_main_coordinates)[translation_direction].size() ;
   for (size_t i=0;i<ng;++i)
     (*global_main_coordinates)[translation_direction](i) += 
     		translation_magnitude ;
		
   size_t nl = (*local_main_coordinates)[translation_direction].size() ;
   for (size_t i=0;i<nl;++i)
     (*local_main_coordinates)[translation_direction](i) += 
     		translation_magnitude ;
		
   min_coordinate_on_current_processor(translation_direction) += 
     		translation_magnitude ;	
   max_coordinate_on_current_processor(translation_direction) += 
     		translation_magnitude ;						
   min_coordinate_with_halozone(translation_direction) += 
     		translation_magnitude ;	     
   max_coordinate_with_halozone(translation_direction) += 
     		translation_magnitude ;	
		
   translation_distance += translation_magnitude ;           
}




//----------------------------------------------------------------------
size_t
FV_Mesh:: get_translation_direction( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_translation_direction" ) ;

   return( translation_direction ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_translation_magnitude( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_translation_magnitude" ) ;

   return( translation_magnitude ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_translation_distance( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_translation_distance" ) ;

   return( translation_distance ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_translation_active( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_translation_active" ) ;

   return( translation_projection ) ;
}




//----------------------------------------------------------------------
size_t
FV_Mesh:: get_periodic_flow_direction( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_periodic_flow_direction" ) ;

   return( periodic_flow_direction ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_periodic_pressure_drop( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_periodic_pressure_drop" ) ;

   return( *periodic_pressure_drop ) ;
}




//----------------------------------------------------------------------
void
FV_Mesh:: set_periodic_pressure_drop( double const& ppd )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: set_periodic_pressure_drop" ) ;

   *periodic_pressure_drop = ppd ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_periodic_flow_rate( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_periodic_flow_rate" ) ;

   return( *periodic_flow_rate ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_periodic_pressure_drop( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_periodic_pressure_drop" ) ;

   // Return whether a pressure drop is imposed in a periodic direction
   // True if is_periodic_flow_rate() is true
   return( periodic_pressure_drop != NULL ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_periodic_flow_rate( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_periodic_flow_rate" ) ;

   // Return whether a flow rate is imposed in a periodic direction 
   return( periodic_flow_rate != NULL ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_periodic_flow( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_periodic_flow" ) ;

   return( periodic_flow_rate != NULL || periodic_pressure_drop != NULL ) ;
}




//----------------------------------------------------------------------
bool
FV_Mesh:: is_periodic_domain( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: is_periodic_domain" ) ;
 
   // Return whether the domain is periodic is at least one direction
   bool result = false ;
   for (size_t i=0;i<periodic->size() && !result;++i) result = (*periodic)(i);
      
   return( result ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_main_domain_min_coordinate( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_main_domain_min_coordinate" ) ;

   return( main_domain_min_coordinates(direction) ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_main_domain_max_coordinate( size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_main_domain_max_coordinate" ) ;

   return( main_domain_max_coordinates(direction) ) ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_main_domain_boundary_perp_to_direction_measure( 
	size_t direction ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_Mesh:: get_main_domain_boundary_perp_to_direction_measure" ) ;
   
   double result = 0. ;
   size_t dim = global_main_coordinates->size();
   
   if ( dim == 2 )
   {
     switch( direction )
     {
       case 0: 
         result = main_domain_max_coordinates(1)
	 	- main_domain_min_coordinates(1) ;
         break;
       case 1: 
         result = main_domain_max_coordinates(0)
	 	- main_domain_min_coordinates(0) ;
         break;       
     }
   }
   else
   {
     switch( direction )
     {
       case 0: 
         result = ( main_domain_max_coordinates(1)
	 	- main_domain_min_coordinates(1) ) *
		( main_domain_max_coordinates(2)
	 	- main_domain_min_coordinates(2) ) ;
         break;
       case 1: 
         result = ( main_domain_max_coordinates(0)
	 	- main_domain_min_coordinates(0) ) *
		( main_domain_max_coordinates(2)
	 	- main_domain_min_coordinates(2) ) ;
         break;       
       case 2: 
         result = ( main_domain_max_coordinates(1)
	 	- main_domain_min_coordinates(1) ) *
		( main_domain_max_coordinates(0)
	 	- main_domain_min_coordinates(0) ) ;
         break;
     }   
   } 
   
   return( result ) ;
}




//---------------------------------------------------------------------------
bool
FV_Mesh:: between( doubleVector const* mesh, double const& x, size_t& i0 )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Utils:: between" ) ;
   
   bool in = false ;

   if ( x >= (*mesh)(0) && x <= (*mesh)(mesh->size()-1) )
   {
     i0 = 0 ;
     size_t i1 = mesh->size()-1, ibar = 0 ;
     while ( i1-i0 != 1 )
     {
       ibar = size_t( 0.5 * ( i0 + i1 ) );
       if ( x > (*mesh)(ibar) ) i0 = ibar ;
       else i1 = ibar ; 
     }
     in = true ;    
   }
  
   return in ;      
}




//---------------------------------------------------------------------------
bool
FV_Mesh:: between_subinterval( doubleVector const* mesh,
       		double const& x, size_t& i0,
		const size_t& i0_init, const size_t& i1_init )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Utils:: between_subinterval" ) ;
   
   bool in = false ;

   if ( x >= (*mesh)(i0_init) && x <= (*mesh)(i1_init) )
   {
     i0 = i0_init ;
     size_t i1 = i1_init, ibar = 0 ;
     while ( i1-i0 != 1 )
     {
       ibar = size_t( 0.5 * ( i0 + i1 ) );
       if ( x > (*mesh)(ibar) ) i0 = ibar ;
       else i1 = ibar ; 
     }
     in = true ;    
   }
  
   return in ;      
}




//---------------------------------------------------------------------------
size_t
FV_Mesh:: max_index( doubleVector const* mesh, double const& x )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Utils:: max_index" ) ;
   
   size_t imax = 0 ;

   if ( x > (*mesh)(mesh->size()-1) ) imax = mesh->size()-1 ;
   else if ( between( mesh, x, imax ) ) ++imax;

   return imax ;     
} 




//---------------------------------------------------------------------------
size_t
FV_Mesh:: min_index( doubleVector const* mesh, double const& x )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Utils:: min_index" ) ;
   
   size_t imin = mesh->size()-1 ;

   if ( x < (*mesh)(0) ) imin = 0 ;
   else if ( between( mesh, x, imin ) ) {}

   return imin ;      
} 




//---------------------------------------------------------------------------
size_t_vector
FV_Mesh:: subvector_sizes_sequence( size_t const& ntotal,
      		size_t const& n )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Utils:: subvector_sizes_sequence" ) ;

   size_t subsizemin = ntotal / n ;
   double res = double(ntotal) / double(n) - subsizemin ;
   size_t_vector sizes_sequence( n );
   
   // Exact even decomposition
   if ( fabs(res) < 1.e-12 )
     for (size_t i=0;i<n;++i) sizes_sequence(i) = subsizemin ;
   else
   {    
     size_t nmin = size_t( n * ( 1. - res + 1.e-12 ) ); // add 1.e-12 to avoid
     	// rounding errors
     size_t sum = 0 ;
     for (size_t i=0;i<nmin;++i) 
     {
       sizes_sequence(i) = subsizemin ;
       sum += sizes_sequence(i) ;
     }
     for (size_t i=nmin;i<n-1;++i) 
     {
       sizes_sequence(i) = subsizemin + 1 ;
       sum += sizes_sequence(i) ;
     }     
     sizes_sequence(n-1) = ntotal - sum ; 
   }
   
   return sizes_sequence;   
} 




//----------------------------------------------------------------------
void
FV_Mesh:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "FV_Mesh" ) ;
   
   writer->add_entry( "translation_active", MAC_Bool::create( 0, 
   	translation_projection ) ) ;
   writer->add_entry( "translation_direction", MAC_Int::create( 0, 
   	translation_direction ) ) ;
   writer->add_entry( "translation_magnitude", MAC_Double::create( 0, 
   	translation_magnitude ) ) ;	
   writer->add_entry( "translation_distance", 
   	MAC_Double::create( 0, translation_distance ) ) ;
   
   writer->finalize_object() ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
FV_Mesh:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "FV_Mesh" ) ;

   // Retrieving stored datas :
   bool read_translation_active = 
   	reader->data_of_entry( "translation_active" )->to_bool() ;
   int read_translation_direction = 
   	reader->data_of_entry( "translation_direction" )->to_int() ;
   double read_translation_magnitude = 
   	reader->data_of_entry( "translation_magnitude" )->to_double() ;

   // Does some checks
   MAC_ASSERT( read_translation_active == translation_projection ) ;
   MAC_ASSERT( read_translation_direction == (int) translation_direction ) ;
   MAC_ASSERT( fabs( read_translation_magnitude - translation_magnitude )
   	< 1e-12 ) ;

   translation_distance = 
   	reader->data_of_entry( "translation_distance" )->to_double() ;
	
   if ( fabs( translation_distance ) > 1e-12 )
   {
     if ( MAC_Exec::communicator()->rank() == 0 ) 
     {
        FV::out() << "*** FV_Mesh translation-projection" << endl << endl;
        FV::out() << "      Initial translation distance = " << 
		translation_distance << endl << endl;
     }             
     size_t ng = (*global_main_coordinates)[translation_direction].size() ;
     for (size_t i=0;i<ng;++i)
       (*global_main_coordinates)[translation_direction](i) += 
     		translation_distance ;
		
     size_t nl = (*local_main_coordinates)[translation_direction].size() ;
     for (size_t i=0;i<nl;++i)
       (*local_main_coordinates)[translation_direction](i) += 
     		translation_distance ;
		
     min_coordinate_on_current_processor(translation_direction) += 
     		translation_distance ;	
     max_coordinate_on_current_processor(translation_direction) += 
     		translation_distance ;						
     min_coordinate_with_halozone(translation_direction) += 
     		translation_distance ;	     
     max_coordinate_with_halozone(translation_direction) += 
     		translation_distance ;	   
   }

   reader->end_object_retrieval() ;

   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//----------------------------------------------------------------------
boolVector const*
FV_Mesh:: get_periodic_directions( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_periodic_directions" ) ;
   MAC_CHECK_PRE( periodic != NULL ) ;

   boolVector const* result = periodic ;
   
   return result ;
}




//----------------------------------------------------------------------
int const*
FV_Mesh:: get_MPI_coordinates( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_MPI_coordinates" ) ;

   int const* result = MPI_coordinates ;
   
   return result ;
}





//----------------------------------------------------------------------
int const*
FV_Mesh:: get_domain_decomposition( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_domain_decomposition" ) ;
      
   int const* result = domain_decomposition ;
   
   return result ;
}




//----------------------------------------------------------------------
list<size_t> const*
FV_Mesh:: get_MPI_neighbors_ranks( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_MPI_neighbors_ranks" ) ;

   list<size_t> const* result = &MPI_neighbors_World ;
   
   return result ;
}




//----------------------------------------------------------------------
list<size_t> const*
FV_Mesh:: get_MPI_periodic_neighbors_ranks( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_MPI_periodic_neighbors_ranks" ) ;

   list<size_t> const* result = &MPI_periodic_neighbors_World ;
   
   return result ;
}




//----------------------------------------------------------------------
double
FV_Mesh:: get_main_boundary_measure( string const& boundary_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_main_boundary_measure" ) ;
   MAC_CHECK( FV_DomainBuilder:: is_main_color( 
   	FV_DomainBuilder:: get_color_number( boundary_name ) ) ) ;

   double result = 0. ;
   size_t boundary_color = 
   	FV_DomainBuilder:: get_color_number( boundary_name ) ;
   bool is_3D = ( global_main_coordinates->size() == 3 );
   
   switch ( boundary_color )
   {
     case FV_BC_LEFT:
       result = ( main_domain_max_coordinates(1) -
       			main_domain_min_coordinates(1) );
       if ( is_3D )
         result *= ( main_domain_max_coordinates(2) -
       			main_domain_min_coordinates(2) );
       break;
       
     case FV_BC_RIGHT:
       result = ( main_domain_max_coordinates(1) -
       			main_domain_min_coordinates(1) );
       if ( is_3D )
         result *= ( main_domain_max_coordinates(2) -
       			main_domain_min_coordinates(2) );
       break;
       
     case FV_BC_BOTTOM:
       result = ( main_domain_max_coordinates(0) -
       			main_domain_min_coordinates(0) );
       if ( is_3D )
         result *= ( main_domain_max_coordinates(2) -
       			main_domain_min_coordinates(2) );
       break;
       
     case FV_BC_TOP:
       result = ( main_domain_max_coordinates(0) -
       			main_domain_min_coordinates(0) );
       if ( is_3D )
         result *= ( main_domain_max_coordinates(2) -
       			main_domain_min_coordinates(2) );
       break;
       
     case FV_BC_BEHIND:
       result = ( main_domain_max_coordinates(1) -
       			main_domain_min_coordinates(1) )
		* ( main_domain_max_coordinates(0) -
       			main_domain_min_coordinates(0) );
       break;
       
     case FV_BC_FRONT:
       result = ( main_domain_max_coordinates(1) -
       			main_domain_min_coordinates(1) )
		* ( main_domain_max_coordinates(0) -
       			main_domain_min_coordinates(0) );
       break;                          
   }
     
   return result ;   
}




//internal--------------------------------------------------------------
void
FV_Mesh_ERROR:: n1( MAC_ModuleExplorer const* exp )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "empty meshing encountered." << endl << "MODULE FV_Mesh required";
   MAC_Error::object()->raise_module_error( exp, mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_Mesh_ERROR:: n2( MAC_ModuleExplorer const* exp )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "splitting strategy missing." << endl << 
   	"MODULE splitting_strategy required";
   MAC_Error::object()->raise_module_error( exp, mesg.str() ) ;
}




//internal--------------------------------------------------------------
void
FV_Mesh_ERROR:: n3( MAC_ModuleExplorer const* exp )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Periodic flow features are wrong. Possible reasons are:" << endl << 
   	"   * periodic flow direction is not a periodic direction" 
	<< endl <<
	"   * pressure drop and flow rate are both specified";  
   MAC_Error::object()->raise_module_error( exp, mesg.str() ) ;
}
