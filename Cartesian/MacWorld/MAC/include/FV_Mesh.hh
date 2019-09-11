#ifndef FV_MESH_HH
#define FV_MESH_HH

#include <mpi.h>
#include <MAC_Object.hh>
#include <MAC_Data.hh>
#include <doubleVector.hh>
#include <stringVector.hh>
#include <size_t_vector.hh>
#include <boolVector.hh>
#include <string>
#include <vector>
#include <list>
using std::vector;
using std::string;
using std::list;

class MAC_Module ;
class MAC_ModuleExplorer ;


/** @brief The Class FV_Mesh.

Structured rectangular or cubic mesh for Finite Volume FV scheme.

@author A. Wachs - Particulate flow project 2010-2012 */


class FV_Mesh : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Attributes

      static stringVector* directionName;  


   //-- Instance delivery and initialization
      
      static FV_Mesh* create( MAC_Object* a_owner,
       	MAC_ModuleExplorer const* exp,
	size_t dim ) ;

      
   //-- Characteristics

      // name of `self'
      std::string const& register_name( void ) const ;

      
   //-- Input - Output
      
      void print( std::ostream& os, size_t indent_width ) const ;

      
   //-- Post processing

      void write_grid( MAC_Module* base, bool parallel ) const;


   //-- Persistence   

      virtual void save_state( MAC_ObjectWriter* writer ) const ;

      virtual void restore_state( MAC_ObjectReader* reader ) ;


   //-- Cell balance
      
      size_t point_owner( MAC_Data const* formula, size_t const& nb_ranks ) 
      	const ;


   //-- Access
      
      size_t nb_space_dimensions( void ) const ; 
      
      vector< doubleVector > const* get_global_main_coordinates( void ) const;
      
      size_t_vector const* get_global_max_index( void ) const;
      
      size_t_vector const* get_global_min_index_in_domain( void ) const;
      
      size_t_vector const* get_global_max_index_in_domain( void ) const;

      vector< doubleVector > const* get_local_main_coordinates( void ) const ;
      
      size_t_vector const* get_local_max_index_in_global( void ) const ;
      
      size_t_vector const* get_local_min_index_in_global( void ) const ;  
      
      size_t_vector const* get_local_max_index_in_global_on_current_proc( 
      	void ) const ;
      
      size_t_vector const* get_local_min_index_in_global_on_current_proc( 
      	void ) const ;  

      vector< boolVector > const* get_nodes_owner( void ) const ; 
      
      size_t get_local_nb_points( size_t direction ) const ; 
      
      size_t get_local_nb_points_on_current_proc( size_t direction ) const ; 
      
      size_t get_security_bandwidth( void ) const ; 
      
      double get_smallest_grid_size( void ) const ;
      
      double get_smallest_constant_grid_size( void ) const ;      
      
      size_t get_translation_direction( void ) const ;
      
      double get_translation_magnitude( void ) const ;
      
      double get_translation_distance( void ) const ;
      
      boolVector const* get_periodic_directions( void ) const ; 
      
      size_t get_periodic_flow_direction( void ) const ;

      double get_periodic_pressure_drop( void ) const ;

      double get_periodic_flow_rate( void ) const ;
      
      bool is_periodic_pressure_drop( void ) const ;
      
      bool is_periodic_flow_rate( void ) const ; 
      
      bool is_periodic_flow( void ) const ;
      
      bool is_periodic_domain( void ) const ;        
      
      int const* get_MPI_coordinates( void ) const ;

      int const* get_domain_decomposition( void ) const ;

      list<size_t> const* get_MPI_neighbors_ranks( void ) const ;
      
      list<size_t> const* get_MPI_periodic_neighbors_ranks( void ) const ;

      
   //-- Geometry
   
      bool is_in_domain_on_current_processor( double const& x, 
      	double const& y ) const ;
	
      bool is_in_domain_on_current_processor( double const& x, 
      	double const& y, double const& z ) const ;
	
      bool is_in_domain_on_current_processor( double const& coor, 
      	size_t direction, double const& tol ) const ;	
	
      bool is_in_domain_with_halozone( double const& x, double const& y )
      	const ;
	
      bool is_in_domain_with_halozone( double const& x, double const& y,
      	double const& z ) const ;
	
      bool is_in_domain_with_halozone_plus_ext( double const& x, 
      	double const& y, double const& ext ) const ;
	
      bool is_in_domain_with_halozone_plus_ext( double const& x, 
      	double const& y, double const& z, double const& ext ) const ;
	
      bool is_in_main_domain( double const& coor, size_t direction,
      	double const& tol ) const ;	
	
      double get_min_coordinate_on_current_processor( size_t direction ) 
      	const ;
	
      double get_min_coordinate_with_halozone( size_t direction ) const ;
      
      double get_max_coordinate_on_current_processor( size_t direction ) 
      	const ;
	
      double get_max_coordinate_with_halozone( size_t direction ) const ;
      
      double get_main_domain_min_coordinate( size_t direction ) const ;
      
      double get_main_domain_max_coordinate( size_t direction ) const ; 
      
      double get_main_domain_boundary_perp_to_direction_measure( 
      	size_t direction ) const ; 
	
      double get_main_boundary_measure( string const& boundary_name ) const ;


   //-- Translation
   
      void translation( void ) ;
      
      bool is_translation_active( void ) const ;


   //-- Modification
   
      void set_periodic_pressure_drop( double const& ppd ) ;
   
   
   //-- Utilities
   
      /** @brief Given a sorted vector of double, returns whether the input
      value lies in the interval defined by the vector and if yes, the index i0
      which implies that the input value is in the interval x[i0]:x[i0+1]
      @param mesh the sorted vector (generally a 1D mesh)
      @param x the input value
      @param i0 the index */        
      static bool between( doubleVector const* mesh, double const& x, 
       	size_t& i0 ) ;
	
      /** @brief Given a sorted vector of double, returns whether the input
      value lies in the sub-interval defined by mesh[i0_init] and mesh[i1_init]
      and if yes, the index i0 which implies that the input value is in the 
      interval x[i0]:x[i0+1]
      @param mesh the sorted vector (generally a 1D mesh)
      @param x the input value
      @param i0 the index 
      @param i0_init the sub-interval lower index 
      @param i1_init the sub-interval upper index */        
      static bool between_subinterval( doubleVector const* mesh,
       		double const& x, size_t& i0,
		const size_t& i0_init, const size_t& i1_init ) ;
		
      /** @brief Given a sorted vector of double, returns the smallest
      index imax such that mesh[imax] > x; in case x is larger than the largest
      value of mesh, the method returns the largest index of mesh 
      @param mesh the sorted vector (generally a 1D mesh)
      @param x the input value */        
      static size_t max_index( doubleVector const* mesh, double const& x ) ;
       
      /** @brief Given a sorted vector of double, returns the largest
      index imin such that mesh[imin] < x; in case x is smaller than the
      smallest value of mesh, the method returns 0, the smallest index of mesh 
      @param mesh the sorted vector (generally a 1D mesh)
      @param x the input value */        
      static size_t min_index( doubleVector const* mesh, double const& x ) ;
	
      /** brief Divide a vector in n sub-vectors whose size varies by 1 at most
      (takes the vector size as entry and returns the sequence of subvector 
      sizes)
      @param ntotal size of the total vector
      @param n number of subvectors */
      static size_t_vector subvector_sizes_sequence( size_t const& ntotal,
      		size_t const& n) ; 
		
      		            
            
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      FV_Mesh( void ) ;
     ~FV_Mesh( void ) ;
      FV_Mesh( FV_Mesh const& other ) ;
      FV_Mesh& operator=( FV_Mesh const& other ) ;

      FV_Mesh( MAC_Object* a_owner,
       	MAC_ModuleExplorer const* exp,
	size_t dim ) ;
      
   //-- Attributes
   
      vector< doubleVector >* global_main_coordinates; 
      size_t global_number_of_cells;
      size_t_vector global_max_index;
      size_t_vector global_min_index_in_domain;
      size_t_vector global_max_index_in_domain;            
      vector< doubleVector >* local_main_coordinates;
      vector< boolVector >* on_current_processor;       
      size_t local_number_of_cells;
      size_t_vector local_max_index_in_global;
      size_t_vector local_min_index_in_global;       
      size_t_vector local_max_index_in_global_on_current_proc;
      size_t_vector local_min_index_in_global_on_current_proc;       
      size_t security_bandwidth;
      doubleVector min_coordinate_on_current_processor;
      doubleVector max_coordinate_on_current_processor;
      doubleVector min_coordinate_with_halozone;
      doubleVector max_coordinate_with_halozone;
      doubleVector main_domain_min_coordinates;
      doubleVector main_domain_max_coordinates;      
      double smallest_grid_size ;
      size_t translation_direction;
      double translation_magnitude;
      double translation_distance;
      bool translation_projection; 
      boolVector* periodic; 
      size_t periodic_flow_direction ;
      double* periodic_pressure_drop ;
      double* periodic_flow_rate ;
      MPI_Comm* structured_cartesian_Comm ;
      int* MPI_coordinates ; 
      list<size_t> MPI_neighbors_World ;
      list<size_t> MPI_periodic_neighbors_World ; 
      int* domain_decomposition ;
} ;

// Inline methods
//----------------------------------------------------------------------
inline vector< doubleVector > const* 
FV_Mesh:: get_local_main_coordinates( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_Mesh:: get_local_main_coordinates" ) ;
   MAC_CHECK( local_main_coordinates != 0 ) ;

   vector< doubleVector > const* result = local_main_coordinates;

   return( result ) ;
}  

#endif

