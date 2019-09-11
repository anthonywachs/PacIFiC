#ifndef FV_BOUNDARY_CONDITION_HH
#define FV_BOUNDARY_CONDITION_HH

#include <MAC_Object.hh>
#include <FV_DiscreteField.hh>
#include <string>
#include <vector>
#include <list>
using std::list ;
using std::string ;
using std::pair ;

class MAC_DataWithContext ;
class stringVector ;
class boolVector ;
class MAC_Context ;


/** @brief The Class FV_BoundaryCondition.

Define a boundary condition for a scalar or vectorial field

@author A. Wachs - Particulate flow project 2010-2012 */

class FV_BoundaryCondition : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FV_BoundaryCondition* create( MAC_Object* a_owner, 
 		size_t number_of_components, size_t color ) ; 
		
      FV_BoundaryCondition* create_clone( MAC_Object* a_owner ) const ; 
		
		
   //-- Modification
   
      void add_MacTriplet( size_t component, FV_TRIPLET const& mt ) ; 
      
      void add_MacTriplet( size_t component, size_t i, size_t j, size_t k ) ;
	
      void read_dirichlet_BC( MAC_ModuleExplorer* sse,
      	size_t nb_space_dimensions, string const& field_name ) ;
      
      void set_imposed_DOF_values( FV_DiscreteField* ff ) ; 
      
      double get_imposed_DOF_values( size_t & component ) const;
      
      void set_free_DOF_values( FV_DiscreteField* ff, size_t level ) ; 
      
      void set_free_DOF_values( FV_DiscreteField* ff, size_t component, 
      	size_t level ) ;            
      
      void set_unknown_on_BC( size_t component ) ; 
      
      void set_periodic( void ) ;
      
      void set_none( void ) ;            
      
      void set_shift_MacTriplet( size_t component, 
      	int i, int j, int k, double coefficient = 1. ) ;
	
      void set_shift_MacTriplet( size_t component, 
      	FV_SHIFT_TRIPLET const* mst, double coefficient = 1. ) ;	
	
      void add_shift_MacTriplet( size_t component, 
      	int i, int j, int k, double coefficient = 1. ) ;
	
      void add_shift_MacTriplet( size_t component, 
      	FV_SHIFT_TRIPLET const* mst, double coefficient = 1. ) ;	


   //-- Input - Output
      
      void print( std::ostream& os, size_t indent_width,
      	size_t nb_space_dimensions ) const;

	
   //-- Access to data
   
      bool is_dirichlet( size_t component ) const ;
      
      bool is_neumann( size_t component ) const ;
      
      bool is_periodic( size_t component ) const ; 
      
      bool is_none( size_t component ) const ;            
      
      size_t get_color_ID( void ) const ; 
      
      bool has_unknown( size_t component ) const ;
      
      bool has_DOF_on_proc( void ) const ;  
      
         
   //-- Utilities
   
      double compute_boundary_cell_centered_DOF_integral( 
      	FV_DiscreteField const* ff, size_t component, size_t level ) const ;
	
      double compute_boundary_mean_normal_derivative( 
      	FV_DiscreteField const* ff, size_t component, size_t level ) const ;	 
            
   
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      FV_BoundaryCondition( void ) ;
     ~FV_BoundaryCondition( void ) ;
      FV_BoundaryCondition( FV_BoundaryCondition const& other ) ;
      FV_BoundaryCondition& operator=( FV_BoundaryCondition const& other ) ;
     
      FV_BoundaryCondition( MAC_Object* a_owner, 
 		size_t number_of_components, size_t color ) ;

      void set_BC_type( string const& type_of_boundary_condition, 
      	size_t component ) ;
	
   //-- Attributes

      vector< list<FV_TRIPLET> > nodes_localStructuredNumbering;
      stringVector* BC_TYPE ;
      boolVector* UNKNOWN_ON_BC ;
      vector< list< pair<double,FV_SHIFT_TRIPLET> > > SHIFT_UPDATE_BC_VALUES ;
      vector< MAC_DataWithContext const* > BC_VALUES_PER_COMP;
      MAC_DataWithContext const* BC_VALUES_ALL_COMPS;
      size_t NB_COMPS ;
      size_t BC_COLOR ;
      MAC_Context* CTX ; 
      MAC_DoubleVector* COORDS ;
      boolVector* READ; 
} ;

#endif
