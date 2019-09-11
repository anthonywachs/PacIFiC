#ifndef FV_DISCRETE_FIELD_CENTERED_HH
#define FV_DISCRETE_FIELD_CENTERED_HH

#include <FV_DiscreteField.hh>
#include <string>
#include <boolVector.hh>
#include <size_t_vector.hh>
#include <doubleVector.hh>
#include <doubleArray3D.hh>
#include <vector>
using std::vector ;

class MAC_ModuleExplorer ;
class MAC_Module ;
class stringVector ;
class LA_SeqVector ;
class FV_Mesh ;
class intArray3D ;


/** @brief The Class FV_DiscreteField_Centered.

Discrete field of the centered FV/Finite Volume type.

@author A. Wachs - Particulate flow project 2010-2012 */

class FV_DiscreteField_Centered : public FV_DiscreteField
{
   public: //-----------------------------------------------------------

   //-- DOFs status

      virtual bool is_global_triplet( int i, int j, int k,
      	size_t component ) const ;

      virtual bool is_global_triplet_local_DOF( size_t i, size_t j, size_t k,
      	size_t comp ) const ;

      virtual bool DOF_is_unknown_handled_by_proc( size_t i, size_t j, size_t k,
      	size_t component ) const ;	

      virtual bool DOF_has_imposed_Dirichlet_value( size_t i, size_t j, 
      	size_t k, size_t component ) const ;	

      virtual bool DOF_on_BC( size_t i, size_t j, size_t k, 
      	size_t component ) const  ; 	

      virtual double get_DOF_coordinate( size_t i, size_t component, 
      	size_t direction ) const ;

      virtual double get_DOF_coordinate_Assembling( int i, size_t component, 
      	size_t direction ) const ;
	
      virtual double get_cell_size( size_t i, size_t component, 
      	size_t direction ) const ;
	
      virtual double get_face_perp_to_direction_measure( size_t i, size_t j, 
      	size_t k, size_t component, size_t direction ) const ;
	
      virtual double get_cell_measure( size_t i, size_t j, size_t k,
      	size_t component ) const ;
	
      virtual bool DOF_offset( int &i, int &j, int &k,
        size_t_vector center, size_t_vector stencil, vector<double> &offset, 
        size_t component ) const ;	

      virtual doubleVector const* get_DOF_coordinates_vector( size_t component, 
      	size_t direction ) const ;	
		
	
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width,
      	bool b_values = false ) const ;

      virtual void compute_normLinf( double time ) const ;


   //-- Modification
      virtual void set_postprocessing_options( 
      		std::string const& a_location,
		std::string const& a_paraview_fname );

      virtual void copy_DOFs_value( size_t source_level, size_t target_level ) ;
      
      virtual void add_to_DOFs_value( size_t component, size_t level,
            double const& val ) ;

   //-- Initialization
   
      virtual void build_BCs( MAC_ModuleExplorer const* exp, 
      	FV_DomainBuilder const* DB ) ;      

      
   //-- Post processing

      virtual void write_field(MAC_Module* point_data,
	MAC_Module* cell_data) const ;

      virtual double interpolate_values_at_nodes( 
      		size_t i, size_t shift_i,
		size_t j, size_t shift_j,
		size_t k, size_t shift_k,
		size_t component,
               	size_t level ) const ;


   //-- Field & System numbering
   
      size_t get_min_index_unknown_handled_by_proc( size_t component, 
      	size_t direction ) const ;
				
      size_t get_max_index_unknown_handled_by_proc( size_t component, 
      	size_t direction ) const ;
	
      size_t get_min_index_unknown_on_proc( size_t component, 
      	size_t direction ) const ;
	
      size_t get_max_index_unknown_on_proc( size_t component, 
      	size_t direction ) const ;	
	
      size_t get_local_nb_dof( size_t component, 
      	size_t direction ) const ;
	
      /* @brief Interpolate staggered field on centered Control Volume */
      virtual FV_SHIFT_TRIPLET	shift_staggeredToCentered( void ) const ;

	
   //-- Translation-projection
   
      virtual void create_transproj_interpolation( void ) ;
      
      virtual void translation_projection( size_t const& level, 
      	size_t const& temporary_level, bool translate_mesh = true,
	doubleVector const* outOfDomain_values = NULL ) ; 
	
      virtual void restore_translated_field_mesh( void ) ;


   //-- Utilities
   
      virtual double compute_boundary_cell_centered_DOF_integral( 
      	size_t component, size_t level, std::string const& boundary_name ) 
	const ;


   //-- Spatial discretization and matrix assembly 
   
      void assemble_pDivv_matrix( FV_DiscreteField const* UU,
	double const& coef, LA_Matrix *MAT, LA_Vector *VEC_rhs ) const;
	
      void assemble_advection_Upwind( 
      	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const;
	
      void assemble_advection_TVD( 
      	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const;
        
                	
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      FV_DiscreteField_Centered( void ) ;
     ~FV_DiscreteField_Centered( void ) ;
      FV_DiscreteField_Centered( FV_DiscreteField_Centered const& other ) ;
      FV_DiscreteField_Centered& operator=( 
      	FV_DiscreteField_Centered const& other ) ;

      FV_DiscreteField_Centered( MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type,
		size_t a_nb_components,
		size_t a_depth ) ;
		
      FV_DiscreteField_Centered( std::string const& a_type ) ;	
      
      FV_DiscreteField_Centered( MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type ) ;	      
      
      FV_DiscreteField* create_replica( MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type,
		size_t a_nb_components,
		size_t a_depth ) const; 

      FV_DiscreteField* create_clone_replica( 
      		MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type ) const ; 
   

      virtual void set_nb_dof_post_processing( void ) ;
      
      bool check_invariant_cell_features( size_t i, size_t direction  ) const ;

      void translate_field_mesh( const size_t& trans_dir, 
      	const double &trans_dist ) ;
	
      void one_DOF_pDivv( FV_DiscreteField const* UU,
	LA_Matrix *MAT, LA_Vector *VEC_rhs, 
	int i, int j, int k, size_t component, 
	size_t center_pos_in_matrix, double ai ) const;

		   		     	      
   //-- Attributes

      static FV_DiscreteField_Centered const* PROTOTYPE ;
            
} ;


#endif

