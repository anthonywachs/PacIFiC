#ifndef FV_DISCRETE_FIELD_HH
#define FV_DISCRETE_FIELD_HH

#include <MAC_Module.hh>
#include <MAC_Object.hh>
#include <string>
#include <boolVector.hh>
#include <size_t_array2D.hh>
#include <doubleArray3D.hh>
#include <intArray3D.hh>
#include <boolArray3D.hh>
#include <doubleVector.hh>
#include <vector>
#include <list>
#include <geomVector.hh>

using std::vector ;
using std::list ;
using std::pair ;

class MAC_ModuleExplorer ;
class stringVector ;
class LA_SeqVector ;
class FV_Mesh ;
class FV_BoundaryCondition ;
class FV_DomainBuilder ;
class FV_TimeIterator ;
class MAC_ObjectRegister ;
class MAC_Context ;
class MAC_DoubleVector ;
class MAC_DataWithContext ;
class MAC_Context ;
class LA_Matrix ;
class LA_Vector ;


/** @brief The Structure FV_TRIPLET.
FV Structured numbering
@author A. Wachs - Particulate flow project 2010-2012 */
struct FV_TRIPLET
{
   size_t i;
   size_t j;
   size_t k;
   
   size_t& index( size_t direction )
   {
     size_t* pIndex = 0 ;
     switch( direction )
     {
       case 0: pIndex = &i; break;
       case 1: pIndex = &j; break;
       case 2: pIndex = &k; break;
     }
     return ( *pIndex ) ;
   } ;
   
   size_t get_index( size_t direction ) const
   {
     size_t index_ = 0;
     switch( direction )
     {
       case 0: index_=i; break;
       case 1: index_=j; break;
       case 2: index_=k; break;
     }
     return ( index_ ) ;
   } ;  
};


/** @brief The Structure FV_SHIFT_TRIPLET.
FV Structured numbering shift
@author A. Wachs - Particulate flow project 2010-2012 */
struct FV_SHIFT_TRIPLET
{
   int i;
   int j;
   int k;
   
   int& index( size_t direction )
   {
     int* pIndex = 0 ;
     switch( direction )
     {
       case 0: pIndex = &i; break;
       case 1: pIndex = &j; break;
       case 2: pIndex = &k; break;
     }
     return ( *pIndex ) ;
   } ;
   
   int get_index( size_t direction ) const
   {
     int index_ = 0;
     switch( direction )
     {
       case 0: index_=i; break;
       case 1: index_=j; break;
       case 2: index_=k; break;
     }
     return ( index_ ) ;
   } ;    
};


/** @brief The Class FV_DiscreteField.

Discrete field of the FV/Finite Volume type, either centered or staggered.

@author A. Wachs - Particulate flow project 2010-2012 */

class FV_DiscreteField : public MAC_Object
{
   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization

      static FV_DiscreteField* create( MAC_Object* a_owner,
      		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type,
		size_t a_nb_components,
		size_t a_depth ) ;
		
      FV_DiscreteField* create_clone( MAC_Object* a_owner,
      		std::string const& name_of_new_field ) const ;		

      static size_t nb_objects( void ) ;
      
      virtual void build_BCs( MAC_ModuleExplorer const* exp, 
      	FV_DomainBuilder const* DB ) = 0 ;
	
      void build_field_numbering( void ) ;
      
      void initialize_DOFs( MAC_ModuleExplorer const* exp ) ;
      
      void allocate_DLMFDconstrainedDOFs( void ) ;	

      
   //-- Identification

      size_t id_number( void ) const ;

      std::string const& name( void ) const ;
      
      std::string const& discretization_type( void ) const ;
      
      std::string const& paraview_location( void ) const ;
      
      
   //-- DOFs status

      size_t nb_cells( void ) const ;

      size_t nb_components( void ) const ;

      size_t storage_depth( void ) const ;
      
      virtual bool is_global_triplet( int i, int j, int k,
      	size_t component ) const = 0 ;

      virtual bool is_global_triplet_local_DOF( size_t i, size_t j, size_t k,
      	size_t component ) const = 0 ;
	
      bool DOF_is_unknown( size_t i, size_t j, size_t k,
      	size_t component ) const ;

      virtual bool DOF_is_unknown_handled_by_proc( size_t i, size_t j, size_t k,
      	size_t component ) const = 0 ;	
      
      bool DOF_in_domain( int i, int j, int k, size_t component ) const ;
      
      bool DOF_on_proc( int i, int j, int k, size_t component ) const ;      	
	
      virtual bool DOF_has_imposed_Dirichlet_value( size_t i, size_t j, 
      	size_t k, size_t component ) const = 0 ;
	
      virtual bool DOF_on_BC( size_t i, size_t j, size_t k, size_t component ) 
      	const = 0 ; 		

      double DOF_value( size_t i, size_t j, size_t k,
      	size_t component, size_t level ) const ; 
	
      int DOF_color( size_t i, size_t j, size_t k,
      	size_t component ) const ; 
	
      virtual double get_DOF_coordinate( size_t i, size_t component, 
      	size_t direction ) const = 0 ;
	
      virtual double get_DOF_coordinate_Assembling( int i, size_t component, 
      	size_t direction ) const = 0 ;
		
      virtual double get_cell_size( size_t i, size_t component, 
      	size_t direction ) const = 0 ;
	
      virtual double get_face_perp_to_direction_measure( size_t i, size_t j, 
      	size_t k, size_t component, size_t direction ) const = 0 ;
	
      virtual double get_cell_measure( size_t i, size_t j, size_t k,
      	size_t component ) const = 0 ;
	
      virtual bool DOF_offset( int &i, int &j, int &k,
        size_t_vector center, size_t_vector stencil, vector<double> &offset, 
        size_t component ) const = 0 ;
	
      virtual doubleVector const* get_DOF_coordinates_vector( 
      	size_t component, size_t direction ) const = 0 ;
      
      void extract_unknown_DOFs_value( size_t level, 
	LA_SeqVector* vec ) const ;
	
      size_t DOF_local_number( size_t i, size_t j, size_t k,
      	size_t component ) const ;
	
      size_t DOF_global_number( size_t i, size_t j, size_t k,
      	size_t component ) const ;
	
      // Here DLMFDconstrainedDOFs are define as a Boolean
      // bool DOF_is_constrained( size_t i, size_t j, size_t k,
      // 	size_t component ) const ;
      
      // Here DLMFDconstrainedDOFs are define as a integer
      int DOF_is_constrained( size_t i, size_t j, size_t k,
      	size_t component ) const ;
	
      size_t global_index_from_local( size_t i, size_t component, 
      	size_t direction ) const ;


   //-- Access to other data

      FV_Mesh const* primary_grid( void ) const ; 
      
      bool all_BCs_nonDirichlet( size_t component ) const ; 
      
      list< pair<size_t,double> > main_geometric_boundaries_on_proc( void ) 
      	const ;
      
      /** @brief Interpolate field values in 3D
      @param X_coordinate x coordinate
      @param Y_coordinate y coordinate
      @param Z_coordinate z coordinate
      @param component field component
      @param level component storage level*/
      double interpolateFieldValues( const double &X_coordinate,
	const double &Y_coordinate, const double &Z_coordinate,
	size_t component, size_t level ) const ;
	
      /**
        @brief Reconstruct field values in 3D
        @param X_coordinate x coordinate
        @param Y_coordinate y coordinate
        @param Z_coordinate z coordinate
        @param component field component
        @param level component storage level
        @param kernelWidth length scale of Gaussian function
        @param dp Particle diameter
      */
      double GaussianReconstructionFieldValues(
        const double &X_coordinate,
        const double &Y_coordinate,
        const double &Z_coordinate,
        size_t component,
        size_t level,
        const double kernelWidth,
        const double dp ) const ;    
	
      /**
        @brief Reconstruct field values in 3D using a corrective Kernel
        @param X_coordinate x coordinate
        @param Y_coordinate y coordinate
        @param Z_coordinate z coordinate
        @param component field component
        @param level component storage level
        @param kernelWidth length scale of the corrective kernel
        @param dp Particle diameter
      */
      double CorrectiveKernelAverageFieldValues(
        const double &X_coordinate,
        const double &Y_coordinate,
        const double &Z_coordinate,
        size_t component,
        size_t level,
        const double kernelWidth,
        const double dp ) const ;
		
      /** @brief Interpolate field values in 2D
      @param X_coordinate x coordinate
      @param Y_coordinate y coordinate
      @param component field component
      @param level component storage level*/
      double interpolateFieldValues( const double &X_coordinate,
          const double &Y_coordinate, size_t component, size_t level ) const ;

      /** @brief Compute gradient in 3D, interpolating surounding nodes 
      field values for each face of the owner cell
      @param X_coordinate x coordinate
      @param Y_coordinate y coordinate
      @param Z_coordinate z coordinate
      @param component field component
      @param level component storage level*/
      geomVector interpolateGradient( const double &X_coordinate,
	const double &Y_coordinate, const double &Z_coordinate, 
	size_t component, size_t level ) const ;	
	
      /** @brief Compute laplacian in 3D, interpolating surounding nodes 
      field values for each face of the owner cell
      @param X_coordinate x coordinate
      @param Y_coordinate y coordinate
      @param Z_coordinate z coordinate
      @param component field component
      @param level component storage level*/
      double interpolateLap( const double &X_coordinate,
	const double &Y_coordinate, const double &Z_coordinate, 
	size_t component, size_t level ) const ;
	
      /**
        @brief Compute gradient in 3D, using a corrective Kernel
        adapted to the Cappecelatro method
        @param X_coordinate x coordinate
        @param Y_coordinate y coordinate
        @param Z_coordinate z coordinate
        @param component field component
        @param level component storage level
        @param kernelWidth length scale of Gaussian function
        @param dp Particle diameter
      */
      geomVector CorrectiveKernelGradientFieldValues(
        const double &X_coordinate,
        const double &Y_coordinate,
        const double &Z_coordinate,
        size_t component,
        size_t level,
        const double kernelWidth,
        const double dp ) const ;

      /**
        @brief Compute gradient in 3D, using two-cubes in both
	sides (not used for the moment)
        @param X_coordinate x coordinate
        @param Y_coordinate y coordinate
        @param Z_coordinate z coordinate
        @param component field component
        @param level component storage level
        @param kernelWidth length scale of Gaussian function
        @param dp Particle diameter
      */
      geomVector interpolateGradientWithKernel(
        const double &X_coordinate,
        const double &Y_coordinate,
        const double &Z_coordinate,
        size_t component,
        size_t level,
        const double kernelWidth,
        const double dp ) const ;
		
	
      /**
        @brief Compute Laplacian in 3D, equivalent of "interpolateGradient
	WithKernel" method, used to compute explicit viscouse drag
        @param X_coordinate x coordinate
        @param Y_coordinate y coordinate
        @param Z_coordinate z coordinate
        @param component field component
        @param level component storage level
        @param kernelWidth length scale of Gaussian function
        @param dp Particle diameter
      */
      double interpolateLaplacianWithKernel(
        const double &X_coordinate,
        const double &Y_coordinate,
        const double &Z_coordinate,
        size_t component,
        size_t level,
        const double kernelWidth,
        const double dp ) const ;
	
      /** @brief Compute gradient in 2D, interpolating surounding nodes 
      field values for each edge of the owner cell
      @param X_coordinate x coordinate
      @param Y_coordinate y coordinate
      @param component field component
      @param level component storage level*/
      geomVector interpolateGradient( const double &X_coordinate,
	const double &Y_coordinate, size_t component, size_t level ) const ;

      geomVector cellsize2D( const double &X_coordinate,
	const double &Y_coordinate, size_t component, size_t level ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width,
      	bool b_values = false ) const = 0;
	
      void out_endOfBuilding( std::ostream& os, size_t indent_width,
      	size_t rank = 0 ) const ;
	
      virtual void compute_normLinf( double time ) const = 0;

	
   //-- Modification
   
      void set_imposed_DOF_values( void ) ;

      void set_neumann_DOF_values( void ) ;
      
      void set_DOF_value( size_t i, size_t j, size_t k,
      	size_t component, size_t level, double val ) ;
      
      void set_DOFs_value( size_t component, size_t level, double value ) ;

      void add_value_to_DOF( size_t i, size_t j, size_t k,
      	size_t component, size_t level, double value ) ;

      void update_free_DOFs_value( size_t level,
	LA_SeqVector const* vec ) ;
	
      void add_to_free_DOFs_value( size_t level,
	LA_SeqVector const* vec ) ;	
	
      void set_BC_values_modif_status( bool const& allowed ) ;
      
      virtual void set_postprocessing_options( 
      		std::string const& a_location,
		std::string const& a_paraview_fname ) = 0;

      virtual void copy_DOFs_value( size_t source_level, size_t target_level ) ;
      
      virtual void add_to_DOFs_value( size_t component, size_t level,
      		double const& val ) ;      
      
      void initialize_DLMFDconstrainedDOFs( void ) ;

      // Here DLMFDconstrainedDOFs are define as a Boolean
      //void set_DOF_constrained( size_t i, size_t j, size_t k,
      //	size_t component ) ;      
  
      // Here DLMFDconstrainedDOFs are define as a Boolean
      void set_DOF_constrained( size_t i, size_t j, size_t k,
            size_t component, bool b_has_extended_constrained = 0 ) ;      
    

            
   //-- Post processing
   
      virtual void write_field( MAC_Module* point_data,
	MAC_Module* cell_data ) const = 0;
			       
      virtual double interpolate_values_at_nodes( 
      		size_t i, size_t shift_i,
		size_t j, size_t shift_j,
		size_t k, size_t shift_k,
		size_t component,
		size_t level ) const = 0;
			       

   //-- Persistence   

      virtual void save_state( MAC_ObjectWriter* writer ) const ;

      virtual void restore_state( MAC_ObjectReader* reader ) ;
      
      static void read_state_nonrestored( MAC_ObjectReader* reader ) ;      


   //-- Field & System numbering
   
      void build_system_numbering( size_t_vector* idx_locs, 
      	size_t_vector* idx_globs ) const ;
	
      size_t nb_global_unknowns( void ) const ;	

      size_t nb_local_unknowns( void ) const ;
      
      size_t nb_local_unknowns_handled_by_proc( void ) const ;
      
      virtual size_t get_min_index_unknown_handled_by_proc( size_t component, 
      	size_t direction ) const = 0 ;
	
      virtual size_t get_max_index_unknown_handled_by_proc( size_t component, 
      	size_t direction ) const = 0 ;
	
      virtual size_t get_min_index_unknown_on_proc( size_t component, 
      	size_t direction ) const = 0 ;
	
      virtual size_t get_max_index_unknown_on_proc( size_t component, 
      	size_t direction ) const = 0 ;	
	
      virtual size_t get_local_nb_dof( size_t component, 
      	size_t direction ) const = 0 ;		

      void check_field_numbering( void ) const ;
      
      	
   //-- Interpolation tools

      /* @brief Interpolate staggered field on centered Control Volume */
      virtual FV_SHIFT_TRIPLET	shift_staggeredToCentered( void ) const ;
      
      /* @brief Interpolate vertex field on staggered Control Volume */
      virtual FV_SHIFT_TRIPLET	shift_vertexToStaggered( void ) const ;      
			       
      /* @brief Interpolate one component of a staggered field
      	on another component of a staggered Control Volume */
      virtual FV_SHIFT_TRIPLET	shift_staggeredToStaggered( 
      			size_t component ) const ; 
			
      /* @brief Interpolate vorticity field on staggered Control Volume */
      virtual FV_SHIFT_TRIPLET	shift_vorticityToStaggered( 
      			size_t component ) const ; 
			
      /* @brief Interpolate staggered field on tensor Control Volume */
      virtual FV_SHIFT_TRIPLET	shift_staggeredToTensor( 
      			size_t component ) const ; 
			
      /* @brief Interpolate tensor field on tensor Control Volume */
      virtual FV_SHIFT_TRIPLET	shift_tensorToTensor( 
      			size_t component ) const ; 
			
      /* @brief Interpolate comp1 on comp2 in 2D */
      virtual double interpolateOneCompOnAnotherComp(
		size_t i, size_t j,
		size_t comp1, size_t comp2, size_t level,
		FV_SHIFT_TRIPLET shift ) const ;

      /* @brief Interpolate comp1 on comp2 in 3D */
      virtual double interpolateOneCompOnAnotherComp(
		size_t i, size_t j, size_t k,
		size_t comp1, size_t comp2, size_t level,
		FV_SHIFT_TRIPLET shift ) const ;


   //-- Translation-projection
   
      virtual void create_transproj_interpolation( void ) = 0;
      
      virtual void translation_projection( size_t const& level, 
      	size_t const& temporary_level, bool translate_mesh = true,
	doubleVector const* outOfDomain_values = NULL ) = 0 ; 
	
      virtual void restore_translated_field_mesh( void ) = 0 ; 
      
      void check_field_primary_meshes_coincide( bool const &force, 
      	std::ostream& os, size_t indent_width ) ; 
	
	            	     
   //-- Utilities
   
      virtual double compute_boundary_cell_centered_DOF_integral( 
      	size_t component, size_t level, std::string const& boundary_name ) 
	const ;

      virtual double compute_boundary_mean_normal_derivative( 
      	size_t component, size_t level, std::string const& boundary_name ) 
	const ;
      
      virtual double compute_grad_at_cell_center(
      		size_t i, size_t j,
		FV_SHIFT_TRIPLET shift,
		size_t component, size_t direction,
		double dxC, double dyC ) ;      
      
      virtual double compute_grad_at_cell_center(
      		size_t i, size_t j, size_t k,
		FV_SHIFT_TRIPLET shift,
		size_t component, size_t direction,
		double dxC, double dyC, double dzC ) ; 
		
      void synchronize_field( size_t level ) ;

      virtual void synchronize_extended_DOF_constrained( 
		size_t level, 
		list<FV_TRIPLET> &ijk, 
		list<size_t> &comps );  
      
      void set_HaloZone_Features( size_t level ) ;      

      void set_HaloZone_Features_Periodic( void ) ;
      
      FV_BoundaryCondition const* get_BC( size_t color ) const;      


   //-- Spatial discretization and matrix assembly
   
      void assemble_constantcoef_laplacian_matrix( double const& coef_lap,
      		LA_Matrix *MAT, LA_Vector *VEC_rhs,
		bool const& rescale = false ) const;
		
      void assemble_mass_matrix( double const& coef, 
      		LA_Matrix *MAT, double const& power_index = 1 ) const;	

      void assemble_mass_vector( double const& coef, 
      		LA_Vector *VEC, double const& power_index = 1 ) const;		

      virtual void assemble_pDivv_matrix( FV_DiscreteField const* UU,
	double const& coef, LA_Matrix *MAT, LA_Vector *VEC_rhs ) const;	
	
      virtual void assemble_tauGradv_tensor_divergence_matrix( 
      	FV_DiscreteField const* DD,
	double const& coef, LA_Matrix *MAT ) const;		
		
      /** @brief Compute the slope limiter SuperBee for the TVD scheme
      @param theta the slope */
      static double SuperBee_phi( double const& theta );
      
      virtual double compute_CFL( FV_TimeIterator const* t_it,
	size_t level ) const;
	
      virtual void assemble_advection_Upwind( 
      	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const;
	
      virtual void assemble_advection_TVD( 
      	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const;		


   protected: //--------------------------------------------------------

      virtual ~FV_DiscreteField( void ) ;
      
      FV_DiscreteField( MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type ) ;
      
      FV_DiscreteField( std::string const& a_type ) ;
      
      FV_DiscreteField( MAC_Object* a_owner,
      		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type,
		size_t a_nb_components,
		size_t a_depth ) ;      
      
      virtual FV_DiscreteField* create_clone_replica( 
      		MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type ) const = 0 ; 
		
      virtual FV_DiscreteField* create_replica( MAC_Object* a_owner,
		FV_Mesh const* a_primary_mesh,
		std::string const& a_name,
		std::string const& a_type,
		size_t a_nb_components,
		size_t a_depth ) const = 0 ; 		
		
      bool is_a_prototype( void ) const ;
      
      void read_BCs( MAC_ModuleExplorer const* exp, 
      	FV_DomainBuilder const* DB ) ;
		
      void set_DOF_colors( size_t component ) ;
	
      void set_PERIODIC_default( size_t component, 
      	size_t_vector const* periodic_depth = NULL ) ;	

      void set_DOF_status( size_t component ) ;  
      
      void print_BC_FieldValues( std::ostream& os, size_t indent_width,
      	bool b_values = false ) const ;      
            		
      void set_freeDOFonBC_update_features( void ) ;
      
      
   //-- Attributes

      FV_Mesh const* PRIMARY_GRID ;
      std::string const FNAME ; 
      std::string const FDISCRETIZATION ;      
      size_t const ID ;

      size_t DIM ;
      size_t NB_COMPS ;
      size_t NB_CELLS ;
      size_t STO_DEPTH ; 
       
      vector< vector< doubleArray3D > >* VALUES ;
      list< FV_BoundaryCondition* > SET_OF_BCS ;
      vector< FV_BoundaryCondition* >* V_BCS ;      
      
      vector< intArray3D >* UNK_LOCAL_NUMBERING ;
      vector< intArray3D >* UNK_GLOBAL_NUMBERING ;
      size_t NB_LOCAL_UNKNOWNS ;
      size_t NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC ;
      size_t NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_BUFFERZONE ; 
      size_t NB_LOCAL_UNKNOWNS_HANDLED_BY_PROC_IN_PERIODIC_BUFFERZONE ;
      size_t NB_LOCAL_DOF ; 
      size_t NB_GLOBAL_UNKNOWNS ; 
      
      MAC_DataWithContext const* INITIALIZER;
      MAC_Context* CTX ; 
      MAC_DoubleVector* COORDS ;

      vector< size_t_vector >* global_max_index;
      vector< size_t_vector >* local_max_index_in_global;
      vector< size_t_vector >* local_min_index_in_global;
      vector< size_t_vector >* local_dof_number;
      vector< size_t_vector >* max_index_unknown_handled_by_proc;
      vector< size_t_vector >* min_index_unknown_handled_by_proc;
      vector< size_t_vector >* max_index_unknown_on_proc;
      vector< size_t_vector >* min_index_unknown_on_proc;      
      
      
      vector< vector< doubleVector > >* global_main_coordinates;       
      vector< vector< doubleVector > >* local_main_coordinates;
      vector< vector< doubleVector > >* local_cell_size;            

      vector< vector< size_t_vector > >* on_current_processor; /* 0 = on proc,
      	1 = bufferzone, 2 = halozone, 3 = bufferzone periodic, 4 = halozone
	periodic, 5 = bufferzone standard + periodic */
      
      vector< intArray3D >* DOFcolors;
      vector< intArray3D >* DOFstatus; 
      //vector< boolArray3D >* DLMFDconstrainedDOFs;
      vector< intArray3D >* DLMFDconstrainedDOFs;     
      
      boolVector* PERIODIC;      
      vector< FV_SHIFT_TRIPLET >* PERIODIC_SHIFT;
            
      std::string PARAVIEW_FNAME ; 
      std::string LOCATION ;      
      size_t NB_DOF_POSTPROCESSING_PER_COMP ; 
      
      bool ALL_COMPS_SAME_LOCATION ;
      bool SET_BC_VALUES_ALLOWED ; 
      
      vector< size_t_vector >* transproj_interpolation ;      

      list< vector< double* > > halozone_values_list ;
      list< vector< double* > > bufferzone_values ; 
      size_t synchronized_level ;     
      
   private: //----------------------------------------------------------

      FV_DiscreteField( void ) ;
      FV_DiscreteField( FV_DiscreteField const& other ) ;
      FV_DiscreteField& operator=( FV_DiscreteField const& other ) ;

      static MAC_ObjectRegister* plugins_map( void ) ;    

      virtual void set_nb_dof_post_processing( void ) = 0 ;

      void set_min_max_indices_unknown_handled_by_proc( 
      		size_t i, size_t j, size_t k, size_t component ) ;
		
      void set_min_max_indices_unknown_on_proc( 
      		size_t i, size_t j, size_t k, size_t component ) ;		

      virtual void translate_field_mesh( const size_t& trans_dir, 
      	const double &trans_dist ) = 0 ;
	
      void build_periodic_numbering( void ) ;
      
      double one_DOF_laplacian( LA_Matrix *MAT, LA_Vector *VEC_rhs, 
      		bool const& rescale,
          	int i, int j, int k,
          	size_t component, bool center_unknown_on_BC,
          	size_t center_pos_in_matrix, double ai ) const;	

   //-- Class attributes

      static size_t NB_INSTANCES ;
      bool const IS_PROTO ;             
} ;

const size_t OUT_OF_DOMAIN_PROJTRANS = 88488848 ;

struct FV_DiscreteField_ERROR
{
   static void n1( std::string const& name ) ;
   static void n2( std::string const& name ) ;
   static void n3( std::string const& name ) ;
   static void n4( std::string const& name ) ;         
} ;

#endif

