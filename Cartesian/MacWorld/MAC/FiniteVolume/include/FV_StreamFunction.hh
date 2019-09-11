#ifndef FV_STREAMFUNCTION_HH
#define FV_STREAMFUNCTION_HH

#include <FV_OneStepIteration.hh>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using std::list ;
using std::pair ;
using std::string ;

class MY_Communicator ;
class FV_DiscreteField ;
class LA_Matrix ;
class LA_Vector ;
class LA_SeqVector ;
class LA_Scatter ;
class LA_Solver ;
class FV_SystemNumbering ;


/** @brief The Class FV_StreamFunction.

Server for the resolution of the 2D stream function equation 
-lap(sf) = vorticity.

@author A. Wachs - Particulate flow project 2010-2012 */

class FV_StreamFunction : public FV_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @name Substeps of the step by step progression */
      //@{
      /** @brief Tasks performed at initialization of the algorithm, before
      starting the time stepping loop
      @param t_it time iterator */
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename ) ;
      
      /** @brief Perform one time step
      @param t_it time iterator */      
      virtual void do_one_inner_iteration( FV_TimeIterator const* t_it ) ;
      
      /** @brief Tasks performed at initialization of each time step 
      @param t_it time iterator */       
      virtual void do_before_inner_iterations_stage( 
      	FV_TimeIterator const* t_it );
      
      /** @brief Tasks performed after of each time step 
      @param t_it time iterator */       
      virtual void do_after_inner_iterations_stage( 
      	FV_TimeIterator const* t_it );
      
      /** @brief Tasks performed at the end of the time stepping loop */      
      virtual void do_after_time_stepping( void );
      
      /** @brief Save additional data than fields 
      @param t_it time iterator 
      @param cycleNumber cycle number */      
      virtual void do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber );      	
      //@}
      
   //-- Additional post-processing

      virtual void do_more_post_processing( FV_DomainAndFields* dom,
      		MAC_ModuleExplorer const* exp ) ;
            
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */          
      ~FV_StreamFunction( void ) ;
     
      /** @brief Copy constructor */      
      FV_StreamFunction( FV_StreamFunction const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */        
      FV_StreamFunction& operator=( FV_StreamFunction const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object   
      @param dom mesh and fields
      @param exp to read the data file */                 
      FV_StreamFunction( MAC_Object* a_owner, 
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */      
      FV_StreamFunction( void ) ;

      /** @brief Create a clone
      @param a_owner the MAC-based object
      @param dom mesh and fields
      @param exp to read the data file */
      virtual FV_StreamFunction* create_replica( 
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;
      //@}


   //-- Basic discrete system building
   
      /** @name Basic discrete system building */
      //@{
      /** @brief Assemble laplacian */
      void assemble_StreamFunctionLaplacian_matrix( void ) ;

      /** @brief Assemble one item of laplacian matrix */
      double one_DOF_StreamFunctionLaplacian( 
	int i, int j, 
	size_t center_pos_in_matrix, 
	double ai ) ;
	
      /** @brief Assemble vorticity rhs */
      void assemble_StreamFunctionRHS( void ) ;			         
      //@}

      
   private: //----------------------------------------------------------------
      
   //-- Class attributes

      static FV_StreamFunction const* PROTOTYPE ;

   //-- Attributes

      // Fields
      FV_DiscreteField const* UU;
      FV_DiscreteField* SF; 
      size_t velocity_level ;
      string velocity_name ;
      string streamfunction_name ;             

      // MPI
      size_t nb_procs;
      size_t my_rank;
      size_t is_master;
      MAC_Communicator const* macCOMM; 
      
      // Matrix system
      FV_SystemNumbering* SF_NUM ;
      LA_Vector * VEC_SF ;
      LA_SeqVector * SF_LOC ; 
      LA_Matrix * MAT_SF_Laplacian ;
      LA_Vector * VEC_SF_rhs ; 
      LA_Solver* SOLVER_SF ;
      list< pair< size_t, double > > SF_DirichletBCValues ;
} ;

#endif
