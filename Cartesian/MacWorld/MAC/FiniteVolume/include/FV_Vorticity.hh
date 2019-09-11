#ifndef FV_VORTICITY_HH
#define FV_VORTICITY_HH

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
class LA_Vector ;
class LA_DistVector ;
class LA_SeqVector ;
class LA_Scatter ;
class LA_Matrix ;
class FV_SystemNumbering ;


/** @brief The Class FV_Vorticity.

Server for the computation and post-processing of the vorticity field.

@author A. Wachs - Particulate flow project 2010-2012 */

class FV_Vorticity : public FV_OneStepIteration
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
      ~FV_Vorticity( void ) ;
     
      /** @brief Copy constructor */      
      FV_Vorticity( FV_Vorticity const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */        
      FV_Vorticity& operator=( FV_Vorticity const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object   
      @param dom mesh and fields
      @param exp to read the data file */                 
      FV_Vorticity( MAC_Object* a_owner, 
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */      
      FV_Vorticity( void ) ;

      /** @brief Create a clone
      @param a_owner the MAC-based object
      @param dom mesh and fields
      @param exp to read the data file */
      virtual FV_Vorticity* create_replica( 
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;
      //@}


   //-- Utilities
   
      /** @name Utilities */
      //@{
      /** @brief Prepare Computation of the vorticity field */
      void prepare_vorticity_computation( void ) ;
      
      /** @brief Compute the vorticity field in 2D */
      void compute_vorticity_2D( void ) ;
      
      /** @brief Compute the vorticity field in 3D */
      void compute_vorticity_3D( void ) ;      		         
      //@}

      
   private: //----------------------------------------------------------------
      
   //-- Class attributes

      static FV_Vorticity const* PROTOTYPE ;

   //-- Attributes
      
      bool b_at_each_TS;

      // Fields
      FV_DiscreteField const* UU;
      FV_DiscreteField* OM; 
      size_t velocity_level ;
      string velocity_name ;       

      // MPI
      size_t nb_procs;
      size_t my_rank;
      size_t is_master;
      MAC_Communicator const* macCOMM; 
      
      // Matrix system
      FV_SystemNumbering* OM_NUM ;
      LA_Vector * VEC_OM ;
      LA_SeqVector * OM_LOC ;
      LA_Matrix * MAT_OM ;       
} ;

#endif
