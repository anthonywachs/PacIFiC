#ifndef _FS_INTERFACE_TEST__
#define _FS_INTERFACE_TEST__

#include <FV_OneStepIteration.hh>
#include <FV_DiscreteField.hh>
#include <string>
#include <sstream>
using std::string;
using std::istringstream ;
class MAC_Communicator ;
class MAC_ListIdentity ;
class FS_SolidPlugIn ;
class DS_AllRigidBodies ;


/** @brief The Class FS_Interface_Test.

Test of the interface between the fluid solver and Grains3D. 

@author A. Wachs - Pacific project 2021 */

class FS_Interface_Test : public FV_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @name Substeps of the step by step progression */
      //@{
      /** @brief Tasks performed at initialization of the algorithm, before
      starting the time stepping loop
      @param t_it time iterator */
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it, 
      	string const& basename ) ;
      
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
      //@}      
            

   protected: //--------------------------------------------------------------


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */          
      ~FS_Interface_Test( void ) ;
     
      /** @brief Copy constructor */      
      FS_Interface_Test( FS_Interface_Test const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */        
      FS_Interface_Test& operator=( 
      	FS_Interface_Test const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the PEL-based object   
      @param dom mesh and fields
      @param exp to read the data file */                 
      FS_Interface_Test( MAC_Object* a_owner, 
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */      
      FS_Interface_Test( void ) ;

      /** @brief Create a clone
      @param a_owner the PEL-based object
      @param dom mesh and fields
      @param exp to read the data file */
      virtual FS_Interface_Test* create_replica( 
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;
      //@}


   private: //----------------------------------------------------------------
      
   //-- Class attributes

      static FS_Interface_Test const* PROTOTYPE ;

   //-- Attributes   

      FV_DiscreteField* PP;
      FV_DiscreteField* UU;      
      size_t dimension;
       
      // MPI parameters
      size_t nb_ranks;
      size_t my_rank;
      size_t is_master;
      MAC_Communicator const* macCOMM; 
      
      // Fluid-solid interface
      string solidSolverType;
      FS_SolidPlugIn* solidSolver;
      bool b_solidSolver_parallel;
      string solidSolver_insertionFile;
      string solidSolver_simulationFile;
      istringstream* solidFluid_transferStream; 
      DS_AllRigidBodies* allrigidbodies; 
      bool b_particles_as_fixed_obstacles;  
                          
} ;

#endif
