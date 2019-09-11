#ifndef REG_HEAT_TRANSFER_HH
#define REG_HEAT_TRANSFER_HH

#include <MAC_Object.hh>
#include <FV_DiscreteField.hh>
#include <VPH_HeatTransferSystem.hh>
#include <utility>
using namespace std;


class MAC_ModuleExplorer ;
class MAC_ListIdentity ;
class MAC_TimeIterator ;
class FV_DiscreteField ;
class FV_DomainAndFields ;
class FV_TimeIterator ;
class MAC_Communicator ;


/** @brief The Structure NavierStokes2Temperature.

Input data to be transferred from Navier Stokes solver to heat transfer solver.

@author A. Wachs - Pacific project 2018-2019 */
struct NavierStokes2Temperature
{
  double density_ ;
  bool b_restart_ ;      
  FV_DomainAndFields const* dom_ ; 
  FV_DiscreteField const* UU_ ;
  size_t levelAdvectingVelocity_ ; 
  double imposed_CFL_ ;
  string resultsDirectory_ ; 
};



/** @brief The Class VPH_HeatTransfer.

Solver for the time-dependent advection diffusion heat transfer problem. 

@author A. Wachs - Pacific project 2018-2019 */

class VPH_HeatTransfer : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */      
      VPH_HeatTransfer( void ) ;

      /** @brief Destructor */       
      virtual ~VPH_HeatTransfer( void ) ;	

      /** @brief Copy constructor */       
      VPH_HeatTransfer( 
      	VPH_HeatTransfer const& other ) ;

      /** @brief Operator == 
      @param other the right hand side */   
      VPH_HeatTransfer& operator=( 
      	VPH_HeatTransfer const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param fromNS structure containing input data from the NS solver */
      VPH_HeatTransfer ( MAC_Object* a_owner,
          MAC_ModuleExplorer const* exp,
          struct NavierStokes2Temperature const& fromNS );      
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of VPH_HeatTransfer
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param fromNS structure containing input data from the NS solver */
      static VPH_HeatTransfer* create( MAC_Object* a_owner,
          MAC_ModuleExplorer const* exp,
          struct NavierStokes2Temperature const& fromNS ) ;


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
      
      /** @brief Save additional data than fields 
      @param t_it time iterator       
      @param cycleNumber cycle number */      
      virtual void do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber  );
	
      // Save other data than FV fields for restart 
      virtual void do_additional_save_for_restart( 
      	FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, string const& basename ) ;
      //@}      


   //-- Persistence

      // Add objects to be stored for restart
      virtual void add_storable_objects( MAC_ListIdentity* list ) ;
      
            
   protected: //--------------------------------------------------------


   private: //----------------------------------------------------------
   
   //-- Utilities

      /** @name Utilities */
      //@{		      
      /** @ brief Reload other data than FV fields for restart  */
      void do_additional_reload( string const& basename ) ;           
      //@}     
   

   private: //----------------------------------------------------------
      
   //-- Attributes

      // Fields
      FV_DiscreteField* TT ; 
      FV_DiscreteField const* UU ;

      VPH_HeatTransferSystem* GLOBAL_EQ ;

      // MPI parameters
      size_t nb_ranks;
      size_t my_rank;
      size_t is_master;
      MAC_Communicator const* macCOMM;
      
      // Parameters
      double density;
      double heat_capacity;
      double thermal_conductivity;

      // Numerical parameters
      double imposed_CFL ;
      string AdvectionScheme ;
      string resultsDirectory;      
      size_t DiffusionTimeAccuracy ;
      size_t AdvectionTimeAccuracy ;
      size_t levelAdvectingVelocity ;

      // Restart
      bool b_restart ;      
      
} ; 

#endif
