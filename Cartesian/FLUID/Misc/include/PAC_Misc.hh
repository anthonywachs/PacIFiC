#ifndef PAC_MISC_HH
#define PAC_MISC_HH

#include <MAC.hh>
#include <iostream>
#include <fstream>
#include <string>
using std::cout;
using std::endl;
using std::string;

class FV_DiscreteField;


/** @brief The Class PAC_Misc.

Miscellaneous routines for all solvers.

@author A. Wachs - Pacific project 2018-2019 */

class PAC_Misc
{
  public: //-----------------------------------------------------------------

   //-- Utilities

      /** @name Utilities */
      //@{      
      /** @brief Compute flow rate
      @param FF field
      @param level field level      
      @param boundary_name boundary name */
      static double compute_flow_rate( FV_DiscreteField const* FF,
      	size_t level, string const& boundary_name );
	
      /** @ brief Clear results and restart directories in case of a new run */
      static void clearAllFiles( string const& resultsDirectory,
      	string const& savingsDirectory,
	const size_t& process_rank ) ;	
      //@}
      	

  protected: //--------------------------------------------------------------


  private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */          
      ~PAC_Misc( void ) {}

      /** @brief Copy constructor */      
      PAC_Misc( PAC_Misc const& other ) {}

      /** @brief Constructor without argument */      
      PAC_Misc( void ) {}
      //@}

} ;

#endif
