#ifndef FV_HH
#define FV_HH

#include <MAC.hh>
#include <iostream>


/** @brief The Class FV.

Utilities for FV/Finite Volume scheme.

@author A. Wachs - Particulate flow project 2011-2013 */

class FV
{
   public: //-----------------------------------------------------------------

   //-- Output methods

      /** @name Output methods */
      //@{         
      /** @brief Standard output */
      static std::ostream& out(	void ) ;	   
      //@} 
            

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */          
      ~FV( void ) {}
     
      /** @brief Copy constructor */      
      FV( FV const& other ) {}

      /** @brief Constructor without argument */      
      FV( void ) {}
      //@}      
            
} ;

#endif
