#ifndef _DS_3DCYLINDER__
#define _DS_3DCYLINDER__

#include <DS_RigidBody.hh>
#include <string>
using std::string;


/** @brief The class DS_3Dcylinder.

A moving or stationary rigid 3D cylinder of axisymmetric cross-section in the 
Direction Splitting solver. 
 
@author A. Wachs - Pacific project 2021 */

class DS_3Dcylinder: public DS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_3Dcylinder(); 
      
      /** @brief Constructor with arguments 
      @param pgrb pointer to the corresponding geometric rigid body */
      DS_3Dcylinder( FS_RigidBody* pgrb );      	   

      /** @brief Destructor */
      ~DS_3Dcylinder();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates the sphere features from its corresponding 
      geometric rigid body */
      void update();     
      //@}
      

   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the attributes in a stream 
      @param out output stream 
      @param indent_width indentation width */
      void display( ostream& out, size_t const& indent_width ) const;    
      
      /** @brief Computes the hydrodynamic force and torque and store the values
      in the corresponding geometric sphere
      @param PP the pressure field
      @param UU the velocity field */
      void compute_hydro_force_torque( FV_DiscreteField const* PP,
	FV_DiscreteField const* UU );       
      //@} 


   protected: //--------------------------------------------------------------

   //-- Attributes  

      /**@name Parameters */
      //@{ 
      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_3Dcylinder object */
      DS_3Dcylinder( DS_3Dcylinder const& copy );
      //@}    
};

#endif
