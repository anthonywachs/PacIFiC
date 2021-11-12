#ifndef _DS_RIGIDBODY_BUILDERFACTORY__
#define _DS_RIGIDBODY_BUILDERFACTORY__

#include <iostream>
using std::istream;
class FS_RigidBody;
class DS_RigidBody;


/** @brief The Class DS_RigidBody_BuilderFactory.

Rigid body builder factory, creates an instance of rigid body using a unsigned
integer (numebr of corners for a polyhedron or a code for other shapes) as 
the main parameter.

@author A. Wachs - Pacific project 2021 */

class DS_RigidBody_BuilderFactory
{
   public: //-----------------------------------------------------------------
   
   //-- Static methods

      /** @name Static methods */
      //@{
      /** @brief Creates a Direction Splitting rigid body
      @param pgrb pointer to the corresponding geometric rigid body */
      static DS_RigidBody* create( FS_RigidBody* pgrb );

	
   protected: //--------------------------------------------------------------
	

   private: //----------------------------------------------------------------
   
   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_RigidBody_BuilderFactory();

      /** @brief Destructor */
      ~DS_RigidBody_BuilderFactory();
      //@}   
};

#endif
