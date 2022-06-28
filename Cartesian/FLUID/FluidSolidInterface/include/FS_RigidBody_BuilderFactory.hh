#ifndef _FS_RIGIDBODY_BUILDERFACTORY__
#define _FS_RIGIDBODY_BUILDERFACTORY__

#include <iostream>
using std::istream;
class FS_RigidBody;


/** @brief The Class FS_RigidBody_BuilderFactory.

Rigid body builder factory, creates an instance of rigid body using a unsigned
integer (numebr of corners for a polyhedron or a code for other shapes) as 
the main parameter.

@author A. Wachs - Pacific project 2021 */

class FS_RigidBody_BuilderFactory
{
   public: //-----------------------------------------------------------------
   
   //-- Static methods

      /** @name Static methods */
      //@{
      /** @brief Creates a rigid body
      @param dimens number of space dimensions  
      @param in input stream    
      @param ncorners_ number of corners or a code 
      @param id_ identification number */
      static FS_RigidBody* create( size_t& dimens, istream& in, 
      	size_t& ncorners_, size_t& id_ );

	
   protected: //--------------------------------------------------------------
	

   private: //----------------------------------------------------------------
   
   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_RigidBody_BuilderFactory();

      /** @brief Destructor */
      ~FS_RigidBody_BuilderFactory();
      //@}   
};

#endif
