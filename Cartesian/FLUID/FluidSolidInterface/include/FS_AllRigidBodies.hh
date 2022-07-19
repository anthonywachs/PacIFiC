#ifndef _FS_ALLRIGIDBODIES__
#define _FS_ALLRIGIDBODIES__

#include <geomVector.hh>
#include <vector>
#include <iostream>
#include <sstream>
using std::vector;
using std::istream ;
using std::ostream ;
using std::istringstream ;
class FS_RigidBody;


/** @brief The class FS_AllRigidBodies.

The array of all rigid bodies. 
 
@author A. Wachs - Pacific project 2021 */

class FS_AllRigidBodies
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_AllRigidBodies(); 
      
      /** @brief Constructor with arguments
      @param dimens number of space dimensions 
      @param in input stream where features of rigid bodies are read 
      @param b_particles_as_fixed_obstacles treat all rigid bodies as fixed 
      obstacles */
      FS_AllRigidBodies( size_t& dimens, istream& in, 
      	bool const& b_particles_as_fixed_obstacles );       	   

      /** @brief Destructor */
      ~FS_AllRigidBodies();
      //@}


   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns the total number of rigid bodies */
      size_t get_number_rigid_bodies() const;
      
      /** @brief Returns the number of particles */
      size_t get_number_particles() const; 
      
      /** @brief Returns a pointer to the ith geometric rigid body
      @param i rigid body number */
      FS_RigidBody* get_ptr_rigid_body( size_t i ) const;             
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates all rigid bodies
      @param in input stream where features of rigid bodies are read */
      void update( istream& in ); 
      
      /** @brief Sets the translational and angular velocity of all rigid bodies
      to zero */
      void nullify_velocity();           
      //@}  

      
   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the attributes in a stream 
      @param out output stream 
      @param indent_width indentation width */
      void display( ostream& out, size_t const& indent_width ) const;
      //@}
      
      
   protected: //--------------------------------------------------------------


   private: //----------------------------------------------------------------

   //-- Attributes  

      /**@name Parameters */
      //@{
      size_t m_space_dimension; /**< Space dimension */
      size_t m_npart; /**< number of particles */
      size_t m_nrb; /**< total number of rigid bodies = number of 
      	particles + number of obstacles, npart first rigid bodies are always 
	particles while ( m_nrb - m_npart ) last rigid bodies are obstacles */
      vector<FS_RigidBody*> m_allrigidbodies; /**< the vector of all rigid 
    	bodies */        
      //@}  


   //-- Constructors & Destructor
  
      /**@name Constructors  */
      //@{
      /** @brief Copy constructor
      @param copy copied FS_AllRigidBodies object */
      FS_AllRigidBodies( FS_AllRigidBodies const& copy );
      //@}    
};

#endif
