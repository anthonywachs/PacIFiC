#ifndef _DS_ALLRIGIDBODIES__
#define _DS_ALLRIGIDBODIES__

#include <geomVector.hh>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
using std::vector;
using std::istream ;
using std::ostream ;
using std::istringstream ;
using std::string;
class DS_RigidBody;
class FS_AllRigidBodies;
class FV_DiscreteField;


/** @brief The class DS_AllRigidBodies.

The array of all rigid bodies in the Direction Splitting solver. 
 
@author A. Wachs - Pacific project 2021 */

class DS_AllRigidBodies
{
   public: //-----------------------------------------------------------------
  
   //-- Constructors & Destructor
  
      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_AllRigidBodies();
      
      /** @brief Constructor with arguments
      @param dimens number of space dimensions
      @param in input stream where features of rigid bodies are read
      @param b_particles_as_fixed_obstacles treat all rigid bodies as fixed 
      obstacles */
      DS_AllRigidBodies( size_t& dimens, istream& in, 
      	bool const& b_particles_as_fixed_obstacles );       	   

      /** @brief Destructor */
      ~DS_AllRigidBodies();
      //@}

      
   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns the total number of rigid bodies */
      size_t get_number_rigid_bodies() const;
      
      /** @brief Returns the number of particles */
      size_t get_number_particles() const;
      
      /** @brief Returns a constant pointer to the FS_AllRigidBodies object that
      contains the vector of all corresponding geometric rigid bodies */
      FS_AllRigidBodies const* get_ptr_FS_AllRigidBodies() const;
      
      /** @brief Returns a const pointer to the ith DS rigid body */
      DS_RigidBody const* get_ptr_rigid_body( size_t i ) const; 
      
      /** @brief Returns a pointer to the ith DS rigid body */
      DS_RigidBody* get_ptr_rigid_body( size_t i );                        
      //@}
      
      
   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates all rigid bodies
      @param in input stream where features of rigid bodies are read */
      void update( istream& in );
      //@}      


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the geometric attributes of the rigid bodies in a stream
      @param out output stream 
      @param indent_width indentation width */
      void display_geometric( ostream& out, size_t const& indent_width ) const;
      
      /** @brief Writes the attributes of the rigid bodies in a stream 
      @param out output stream 
      @param indent_width indentation width */
      void display( ostream& out, size_t const& indent_width ) const; 
      
      /** @brief Computes the hydrodynamic force and torque acting on each rigid
      body and store the values in its corresponding geometric rigid body
      @param PP the pressure field
      @param UU the velocity field */
      void compute_hydro_force_torque( FV_DiscreteField const* PP,
	FV_DiscreteField const* UU );           
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
      vector<DS_RigidBody*> m_allDSrigidbodies; /**< the vector of all 
    	Direction Splitting rigid bodies */
      FS_AllRigidBodies* m_FSallrigidbodies; /**< the pointer to the 
    	FS_AllRigidBodies object that contains the vector of all 
    	corresponding geometric rigid bodies */	        
      //@}  


    //-- Constructors & Destructor  
      
       /**@name Constructors & Destructor */
       //@{
       /** @brief Copy constructor
       @param copy copied DS_AllRigidBodies object */
       DS_AllRigidBodies( DS_AllRigidBodies const& copy );
       //@}    
};

#endif
