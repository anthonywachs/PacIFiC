#ifndef _CONTACTBUILDERFACTORY_HH_
#define _CONTACTBUILDERFACTORY_HH_

#include "ReaderXML.hh"
#include <iostream>
#include <map>
#include <string>
using namespace std;

class ContactForceModel;
class Particle;
class Obstacle;


/** @brief The class ContactBuilderFactory.

    Creates and stores tables (couple of materials)-(contact force model) and
    returns the contact force model for a given couple of materials.

    @author Grains Project - IFP - 2007 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ContactBuilderFactory
{
  public:
    /** @name Static methods */
    //@{
    /** @brief Returns the contact associated to a couple of materials
    @param matA material A
    @param matB material B */
    static ContactForceModel* contactForceModel( string const& matA,
	string const& matB );

    /** @brief Creates the table (couple of materials)-(contact force model)
    from an XML node
    @param root XML node */
    static void define( DOMNode* root );

    /** @brief Adds a material to the simulation
    @param material material name
    @param matForObstacle whether the material is for an obstacle */
    static void defineMaterial( string const& material,
	bool const& matForObstacle );

    /** @brief Reloads and recreates contact force models and tables
    @param file input stream */
    static void reload( istream& file );

    /** @brief Makes sure that materials for obstacles only are identified after
    a reload as these data are not stored in the reload file
    @param referenceParticles reference particles */
    static void set_materialsForObstaclesOnly_reload(
   	vector<Particle*> const* referenceParticles );

    /** @brief Writes the whole contact force model data in a stream
    @param file output stream
    @param contactFile root name of the contact reload file
    @param rank rang du processus */
    static void save( ostream& file, string const& contactFile,
    	int const& rank );

    /** @brief Returns the name of the material given its ID number
    @param num ID number in the map(string,int) */
    static string materialName( int num );

    /** @brief Returns the ID number of the material given its name
    @param mat material name in the map(string,int) */
    static int materialIDnumber( string const& mat );

    /** @brief Checks that all couples of materials in which one of the
    materials is assigned to a particle has an associated contact force model.
    If not, the two material names are copied in matA and matB
    @param matA material A
    @param matB material B */
    static bool checkContactForceModelsExist( string& matA, string& matB );

    /** @brief Deletes all contact force models */
    static void eraseAllContactForceModels();
    //@}


  private:
    /** @name Enumeration */
    //@{
    enum ContactType {
       /// Hookean spring, dashpot and Coulomb friction
       Hooke,
       /// Memory model: hookean spring and dashpot in all directions; Coulomb
       /// friction in tangential and rotational directions
       HookeMemory,
       /// Hertz spring, dashpot and Coulomb friction
       Hertz,
       /// Memory model: Hertz non-linear spring and dashpot in all directions; 
       /// Coulomb friction in tangential and rotational directions
       HertzMemory              
    };
    //@}


    /** @name Structure */
    //@{
    struct ContactFeatures
    {
      ContactType name;
      map<string,double> values;
    };
    //@}


    /** @name Methods */
    //@{
    /** @brief Reads and returns the contact parameters of the contact force
    model from an XML node
    @param root XML node */
    static pair<ContactBuilderFactory::ContactFeatures,ContactForceModel*>
  	defineParameters( DOMNode *root ) ;
    //@}


    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    ContactBuilderFactory() {};

    /** @brief Destructor (forbidden) */
    ~ContactBuilderFactory() {};
    //@}


    /** @name Parameters */
    //@{
    static map<string,int> m_materials; /**< ID number of materials in the table
    	(material names)-(material ID numbers) */
    static map<string,bool> m_materialsForObstaclesOnly; /**< table that defines
    	whether a material is associated to obstacles only */
    static map<int,ContactBuilderFactory::ContactFeatures>
  	m_contactParametres; /**< (contact ID)-(contact features) table */
    static map<int,ContactForceModel*> m_contactForceModels; /**<
    	(contact ID)-(contact force model) table */
    static int m_value; /**< variable to define material ID number */
    //@}
};

#endif
