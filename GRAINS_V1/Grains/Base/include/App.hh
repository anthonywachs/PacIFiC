#ifndef _APP_HH_
#define _APP_HH_

#include "Particle.hh"
#include "SimpleObstacle.hh"
#include "Basic.hh"
#include "Error.hh"
#include <list>
using namespace std;


/** @brief The class App.

    Used to compute forces and torques exerted on rigid bodies.

    @author Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class App
{
  public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    App();

    /** @brief Destructor */
    virtual ~App();
    //@}


    /**@name Virtual methods */
    //@{
    /** @brief Computes forces and torques exerted on rigid bodies
    @param time physical time
    @param dt time step magnitude 
    @param particles active particles */
    virtual void ComputeForces( double time, double dt,
    	list<Particle*> const* particles ) = 0;  
    //@}

 
    /** @name Methods */
    //@{  
    /** @brief Sets the name of the application
    @param name_ app name */
    void setName( string const& name_ );
  
    /** @brief Returns the name of the application */
    string getName() const;
    
    /** @brief Returns whether the name of the application matches an input name
    @param name_ name */
    bool isName( string const& name_ );         
    //@}


    /** @name Static methods */
    //@{
    /** @brief Sets domain dimensions
    @param xmax maximum coordinate in x 
    @param ymax maximum coordinate in y
    @param zmax maximum coordinate in z
    @param ox x-coordinate of origin
    @param oy y-coordinate of origin 
    @param oz z-coordinate of origin */
    static void set_dimensions( double xmax, double ymax, double zmax,
  	double ox, double oy, double oz );

    /** @brief Sets domain periodicity
    @param vper vector of periodicity (3 booleans) */
    static void set_periodicity( vector<bool> const& vper );
  
    /** @brief Sets local domain size
    @param lx_ length in x 
    @param ly_ length in y
    @param lz_ length in z */
    static void set_local_domain_size( double lx_, double ly_, double lz_ );
  
    /** @brief Sets local domain origin
    @param nprocsdir number of subdomains in each direction
    @param MPIcoords subdomain coordinates in the MPI Cartesian topology */
    static void set_local_domain_origin( int const* nprocsdir, 
    	int const* MPIcoords );   
  
    /** @brief Gets local domain origin 
    @param x x-coordinate of local domain origin
    @param y y-coordinate of local domain origin
    @param z z-coordinate of local domain origin */
    static void get_local_domain_origin( double& x, double& y, double& z );   

    /** @brief Gets global domain origin
    @param x x-coordinate of global domain origin
    @param y y-coordinate of global domain origin
    @param z z-coordinate of global domain origin */
    static void get_origin( double& x, double& y, double& z );   

    /** @brief Gets local domain size
    @param lx length in x 
    @param ly length in y
    @param lz length in z */
    static void get_local_domain_size( double& lx, double& ly, double& lz );   

    /** @brief Gets global domain size
    @param lx length in x 
    @param ly length in y
    @param lz length in z */
    static void get_size( double& lx, double& ly, double& lz );   
  
    /** @brief Returns whether a point belongs to the global domain
    @param position point */
    static bool isInDomain( Point3 const* position );
    
    /** @brief Returns whether a point belongs to the global domain in a given 
    direction
    @param position point 
    @param dir direction */
    static bool isInDomain( Point3 const* position, size_t const& dir );    
  
    /** @brief Returns whether a point belongs to the local domain
    @param position point */
    static bool isInLocalDomain( Point3 const* position ); 
    
    /** @brief Writes the domain features in an output stream
    @param output output stream 
    @param oshift empty string to shift the output */
    static void output_domain_features( ostream& output, string const& oshift );

    /** @brief Returns whether the domain is periodic in a direction
    @param dir direction */
    static bool isPeriodic( size_t dir ); 
    //@}  

  
  protected:
    /**@name Parameters */
    //@{  
    string m_name; /**< Application name */
    //@}

    /**@name Parameters Static */
    //@{  
    static Vector3 m_domain_global_size; /**< Domain global size */     
    static Vector3 m_domain_local_size; /**< Domain local size */   
    static Point3 m_domain_global_origin; /**< Domain global origin */  
    static Point3 m_domain_local_origin; /**< Domain local origin */
    static Point3 m_domain_global_max; /**< Domain global max point */    
    static Point3 m_domain_local_max; /**< Domain local max point */       
    static vector<bool> m_domain_global_periodicity; /**< vector of domain 
    	periodicity (3 booleans) */ 
    static bool m_domain_global_periodic; /**< true if the domain is periodic 
    	in at least one direction */
    static vector<Vector3> m_domain_global_periodic_vectors; /**< 
    	periodic vectors as a function of the geographic position */
    static vector< vector<int> > m_periodic_vector_indices; /**< periodic vector
    	indices (geographic positions) as a function of the geographic position.
	For instance, GEOPOS_NORTH_WEST has GEOPOS_NORTH, GEOPOS_WEST and 
	GEOPOS_NORTH_WEST */
    //@}
    

  private:
    /** @name Static methods */
    //@{
    /** @brief Writes all periodic vectors for a given geographic position in 
    an output stream
    @param output output stream 
    @param oshift empty string to shift the output
    @param geoloc_ geographic position */
    static void output_periodic_vectors_per_geopos( ostream& output, 
    	string const& oshift, GeoPosition geoloc_ );
    //@}      
};

#endif

     
