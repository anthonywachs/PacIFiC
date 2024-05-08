#ifndef _GRAINSEXEC_HH_
#define _GRAINSEXEC_HH_

#include <mpi.h>
#include "Vector3.hh"
#include "Point3.hh"
#include "Matrix.hh"
#include "VertexBase.hh"
#include "IndexArray.hh"
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <set>
using std::string;
using std::ios_base;
using std::list;
using std::vector;
using std::set;
using namespace solid;

class GrainsMPIWrapper;
class App;


/** @brief Insertion window type */
enum WindowType
{
  WINDOW_BOX, /**< Box */
  WINDOW_CYLINDER, /**< Cylinder */
  WINDOW_ANNULUS, /**< Annulus */
  WINDOW_LINE, /**< Line */
  WINDOW_NONE /**< unknown */
};


/** @brief Insertion window feature */
struct Window
{
  WindowType ftype; /**< Window type */
  Point3 ptA; /**< Box 1st corner or center of lower disk of the cylinder */
  Point3 ptB; /**< Box 2nd corner */
  double radius; /**< Cylinder radius */
  double radius_int; /**< Inner cylinder radius in case of annulus */
  double height; /**< Cylinder height */
  Direction axisdir; /**< Cylinder axis direction */
};


/** @brief Structured array of positions */
struct StructArrayInsertion
{
  struct Window box; /**< Window */
  size_t NX; /**< Number of positions in x */
  size_t NY; /**< Number of positions in y */
  size_t NZ; /**< Number of positions in z */
};


/** @brief The class GrainsExec.

    Fully static class that manages global variables.

    @author A. WACHS - Institut Francais du Petrole - 2010 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class GrainsExec
{
  public:
    /** @name Get methods */
    //@{
    /** @brief Returns a pointer to the MPI wrapper */
    static GrainsMPIWrapper* getComm();

    /** @brief Returns the list of applications */
    static list<App*> get_listApp();

    /** @brief Returns the total number of particles in the physical system 
    (i.e. on all subdomains/processes), i.e. sum of total number of active 
    particles with tag 0 or 1 and inactive particles */
    static size_t getTotalNumberPhysicalParticles();
    //@}


    /** @name Set methods */
    //@{
    /** @brief Sets the MPI wrapper
    @param wrapper_ wrapper */
    static void setComm( GrainsMPIWrapper* wrapper_ );

    /** @brief Sets the list of applications
    @param allApp_ list of pointers to applications */
    static void set_listApp( list<App*> allApp_ );

    /** @brief Sets the total number of particles in the physical system 
    (i.e. on all subdomains/processes), i.e. sum of total number of active 
    particles with tag 0 or 1 and inactive particles
    @param nb_ total number of physical particles */
    static void setTotalNumberPhysicalParticles( size_t const& nb_ );
    //@}


    /** @name Useful methods */
    //@{
    /** @brief Writes a float number with a prescribed number of digits in a
    string
    @param figure the float number
    @param size number of digits */
    static string doubleToString( double const& figure, int const& size );

    /** @brief Writes a float number with a prescribed format and a prescribed
    number of digits after the decimal point in a string
    @param format the format
    @param digits number of digits after the decimal point
    @param number the float number */
    static string doubleToString( ios_base::fmtflags format, int digits,
      	double const& number );

    /** @brief Writes an integer in a string
    @param figure the integer */
    static string intToString( int const& figure );

    /** @brief Checks the last output time in a file and deletes all subsequent
    lines, i.e., each line whose time is after the current time
    @param filename file name
    @param current_time current physical time */
    static void checkTime_outputFile( string const& filename,
	const double& current_time );

    /** @brief Returns the root of a complete file name. Example: if input is
    Titi/tutu/toto, returns Titi/tutu
    @param FileName complete file name */
    static string extractRoot( string const& FileName );

    /** @brief Returns the file name of a complete file name. Example: if input
    is Titi/tutu/toto, returns toto
    @param FileName complete file name */
    static string extractFileName( string const& FileName );

    /** @brief Checks that all reload files are in the same directory (primarily
    checks that files for polyhedrons et polygons are there) */
    static void checkAllFilesForReload();

    /** @brief Returns the reload file name from the restart time table
    @param rootName root name of the restart time table file
    @param RFTable_ext extension of the restart time table file (in general
    	"_RFTable.txt") */
    static string restartFileName_AorB( string const& rootName,
  	string const& RFTable_ext );

    /** @brief Returns a random rotation matrix
    @param dim number of space dimensions */
    static Matrix RandomRotationMatrix( size_t dim );

    /** @brief Returns a random unit vector
    @param dim number of space dimensions */
    static Vector3 RandomUnitVector( size_t dim );

    /** @brief Returns whether a sphere is fully in, fully out or intersects
    an axis-aligned cylinder. Returned values are 0, 1 and 2 respectively
    @param SphereCenter sphere center
    @param SphereRadius sphere radius
    @param CylBottomCentre center of lower disk of the cylinder
    @param CylRadius cylinder radius
    @param CylHeight cylinder height
    @param CylAxisDir cylinder axis direction
    @param tol tolerance for the test fully in */
    static size_t AACylinderSphereIntersection( Point3 const& SphereCenter,
    	double const& SphereRadius,
	Point3 const& CylBottomCentre,
	double const& CylRadius,
	double const& CylHeight,
	size_t const& CylAxisDir,
	double const& tol = EPSILON );

    /** @brief Returns whether a point belongs to an axis-aligned cylinder.
    @param pt the point
    @param CylBottomCentre center of lower disk of the cylinder
    @param CylRadius cylinder radius
    @param CylHeight cylinder height
    @param CylAxisDir cylinder axis direction
    @param tol tolerance */
    static bool isPointInAACylinder( Point3 const& pt,
	Point3 const& CylBottomCentre,
	double const& CylRadius,
	double const& CylHeight,
	size_t const& CylAxisDir,
	double const& tol = EPSILON );
	
    /** @brief Computes and returns the 4 x 4 determinant of 4 points
    @param p1 point 1    
    @param p2 point 2    
    @param p3 point 3 
    @param p4 point 4 */
    static double PointDeterm4by4( Point3 const& p1, Point3 const& p2,
    	Point3 const& p3, Point3 const& p4 );
	
    /** @brief Returns whether a point lies inside a tetrahedron
    @param p1 point 1    
    @param p2 point 2    
    @param p3 point 3 
    @param p4 point 4 
    @param p point to check 
    @param check additional geometric checks if true, default value is false */
    static bool isPointInTetrahedron( Point3 const& p1, Point3 const& p2,
    	Point3 const& p3, Point3 const& p4, Point3 const& p, 
	bool check = false );
	
    /** @brief Returns the full result file name
    @param rootname root file name 
    @param addrank add rank number */
    static string fullResultFileName( string const& rootname,
    	bool addrank = true ); 	
    //@}


    /** @name Garbage collector for shape description */
    //@{
    /** @brief Frees the content of the garbage collector */
    static void GarbageCollector();

    /** @brief Adds a list of reference points of a polytope
    @param refPB list of reference points used by the polytopes
    @param refVB pointer to the list of reference points and usd to access these
    	reference points */
    static void addOnePolytopeRefPointBase( Point3* refPB, VertexBase* refVB );

    /** @brief Adds a description of vertices of a polytope
    @param idar description of vertices */
    static void addOnePolytopeNodeNeighbors( IndexArray* idar );

    /** @brief Adds an array of indices of vertex of a polytope
    @param idar array of indices of vertex */
    static void addOnePolytopeNodeIndex( IndexArray* idar );

    /** @brief Adds a face connectivity of a polyhedron
    @param faceCon face connectivity */
    static void addOnePolyhedronFaceConnectivity(
  	vector< vector<int> >* faceCon );
    //@}


    /** @name Used memory tracing */
    //@{
    /** @brief Returns memory used by this process */
    static size_t used_memory( void );

    /** @brief Writes memory used by this process in a stream
    @param os output stream
    @param memory used memory */
    static void display_memory( ostream& os, size_t memory );
    //@}


    /** @name Parameters */
    //@{
    static bool m_MPI; /**< whether the computation is serial or MPI */
    static string m_TIScheme; /**< Time integration scheme type */
    static bool m_SphereAsPolyParaview; /**< in Paraview, true if spheres are
    	faceted as polyhedrons, false if post-processed as a vectorial field */
    static int m_MPI_verbose; /**< MPI verbosity level, 2 levels: 0=none,
    	1=particles */
    static bool m_isReloaded; /**< whether the simulation starts from a reload 
    	state */ 
    static string m_ReloadType; /**< Reload type: "new" for a new simulation and
    	 "same" for an on-going simulation */
    static Vector3 m_vgravity; /**< gravity vector */
    static Vector3* m_translationParaviewPostProcessing; /**< in case of a
    	coupling with the fluid in the translation-projection mode; if not
	active, pointer is NULL *? */
    static bool m_periodic; /**< Is the domain periodic ? */
    static bool m_isGrainsCompFeatures; /**< Is the Grains simulation of type
    	CompFeatures ? */
    static bool m_isGrainsPorosity; /**< Is the Grains simulation of type
    	Porosity ? */
    static string m_ReloadDirectory; /**< Directory where reload files were
    	read */
    static string m_SaveDirectory; /**< Directory where reload files are
    	written */
    static bool m_SaveMPIInASingleFile; /**< In MPI, true if all processes write
    	the particle data in a single file, false if each process writes its own
	data */
    static bool m_ReadMPIInASingleFile; /**< In MPI, true if all processes read
    	the particle data in a single file, false if each process reads its own
	data */		
    static set<string> m_additionalDataFiles; /**< additional files for reload
    	(primarily files for polyhedrons et polygons) */
    static bool m_writingModeHybrid; /**< Is writing mode hybrid, i.e., a text
    	header for reference particles and obstacles, and a binary file for
	particles */
    static bool m_readingModeHybrid; /**< Is reading mode hybrid, i.e., a text
    	header for reference particles and obstacles, and a binary file for
	particles */	
    static string m_GRAINS_HOME; /**< Main Grains directory */
    static string m_reloadFile_suffix; /**< Reload file suffix (A or B) */
    static bool m_exception_Contact; /**< Contact exception */
    static bool m_exception_Displacement; /**< Displacement exception */
    static bool m_exception_Simulation; /**< Simulation exception */
    static string m_shift0; /**< string of 0 blank space */
    static string m_shift1; /**< string of 1 blank space */
    static string m_shift2; /**< string of 2 blank spaces */
    static string m_shift3; /**< string of 3 blank spaces */
    static string m_shift6; /**< string of 6 blank spaces */
    static string m_shift9; /**< string of 9 blank spaces */
    static string m_shift12; /**< string of 12 blank spaces */
    static string m_shift15; /**< string of 15 blank spaces */
    static bool m_output_data_at_this_time; /**< writes data in
    	result files at this time */
    static string m_inputFile; /**< Grains3D major input file */
    static int m_return_syscmd; /**< Returned value of system command */
    static bool m_preCollision_cyl; /** precollision test w/ bounding cyls **/
    static Point3 m_defaultInactivePos; /**< Default position of inactive 
    	particles */    
    //@}


  private:
    /** @name Parameters */
    //@{
    static GrainsMPIWrapper* m_wrapper; /**< Wrapper MPI */
    static list<App*> m_allApp; /**< List of all applications used in
    	the simulation (the 1st application is contact detection, i.e., the
	LinkedCell) */
    static size_t m_total_nb_physical_particles; /**< total number of particles
    	in the physical system (i.e. on all subdomains/processes), i.e. sum of 
	total number of active particles with tag 0 or 1 and inactive 
	particles */
    static list< pair<Point3*,VertexBase *> > m_allPolytopeRefPointBase; /**<
  	list of reference points used by the polytopes */
    static list<IndexArray*> m_allPolytopeNodeNeighbors; /**< list of
    	descriptions of vertices of polytopes (pointers to arrays of
	IndexArray) */
    static list<IndexArray*> m_allPolytopeNodesIndex; /**< list of arrays of
  	indices of vertex used by the polytopes */
    static list<vector< vector<int> >*> m_allPolyhedronFacesConnectivity; /**<
  	list of face connectivity in the polyhedrons */
    //@}


    /** @name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    GrainsExec();

    /** @brief Destructor */
    ~GrainsExec();
    //@}
};

#endif
