#include "App.hh"
#include "AllComponents.hh"
#include "GrainsExec.hh"

double App::m_domain_global_size_X = 0.;
double App::m_domain_global_size_Y = 0.;
double App::m_domain_global_size_Z = 0.;
double App::m_domain_local_size_X = 0.;
double App::m_domain_local_size_Y = 0.;
double App::m_domain_local_size_Z = 0.;
Point3 App::m_domain_local_origin;
Point3 App::m_domain_global_origin;
vector<bool> App::m_domain_global_periodicity;
bool App::m_domain_global_periodic = false;
vector<Vector3> App::m_domain_global_periodic_vectors;
vector< vector<int> > App::m_periodic_vector_indices;


// ----------------------------------------------------------------------------
// Default constructor
App::App()
{}




// ----------------------------------------------------------------------------
// Destructeur
App::~App()
{}




// ----------------------------------------------------------------------------
// Sets domain dimensions
void App::set_dimensions( double xmax, double ymax, double zmax,
  	double ox, double oy, double oz )
{
  App::m_domain_global_origin[0] = ox;
  App::m_domain_global_origin[1] = oy;
  App::m_domain_global_origin[2] = oz;
  App::m_domain_local_origin[0] = ox;
  App::m_domain_local_origin[1] = oy;
  App::m_domain_local_origin[2] = oz; 
  App::m_domain_global_size_X = xmax - m_domain_global_origin[0];
  App::m_domain_global_size_Y = ymax - m_domain_global_origin[1];
  App::m_domain_global_size_Z = zmax - m_domain_global_origin[2];
  App::m_domain_local_size_X = xmax - m_domain_global_origin[0];
  App::m_domain_local_size_Y = ymax - m_domain_global_origin[1];
  App::m_domain_local_size_Z = zmax - m_domain_global_origin[2];
  assert( App::m_domain_global_size_X >= 0. && App::m_domain_global_size_Y >= 0.
 	&& App::m_domain_global_size_Z >= 0. );   
}




// ----------------------------------------------------------------------------
// Sets domain periodicity
void App::set_periodicity( vector<bool> const& vper )
{
  m_domain_global_periodicity = vper;
  if ( vper[X] || vper[Y] || vper[Z] ) m_domain_global_periodic = true;
  else m_domain_global_periodic = false;
  
  if ( m_domain_global_periodic_vectors.empty() )
  {
    m_domain_global_periodic_vectors.reserve( 27 );
    for (int i=0;i<27;++i) 
      m_domain_global_periodic_vectors.push_back( Vector3Nul );
      
    // Single periodicity
    // West East
    if ( m_domain_global_periodicity[X] )
    {
      m_domain_global_periodic_vectors[GEOPOS_WEST][X] = 
      	m_domain_global_size_X ;  
      m_domain_global_periodic_vectors[GEOPOS_EAST][X] = 
      	- m_domain_global_size_X ;
    }

    // South North
    if ( m_domain_global_periodicity[Y] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH][Y] = 
      	m_domain_global_size_Y ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH][Y] = 
      	- m_domain_global_size_Y ;
    }

    // Behind Front
    if ( m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_BEHIND][Z] = 
      	m_domain_global_size_Z ;
      m_domain_global_periodic_vectors[GEOPOS_FRONT][Z] = 
      	- m_domain_global_size_Z ;
    }


    // Double periodicity    
    // South North West East
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Y] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST][X] = 
      	m_domain_global_size_X ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST][Y] = 
      	m_domain_global_size_Y ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST][X] = 
      	- m_domain_global_size_X ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST][Y] = 
      	m_domain_global_size_Y ;        
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST][X] = 
      	m_domain_global_size_X ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST][Y] = 
      	- m_domain_global_size_Y ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST][X] = 
      	- m_domain_global_size_X ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST][Y] = 
      	- m_domain_global_size_Y ;
    }
    
    // South North Behind Front
    if ( m_domain_global_periodicity[Y] && m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_BEHIND][Y] = 
      	m_domain_global_size_Y ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_BEHIND][Z] = 
      	m_domain_global_size_Z ;	     
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_FRONT][Y] = 
      	m_domain_global_size_Y ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_FRONT][Z] = 
      	- m_domain_global_size_Z ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_BEHIND][Y] = 
      	- m_domain_global_size_Y ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_BEHIND][Z] = 
      	m_domain_global_size_Z ;     
      m_domain_global_periodic_vectors[GEOPOS_NORTH_FRONT][Y] = 
      	- m_domain_global_size_Y ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_FRONT][Z] = 
      	- m_domain_global_size_Z ;
    }      

    // West East Behind Front
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_WEST_BEHIND][X] = 
      	m_domain_global_size_X ;    
      m_domain_global_periodic_vectors[GEOPOS_WEST_BEHIND][Z] = 
      	m_domain_global_size_Z ;    
      m_domain_global_periodic_vectors[GEOPOS_WEST_FRONT][X] = 
      	m_domain_global_size_X ;    
      m_domain_global_periodic_vectors[GEOPOS_WEST_FRONT][Z] = 
      	- m_domain_global_size_Z ;   
      m_domain_global_periodic_vectors[GEOPOS_EAST_BEHIND][X] = 
      	- m_domain_global_size_X ;    
      m_domain_global_periodic_vectors[GEOPOS_EAST_BEHIND][Z] = 
      	m_domain_global_size_Z ;
      m_domain_global_periodic_vectors[GEOPOS_EAST_FRONT][X] = 
      	- m_domain_global_size_X ;    
      m_domain_global_periodic_vectors[GEOPOS_EAST_FRONT][Z] = 
      	- m_domain_global_size_Z ; 
    }   

    
    // Triple periodicity
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Y] 
    	&& m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_BEHIND][X] = 
      	m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_BEHIND][Y] = 
      	m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_BEHIND][Z] = 
      	m_domain_global_size_Z ;  		     
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_FRONT][X] = 
      	m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_FRONT][Y] = 
      	m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_FRONT][Z] = 
      	- m_domain_global_size_Z ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_BEHIND][X] = 
      	m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_BEHIND][Y] = 
      	- m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_BEHIND][Z] = 
      	m_domain_global_size_Z ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_FRONT][X] = 
      	m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_FRONT][Y] = 
      	- m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_FRONT][Z] = 
      	- m_domain_global_size_Z ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_BEHIND][X] = 
      	- m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_BEHIND][Y] = 
      	m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_BEHIND][Z] = 
      	m_domain_global_size_Z ;  
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_FRONT][X] = 
      	- m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_FRONT][Y] = 
      	m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_FRONT][Z] = 
      	- m_domain_global_size_Z ;  
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_BEHIND][X] = 
      	- m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_BEHIND][Y] = 
      	- m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_BEHIND][Z] = 
      	m_domain_global_size_Z ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_FRONT][X] = 
      	- m_domain_global_size_X ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_FRONT][Y] = 
      	- m_domain_global_size_Y ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_FRONT][Z] = 
      	- m_domain_global_size_Z ;
    }            
  }
  
  if ( m_periodic_vector_indices.empty() )
  {
    m_periodic_vector_indices.reserve(30);
    vector<int> emptyVECINT;
    for (int i=0;i<30;++i) 
      m_periodic_vector_indices.push_back(emptyVECINT);
    vector<int>* work = NULL;
  
    // NORTH => NORTH
    work = new vector<int>(1,0);
    (*work)[0] = GEOPOS_NORTH;
    m_periodic_vector_indices[GEOPOS_NORTH] = *work;
    work->clear();
    delete work;
  
    // NORTH_EAST => NORTH, EAST, NORTH_EAST
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_EAST;  
    (*work)[2] = GEOPOS_NORTH_EAST;
    m_periodic_vector_indices[GEOPOS_NORTH_EAST] = *work;    
    work->clear();
    delete work;
  
    // NORTH_WEST => NORTH, WEST, NORTH_WEST
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_WEST;  
    (*work)[2] = GEOPOS_NORTH_WEST;
    m_periodic_vector_indices[GEOPOS_NORTH_WEST] = *work;    
    work->clear();
    delete work;    

    // NORTH_FRONT => NORTH, FRONT, NORTH_FRONT
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_FRONT;  
    (*work)[2] = GEOPOS_NORTH_FRONT;
    m_periodic_vector_indices[GEOPOS_NORTH_FRONT] = *work;    
    work->clear();
    delete work;  

    // NORTH_BEHIND => NORTH, BEHIND, NORTH_BEHIND
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_BEHIND;  
    (*work)[2] = GEOPOS_NORTH_BEHIND;
    m_periodic_vector_indices[GEOPOS_NORTH_BEHIND] = *work;    
    work->clear();
    delete work; 

    // NORTH_EAST_FRONT => NORTH, EAST, FRONT, EAST_FRONT, NORTH_EAST, 
    // NORTH_FRONT, NORTH_EAST_FRONT
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_EAST;  
    (*work)[2] = GEOPOS_FRONT;
    (*work)[3] = GEOPOS_EAST_FRONT;  
    (*work)[4] = GEOPOS_NORTH_EAST;  
    (*work)[5] = GEOPOS_NORTH_FRONT;  
    (*work)[6] = GEOPOS_NORTH_EAST_FRONT;  
    m_periodic_vector_indices[GEOPOS_NORTH_EAST_FRONT] = *work;    
    work->clear();
    delete work; 

    // NORTH_EAST_BEHIND => NORTH, EAST, BEHIND, EAST_BEHIND, NORTH_EAST, 
    // NORTH_BEHIND, NORTH_EAST_BEHIND
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_EAST;  
    (*work)[2] = GEOPOS_BEHIND;
    (*work)[3] = GEOPOS_EAST_BEHIND;  
    (*work)[4] = GEOPOS_NORTH_EAST;  
    (*work)[5] = GEOPOS_NORTH_BEHIND;  
    (*work)[6] = GEOPOS_NORTH_EAST_BEHIND;  
    m_periodic_vector_indices[GEOPOS_NORTH_EAST_BEHIND] = *work;    
    work->clear();
    delete work; 

    // NORTH_WEST_FRONT => NORTH, WEST, FRONT, WEST_FRONT, NORTH_WEST, 
    // NORTH_FRONT, NORTH_WEST_FRONT
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_WEST;  
    (*work)[2] = GEOPOS_FRONT;
    (*work)[3] = GEOPOS_WEST_FRONT;  
    (*work)[4] = GEOPOS_NORTH_WEST;  
    (*work)[5] = GEOPOS_NORTH_FRONT;  
    (*work)[6] = GEOPOS_NORTH_WEST_FRONT;  
    m_periodic_vector_indices[GEOPOS_NORTH_WEST_FRONT] = *work;    
    work->clear();
    delete work; 

    // NORTH_WEST_BEHIND => NORTH, WEST, BEHIND, WEST_BEHIND, NORTH_WEST, 
    // NORTH_BEHIND, NORTH_WEST_BEHIND
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_WEST;  
    (*work)[2] = GEOPOS_BEHIND;
    (*work)[3] = GEOPOS_WEST_BEHIND;  
    (*work)[4] = GEOPOS_NORTH_WEST;  
    (*work)[5] = GEOPOS_NORTH_BEHIND;  
    (*work)[6] = GEOPOS_NORTH_WEST_BEHIND;  
    m_periodic_vector_indices[GEOPOS_NORTH_WEST_BEHIND] = *work;    
    work->clear();
    delete work; 
 
    // SOUTH => SOUTH
    work = new vector<int>(1,0);
    (*work)[0] = GEOPOS_SOUTH;
    m_periodic_vector_indices[GEOPOS_SOUTH] = *work;
    work->clear();
    delete work;
  
    // SOUTH_EAST => SOUTH, EAST, SOUTH_EAST
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_EAST;  
    (*work)[2] = GEOPOS_SOUTH_EAST;
    m_periodic_vector_indices[GEOPOS_SOUTH_EAST] = *work;    
    work->clear();
    delete work;
  
    // SOUTH_WEST => SOUTH, WEST, SOUTH_WEST
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_WEST;  
    (*work)[2] = GEOPOS_SOUTH_WEST;
    m_periodic_vector_indices[GEOPOS_SOUTH_WEST] = *work;    
    work->clear();
    delete work;    

    // SOUTH_FRONT => SOUTH, FRONT, SOUTH_FRONT
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_FRONT;  
    (*work)[2] = GEOPOS_SOUTH_FRONT;
    m_periodic_vector_indices[GEOPOS_SOUTH_FRONT] = *work;    
    work->clear();
    delete work;  

    // SOUTH_BEHIND => SOUTH, BEHIND, SOUTH_BEHIND
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_BEHIND;  
    (*work)[2] = GEOPOS_SOUTH_BEHIND;
    m_periodic_vector_indices[GEOPOS_SOUTH_BEHIND] = *work;    
    work->clear();
    delete work; 

    // SOUTH_EAST_FRONT => SOUTH, EAST, FRONT, EAST_FRONT, SOUTH_EAST, 
    // SOUTH_FRONT, SOUTH_EAST_FRONT
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_EAST;  
    (*work)[2] = GEOPOS_FRONT;
    (*work)[3] = GEOPOS_EAST_FRONT;  
    (*work)[4] = GEOPOS_SOUTH_EAST;  
    (*work)[5] = GEOPOS_SOUTH_FRONT;  
    (*work)[6] = GEOPOS_SOUTH_EAST_FRONT;  
    m_periodic_vector_indices[GEOPOS_SOUTH_EAST_FRONT] = *work;    
    work->clear();
    delete work; 

    // SOUTH_EAST_BEHIND => SOUTH, EAST, BEHIND, EAST_BEHIND, SOUTH_EAST, 
    // SOUTH_BEHIND, SOUTH_EAST_BEHIND
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_EAST;  
    (*work)[2] = GEOPOS_BEHIND;
    (*work)[3] = GEOPOS_EAST_BEHIND;  
    (*work)[4] = GEOPOS_SOUTH_EAST;  
    (*work)[5] = GEOPOS_SOUTH_BEHIND;  
    (*work)[6] = GEOPOS_SOUTH_EAST_BEHIND;  
    m_periodic_vector_indices[GEOPOS_SOUTH_EAST_BEHIND] = *work;    
    work->clear();
    delete work; 

    // SOUTH_WEST_FRONT => SOUTH, WEST, FRONT, WEST_FRONT, SOUTH_WEST, SOUTH_FRONT,
    // SOUTH_WEST_FRONT
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_WEST;  
    (*work)[2] = GEOPOS_FRONT;
    (*work)[3] = GEOPOS_WEST_FRONT;  
    (*work)[4] = GEOPOS_SOUTH_WEST;  
    (*work)[5] = GEOPOS_SOUTH_FRONT;  
    (*work)[6] = GEOPOS_SOUTH_WEST_FRONT;  
    m_periodic_vector_indices[GEOPOS_SOUTH_WEST_FRONT] = *work;    
    work->clear();
    delete work; 

    // SOUTH_WEST_BEHIND => SOUTH, WEST, BEHIND, WEST_BEHIND, SOUTH_WEST, 
    // SOUTH_BEHIND, SOUTH_WEST_BEHIND
    work = new vector<int>(7,0);
    (*work)[0] = GEOPOS_SOUTH;
    (*work)[1] = GEOPOS_WEST;  
    (*work)[2] = GEOPOS_BEHIND;
    (*work)[3] = GEOPOS_WEST_BEHIND;  
    (*work)[4] = GEOPOS_SOUTH_WEST;  
    (*work)[5] = GEOPOS_SOUTH_BEHIND;  
    (*work)[6] = GEOPOS_SOUTH_WEST_BEHIND;  
    m_periodic_vector_indices[GEOPOS_SOUTH_WEST_BEHIND] = *work;    
    work->clear();
    delete work; 
  
    // EAST => EAST
    work = new vector<int>(1,0);
    (*work)[0] = GEOPOS_EAST;
    m_periodic_vector_indices[GEOPOS_EAST] = *work;
    work->clear();
    delete work;  

    // WEST => WEST
    work = new vector<int>(1,0);
    (*work)[0] = GEOPOS_WEST;
    m_periodic_vector_indices[GEOPOS_WEST] = *work;
    work->clear();
    delete work;  

    // EAST_FRONT => EAST, FRONT, EAST_FRONT
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_EAST;
    (*work)[1] = GEOPOS_FRONT;  
    (*work)[2] = GEOPOS_EAST_FRONT;
    m_periodic_vector_indices[GEOPOS_EAST_FRONT] = *work;    
    work->clear();
    delete work;  

    // EAST_BEHIND => EAST, BEHIND, EAST_BEHIND
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_EAST;
    (*work)[1] = GEOPOS_BEHIND;  
    (*work)[2] = GEOPOS_EAST_BEHIND;
    m_periodic_vector_indices[GEOPOS_EAST_BEHIND] = *work;    
    work->clear();
    delete work; 

    // WEST_FRONT => WEST, FRONT, WEST_FRONT
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_WEST;
    (*work)[1] = GEOPOS_FRONT;  
    (*work)[2] = GEOPOS_WEST_FRONT;
    m_periodic_vector_indices[GEOPOS_WEST_FRONT] = *work;    
    work->clear();
    delete work;  

    // WEST_BEHIND => WEST, BEHIND, WEST_BEHIND
    work = new vector<int>(3,0);
    (*work)[0] = GEOPOS_WEST;
    (*work)[1] = GEOPOS_BEHIND;  
    (*work)[2] = GEOPOS_WEST_BEHIND;
    m_periodic_vector_indices[GEOPOS_WEST_BEHIND] = *work;    
    work->clear();
    delete work; 
   
    // FRONT => FRONT
    work = new vector<int>(1,0);
    (*work)[0] = GEOPOS_FRONT;
    m_periodic_vector_indices[GEOPOS_FRONT] = *work;
    work->clear();
    delete work;  

    // BEHIND => BEHIND
    work = new vector<int>(1,0);
    (*work)[0] = GEOPOS_BEHIND;
    m_periodic_vector_indices[GEOPOS_BEHIND] = *work;
    work->clear();
    delete work; 
    
    // Special cases for serial periodicity with a single cell in the main 
    // domain in the periodic direction(s)
    // EASTWEST => EAST, WEST
    work = new vector<int>(2,0);
    (*work)[0] = GEOPOS_EAST;
    (*work)[1] = GEOPOS_WEST;    
    m_periodic_vector_indices[GEOPOS_EASTWEST] = *work;
    work->clear();
    delete work;
    
    // NORTHSOUTH => NORTH, SOUTH
    work = new vector<int>(2,0);
    (*work)[0] = GEOPOS_NORTH;
    (*work)[1] = GEOPOS_SOUTH;    
    m_periodic_vector_indices[GEOPOS_NORTHSOUTH] = *work;
    work->clear();
    delete work;
    
    // GEOPOS_EASTWESTNORTHSOUTH => EAST, WEST, NORTH, NORTH_EAST, NORTH_WEST,
    // SOUTH_EAST, SOUTH_WEST
    work = new vector<int>(8,0);
    (*work)[0] = GEOPOS_EAST;
    (*work)[1] = GEOPOS_WEST;  
    (*work)[2] = GEOPOS_NORTH;
    (*work)[3] = GEOPOS_SOUTH;  
    (*work)[4] = GEOPOS_NORTH_EAST;  
    (*work)[5] = GEOPOS_NORTH_WEST;  
    (*work)[6] = GEOPOS_SOUTH_EAST; 
    (*work)[7] = GEOPOS_SOUTH_WEST;     
    m_periodic_vector_indices[GEOPOS_EASTWESTNORTHSOUTH] = *work;
    work->clear();
    delete work;               
  } 
}




// ----------------------------------------------------------------------------
// Sets local domain size
void App::set_local_domain_size( double lx_, double ly_, double lz_ )
{
  App::m_domain_local_size_X = lx_;
  App::m_domain_local_size_Y = ly_;
  App::m_domain_local_size_Z = lz_;
}




// ----------------------------------------------------------------------------
// Sets local domain origin
void App::set_local_domain_origin( int const* nprocsdir, int const* MPIcoords )
{
  m_domain_local_origin[0] = m_domain_global_origin[0] 
  	+ MPIcoords[0] * m_domain_local_size_X ;
  m_domain_local_origin[1] = m_domain_global_origin[1] 
  	+ MPIcoords[1] * m_domain_local_size_Y ;  
  m_domain_local_origin[2] = m_domain_global_origin[2] 
  	+ MPIcoords[2] * m_domain_local_size_Z ;  
}




// ----------------------------------------------------------------------------
// Gets local domain origin
void App::get_local_domain_origin( double& x, double& y, double& z ) 
{
  x = m_domain_local_origin[0] ;
  y = m_domain_local_origin[1] ;  
  z = m_domain_local_origin[2] ;
}




// ----------------------------------------------------------------------------
// Gets global domain origin
void App::get_origin( double& x, double& y, double& z ) 
{
  x = m_domain_global_origin[0] ;
  y = m_domain_global_origin[1] ;  
  z = m_domain_global_origin[2] ;
}




// ----------------------------------------------------------------------------
// Gets local domain size
void App::get_local_domain_size( double& lx, double& ly, double& lz ) 
{
  lx = m_domain_local_size_X ;
  ly = m_domain_local_size_Y ;  
  lz = m_domain_local_size_Z ;
}




// ----------------------------------------------------------------------------
// Gets global domain size
void App::get_size( double& lx, double& ly, double& lz ) 
{
  lx = m_domain_global_size_X ;
  ly = m_domain_global_size_Y ;  
  lz = m_domain_global_size_Z ;
}




// ----------------------------------------------------------------------------
// Returns whether a point belongs to the global doamin
bool App::isInDomain( Point3 const* position )
{
  bool isIn = true;
  
  if ( (*position)[0] < m_domain_global_origin[0] 
  	|| (*position)[0] > m_domain_global_origin[0] + m_domain_global_size_X 
  	|| (*position)[1] < m_domain_global_origin[1] 
	|| (*position)[1] > m_domain_global_origin[1] + m_domain_global_size_Y
  	|| (*position)[2] < m_domain_global_origin[2] 
	|| (*position)[2] > m_domain_global_origin[2] + m_domain_global_size_Z )
    isIn = false;
  
  return ( isIn );
}  




// ----------------------------------------------------------------------------
// Returns whether a point belongs to the local domain
bool App::isInLocalDomain( Point3 const* position )
{
  bool isIn = true;
  
  if ( (*position)[0] < m_domain_local_origin[0] 
  	|| (*position)[0] > m_domain_local_origin[0] + m_domain_local_size_X 
  	|| (*position)[1] < m_domain_local_origin[1] 
	|| (*position)[1] > m_domain_local_origin[1] + m_domain_local_size_Y
  	|| (*position)[2] < m_domain_local_origin[2] 
	|| (*position)[2] > m_domain_local_origin[2] + m_domain_local_size_Z ) 
    isIn = false;
  
  return ( isIn );
} 




// ----------------------------------------------------------------------------
// Tells if an app corresponding to the argument exists or not
bool App::isName( string const& name_ )
{
  return ( m_name == name_ ) ;
}




// ----------------------------------------------------------------------------
// Sets the name of the application
void App::setName( string const& name_ ) 
{ 
  m_name = name_; 
}




// ----------------------------------------------------------------------------
// Returns the name of the application
string App::getName() const 
{ 
  return ( m_name ); 
} 




// ----------------------------------------------------------------------------
// Writes the domain features in an output stream
void App::output_domain_features( ostream& output, string const& oshift )
{
  output << oshift << "Global domain size  = " << m_domain_global_size_X 
  	<< " x " << m_domain_global_size_Y 
	<< " x " << m_domain_global_size_Z << endl;
  output << oshift << "Periodicity = ";
  if ( !m_domain_global_periodic )
    output << "No" << endl;
  else
  {
    if ( m_domain_global_periodicity[X] ) output << "X ";
    if ( m_domain_global_periodicity[Y] ) output << "Y ";    
    if ( m_domain_global_periodicity[Z] ) output << "Z ";
    output << endl;

    output << oshift << "Global periodic vectors" << endl;
    for (size_t i=0;i<27;++i)
      output << oshift << GrainsExec::m_shift3 << 
      	Cell::getGeoPositionName( int(i) ) << " = " 
	<< m_domain_global_periodic_vectors[i] << endl; 

    output << oshift << "Periodic vectors per geographic position" << endl;

    // Single periodicity
    // West East
    if ( m_domain_global_periodicity[X] )
    {
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_WEST );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_EAST );      
    }
    
    if ( m_domain_global_periodicity[Y] )
    {
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_SOUTH );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_NORTH );  
    }
    
    if ( m_domain_global_periodicity[Z] )
    {
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_FRONT );  
    }
    
    // Double periodicity    
    // South North West East    
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Y] )
    {
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_SOUTH_WEST );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_SOUTH_EAST );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_NORTH_WEST );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_NORTH_EAST );
    }
    
    // South North BEHIND FRONT
    if ( m_domain_global_periodicity[Y] && m_domain_global_periodicity[Z] )
    {
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_SOUTH_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_SOUTH_FRONT );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_NORTH_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_NORTH_FRONT );
    }      

    // West East BEHIND FRONT
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Z] )
    {
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_WEST_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_WEST_FRONT );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_EAST_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, GEOPOS_EAST_FRONT );
    }   
    
    // Triple periodicity
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Y] 
    	&& m_domain_global_periodicity[Z] )
    {
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_SOUTH_WEST_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_SOUTH_WEST_FRONT );
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_NORTH_WEST_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_NORTH_WEST_FRONT );
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_SOUTH_EAST_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_SOUTH_EAST_FRONT );
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_NORTH_EAST_BEHIND );
      output_periodic_vectors_per_geopos( output, oshift, 
      	GEOPOS_NORTH_EAST_FRONT );
    }                
  }       	
  output << oshift << "Local domain size  = " << m_domain_local_size_X 
  	<< " x " << m_domain_local_size_Y 
	<< " x " << m_domain_local_size_Z << endl;
  output << oshift << "Global origin = " << m_domain_global_origin << endl;
  output << oshift << "Local origin = " << m_domain_local_origin << endl;  
}    




// ----------------------------------------------------------------------------
// Writes all periodic vectors for a given geographic position in an output 
// stream
void App::output_periodic_vectors_per_geopos( ostream& output, 
	string const& oshift, GeoPosition geoloc_ )
{
  output << oshift << GrainsExec::m_shift3 << 
      	Cell::getGeoPositionName( geoloc_ ) << endl;
  for ( size_t i=0;i<m_periodic_vector_indices[geoloc_].size();++i)
    output << oshift << GrainsExec::m_shift6 << 
	m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc_][i]] << endl;  
}




// ----------------------------------------------------------------------------
// Returns whether the domain is periodic in a direction
bool App::isPeriodic( size_t dir )
{
  return ( m_domain_global_periodicity[dir] );
}
