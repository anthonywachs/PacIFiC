#include "App.hh"
#include "AllComponents.hh"
#include "GrainsExec.hh"

Vector3 App::m_domain_global_size;
Vector3 App::m_domain_local_size;
Point3 App::m_domain_global_origin;
Point3 App::m_domain_local_origin;
Point3 App::m_domain_global_max;
Point3 App::m_domain_local_max;
vector<bool> App::m_domain_global_periodicity;
bool App::m_domain_global_periodic = false;
vector<Vector3> App::m_domain_global_periodic_vectors;
vector< vector<int> > App::m_periodic_vector_indices;


// ----------------------------------------------------------------------------
// Default constructor
App::App()
{}




// ----------------------------------------------------------------------------
// Destructor
App::~App()
{}




// ----------------------------------------------------------------------------
// Sets domain dimensions
void App::set_dimensions( double xmax, double ymax, double zmax,
  	double ox, double oy, double oz )
{
  m_domain_global_origin[X] = ox;
  m_domain_global_origin[Y] = oy;
  m_domain_global_origin[Z] = oz;
  
  m_domain_global_max[X] = xmax;
  m_domain_global_max[Y] = ymax;  
  m_domain_global_max[Z] = zmax;    
   
  for (size_t i=0;i<3;++i)
  {
    m_domain_local_origin[i] = m_domain_global_origin[i];
    m_domain_local_max[i] = m_domain_global_max[i];
    m_domain_global_size[i] = m_domain_global_max[i] 
    	- m_domain_global_origin[i];
    m_domain_local_size[i] = m_domain_global_size[i];
  }

  assert( App::m_domain_global_size[X] >= 0. 
  	&& App::m_domain_global_size[Y] >= 0.
 	&& App::m_domain_global_size[Z] >= 0. );   
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
      m_domain_global_periodic_vectors.push_back( Vector3Null );
      
    // Single periodicity
    // West East
    if ( m_domain_global_periodicity[X] )
    {
      m_domain_global_periodic_vectors[GEOPOS_WEST][X] = 
      	m_domain_global_size[X] ;  
      m_domain_global_periodic_vectors[GEOPOS_EAST][X] = 
      	- m_domain_global_size[X] ;
    }

    // South North
    if ( m_domain_global_periodicity[Y] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH][Y] = 
      	m_domain_global_size[Y] ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH][Y] = 
      	- m_domain_global_size[Y] ;
    }

    // Behind Front
    if ( m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_BEHIND][Z] = 
      	m_domain_global_size[Z] ;
      m_domain_global_periodic_vectors[GEOPOS_FRONT][Z] = 
      	- m_domain_global_size[Z] ;
    }


    // Double periodicity    
    // South North West East
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Y] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST][X] = 
      	m_domain_global_size[X] ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST][Y] = 
      	m_domain_global_size[Y] ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST][X] = 
      	- m_domain_global_size[X] ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST][Y] = 
      	m_domain_global_size[Y] ;        
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST][X] = 
      	m_domain_global_size[X] ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST][Y] = 
      	- m_domain_global_size[Y] ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST][X] = 
      	- m_domain_global_size[X] ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST][Y] = 
      	- m_domain_global_size[Y] ;
    }
    
    // South North Behind Front
    if ( m_domain_global_periodicity[Y] && m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_BEHIND][Y] = 
      	m_domain_global_size[Y] ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_BEHIND][Z] = 
      	m_domain_global_size[Z] ;	     
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_FRONT][Y] = 
      	m_domain_global_size[Y] ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_FRONT][Z] = 
      	- m_domain_global_size[Z] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_BEHIND][Y] = 
      	- m_domain_global_size[Y] ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_BEHIND][Z] = 
      	m_domain_global_size[Z] ;     
      m_domain_global_periodic_vectors[GEOPOS_NORTH_FRONT][Y] = 
      	- m_domain_global_size[Y] ;    
      m_domain_global_periodic_vectors[GEOPOS_NORTH_FRONT][Z] = 
      	- m_domain_global_size[Z] ;
    }      

    // West East Behind Front
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_WEST_BEHIND][X] = 
      	m_domain_global_size[X] ;    
      m_domain_global_periodic_vectors[GEOPOS_WEST_BEHIND][Z] = 
      	m_domain_global_size[Z] ;    
      m_domain_global_periodic_vectors[GEOPOS_WEST_FRONT][X] = 
      	m_domain_global_size[X] ;    
      m_domain_global_periodic_vectors[GEOPOS_WEST_FRONT][Z] = 
      	- m_domain_global_size[Z] ;   
      m_domain_global_periodic_vectors[GEOPOS_EAST_BEHIND][X] = 
      	- m_domain_global_size[X] ;    
      m_domain_global_periodic_vectors[GEOPOS_EAST_BEHIND][Z] = 
      	m_domain_global_size[Z] ;
      m_domain_global_periodic_vectors[GEOPOS_EAST_FRONT][X] = 
      	- m_domain_global_size[X] ;    
      m_domain_global_periodic_vectors[GEOPOS_EAST_FRONT][Z] = 
      	- m_domain_global_size[Z] ; 
    }   

    
    // Triple periodicity
    if ( m_domain_global_periodicity[X] && m_domain_global_periodicity[Y] 
    	&& m_domain_global_periodicity[Z] )
    {
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_BEHIND][X] = 
      	m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_BEHIND][Y] = 
      	m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_BEHIND][Z] = 
      	m_domain_global_size[Z] ;  		     
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_FRONT][X] = 
      	m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_FRONT][Y] = 
      	m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_FRONT][Z] = 
      	- m_domain_global_size[Z] ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_BEHIND][X] = 
      	m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_BEHIND][Y] = 
      	- m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_BEHIND][Z] = 
      	m_domain_global_size[Z] ;
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_FRONT][X] = 
      	m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_FRONT][Y] = 
      	- m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_FRONT][Z] = 
      	- m_domain_global_size[Z] ;
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_BEHIND][X] = 
      	- m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_BEHIND][Y] = 
      	m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_BEHIND][Z] = 
      	m_domain_global_size[Z] ;  
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_FRONT][X] = 
      	- m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_FRONT][Y] = 
      	m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_FRONT][Z] = 
      	- m_domain_global_size[Z] ;  
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_BEHIND][X] = 
      	- m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_BEHIND][Y] = 
      	- m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_BEHIND][Z] = 
      	m_domain_global_size[Z] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_FRONT][X] = 
      	- m_domain_global_size[X] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_FRONT][Y] = 
      	- m_domain_global_size[Y] ; 
      m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_FRONT][Z] = 
      	- m_domain_global_size[Z] ;
    }            
  }
  
  if ( m_periodic_vector_indices.empty() )
  {
    m_periodic_vector_indices.reserve(30);
    vector<int> emptyVECINT;
    list<int> indices;    
    for (int i=0;i<30;++i) 
      m_periodic_vector_indices.push_back( emptyVECINT );
  
    // NORTH => NORTH
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH] = w;
      indices.clear();
    }  
  
    // NORTH_EAST => NORTH, EAST, NORTH_EAST
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_EAST );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_EAST] = w;
      indices.clear();
    }        
  
    // NORTH_WEST => NORTH, WEST, NORTH_WEST
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_WEST );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_WEST] = w;
      indices.clear();
    }  

    // NORTH_FRONT => NORTH, FRONT, NORTH_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_FRONT );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_FRONT] = w;
      indices.clear();
    } 

    // NORTH_BEHIND => NORTH, BEHIND, NORTH_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_BEHIND );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_BEHIND] = w;
      indices.clear();
    }

    // NORTH_EAST_FRONT => NORTH, EAST, FRONT, EAST_FRONT, NORTH_EAST, 
    // NORTH_FRONT, NORTH_EAST_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_EAST_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_FRONT] 
    	!= Vector3Null ) indices.push_back( GEOPOS_NORTH_EAST_FRONT ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_EAST_FRONT] = w;
      indices.clear();
    }

    // NORTH_EAST_BEHIND => NORTH, EAST, BEHIND, EAST_BEHIND, NORTH_EAST, 
    // NORTH_BEHIND, NORTH_EAST_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_EAST_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST_BEHIND] 
    	!= Vector3Null ) indices.push_back( GEOPOS_NORTH_EAST_BEHIND ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_EAST_BEHIND] = w;
      indices.clear();
    }

    // NORTH_WEST_FRONT => NORTH, WEST, FRONT, WEST_FRONT, NORTH_WEST, 
    // NORTH_FRONT, NORTH_WEST_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_WEST_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_FRONT] 
    	!= Vector3Null ) indices.push_back( GEOPOS_NORTH_WEST_FRONT ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_WEST_FRONT] = w;
      indices.clear();
    }

    // NORTH_WEST_BEHIND => NORTH, WEST, BEHIND, WEST_BEHIND, NORTH_WEST, 
    // NORTH_BEHIND, NORTH_WEST_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_WEST_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST_BEHIND] 
    	!= Vector3Null ) indices.push_back( GEOPOS_NORTH_WEST_BEHIND ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTH_WEST_BEHIND] = w;
      indices.clear();
    } 
 
    // SOUTH => SOUTH
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH] = w;
      indices.clear();
    }        
  
    // SOUTH_EAST => SOUTH, EAST, SOUTH_EAST
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_EAST );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_EAST] = w;
      indices.clear();
    }
  
    // SOUTH_WEST => SOUTH, WEST, SOUTH_WEST
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_WEST );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_WEST] = w;
      indices.clear();
    }   

    // SOUTH_FRONT => SOUTH, FRONT, SOUTH_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_FRONT );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_FRONT] = w;
      indices.clear();
    }

    // SOUTH_BEHIND => SOUTH, BEHIND, SOUTH_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_BEHIND );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_BEHIND] = w;
      indices.clear();
    }

    // SOUTH_EAST_FRONT => SOUTH, EAST, FRONT, EAST_FRONT, SOUTH_EAST, 
    // SOUTH_FRONT, SOUTH_EAST_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_EAST_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_FRONT] 
    	!= Vector3Null ) indices.push_back( GEOPOS_SOUTH_EAST_FRONT ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_EAST_FRONT] = w;
      indices.clear();
    } 

    // SOUTH_EAST_BEHIND => SOUTH, EAST, BEHIND, EAST_BEHIND, SOUTH_EAST, 
    // SOUTH_BEHIND, SOUTH_EAST_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_EAST_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST_BEHIND] 
    	!= Vector3Null ) indices.push_back( GEOPOS_SOUTH_EAST_BEHIND ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_EAST_BEHIND] = w;
      indices.clear();
    }

    // SOUTH_WEST_FRONT => SOUTH, WEST, FRONT, WEST_FRONT, SOUTH_WEST, 
    // SOUTH_FRONT, SOUTH_WEST_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_WEST_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_FRONT ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_FRONT] 
    	!= Vector3Null ) indices.push_back( GEOPOS_SOUTH_WEST_FRONT ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_WEST_FRONT] = w;
      indices.clear();
    }

    // SOUTH_WEST_BEHIND => SOUTH, WEST, BEHIND, WEST_BEHIND, SOUTH_WEST, 
    // SOUTH_BEHIND, SOUTH_WEST_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_WEST_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_BEHIND ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST_BEHIND] 
    	!= Vector3Null ) indices.push_back( GEOPOS_SOUTH_WEST_BEHIND ); 
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_SOUTH_WEST_BEHIND] = w;
      indices.clear();
    }
  
    // EAST => EAST
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_EAST] = w;
      indices.clear();
    } 

    // WEST => WEST
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_WEST] = w;
      indices.clear();
    } 

    // EAST_FRONT => EAST, FRONT, EAST_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_EAST_FRONT );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_EAST_FRONT] = w;
      indices.clear();
    }  

    // EAST_BEHIND => EAST, BEHIND, EAST_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_EAST_BEHIND );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_EAST_BEHIND] = w;
      indices.clear();
    } 

    // WEST_FRONT => WEST, FRONT, WEST_FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_WEST_FRONT );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_WEST_FRONT] = w;
      indices.clear();
    }  

    // WEST_BEHIND => WEST, BEHIND, WEST_BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_WEST_BEHIND );      
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_WEST_BEHIND] = w;
      indices.clear();
    }
   
    // FRONT => FRONT
    if ( m_domain_global_periodic_vectors[GEOPOS_FRONT] != Vector3Null )
      indices.push_back( GEOPOS_FRONT );
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_FRONT] = w;
      indices.clear();
    }   

    // BEHIND => BEHIND
    if ( m_domain_global_periodic_vectors[GEOPOS_BEHIND] != Vector3Null )
      indices.push_back( GEOPOS_BEHIND );
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_BEHIND] = w;
      indices.clear();
    }
 
    
    // Special cases for serial periodicity with a single cell in the main 
    // domain in the periodic direction(s)
    // EASTWEST => EAST, WEST
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );          
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_EASTWEST] = w;
      indices.clear();
    }
    
    // NORTHSOUTH => NORTH, SOUTH
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );          
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_NORTHSOUTH] = w;
      indices.clear();
    }
    
    // GEOPOS_EASTWESTNORTHSOUTH => EAST, WEST, NORTH, NORTH_EAST, NORTH_WEST,
    // SOUTH_EAST, SOUTH_WEST
    if ( m_domain_global_periodic_vectors[GEOPOS_EAST] != Vector3Null )
      indices.push_back( GEOPOS_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_WEST] != Vector3Null )
      indices.push_back( GEOPOS_WEST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH] != Vector3Null )
      indices.push_back( GEOPOS_NORTH ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_NORTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_NORTH_WEST ); 
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_EAST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_EAST );    
    if ( m_domain_global_periodic_vectors[GEOPOS_SOUTH_WEST] != Vector3Null )
      indices.push_back( GEOPOS_SOUTH_WEST ); 
    if ( indices.size() )
    {
      vector<int> w( begin(indices), end(indices) );
      m_periodic_vector_indices[GEOPOS_EASTWESTNORTHSOUTH] = w;
      indices.clear();
    }               
  } 
}




// ----------------------------------------------------------------------------
// Sets local domain size
void App::set_local_domain_size( double lx_, double ly_, double lz_ )
{
  m_domain_local_size[X] = lx_;
  m_domain_local_size[Y] = ly_;
  m_domain_local_size[Z] = lz_;
}




// ----------------------------------------------------------------------------
// Sets local domain origin
void App::set_local_domain_origin( int const* nprocsdir, int const* MPIcoords )
{
  for (size_t i=0;i<3;++i)
  {
    m_domain_local_origin[i] = m_domain_global_origin[i] 
  	+ MPIcoords[i] * m_domain_local_size[i] ;
    m_domain_local_max[i] = m_domain_local_origin[i] + m_domain_local_size[i];
  } 
}




// ----------------------------------------------------------------------------
// Gets local domain origin
void App::get_local_domain_origin( double& x, double& y, double& z ) 
{
  x = m_domain_local_origin[X] ;
  y = m_domain_local_origin[Y] ;  
  z = m_domain_local_origin[Z] ;
}




// ----------------------------------------------------------------------------
// Gets global domain origin
void App::get_origin( double& x, double& y, double& z ) 
{
  x = m_domain_global_origin[X] ;
  y = m_domain_global_origin[Y] ;  
  z = m_domain_global_origin[Z] ;
}




// ----------------------------------------------------------------------------
// Gets local domain size
void App::get_local_domain_size( double& lx, double& ly, double& lz ) 
{
  lx = m_domain_local_size[X] ;
  ly = m_domain_local_size[Y] ;  
  lz = m_domain_local_size[Z] ;
}




// ----------------------------------------------------------------------------
// Gets global domain size
void App::get_size( double& lx, double& ly, double& lz ) 
{
  lx = m_domain_global_size[X] ;
  ly = m_domain_global_size[Y] ;  
  lz = m_domain_global_size[Z] ;
}




// ----------------------------------------------------------------------------
// Returns whether a point belongs to the global doamin
bool App::isInDomain( Point3 const* position )
{
  bool isIn = true;
  
  if ( (*position)[X] < m_domain_global_origin[X] 
  	|| (*position)[X] > m_domain_global_origin[X] 
		+ m_domain_global_size[X] 
  	|| (*position)[Y] < m_domain_global_origin[Y] 
	|| (*position)[Y] > m_domain_global_origin[Y] 
		+ m_domain_global_size[Y]
  	|| (*position)[Z] < m_domain_global_origin[Z] 
	|| (*position)[Z] > m_domain_global_origin[Z] 
		+ m_domain_global_size[Z] )
    isIn = false;
  
  return ( isIn );
} 



 
// ----------------------------------------------------------------------------
// Returns whether a point belongs to the global domain in a given direction
bool App::isInDomain( Point3 const* position, size_t const& dir )
{
  bool isIn = true;
  
  switch ( dir )
  {
    case 0:
      if ( (*position)[X] < m_domain_global_origin[X] 
  	|| (*position)[X] > m_domain_global_origin[X] 
		+ m_domain_global_size[X] )
      isIn = false; 
      break;
      
    case 1:
      if ( (*position)[Y] < m_domain_global_origin[Y] 
	|| (*position)[Y] > m_domain_global_origin[Y] 
		+ m_domain_global_size[Y] )
      isIn = false; 
      break;
      
    default:
      if ( (*position)[Z] < m_domain_global_origin[Z] 
	|| (*position)[Z] > m_domain_global_origin[Z] 
		+ m_domain_global_size[Z] )
      isIn = false;                        
      break;  
  }
  
  return ( isIn );
}  



// ----------------------------------------------------------------------------
// Returns whether a point belongs to the local domain
bool App::isInLocalDomain( Point3 const* position )
{
  bool isIn = true;
  
  if ( (*position)[X] < m_domain_local_origin[X] 
  	|| (*position)[X] > m_domain_local_origin[X] + m_domain_local_size[X] 
  	|| (*position)[Y] < m_domain_local_origin[Y] 
	|| (*position)[Y] > m_domain_local_origin[Y] + m_domain_local_size[Y]
  	|| (*position)[Z] < m_domain_local_origin[Z] 
	|| (*position)[Z] > m_domain_local_origin[Z] + m_domain_local_size[Z] ) 
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
  output << oshift << "Global domain size  = " << m_domain_global_size[X] 
  	<< " x " << m_domain_global_size[Y] 
	<< " x " << m_domain_global_size[Z] << endl;
  output << oshift << "Periodicity = ";
  if ( !m_domain_global_periodic )
    output << "No" << endl;
  else
  {
    if ( m_domain_global_periodicity[X] ) output << "X ";
    if ( m_domain_global_periodicity[Y] ) output << "Y ";    
    if ( m_domain_global_periodicity[Z] ) output << "Z ";
    output << endl;
    
    if ( GrainsExec::m_partialPer_is_active )
    {
      PartialPeriodicity const* pp = GrainsExec::getPartialPeriodicity(); 
      output << oshift << "Partial Periodicity = " << 
      	( pp->dir == X ? "X" : ( pp->dir == Y ? "Y" : "Z" ) ) << 
	( pp->comp == LLO_LARGER ? " > " : " < " ) <<
	pp->limit << endl;        
    }

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
  output << oshift << "Local domain size  = " << m_domain_local_size[X] 
  	<< " x " << m_domain_local_size[Y] 
	<< " x " << m_domain_local_size[Z] << endl;
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
