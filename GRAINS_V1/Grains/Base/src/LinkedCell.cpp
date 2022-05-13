#include "GrainsMPIWrapper.hh"
#include "AllComponents.hh"
#include "MPINeighbors.hh"
#include "LinkedCell.hh"  
#include "Cell.hh"
#include "RigidBodyWithCrust.hh"
#include "Box.hh"
#include "Grains.hh"
#include "GrainsExec.hh"
#include <algorithm>


// ----------------------------------------------------------------------------
// Default constructor
LinkedCell::LinkedCell() 
  : AppCollision()
  , m_nb( 0 )
  , m_nbi( 0 )
  , m_nbj( 0 )
  , m_nbk( 0 )
  , m_cellsize_X( 0. )
  , m_cellsize_Y( 0. )
  , m_cellsize_Z( 0. )
  , m_LC_local_xmin( 0. )
  , m_LC_local_ymin( 0. )
  , m_LC_local_zmin( 0. )
  , m_LC_local_xmax( 0. )
  , m_LC_local_ymax( 0. )
  , m_LC_local_zmax( 0. )
  , m_extendedBBox( NULL )
{}




// ----------------------------------------------------------------------------
// Destructor
LinkedCell::~LinkedCell()
{
  vector<Cell*>::iterator cell_;
  for (cell_=m_allcells.begin(); cell_!=m_allcells.end(); cell_++)
    delete *cell_;
  delete m_extendedBBox;
}




// ----------------------------------------------------------------------------
// Sets the linked cell grid in serial mode
void LinkedCell::set( double cellsize_, string const& oshift )
{
  m_LC_global_origin = m_domain_global_origin;
  m_LC_local_origin = m_domain_local_origin;

  // Number of cells and cell edge length in each direction and 
  // Default in 1 unique cell if cellsize_ is zero 
  if ( cellsize_ > 1.e-10 ) 
  {
    m_nbi = (int)( ( App::m_domain_local_size_X + EPSILON ) / cellsize_);
    m_cellsize_X = App::m_domain_local_size_X / m_nbi ;
    m_nbj = (int)( ( App::m_domain_local_size_Y + EPSILON ) / cellsize_);
    m_cellsize_Y = App::m_domain_local_size_Y / m_nbj ;    
    m_nbk = (int)( ( App::m_domain_local_size_Z + EPSILON ) / cellsize_);
    m_cellsize_Z = App::m_domain_local_size_Z / m_nbk ; 
    if ( !m_nbk ) 
    {
      m_nbk = 1;
      m_cellsize_Z = m_cellsize_Y;
    }
    
    // Periodicity
    if ( m_domain_global_periodicity[X] ) 
    {
      m_LC_global_origin.Move( - m_cellsize_X, 0., 0. );
      m_LC_local_origin.Move( - m_cellsize_X, 0., 0. );
      m_nbi += 2;
    }
    if ( m_domain_global_periodicity[Y] ) 
    {
      m_LC_global_origin.Move( 0., - m_cellsize_Y, 0. );
      m_LC_local_origin.Move( 0., - m_cellsize_Y, 0. ); 
      m_nbj += 2;     
    }    
    if ( m_domain_global_periodicity[Z] ) 
    {
      m_LC_global_origin.Move( 0., 0., - m_cellsize_Z );
      m_LC_local_origin.Move( 0., 0., - m_cellsize_Z ); 
      m_nbk += 2;     
    }      	                    
  } 
  else  m_nbi = m_nbj = m_nbk = 1;

  m_nb = m_nbi * m_nbj * m_nbk;

  cout << oshift << "Number of cells = " << m_nbi << " " << m_nbj << " " 
  	<< m_nbk << " = " << m_nbi * m_nbj * m_nbk << endl;
  cout << oshift << "Cell size  = " << m_cellsize_X << " x " << m_cellsize_Y <<
  	" x " << m_cellsize_Z << endl;
  cout << oshift << "Global origin = " << m_LC_global_origin << endl;
  cout << oshift << "Local origin = " << m_LC_global_origin << endl;
    
  m_LC_local_xmin = m_LC_local_origin[0];
  m_LC_local_ymin = m_LC_local_origin[1]; 
  m_LC_local_zmin = m_LC_local_origin[2];
  m_LC_local_xmax = m_LC_local_xmin + m_nbi * m_cellsize_X;
  m_LC_local_ymax = m_LC_local_ymin + m_nbj * m_cellsize_Y;
  m_LC_local_zmax = m_LC_local_zmin + m_nbk * m_cellsize_Z;
  m_extendedBBox = new BBox(
  	Point3( m_LC_local_xmin - 0.5 * m_cellsize_X, 
		m_LC_local_ymin - 0.5 * m_cellsize_Y,
		m_LC_local_zmin - 0.5 * m_cellsize_Z ),
	Point3( m_LC_local_xmax + 0.5 * m_cellsize_X,
		m_LC_local_ymax + 0.5 * m_cellsize_Y,
		m_LC_local_zmax + 0.5 * m_cellsize_Z ) );
     	
  cout << oshift << "Local size  = " << m_LC_local_xmax - m_LC_local_xmin 
  	<< " x " << m_LC_local_ymax - m_LC_local_ymin 
	<< " x " << m_LC_local_zmax - m_LC_local_zmin << endl;
  
  // Cells construction 
  m_allcells.reserve( m_nb );
  Cell::setNbCellsPerDirection( m_nbi, m_nbj, m_nbk );
  for (int j=0; j<m_nbj; j++) 
    for (int k=0; k<m_nbk; k++)
      for (int i=0; i<m_nbi; i++)
        m_allcells.push_back( new Cell( getCellNumber( i, j, k ), 
                i, j, k, m_LC_local_origin, 
		m_cellsize_X, m_cellsize_Y, m_cellsize_Z, 
                m_LC_local_xmax, m_LC_local_ymax, m_LC_local_zmax ) );

  // Sets the the list of neighboring cells over which broad phase
  // contact detection is performed
  setCellContactNeighborhood(); 
  
  // If periodic, assign tag and geographic position to the cells
  if ( m_domain_global_periodicity[X] ) 
  { 
    for (int j=0; j<m_nbj; j++) 
      for (int k=0; k<m_nbk; k++) 
      {
        // Set tag = 2 to 1st and last row
	m_allcells[getCellNumber( 0, j, k )]->m_tag = 2;
	m_allcells[getCellNumber( m_nbi - 1, j, k )]->m_tag = 2;
	
	// Set tag = 1 and geopos = GEOPOS_WEST in 2nd row
	m_allcells[getCellNumber( 1, j, k )]->m_tag = 1;
	m_allcells[getCellNumber( 1, j, k )]->m_GeoPosCell = GEOPOS_WEST;
	
	// Set tag = 1 and geopos = GEOPOS_EAST in penultimate row
	m_allcells[getCellNumber( m_nbi - 2, j, k )]->m_tag = 1;
	m_allcells[getCellNumber( m_nbi - 2, j, k )]->m_GeoPosCell = 
		GEOPOS_EAST;
      }	     
  }
  
  if ( m_domain_global_periodicity[Y] ) 
  { 
    for (int i=0; i<m_nbi; i++) 
      for (int k=0; k<m_nbk; k++) 
      {
        // Set tag = 2 and geopos = GEOPOS_NONE to 1st and last row
	m_allcells[getCellNumber( i, 0, k )]->m_tag = 2;
	m_allcells[getCellNumber( i, m_nbj - 1, k )]->m_tag = 2;
	m_allcells[getCellNumber( i, 0, k )]->m_GeoPosCell = GEOPOS_NONE ;
	m_allcells[getCellNumber( i, m_nbj - 1, k )]->m_GeoPosCell = 
		GEOPOS_NONE;	
	
	// Set tag = 1 and geopos += GEOPOS_SOUTH in 2nd row
	if ( m_allcells[getCellNumber( i, 1, k )]->m_tag != 2 )
	{
	  m_allcells[getCellNumber( i, 1, k )]->m_tag = 1;
	  switch( int(m_allcells[getCellNumber( i, 1, k )]->m_GeoPosCell) )
	  {
	    case GEOPOS_NONE:
	      m_allcells[getCellNumber( i, 1, k )]->m_GeoPosCell = GEOPOS_SOUTH;
	      break;
	    case GEOPOS_WEST:
	      m_allcells[getCellNumber( i, 1, k )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_WEST;
	      break;
	    case GEOPOS_EAST:
	      m_allcells[getCellNumber( i, 1, k )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_EAST;
	      break;	      
	  }
	}
	
	// Set tag = 1 and geopos += GEOPOS_NORTH in penultimate row
	if ( m_allcells[getCellNumber( i, m_nbj - 2, k )]->m_tag != 2 )
	{
	  m_allcells[getCellNumber( i, m_nbj - 2, k )]->m_tag = 1;
	  switch( int(m_allcells[getCellNumber( i, m_nbj - 2, k )]
	  	->m_GeoPosCell) )
	  {
	    case GEOPOS_NONE:
	      m_allcells[getCellNumber( i, m_nbj - 2, k )]->m_GeoPosCell = 
	      	GEOPOS_NORTH;
	      break;
	    case GEOPOS_WEST:
	      m_allcells[getCellNumber( i, m_nbj - 2, k )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_WEST;
	      break;
	    case GEOPOS_EAST:
	      m_allcells[getCellNumber( i, m_nbj - 2, k )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_EAST;
	      break;	      
	  }
	}
      }	     
  }
  
  if ( m_domain_global_periodicity[Z] ) 
  {
    for (int i=0; i<m_nbi; i++) 
      for (int j=0; j<m_nbj; j++) 
      {
        // Set tag = 2 and geopos = GEOPOS_NONE to 1st and last row
	m_allcells[getCellNumber( i, j, 0 )]->m_tag = 2;
	m_allcells[getCellNumber( i, j, m_nbk - 1 )]->m_tag = 2;
	m_allcells[getCellNumber( i, j, 0 )]->m_GeoPosCell = GEOPOS_NONE ;
	m_allcells[getCellNumber( i, j, m_nbk - 1 )]->m_GeoPosCell = 
		GEOPOS_NONE;      

	// Set tag = 1 and geopos += GEOPOS_BEHIND in 2nd row
	if ( m_allcells[getCellNumber( i, j, 1 )]->m_tag != 2 )
	{
	  m_allcells[getCellNumber( i, j, 1 )]->m_tag = 1;
	  switch( int(m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell) )
	  {
	    case GEOPOS_NONE:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_BEHIND;
	      break;
	    case GEOPOS_WEST:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_WEST_BEHIND;
	      break;
	    case GEOPOS_EAST:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_EAST_BEHIND;
	      break;
	    case GEOPOS_SOUTH:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_BEHIND;
	      break;	      
	    case GEOPOS_NORTH:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_BEHIND;
	      break;	      
	    case GEOPOS_SOUTH_WEST:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_WEST_BEHIND;
	      break;
	    case GEOPOS_SOUTH_EAST:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_EAST_BEHIND;
	      break;
	    case GEOPOS_NORTH_WEST:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_WEST_BEHIND;
	      break;
	    case GEOPOS_NORTH_EAST:
	      m_allcells[getCellNumber( i, j, 1 )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_EAST_BEHIND;
	      break;	      	      	      	      
	  }
	}
	
	// Set tag = 1 and geopos += GEOPOS_FRONT in penultimate row
	if ( m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_tag != 2 )
	{
	  m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_tag = 1;
	  switch( int(m_allcells[getCellNumber( i, j, m_nbk - 2 )]
	  	->m_GeoPosCell) )
	  {
	    case GEOPOS_NONE:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_FRONT;
	      break;
	    case GEOPOS_WEST:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_WEST_FRONT;
	      break;
	    case GEOPOS_EAST:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_EAST_FRONT;
	      break;
	    case GEOPOS_SOUTH:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_FRONT;
	      break;	      
	    case GEOPOS_NORTH:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_FRONT;
	      break;	      
	    case GEOPOS_SOUTH_WEST:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_WEST_FRONT;
	      break;
	    case GEOPOS_SOUTH_EAST:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_SOUTH_EAST_FRONT;
	      break;
	    case GEOPOS_NORTH_WEST:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_WEST_FRONT;
	      break;
	    case GEOPOS_NORTH_EAST:
	      m_allcells[getCellNumber( i, j, m_nbk - 2 )]->m_GeoPosCell = 
	      	GEOPOS_NORTH_EAST_FRONT;
	      break;	      	      	      	      
	  }
	}	
      }  
  }
  
  // list of buffer cells (i.e. tag = 1)
  for (int j=0; j<m_nbj; j++) 
    for (int k=0; k<m_nbk; k++)
      for (int i=0; i<m_nbi; i++)
        if ( m_allcells[getCellNumber( i, j, k )]->m_tag == 1 )
	  m_buffer_cells.push_back( m_allcells[getCellNumber( i, j, k )] );  
}




// ----------------------------------------------------------------------------
// Sets the linked cell grid in parallel mode
void LinkedCell::set( double cellsize_, int const* nprocsdir,
	int const* MPIcoords, MPINeighbors const* voisins, 
	int const* MPIperiod, string const& oshift )
{
  m_LC_global_origin = m_domain_global_origin;
  m_LC_local_origin = m_domain_local_origin;
  
  // Number of cells and cell edge length in each direction and 
  // Default in 1 unique cell if cellsize_ is zero 
  if ( cellsize_ > 1.e-10 ) 
  {
    m_nbi = (int)( ( App::m_domain_local_size_X + EPSILON ) / cellsize_);
    m_cellsize_X = App::m_domain_local_size_X / m_nbi ;
    m_nbj = (int)( ( App::m_domain_local_size_Y + EPSILON ) / cellsize_);
    m_cellsize_Y = App::m_domain_local_size_Y / m_nbj ;    
    m_nbk = (int)( ( App::m_domain_local_size_Z + EPSILON ) / cellsize_);
    m_cellsize_Z = App::m_domain_local_size_Z / m_nbk ; 
    if ( !m_nbk ) 
    {
      m_nbk = 1;
      m_cellsize_Z = m_cellsize_Y;
    }	  
    
    // Add cells in halo zones
    int suppX = 0, suppY = 0, suppZ = 0;
    if ( MPIcoords[0] != 0 || MPIperiod[0] ) suppX++;
    if ( MPIcoords[0] != nprocsdir[0] - 1 ||  MPIperiod[0] ) suppX++;
    m_nbi += suppX;
    if ( MPIcoords[1] != 0 ||  MPIperiod[1] ) suppY++;
    if ( MPIcoords[1] != nprocsdir[1] - 1 ||  MPIperiod[1] ) suppY++;
    m_nbj += suppY;
    if ( MPIcoords[2] != 0 ||  MPIperiod[2] ) suppZ++;
    if ( MPIcoords[2] != nprocsdir[2] - 1 ||  MPIperiod[2] ) suppZ++;
    m_nbk += suppZ;   
  
    // Local origin of the linked cell grid
    if ( voisins->rank( -1, 0, 0 ) != -1 ) m_LC_local_origin.Move(
    	- m_cellsize_X, 0., 0. );
    if ( voisins->rank( 0, -1, 0 ) != -1 ) m_LC_local_origin.Move(
    	0., - m_cellsize_Y, 0. );
    if ( voisins->rank( 0, 0, -1 ) != -1 ) m_LC_local_origin.Move(
    	0., 0., - m_cellsize_Z );
	
    // Periodicity
    if ( MPIperiod[0] ) m_LC_global_origin.Move( - m_cellsize_X, 0., 0. );
    if ( MPIperiod[1] )	m_LC_global_origin.Move( 0., - m_cellsize_Y, 0. );
    if ( MPIperiod[2] )	m_LC_global_origin.Move( 0., 0., - m_cellsize_Z );
  } 
  else  m_nbi = m_nbj = m_nbk = 1;

  m_nb = m_nbi * m_nbj * m_nbk;

  m_LC_local_xmin = m_LC_local_origin[0];
  m_LC_local_ymin = m_LC_local_origin[1]; 
  m_LC_local_zmin = m_LC_local_origin[2];
  m_LC_local_xmax = m_LC_local_xmin + m_nbi * m_cellsize_X;
  m_LC_local_ymax = m_LC_local_ymin + m_nbj * m_cellsize_Y;
  m_LC_local_zmax = m_LC_local_zmin + m_nbk * m_cellsize_Z;
  m_extendedBBox = new BBox(
  	Point3( m_LC_local_xmin - 0.5 * m_cellsize_X,
		m_LC_local_ymin - 0.5 * m_cellsize_Y,
		m_LC_local_zmin - 0.5 * m_cellsize_Z ),
	Point3( m_LC_local_xmax + 0.5 * m_cellsize_X,
		m_LC_local_ymax + 0.5 * m_cellsize_Y,
		m_LC_local_zmax + 0.5 * m_cellsize_Z ) );

  if ( voisins->rank( 0, 0, 0 ) == 0 )
  {  
    cout << "Linked-cell grid on proc 0" << endl;
    cout << "   Number of cells = " << m_nbi << " " << m_nbj << " " 
    	<< m_nbk << " = " << m_nbi * m_nbj * m_nbk << endl;
    cout << "   Cell size = " << m_cellsize_X << " x " << m_cellsize_Y << 
  	" x " << m_cellsize_Z << endl;
    cout << "   Global origin = " << m_LC_global_origin[X] << " " << 
    	m_LC_global_origin[Y] << " " <<
	m_LC_global_origin[Z] << endl;
    cout << "   Local origin = " << m_LC_local_origin[X] << " " << 
    	m_LC_local_origin[Y] << " " <<
	m_LC_local_origin[Z] << endl;	
    cout << "   Max coordinates of local grid = " << m_LC_local_xmax << " " <<
    	 m_LC_local_ymax << " " << m_LC_local_zmax << endl << endl;
  }

  // Cells construction
  int tag = 0;
  GeoPosition geoLoc = GEOPOS_NONE;
  m_allcells.reserve( m_nb );
  Cell::setNbCellsPerDirection( m_nbi, m_nbj, m_nbk );
  for (int j=0; j<m_nbj; j++) 
  {
    for (int k=0; k<m_nbk; k++) 
    {
      for (int i=0; i<m_nbi; i++) 
      {
        // General case
	geoLoc = GEOPOS_NONE;
	
	if ( i == 0 ) tag = 2;
	
	else if ( i == m_nbi - 1 ) tag = 2;
	
	else if ( i == 1 )
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == 1 || j == 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_SOUTH_WEST;	    
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 2 || j == m_nbj - 3 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_NORTH_WEST;	    
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	       geoLoc = GEOPOS_NORTH_WEST;
	      }
	    }
	  }	  
	  else if ( j == m_nbj - 1 ) tag = 2;	  
	  else 
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_WEST;	    
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_WEST;
	      }
	    }
	  }
	}
	
	else if ( i == m_nbi - 2 )
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == 1 || j == 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_SOUTH_EAST;	    
	    }
	    // 3D geometry
	    else
	    {	    
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 2 || j == m_nbj - 3 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_NORTH_EAST;	    
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST;
	      }
	    }
	  }	  
	  else if ( j == m_nbj - 1 ) tag = 2;	  
	  else 
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_EAST;	    
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_EAST;
	      }
	    }
	  }
	}

	else if ( i == 2 )
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == 1 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_SOUTH_WEST;	    
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST;
	      }
	    }
	  }
	  else if ( j == 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }	  
	  else if ( j == m_nbj - 3 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {	    
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_NORTH_WEST;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_WEST;
	      }
	    }
	  }	  	  
	  else if ( j == m_nbj - 1 ) tag = 2;	  
	  else 
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {	    
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_WEST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_WEST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	}

	else if ( i == m_nbi - 3 )
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == 1 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_SOUTH_EAST;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST;
	      }
	    }
	  }
	  else if ( j == 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }	  
	  else if ( j == m_nbj - 3 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_NORTH_EAST;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_EAST;
	      }
	    }
	  }	  	  
	  else if ( j == m_nbj - 1 ) tag = 2;	  
	  else 
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_EAST_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_EAST_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	}
		
	else
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == m_nbj - 1 ) tag = 2;
	  else if ( j == 1 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_SOUTH;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH;
	      }
	    }  
	  }
	  else if ( j == 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }  
	  }
	  else if ( j == m_nbj - 3 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_BEHIND;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }  
	  }	  	  
	  else if ( j == m_nbj - 2 )
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = GEOPOS_NORTH;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_BEHIND;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH_FRONT;	    
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_NORTH;
	      }
	    }  
	  }
	  else
	  {
	    // 2D geometry
	    if ( m_nbk == 1 )
	    {
	      tag = 0;  
	    }
	    // 3D geometry
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else if ( k == 1 ) 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_BEHIND;
	      }
	      else if ( k == m_nbk - 2 ) 
	      {
	        tag = 1;
	        geoLoc = GEOPOS_FRONT;
	      }
	      else tag = 0;
	    }	  
	  } 	
	}
	m_allcells.push_back( new Cell( getCellNumber( i, j, k ), 
		i, j, k, m_LC_local_origin, 
		m_cellsize_X, m_cellsize_Y, m_cellsize_Z, 		
		m_LC_local_xmax, m_LC_local_ymax, m_LC_local_zmax, 
		tag, geoLoc ) );
      }
    }
  }  

  // Special cases in each direction
  // No neighbor to the left
  if ( voisins->rank( -1, 0, 0 ) == -1 )
  {
    for (int i=0; i<3; i++) 
      for (int j=0; j<m_nbj; j++) 
        for (int k=0; k<m_nbk; k++)
	{ 
	  getCell( i, j, k )->m_tag = getCell( 3, j, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell = 
	  	getCell( 3, j, k )->m_GeoPosCell;
	}
  }
  
  // No neighbor to the right
  if ( voisins->rank( 1, 0, 0 ) == -1 )
  {
    for (int i=m_nbi-3; i<m_nbi; i++) 
      for (int j=0; j<m_nbj; j++) 
        for (int k=0; k<m_nbk; k++)
	{ 
	  getCell( i, j, k )->m_tag = getCell( m_nbi-4, j, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell = 
	  	getCell( m_nbi-4, j, k )->m_GeoPosCell;
	}	  
  }  

  // No neighbor at the bottom
  if ( voisins->rank( 0, -1, 0 ) == -1 )
  {
    for (int j=0; j<3; j++) 
      for (int i=0; i<m_nbi; i++) 
        for (int k=0; k<m_nbk; k++) 
	{
	  getCell( i, j, k )->m_tag = getCell( i, 3, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell = 
	  	getCell( i, 3, k )->m_GeoPosCell;
	}	  
  }
  
  // No neighbor at the top
  if ( voisins->rank( 0, 1, 0 ) == -1 )
  {
    for (int j=m_nbj-3; j<m_nbj; j++)  
      for (int i=0; i<m_nbi; i++)
        for (int k=0; k<m_nbk; k++)
	{ 
	  getCell( i, j, k )->m_tag = getCell( i, m_nbj-4, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell = 
	  	getCell( i, m_nbj-4, k )->m_GeoPosCell;
	}	  
  }  

  // 3D geometry
  if ( m_nbk > 1 )
  {
    // No neighbor behind
    if ( voisins->rank( 0, 0, -1 ) == -1 )
    {
      for (int k=0; k<3; k++) 
        for (int i=0; i<m_nbi; i++) 
          for (int j=0; j<m_nbj; j++)
	  {
	    getCell( i, j, k )->m_tag = getCell( i, j, 3 )->m_tag;
	    getCell( i, j, k )->m_GeoPosCell = 
	    	getCell( i, j, 3 )->m_GeoPosCell;
	  }	  
    }
  
    // No neighbor at the front
    if ( voisins->rank( 0, 0, 1 ) == -1 )
    {
      for (int k=m_nbk-3; k<m_nbk; k++)  
        for (int i=0; i<m_nbi; i++)
          for (int j=0; j<m_nbj; j++) 
	  {
	    getCell( i, j, k )->m_tag = getCell( i, j, m_nbk-4 )->m_tag;
	    getCell( i, j, k )->m_GeoPosCell = 
	    	getCell( i, j, m_nbk-4 )->m_GeoPosCell;
	  }	  
    }  
  }
  
  // Sets the the list of neighboring cells over which broad phase
  // contact detection is performed
  setCellContactNeighborhood(); 
}




// ----------------------------------------------------------------------------
// Sets the list of neighboring cells over which broad phase contact detection 
// is performed (3 cells above, 1 cell to the right, 9 cells behind)
void LinkedCell::setCellContactNeighborhood()
{
  // Contact neighborhood
  // 3 cells above        => i-1<->i+1, j,   k+1
  // 1 cell to the right  => i+1,       j,   k
  // 9 cells behind       => i-1<->i+1, j+1, k+1<->k-1
  Cell *cell_ = NULL;
  Cell *neighborc = NULL;
  int icel, jcel, kcel;
  int cel[3][13] = {{-1,  0, +1, +1, -1,  0, +1, -1,  0, +1, -1,  0, +1},
		    { 0,  0,  0,  0, +1, +1, +1, +1, +1, +1, +1, +1, +1},
		    {+1, +1, +1,  0, +1, +1, +1,  0,  0,  0, -1, -1, -1}};
  
  for (int i=0; i<m_nb; i++) 
  {
    cell_ = m_allcells[i];
    icel = (*cell_)[X];
    jcel = (*cell_)[Y];
    kcel = (*cell_)[Z];

    for (int j=0; j<13; j++) 
    {
      neighborc = getCell( icel+cel[X][j], jcel+cel[Y][j], kcel+cel[Z][j] );
      if ( neighborc ) cell_->addNeighboringCellContact( neighborc );
    }
  } 
}




// ----------------------------------------------------------------------------
// Sets the list of all neighboring cells
void LinkedCell::setCellCompleteNeighborhood()
{
  Cell *cell_ = NULL;
  Cell *neighborc = NULL;
  int icel, jcel, kcel;
  for (int i=0; i<m_nb; i++) 
  {
    cell_ = m_allcells[i];
    icel = (*cell_)[X];
    jcel = (*cell_)[Y];
    kcel = (*cell_)[Z];
    
    for (int k=-1;k<2;++k)
      for (int l=-1;l<2;++l)      
        for (int m=-1;m<2;++m)
          if ( k || l || m )
          {
            neighborc = getCell( icel+k, jcel+l, kcel+m );
            if ( neighborc ) cell_->addNeighboringCell( neighborc );
          } 
  }   
}




// ----------------------------------------------------------------------------
// Computes forces and torques exerted on rigid bodies
void LinkedCell::ComputeForces( double time, double dt,
  	list<Particle*> const* particles )
{
  Particle* reference;
  Cell* cell_;
  list<Particle*> neighborparticles;

  // Particle-particle contacts
  Point3 centre;
  list<Particle*>::const_iterator particle;
  list<Particle*>::iterator neighborp;
  list<Cell*>::iterator around;
  int id[3];

  for( particle=particles->begin(); particle!=particles->end();
       particle++)
  {
    reference = *particle;

    // Search for neighboring particle
    // In the local cell: we only detect collisions with neighboring particles
    // that are located beyond the current particle in the list of particles in
    // the local cell. This enables to detect only once a contact between 2
    // particles  
    // In neighboring cells: the contact neighborood (3 cells above, 1 cell to 
    // the right, 9 cells behind) guarantees that contacts between particles
    // are only detected once
    centre = *(*particle)->getPosition();
    Cell::GetCell( centre, id );
    cell_ = getCell( id[X], id[Y], id[Z] );
    neighborp = find( cell_->m_particles.begin(), cell_->m_particles.end(),
		   reference );

    neighborp++;
    for( ; neighborp!=cell_->m_particles.end(); neighborp++ )
      neighborparticles.push_back(*neighborp);

    for( around=cell_->m_neighborsContact.begin();
    	 around!=cell_->m_neighborsContact.end(); around++ ) 
      for( neighborp=(*around)->m_particles.begin();
           neighborp!=(*around)->m_particles.end(); neighborp++ )
        neighborparticles.push_back(*neighborp);

    // Narrow phase contact detection of a reference particle with neighboring 
    // particles and force computation
    // In case of a composite particle, we swap neighbor and reference such that
    // the calling object is always the composite particle
    for( neighborp=neighborparticles.begin(); 
    	neighborp!=neighborparticles.end(); neighborp++ )
      if( (*neighborp)->isCompositeParticle() )
        (*neighborp)->InterAction( reference, dt, time, this );
      else
        reference->InterAction( *neighborp, dt, time, this );

    neighborparticles.clear();

    // Narrow phase contact detection of a reference particle with neighboring 
    // obstacles and force computation
    // In case of a composite particle, we swap neighbor and reference such that
    // the calling object is always the composite particle
    list<SimpleObstacle*>::iterator myObs;
    for( myObs=cell_->m_obstacles.begin();
         myObs!=cell_->m_obstacles.end(); myObs++ )
    {
      if ( reference->isCompositeParticle() )
        reference->InterAction( *myObs, dt, time, this );
      else  
        (*myObs)->InterAction( reference, dt, time, this );
    }
  }
}




// ----------------------------------------------------------------------------
// Returns a pointer to the cell given its ijk indexing
Cell* LinkedCell::getCell( int i, int j, int k ) const
{
  Cell* cell_ = NULL;
  bool valid = -1 < i && i < m_nbi && -1 < j && j < m_nbj 
  	&& -1 < k && k < m_nbk;  
  if ( valid ) cell_ =  m_allcells[getCellNumber( i, j, k )];

  return ( cell_ );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the cell that contains a point
Cell* LinkedCell::getCell( const Point3 &position ) const
{
  int id[3];
  Cell::GetCell( position, id );
  
  return ( getCell( id[X], id[Y], id[Z] ) );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the cell given its ijk indexing
int LinkedCell::getCellNumber( int i, int j, int k ) const
{
  return ( j * m_nbk * m_nbi + k * m_nbi + i ) ;
}  




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// using the method Component::isContact
bool LinkedCell::isContact( Particle const* particle ) const
{
  bool contact = false;

  Point3 centre = *(particle->getPosition());
  int id[3];
  Cell::GetCell( centre, id );
  Cell* cell_ = getCell( id[X], id[Y], id[Z] );
  Cell* neighborc = NULL ;  
  
  if ( cell_ )
  {
    // In the local cell
    contact = cell_->isContact( particle );
      
    // In the neighboring cells  
    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)      
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    neighborc = getCell( id[X]+k, id[Y]+l, id[Z]+m );
            if ( neighborc ) contact = neighborc->isContact( particle );
	  } 
  }
  
  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// using the method Component::isContactWithCrust
bool LinkedCell::isContactWithCrust( Particle const* particle ) const
{
  bool contact = false;

  Point3 centre = *(particle->getPosition());
  int id[3];
  Cell::GetCell( centre, id );
  Cell* cell_ = getCell( id[X], id[Y], id[Z] );
  Cell* neighborc = NULL ;  

  if ( cell_ )
  {
    // In the local cell
    contact = cell_->isContactWithCrust( particle );

    // In the neighboring cells  
    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)      
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    neighborc = getCell( id[X]+k, id[Y]+l, id[Z]+m );
            if ( neighborc ) contact = neighborc->isContactWithCrust( 
	    	particle );
	  }  
  }
  
  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is close to another component
// using the method Component::isClose
bool LinkedCell::isClose( Particle const* particle ) const
{
  bool contact = false;

  Point3 centre = *(particle->getPosition());
  int id[3];
  Cell::GetCell( centre, id );
  Cell* cell_ = getCell( id[X], id[Y], id[Z] );
  Cell* neighborc = NULL ; 
  
  if ( cell_ )
  {
    // In the local cell
    contact = cell_->isClose( particle );

    // In the neighboring cells  
    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)      
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    neighborc = getCell( id[X]+k, id[Y]+l, id[Z]+m );
            if ( neighborc ) contact = neighborc->isClose( particle );
	  }
      
  }

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is close to another component
// using the method Component::isCloseWithCrust
bool LinkedCell::isCloseWithCrust( Particle const* particle ) const
{
  bool contact = false;

  Point3 centre = *(particle->getPosition());
  int id[3];
  Cell::GetCell( centre, id );
  Cell* cell_ = getCell( id[X], id[Y], id[Z] );
  Cell* neighborc = NULL ; 
  
  if ( cell_ )
  {
    // In the local cell
    contact = cell_->isCloseWithCrust( particle );

    // In the neighboring cells  
    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)      
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    neighborc = getCell( id[X]+k, id[Y]+l, id[Z]+m );
            if ( neighborc ) contact = neighborc->isCloseWithCrust( particle );
	  }      
  }
  
  return ( contact );
}




// ----------------------------------------------------------------------------
// Links a particle with the linked cell grid without checking if the particle 
// overlaps with another rigid body
void LinkedCell::Link( Particle* particle )
{
  Point3 centre = *(particle->getPosition());
  int id[3];
  Cell::GetCell( centre, id );
  Cell* cell_ = getCell( id[X], id[Y], id[Z] );
  particle->setCellTagGeoPosition( cell_, cell_->m_tag, cell_->m_GeoPosCell );
  particle->setCellTagGeoPosition_nm1( cell_, cell_->m_tag, 
  	cell_->m_GeoPosCell );   
  	
  if ( !cell_->contains( particle ) ) cell_->add( particle );
}




// ----------------------------------------------------------------------------
// Links an obstacle with the linked cell grid
void LinkedCell::Link( Obstacle* obstacle )
{  
  // We search intersection between the obstacle and twice expanded cells
  // i.e. cells expanded by a least the maximum circumscribed radius of the
  // largest particle in the simulation, hence guaranteeing that no collision
  // between particles and the obstacle is missed 

  AppCollision::Link(obstacle);
  list<SimpleObstacle*> list_obstacles = obstacle->getObstacles();
  list<SimpleObstacle*>::iterator myObs;
  Cell* cell_ = NULL;
  Transform CelPosition;
  double alpha=2.;  
  Point3 const* cg = NULL;
  
  for (myObs=list_obstacles.begin();myObs!=list_obstacles.end();myObs++)
  {
    RigidBody const* obstacleRigidBody = (*myObs)->getRigidBody();
    BBox obsBox = obstacleRigidBody->BoxRigidBody();
    if ( (*myObs)->getMaterial() == "periode" ) alpha = 3.1;
    else alpha = 2.;
    Convex* convexCell = new Box( alpha * m_cellsize_X, alpha * m_cellsize_Y,
    	alpha * m_cellsize_Z);
    RigidBody CelRigidBody( convexCell, CelPosition );
    
    // Intersection géométrique de la cellule avec l'obstacle
    // Intersection of the cell with the obstacle
    for (int i=0; i<m_nb; i++) 
    {
      cell_ = m_allcells[i];
      cg = cell_->getCentre();      
      if ( obsBox.InZone( cg, 0.5 * alpha * m_cellsize_X, 
      		0.5 * alpha * m_cellsize_Y,
      		0.5 * alpha * m_cellsize_Z ) )
      {
        CelRigidBody.setOrigin( (*cg)[X], (*cg)[Y], (*cg)[Z] );
	if ( CelRigidBody.isContact( *obstacleRigidBody ) ) 
        {
          cell_->addObstacle( *myObs );
          (*myObs)->add( cell_ );
        }
      }
    }

    // Rem: we do not explicitly destroy the convex convexCell because the 
    // destructor of CelRigidBody, object of type RigidBody, takes care of it
    // (cf RigidBody.cpp)
  }   
}




// ----------------------------------------------------------------------------
// Updates links between particles & obstacles and the linked cell grid
void LinkedCell::LinkUpdate( double time, double dt,
  	list<Particle*>* particles ) throw(SimulationError)
{
  // If the particle is not active anymore, we remove it from the cell it
  // belongs to and from the list of active particles
  list<Particle*>::iterator particle;
  for (particle=particles->begin(); particle!=particles->end(); ) 
  {
    switch ( (*particle)->getActivity() ) 
    {
      case COMPUTE:
        particle++;
        break;

      default:
        (*particle)->getCell()->remove( *particle );
        particle = particles->erase( particle );
        break;
    }
  }

  // Update obstacles in case they move
  list<SimpleObstacle*>::iterator myObs;
  for (myObs=m_allObstacles.begin();myObs!=m_allObstacles.end();myObs++)
    if ( (*myObs)->hasMoved() )
    {
      // Check whether the obstacle intersects the local linked cell grid
      if ( intersect( *(*myObs)->getObstacleBox() , *m_extendedBBox ) )
      {
	// If the obstacle has not been linked yet, we perform a classic Link
	// Otherwise we perform a LinkUpdate
	if ( (*myObs)->getInCells()->empty() ) Link( *myObs );
	else LinkUpdate( time, dt, *myObs );
      }	
    } 

  // Update active particles
  for (particle=particles->begin(); particle!=particles->end(); 
	 particle++) 
    LinkUpdateActiveParticle( *particle );
}




// ----------------------------------------------------------------------------
// Updates the link of an active particle and the linked cell grid
void LinkedCell::LinkUpdateActiveParticle( Particle* particle ) 
	throw (SimulationError)
{
  if ( particle->getActivity() != COMPUTE )
  {
      cout << "\nParticle not active " << particle->getID() << endl;
      cout << "            " << *particle->getPosition() << endl;
      GrainsExec::m_exception_Simulation = true;
      throw(SimulationError("LinkedCell::LinkUpdateActiveParticle"));  
  }

  // Copies the cell the particle belonged to, the particle tag and 
  // the geographic location of the particle from current time to previous 
  // time
  particle->copyCellTagGeoPosition_n_to_nm1();	

  // Cell the particle belongs to at the current discrete time
  Point3 centre = *(particle->getPosition());
  int id[3];  
  Cell::GetCell( centre, id );
  Cell* cellNew = getCell( id[X], id[Y], id[Z] );

  if ( cellNew == NULL ) 
  {
    cout << "\nParticle " << particle->getID()       << endl;
    cout << "            " << *particle->getPosition() << endl;
    GrainsExec::m_exception_Simulation = true;    
    throw(SimulationError("LinkedCell::LinkUpdateActiveParticle"));
  }

  // Update if the particle has moved to a different cell 
  Cell* cellNm1 = particle->getCellNm1();   
  if ( cellNew != cellNm1 ) 
  {
    cellNew->add( particle );
    cellNm1->remove( particle ); 
  }
  
  // Set the cell the particle belongs to, the particle tag and 
  // the geographic location of the particle at the current time
  particle->setCellTagGeoPosition( cellNew, cellNew->m_tag, 
  	cellNew->m_GeoPosCell );      
}




// ----------------------------------------------------------------------------
// Updates the link between the cells and a simple obstacle
void LinkedCell::LinkUpdate( double time, double dt, SimpleObstacle *myObs )
{
  if ( myObs->performLinkUpdate() )
  {
    const RigidBody* obstacleRigidBody = myObs->getRigidBody();
    BBox obsBox = obstacleRigidBody->BoxRigidBody(); 
    Cell* cell_ = NULL;
    Vector3 deplMax;
    Point3 const* cg = NULL;    
    
    // le coefficient 1.2 donne une marge d'erreur de 20%, ce qui signifie
    // qu'on suppose que le vecteur vitesse de l'obstacle sur les n pas de time
    // suivants ne varie pas de plus de 20%
    // Attention: rien dans le code verifie cette hypothèse !!    
    int updateFreq = myObs->getObstacleLinkedCellUpdateFrequency();
    double coefApprox = updateFreq == 1 ? 1. : 1.2;
    deplMax = myObs->vitesseMaxPerDirection() * coefApprox
   	* updateFreq * dt;
        
    Vector3 CelExtent( m_cellsize_X + deplMax[X], m_cellsize_Y + deplMax[Y],
    	m_cellsize_Z + deplMax[Z] );
	
    myObs->resetInCells();  
    for (int i=0; i<m_nb; i++) 
    {
      cell_ = m_allcells[i];
      cg = cell_->getCentre();      
      if ( obsBox.InZone( cg, CelExtent[X], CelExtent[Y], CelExtent[Z] ) )
      {
        cell_->addObstacle( myObs );
        myObs->add( cell_ );
      }
    }    
  } 
  
//   const list<Cell*>* voisinageCourant = myObs->getInCells();
//   list<Cell*> voisinageEtendu = *voisinageCourant;
//   const list<Cell*>*	celluleVoisinageComplet = NULL;
//   list<Cell*>::const_iterator icellule,icelluleVoisine;
//   list<Cell*>::iterator il;  
//   const RigidBody* obstacleRigidBody = myObs->getRigidBody();
//   BBox obsBox = obstacleRigidBody->BoxRigidBody(); 
//   Cell* cellule = NULL;
//   Transform CelPosition;
//   double alpha = 2.;  
//   Convex* convexCell = new Box(alpha*m_cellsize_X,alpha*m_cellsize_Y,
//     	alpha*m_cellsize_Z);
//   RigidBody CelRigidBody(convexCell,CelPosition);
// 
//   // Voisinage etendu  
//   for (icellule=voisinageCourant->begin();icellule!=voisinageCourant->end();
//   	icellule++)
//   {
//     celluleVoisinageComplet = (*icellule)->getCompleteNeighborhood();
//     for (icelluleVoisine=celluleVoisinageComplet->begin();
//     	icelluleVoisine!=celluleVoisinageComplet->end();icelluleVoisine++)
//       voisinageEtendu.push_back(*icelluleVoisine);
//   } 
//   voisinageEtendu.sort();
//   voisinageEtendu.unique();   
// 
//   // Intersection géométrique de la cellule avec l'obstacle 
//   // dans le voisinage etendu
//   myObs->resetInCells();
//   for (il=voisinageEtendu.begin();il!=voisinageEtendu.end();il++)
//   {
//     cellule = *il;
//     Point3 cg = cellule->Gravite();      
//     if (obsBox.InZone(cg,m_cellsize_X,m_cellsize_Y,m_cellsize_Z))
//     {
//       CelRigidBody.setOrigin(cg[X],cg[Y],cg[Z]);
//       if (CelRigidBody.isContact(*obstacleRigidBody)) 
//       {
//         cellule->addObstacle(myObs);
//         myObs->add(cellule);
//       }
//     }
//   }     
}




// ----------------------------------------------------------------------------
// Removes a particle from the linked cell grid
void LinkedCell::remove( Particle* particle )
{
  Point3 centre = *(particle->getPosition());
  int id[3];
  Cell::GetCell( centre, id );
  Cell* cell_ = getCell( id[X], id[Y], id[Z] );
  if ( cell_ ) cell_->remove( particle );
}




// ----------------------------------------------------------------------------
// Removes an obstacle from the linked cell grid
void LinkedCell::remove( SimpleObstacle* obs )
{
  AppCollision::remove( obs );
  obs->resetInCells();  
}


 
  
// ----------------------------------------------------------------------------
// Output operator
ostream& operator <<( ostream& f, LinkedCell const& LC )
{
  vector<Cell*>::const_iterator iv; 
  int nbCelSansObstacle = 0;
  for (iv=LC.m_allcells.begin();iv!=LC.m_allcells.end();iv++)
    if ( (*iv)->numberOfObstacles() == 0 ) nbCelSansObstacle++; 

  f << "Total number of cells = " << LC.m_nb << endl;
  f << "Number of cells without any obstacle in their neighborhood = " << 
  	nbCelSansObstacle << endl;  
  f << "Number of cells in X x Y x Z = " << LC.m_nbi << " x " << LC.m_nbj 
  	<< " x " << LC.m_nbk << endl;   
  f << "Cell size in X x Y x Z = " << LC.m_cellsize_X << " x " 
  	<< LC.m_cellsize_Y << " x " << LC.m_cellsize_Z << endl;     
  f << "Local origin = " << LC.m_LC_local_origin;
  f << "Max coordinates of local grid = " << LC.m_LC_local_xmax << " " 
  	<< LC.m_LC_local_ymax << " " << LC.m_LC_local_zmax << endl;
  f << "CELLS" << endl;
  for (iv=LC.m_allcells.begin();iv!=LC.m_allcells.end();iv++)
    f << *(*iv) << endl << endl;
  
  return ( f );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the linked cell grid
bool LinkedCell::isInLinkedCell( Point3 const& position ) const
{
  bool isIn = true;
  
  if ( position[0] < m_LC_local_xmin || position[0] > m_LC_local_xmax
  	|| position[1] < m_LC_local_ymin || position[1] > m_LC_local_ymax
  	|| position[2] < m_LC_local_zmin || position[2] > m_LC_local_zmax ) 
    isIn = false;
  
  return ( isIn );
}  




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the linked cell grid
bool LinkedCell::isInLinkedCell( double const& gx, double const& gy,
	double const& gz ) const
{
  bool isIn = true;
  
  if ( gx < m_LC_local_xmin || gx > m_LC_local_xmax
  	|| gy < m_LC_local_ymin || gy > m_LC_local_ymax
  	|| gz < m_LC_local_zmin || gz > m_LC_local_zmax ) 
    isIn = false;
  
  return ( isIn );
}  




// ----------------------------------------------------------------------------
// Updates active particles with an interior tag (tag=0)
void LinkedCell::updateInteriorTag( double time,
	list<Particle*>* particles,  
	list<Particle*>* particlesHalozone,
	GrainsMPIWrapper const* wrapper )
{
  list<Particle*>::iterator particle;
  int tag_nm1,tag;
  Cell* current_cell = NULL;
  
  for (particle=particles->begin(); particle!=particles->end(); 
	 particle++) 
  {
    tag_nm1 = (*particle)->getTag();
    current_cell = getCell( *(*particle)->getPosition() );
    if ( tag_nm1 == 0 )
    {
      tag = (*particle)->setTag( current_cell->m_tag );
    
      // Interior to Halozone (0 -> 1)
      if ( tag == 1 )
      {
        if ( GrainsExec::m_MPI_verbose )
	{
	  ostringstream oss;
	  oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT ) <<
		" Interior to Halozone (0 -> 1)               Id = " << 
      		(*particle)->getID() << " " << *(*particle)->getPosition();
	  GrainsMPIWrapper::addToMPIString( oss.str() );
	}
        (*particle)->setGeoPosition( current_cell->m_GeoPosCell );
	particlesHalozone->push_back( *particle );  
      }
    }
  }
} 




// ----------------------------------------------------------------------------
// Updates active particles with a halo tag (tag=1) or a clone tag (tag=2)
void LinkedCell::updateHalozoneCloneTag( double time,
	list<Particle*>* particlesHalozone,
	list<Particle*>* particlesClones,
	GrainsMPIWrapper const* wrapper )
{
  list<Particle*>::iterator particle;
  int tag;
  
  for (particle=particlesHalozone->begin(); 
  	particle!=particlesHalozone->end(); ) 
  {
    tag = (*particle)->setTag( 
    	getCell( *(*particle)->getPosition() )->m_tag );
    
    // Halozone to Clone (1 -> 2)
    if ( tag == 2 )
    {
      if ( GrainsExec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT ) <<
      		" Halozone to Clone (1 -> 2)                  Id = " << 
      		(*particle)->getID() << " " << *(*particle)->getPosition();
        GrainsMPIWrapper::addToMPIString( oss.str() );
      }	
      particlesClones->push_back( *particle );
      (*particle)->setGeoPosition( GEOPOS_NONE );       
      particle = particlesHalozone->erase( particle );
    }
    // Halozone to Interior (1 -> 0)
    else if (tag==0)
    {
      if ( GrainsExec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT ) <<
      		" Halozone to Interior (1 -> 0)               Id = " << 
      		(*particle)->getID() << " " << *(*particle)->getPosition();
        GrainsMPIWrapper::addToMPIString( oss.str() );
      }
      (*particle)->setGeoPosition( GEOPOS_NONE );        
      particle = particlesHalozone->erase( particle );      
    } 
    else particle++;           
  }  

  for (particle=particlesClones->begin(); 
  	particle!=particlesClones->end(); ) 
  {
    tag = (*particle)->setTag( 
    	getCell( *(*particle)->getPosition() )->m_tag );
    
    // Clone to Halozone (2 -> 1)
    if ( tag == 1 )
    {
      if ( GrainsExec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT ) <<
      		" Clone to Halozone (2 -> 1)                  Id = " << 
      		(*particle)->getID() << " " << *(*particle)->getPosition();
        GrainsMPIWrapper::addToMPIString( oss.str() );
      }
      (*particle)->setGeoPosition( 
      	getCell( *(*particle)->getPosition() )->m_GeoPosCell );
      particlesHalozone->push_back( *particle );
      particle = particlesClones->erase( particle );  
    }
    else particle++;           
  }    
} 




// ----------------------------------------------------------------------------
// Removes clone particles that exited the local linked cell grid
void LinkedCell::DestroyOutOfDomainClones( double time,
	list<Particle*>* particles,
	list<Particle*>* particlesClones,
	GrainsMPIWrapper const* wrapper )
{
  list<Particle*>::iterator particle;
  list<App*>::iterator app;
  Particle *pdestroy=NULL;
  
  for (particle=particlesClones->begin(); particle!=particlesClones->end();)
  {
    if ( !isInLinkedCell( *(*particle)->getPosition() ) )
    {
      if ( GrainsExec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT ) <<
      		" Destroy clone                               Id = " << 
      		(*particle)->getID() << " " << *(*particle)->getPosition();
        GrainsMPIWrapper::addToMPIString( oss.str() );
      }
      pdestroy = *particle;
	
      // Removes the clone particle from the last cell it belonged to
      pdestroy->getCellNm1()->remove( pdestroy );
 
      // Removes the clone particle from the list of active particles
      removeParticleFromList( *particles, pdestroy );
	  
      // Destroy the clone particle
      delete pdestroy;
      
      // Removes the clone particle from the list of active clone particles
      particle = particlesClones->erase( particle );
    }
    else particle++;
  }   

}








// ----------------------------------------------------------------------------
// Returns the cell edge length in a direction
double LinkedCell::getCellSize( int const& dir ) const
{
  double size = 0.;
  
  switch(dir)
  {
    case 0: size = m_cellsize_X;
      break;
    case 1: size = m_cellsize_Y;
      break;
    case 2: size = m_cellsize_Z;
      break;
  }
  
  return ( size );
}      
  



// ----------------------------------------------------------------------------
// Returns a pointer to the vector of all cells of the linked cell grid
vector<Cell*> const* LinkedCell::getAllCells() const
{
  return ( &m_allcells );
}




// ----------------------------------------------------------------------------
// Attempts to insert a particle in serial mode */
bool LinkedCell::insertParticleSerial( Particle* particle, 
	list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	vector<Particle*> const* ReferenceParticles,
	bool const& periodic, bool const& force_insertion )
{
  bool insert = false, contact = true;
  GeoPosition geoloc = GEOPOS_NONE; 
  
  // If insertion is not forced, check contacts with other particles and 
  // obstacles
  if ( !force_insertion )
  {
    // Check with master particle
    contact = isContactWithCrust( particle );

    // In case of a periodic domain, check if periodic clones have contacts
    if ( periodic && !contact )
    {
      // Get the cell where the master particle is located
      Point3 centre = *(particle->getPosition());
      int id[3];
      Cell::GetCell( centre, id );
      Cell* cell_ = getCell( id[X], id[Y], id[Z] );
      geoloc = cell_->m_GeoPosCell;

      // Loop over the domain periodic vectors for this geographic position
      for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size() && 
      	!contact;++i)
      {
        // Translate particle by periodic vector i
	particle->Translate( m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );

	// Check contact of the periodic clone
	contact = isContactWithCrust( particle );
	
	// If no contact, translate back to original position
        if ( !contact )
	  particle->Translate( - m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );
      }      
    }  
  }
  
  // If no contact or force insertion, insert particle, create periodic clones 
  // and insert them
  if ( !contact || force_insertion )
  {
    insert = true;      
      
    // Link master particle
    Link( particle );
      
    // Create periodic clones and insert them
    if ( periodic )
    {
      Particle* clone = NULL;
      
      // Loop over the domain periodic vectors for this geographic position
      for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size();++i)
      {
        // Create periodic clone
	clone = particle->createCloneCopy( particle->getID(), 
		(*ReferenceParticles)[particle->getGeometricType()], 
		*(particle->getTranslationalVelocity()),
		*(particle->getQuaternionRotation()),	 
		*(particle->getAngularVelocity()),	 
		*(particle->getRigidBody()->getTransform()),
		COMPUTE );		
	
	// Translate to its periodic position
	clone->Translate( m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );

        // Link periodic clone
        Link( clone );
		
	// Insert into periodic clone map
	particlesPeriodicClones->insert( 
		pair<int,Particle*>( particle->getID(), clone ) );
	
	// Insert into active particle list
	particles->push_back( clone );
      }      
    }
  }
  
  return ( insert ); 
}




// ----------------------------------------------------------------------------
// Updates periodic clones and destroy those out of the linked cell grid in 
// serial mode
void LinkedCell::updateDestroyPeriodicClones( list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones )
{
  list<Cell*>::iterator il;
  list<Particle*>::iterator particle;
  GeoPosition geoloc = GEOPOS_NONE; 
  multimap<int,Particle*>::iterator imm;
  size_t ncid = 0;
  pair < multimap<int,Particle*>::iterator, 
  	multimap<int,Particle*>::iterator > crange;
  int particleID = 0; 
  bool found = false;
  Particle* periodic_clone = NULL ; 
  	 
  // Loop over all buffer cells to update periodic clones
  for (il=m_buffer_cells.begin();il!=m_buffer_cells.end();il++)
  {
    // Get the cell geographic position
    geoloc = (*il)->m_GeoPosCell;
    
    // For all particles in the cell
    for (particle=(*il)->m_particles.begin();particle!=(*il)->m_particles.end();
    	particle++)
    {
      particleID = (*particle)->getID();
      ncid = particlesPeriodicClones->count( particleID );
      
      // Case 1: 1 periodic clone
      if ( ncid == 1 )
      {
        // Find the periodic clone
	imm = particlesPeriodicClones->find( particleID );

	// Update the periodic clone features
	imm->second->setTransform( 
		*((*particle)->getRigidBody()->getTransform()) );
	imm->second->Translate( m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][0]] );
	imm->second->setTranslationalVelocity(
		*((*particle)->getTranslationalVelocity()) );
	imm->second->setQuaternionRotation( 
		*((*particle)->getQuaternionRotation()) );
	imm->second->setAngularVelocity( *((*particle)->getAngularVelocity()) );
      }
      else
      // Case 2: multiple periodic clones
      {      
        crange = particlesPeriodicClones->equal_range( particleID );
                 	
	// Loop over the domain periodic vectors for this geographic position
        for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size();++i)
        {
          found = false;
	  
	  // Find the periodic clone assuming that the master particle has not 
	  // moved by more than twice its crust thickness (an assumption that
	  // is always verified in Grains3D)	  
	  for (imm=crange.first; imm!=crange.second && !found; )
            if ( imm->second->getPosition()->DistanceTo( 
	    	*((*particle)->getPosition()) 
		+ m_domain_global_periodic_vectors[
			m_periodic_vector_indices[geoloc][i]] ) < 
		2. * (*particle)->getCrustThickness() )
	    {
	      periodic_clone = imm->second;
	      found = true;
	    }
	    else imm++;
	    
	  // Update the periodic clone features
	  periodic_clone->setTransform( 
		*((*particle)->getRigidBody()->getTransform()) );
	  periodic_clone->Translate( m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );
	  periodic_clone->setTranslationalVelocity(
		*((*particle)->getTranslationalVelocity()) );
	  periodic_clone->setQuaternionRotation( 
		*((*particle)->getQuaternionRotation()) );
	  periodic_clone->setAngularVelocity( 
	  	*((*particle)->getAngularVelocity()) );	   
        }      
      }     
    }
  }
  
  // Destroy periodic clones that are out of the linked cell grid
  Particle *pdestroy = NULL;
  for (imm=particlesPeriodicClones->begin();
  	imm!=particlesPeriodicClones->end();)
    if ( isInLinkedCell( *(imm->second->getPosition()) ) )
      imm++;
    else
    {
      pdestroy = imm->second;

      // Suppress the periodic clone from the cell it belonged to before exiting
      // the linked cell grid
      // Note: at that point, we have not done LinkUpdate yet, so the the cell 
      // it belonged to before exiting the linked cell grid is getCell(), not 
      // getCellNm1()
      pdestroy->getCell()->remove( pdestroy );
      
      // Remove the periodic particle from the list of active particles
      removeParticleFromList( *particles, pdestroy );
	  
      // Destroy the clone particle
      delete pdestroy;
      
      // Removes the periodic clone particle from the map of periodic clones
      imm = particlesPeriodicClones->erase( imm );          
    }
}




// ----------------------------------------------------------------------------
// Creates/destroys periodic clones after LinkUpdate in serial mode
void LinkedCell::createDestroyPeriodicClones( list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	vector<Particle*> const* ReferenceParticles )
{
  int tag = 0, tag_nm1 = 0;
  GeoPosition geoloc = GEOPOS_NONE, geoloc_nm1 = GEOPOS_NONE;
  multimap<int,Particle*>::iterator imm;
  pair < multimap<int,Particle*>::iterator, 
  	multimap<int,Particle*>::iterator > crange;
  int particleID = 0;
  bool found = false; 
  Particle* periodic_clone = NULL;
      
  for (list<Particle*>::iterator particle=particles->begin();
	particle!=particles->end();particle++)
  {
    tag = (*particle)->getTag();
    tag_nm1 = (*particle)->getTagNm1();
    
    switch( tag_nm1 )
    {
      case 0:
        // Particle moved from interior to buffer: create new periodic clones
        if ( tag == 1 )
        {  
          geoloc = (*particle)->getGeoPosition();
      
          // Loop over the domain periodic vectors for this geographic position
          for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size();++i)
          {
	    // Create periodic clone
	    periodic_clone = (*particle)->createCloneCopy( (*particle)->getID(),
		(*ReferenceParticles)[(*particle)->getGeometricType()], 
		*((*particle)->getTranslationalVelocity()),
		*((*particle)->getQuaternionRotation()),	 
		*((*particle)->getAngularVelocity()),	 
		*((*particle)->getRigidBody()->getTransform()),
		COMPUTE );
	
	    // Translate to its periodic position
	    periodic_clone->Translate( m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );

            // Link periodic clone
            Link( periodic_clone );
		
	    // Insert into periodic clone map
	    particlesPeriodicClones->insert( 
		pair<int,Particle*>( (*particle)->getID(), periodic_clone ) );
	
	    // Insert into active particle list
	    particles->push_back( periodic_clone );
          }	  
        }
	break;
	
      case 1:
        // If tag == 0: buffer to interior, nothing to do as periodic clones 
        // have exited the linked cell grid and were destroyed 
        // by updateDestroyPeriodicClones
      
        // Particle moved from buffer to halozone: add particle to the periodic
        // clone multimap
        if ( tag == 2 )
        {
	  particlesPeriodicClones->insert( 
		pair<int,Particle*>( (*particle)->getID(), *particle ) );
        }
	else if ( tag == 1 )
	{
	  geoloc = (*particle)->getGeoPosition();
	  geoloc_nm1 = (*particle)->getGeoPositionNm1();
	  
	  // If change of geographic position, search in periodic clone multimap
	  // that all periodic clones exist and if a periodic clone does not
	  // exist, create it
	  if ( geoloc != geoloc_nm1 )
	  {
            particleID = (*particle)->getID();
	    crange = particlesPeriodicClones->equal_range( particleID );
                 	
	    // Loop over the domain periodic vectors for this geographic 
	    // position
            for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size();++i)
            {
              found = false;
	  
	      // Find the periodic clone (same assumption as in 
	      // updateDestroyPeriodicClones)	  
	      for (imm=crange.first; imm!=crange.second && !found; )
                if ( imm->second->getPosition()->DistanceTo( 
	    		*((*particle)->getPosition()) 
			+ m_domain_global_periodic_vectors[
				m_periodic_vector_indices[geoloc][i]] ) < 
			2. * (*particle)->getCrustThickness() )
	          found = true;
	        else imm++;
		
	      // If the periodic clone is not found, create it
	      if ( !found )
	      {
	        // Create periodic clone
	        periodic_clone = (*particle)->createCloneCopy( 
			(*particle)->getID(), 
			(*ReferenceParticles)[(*particle)->getGeometricType()], 
			*((*particle)->getTranslationalVelocity()),
			*((*particle)->getQuaternionRotation()),	 
			*((*particle)->getAngularVelocity()),	 
			*((*particle)->getRigidBody()->getTransform()),
			COMPUTE );
	
	        // Translate to its periodic position
	        periodic_clone->Translate( m_domain_global_periodic_vectors[
			m_periodic_vector_indices[geoloc][i]] );

                // Link periodic clone
                Link( periodic_clone );
		
	        // Insert into periodic clone map
	        particlesPeriodicClones->insert( 
			pair<int,Particle*>( (*particle)->getID(), 
				periodic_clone ) );
	
	        // Insert into active particle list
	        particles->push_back( periodic_clone );	      
	      }
	    }		  
	  }
	}
	break;
	
      default: // i.e. 2
        // Particle moved from halozone to buffer: remove particle from the 
        // periodic clone multimap
        if ( tag == 1 )
        {
          particleID = (*particle)->getID(); 	  
	  crange = particlesPeriodicClones->equal_range( particleID );
	  found = false;
          for (imm=crange.first; imm!=crange.second && !found; )
            if ( imm->second == *particle )
	    {
	      imm = particlesPeriodicClones->erase( imm );
	      found = true;
	    }
	    else imm++;
	}	      
        break;
    }                        
  }
}
