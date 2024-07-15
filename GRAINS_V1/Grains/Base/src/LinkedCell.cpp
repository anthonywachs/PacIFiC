#include "GrainsMPIWrapper.hh"
#include "AllComponents.hh"
#include "GrainsBuilderFactory.hh"
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
size_t LinkedCell::set( double cellsize_, string const& oshift )
{
  size_t error = 0;
  
  // In 2D, we set the domain size in Z to 2*EPSILON and the origin in Z 
  // to -EPSILON such that the Z position of the cells is always 0
  if ( GrainsBuilderFactory::getContext() == DIM_2 )
  {
    m_domain_global_origin[Z] = - EPSILON;
    m_domain_local_origin[Z] = - EPSILON; 
    m_domain_global_size_Z = 2. * EPSILON;
    m_domain_local_size_Z = 2. * EPSILON;
    m_domain_global_periodicity[Z] = false;
  }
  
  m_LC_global_origin = m_domain_global_origin;
  m_LC_global_max = m_domain_global_origin;
  m_LC_global_max[X] += m_domain_global_size_X;
  m_LC_global_max[Y] += m_domain_global_size_Y;  
  m_LC_global_max[Z] += m_domain_global_size_Z;  
  m_LC_local_origin = m_domain_local_origin;

  // Number of cells and cell edge length in each direction and
  // Default is 1 unique cell if cellsize_ is zero
  if ( cellsize_ > EPSILON )
  {
    m_nbi = (int)( ( App::m_domain_local_size_X + EPSILON ) / cellsize_ );
    if ( !m_nbi ) m_nbi = 1;
    m_cellsize_X = App::m_domain_local_size_X / m_nbi ;
    m_nbj = (int)( ( App::m_domain_local_size_Y + EPSILON ) / cellsize_ );
    if ( !m_nbj ) m_nbj = 1;    
    m_cellsize_Y = App::m_domain_local_size_Y / m_nbj ;
    if ( GrainsBuilderFactory::getContext() == DIM_2 )
    {
      m_nbk = 1;
      m_cellsize_Z = App::m_domain_local_size_Z;
    }
    else
    {
      m_nbk = (int)( ( App::m_domain_local_size_Z + EPSILON ) / cellsize_);
      if ( !m_nbk ) m_nbk = 1;
      m_cellsize_Z = App::m_domain_local_size_Z / m_nbk ;
    }
  }
  else  
  {
    m_nbi = m_nbj = m_nbk = 1;
    m_cellsize_X = App::m_domain_local_size_X;
    m_cellsize_Y = App::m_domain_local_size_Y;    
    m_cellsize_Z = App::m_domain_local_size_Z;
  }
  
  // Periodicity
  if ( m_domain_global_periodicity[X] )
  {
    m_LC_global_origin.Move( - m_cellsize_X, 0., 0. );
    m_LC_global_max.Move( m_cellsize_X, 0., 0. );
    m_LC_local_origin.Move( - m_cellsize_X, 0., 0. );
    m_nbi += 2;
  }
  if ( m_domain_global_periodicity[Y] )
  {
    m_LC_global_origin.Move( 0., - m_cellsize_Y, 0. );
    m_LC_global_max.Move( 0., m_cellsize_Y, 0. );    
    m_LC_local_origin.Move( 0., - m_cellsize_Y, 0. );
    m_nbj += 2;
  }
  if ( m_domain_global_periodicity[Z] )
  {
    m_LC_global_origin.Move( 0., 0., - m_cellsize_Z );
    m_LC_global_max.Move( 0., 0., m_cellsize_Z );     
    m_LC_local_origin.Move( 0., 0., - m_cellsize_Z );
    m_nbk += 2;
  }      

  m_nb = m_nbi * m_nbj * m_nbk;

  cout << oshift << "Number of cells = " << m_nbi << " " << m_nbj << " "
  	<< m_nbk << " = " << m_nbi * m_nbj * m_nbk << endl;
  cout << oshift << "Cell size  = " << m_cellsize_X << " x " << m_cellsize_Y <<
  	" x " << m_cellsize_Z << endl;
  cout << oshift << "Global origin = " << m_LC_global_origin << endl;
  cout << oshift << "Global max = " << m_LC_global_max << endl;  
  cout << oshift << "Local origin = " << m_LC_global_origin << endl;

  m_LC_local_origin[X] = m_LC_local_origin[0];
  m_LC_local_origin[Y] = m_LC_local_origin[1];
  m_LC_local_origin[Z] = m_LC_local_origin[2];
  m_LC_local_max[X] = m_LC_local_origin[X] + m_nbi * m_cellsize_X;
  m_LC_local_max[Y] = m_LC_local_origin[Y] + m_nbj * m_cellsize_Y;
  m_LC_local_max[Z] = m_LC_local_origin[Z] + m_nbk * m_cellsize_Z;
  m_extendedBBox = new BBox(
  	Point3( m_LC_local_origin[X] - 0.5 * m_cellsize_X,
		m_LC_local_origin[Y] - 0.5 * m_cellsize_Y,
		m_LC_local_origin[Z] - 0.5 * m_cellsize_Z ),
	Point3( m_LC_local_max[X] + 0.5 * m_cellsize_X,
		m_LC_local_max[Y] + 0.5 * m_cellsize_Y,
		m_LC_local_max[Z] + 0.5 * m_cellsize_Z ) );

  cout << oshift << "Local size  = " << m_LC_local_max[X] - m_LC_local_origin[X]
  	<< " x " << m_LC_local_max[Y] - m_LC_local_origin[Y]
	<< " x " << m_LC_local_max[Z] - m_LC_local_origin[Z] << endl;

  // Cells construction
  m_allcells.reserve( m_nb );
  Cell::setNbCellsPerDirection( m_nbi, m_nbj, m_nbk );
  for (int j=0; j<m_nbj; j++)
    for (int k=0; k<m_nbk; k++)
      for (int i=0; i<m_nbi; i++)
        m_allcells.push_back( new Cell( getCellNumber( i, j, k ),
                i, j, k, m_LC_local_origin,
		m_cellsize_X, m_cellsize_Y, m_cellsize_Z,
                m_LC_local_max[X], m_LC_local_max[Y], m_LC_local_max[Z] ) );

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

	if ( m_nbi == 3 )
	// Special case of a single cell in the main domain
	{
	  // Set tag = 1 and geopos = GEOPOS_EASTWEST in 2nd row
	  m_allcells[getCellNumber( 1, j, k )]->m_tag = 1;
	  m_allcells[getCellNumber( 1, j, k )]->m_GeoPosCell = GEOPOS_EASTWEST;	
	}
	else
	// General case
	{
	  // Set tag = 1 and geopos = GEOPOS_WEST in 2nd row
	  m_allcells[getCellNumber( 1, j, k )]->m_tag = 1;
	  m_allcells[getCellNumber( 1, j, k )]->m_GeoPosCell = GEOPOS_WEST;

	  // Set tag = 1 and geopos = GEOPOS_EAST in penultimate row
	  m_allcells[getCellNumber( m_nbi - 2, j, k )]->m_tag = 1;
	  m_allcells[getCellNumber( m_nbi - 2, j, k )]->m_GeoPosCell =
		GEOPOS_EAST;
	}
      }
  }

  if ( m_domain_global_periodicity[Y] )
  {
    if ( !m_domain_global_periodicity[X] || 
    	 ( m_domain_global_periodicity[X] && m_nbi > 3 && m_nbj > 3 ) )
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
	
	if ( m_nbj == 3 )
	// Special case of a single cell in the main domain
	{
	  // Set tag = 1 and geopos = GEOPOS_NORTHSOUTH in 2nd row
	  m_allcells[getCellNumber( i, 1, k )]->m_tag = 1;
	  m_allcells[getCellNumber( i, 1, k )]->m_GeoPosCell = 
	  	GEOPOS_NORTHSOUTH;	
	}
	else
	// General case
	{	
	  // Set tag = 1 and geopos += GEOPOS_SOUTH in 2nd row
	  if ( m_allcells[getCellNumber( i, 1, k )]->m_tag != 2 )
	  {
	    m_allcells[getCellNumber( i, 1, k )]->m_tag = 1;
	    switch( int(m_allcells[getCellNumber( i, 1, k )]->m_GeoPosCell) )
	    {
	      case GEOPOS_NONE:
	        m_allcells[getCellNumber( i, 1, k )]->m_GeoPosCell = 
			GEOPOS_SOUTH;
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
    }
    else if ( m_domain_global_periodicity[X] && m_nbi == 3 && m_nbj == 3 )
    {
      for (int i=0; i<m_nbi; i++)
        for (int j=0; j<m_nbj; j++)
          for (int k=0; k<m_nbk; k++)
	  {
	    m_allcells[getCellNumber( i, j, k )]->m_tag = 2;
	    m_allcells[getCellNumber( i, j, k )]->m_GeoPosCell = GEOPOS_NONE ;
	  }
      
      for (int k=0; k<m_nbk; k++)
      {
        m_allcells[getCellNumber( 1, 1, k )]->m_tag = 1;
	m_allcells[getCellNumber( 1, 1, k )]->m_GeoPosCell = 
		GEOPOS_EASTWESTNORTHSOUTH ;
      }     
    }
    else 
    {
      cout << "Serial bi-periodicity in XY with a single cell in the"
      	<< " main domain in one direction and more than one cell in the other"
	<< " direction is not handled" << endl;
      cout << "This configuration has not been implemented yet" << endl;
      error = 1;
    }
  }

  if ( !error )
  {
    if ( m_domain_global_periodicity[Z] )
    {
      if ( m_nbk > 3 )
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
      else
      {
        cout << "Serial periodicity in Z with a single cell in the"
      	<< " main domain in the Z direction is not handled" << endl;
        cout << "This configuration has not been implemented yet" << endl;
        error = 1;
      }
    }
  }      

  // list of buffer cells (i.e. tag = 1)
  if ( !error )
    for (int j=0; j<m_nbj; j++)
      for (int k=0; k<m_nbk; k++)
        for (int i=0; i<m_nbi; i++)
          if ( m_allcells[getCellNumber( i, j, k )]->m_tag == 1 )
	    m_buffer_cells.push_back( m_allcells[getCellNumber( i, j, k )] );
	    
  return ( error );	    
}




// ----------------------------------------------------------------------------
// Sets the linked cell grid in parallel mode
size_t LinkedCell::set( double cellsize_, int const* nprocsdir,
	int const* MPIcoords, MPINeighbors const* voisins,
	string const& oshift )
{
  size_t error = 0;

  // In 2D, we set the domain size in Z to 2*EPSILON and the origin in Z 
  // to -EPSILON such that the Z position of the cells is always 0
  if ( GrainsBuilderFactory::getContext() == DIM_2 )
  {
    m_domain_global_origin[Z] = - EPSILON;
    m_domain_local_origin[Z] = - EPSILON; 
    m_domain_global_size_Z = 2. * EPSILON;
    m_domain_local_size_Z = 2. * EPSILON;
    m_domain_global_periodicity[Z] = false;
  }

  m_LC_global_origin = m_domain_global_origin;
  m_LC_global_max = m_domain_global_origin;
  m_LC_global_max[X] += m_domain_global_size_X;
  m_LC_global_max[Y] += m_domain_global_size_Y;  
  m_LC_global_max[Z] += m_domain_global_size_Z;  
  m_LC_local_origin = m_domain_local_origin;

  // Number of cells and cell edge length in each direction and
  // Default is 1 unique cell if cellsize_ is zero
  if ( cellsize_ > EPSILON )
  {
    m_nbi = (int)( ( App::m_domain_local_size_X + EPSILON ) / cellsize_);
    if ( !m_nbi ) m_nbi = 1;
    m_cellsize_X = App::m_domain_local_size_X / m_nbi ;
    m_nbj = (int)( ( App::m_domain_local_size_Y + EPSILON ) / cellsize_);
    if ( !m_nbj ) m_nbj = 1;
    m_cellsize_Y = App::m_domain_local_size_Y / m_nbj ;
    if ( GrainsBuilderFactory::getContext() == DIM_2 )
    {
      m_nbk = 1;
      m_cellsize_Z = App::m_domain_local_size_Z;
    } 
    else
    {
      m_nbk = (int)( ( App::m_domain_local_size_Z + EPSILON ) / cellsize_);
      if ( !m_nbk ) m_nbk = 1;
      m_cellsize_Z = App::m_domain_local_size_Z / m_nbk ;    
    }   
  }
  else  
  {
    m_nbi = m_nbj = m_nbk = 1;
    m_cellsize_X = App::m_domain_local_size_X;
    m_cellsize_Y = App::m_domain_local_size_Y;    
    m_cellsize_Z = App::m_domain_local_size_Z;
  }

  // Add cells in clone zones
  int suppX = 0, suppY = 0, suppZ = 0;
  if ( MPIcoords[0] != 0 || m_domain_global_periodicity[0] ) suppX++;
  if ( MPIcoords[0] != nprocsdir[0] - 1 ||  m_domain_global_periodicity[0] ) 
    suppX++;
  m_nbi += suppX;
  if ( MPIcoords[1] != 0 ||  m_domain_global_periodicity[1] ) suppY++;
  if ( MPIcoords[1] != nprocsdir[1] - 1 ||  m_domain_global_periodicity[1] ) 
    suppY++;
  m_nbj += suppY;
  if ( MPIcoords[2] != 0 ||  m_domain_global_periodicity[2] ) suppZ++;
  if ( MPIcoords[2] != nprocsdir[2] - 1 ||  m_domain_global_periodicity[2] ) 
    suppZ++;
  m_nbk += suppZ;

  // Local origin of the linked cell grid
  if ( voisins->rank( -1, 0, 0 ) != -1 ) m_LC_local_origin.Move(
    	- m_cellsize_X, 0., 0. );
  if ( voisins->rank( 0, -1, 0 ) != -1 ) m_LC_local_origin.Move(
    	0., - m_cellsize_Y, 0. );
  if ( voisins->rank( 0, 0, -1 ) != -1 ) m_LC_local_origin.Move(
    	0., 0., - m_cellsize_Z );

  // Periodicity
  if ( m_domain_global_periodicity[0] )
  { 
    m_LC_global_origin.Move( - m_cellsize_X, 0., 0. );
    m_LC_global_max.Move( m_cellsize_X, 0., 0. );    
  }
  if ( m_domain_global_periodicity[1] )
  {	
    m_LC_global_origin.Move( 0., - m_cellsize_Y, 0. );
    m_LC_global_max.Move( 0., m_cellsize_Y, 0. );    
  }
  if ( m_domain_global_periodicity[2] )	
  {
    m_LC_global_origin.Move( 0., 0., - m_cellsize_Z );
    m_LC_global_max.Move( 0., 0., m_cellsize_Z );    
  }

  m_nb = m_nbi * m_nbj * m_nbk;

  if ( voisins->rank( 0, 0, 0 ) == 0 )
  {
    cout << oshift << "Linked-cell grid on proc 0" << endl;
    cout << oshift << "Number of cells = " << m_nbi << " " << m_nbj << " "
    	<< m_nbk << " = " << m_nbi * m_nbj * m_nbk << endl;
    cout << oshift << "Cell size = " << m_cellsize_X << " x " << m_cellsize_Y <<
  	" x " << m_cellsize_Z << endl;
    cout << oshift << "Global origin = " << m_LC_global_origin << endl;
    cout << oshift << "Global max = " << m_LC_global_max << endl;    
    cout << oshift << "Local origin = " << m_LC_global_origin << endl;
  }

  m_LC_local_origin[X] = m_LC_local_origin[0];
  m_LC_local_origin[Y] = m_LC_local_origin[1];
  m_LC_local_origin[Z] = m_LC_local_origin[2];
  m_LC_local_max[X] = m_LC_local_origin[X] + m_nbi * m_cellsize_X;
  m_LC_local_max[Y] = m_LC_local_origin[Y] + m_nbj * m_cellsize_Y;
  m_LC_local_max[Z] = m_LC_local_origin[Z] + m_nbk * m_cellsize_Z;
  m_extendedBBox = new BBox(
  	Point3( m_LC_local_origin[X] - 0.5 * m_cellsize_X,
		m_LC_local_origin[Y] - 0.5 * m_cellsize_Y,
		m_LC_local_origin[Z] - 0.5 * m_cellsize_Z ),
	Point3( m_LC_local_max[X] + 0.5 * m_cellsize_X,
		m_LC_local_max[Y] + 0.5 * m_cellsize_Y,
		m_LC_local_max[Z] + 0.5 * m_cellsize_Z ) );

  if ( voisins->rank( 0, 0, 0 ) == 0 )
  {
    cout << oshift << "Local size  = " 
    	<< m_LC_local_max[X] - m_LC_local_origin[X]
  	<< " x " << m_LC_local_max[Y] - m_LC_local_origin[Y]
	<< " x " << m_LC_local_max[Z] - m_LC_local_origin[Z] << endl;
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
	      else
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_WEST;
	      }
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
	      else
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH_EAST;
	      }
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
	      else if ( k == 1  )
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
	      else
	      {
	        tag = 1;
	        geoLoc = GEOPOS_EAST;
	      }
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
	      else
	      {
	        tag = 1;
	        geoLoc = GEOPOS_SOUTH;
	      }
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
		m_LC_local_max[X], m_LC_local_max[Y], m_LC_local_max[Z],
		tag, geoLoc ) );
      }
    }
  }

  // Special cases in each direction
  // No neighbor to the left
  if ( voisins->rank( -1, 0, 0 ) == -1 )
  {
    for (int i=0; i<2; i++)
      for (int j=0; j<m_nbj; j++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCell( i, j, k )->m_tag = getCell( 2, j, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell =
	  	getCell( 2, j, k )->m_GeoPosCell;
	}
  }

  // No neighbor to the right
  if ( voisins->rank( 1, 0, 0 ) == -1 )
  {
    for (int i=m_nbi-2; i<m_nbi; i++)
      for (int j=0; j<m_nbj; j++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCell( i, j, k )->m_tag = getCell( m_nbi-3, j, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell =
	  	getCell( m_nbi-3, j, k )->m_GeoPosCell;
	}
  }

  // No neighbor at the bottom
  if ( voisins->rank( 0, -1, 0 ) == -1 )
  {
    for (int j=0; j<2; j++)
      for (int i=0; i<m_nbi; i++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCell( i, j, k )->m_tag = getCell( i, 2, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell =
	  	getCell( i, 2, k )->m_GeoPosCell;
	}
  }

  // No neighbor at the top
  if ( voisins->rank( 0, 1, 0 ) == -1 )
  {
    for (int j=m_nbj-2; j<m_nbj; j++)
      for (int i=0; i<m_nbi; i++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCell( i, j, k )->m_tag = getCell( i, m_nbj-3, k )->m_tag;
	  getCell( i, j, k )->m_GeoPosCell =
	  	getCell( i, m_nbj-3, k )->m_GeoPosCell;
	}
  }

  // 3D geometry
  if ( m_nbk > 1 )
  {
    // No neighbor behind
    if ( voisins->rank( 0, 0, -1 ) == -1 )
    {
      for (int k=0; k<2; k++)
        for (int i=0; i<m_nbi; i++)
          for (int j=0; j<m_nbj; j++)
	  {
	    getCell( i, j, k )->m_tag = getCell( i, j, 2 )->m_tag;
	    getCell( i, j, k )->m_GeoPosCell =
	    	getCell( i, j, 2 )->m_GeoPosCell;
	  }
    }

    // No neighbor at the front
    if ( voisins->rank( 0, 0, 1 ) == -1 )
    {
      for (int k=m_nbk-2; k<m_nbk; k++)
        for (int i=0; i<m_nbi; i++)
          for (int j=0; j<m_nbj; j++)
	  {
	    getCell( i, j, k )->m_tag = getCell( i, j, m_nbk-3 )->m_tag;
	    getCell( i, j, k )->m_GeoPosCell =
	    	getCell( i, j, m_nbk-3 )->m_GeoPosCell;
	  }
    }
  }

  // Sets the the list of neighboring cells over which broad phase
  // contact detection is performed
  setCellContactNeighborhood();
  
  return ( error );  
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
  list<SimpleObstacle*>::iterator myObs;

  // Particle-particle contacts
  Point3 centre;
  list<Particle*>::const_iterator particle;
  list<Particle*>::iterator neighborp;
  list<Cell*>::iterator around;
  int id[3];

  for (particle=particles->cbegin(); particle!=particles->cend(); particle++)
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
    {
      if( (*neighborp)->isCompositeParticle() )
        (*neighborp)->InterAction( reference, dt, time, this );
      else
        reference->InterAction( *neighborp, dt, time, this );
    }

    neighborparticles.clear();

    // Narrow phase contact detection of a reference particle with neighboring
    // obstacles and force computation
    // In case of a composite particle, we swap neighbor and reference such that
    // the calling object is always the composite particle
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
// Returns a list of pointers to the cell that contains a point
// and the neighboring cells to that cell
list<Cell*> LinkedCell::getCellAndCellNeighborhood( Point3 const& position ) 
	const
{
  list<Cell*> cells;
  Cell* neighborc = NULL;  
  
  int id[3];
  Cell::GetCell( position, id );
  cells.push_back( getCell( id[X], id[Y], id[Z] ) );
  
  for (int k=-1;k<2;++k)
    for (int l=-1;l<2;++l)
      for (int m=-1;m<2;++m)
        if ( k || l || m )
        {
          neighborc = getCell( id[X]+k, id[Y]+l, id[Z]+m );
          if ( neighborc ) cells.push_back( neighborc );
        }  
    
  return ( cells ) ;
}



// ----------------------------------------------------------------------------
// Returns the cell number given its ijk indexing
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

  AppCollision::Link( obstacle );
  
  list<SimpleObstacle*>::iterator myObs;
  Cell* cell_ = NULL;
  Transform cellPosition;
  double alpha = 2.;
  Point3 const* cg = NULL;
  bool add = false;

  for (myObs=m_allSimpleObstacles.begin();myObs!=m_allSimpleObstacles.end();
  	myObs++)
  {
    RigidBody const* obstacleRigidBody = (*myObs)->getRigidBody();
    BBox const* obstacleBBox = (*myObs)->getObstacleBox();
    Vector3 cellBoxExtension( 0.5 * alpha * m_cellsize_X, 
    	0.5 * alpha * m_cellsize_Y,
	0.5 * alpha * m_cellsize_Z );    
    Convex* cellBox = new Box( 2. * cellBoxExtension[X], 
    	2. * cellBoxExtension[Y],
    	2. * cellBoxExtension[Z] );
    RigidBody cellBoxRigidBody( cellBox, cellPosition );

    // Intersection of the cell with the obstacle
    for (int i=0; i<m_nb; i++)
    {
      cell_ = m_allcells[i];
      cg = cell_->getCentre();
      add = false;
      if ( obstacleBBox->InZone( cg, cellBoxExtension[X], cellBoxExtension[Y],
      		cellBoxExtension[Z] ) )
      {
        if ( (*myObs)->isSTLObstacle() ) add = true; // Temporary, TO DO
	else 
	{
	  cellBoxRigidBody.setOrigin( (*cg)[X], (*cg)[Y], (*cg)[Z] );
	  add = cellBoxRigidBody.isContact( *obstacleRigidBody );	
	}
	
	if ( add )
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
  	list<Particle*>* particles )
{
  try
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
    
    // Obstacle periodicity
    if ( m_domain_global_periodic ) m_obstacles->periodicity( this );

    // Update obstacles in case they move
    list<SimpleObstacle*>::iterator myObs;
    for (myObs=m_allSimpleObstacles.begin();myObs!=m_allSimpleObstacles.end();myObs++)
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
  catch (const SimulationError&) {
    throw SimulationError();
  }
}




// ----------------------------------------------------------------------------
// Updates the link of an active particle and the linked cell grid
void LinkedCell::LinkUpdateActiveParticle( Particle* particle )
{
  try
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
  catch (const SimulationError&) {
    throw SimulationError();
  }
}




// ----------------------------------------------------------------------------
// Updates the link between the cells and a simple obstacle
void LinkedCell::LinkUpdate( double time, double dt, SimpleObstacle *myObs )
{
  // We search intersection between the obstacle AABBox and twice expanded cells
  // i.e. cells expanded by a least the maximum circumscribed radius of the
  // largest particle in the simulation, hence guaranteeing that no collision
  // between particles and the obstacle is missed
  // This method is highly sub-optimal for the following two reasons:
  // 1) the whole linked cell grid is searched
  // 2) we use the AABBox of the obstacle such that the geometric intersection
  // test is faster than relying on GJK, leading to unnecessary cells whenever
  // the obstacle is not "reasonably" aligned with the coordinate axis

  if ( myObs->performLinkUpdate() )
  {
    BBox const* obstacleBBox = myObs->getObstacleBox();
    Cell* cell_ = NULL;
    Point3 const* cg = NULL;
    double alpha = 2.;    
    Vector3 cellBoxExtension( 0.5 * alpha * m_cellsize_X, 
    	0.5 * alpha * m_cellsize_Y,
	0.5 * alpha * m_cellsize_Z );

    myObs->resetInCells();
    for (int i=0; i<m_nb; i++)
    {
      cell_ = m_allcells[i];
      cg = cell_->getCentre();
      if ( obstacleBBox->InZone( cg, cellBoxExtension[X], cellBoxExtension[Y],
      		cellBoxExtension[Z] ) )
      {
        cell_->addObstacle( myObs );
        myObs->add( cell_ );
      }
    }
  }
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
  f << "Max coordinates of local grid = " << LC.m_LC_local_max[X] << " "
  	<< LC.m_LC_local_max[Y] << " " << LC.m_LC_local_max[Z] << endl;
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

  if ( position[0] < m_LC_local_origin[X] || position[0] > m_LC_local_max[X]
    || position[1] < m_LC_local_origin[Y] || position[1] > m_LC_local_max[Y]
    || position[2] < m_LC_local_origin[Z] || position[2] > m_LC_local_max[Z] )
    isIn = false;

  return ( isIn );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the linked cell grid
bool LinkedCell::isInLinkedCell( double const& gx, double const& gy,
	double const& gz ) const
{
  bool isIn = true;

  if ( gx < m_LC_local_origin[X] || gx > m_LC_local_max[X]
  	|| gy < m_LC_local_origin[Y] || gy > m_LC_local_max[Y]
  	|| gz < m_LC_local_origin[Z] || gz > m_LC_local_max[Z] )
    isIn = false;

  return ( isIn );
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
      		" Destroy clone               Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
        GrainsMPIWrapper::addToMPIString( oss.str() );
      }
      pdestroy = *particle;

      // Removes the clone particle from the last cell it belonged to
      // Note: we user getCell because we have moved the particle already
      // but did not execute LinkUpdate yet (see the sequence of steps in 
      // GrainsMPI::Simulation) 
      pdestroy->getCell()->remove( pdestroy );

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
  bool insert = false, contact = false;
  GeoPosition geoloc = GEOPOS_NONE;

  // If insertion is not forced, check contacts with other particles and
  // obstacles

  // Check with master particle
  if ( !force_insertion ) contact = isContactWithCrust( particle );

  // In case of a periodic domain, we need the geoloc of the cell the particle
  // belongs to
  if ( periodic && !contact )
  {
    // Get the cell where the master particle is located
    Point3 centre = *(particle->getPosition());
    int id[3];
    Cell::GetCell( centre, id );
    Cell* cell_ = getCell( id[X], id[Y], id[Z] );
    geoloc = cell_->m_GeoPosCell;
  }

  // In case of a periodic domain, check if periodic clones have contacts
  if ( periodic && !force_insertion && !contact )
  {  
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
		COMPUTE, particle->getContactMap() );

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
// Attempts to insert a particle in parallel mode
pair<bool,bool> LinkedCell::insertParticleParallel( double time,
	Particle* particle, 
	list<Particle*>* particles,
	list<Particle*>* particlesClones,
	vector<Particle*> const* ReferenceParticles,
	bool const& periodic,
    	bool const& force_insertion,
	GrainsMPIWrapper const* wrapper )
{
  GeoPosition geoloc = GEOPOS_NONE;
  int source = 0;
  bool contact = false;

  // First is "Is in LinkedCell?" and second is "Is contact?"
  pair<bool,bool> insert(false,false);

  // Check whether the particle belongs to the local Linked Cell
  insert.first = isInLinkedCell( *(particle->getPosition()) );
  
  // Check contacts
  if ( !force_insertion )
  {
    if ( insert.first ) insert.second = isContactWithCrust( particle );  
    insert.second = wrapper->max_INT( insert.second );     
  }
  
  // Periodic clones
  if ( !insert.second && periodic ) 
  {
    // The particle is only in the local domain of a single proc
    if ( insert.first && isInLocalDomain( particle->getPosition() ) )
    {
      // Get the cell where the master particle is located
      Point3 centre = *(particle->getPosition());
      int id[3];
      Cell::GetCell( centre, id );
      Cell* cell_ = getCell( id[X], id[Y], id[Z] );
      geoloc = cell_->m_GeoPosCell;
      source = wrapper->get_rank();    
    }
    
    // Broadcast the geoloc
    source = wrapper->max_INT( source ); 
    geoloc = GeoPosition( wrapper->Broadcast_INT( int(geoloc), source ) );
    
    // Loop over the domain periodic vectors for this geographic position
    for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size() &&
	!contact;++i)
    {
      // Translate particle by periodic vector i
      particle->Translate( m_domain_global_periodic_vectors[
	m_periodic_vector_indices[geoloc][i]] );

      // Check contact of the periodic clone
      if ( isInLinkedCell( *(particle->getPosition()) ) )
        contact = isContactWithCrust( particle );

      // If no contact, translate back to original position
      if ( !contact )
	particle->Translate( - m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );
    }    
    contact = wrapper->max_INT( contact ); 
    if ( contact ) insert.second = true;
    
    // If no contact for periodic clones, create them in local domains
    // where they exist
    if ( !contact )
    {
      Particle* clone = NULL;

      // Loop over the domain periodic vectors for this geographic position
      for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size();++i)
      {
        // Translate particle by periodic vector i
        particle->Translate( m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );
		
        if ( isInLinkedCell( *(particle->getPosition()) ) )
	{
          if ( GrainsExec::m_MPI_verbose )
          {
            ostringstream oss;
            oss << "   t=" << GrainsExec::doubleToString(time,TIMEFORMAT)
		<< " Create Clone                Id = " 
		<< particle->getID()
		<< " Type = " << particle->getGeometricType() << " " 
		<< (*(particle->getPosition()))[X] << " " 
		<< (*(particle->getPosition()))[Y] << " " 
		<< (*(particle->getPosition()))[Z]
		<< endl;
            GrainsMPIWrapper::addToMPIString(oss.str());
          }

          // Create periodic clone
	  clone = particle->createCloneCopy( particle->getID(),
		(*ReferenceParticles)[particle->getGeometricType()],
		*(particle->getTranslationalVelocity()),
		*(particle->getQuaternionRotation()),
		*(particle->getAngularVelocity()),
		*(particle->getRigidBody()->getTransform()),
		COMPUTE, particle->getContactMap() );

          // Link periodic clone
          Link( clone );

	  // Insert into active particle list and clone particle list
	  particles->push_back( clone );
	  particlesClones->push_back( clone );
	}
	
        // Translate back to original position
        particle->Translate( - m_domain_global_periodic_vectors[
		m_periodic_vector_indices[geoloc][i]] );	
      }
    }    
  } 
      
  // If contact is false and inLinkCell is true, link particle 
  if ( !insert.second && insert.first ) Link( particle ); 
  
  return ( insert );  
}




// ----------------------------------------------------------------------------
// Updates periodic clones and destroy those out of the linked cell grid in
// serial mode
void LinkedCell::updateDestroyPeriodicClones( list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	bool destroyOnly )
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
  if ( !destroyOnly )
  {
    for (il=m_buffer_cells.begin();il!=m_buffer_cells.end();il++)
    {
      // Get the cell geographic position
      geoloc = (*il)->m_GeoPosCell;

      // For all particles in the cell
      for (particle=(*il)->m_particles.begin();
      	particle!=(*il)->m_particles.end(); particle++)
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
	  imm->second->setAngularVelocity( 
	  	*((*particle)->getAngularVelocity()) );
	  imm->second->setContactMap( *((*particle)->getContactMap()) ); 
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
	    periodic_clone->setContactMap( *((*particle)->getContactMap()) ); 	
          }
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
		COMPUTE, (*particle)->getContactMap() );

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

        // Particle moved from buffer zone to clone zone: add particle to the 
	// periodic clone multimap
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
			COMPUTE, (*particle)->getContactMap() );

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
        // Particle moved from clone zone to buffer zone: remove particle from 
	// the periodic clone multimap
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

  // Special of periodicity with a single cell in the main domain in that 
  // direction
  // Some particles that moved from tag 2 to tag 1 may not have all their
  // periodic clones in the system depending on the order of particles
  // in the list of active particles
  // Consequently, we check whether all periodic clones of each particle 
  // tagged 1 exists and if not we create them
  // Important: this problem does not arise as soon as there are 2 cells in 
  // the main domain in that direction
  if ( ( m_nbi == 3 && m_domain_global_periodicity[X] )
  	|| ( m_nbj == 3 && m_domain_global_periodicity[Y] )
	|| ( m_nbk == 3 && m_domain_global_periodicity[Z] ) )
  {
    for (list<Particle*>::iterator particle=particles->begin();
	particle!=particles->end();particle++)
    {
      tag = (*particle)->getTag();
      tag_nm1 = (*particle)->getTagNm1();

      if ( tag == 1 && tag_nm1 == 2 )
      {	
        geoloc = (*particle)->getGeoPosition();
        crange = particlesPeriodicClones->equal_range( particleID );
      
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
		COMPUTE, (*particle)->getContactMap() );

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
      }
    }
  }
}




// ----------------------------------------------------------------------------
// Checks periodic clones in serial mode when a simulation is reloaded
void LinkedCell::checkPeriodicClonesReload( list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	vector<Particle*> const* ReferenceParticles, 
	double const& time )
{
  GeoPosition geoloc = GEOPOS_NONE;
  list<Cell*>::iterator il;
  list<Particle*>::iterator particle;
  multimap<int,Particle*>::iterator imm;
  pair < multimap<int,Particle*>::iterator,
  	multimap<int,Particle*>::iterator > crange;
  int particleID = 0;
  bool found = false;
  Particle* periodic_clone = NULL;
  size_t counter = 0;
  Point3 gc;  

  // Destroy periodic clones that are out of the linked cell grid
  updateDestroyPeriodicClones( particles, particlesPeriodicClones, true );

  // Create periodic clones that do not yet exist
  for (il=m_buffer_cells.begin();il!=m_buffer_cells.end();il++)
  {
    // Get the cell geographic position
    geoloc = (*il)->m_GeoPosCell;

    // For all particles in the cell
    for (particle=(*il)->m_particles.begin();particle!=(*il)->m_particles.end();
    	particle++)
    {
      particleID = (*particle)->getID();
      crange = particlesPeriodicClones->equal_range( particleID );

      // Loop over the domain periodic vectors for this geographic position
      for ( size_t i=0;i<m_periodic_vector_indices[geoloc].size();++i)
      {
        found = false;

	// Try to find the periodic clone
	for (imm=crange.first; imm!=crange.second && !found; )
          if ( imm->second->getPosition()->DistanceTo(
	    	*((*particle)->getPosition())
		+ m_domain_global_periodic_vectors[
			m_periodic_vector_indices[geoloc][i]] ) <
		2. * (*particle)->getCrustThickness() ) found = true;
	  else imm++; 
	  
	// Create the periodic clone if it was not found
	if ( !found )
	{
          ++counter;
	  
	  // Create periodic clone
	  periodic_clone = (*particle)->createCloneCopy(
		(*particle)->getID(),
		(*ReferenceParticles)[(*particle)->getGeometricType()],
		*((*particle)->getTranslationalVelocity()),
		*((*particle)->getQuaternionRotation()),
		*((*particle)->getAngularVelocity()),
		*((*particle)->getRigidBody()->getTransform()),
		COMPUTE, (*particle)->getContactMap() );

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
    }
  }
  
  if ( counter )
    cout << "Creation of additional periodic clones (different"
	<< " link cell grid or clones not stored in the reload file)" 
	<< endl; 
}




// ----------------------------------------------------------------------------
// Returns an array of point coordinates of the local grid in a direction
vector<double> LinkedCell::local_coordinates( size_t const& dir ) const
{
  int npts = ( dir == 0 ? m_nbi + 1 : ( dir == 1 ? m_nbj + 1 : m_nbk + 1 ) );
  vector<double> coord( npts, 0. );
  double step = ( dir == 0 ? m_cellsize_X : 
  	( dir == 1 ? m_cellsize_Y : m_cellsize_Z ) );
  
  for (int i=0;i<npts;i++) 
    coord[i] = m_LC_local_origin[dir] + double(i) * step;
  
  return ( coord );
}




// ----------------------------------------------------------------------------
// Returns an array of point coordinates of the global grid in a direction
vector<double> LinkedCell::global_coordinates( size_t const& dir ) const
{
  int npts = ( dir == 0 ? 
	int( ( m_LC_global_max[X] - m_LC_global_origin[X] + EPSILON ) 
	/ m_cellsize_X ) + 1 : ( dir == 1 ? 
	int( ( m_LC_global_max[Y] - m_LC_global_origin[Y] + EPSILON ) 
	/ m_cellsize_Y ) + 1 : 
	int( ( m_LC_global_max[Z] - m_LC_global_origin[Z] + EPSILON ) 
	/ m_cellsize_Z ) + 1 ) );
  vector<double> coord( npts, 0. );
  double step = ( dir == 0 ? m_cellsize_X : 
  	( dir == 1 ? m_cellsize_Y : m_cellsize_Z ) );
  
  for (int i=0;i<npts;i++) 
    coord[i] = m_LC_global_origin[dir] + double(i) * step;
  
  return ( coord );
}




// ----------------------------------------------------------------------------
// Checks that none of the structured array positions is exactly 
// at a limit of the linked cell grid, otherwise shift by 1e-12 
void LinkedCell::checkStructuredArrayPositionsMPI( struct StructArrayInsertion* 
    	InsertionArray, GrainsMPIWrapper const* wrapper ) const
{
  Point3 position;
  double geoshift = 1.e-12;      
  double deltax = ( InsertionArray->box.ptB[X]
  	- InsertionArray->box.ptA[X] ) / double(InsertionArray->NX) ;
  double deltay = ( InsertionArray->box.ptB[Y]
  	- InsertionArray->box.ptA[Y] ) / double(InsertionArray->NY) ;
  double deltaz = ( InsertionArray->box.ptB[Z]
  	- InsertionArray->box.ptA[Z] ) / double(InsertionArray->NZ) ;
  vector<size_t> coorMatchLocLim( 3, 0 );
  size_t k, l, m;
  
  for (k=0;k<InsertionArray->NX && !coorMatchLocLim[X];++k)
  {
    position[X] = InsertionArray->box.ptA[X] + ( double(k) + 0.5 ) * deltax;
    if ( fabs( position[X] - m_LC_local_origin[X] ) < geoshift
    	|| fabs( position[X] - m_LC_local_max[X] ) < geoshift ) 
      coorMatchLocLim[X] = 1;
  }
  
  for (l=0;l<InsertionArray->NY && !coorMatchLocLim[Y];++l)
  {
    position[Y] = InsertionArray->box.ptA[Y] + ( double(l) + 0.5 ) * deltay;
    if ( fabs( position[Y] - m_LC_local_origin[Y] ) < geoshift
    	|| fabs( position[Y] - m_LC_local_max[Y] ) < geoshift ) 
      coorMatchLocLim[Y] = 1;
  }  

  for (m=0;m<InsertionArray->NZ && !coorMatchLocLim[Z];++m)
  {
    position[Z] = InsertionArray->box.ptA[Z] + ( double(m) + 0.5 ) * deltaz;
    if ( fabs( position[Z] - m_LC_local_origin[Z] ) < geoshift
    	|| fabs( position[Z] - m_LC_local_max[Z] ) < geoshift ) 
      coorMatchLocLim[Z] = 1;
  }
  
  coorMatchLocLim[X] = wrapper->max_UNSIGNED_INT( coorMatchLocLim[X] );
  if ( coorMatchLocLim[X] ) InsertionArray->box.ptA[X] += geoshift;
  coorMatchLocLim[Y] = wrapper->max_UNSIGNED_INT( coorMatchLocLim[Y] );
  if ( coorMatchLocLim[Y] ) InsertionArray->box.ptA[Y] += geoshift;  
  coorMatchLocLim[Z] = wrapper->max_UNSIGNED_INT( coorMatchLocLim[Z] );
  if ( coorMatchLocLim[Z] ) InsertionArray->box.ptA[Z] += geoshift;
  
  if ( ( coorMatchLocLim[X] || coorMatchLocLim[Y] || coorMatchLocLim[Z] )
  	&& wrapper->get_rank() == 0 )
  {
    cout << endl << "Warning: Structured array positions: some coordinates"
    	<< " exactly match local linked cell grid limits in the following "
	<< "directions:" << endl;
    for (size_t i=0;i<3;++i)
      if ( coorMatchLocLim[i] )
        cout << "   * " << ( i == 0 ? "X" : i == 1 ? "Y" : "Z" ) << 
      	" automatic translation of " << geoshift << endl;	
  }
}
