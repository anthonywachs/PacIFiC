#include "MPINeighbors.hh"
#include "GrainsMPIWrapper.hh"


// ----------------------------------------------------------------------------
// Default constructor
MPINeighbors::MPINeighbors()  
{}




// ----------------------------------------------------------------------------
// Constructor with input parameters 
MPINeighbors::MPINeighbors( MPI_Comm &commgrainsMPI_3D, int const *center_coords,
  	int const* nprocsdir, int const* period )  
{
  struct MPICartInfos empty;
  int i, j, k, ii, jj, kk, l, pos, rank_;
  int *coords_ = new int[3];
  
  m_data.reserve(27); 
  for (i=0;i<27;++i) m_data.push_back(empty);
  for (i=0;i<27;++i) m_data[i].coords = new int[3];  

  m_nneighbors = 0;
  for (ii=0;ii<3;ii++)
    for (jj=0;jj<3;jj++)
      for (kk=0;kk<3;kk++)
      {
        i = ii-1;
	j = jj-1;
	k = kk-1;
	pos = 9 * ( k + 1 ) + 3 * ( j + 1 ) + i + 1;
	
	coords_[0] = center_coords[0] + i;
	coords_[1] = center_coords[1] + j;	
	coords_[2] = center_coords[2] + k;
	for (l=0;l<3;++l) 
	  if ( period[l] )
	  {
	    if ( coords_[l] < 0 ) coords_[l] += nprocsdir[l] ;
	    if ( coords_[l] > nprocsdir[l] - 1 ) coords_[l] -= nprocsdir[l] ;
	  }	
	
	if ( ( coords_[0] < 0 || coords_[0] > nprocsdir[0] - 1 )
		|| ( coords_[1] < 0 || coords_[1] > nprocsdir[1] - 1 )
		|| ( coords_[2] < 0 || coords_[2] > nprocsdir[2] - 1 ) )
	{
	  m_data[pos].rank = -1;
	  for (l=0;l<3;++l) m_data[pos].coords[l] = -1; 
	}
	else
	{
	  for (l=0;l<3;++l) m_data[pos].coords[l] = coords_[l];
	  MPI_Cart_rank( commgrainsMPI_3D, coords_, &rank_ );
          m_data[pos].rank = rank_;
	  if ( i != 0 || j != 0 || k != 0 ) 
	  {
	    m_rankNeighborsOnly.push_back( rank_ );
	    m_geolocNeighborsOnly.push_back( 
	    	GrainsMPIWrapper::getGeoPosition( i, j, k ) );
	  }
	  ++m_nneighbors;
	}
      }

  delete [] coords_;
}




// ----------------------------------------------------------------------------
// Destructor
MPINeighbors::~MPINeighbors()
{
  for (int i=0;i<27;++i) delete [] m_data[i].coords;
  m_data.clear();
}




// ----------------------------------------------------------------------------
// Returns the MPI rank of a neighbor based on a relative position (-1,0 ou +1) 
int MPINeighbors::rank( int i, int j, int k ) const
{
  return ( m_data[ 9 * ( k + 1 ) + 3 * ( j + 1 ) + i + 1 ].rank );
}




// ----------------------------------------------------------------------------
// Returns the coordinates in the Cartesian MPI topology of a 
// neighbor based on a relative position (-1,0 ou +1) 
int const* MPINeighbors::coordinates( int i, int j, int k ) const
{
  return ( m_data[ 9 * ( k + 1 ) + 3 * ( j + 1 ) + i + 1 ].coords );
}




// ----------------------------------------------------------------------------
// Display
ostream& operator <<( ostream &os, MPINeighbors const& M )
{
  int const* coords_;  
  
  os << "Number of true neighbors including itself = " << M.m_nneighbors 
  	<< endl;
  os << "Number of true neighbors without itself = " << 
  	M.m_rankNeighborsOnly.size() << endl;	
  for (int i=-1;i<2;i++)
    for (int j=-1;j<2;j++)
      for (int k=-1;k<2;k++)
      {  
        os << "Neighbor (" << i << "," << j << "," << k << ")" << endl;
	coords_ = M.coordinates( i, j, k );
	os << "Position in MPI topology = " << coords_[0] << " " <<
		coords_[1] << " " << coords_[2] << endl;
	os << "Rank = " << M.rank( i, j, k ) << endl;
	os << endl;
      }
      
  return ( os );

}




// ----------------------------------------------------------------------------
// Returns the MPI rank of all neighbors (including the process
// itself) in a array of integers
int* MPINeighbors::rankMPINeighbors() const
{
  vector<struct MPICartInfos>::const_iterator iv;
  int i=0;
  int* rv = new int[m_nneighbors];

  for (iv=m_data.begin();iv!=m_data.end();iv++)
    if ( iv->rank != -1 )
    {
      rv[i]=iv->rank;
      ++i;
    }

  return ( rv );

}




// ----------------------------------------------------------------------------
// Returns the MPI rank of all neighbors (excluding the process
// itself) in a list of integers 
list<int> const* MPINeighbors::rankMPINeighborsOnly() const
{
  return ( &m_rankNeighborsOnly );
}




// ----------------------------------------------------------------------------
// Returns the MPI geolocalisation of all neighbors (excluding the 
// process itself) in a list of GeoPosition
list<GeoPosition> const* MPINeighbors::geolocMPINeighborsOnly() const
{
  return ( &m_geolocNeighborsOnly );
}




// ----------------------------------------------------------------------------
// Returns the number of neighbors including the process itself
int MPINeighbors::nbMPINeighbors() const 
{ 
  return ( m_nneighbors );
}
