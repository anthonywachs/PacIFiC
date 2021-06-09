#include "Voisins.hh"
#include "MPIWrapperGrains.hh"


//-----------------------------------------------------------------------------
// Constructeur par defaut
Voisins::Voisins()  
{}




//-----------------------------------------------------------------------------
// Constructeur avec arguments
Voisins::Voisins( MPI_Comm &commgrainsMPI_3D, int const *center_coords,
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
	    	MPIWrapperGrains::getMPIGeoLocalisation( i, j, k ) );
	  }
	  ++m_nneighbors;
	}
      }

  delete [] coords_;

}




//-----------------------------------------------------------------------------
// Destructeur
Voisins::~Voisins()
{
  for (int i=0;i<27;++i) delete [] m_data[i].coords;
  m_data.clear();
}




//-----------------------------------------------------------------------------
// Rang
int Voisins::rank( int i, int j, int k ) const
{
  return m_data[ 9 * ( k + 1 ) + 3 * ( j + 1 ) + i + 1 ].rank;
}




//-----------------------------------------------------------------------------
// Coordinates
int const* Voisins::coordinates( int i, int j, int k ) const
{
  return m_data[ 9 * ( k + 1 ) + 3 * ( j + 1 ) + i + 1 ].coords;
}




//-----------------------------------------------------------------------------
// Display
ostream& operator <<( ostream &os, Voisins const& M )
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
      
  return os;

}




//-----------------------------------------------------------------------------
// Renvoie le rang des voisins (avec le processeur lui même) 
// dans un tableau d'entiers
int* Voisins::rangVoisins() const
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

  return rv;

}




//-----------------------------------------------------------------------------
// Renvoie le rang des voisins (sans le processeur lui même) 
list<int> const* Voisins::rangVoisinsSeuls() const
{
  return &m_rankNeighborsOnly;
}




//-----------------------------------------------------------------------------
// Renvoie la geolocalisation des voisins (sans le processeur lui même) 
list<MPIGeoLocalisation> const* Voisins::geolocVoisinsSeuls() const
{
  return &m_geolocNeighborsOnly;
}
