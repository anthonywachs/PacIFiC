#include "MPIWrapperGrains.hh"
#include "EnsComposant.H"
#include "Voisins.hh"
#include "LinkedCell.H"
#include "Cellule.H"
#include "FormeVdW.H"
#include "Box.H"
#include "Grains.H"
#include "Grains_Exec.hh"
#include <algorithm>


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Constructeur avec definition de la cellule unitaire. */
LinkedCell::LinkedCell()
  : AppSec(),
  m_nb( 0 ),
  m_nbi( 0 ),
  m_nbj( 0 ),
  m_nbk( 0 ),
  m_arete_X( 0. ),
  m_arete_Y( 0. ),
  m_arete_Z( 0. ),
  m_xmin( 0. ),
  m_ymin( 0. ),
  m_zmin( 0. ),
  m_xmax( 0. ),
  m_ymax( 0. ),
  m_zmax( 0. ),
  m_paveEtendu( NULL )
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
LinkedCell::~LinkedCell()
{
  vector<Cellule*>::iterator cellule;
  for (cellule=m_cellules.begin(); cellule!=m_cellules.end(); cellule++)
    delete *cellule;
  delete m_paveEtendu;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition de l'espace du Linked_Cell.
void LinkedCell::define( double arete )
{
  m_LC_origine_globale = m_origine_globale;
  m_LC_origine_locale = m_origine_locale;

  // Nombre de cellules dans les 3 directions
  // WRN : Protection si pas de particules presentes
  if ( arete > 1.e-10 )
  {
    m_nbi = (int)( ( App::m_lX + EPS ) / arete);
    m_nbi = m_nbi == 0 ? 1 : m_nbi;
    m_arete_X = App::m_lX / m_nbi ;
    m_nbj = (int)( ( App::m_lY + EPS ) / arete);
    m_nbj = m_nbj == 0 ? 1 : m_nbj;
    m_arete_Y = App::m_lY / m_nbj ;
    m_nbk = (int)( ( App::m_lZ + EPS ) / arete);
    m_nbk = m_nbk == 0 ? 1 : m_nbk;
    m_arete_Z = App::m_lZ / m_nbk ;
    if ( !m_nbk )
    {
      m_nbk = 1;
      m_arete_Z = m_arete_Y;
    }
  }
  else  m_nbi = m_nbj = m_nbk = 1;

  m_nb = m_nbi * m_nbj * m_nbk;

  cout << "Maillage Linked-cell" << endl;
  cout << "   Nombre de cellules = " << m_nbi << " " << m_nbj << " "
  	<< m_nbk << " = " << m_nbi * m_nbj * m_nbk << endl;
  cout << "   Taille des cellules = " << m_arete_X << " x " << m_arete_Y <<
  	" x " << m_arete_Z << endl << endl;

  m_xmin = m_LC_origine_locale[0];
  m_ymin = m_LC_origine_locale[1];
  m_zmin = m_LC_origine_locale[2];
  m_xmax = m_xmin + m_nbi * m_arete_X;
  m_ymax = m_ymin + m_nbj * m_arete_Y;
  m_zmax = m_zmin + m_nbk * m_arete_Z;
  m_paveEtendu = new BBox(
  	Point( m_xmin - 0.5 * m_arete_X,
		m_ymin - 0.5 * m_arete_Y,
		m_zmin - 0.5 * m_arete_Z ),
	Point( m_xmax + 0.5 * m_arete_X,
		m_ymax + 0.5 * m_arete_Y,
		m_zmax + 0.5 * m_arete_Z ) );

  // Construction des cellules
  m_cellules.reserve( m_nb );
  Cellule::SetEspace( m_nbi, m_nbj, m_nbk );
  for (int j=0; j<m_nbj; j++)
    for (int k=0; k<m_nbk; k++)
      for (int i=0; i<m_nbi; i++)
        m_cellules.push_back( new Cellule( getCelluleNumber( i, j, k ),
                i, j, k, m_arete_X, m_arete_Y, m_arete_Z,
                m_xmax, m_ymax, m_zmax ) );

  // Affectation du voisinage
  affecteVoisinage();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition de l'espace du Linked_Cell.
void LinkedCell::define( double arete, int const* nprocsdir,
	int const* MPIcoords, Voisins const* voisins, int const* MPIperiod )
{
  m_LC_origine_globale = m_origine_globale;
  m_LC_origine_locale = m_origine_locale;

  // Nombre de cellules dans les 3 directions
  // WRN : Protection si pas de particules presentes
  if (arete > 1.e-10)
  {
    m_nbi = (int)( ( App::m_lX + EPS ) / arete);
    m_nbi = m_nbi == 0 ? 1 : m_nbi;
    m_arete_X = App::m_lX / m_nbi ;
    m_nbj = (int)( ( App::m_lY + EPS ) / arete);
    m_nbj = m_nbj == 0 ? 1 : m_nbj;
    m_arete_Y = App::m_lY / m_nbj ;
    m_nbk = (int)( ( App::m_lZ + EPS ) / arete);
    m_nbk = m_nbk == 0 ? 1 : m_nbk;
    m_arete_Z = App::m_lZ / m_nbk ;
    if ( !m_nbk )
    {
      m_nbk = 1;
      m_arete_Z = m_arete_Y;
    }

    // Ajout de cellules pour les zones de recouvrement
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

    // Origine locale du LinkedCell
    if ( voisins->rank( -1, 0, 0 ) != -1 ) m_LC_origine_locale.Deplacer(
    	- m_arete_X, 0., 0. );
    if ( voisins->rank( 0, -1, 0 ) != -1 ) m_LC_origine_locale.Deplacer(
    	0., - m_arete_Y, 0. );
    if ( voisins->rank( 0, 0, -1 ) != -1 ) m_LC_origine_locale.Deplacer(
    	0., 0., - m_arete_Z );

    // Periodicite traitee dans la topologie MPI
    if ( MPIperiod[0] ) m_LC_origine_globale.Deplacer( - m_arete_X, 0., 0. );
    if ( MPIperiod[1] )	m_LC_origine_globale.Deplacer( 0., - m_arete_Y, 0. );
    if ( MPIperiod[2] )	m_LC_origine_globale.Deplacer( 0., 0., - m_arete_Z );
  }
  else  m_nbi = m_nbj = m_nbk = 1;

  m_nb = m_nbi * m_nbj * m_nbk;

  m_xmin = m_LC_origine_locale[0];
  m_ymin = m_LC_origine_locale[1];
  m_zmin = m_LC_origine_locale[2];
  m_xmax = m_xmin + m_nbi * m_arete_X;
  m_ymax = m_ymin + m_nbj * m_arete_Y;
  m_zmax = m_zmin + m_nbk * m_arete_Z;
  m_paveEtendu = new BBox(
  	Point( m_xmin - 0.5 * m_arete_X,
		m_ymin - 0.5 * m_arete_Y,
		m_zmin - 0.5 * m_arete_Z ),
	Point( m_xmax + 0.5 * m_arete_X,
		m_ymax + 0.5 * m_arete_Y,
		m_zmax + 0.5 * m_arete_Z ) );

  if ( voisins->rank( 0, 0, 0 ) == 0 )
  {
    cout << "Maillage Linked-cell sur le proc 0" << endl;
    cout << "   Nombre de cellules = " << m_nbi << " " << m_nbj << " "
    	<< m_nbk << " = " << m_nbi * m_nbj * m_nbk << endl;
    cout << "   Taille des cellules = " << m_arete_X << " x " << m_arete_Y <<
  	" x " << m_arete_Z << endl;
    cout << "   Origine globale = " << m_LC_origine_globale[X] << " " <<
    	m_LC_origine_globale[Y] << " " <<
	m_LC_origine_globale[Z] << endl;
    cout << "   Origine locale = " << m_LC_origine_locale[X] << " " <<
    	m_LC_origine_locale[Y] << " " <<
	m_LC_origine_locale[Z] << endl;
    cout << "   Limite locale du domaine = " << m_xmax << " " <<
    	 m_ymax << " " << m_zmax << endl << endl;
  }

  // Construction des cellules
  int tag = 0;
  MPIGeoLocalisation geoLoc = MPIGEO_NONE;
  m_cellules.reserve( m_nb );
  Cellule::SetEspace( m_nbi, m_nbj, m_nbk );
  for (int j=0; j<m_nbj; j++)
  {
    for (int k=0; k<m_nbk; k++)
    {
      for (int i=0; i<m_nbi; i++)
      {
        // Cas general
	geoLoc = MPIGEO_NONE;

	if ( i == 0 ) tag = 2;

	else if ( i == m_nbi - 1 ) tag = 2;

	else if ( i == 1 )
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == 1 || j == 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_SOUTH_WEST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 2 || j == m_nbj - 3 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_NORTH_WEST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_WEST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	       geoLoc = MPIGEO_NORTH_WEST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 1 ) tag = 2;
	  else
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_WEST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_WEST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_WEST;
	      }
	    }
	  }
	}

	else if ( i == m_nbi - 2 )
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == 1 || j == 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_SOUTH_EAST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 2 || j == m_nbj - 3 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_NORTH_EAST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 1 ) tag = 2;
	  else
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_EAST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_EAST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_EAST;
	      }
	    }
	  }
	}

	else if ( i == 2 )
	{
	  if ( j == 0 ) tag = 2;
	  else if ( j == 1 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_SOUTH_WEST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST;
	      }
	    }
	  }
	  else if ( j == 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_WEST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 3 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_WEST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_NORTH_WEST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_WEST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_WEST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 1 ) tag = 2;
	  else
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_WEST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_WEST_TOP;
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
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_SOUTH_EAST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST;
	      }
	    }
	  }
	  else if ( j == 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_EAST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 3 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_NORTH_EAST;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_EAST;
	      }
	    }
	  }
	  else if ( j == m_nbj - 1 ) tag = 2;
	  else
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_EAST_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_EAST_TOP;
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
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_SOUTH;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH;
	      }
	    }
	  }
	  else if ( j == 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_SOUTH_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 3 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else tag = 0;
	    }
	  }
	  else if ( j == m_nbj - 2 )
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 1;
	      geoLoc = MPIGEO_NORTH;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == 1 || k == 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 || k == m_nbk - 3 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH_TOP;
	      }
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else
	      {
	        tag = 1;
	        geoLoc = MPIGEO_NORTH;
	      }
	    }
	  }
	  else
	  {
	    // Geometrie 2D
	    if ( m_nbk == 1 )
	    {
	      tag = 0;
	    }
	    // Geometrie 3D
	    else
	    {
	      if ( k == 0 ) tag = 2;
	      else if ( k == m_nbk - 1 ) tag = 2;
	      else if ( k == 1 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_BOTTOM;
	      }
	      else if ( k == m_nbk - 2 )
	      {
	        tag = 1;
	        geoLoc = MPIGEO_TOP;
	      }
	      else tag = 0;
	    }
	  }
	}
	m_cellules.push_back( new Cellule( getCelluleNumber( i, j, k ),
		i, j, k, m_arete_X, m_arete_Y, m_arete_Z, tag,
		m_LC_origine_locale,
		m_xmax, m_ymax, m_zmax, geoLoc ) );
      }
    }
  }

  // Cas particuliers
  // Amont X
  if ( voisins->rank( -1, 0, 0 ) == -1 )
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<m_nbj; j++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCellule( i, j, k )->m_tag = getCellule( 3, j, k )->m_tag;
	  getCellule( i, j, k )->m_GeoLocCell =
	  	getCellule( 3, j, k )->m_GeoLocCell;
	}
  }

  // Aval X
  if ( voisins->rank( 1, 0, 0 ) == -1 )
  {
    for (int i=m_nbi-3; i<m_nbi; i++)
      for (int j=0; j<m_nbj; j++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCellule( i, j, k )->m_tag = getCellule( m_nbi-4, j, k )->m_tag;
	  getCellule( i, j, k )->m_GeoLocCell =
	  	getCellule( m_nbi-4, j, k )->m_GeoLocCell;
	}
  }

  // Amont Y
  if ( voisins->rank( 0, -1, 0 ) == -1 )
  {
    for (int j=0; j<3; j++)
      for (int i=0; i<m_nbi; i++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCellule( i, j, k )->m_tag = getCellule( i, 3, k )->m_tag;
	  getCellule( i, j, k )->m_GeoLocCell =
	  	getCellule( i, 3, k )->m_GeoLocCell;
	}
  }

  // Aval Y
  if ( voisins->rank( 0, 1, 0 ) == -1 )
  {
    for (int j=m_nbj-3; j<m_nbj; j++)
      for (int i=0; i<m_nbi; i++)
        for (int k=0; k<m_nbk; k++)
	{
	  getCellule( i, j, k )->m_tag = getCellule( i, m_nbj-4, k )->m_tag;
	  getCellule( i, j, k )->m_GeoLocCell =
	  	getCellule( i, m_nbj-4, k )->m_GeoLocCell;
	}
  }

  // Geometrie 3D
  if ( m_nbk > 1 )
  {
    // Amont Z
    if ( voisins->rank( 0, 0, -1 ) == -1 )
    {
      for (int k=0; k<3; k++)
        for (int i=0; i<m_nbi; i++)
          for (int j=0; j<m_nbj; j++)
	  {
	    getCellule( i, j, k )->m_tag = getCellule( i, j, 3 )->m_tag;
	    getCellule( i, j, k )->m_GeoLocCell =
	    	getCellule( i, j, 3 )->m_GeoLocCell;
	  }
    }

    // Aval Z
    if ( voisins->rank( 0, 0, 1 ) == -1 )
    {
      for (int k=m_nbk-3; k<m_nbk; k++)
        for (int i=0; i<m_nbi; i++)
          for (int j=0; j<m_nbj; j++)
	  {
	    getCellule( i, j, k )->m_tag = getCellule( i, j, m_nbk-4 )->m_tag;
	    getCellule( i, j, k )->m_GeoLocCell =
	    	getCellule( i, j, m_nbk-4 )->m_GeoLocCell;
	  }
    }
  }

  // Affectation du voisinage
  affecteVoisinage();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// D�finition le voisinage des cellules
void LinkedCell::affecteVoisinage()
{
  // Affectation du voisinage pour les contacts
  // 3 cellules au-dessus => i-1<->i+1, j,   k+1
  // 1 cellule  a droite  => i+1,       j,   k
  // 9 cellules arriere   => i-1<->i+1, j+1, k+1<-->k-1
  Cellule *cellule = NULL;
  Cellule *voisine = NULL;
  int icel, jcel, kcel;
  int cel[3][13] = {{-1,  0, +1, +1, -1,  0, +1, -1,  0, +1, -1,  0, +1},
		    //{ 0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		    { 0,  0,  0,  0, +1, +1, +1, +1, +1, +1, +1, +1, +1},
		    {+1, +1, +1,  0, +1, +1, +1,  0,  0,  0, -1, -1, -1}};

  for (int i=0; i<m_nb; i++)
  {
    cellule = m_cellules[i];
    icel = (*cellule)[X];
    jcel = (*cellule)[Y];
    kcel = (*cellule)[Z];

    for (int j=0; j<13; j++)
    {
      voisine = getCellule( icel+cel[X][j], jcel+cel[Y][j], kcel+cel[Z][j] );
      if ( voisine ) cellule->addVoisineContact( voisine );
    }
  }

//   // Affectation du voisinage complet
//   for (int i=0; i<m_nb; i++)
//   {
//     cellule = m_cellules[i];
//     icel = (*cellule)[X];
//     jcel = (*cellule)[Y];
//     kcel = (*cellule)[Z];
//
//     for (int k=-1;k<2;++k)
//       for (int l=-1;l<2;++l)
//         for (int m=-1;m<2;++m)
//           if ( k || l || m )
// 	  {
// 	    voisine = getCellule( icel+k, jcel+l, kcel+m );
//             if ( voisine ) cellule->addVoisine( voisine );
// 	  }
//   }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation du voisinage complet
// D. RAKOTONIRINA - Dec 2014 - Creation
void LinkedCell::affecteVoisinageComplet()
{
  Cellule *cellule = NULL;
  Cellule *voisine = NULL;
  int icel, jcel, kcel;
  for (int i=0; i<m_nb; i++)
  {
    cellule = m_cellules[i];
    icel = (*cellule)[X];
    jcel = (*cellule)[Y];
    kcel = (*cellule)[Z];

    for (int k=-1;k<2;++k)
      for (int l=-1;l<2;++l)
        for (int m=-1;m<2;++m)
          if ( k || l || m )
          {
            voisine = getCellule( icel+k, jcel+l, kcel+m );
            if ( voisine ) cellule->addVoisine( voisine );
          }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces d'interaction.
// D. RAKOTONIRINA - Juil 2014 - Modification
void LinkedCell::CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules )
{
  Particule* reference;
  Cellule* cellule;
  list<Particule*> voisines;

  // Contact entre particules
  Point centre;
  list<Particule*>::const_iterator particule;
  list<Particule*>::iterator voisine;
  list<Cellule*>::iterator around;
  int id[3];

  for( particule=particules->begin(); particule!=particules->end();
       particule++)
  {
    reference = *particule;

    // Recherche des particules voisines
    // Dans la cellule locale:
    // On recherche les voisines situees au dela dans la liste
    // Cela permet de ne considerer un contact entre 2 particules
    // Exemple: soit 4 particules dans la cellule aux positions 1,2,3,4
    // Si on s'interesse � la particule 3, le seul contact possible est 3-4
    // Puis quand on considere 4 plus tard dans la boucle, plus aucun contact
    // possible n'est recherche et donc le couple 3-4 = 4-3 n'est pris en compte
    // qu'une fois
    // Dans les cellules adjacentes:
    // Grace au voisinage de contact (9 derriere, 3 au-dessus, 1 droite), les
    // couples de contact possibles ne sont egalement compt�s qu'une seule fois
    centre = *(*particule)->getPosition();
    Cellule::GetCellule( centre, id );
    cellule = getCellule( id[X], id[Y], id[Z] );
    voisine = find( cellule->m_particules.begin(), cellule->m_particules.end(),
		   reference );

    voisine++;
    for( ; voisine!=cellule->m_particules.end(); voisine++ )
      voisines.push_back(*voisine);

    for( around=cellule->m_voisinesContact.begin();
    	 around!=cellule->m_voisinesContact.end(); around++ ) {
      for( voisine=(*around)->m_particules.begin();
           voisine!=(*around)->m_particules.end(); voisine++ ) {
        voisines.push_back(*voisine);
      }
    }

    // Evaluation des contacts avec les particules voisines
    // Inversion de l'appel dans le cas de particule composite
    for( voisine=voisines.begin(); voisine!=voisines.end(); voisine++ )
      if( (*voisine)->isCompParticule() )
        (*voisine)->InterAction( reference, dt, time, this );
      else
        reference->InterAction( *voisine, dt, time, this );

    voisines.clear();

    // Evaluation des contacts avec les obstacles
    list<MonObstacle*>::iterator myObs;
    for( myObs=cellule->m_obstacles.begin();
         myObs!=cellule->m_obstacles.end(); myObs++ )
    {
      // Inversion de l'appel dans le cas de particule composite
      // ObstaclePeriodique ne fait pas l'appel
      if ( reference->isCompParticule() && !(*myObs)->isObsPeriodique() )
        reference->InterAction( *myObs, dt, time, this );
      else
        (*myObs)->InterAction( reference, dt, time, this );
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces d'interaction pour l'initialisation des conditions de
// cohesion. ( based on LinkedCell::CalculerForces )
void LinkedCell::InitializeCohesiveForces( Scalar time, Scalar dt,
    list<Particule*> const* particules )
{
  Particule* reference;
  Cellule* cellule;
  list<Particule*> voisines;

  // Contact entre particules
  Point centre;
  list<Particule*>::const_iterator particule;
  list<Particule*>::iterator voisine;
  list<Cellule*>::iterator around;
  int id[3];

  for (particule=particules->begin(); particule!=particules->end();
       particule++)
  {
    reference = *particule;

    // Recherche des particules voisines
    // Dans la cellule locale:
    // On recherche les voisines situees au dela dans la liste
    // Cela permet de ne considerer un contact entre 2 particules
    // Exemple: soit 4 particules dans la cellule aux positions 1,2,3,4
    // Si on s'interesse � la particule 3, le seul contact possible est 3-4
    // Puis quand on considere 4 plus tard dans la boucle, plus aucun contact
    // possible n'est recherche et donc le couple 3-4 = 4-3 n'est pris en compte
    // qu'une fois
    // Dans les cellules adjacentes:
    // Grace au voisinage de contact (9 derriere, 3 au-dessus, 1 droite), les
    // couples de contact possibles ne sont egalement compt�s qu'une seule fois
    centre = *(*particule)->getPosition();
    Cellule::GetCellule( centre, id );
    cellule = getCellule( id[X], id[Y], id[Z] );
    voisine = find( cellule->m_particules.begin(), cellule->m_particules.end(),
        reference );

    voisine++;
    for ( ; voisine!=cellule->m_particules.end(); voisine++)
      voisines.push_back(*voisine);

    for( around=cellule->m_voisinesContact.begin();
         around!=cellule->m_voisinesContact.end(); around++)
      for( voisine=(*around)->m_particules.begin();
           voisine!=(*around)->m_particules.end(); voisine++)
        voisines.push_back(*voisine);

    // Evaluation des contacts avec les particules voisines
    // Inversion de l'appel dans le cas de particule composite
    for (voisine=voisines.begin(); voisine!=voisines.end(); voisine++)
    {
      if( (*voisine)->isCompParticule() )
        (*voisine)->InterActionCohesiveInit( reference, dt, time, this );
      else
        reference->InterActionCohesiveInit( *voisine, dt, time, this );
    }

    voisines.clear();

    // A FAIRE !!!

/*     // Evaluation des contacts avec les obstacles si cohesion avec obstacles
     list<MonObstacle*>::iterator myObs;
     for( myObs=cellule->m_obstacles.begin();
          myObs!=cellule->m_obstacles.end(); myObs++ )
     {
       // Inversion de l'appel dans le cas de particule composite
       // ObstaclePeriodique ne fait pas l'appel
       if ( reference->isCompParticule() && !(*myObs)->isObsPeriodique() )
         reference->InterActionBIS( *myObs, dt, time, this );
       else
         (*myObs)->InterActionBIS( reference, dt, time, this );
     }
*/
  }
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces d'interaction pour les clones periodiques
void LinkedCell::CalculerForcesClonesPeriodiques( Scalar time, Scalar dt,
  	list<Particule*>* particulesClonesPeriodiques )
{
  Particule* reference;
  Cellule* cellule;
  list<Particule*> voisines;

  // Contact entre particules
  Point centre;
  list<Particule*>::iterator particule;
  list<Particule*>::iterator voisine;
  list<Cellule*>::iterator   around;
  int id[3];

  for (particule=particulesClonesPeriodiques->begin();
  	particule!=particulesClonesPeriodiques->end();
       	particule++)
  {
    reference = *particule;

    // Recherche des particules voisines
    // Dans la cellule locale:
    // On recherche les voisines situees au dela dans la liste
    // Cela permet de ne considerer un contact entre 2 particules
    // Exemple: soit 4 particules dans la cellule aux positions 1,2,3,4
    // Si on s'interesse � la particule 3, le seul contact possible est 3-4
    // Puis quand on considere 4 plus tard dans la boucle, plus aucun contact
    // possible n'est recherche et donc le couple 3-4 = 4-3 n'est pris en compte
    // qu'une fois
    // Dans les cellules adjacentes:
    // Grace au voisinage de contact (9 derriere, 3 au-dessus, 1 droite), les
    // couples de contact possibles ne sont egalement compt�s qu'une seule fois
    centre = *(*particule)->getPosition();
    Cellule::GetCellule( centre, id );
    cellule = getCellule( id[X], id[Y], id[Z] );
    voisine = find( cellule->m_particules.begin(), cellule->m_particules.end(),
		   reference );

    voisine++;
    for ( ; voisine!=cellule->m_particules.end(); voisine++)
      voisines.push_back(*voisine);

    for (around=cellule->m_voisinesContact.begin();
    	around!=cellule->m_voisinesContact.end(); around++) {
      for (voisine=(*around)->m_particules.begin();
	   voisine!=(*around)->m_particules.end(); voisine++) {
	voisines.push_back(*voisine);
      }
    }

    // Evaluation des contacts avec les particules voisines
    // Inversion de l'appel dans le cas de particule composite
    for (voisine=voisines.begin(); voisine!=voisines.end(); voisine++)
      if ( (*voisine)->isCompParticule() )
	(*voisine)->InterAction( reference, dt, time, this );
      else
	reference->InterAction( *voisine, dt, time, this );

    voisines.clear();

    // Pas d'evaluation des contacts avec les obstacles
    // car la particule est un clone periodique: on n'evalue pas les contacts
    // avec les obstacles car ces efforts sont transmis � la particule
    // reference et ainsi les contacts seraient compt�s 2 fois
  }
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces d'interaction.
// D. RAKOTONIRINA - Juil 2014 - Modification
list<struct PointForcePostProcessing>* LinkedCell::CalculerForcesPostProcessing(
	list<Particule*> const* particules, Scalar dt ) const
{
  Particule* reference;
  Cellule* cellule;
  list<Particule*> voisines;
  list<PointForcePostProcessing>* listOfContacts =
  	new list<PointForcePostProcessing>;

  // Contact entre particules
  Point centre;
  list<Particule*>::const_iterator particule;
  list<Particule*>::iterator voisine;
  list<Cellule*>::iterator   around;
  int id[3];

  for (particule=particules->begin(); particule!=particules->end();
       particule++) {
    reference = *particule;
    centre = *(*particule)->getPosition();
    Cellule::GetCellule( centre, id );
    cellule = getCellule( id[X], id[Y], id[Z] );
    voisine = find( cellule->m_particules.begin(), cellule->m_particules.end(),
		   reference );

    voisine++;
    for ( ; voisine!=cellule->m_particules.end(); voisine++)
      voisines.push_back(*voisine);

    for (around=cellule->m_voisinesContact.begin();
    	around!=cellule->m_voisinesContact.end(); around++) {
      for (voisine=(*around)->m_particules.begin();
	   voisine!=(*around)->m_particules.end(); voisine++) {
	voisines.push_back(*voisine);
      }
    }

    // Evaluation des contacts avec les particules voisines
    // Inversion de l'appel dans le cas de particule composite
    for ( voisine=voisines.begin(); voisine!=voisines.end(); voisine++ )
      if ( (*voisine)->isCompParticule() )
	(*voisine)->InterActionPostProcessing( reference, dt, listOfContacts );
      else
	reference->InterActionPostProcessing( *voisine, dt, listOfContacts );

    voisines.clear();

    // Evaluation des contacts avec les obstacles
    // Inversion de l'appel dans le cas de particule composite
    list<MonObstacle*>::iterator myObs;
    for (myObs=cellule->m_obstacles.begin();
      myObs!=cellule->m_obstacles.end();myObs++)
      if ( reference->isCompParticule() && !(*myObs)->isObsPeriodique() )
        reference->InterActionPostProcessing( *myObs, dt, listOfContacts );
      else
        (*myObs)->InterActionPostProcessing( reference, dt, listOfContacts );
  }

  return listOfContacts;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reference de la cellule liee a la position indiquee.
Cellule* LinkedCell::getCellule( int i, int j, int k ) const
{
  Cellule* cellule = NULL;
  bool valid = -1 < i && i < m_nbi && -1 < j && j < m_nbj
  	&& -1 < k && k < m_nbk;
  if ( valid ) cellule =  m_cellules[getCelluleNumber( i, j, k )];
  return cellule;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cellule associee a la position designee
Cellule* LinkedCell::getCellule( const Point &position ) const
{
  int id[3];
  Cellule::GetCellule( position, id );
  return getCellule( id[X], id[Y], id[Z] );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Numero de la cellule liee a la position indiquee.
int LinkedCell::getCelluleNumber( int i, int j, int k ) const
{
  return ( j * m_nbk * m_nbi + k * m_nbi + i ) ;
}




// ----------------------------------------------------------------------------
// La particule indiquee est elle en contact avec une autre particule
// G.FERRER - Aout.2002 - Creation
bool LinkedCell::isContact( const Particule* reference ) const
{
  bool contact = false;

  Point centre = *reference->getPosition();
  int id[3];
  Cellule::GetCellule( centre, id );
  Cellule* cellule = getCellule( id[X], id[Y], id[Z] );
  Cellule* voisine = NULL ;

  if ( cellule )
  {
    // Verification avec les particules de la cellule
    // Plus forte probabilite de contact => verification 1
    contact = cellule->isContact( reference );

//     // Verification avec l'ensemble des cellules voisines
//     list<Cellule*>::iterator icell;
//     for (icell=cellule->m_allVoisines.begin();
//     	icell!=cellule->m_allVoisines.end() && !contact; icell++)
//       contact = (*icell)->isContact( reference );

    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    voisine = getCellule( id[X]+k, id[Y]+l, id[Z]+m );
            if ( voisine ) contact = voisine->isContact( reference );
	  }
  }

  return contact;
}




// ----------------------------------------------------------------------------
// La particule indiquee est elle en contact avec une autre particule, variante
// VdW
// AWACHS - Janv.2011 - Creation
bool LinkedCell::isContactVdW( const Particule* reference ) const
{
  bool contact = false;

  Point centre = *reference->getPosition();
  int id[3];
  Cellule::GetCellule( centre, id );
  Cellule* cellule = getCellule( id[X], id[Y], id[Z] );
  Cellule* voisine = NULL ;

  if ( cellule )
  {
    // Verification avec les particules de la cellule
    // Plus forte probabilite de contact => verification 1
    contact = cellule->isContactVdW( reference );

    // Verification avec l'ensemble des cellules voisines
//     list<Cellule*>::iterator icell;
//     for (icell=cellule->m_allVoisines.begin();icell!=cellule->m_allVoisines.end()
//   	&& !contact; icell++)
//       contact = (*icell)->isContactVdW( reference );

    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    voisine = getCellule( id[X]+k, id[Y]+l, id[Z]+m );
            if ( voisine ) contact = voisine->isContactVdW( reference );
	  }
  }
  return contact;
}




// ----------------------------------------------------------------------------
// La particule indiquee est elle proche d'une autre particule
// AWACHS - Janv.2011 - Creation
bool LinkedCell::isProche( const Particule* reference ) const
{
  bool contact = false;

  Point centre = *reference->getPosition();
  int id[3];
  Cellule::GetCellule( centre, id );
  Cellule* cellule = getCellule( id[X], id[Y], id[Z] );
  Cellule* voisine = NULL ;

  if ( cellule )
  {
    // Verification avec les particules de la cellule
    // Plus forte probabilite de contact => verification 1
    contact = cellule->isProche( reference );

    // Verification avec l'ensemble des cellules voisines
//     list<Cellule*>::iterator icell;
//     for (icell=cellule->m_allVoisines.begin();icell!=cellule->m_allVoisines.end()
//   	&& !contact; icell++)
//       contact = (*icell)->isProche( reference );
    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    voisine = getCellule( id[X]+k, id[Y]+l, id[Z]+m );
            if ( voisine ) contact = voisine->isProche( reference );
	  }

  }

  return contact;
}




// ----------------------------------------------------------------------------
// La particule indiquee est elle proche d'une autre particule, variante
// VdW
// AWACHS - Janv.2011 - Creation
bool LinkedCell::isProcheVdW( const Particule* reference ) const
{
  bool contact = false;

  Point centre = *reference->getPosition();
  int id[3];
  Cellule::GetCellule( centre, id );
  Cellule* cellule = getCellule( id[X], id[Y], id[Z] );
  Cellule* voisine = NULL ;

  if ( cellule )
  {
    // Verification avec les particules de la cellule
    // Plus forte probabilite de contact => verification 1
    contact = cellule->isProcheVdW( reference );

    // Verification avec l'ensemble des cellules voisines
//     list<Cellule*>::iterator icell;
//     for (icell=cellule->m_allVoisines.begin();icell!=cellule->m_allVoisines.end()
//   	&& !contact; icell++)
//       contact = (*icell)->isProcheVdW( reference );
    for (int k=-1;k<2 && !contact;++k)
      for (int l=-1;l<2 && !contact;++l)
        for (int m=-1;m<2 && !contact;++m)
          if ( k || l || m )
	  {
	    voisine = getCellule( id[X]+k, id[Y]+l, id[Z]+m );
            if ( voisine ) contact = voisine->isProcheVdW( reference );
	  }
  }

  return contact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Association de la particule avec l'algorithme.
void LinkedCell::Link( Particule* particule )
{
  Point centre = *particule->getPosition();
  int id[3];
  Cellule::GetCellule( centre, id );
  Cellule* cellule = getCellule( id[X], id[Y], id[Z] );
  particule->setCelluleNm1( cellule );
  particule->setTag( cellule->m_tag );
  particule->setGeoLocalisation( cellule->m_GeoLocCell );

  if ( !cellule->Contient(particule) ) cellule->add( particule );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Association de l'obstacle avec l'algorithme.
void LinkedCell::Link( Obstacle* obstacle )
{
  AppSec::Link(obstacle);
  list<MonObstacle*> list_obstacles = obstacle->getObstacles();
  list<MonObstacle*>::iterator myObs;
  Cellule* cellule = NULL;
  Transform CelPosition;
  double alpha=2.;
  Point const* cg = NULL;

  // Sachant que les cellules ont une taille de 2*rayon_max, pour attacher
  // un obstacle � une cellule, on construit une boite de taille 2 fois
  // la taille de la cellule, soit cette boite est plus grande de rayon_max
  // dans toutes les directions. Ainsi, si une particule appartient � une
  // cellule, le point de sa surface le plus eloigne de la cellule est au max
  // a une distance rayon_max, ce qui garantit qu'aucune collision particule-
  // obstacle n'est oubliee.
  // Dans le cas d'obstacle periodique, on augmente la taille � 3.1 fois
  // la taille de la cellule
  for (myObs=list_obstacles.begin();myObs!=list_obstacles.end();myObs++)
  {
    const Forme* obstacleForme = (*myObs)->getForme();
    BBox obsBox = obstacleForme->BoxForme();
    if ( (*myObs)->materiau() == "periode" ) alpha = 3.1;
    else alpha = 2.;
    Convex* convexCellule = new Box( alpha * m_arete_X, alpha * m_arete_Y,
    	alpha * m_arete_Z);
    Forme CelForme( convexCellule, CelPosition );

    // Intersection g�om�trique de la cellule avec l'obstacle
    for (int i=0; i<m_nb; i++)
    {
      cellule = m_cellules[i];
      cg = cellule->getCentre();
      if ( obsBox.InZone( cg, 0.5 * alpha * m_arete_X, 0.5 * alpha * m_arete_Y,
      		0.5 * alpha * m_arete_Z ) )
      {
        CelForme.setOrigin( (*cg)[X], (*cg)[Y], (*cg)[Z] );
	if ( CelForme.isContact( *obstacleForme ) )
        {
          cellule->addObstacle( *myObs );
          (*myObs)->add( cellule );
        }
      }
    }

    // Rem: il ne faut pas detruire le convex convexCellule car le destructeur
    // de CelForme, objet de type Forme, s'en charge (cf Forme.cpp)
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Actualisation des associations des particules avec l'algorithme.
void LinkedCell::LinkUpdate( Scalar time, Scalar dt,
  	list<Particule*>* particules) throw(ErreurSimulation)
{
  // Si la particule n'est plus active, on la supprime de l'algorithme
  list<Particule*>::iterator particule;
  for (particule=particules->begin(); particule!=particules->end(); )
  {
    switch ( (*particule)->getActivity() )
    {
      case COMPUTE:
        particule++;
        break;

      default:
        (*particule)->getCelluleNm1()->remove( *particule );
        particule = particules->erase( particule );
        break;
    }
  }

  // Les obstacles etant traites dans le LinkedCell, on traite la mise � jour
  // de leurs liens avec les cellules
  list<MonObstacle*>::iterator myObs;
  for (myObs=m_allObstacles.begin();myObs!=m_allObstacles.end();myObs++)
    // Est ce que l'obstacle a boug�: si oui mise � jour, sinon ne rien faire
    if ( (*myObs)->hasMoved() )
    {
      // Est ce que l'obstacle a potentiellement une intersection avec le
      // LinkedCell ?
      if ( intersect( *(*myObs)->getObstacleBox() , *m_paveEtendu ) )
      {
	// Si l'obstacle n'a pas encore ete link� sur ce LinkedCell (sous-
	// domaine ou processeur) => Link classique
	// sinon => Update optimis� sur les cellules voisines
	if ( (*myObs)->getInCells()->empty() ) Link( *myObs );
	else LinkUpdate( time, dt, *myObs );
      }
    }

  for (particule=particules->begin(); particule!=particules->end();
	 particule++)
    LinkUpdateActiveParticule( *particule );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Actualisation d'une particule active avec l'algorithme
// A.WACHS - Oct.2009 - Creation
void LinkedCell::LinkUpdateActiveParticule( Particule* particule )
	throw (ErreurSimulation)
{
  if ( particule->getActivity() != COMPUTE )
  {
      cout << "\nParticule non active " << particule->getID() << endl;
      cout << "            " << *particule->getPosition() << endl;
      Grains_Exec::m_exception_Simulation = true;
      throw(ErreurSimulation("LinkedCell::LinkUpdateActiveParticule"));
  }

  // Cellule dans laquelle se trouvait la particule au temps pr�c�dent
  Cellule* celluleNm1 = particule->getCelluleNm1();

  // Cellule actuelle de la particule
  Point centre = *particule->getPosition();
  int id[3];
  Cellule::GetCellule( centre, id );
  Cellule* celluleNew = getCellule( id[X], id[Y], id[Z] );
  particule->setCelluleNm1( celluleNew );

  if ( celluleNew == NULL )
  {
    cout << "\nParticule " << particule->getID()       << endl;
    cout << "            " << *particule->getPosition() << endl;
    Grains_Exec::m_exception_Simulation = true;
    throw(ErreurSimulation("LinkedCell::LinkUpdateActiveParticule"));
  }

  // Si la particule n'est pas dans la cellule actuelle,
  // on l'installe & on la supprime de sa cellule precedente
  if ( celluleNew != celluleNm1 )
  {
    celluleNew->add( particule );
    celluleNm1->remove( particule );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise a jour de l'association de l'obstacle avec l'algorithme.
// A.WACHS - Aout.2009 - Creation
void LinkedCell::LinkUpdate( Scalar time, Scalar dt, MonObstacle *myObs )
{
  if ( myObs->performLinkUpdate() )
  {
    const Forme* obstacleForme = myObs->getForme();
    BBox obsBox = obstacleForme->BoxForme();
    Cellule* cellule = NULL;
    Vecteur deplMax;
    Point const* cg = NULL;

    // le coefficient 1.2 donne une marge d'erreur de 20%, ce qui signifie
    // qu'on suppose que le vecteur vitesse de l'obstacle sur les n pas de temps
    // suivants ne varie pas de plus de 20%
    // Attention: rien dans le code verifie cette hypoth�se !!
    int updateFreq = myObs->getObstacleLinkedCellUpdateFrequency();
    double coefApprox = updateFreq == 1 ? 1. : 1.2;
    deplMax = myObs->vitesseMaxPerDirection() * coefApprox
   	* updateFreq * dt;

    Vecteur CelExtent( m_arete_X + deplMax[X], m_arete_Y + deplMax[Y],
    	m_arete_Z + deplMax[Z] );

    myObs->resetInCells();
    for (int i=0; i<m_nb; i++)
    {
      cellule = m_cellules[i];
      cg = cellule->getCentre();
      if ( obsBox.InZone( cg, CelExtent[X], CelExtent[Y], CelExtent[Z] ) )
      {
        cellule->addObstacle( myObs );
        myObs->add( cellule );
      }
    }
  }

//   const list<Cellule*>* voisinageCourant = myObs->getInCells();
//   list<Cellule*> voisinageEtendu = *voisinageCourant;
//   const list<Cellule*>*	celluleVoisinageComplet = NULL;
//   list<Cellule*>::const_iterator icellule,icelluleVoisine;
//   list<Cellule*>::iterator il;
//   const Forme* obstacleForme = myObs->getForme();
//   BBox obsBox = obstacleForme->BoxForme();
//   Cellule* cellule = NULL;
//   Transform CelPosition;
//   double alpha = 2.;
//   Convex* convexCellule = new Box(alpha*m_arete_X,alpha*m_arete_Y,
//     	alpha*m_arete_Z);
//   Forme CelForme(convexCellule,CelPosition);
//
//   // Voisinage etendu
//   for (icellule=voisinageCourant->begin();icellule!=voisinageCourant->end();
//   	icellule++)
//   {
//     celluleVoisinageComplet = (*icellule)->getVoisinageComplet();
//     for (icelluleVoisine=celluleVoisinageComplet->begin();
//     	icelluleVoisine!=celluleVoisinageComplet->end();icelluleVoisine++)
//       voisinageEtendu.push_back(*icelluleVoisine);
//   }
//   voisinageEtendu.sort();
//   voisinageEtendu.unique();
//
//   // Intersection g�om�trique de la cellule avec l'obstacle
//   // dans le voisinage etendu
//   myObs->resetInCells();
//   for (il=voisinageEtendu.begin();il!=voisinageEtendu.end();il++)
//   {
//     cellule = *il;
//     Point cg = cellule->Gravite();
//     if (obsBox.InZone(cg,m_arete_X,m_arete_Y,m_arete_Z))
//     {
//       CelForme.setOrigin(cg[X],cg[Y],cg[Z]);
//       if (CelForme.isContact(*obstacleForme))
//       {
//         cellule->addObstacle(myObs);
//         myObs->add(cellule);
//       }
//     }
//   }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Suppression de la particule de l'APPlication
void LinkedCell::remove( Particule* particule )
{
  Point centre = *particule->getPosition();
  int id[3];
  Cellule::GetCellule( centre, id );
  Cellule* cellule = getCellule( id[X], id[Y], id[Z] );
  if ( cellule ) cellule->remove( particule );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Suppression de l'obstacle de l'APPlication
void LinkedCell::remove( MonObstacle* obs )
{
  AppSec::remove( obs );
  obs->resetInCells();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Transfert de la particule dans le LinkedCell lie a son deplacement
void LinkedCell::shift( Particule& particule,
	const Vecteur& vecteur )
{
  Point ancien  = *particule.getPosition();
  Point nouveau = ancien + vecteur;

  Cellule* celluleA = getCellule( ancien );
  Cellule* celluleB = getCellule( nouveau );
  if ( celluleA != celluleB )
  {
    celluleA->remove( &particule );
    celluleB->add( &particule );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Operateur <<
ostream& operator <<( ostream &f, const LinkedCell &LC )
{
  vector<Cellule*>::const_iterator iv;
  int nbCelSansObstacle = 0;
  for (iv=LC.m_cellules.begin();iv!=LC.m_cellules.end();iv++)
    if ( (*iv)->nombreObstacles() == 0 ) nbCelSansObstacle++;

  f << "Nb total de cellules = " << LC.m_nb << endl;
  f << "Nb de cellules sans obstacles dans leur voisinage = " <<
  	nbCelSansObstacle << endl;
  f << "Nb de cellules en X x Y x Z = " << LC.m_nbi << " x " << LC.m_nbj
  	<< " x " << LC.m_nbk << endl;
  f << "Taille des cellules en X x Y x Z = " << LC.m_arete_X << " x "
  	<< LC.m_arete_Y << " x " << LC.m_arete_Z << endl;
  f << "Origine locale de la grille = " << LC.m_LC_origine_locale;
  f << "Limite de la grille = " << LC.m_xmax << " " << LC.m_ymax << " " <<
  	LC.m_zmax << endl;
  f << "CELLULES" << endl;
  for (iv=LC.m_cellules.begin();iv!=LC.m_cellules.end();iv++)
    f << *(*iv) << endl << endl;

  return f;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Appartenance du centre de gravite d'une particule au Linked Cell
bool LinkedCell::isInLinkedCell( const Point &position ) const
{
  bool isIn = true;

  if ( position[0] < m_xmin || position[0] > m_xmax
  	|| position[1] < m_ymin || position[1] > m_ymax
  	|| position[2] < m_zmin || position[2] > m_zmax ) isIn = false;

  return isIn;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Appartenance du centre de gravite d'une particule au Linked Cell
bool LinkedCell::isInLinkedCell( const double &gx, const double &gy,
	const double &gz ) const
{
  bool isIn = true;

  if ( gx < m_xmin || gx > m_xmax
  	|| gy < m_ymin || gy > m_ymax
  	|| gz < m_zmin || gz > m_zmax ) isIn = false;

  return isIn;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise a jour des tags des particules interieures
void LinkedCell::updateInteriorTag( Scalar time,
	list<Particule*>* ENSparticules,
	list<Particule*>* particulesHalozone,
	MPIWrapperGrains const* wrapper )
{
  list<Particule*>::iterator particule;
  int tag_nm1,tag;
  Cellule* current_cell = NULL;

  for (particule=ENSparticules->begin(); particule!=ENSparticules->end();
	 particule++)
  {
    tag_nm1 = (*particule)->getTag();
    current_cell = getCellule( *(*particule)->getPosition() );
    if ( tag_nm1 == 0 )
    {
      tag = (*particule)->setTag( current_cell->m_tag );

      // Interior to Halozone (0 -> 1)
      if ( tag == 1 )
      {
        if ( Grains_Exec::m_MPI_verbose )
	{
	  ostringstream oss;
	  oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT ) <<
		" Interior to Halozone (0 -> 1)               Id = " <<
      		(*particule)->getID() << " " << *(*particule)->getPosition();
	  MPIWrapperGrains::addToMPIString( oss.str() );
	}
        (*particule)->setGeoLocalisation( current_cell->m_GeoLocCell );
	particulesHalozone->push_back(*particule);
      }
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// La particule est elle en contact avec un obstacle periodique
bool LinkedCell::isContactVdW_ObstaclePeriodique( const Particule* particule )
	const
{
  Cellule* current_cell = getCellule( *particule->getPosition() );
  return ( current_cell->isContactVdW_ObstaclePeriodique( particule ) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise a jour des tags des particules de la zone de recouvrement et des clones
void LinkedCell::updateHalozoneCloneTag( Scalar time,
	list<Particule*>* particulesHalozone,
	list<Particule*>* particulesClones,
	MPIWrapperGrains const* wrapper )
{
  list<Particule*>::iterator particule;
  int tag;

  for (particule=particulesHalozone->begin();
  	particule!=particulesHalozone->end(); )
  {
    tag = (*particule)->setTag(
    	getCellule( *(*particule)->getPosition() )->m_tag );

    // Halozone to Clone (1 -> 2)
    if ( tag == 2 )
    {
      if ( Grains_Exec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT ) <<
      		" Halozone to Clone (1 -> 2)                  Id = " <<
      		(*particule)->getID() << " " << *(*particule)->getPosition();
        MPIWrapperGrains::addToMPIString( oss.str() );
      }
      particulesClones->push_back(*particule);
      (*particule)->setGeoLocalisation( MPIGEO_NONE );
      particule = particulesHalozone->erase( particule );
    }
    // Halozone to Interior (1 -> 0)
    else if (tag==0)
    {
      if ( Grains_Exec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT ) <<
      		" Halozone to Interior (1 -> 0)               Id = " <<
      		(*particule)->getID() << " " << *(*particule)->getPosition();
        MPIWrapperGrains::addToMPIString( oss.str() );
      }
      (*particule)->setGeoLocalisation( MPIGEO_NONE );
      particule = particulesHalozone->erase( particule );
    }
    else particule++;
  }

  for (particule=particulesClones->begin();
  	particule!=particulesClones->end(); )
  {
    tag = (*particule)->setTag(
    	getCellule( *(*particule)->getPosition() )->m_tag );

    // Clone to Halozone (2 -> 1)
    if ( tag == 1 )
    {
      if ( Grains_Exec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT ) <<
      		" Clone to Halozone (2 -> 1)                  Id = " <<
      		(*particule)->getID() << " " << *(*particule)->getPosition();
        MPIWrapperGrains::addToMPIString( oss.str() );
      }
      (*particule)->setGeoLocalisation(
      	getCellule( *(*particule)->getPosition() )->m_GeoLocCell );
      particulesHalozone->push_back(*particule);
      particule = particulesClones->erase( particule );
    }
    else particule++;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Suppression des clones sortis de la grille
void LinkedCell::DestroyOutOfDomainClones( Scalar time,
	list<Particule*>* particulesClones,
	set<Particule*>* particulesReferencesPeriodiques,
	list<Particule*>* ENSparticules,
	MPIWrapperGrains const* wrapper )
{
  list<Particule*>::iterator particule;
  list<App*>::iterator app;
  Particule *pdestroy=NULL;

  for (particule=particulesClones->begin(); particule!=particulesClones->end();)
  {
    if ( !isInLinkedCell( *(*particule)->getPosition() ) )
    {
      if ( Grains_Exec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT ) <<
      		" Destroy clone                               Id = " <<
      		(*particule)->getID() << " " << *(*particule)->getPosition();
        MPIWrapperGrains::addToMPIString( oss.str() );
      }
      pdestroy = *particule;

      // Suppression de la particule dans la derni�re cellule � laquelle
      // elle a appartenu avant de sortir du domaine
      pdestroy->getCelluleNm1()->remove( pdestroy );

      // Suppression des differentes listes
      removeParticuleFromList( *ENSparticules, pdestroy );
      if ( Grains_Exec::m_periodique == true )
        if ( pdestroy->getNombreClonesPeriodiques() )
	  removeParticuleFromSet( *particulesReferencesPeriodiques, pdestroy );

      // Destruction de l'objet point�
      delete pdestroy;

      // Suppression de la liste des clones
      particule = particulesClones->erase( particule );
    }
    else particule++;
  }

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Traitement des particules periodiques: creation et destruction des
// clones periodiques de ces particules
void LinkedCell::LinkUpdateParticulesPeriodiques( Scalar time,
  	list<Particule*>* ENSparticules,
  	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference )
{
  if ( Grains_Exec::m_periodique == true )
  {
    list<MonObstacle*>::iterator myObs;
    list<Particule*>::iterator particule;

    // Destruction des clones multi-periodiques
    for (particule=particulesClonesPeriodiques->begin();
    	particule!=particulesClonesPeriodiques->end(); )
      if ( (*particule)->getNbPeriodes() > 1 )
      {
        Particule *pdestroy=*particule;

        // Suppression de la particule dans le LinkedCell
	// deja realise au pas de temps precedent
	// par EnsComposant::updateClonesPeriodiques

        // Suppression dans la liste des clones de la particule de reference
        pdestroy->getPeriodicReference()->erasePeriodicClone( pdestroy );

        // Destruction de l'objet point�
        delete pdestroy;

        // Suppression de la liste particulesClonesPeriodiques
        particule = particulesClonesPeriodiques->erase( particule );
      }
      else particule++;

    // LinkUpdate des clones periodiques
    for (particule=particulesClonesPeriodiques->begin();
    	particule!=particulesClonesPeriodiques->end(); particule++)
      LinkUpdateActiveParticule( *particule );

    // Creation/destruction des clones uni-periodiques
    for (myObs=m_allObstacles.begin();myObs!=m_allObstacles.end();myObs++)
      (*myObs)->LinkUpdateParticulesPeriodiques( particulesClonesPeriodiques,
    	ParticuleClassesReference, this );

    // Creation des clones multi-periodiques
    for (particule=ENSparticules->begin(); particule!=ENSparticules->end();
	 particule++)
      (*particule)->createMultiPeriodicClones( particulesClonesPeriodiques,
      	ParticuleClassesReference, this );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Traitement des particules periodiques en parallele
void LinkedCell::LinkUpdateParticulesPeriodiques_MPI_Step1(
	Scalar time,
	list<int>& ClonestoDestroy,
  	list<int>& ClonestoParticules,
	list<int>& PartRefPerHalozone,
  	list<int>& PartRefPerOutDomainHalozone,
  	list<int>& InNotRefPerHalozone,
  	list<Particule*>* particulesClonesPeriodiques,
  	set<Particule*>* particulesReferencesPeriodiques,
	list<Particule*>* ENSparticules,
  	list<Particule*>* particulesHalozone )
{
  Particule *pdestroy = NULL;
  int old_tag,tag;

  if ( Grains_Exec::m_periodique == true )
  {
    list<MonObstacle*>::iterator myObs;
    list<Particule*>::iterator particule;

    // LinkUpdate des clones periodiques
    // Si le clone periodique est hors du domaine de calcul -> destruction
    for (particule=particulesClonesPeriodiques->begin();
    	particule!=particulesClonesPeriodiques->end(); )
      if ( !isInLinkedCell( *(*particule)->getPosition() ) )
      {
        if ( Grains_Exec::m_MPI_verbose )
        {
	  ostringstream oss;
          oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT ) <<
      		" Destroy periodic clone         Id = " <<
      		(*particule)->getPeriodicReferenceID() << " " <<
		*(*particule)->getPosition();
          MPIWrapperGrains::addToMPIString( oss.str() );
	}
        pdestroy = *particule;

        // Suppression de la particule dans la derni�re cellule � laquelle
        // elle a appartenu avant de sortir du domaine
        pdestroy->getCelluleNm1()->remove( pdestroy );

        // Destruction de l'objet point�
        delete pdestroy;

        // Suppression de la liste des clones periodiques
        particule = particulesClonesPeriodiques->erase( particule );
      }
      else
      {
	old_tag = (*particule)->getCelluleNm1()->m_tag;
	LinkUpdateActiveParticule( *particule );
	tag = (*particule)->setTag(
		getCellule( *(*particule)->getPosition() )->m_tag );

	if ( tag != old_tag )
	{
          if ( Grains_Exec::m_MPI_verbose )
          {
	    ostringstream oss;
            oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT ) <<
      		" Periodic Clone (" << old_tag << " -> " << tag <<
		")                     Id = " <<
      		(*particule)->getID() << " " <<
		(*particule)->getPeriodicReferenceID() << " " <<
		*(*particule)->getPosition();
            MPIWrapperGrains::addToMPIString( oss.str() );
	  }
	}
	particule++;
      }

    // Gestion des particules de reference des clones uni-periodiques
    for (myObs=m_allObstacles.begin();myObs!=m_allObstacles.end();myObs++)
      (*myObs)->LinkUpdateParticulesPeriodiques_MPI( time,
	ClonestoDestroy, ClonestoParticules,
  	PartRefPerHalozone, PartRefPerOutDomainHalozone,
  	InNotRefPerHalozone, particulesReferencesPeriodiques,
	ENSparticules, particulesHalozone, this );

    // Creation des clones multi-periodiques
    // TO DO
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Contient une paroi periodique
bool LinkedCell::intersectObstaclePeriodique() const
{
  bool b_intersect = false;

  if ( Grains_Exec::m_periodique == true )
  {
    BBox boite(
  	Point( m_xmin - 1.05 * m_arete_X,
		m_ymin - 1.05 * m_arete_Y,
		m_zmin - 1.05 * m_arete_Z ),
	Point( m_xmax + 1.05 * m_arete_X,
		m_ymax + 1.05 * m_arete_Y,
		m_zmax + 1.05 * m_arete_Z ) );

    // Creation/destruction des clones uni-periodiques
    list<MonObstacle*>::const_iterator myObs;
    for (myObs=m_allObstacles.begin();myObs!=m_allObstacles.end() &&
    	!b_intersect;myObs++)
      b_intersect = (*myObs)->obstaclePeridiqueIntersect( boite );
  }

  return b_intersect;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie la taille d'une cellule dans une direction
Scalar LinkedCell::getCelluleSize( int const& dir ) const
{
  Scalar size = 0.;

  switch(dir)
  {
    case 0: size = m_arete_X;
      break;
    case 1: size = m_arete_Y;
      break;
    case 2: size = m_arete_Z;
      break;
  }

  return size;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie un pointeur sur le vecteur de cellules
vector<Cellule*> const* LinkedCell::getAllCellules() const
{
  return &m_cellules;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajout des nouvelles references periodiques
void LinkedCell::addNewPeriodicReference_MPI( Scalar time,
	set<Particule*>* particulesReferencesPeriodiques )
{
  list<MonObstacle*>::const_iterator myObs;
  for (myObs=m_allObstacles.begin();myObs!=m_allObstacles.end();myObs++)
    (*myObs)->addNewPeriodicReference_MPI( time,
    	particulesReferencesPeriodiques );
}
