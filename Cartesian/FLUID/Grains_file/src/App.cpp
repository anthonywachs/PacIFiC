#include "App.H"
#include "EnsComposant.H"

Scalar App::m_dX = 0.;
Scalar App::m_dY = 0.;
Scalar App::m_dZ = 0.;
Scalar App::m_lX = 0.;
Scalar App::m_lY = 0.;
Scalar App::m_lZ = 0.;
Point App::m_origine_locale;
Point App::m_origine_globale;

// ----------------------------------------------------------------------------
// Constructeur
// G.FERRER - Mai .2003 - Creation
App::App()
{
}




// ----------------------------------------------------------------------------
// Destructeur
// G.FERRER - Mai .2003 - Creation
App::~App()
{}




// ----------------------------------------------------------------------------
// Accesseur
// G.FERRER - Mai .2003 - Creation
void App::setD( Scalar xmax, Scalar ymax, Scalar zmax,
  	Scalar ox, Scalar oy, Scalar oz )
{
 App::m_origine_globale[0] = ox;
 App::m_origine_globale[1] = oy;
 App::m_origine_globale[2] = oz;
 App::m_dX = xmax - m_origine_globale[0];
 App::m_dY = ymax - m_origine_globale[1];
 App::m_dZ = zmax - m_origine_globale[2];
 App::m_lX = xmax - m_origine_globale[0];
 App::m_lY = ymax - m_origine_globale[1];
 App::m_lZ = zmax - m_origine_globale[2];
 assert( App::m_dX >= 0. && App::m_dY >= 0. && App::m_dZ >= 0. );   
}




// ----------------------------------------------------------------------------
// Accesseur
// A. WACHS - Mai .2009 - Creation
void App::setDlocale( Scalar lx_, Scalar ly_, Scalar lz_ )
{
 App::m_lX = lx_;
 App::m_lY = ly_;
 App::m_lZ = lz_;
}




// ----------------------------------------------------------------------------
// Accesseur
// G.FERRER - Mai .2009 - Creation
void App::setOriginelocale( int const* nprocsdir, int const* MPIcoords ) 
{
  m_origine_locale[0] = m_origine_globale[0] + MPIcoords[0] * m_lX ;
  m_origine_locale[1] = m_origine_globale[1] + MPIcoords[1] * m_lY ;  
  m_origine_locale[2] = m_origine_globale[2] + MPIcoords[2] * m_lZ ;  
}




// ----------------------------------------------------------------------------
// Accesseur
// M.BERNARD - Mai 2012 - Creation
void App::getOrigineLocale( double& x, double& y, double& z ) 
{
  x = m_origine_locale[0] ;
  y = m_origine_locale[1] ;  
  z = m_origine_locale[2] ;
}




// ----------------------------------------------------------------------------
// Accesseur
// D. RAKOTONIRINA - Avril 2017 - Creation
void App::getOrigineGlobale( double& x, double& y, double& z ) 
{
  x = m_origine_globale[0] ;
  y = m_origine_globale[1] ;  
  z = m_origine_globale[2] ;
}




// ----------------------------------------------------------------------------
// Accesseur
// M.BERNARD - Fevrier 2015 - Creation
void App::getDimensionsLocales( double& lx, double& ly, double& lz ) 
{
  lx = m_lX ;
  ly = m_lY ;  
  lz = m_lZ ;
}




// ----------------------------------------------------------------------------
// Accesseur
// D. RAKOTONIRINA - Avril 2017 - Creation
void App::getDimesionsGlobales( double& lx, double& ly, double& lz ) 
{
  lx = m_dX ;
  ly = m_dY ;  
  lz = m_dZ ;
}




// ----------------------------------------------------------------------------
// Affiche les attributs statiques
// A.WACHS Aout 2005
void App::affiche_attributs_statiques( ostream &f )
{
  f << "Domain dimension = " << m_dX << " x " << m_dY << " x " << m_dZ << endl;
  f << "Local domain dimension = " << m_lX << " x " << m_lY << " x " 
  	<< m_lZ << endl;
  f << "Local origin = " << m_origine_locale;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Appartenance du centre de gravite d'une particule au domaine de calcul
bool App::isInDomain( Point const* position )
{
  bool isIn = true;
  
  if ( (*position)[0] < m_origine_globale[0] 
  	|| (*position)[0] > m_origine_globale[0] + m_dX 
  	|| (*position)[1] < m_origine_globale[1] 
	|| (*position)[1] > m_origine_globale[1] + m_dY
  	|| (*position)[2] < m_origine_globale[2] 
	|| (*position)[2] > m_origine_globale[2] + m_dZ ) 
    isIn = false;
  
  return isIn;
}  




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Appartenance du centre de gravite d'une particule au domaine de calcul
bool App::isInLocalDomain( Point const* position )
{
  bool isIn = true;
  
  if ( (*position)[0] < m_origine_locale[0] 
  	|| (*position)[0] > m_origine_locale[0] + m_lX 
  	|| (*position)[1] < m_origine_locale[1] 
	|| (*position)[1] > m_origine_locale[1] + m_lY
  	|| (*position)[2] < m_origine_locale[2] 
	|| (*position)[2] > m_origine_locale[2] + m_lZ ) 
    isIn = false;
  
  return isIn;
} 




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Tells if an app corresponding to the argument exists or not
bool App::isName( const string name_ )
{
  bool result = false;
  if( m_name == name_ )
    result = true;
    
  return result;
}
