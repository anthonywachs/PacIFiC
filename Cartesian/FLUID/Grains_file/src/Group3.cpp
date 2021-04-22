#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "Group3.H"


namespace solid
{
  // --------------------------------------------------------------------------
  // Constructeur par defaut
  Group3::Group3( const Scalar def )
  {
    comp[X] = comp[Y] = comp[Z] = def;
  }




  // --------------------------------------------------------------------------
  // Constructeur avec initialisation.
  Group3::Group3( const Scalar x, const Scalar y, const Scalar z ) 
  {
    comp[X] = x;
    comp[Y] = y;
    comp[Z] = z;
  }




  // --------------------------------------------------------------------------
  // Constructeur par copie.
  Group3::Group3( const Group3& g )
  {
    comp[X] = g.comp[X];
    comp[Y] = g.comp[Y];
    comp[Z] = g.comp[Z];
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // Destructeur.
  Group3::~Group3() 
  {
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // Renvoi l'element.
  const Scalar* Group3::getValue() const 
  {
    return comp;
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // Renvoi l'element.
  Scalar* Group3::getValue() 
  {
    return comp;
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // Modifie les composantes d'un element. 
  void Group3::setValue( const Scalar x, const Scalar y, const Scalar z )
  {
    comp[X] = x;
    comp[Y] = y;
    comp[Z] = z;
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // Modifie les composantes d'un element. 
  void Group3::setValue( const Scalar* g )
  {
    comp[X] = g[X];
    comp[Y] = g[Y];
    comp[Z] = g[Z];
  }




  // --------------------------------------------------------------------------
  // A.WACHS - Oct. 2009 - Creation
  // Annule toutes les composantes d'un element. 
  void Group3::reset()
  {
    comp[X] = comp[Y] = comp[Z] = 0.;
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // nombre de composantes d'un element. 
  int Group3::size() const 
  {
    return 3;
  }



 
  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // surcharge de l'operateur d'indexation (pour affecter).
  Scalar& Group3::operator[]( const int i )
  {
    return comp[i];
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // surcharge de l'operateur d'indexation (pour affecter).
  const Scalar& Group3::operator[]( const int i ) const 
  {
    return comp[i];
  }




  // --------------------------------------------------------------------------
  // Operateur unaire - 
  Group3 Group3::operator - () const
  {
    return Group3(-comp[X], -comp[Y], -comp[Z]);
  }




  // --------------------------------------------------------------------------
  // Produit scalaire de deux objets de type Group3.
  Scalar Group3::operator * ( const Group3 &g ) const
  {
    return (comp[X]*g.comp[X] + comp[Y]*g.comp[Y] + comp[Z]*g.comp[Z]);
  }




  // --------------------------------------------------------------------------
  // Multiplication par un scalaire
  Group3 Group3::operator * ( Scalar d ) const
  {
    return Group3(comp[X]*d, comp[Y]*d, comp[Z]*d);
  }




  // --------------------------------------------------------------------------
  // Division du Group3 par un scalaire.
  Group3 Group3::operator / ( Scalar d ) const
  {
    return Group3(comp[X]/d, comp[Y]/d, comp[Z]/d);
  }




  // --------------------------------------------------------------------------
  // Addition de deux Group3 composante par composante
  Group3 Group3::operator + ( const Group3 &g2 ) const
  {
    return Group3(comp[X]+g2.comp[X], comp[Y]+g2.comp[Y], comp[Z]+g2.comp[Z]);
  }




  // --------------------------------------------------------------------------
  // Soustraction de deux Group3 composante par composante
  Group3 Group3::operator - ( const Group3 &g2 ) const
  {
    return Group3(comp[X]-g2.comp[X], comp[Y]-g2.comp[Y], comp[Z]-g2.comp[Z]);
  }




  // --------------------------------------------------------------------------
  // Comparaison pour les valeurs
  bool Group3::operator == ( const Group3 &g2 ) const
  {
    return (comp[X]==g2[X] && comp[Y]==g2[Y] && comp[Z]==g2[Z]);
  }




  // --------------------------------------------------------------------------
  // Objets differents par les valeurs */
  bool Group3::operator != ( const Group3 &g2 )
  {
    return !(*this==g2);
  }




  // --------------------------------------------------------------------------
  // Egalite d'affectation d'un vecteur 
  // F.PRADEL - Janvier 2000 - Creation
  Group3& Group3::operator = ( const Group3 &g2 )
  {
    if ( &g2 != this )
    {      
      comp[X] = g2.comp[X];
      comp[Y] = g2.comp[Y];
      comp[Z] = g2.comp[Z];
    }
    return (*this);
  }




  // --------------------------------------------------------------------------
  // Affectation d'une valeur aux composantes du vecteur.
  // G.FERRER - Mars.2000 - Creation
  void Group3::operator = ( Scalar valeur )
  {
    comp[X] = comp[Y] = comp[Z] = valeur;
  }




  // --------------------------------------------------------------------------
  // Multiplication par un scalaire et d'affectation
  Group3& Group3::operator *= ( Scalar d )
  {
    comp[X] *= d;
    comp[Y] *= d;
    comp[Z] *= d;
    return *this;
  }




  // --------------------------------------------------------------------------
  // Division du vecteur courant par un scalaire.
  Group3& Group3::operator /= ( Scalar d )
  {
    comp[X] /= d;
    comp[Y] /= d;
    comp[Z] /= d;
    return *this;
  }




  // --------------------------------------------------------------------------
  // Addition d'un vecteur et d'affectation*/
  Group3& Group3::operator += ( const Group3 &g2 )
  {
    comp[X] += g2.comp[X];
    comp[Y] += g2.comp[Y];
    comp[Z] += g2.comp[Z];
    return *this;  
  }




  // --------------------------------------------------------------------------
  // Soustraction d'un vecteur et d'affectation
  Group3& Group3::operator -= ( const Group3 &g2 )
  {
    comp[X] -= g2.comp[X];
    comp[Y] -= g2.comp[Y];
    comp[Z] -= g2.comp[Z];
    return *this;  
  }




  // --------------------------------------------------------------------------
  // D.PETIT - Juin. 2000 - Creation
  // EfFectue le produit mixte de trois vecteurs
  Scalar triple( const Group3& g1, const Group3& g2, const Group3& g3 )
  {
    return g1.comp[X] * (g2.comp[Y] * g3.comp[Z] - g2.comp[Z] * g3.comp[Y]) +
      g1.comp[Y] * (g2.comp[Z] * g3.comp[X] - g2.comp[X] * g3.comp[Z]) +
      g1.comp[Z] * (g2.comp[X] * g3.comp[Y] - g2.comp[Y] * g3.comp[X]);
  }




  // --------------------------------------------------------------------------
  // Multiplication d'un scalaire par un Group3
  Group3 operator * ( Scalar d, const Group3 &g )
  {
    Group3 result(g);
    result *= d;
    return result;
  }




  // --------------------------------------------------------------------------
  // Operateur d'ecriture.
  ostream &operator << ( ostream &fileOut, const Group3 &g ) 
  {
    fileOut << g[X] << '\t' 
	    << g[Y] << '\t'
	    << g[Z] << '\n';
    return (fileOut);
  }




  // --------------------------------------------------------------------------
  // Operateur de lecture.
  istream &operator >> ( istream &fileIn,Group3 &g ) 
  {
    fileIn >> g[X] 
	   >> g[Y] 
	   >> g[Z];
    return (fileIn);
  }
  


  
  // --------------------------------------------------------------------------
  // Operateur d'ecriture
  // A.WACHS - Janv.2011 - Creation.
  void Group3::writeGroup3( ostream &fileOut ) const
  {
    fileOut << Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	comp[X]) << " " << 
	Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	comp[Y]) << " " << 
	Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	comp[Z]);
  }
  
  

  
  // --------------------------------------------------------------------------
  // Operateur d'ecriture en binaire
  // A.WACHS - Dec.2014 - Creation.
  void Group3::writeGroup3_binary( ostream &fileOut )
  {
    fileOut.write( reinterpret_cast<char*>( &comp[X] ), sizeof(double) );
    fileOut.write( reinterpret_cast<char*>( &comp[Y] ), sizeof(double) );  
    fileOut.write( reinterpret_cast<char*>( &comp[Z] ), sizeof(double) );     
  }
  


  
  // --------------------------------------------------------------------------
  // Operateur de lecture en binaire
  // A.WACHS - Dec.2014 - Creation.
  void Group3::readGroup3_binary( istream &StreamIN )
  {
    StreamIN.read( reinterpret_cast<char*>( &comp[X] ), sizeof(double) );
    StreamIN.read( reinterpret_cast<char*>( &comp[Y] ), sizeof(double) );  
    StreamIN.read( reinterpret_cast<char*>( &comp[Z] ), sizeof(double) );     
  }        
}
