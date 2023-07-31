// D.PETIT - Juil 2000 - Creation
// ============================================================================
#include "Polytope.H"
#include <new>


// ----------------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
// Constructeur avec fichier
Polytope::Polytope( istream &fileIn, int nb_point, VertexBase &ref,
	IndexArray& ia ) : 
  m_base( ref ),
  m_index( ia ) 
{ 
  // Lecture des coordonnees des points
  for(int i=0; i<nb_point; i++) fileIn >> ref[i];
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
// Constructeur par copie
Polytope::Polytope( const Polytope& copie ) : 
  Convex( copie ),
  m_base( copie.m_base ), 
  m_index( copie.m_index )
{
  m_fichPoly = copie.m_fichPoly ; 
}


  

// ----------------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
Polytope::~Polytope()
{
  // m_base et m_index sont detruits par le mecanisme de garbage collector 
  // de Grains_Exec 
}




// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
const Point& Polytope::operator[](int i) const 
{ 
  return m_base[m_index[i]]; 
}




// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
Point& Polytope::operator[](int i) 
{
  return m_base[m_index[i]];
}




// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
int Polytope::numVerts() const 
{ 
  return m_index.size();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant l'enveloppe
vector<Point> Polytope::getEnveloppe() const
{
  vector<Point> enveloppe;
  for (int i=0; i<numVerts(); i++) {
    const Point& p = (*this)[i];
    enveloppe.push_back(p);
  }
  return enveloppe;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code équivalent
int Polytope::getNbCorners() const
{
  return numVerts();
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview 
int Polytope::numberOfPoints_PARAVIEW() const 
{ 
  // centre de gravite + nb de vertex
  return ( 1 + m_index.size() );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview  
void Polytope::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  Point pp, oo;
  
  // Gravity center
  pp = transform(oo);
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  // Corners
  for (int i=0; i<numVerts(); i++) 
  {
    pp = transform((*this)[i]);
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;    
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview  
list<Point> Polytope::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation ) const
{
  list<Point> ParaviewPoints;
  Point pp, oo;
  
  // Gravity center
  pp = transform(oo);
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back(pp);  
  
  // Corners  
  for (int i=0;i<numVerts();++i)
  {
    pp = transform((*this)[i]);
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back(pp);
  }
  
  return ParaviewPoints; 
}
