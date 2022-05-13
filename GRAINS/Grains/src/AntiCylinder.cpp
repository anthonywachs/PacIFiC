#include "AntiCylinder.H"

/**
  @brief The Class AntiCylinder.
  @details An antiCylinder is a cubic obstacle with a cylindrical hole inside
  @author M.BERNARD - April 2015 - Creation
*/

int AntiCylinder::visuNodeNbOnPer = 12;

// ----------------------------------------------------------------------------
AntiCylinder::AntiCylinder(Scalar r , Scalar h ) : 
    m_radius(r), m_halfHeight(h / 2) 
{
}




// ----------------------------------------------------------------------------
// Construction avec un fichier
AntiCylinder::AntiCylinder(istream& fileIn) 
{
  readClass(fileIn);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructor
AntiCylinder::AntiCylinder(DOMNode* root)
{
  m_radius     = ReaderXML::getNodeAttr_Double(root, "Radius" );
  m_halfHeight = ReaderXML::getNodeAttr_Double(root, "Hauteur") / 2.;
}




// ----------------------------------------------------------------------------
AntiCylinder::~AntiCylinder() 
{
}




// ----------------------------------------------------------------------------
// Determine l'inertie et l'inertie inverse d'un cylindre
bool AntiCylinder::BuildInertie( Scalar *inertie,Scalar *inertie_1 ) const 
{
  return true;
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int AntiCylinder:: getShrinkingChoice()const
{
  return(0);
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// fixe le rayon courant 
void AntiCylinder::setShrinkingRadius(Scalar CurrentRadius)
{

}




// ----------------------------------------------------------------------------
//  Determine le rayon circonscrit d'un cylindre
Scalar AntiCylinder::BuildRayonRef() const 
{
  return sqrt(m_radius*m_radius+m_halfHeight*m_halfHeight);
}




// ----------------------------------------------------------------------------
// Construction d'un clone du cylindre
Convex* AntiCylinder::clone() const
{
  return new AntiCylinder(m_radius,2*m_halfHeight);
}




// ----------------------------------------------------------------------------
// Determine le Volume d'un cylindre
Scalar AntiCylinder::getVolume() const
{
  return 2*m_halfHeight*PI*m_radius*m_radius;
}




// ----------------------------------------------------------------------------
Point AntiCylinder::support(const Vecteur& v) const 
{
  cout << "  WARNING, should not go throw ";
  cout << "AntiCylinder::support"<<endl;
  Scalar norm = Norm(v);
  if (norm > EPSILON) {
    Scalar s = sqrt(v[X] * v[X] + v[Z] * v[Z]);
    if (s > EPSILON) {
      Scalar d = m_radius / s;  
      return Point(v[X] * d, v[Y] < 0. ? -m_halfHeight : m_halfHeight, v[Z] * d);
    } else {
      return Point(0., v[Y] < 0. ? -m_halfHeight : m_halfHeight, 0.);
    }
  } else {
    return Point();
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant le Convex pour visu dans une appli externe (Fluide)
// Pour un disque pas de prise en compte de points sommets -> 0
vector<Point> AntiCylinder::getEnveloppe() const
{
  Point point(0.,0.,0.);
  vector<Point> enveloppe(3,point);
  enveloppe[0][Y] = - m_halfHeight;
  enveloppe[1][Y] = - m_halfHeight;
  enveloppe[1][X] = m_radius;  
  enveloppe[2][Y] = m_halfHeight;    
  return enveloppe;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code équivalent
int AntiCylinder::getNbCorners() const
{
  return 1000;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Description des faces
vector< vector<int> > const* AntiCylinder::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return allFaces;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un Cylindre
void AntiCylinder::printClass(ostream &fileOut) const 
{
  fileOut << "*Cylindre\n";
  fileOut << m_radius << '\t' 
          << 2.0*m_halfHeight << '\n';
}



// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture d'un Cylindre
void AntiCylinder::readClass(istream &fileIn) 
{
  fileIn >> m_radius 
         >> m_halfHeight;
  m_halfHeight /= 2.0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de points pour post-processing avec Paraview
int AntiCylinder::numberOfPoints_PARAVIEW() const
{
  return 2*visuNodeNbOnPer+2;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int AntiCylinder::numberOfCells_PARAVIEW() const
{
  return visuNodeNbOnPer;  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview  
void AntiCylinder::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  cout << "  WARNING, should not go throw ";
  cout << "AntiCylinder::write_polygonsPts_PARAVIEW"<<endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview  
list<Point> AntiCylinder::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation ) const
{
  cout << "  WARNING, should not go throw ";
  cout << "AntiCylinder::get_polygonsPts_PARAVIEW"<<endl;
  list<Point> ParaviewPoints;
  return ParaviewPoints; 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview 
void AntiCylinder::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  cout << "  WARNING, should not go throw ";
  cout << "AntiCylinder::write_polygonsStr_PARAVIEW"<<endl;
}
