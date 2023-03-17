/*
  GJK Engine - A Fast and Robust GJK Implementation
  Copyright (C) 1998  Gino van den Bergen

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Library General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Library General Public License for more details.

  You should have received a copy of the GNU Library General Public
  License along with this library; if not, write to the Free
  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  Please send remarks, questions and bug reports to gino@win.tue.nl,
  or write to:
                  Gino van den Bergen
		  Department of Mathematics and Computing Science
		  Eindhoven University of Technology
		  P.O. Box 513, 5600 MB Eindhoven, The Netherlands
*/
// ============================================================================
#include "Grains_Exec.hh"
#include "Polyhedron.H"
#include "EnsComposant.H"
#include "Plan.H"
#include <fstream>
#include <new>
#include <sstream>
#include <vector>
using namespace std;


typedef vector<unsigned int> IndexBuf;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// D.PETIT - Juin. 2000 - Creation
// Construction avec un fichier
Polyhedron::Polyhedron( istream &fileIn, int nb_point, VertexBase& ref,
  	 IndexArray& ia ) :
  Polytope( fileIn, nb_point, ref, ia ),
  m_cobound( NULL ),
  m_curr_vertex( 0 ),
  m_InertiePoly( NULL ),
  m_VolumePoly( 0. )
{
  readface(fileIn);
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Constructeur de copie
Polyhedron::Polyhedron( const Polyhedron& copie ) :
  Polytope( copie )
{
  m_cobound = copie.m_cobound;
  m_curr_vertex = copie.m_curr_vertex;
  m_InertiePoly = NULL;
  m_VolumePoly = copie.m_VolumePoly;
  m_allFaces = copie.m_allFaces;
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Destructeur
Polyhedron::~Polyhedron()
{
  // Le tenseur d'inertie n'est construit que pour les particules de reference
  // de classe, m_InertiePoly est NULL si la particule est active
  if ( m_InertiePoly ) delete m_InertiePoly;

  // m_cobound est detruit par le mecanisme de garbage collector de Grains_Exec
}




// ----------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// creation d'un Polyhedron
Polyhedron* Polyhedron::create( istream &fileIn )
{
  // Lecture du nom puis ouverture du fichier contenant la description
  // du polyedre
  string fichpoly;
  fileIn >> fichpoly;
  fichpoly = Grains_Exec::m_ReloadDirectory + "/"
  	+ Grains_Exec::extractFileName( fichpoly );
  ifstream *PolyIN = new ifstream( fichpoly.c_str() );
  if ( PolyIN->is_open() )
    Grains_Exec::m_additionalDataFiles.insert( fichpoly );

  // Lecture du nb de points definissant le polyedre
  int nb_point;
  *PolyIN >> nb_point;
  *PolyIN >> nb_point;

  // Creation du tableau de points du polyedre
  Point *point = new Point[nb_point];
  VertexBase *vertexbase = new VertexBase( (void *)point );

  // Creation du tableau d'indices des sommets dans le tableau de points
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v );
  delete [] v;

  // Creation du polyedre
  Polyhedron *polyhedron = new Polyhedron( *PolyIN, nb_point, *vertexbase,
  	*ia );
  polyhedron->m_fichPoly = fichpoly;
  PolyIN->close();

  // Les objets m_base, m_cobound, m_index et m_allFaces ne sont crees qu'une
  // seule fois par type de particule, les particules suivantes ne possedent
  // qu'un pointeur sur ces derniers.
  // Un mecanisme de garbage collector a ete mis en place dans la classe
  // Grains_Exec qui conserve les pointeurs et se charge de les detruire
  Grains_Exec::addOnePolytopeRefPointBase( point, vertexbase );
  Grains_Exec::addOnePolytopeNodeNeighbors( polyhedron->m_cobound );
  Grains_Exec::addOnePolytopeNodeIndex( ia );
  Grains_Exec::addOnePolyhedronFaceConnectivity( polyhedron->m_allFaces );

  delete PolyIN;

  return polyhedron;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction d'un Polyhedron
Polyhedron* Polyhedron::create( DOMNode* root )
{
  // Lecture du nom puis ouverture du fichier contenant la description
  // du polyedre
  string fichpoly = ReaderXML::getNodeAttr_String( root, "Name" );
  ifstream PolyIN( fichpoly.c_str() );
  if ( PolyIN.is_open() )
    Grains_Exec::m_additionalDataFiles.insert( fichpoly );

  // Lecture du nb de points definissant le polyedre
  int nb_point;
  PolyIN >> nb_point >> nb_point;

  // Creation du tableau de points du polyedre
  Point *point = new Point[nb_point];
  VertexBase *vertexbase = new VertexBase( (void *)point );

  // Creation du tableau d'indices des sommets dans le tableau de points
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v );
  delete [] v;

  // Creation du polyedre
  Polyhedron *polyhedron = new Polyhedron( PolyIN, nb_point, *vertexbase,
  	*ia );
  polyhedron->m_fichPoly = fichpoly;

  // Les objets m_base, m_cobound, m_index et m_allFaces ne sont crees qu'une
  // seule fois par type de particule, les particules suivantes ne possedent
  // qu'un pointeur sur ces derniers.
  // Un mecanisme de garbage collector a ete mis en place dans la classe
  // Grains_Exec qui conserve les pointeurs et se charge de les detruire
  Grains_Exec::addOnePolytopeRefPointBase( point, vertexbase );
  Grains_Exec::addOnePolytopeNodeNeighbors( polyhedron->m_cobound );
  Grains_Exec::addOnePolytopeNodeIndex( ia );
  Grains_Exec::addOnePolyhedronFaceConnectivity( polyhedron->m_allFaces );

  return polyhedron;
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un Polyhedron
bool Polyhedron::BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const
{
  std::copy(&m_InertiePoly[0], &m_InertiePoly[6], &inertie[0]);

  Scalar determinant = inertie[0]*inertie[3]*inertie[5]
    - inertie[0]*inertie[4]*inertie[4]
    - inertie[5]*inertie[1]*inertie[1]
    - inertie[3]*inertie[2]*inertie[2]
    + 2*inertie[1]*inertie[2]*inertie[4];

  inertie_1[0] = (inertie[3]*inertie[5]-inertie[4]*inertie[4])/determinant;
  inertie_1[1] = (inertie[2]*inertie[4]-inertie[1]*inertie[5])/determinant;
  inertie_1[2] = (inertie[1]*inertie[4]-inertie[2]*inertie[3])/determinant;
  inertie_1[3] = (inertie[0]*inertie[5]-inertie[2]*inertie[2])/determinant;
  inertie_1[4] = (inertie[1]*inertie[2]-inertie[0]*inertie[4])/determinant;
  inertie_1[5] = (inertie[0]*inertie[3]-inertie[1]*inertie[1])/determinant;

  return true;
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Scalar Polyhedron::BuildRayonRef() const
{
  Scalar d , ray = ((*this)[0]) * ((*this)[0]);
  for (int i = 1; i< numVerts(); i++) {
    if ((d =  (*this)[i] * (*this)[i]) > ray)
      ray = d;
  }
  return sqrt(ray);
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int Polyhedron::getShrinkingChoice()const
{
  return(0);
}




// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// Set the shrinking radius value
void Polyhedron:: setShrinkingRadius(Scalar CurrentRadius)
{
}




// -------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
// calcul de l'inertie et du volume d'un tetraedre
// void Polyhedron::CalculVolumeInertie( const Point &H,
// 	const Point &P1, const Point &P2 )
// {
//   // On effectue d'abord un changement de repere
//   // Construction de la matrice de rotation
//   Vecteur x1(1,0,0), y1(0,1,0), z1(0,0,1);
//   Vecteur x2(H), y2(P1 - H);
//   Vecteur z2 = x2^y2;
//
//   Matrix Rot(cos(x1,x2), cos(x1,y2), cos(x1,z2),
// 	     cos(y1,x2), cos(y1,y2), cos(y1,z2),
// 	     cos(z1,x2), cos(z1,y2), cos(z1,z2));
//
//   Matrix tRot = Rot.transpose();
//
//   // on se ramene a la base relative
//   Scalar a = tRot[X] * x2;
//   Scalar b = tRot[Y][X] * P2[X] + tRot[Y][Y] * P2[Y] + tRot[Y][Z] * P2[Z];
//   Scalar c = tRot[Z][X] * P2[X] + tRot[Z][Y] * P2[Y] + tRot[Z][Z] * P2[Z];
//   Scalar d = tRot[Y] * y2;
//
//   Scalar V = a * c * d / 6.0;
//   m_VolumePoly += V;
//   cout << "V = " << V << endl;
//   Scalar inertie[11];
//   Scalar acd_60 = V/10.0;
//   inertie[0] = acd_60*(b*b + b*d + c*c + d*d);
//   inertie[1] = inertie[4] = - acd_60*2*a*(b+d);
//   inertie[2] = inertie[8] = - acd_60*2*c*a;
//   inertie[3] = 0.0;
//   inertie[5] = acd_60*(c*c+6*a*a);
//   inertie[6] = inertie[9] = - acd_60*0.5*c*(2*b+d);
//   inertie[7] = 0.0;
//   inertie[10] = acd_60*(b*b + b*d + 6*a*a + d*d);
//   Matrix I(inertie);
//
//   // changement de base (relative a l'absolue)
//   I = Rot * (I * tRot);
//   //I += GInertie;
//
//   m_InertiePoly[0] += (I.getValue())[0][0];
//   m_InertiePoly[1] += (I.getValue())[0][1];
//   m_InertiePoly[2] += (I.getValue())[0][2];
//   m_InertiePoly[3] += (I.getValue())[1][1];
//   m_InertiePoly[4] += (I.getValue())[1][2];
//   m_InertiePoly[5] += (I.getValue())[2][2];
// }
void Polyhedron::CalculVolumeInertie( const Point &A2,
	const Point &A3, const Point &A4 )
{
  // From Journal of Mathematics and Statistics 1 (1): 8-11, 2004
  // "Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor
  // in Terms of its Vertex Coordinates", F. Tonon

  double x1 = 0., x2 = A2[X], x3 = A3[X], x4 = A4[X],
  	y1 = 0., y2 = A2[Y], y3 = A3[Y], y4 = A4[Y],
	z1 = 0., z2 = A2[Z], z3 = A3[Z], z4 = A4[Z],
	det ;

  det = fabs( ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 )
  	+ ( y2 - y1 ) * ( z3 - z1 ) * (	x4 - x1 )
	+ ( z2 - z1 ) * ( x3 - x1 ) * (	y4 - y1 )
	- ( z2 - z1 ) * ( y3 - y1 ) * (	x4 - x1 )
	- ( x2 - x1 ) * ( z3 - z1 ) * (	y4 - y1 )
	- ( y2 - y1 ) * ( x3 - x1 ) * (	z4 - z1 ) );

  m_VolumePoly += det / 6. ;

  m_InertiePoly[0] += det * ( y1 * y1 + y1 * y2 + y2 * y2
  	+ y1 * y3 + y2 * y3 + y3 * y3
	+ y1 * y4 + y2 * y4 + y3 * y4 + y4 * y4
	+ z1 * z1 + z1 * z2 + z2 * z2
  	+ z1 * z3 + z2 * z3 + z3 * z3
	+ z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4 ) / 60. ;
  m_InertiePoly[1] -= det * ( 2. * x1 * z1 + x2 * z1 + x3 * z1 + x4 * z1
  	+ x1 * z2 + 2. * x2 * z2 + x3 * z2 + x4 * z2
	+ x1 * z3 + x2 * z3 + 2. * x3 * z3 + x4 * z3
	+ x1 * z4 + x2 * z4 + x3 * z4 + 2. * x4 * z4 ) / 120. ;
  m_InertiePoly[2] -= det * ( 2. * x1 * y1 + x2 * y1 + x3 * y1 + x4 * y1
  	+ x1 * y2 + 2. * x2 * y2 + x3 * y2 + x4 * y2
	+ x1 * y3 + x2 * y3 + 2. * x3 * y3 + x4 * y3
	+ x1 * y4 + x2 * y4 + x3 * y4 + 2. * x4 * y4 ) / 120. ;
  m_InertiePoly[3] += det * ( x1 * x1 + x1 * x2 + x2 * x2
  	+ x1 * x3 + x2 * x3 + x3 * x3
	+ x1 * x4 + x2 * x4 + x3 * x4 + x4 * x4
	+ z1 * z1 + z1 * z2 + z2 * z2
  	+ z1 * z3 + z2 * z3 + z3 * z3
	+ z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4 ) / 60. ;
  m_InertiePoly[4] -= det * ( 2. * y1 * z1 + y2 * z1 + y3 * z1 + y4 * z1
  	+ y1 * z2 + 2. * y2 * z2 + y3 * z2 + y4 * z2
	+ y1 * z3 + y2 * z3 + 2. * y3 * z3 + y4 * z3
	+ y1 * z4 + y2 * z4 + y3 * z4 + 2. * y4 * z4 ) / 120. ;
  m_InertiePoly[5] += det * ( x1 * x1 + x1 * x2 + x2 * x2
  	+ x1 * x3 + x2 * x3 + x3 * x3
	+ x1 * x4 + x2 * x4 + x3 * x4 + x4 * x4
	+ y1 * y1 + y1 * y2 + y2 * y2
  	+ y1 * y3 + y2 * y3 + y3 * y3
	+ y1 * y4 + y2 * y4 + y3 * y4 + y4 * y4 ) / 60. ;
}



// -------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
// initialisation de l'inertie et du volume d'un tetraedre
void Polyhedron::Initialisation()
{
  m_VolumePoly  = 0.0;
  m_InertiePoly = new Scalar[6];
  m_InertiePoly[0] = m_InertiePoly[1] = m_InertiePoly[2]
    = m_InertiePoly[3] = m_InertiePoly[4] = m_InertiePoly[5] = 0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Description des faces
vector<vector<int> > const* Polyhedron::getFaces() const
{
  return m_allFaces;
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un Polyhedron
Scalar Polyhedron::getVolume() const
{
  return m_VolumePoly;
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone d'un Polyhedron
Convex* Polyhedron::clone() const
{
  return new Polyhedron(*this);
}




// ----------------------------------------------------------------------
// Fonction support
// D.PETIT - Juil.2000 - Modif
// G.FERRER - Avri.2003 - Modif pour rendre accessible "support 2.0"
// A. WACHS - Aout 2009
// La condition "d <= h" a �t� remplac�e par la condition suivante "d-h < eps"
// pour des prob de precision machine lorsque d = h et evite que le code soit
// machine-dependant; en revanche cela introduit une dependance a eps
// Pour l 'instant, eps = 1.e-16
Point Polyhedron::support( const Vecteur& v ) const
{
  Scalar norm = Norm(v);
  if (norm > EPSILON) {

    if (m_cobound == NULL) {
      return support20(v);
    }

    int last_vertex = -1;
    Scalar h = (*this)[m_curr_vertex] * v, d = 0.;
    for (;;) {
      IndexArray& curr_cobound = m_cobound[m_curr_vertex];
      int i = 0, n = curr_cobound.size();
      while (i != n &&
	     (curr_cobound[i] == last_vertex
//	      || (d = (*this)[curr_cobound[i]] * v) <= h))
	      || (d = (*this)[curr_cobound[i]] * v) - h < 1.e-16 ))
	++i;
      if (i == n) break;
      last_vertex = m_curr_vertex;
      m_curr_vertex = curr_cobound[i];
      h = d;
    }
    return (*this)[m_curr_vertex];
  } else {
    return Point();
  }
}




// ----------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture du polyedre
void Polyhedron::printClass( ostream &fileOut ) const
{
  fileOut << "*Polyhedron\n";
  fileOut << Grains_Exec::extractFileName( m_fichPoly ) << endl;
}




// ----------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture de la sphere
void Polyhedron::readClass( istream &fileIn )
{
  cerr << "Program Error :\n"
       << "Polyhedron::readClass non accessible.\n";
  exit(3);
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction du polyhedron
// void Polyhedron::BuildPolyhedron( const int nbface, const IndexArray *face )
// {
//   m_cobound = new IndexArray[numVerts()];
//
//   IndexBuf* indexBuf = new IndexBuf[numVerts()];
//
//   //On initialise le volume et l'inertie a zero.
//   Initialisation();
//
//   Plan P;
//   Point H;
//   Point G_;// centre de gravite de la face
//   Point P1;
//   Point P2;
//   Point P3;
//   Vecteur u;
//
//   int i, j, k;
//   for(i=0; i<nbface; i++) {
//
//     G_.reset();
//     for(j=0, k=face[i].size()-1; j<face[i].size(); k=j++){
//       indexBuf[face[i][k]].push_back(face[i][j]);
//       G_ += (*this)[face[i][j]];
//     }
//     G_ /= face[i].size();
//
//     // calcul de la projection de l'origine sur la facette.
//     P1 = (*this)[face[i][0]];
//     P2 = (*this)[face[i][1]];
//     P3 = (*this)[face[i][2]];
//     P.setPlan(P1,P2,P3);
//     H = P.Projection(0,0,0);
//     for(j=0, k=face[i].size() - 1; j<face[i].size(); k=j++) {
//       P1 = (*this)[face[i][k]];
//       P2 = (*this)[face[i][j]];
//
//       // verification de l'orientation
//       u = G_;
//       if (triple((P1 - G_),(P2 - G_),u) < 0.){
// 	Point aux = P1;
// 	P1 = P2;
// 	P2 = aux;
//       }
//       // fin de verification
//
//       CalculVolumeInertie(H,P1,P2);
//     }
//   }
//
//   for (i = 0; i < numVerts(); ++i)
//     if (indexBuf[i].size())
//       new(&m_cobound[i]) IndexArray(indexBuf[i].size(), &indexBuf[i][0]);
//
//   m_curr_vertex = 0;
//   while (indexBuf[m_curr_vertex].size() == 0)
//     ++m_curr_vertex;
//
//   delete [] indexBuf;
// }
void Polyhedron::BuildPolyhedron( const int nbface, const IndexArray *face )
{
  m_cobound = new IndexArray[numVerts()];

  IndexBuf* indexBuf = new IndexBuf[numVerts()];

  //On initialise le volume et l'inertie a zero.
  Initialisation();

  Point G_ ;
  int i, j, k;
  for(i=0; i<nbface; i++)
  {
    for(j=0, k=face[i].size()-1; j<face[i].size(); k=j++)
      indexBuf[face[i][k]].push_back(face[i][j]);

    if ( face[i].size() == 3 )
    {
      CalculVolumeInertie( (*this)[face[i][0]], (*this)[face[i][1]],
      	(*this)[face[i][2]]);
    }
    else
    {
      G_.reset();
      for(j=0; j<face[i].size(); ++j) G_ += (*this)[face[i][j]];
      G_ /= face[i].size();
      for(j=0, k=face[i].size() - 1; j<face[i].size(); k=j++)
        CalculVolumeInertie( G_, (*this)[face[i][k]], (*this)[face[i][j]]);
    }
  }

  for (i = 0; i < numVerts(); ++i)
    if (indexBuf[i].size())
      new(&m_cobound[i]) IndexArray(int(indexBuf[i].size()), &indexBuf[i][0]);

  m_curr_vertex = 0;
  while (indexBuf[m_curr_vertex].size() == 0)
    ++m_curr_vertex;

  delete [] indexBuf;
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Lecture des faces
void Polyhedron::readface( istream &fileIn )
{
  char buffer[long_chaine];

  int nbface;
  fileIn >> nbface;
  fileIn.getline(buffer,sizeof(buffer)); // pour aller a la ligne suivante

  IndexArray *face = new IndexArray[nbface];
  IndexBuf facetind;

  // lecture des faces
  for(int i = 0;i < nbface;i++) {
    fileIn.getline(buffer,sizeof(buffer));
    istringstream lign(buffer);
    if (!fileIn.eof()) {
      int integer;
      for(;;){
	if (!lign.eof()) {
	  lign >> integer;
	  facetind.push_back(integer);
	}
	else break;
      }
      new(&face[i]) IndexArray(int(facetind.size()),&facetind[0]);
      facetind.erase(facetind.begin(), facetind.end());
    }
    else break;
  }

  BuildPolyhedron(nbface, face);

  // Stockage de la connectivite des faces
  m_allFaces = new vector< vector<int> >;
  m_allFaces->reserve( nbface );
  for(int i=0;i<nbface;i++)
  {
    vector<int> ww( face[i].size() );
    for (int j=0;j<face[i].size();++j) ww[j] = face[i][j];
    m_allFaces->push_back( ww );
  }

  delete [] face;
}




// ----------------------------------------------------------------------------
// Fonction support sans description des faces.
// Used in orginal SOLID-2.0 software
Point Polyhedron::support20( const Vecteur& v ) const
{
  int c = 0;
  Scalar h = (*this)[0] * v, d;
  for (int i=1; i<numVerts(); ++i) {
    if ((d = (*this)[i] * v) > h) { c=i; h=d; }
  }
  return (*this)[c];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Polyhedron::numberOfCells_PARAVIEW() const
{
  int ncorners = numVerts(), ncells = 0;

  // Tetrahedron or box
  if ( ncorners == 4 || ncorners == 8 ) ncells = 1;
  // Icosahedron
  else if ( ncorners == 12 && m_allFaces->size() == 20 ) ncells = 20;
  // Prism
  else if ( ncorners == 6 && m_allFaces->size() == 5 ) ncells = 1;
  // Octahedron
  else if ( ncorners == 6 && m_allFaces->size() == 8 ) ncells = 8;
   // Dodecahedron
  else if ( ncorners == 20 && m_allFaces->size() == 12 ) ncells = 36;
  // Trancoctahedron
  else if ( ncorners == 24 && m_allFaces->size() == 14 ) ncells = 44;
  // General: not implemented yet !!
  else
  {

  }

  return ( ncells );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview
void Polyhedron::write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
   int ncorners = numVerts();

   // Tetrahedron or box or prism
   if ( ncorners == 4 || ncorners == 8 ||
   	( ncorners == 6 && m_allFaces->size() == 5 ) )
   {
     int count = firstpoint_globalnumber + 1;
     for (int i=0;i<ncorners;++i)
     {
       connectivity.push_back(count);
       ++count;
     }
     last_offset += ncorners;
     offsets.push_back(last_offset);
     if ( ncorners == 4 ) cellstype.push_back(10);
     else if ( ncorners == 8 ) cellstype.push_back(12);
     else cellstype.push_back(13);

     firstpoint_globalnumber += ncorners + 1;
   }
   // Icosahedron or octahedron
   // The icosahedron is split into 20 tetrahedrons using the center of mass and
   // 3 vertices on each triangular face
   // The octahedron is split into 8 tetrahedrons using the center of mass and
   // 3 vertices on each triangular face
   else if ( ( ncorners == 12 && m_allFaces->size() == 20 ) ||
 	  ( ncorners == 6 && m_allFaces->size() == 8 ) )
   {
     size_t nbface = m_allFaces->size();
     for (size_t i=0;i<nbface;++i)
     {
       connectivity.push_back( firstpoint_globalnumber );
       size_t nbv = (*m_allFaces)[i].size();
       for (size_t j=0;j<nbv;++j)
         connectivity.push_back( firstpoint_globalnumber
 		+ (*m_allFaces)[i][j] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);
     }

     firstpoint_globalnumber += ncorners + 1;
   }
   // Dodecahedron
   // Each pentagonal face is split into 3 triangles
   // The dodecahedron is split into 12*3=36 tetrahedrons using the center of
   // mass and each triangular sub-face of each pentagonal face
   else if ( ncorners == 20 && m_allFaces->size() == 12 )
   {
     size_t nbface = m_allFaces->size();
     for (size_t i=0;i<nbface;++i)
     {
       // First tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][1] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][2] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);

       // Second tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][2] + 1 );
        connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][3] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);

       // Third tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][3] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][4] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);
     }

     firstpoint_globalnumber += ncorners + 1;
   }
   // Trancoctahedron
   // Each pentagonal face is split into 3 triangles
   // The trancoctahedron is split into 8*4 + 6*2 = 44 tetrahedrons using the
   // center of mass and each triangular sub-face of each pentagonal face
   else if ( ncorners == 24 && m_allFaces->size() == 14 )
   {
     size_t nbface = m_allFaces->size();
     for (size_t i=0;i<8;++i)
     {
       // First tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][1] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][2] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);

       // Second tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][2] + 1 );
        connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][3] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);

       // Third tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][3] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][4] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);

       // Fourth tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][4] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][5] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);
     }

     for (size_t i=8;i<nbface;++i)
     {
       // First tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][1] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][2] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);

       // Second tetrahedron
       connectivity.push_back( firstpoint_globalnumber );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][2] + 1 );
        connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][3] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                 + (*m_allFaces)[i][0] + 1 );
       last_offset += 4;
       offsets.push_back(last_offset);
       cellstype.push_back(10);
     }

     firstpoint_globalnumber += ncorners + 1;
   }
   // General: not implemented yet !!
   else
   {

   }   
  // int ncorners = numVerts();
  //
  // // Tetraedre ou parallelepipede ou prisme
  // if ( ncorners == 4 || ncorners == 8 ||
  // 	( ncorners == 6 && m_allFaces->size() == 5 ) )
  // {
  //   int count = firstpoint_globalnumber + 1;
  //   for (int i=0;i<ncorners;++i)
  //   {
  //     connectivity.push_back(count);
  //     ++count;
  //   }
  //   last_offset += ncorners;
  //   offsets.push_back(last_offset);
  //   if ( ncorners == 4 ) cellstype.push_back(10);
  //   else if ( ncorners == 8 ) cellstype.push_back(12);
  //   else cellstype.push_back(13);
  //
  //   firstpoint_globalnumber += ncorners + 1;
  // }
  // // Icosaedre
  // // L'icosaedre est decoupe en 20 tetraedres composes du centre de gravite
  // // et des 3 vertex de chaque face
  // else if ( ncorners == 12 && m_allFaces->size() == 20 )
  // {
  //   size_t nbface = m_allFaces->size();
  //   for (size_t i=0;i<nbface;++i)
  //   {
  //     connectivity.push_back( firstpoint_globalnumber );
  //     size_t nbv = (*m_allFaces)[i].size();
  //     for (size_t j=0;j<nbv;++j)
  //       connectivity.push_back( firstpoint_globalnumber
	// 	+ (*m_allFaces)[i][j] + 1 );
  //     last_offset += 4;
  //     offsets.push_back(last_offset);
  //     cellstype.push_back(10);
  //   }
  //
  //   firstpoint_globalnumber += ncorners + 1;
  // }
  // // General: not implemented yet !!
  // else
  // {
  //
  // }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe au format STL pour lien avec openFoam
void Polyhedron::write_convex_STL( ostream &f, const Transform &transform )
	const
{
  int ncorners = numVerts();

  // Parallelepipede
  if ( ncorners == 8 && m_allFaces->size() == 6 )
  {
    Point pp;
    Point GC = transform(pp);
    vector< Point > FC(4,pp);

    for (int i=0;i<int(m_allFaces->size());++i)
    {
      assert( (*m_allFaces)[i].size() == 4 ) ;
      pp = 0.25 * ( (*this)[(*m_allFaces)[i][0]]
      	+ (*this)[(*m_allFaces)[i][1]] + (*this)[(*m_allFaces)[i][2]]
	+ (*this)[(*m_allFaces)[i][3]] ) ;
      Point FaceCenterT = transform(pp);
      Vecteur outward_normal = FaceCenterT - GC;
      outward_normal.normalize();
      for (int j=0;j<4;++j)
        FC[j] = transform((*this)[(*m_allFaces)[i][j]]);

      // On divise la face rectangulaire en 2 triangles
      // Triangle 0
      f << "  facet normal " << outward_normal[X]
      	<< " " << outward_normal[Y]
      	<< " " << outward_normal[Z] << endl;
      f << "    outer loop" << endl;
      f << "      vertex " << FC[0][X] << " " << FC[0][Y] << " "
      	<< FC[0][Z] << endl;
      f << "      vertex " << FC[1][X] << " " << FC[1][Y] << " "
      	<< FC[1][Z] << endl;
      f << "      vertex " << FC[2][X] << " " << FC[2][Y] << " "
      	<< FC[2][Z] << endl;
      f << "    endloop" << endl;
      f << "  endfacet" << endl;
      // Triangle 1
      f << "  facet normal " << outward_normal[X]
      	<< " " << outward_normal[Y]
      	<< " " << outward_normal[Z] << endl;
      f << "    outer loop" << endl;
      f << "      vertex " << FC[0][X] << " " << FC[0][Y] << " "
      	<< FC[0][Z] << endl;
      f << "      vertex " << FC[2][X] << " " << FC[2][Y] << " "
      	<< FC[2][Z] << endl;
      f << "      vertex " << FC[3][X] << " " << FC[3][Y] << " "
      	<< FC[3][Z] << endl;
      f << "    endloop" << endl;
      f << "  endfacet" << endl;
    }
  }
  else
  {
    cout << "Warning for this Convex with ncorners = " << ncorners << " the "
       << "method Polyhedron::write_convex_STL() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
    exit(10);
  }
}
