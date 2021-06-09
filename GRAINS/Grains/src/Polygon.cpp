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
#include "Grains_Exec.hh"
#include "Polygon.H"
#include "EnsComposant.H"
#include "Point.H"
#include "Vecteur.H"
#include <fstream>
#include <sstream>
#include <vector>
typedef std::vector<unsigned int> IndexBuf;



// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction avec un fichier
Polygon::Polygon( istream &fileIn, int nb_point, VertexBase& ref,
  	 IndexArray& ia ) :
  Polytope( fileIn, nb_point, ref, ia ),
  m_cobound( NULL ),
  m_curr_vertex( 0 ), 
  m_InertiePoly( NULL ), 
  m_surface( 0. )
{
  readface( fileIn );
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Constructeur de copie
Polygon::Polygon( const Polygon &copie ) : 
  Polytope( copie )
{
  m_cobound = copie.m_cobound;
  m_curr_vertex = copie.m_curr_vertex;
  m_InertiePoly = NULL;
  m_surface = copie.m_surface;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Polygon::~Polygon() 
{
  // Le tenseur d'inertie n'est construit que pour les particules de reference
  // de classe, m_InertiePoly est NULL si la particule est active
  if ( m_InertiePoly ) delete m_InertiePoly;
  
  // m_cobound est detruit par le mecanisme de garbage collector de Grains_Exec 
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// creation d'un Polygon
Polygon* Polygon::create( istream &fileIn ) 
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
  
  // Creation du tableau de points du polygone    
  Point *point = new Point[nb_point];
  VertexBase *vertexbase = new VertexBase((void *)point);
  
  // Creation du tableau d'indices des sommets dans le tableau de points
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v ); 
  delete [] v; 

  // Creation du polygone  
  Polygon *polygon = new Polygon( *PolyIN, nb_point, *vertexbase, *ia );
  polygon->m_fichPoly = fichpoly;
  PolyIN->close();
  
  // Les objets m_base, m_cobound et m_indexne sont crees qu'une seule fois par
  // type de particule, les particules suivantes ne possedent qu'un pointeur sur
  // ces derniers.
  // Un mecanisme de garbage collector a ete mis en place dans la classe
  // Grains_Exec qui conserve les pointeurs et se charge de les detruire
  Grains_Exec::addOnePolytopeRefPointBase( point, vertexbase );
  Grains_Exec::addOnePolytopeNodeNeighbors( polygon->m_cobound );
  Grains_Exec::addOnePolytopeNodeIndex( ia );  
  
  delete PolyIN;
  
  return polygon;
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int Polygon::getShrinkingChoice()const
{
  return(0);
}




// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// Set the shrinking radius value
void Polygon:: setShrinkingRadius(Scalar CurrentRadius)
{
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation d'un polygone
Polygon* Polygon::create( DOMNode* root ) 
{
  // Lecture du nom puis ouverture du fichier contenant la description 
  // du polygone
  string fichpoly = ReaderXML::getNodeAttr_String( root, "Name" );
  ifstream PolyIN( fichpoly.c_str() );
  if ( PolyIN.is_open() ) 
    Grains_Exec::m_additionalDataFiles.insert( fichpoly );
      
  // Lecture du nb de points definissant le polygone
  int nb_point;
  PolyIN >> nb_point >> nb_point;

  // Creation du tableau de points du polygone
  Point *point = new Point[nb_point];
  VertexBase *vertexbase = new VertexBase((void *)point);

  // Creation du tableau d'indices des sommets dans le tableau de points
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v ); 
  delete [] v; 

  // Creation du polygone
  Polygon *polygon = new Polygon( PolyIN, nb_point, *vertexbase, *ia );
  polygon->m_fichPoly = fichpoly;
  
  // Les objets m_base, m_cobound et m_indexne sont crees qu'une seule fois par
  // type de particule, les particules suivantes ne possedent qu'un pointeur sur
  // ces derniers.
  // Un mecanisme de garbage collector a ete mis en place dans la classe
  // Grains_Exec qui conserve les pointeurs et se charge de les detruire
  Grains_Exec::addOnePolytopeRefPointBase( point, vertexbase );
  Grains_Exec::addOnePolytopeNodeNeighbors( polygon->m_cobound );
  Grains_Exec::addOnePolytopeNodeIndex( ia );  
 
  return polygon;
}




// ----------------------------------------------------------------------
// F.PRADEL - Avril 2003 - Creation
// Lecture des arêtes 
void Polygon::readface( istream &fileIn )
{
  char buffer[long_chaine];
  
  int nbedge;
  fileIn >> nbedge;
  fileIn.getline(buffer,sizeof(buffer)); // pour aller a la ligne suivante
  
  IndexArray *edge = new IndexArray[nbedge];
  IndexBuf facetind;
  
  // lecture des faces
  for(int i = 0;i < nbedge;i++) {
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
      new(&edge[i]) IndexArray(int(facetind.size()),&facetind[0]);
      facetind.erase(facetind.begin(), facetind.end());
    }
    else break;
  }
  
  BuildPolygon(edge[0].size(), edge);
  
  delete [] edge;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un Polygon
bool Polygon::BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const
{
  std::copy(&m_InertiePoly[0], &m_InertiePoly[6], &inertie[0]);

  inertie_1[1] = inertie_1[2] = inertie_1[4] = 0.0;
  inertie_1[0] = 1.0 / inertie[0];
  inertie_1[3] = 1.0 / inertie[3];
  inertie_1[5] = 1.0 / inertie[5];

  return true;
}
  



// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Scalar Polygon::BuildRayonRef() const 
{
  Scalar d , ray = ((*this)[0]) * ((*this)[0]);
  for (int i = 1; i< numVerts(); i++) {
    if ((d =  (*this)[i] * (*this)[i]) > ray)
      ray = d;
  }
  return sqrt(ray);
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone d'un Polygon
Convex* Polygon::clone() const
{
  return new Polygon(*this);
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un Polygon
Scalar Polygon::getVolume() const 
{
  return m_surface;
}




// ----------------------------------------------------------------------------
// Initialisaton de l'inertie et de la surface
// F.PRADEL - Avri.2003 - Creation
void Polygon::Initialisation()
{
  m_surface = 0.;
  m_InertiePoly = new Scalar[6];
  m_InertiePoly[0] = m_InertiePoly[1] = m_InertiePoly[2] = m_InertiePoly[3] = 
    m_InertiePoly[4] = m_InertiePoly[5] = 0.;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Point Polygon::support( const Vecteur& v ) const 
{
  Scalar norm = Norm(v);
  if (norm > EPSILON) {
    const Point ptTmp = (*this)[m_curr_vertex];
    Scalar h = ptTmp * v;
    Scalar d;
    int ni = int(m_curr_vertex) < numVerts()-1 ? m_curr_vertex+1 : 0;
    if ((d = (*this)[ni] * v) > h) { 
      do { 
	h = d; m_curr_vertex = ni;
	if (++ni == numVerts()) ni = 0;
      }
      while ((d = (*this)[ni] * v) > h);
    }
    else {
      ni = m_curr_vertex ? m_curr_vertex-1 : numVerts()-1;
      while ((d = (*this)[ni] * v) > h) {
	h = d; m_curr_vertex = ni;
	if (ni) --ni; else ni = numVerts()-1;
      }
    }  
    return (*this)[m_curr_vertex];
  } else {
    return Point();
  }
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture du polygon
void Polygon::printClass( ostream &fileOut ) const 
{
  fileOut << "*Polygon\n";
  fileOut << Grains_Exec::extractFileName( m_fichPoly ) << endl;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture du polygon
void Polygon::readClass( istream &fileIn ) 
{
  cerr << "Program Error :\n"
       << "Polygon::readClass non accessible.\n";
  exit(3);
}




// ----------------------------------------------------------------------
// F.PRADEL - Juin. 2003 - Creation
// Construction du polyhedron
void Polygon::BuildPolygon( const int nbedge, const IndexArray *edge ) 
{  
  m_cobound = new IndexArray[numVerts()];
  
  IndexBuf* indexBuf = new IndexBuf[numVerts()];
  
  //On initialise le volume et l'inertie a zero.
  Initialisation();
  
  Point G_;// centre de gravite de la face
  Point P1;
  Point P2;
  Vecteur u;
  
  int i;
  //  cout<<" le polygon a "<<nbedge<<" coté"<<endl;
  for (i=0; i<nbedge-1;i++) {
    G_ += (*this)[edge[0][i]];
  }
  G_ /=(1.0*nbedge);
  for(i = 0;i < nbedge-1; i++) { 
    // on prend les deux extrémités de l'arête.  
    P1 = (*this)[edge[0][i]];
    P2 = (*this)[edge[0][i+1]];
    // verification de l'orientation
    
    u = G_;
    if (triple((P1 - G_),(P2 - G_),u) < 0.){
      Point aux = P1;
      P1 = P2;
      P2 = aux;
    }
    // fin de verification

    CalculSurfaceInertie(P1-G_,P2-G_);
  }
  
  for (i = 0; i < numVerts(); ++i) 
    if (indexBuf[i].size()) 
      new(&m_cobound[i]) IndexArray(int(indexBuf[i].size()), &indexBuf[i][0]);   
  
  m_curr_vertex = 0;
  while (indexBuf[m_curr_vertex].size() == 0) ++m_curr_vertex;
  
  delete [] indexBuf;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Avril 2003 - Creation
// calcul de l'inertie et de la surface d'un triangle rectangle
void Polygon::CalculSurfaceInertie( const Point &P, const Point &Q )
{
  // On effectue d'abord un changement de repere
  // Construction de la matrice de rotation
  //  Vecteur x1(1,0,0), y1(0,1,0), z1(0,0,1);
  //attention ici il faut savoir qu'elles sont les 2 dimensions gardées.
	//(X,Y) ou (X,Z)
  Vecteur x2(P), y2(Q - P);
  Point gtemp((P+Q)*1./3.);
  Point Pt(P);
  Point Qt(Q);
  Pt=Pt-gtemp;
  Qt=Qt-gtemp;
  //cout<<"centre de gravité \t" <<gtemp<<endl;
  //cout<<"P\t"<<Pt<<endl;
  //cout<<"Q\t"<<Qt<<endl;
  Vecteur z2 = x2^y2; 
  
  //  Matrix Rot(cos(x1,x2), cos(x1,y2), cos(x1,z2),
  //	     cos(y1,x2), cos(y1,y2), cos(y1,z2),
  //	     cos(z1,x2), cos(z1,y2), cos(z1,z2));
  
  //  Matrix tRot = Rot.transpose();
  
  // on se ramene a la base relative
//  Scalar a = tRot[X] * x2;
 // Scalar b = tRot[Y][X] * P2[X] + tRot[Y][Y] * P2[Y] + tRot[Y][Z] * P2[Z];
  //Scalar c = tRot[Z][X] * P2[X] + tRot[Z][Y] * P2[Y] + tRot[Z][Z] * P2[Z];
  //Scalar d = tRot[Y] * y2;
  // pourquoi utiliser deux façons différentes de l'écrire ?
  
  Scalar S = Norm(z2)/2.0;
  //  cout<<" surface triangle "<<S<<endl;
  m_surface += S;
   
  //  Scalar inertie[11];
  Scalar delta_12 = S/12.0;
  //  inertie[0] = 0.0;
  //  inertie[1] = inertie[4] = 0.0;
  //  inertie[2] = inertie[8] = 0.0;
  //  inertie[3] = 0.0;
  //  inertie[5] = 0.0;
  //  inertie[6] = inertie[9] = 0.0;
  //  inertie[7] = 0.0;
  //  inertie[10] = delta_12*( Pt*Pt+Qt*Qt+gtemp*gtemp)+ S*(gtemp*gtemp );
  //  Matrix I(inertie);

  // changement de base (relative a l'absolue)
  //I = Rot * (I * tRot);
  //I += GInertie;

  //m_InertiePoly[0] += (I.getValue())[0][0];
   //m_InertiePoly[1] += (I.getValue())[0][1];
  //m_InertiePoly[2] += (I.getValue())[0][2];
  //(I.getValue())[1][1];
  //cout<< "Inertie \t"<<m_InertiePoly[3]<<endl;
  //m_InertiePoly[4] += (I.getValue())[1][2];
  //m_InertiePoly[5] += (I.getValue())[2][2];

  // DGR : PLAN 2D actif en XY -> rotation suivant axe Z
  m_InertiePoly[1]  = m_InertiePoly[2] = m_InertiePoly[4] = 0.0;
  m_InertiePoly[0]  = m_InertiePoly[3] = 1.0;
  m_InertiePoly[5] += delta_12*( Pt*Pt+Qt*Qt+gtemp*gtemp)+ S*(gtemp*gtemp );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Polygon::numberOfCells_PARAVIEW() const
{
  return 1;  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview 
void Polygon::write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int ncorners = numVerts();

  int count = firstpoint_globalnumber + 1;
  for (int i=0;i<ncorners;++i)
  {
    connectivity.push_back(count);
    ++count;
  }  
  last_offset += ncorners;    
  offsets.push_back(last_offset);
  cellstype.push_back(7);
  
  firstpoint_globalnumber += ncorners + 1;
}
