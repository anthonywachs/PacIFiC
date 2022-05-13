#include "Disque.H"

int Disque::visuNodeNb = 20;

// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
Disque::Disque(Scalar r) :
  radius(r),
  Shrinking(0),
  initial_radius(r)
{
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
// Constructeur avec fichier
Disque::Disque(istream &fileIn) :
  Shrinking(0),
  initial_radius(0)
{
  readClass(fileIn);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
// M. SULAIMAN MODIF Nov-2015 
Disque::Disque(DOMNode* root) :
  Convex()
{
  initial_radius    = ReaderXML::getNodeAttr_Double(root, "Radius");
  if ( ReaderXML::hasNodeAttr_String( root, "Shrink" ) )
  Shrinking = ReaderXML::getNodeAttr_Int(root,"Shrink");
  radius = initial_radius; 
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
Disque::~Disque() 
{
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
// Determine l'inertie et l'inertie inverse d'un disque
bool Disque::BuildInertie(Scalar *inertie,Scalar *inertie_1) const
{
  /*
  inertie[0]=inertie[1]=inertie[2]=inertie[4]=inertie[5]= 0.0;
  inertie[3] = PI * radius*radius / 2.0;
  */
  // DGR : PLAN 2D actif en XY -> rotation suivant axe Z
  inertie[1] = inertie[2] = inertie[4]= 0.0;
  inertie[0] = inertie[3] = (PI * radius*radius) * radius*radius/4.;
  inertie[5] = (PI * radius*radius) * radius*radius/2.;
  
  inertie_1[1] = inertie_1[2] = inertie_1[4] = 0.0;
  inertie_1[0] = 1.0 / inertie[0];
  inertie_1[3] = 1.0 / inertie[3];
  inertie_1[5] = 1.0 / inertie[5];

  return true;
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
// Determine le rayon circonscrit au disque
Scalar Disque::BuildRayonRef() const 
{
  return radius;
}




// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// return the shrinking choice
int Disque::getShrinkingChoice()const
{
  return Shrinking;
}




// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// Set the shrinking radius value
void Disque::setShrinkingRadius(Scalar CurrentRadius)
{
  radius = CurrentRadius;
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
// Construction d'un clone du disque
Convex* Disque::clone() const
{
  return new Disque(radius);
}




// ----------------------------------------------------------------------------
Point Disque::support(const Vecteur& v) const 
{
  Scalar s = Norm(v);
  if (s > EPSILON) {
    Scalar r = radius / s;
    return Point(v[X] * r, v[Y] * r, v[Z] * r);
  } else 
    return Point();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant le Convex pour visu dans une appli externe (Fluide)
// Pour un disque pas de prise en compte de points sommets -> 0
// G.FERRER - Aout.2004 - Creation
vector<Point> Disque::getEnveloppe() const
{
  vector<Point> enveloppe;
  Point point(0.,0.,0.);
  enveloppe.push_back(point);
  return enveloppe;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code équivalent
int Disque::getNbCorners() const
{
  return 1;
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
// Determine le Volume d'un disque
Scalar Disque::getVolume()const 
{
  return ( PI * radius * radius );
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
// ecriture du disque
void Disque::printClass(ostream &fileOut) const 
{
  fileOut << "*Disque\n";
  if(Shrinking!=1)
    fileOut <<radius<<endl;  
  else
  {
    fileOut <<radius;
    fileOut <<" "<<Shrinking<<'\n'; 
  }
}




// ----------------------------------------------------------------------------
// G.FERRER - Mars.2003 - Creation
// M.SULAIMAN - Nov.2015 - Modification
// lecture du disque
void Disque::readClass(istream &fileIn) 
{
  fileIn >> radius;
  if(Shrinking==1) fileIn >> Shrinking; 
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview 
int Disque::numberOfPoints_PARAVIEW() const 
{ 
  return (visuNodeNb);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview  
// radius replaced by radius
void Disque::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation) const
{
  double pi=acos(-1.);
  double angle = 2.*pi/visuNodeNb ;
  int i;
  Point pp,pptrans;
  
  // Regular points on the surface
  for ( i = 0; i < visuNodeNb ; ++i )
  {
    pp[X] =radius*cos(i*angle);
    pp[Y] =radius*sin(i*angle);
    pptrans = transform(pp);
    if ( translation ) pptrans += *translation;
    f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview  
// radius replaced by radius
list<Point> Disque::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation) const
{
  list<Point> ParaviewPoints;
  double pi=acos(-1.);
  double angle = 2.*pi/visuNodeNb ;
  int i;
  Point pp,pptrans;
  
  // Regular points on the surface
  for ( i = 0; i < visuNodeNb ; ++i )
  {
    pp[X] = radius*cos(i*angle);
    pp[Y] = radius*sin(i*angle);
    pptrans = transform(pp);
    if ( translation ) pptrans += *translation;
    ParaviewPoints.push_back(pptrans);
  } 
  
  return ParaviewPoints;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Disque::numberOfCells_PARAVIEW() const
{
  return 1; 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview 
void Disque::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  int ncorners = visuNodeNb;

  int count=firstpoint_globalnumber;
  for (int i=0;i<ncorners;++i)
  {
    connectivity.push_back(count);
    ++count;
  }  
  last_offset+=ncorners;    
  offsets.push_back(last_offset);
  cellstype.push_back(7);
  
  firstpoint_globalnumber+=ncorners;
}
