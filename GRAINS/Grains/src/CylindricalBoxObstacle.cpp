#include "MPIWrapperGrains.hh"
#include "CylindricalBoxObstacle.H"
#include "Obstacle_BuilderFactory.H"
#include "Box.H"
#include "Cylinder.H"
#include "AntiCylinder.H"
#include "PointC.H"


// ----------------------------------------------------------------------------
// Default constructor
CylindricalBoxObstacle::CylindricalBoxObstacle( const string &s,
    const double &radius_, const double &length_ ) :
  CompObstacle( s ),
  m_length( length_ ),
  m_radius( radius_ ),
  m_transfer2Fluid( true ),
  m_bottomPanel( false ),
  m_topPanel( false )
{
  m_ObstacleType = "CylindricalBox";
  
  Convex *cylinder = new AntiCylinder( m_radius, m_length );

  Transform cylinderTransform;
  Point cylinderCenter;

  cylinderCenter[X] = 0.;
  cylinderCenter[Y] = 0.;
  cylinderCenter[Z] = 0.;
  cylinderTransform.setOrigin( cylinderCenter );

  Matrix cylinderRotation(
      0., 0., -1.,
      1., 0., 0., 
      0., -1., 0. );

  cylinderTransform.setOrigin( cylinderCenter );
  cylinderTransform.setBasis( cylinderRotation );
  
  m_geoFormeVdw = new FormeVdW( cylinder, cylinderTransform );

}




// ----------------------------------------------------------------------------
// Constructeur par decodage d'un noeud XML
CylindricalBoxObstacle::CylindricalBoxObstacle( DOMNode *root ) :
  CompObstacle(),
  m_transfer2Fluid( true ),
  m_bottomPanel( false ),
  m_topPanel( false )
{
  m_ObstacleType = "CylindricalBox";
  
  m_nom = ReaderXML::getNodeAttr_String( root, "Name" );
  
  DOMNode* geometry = ReaderXML::getNode( root, "Geometry" );
  m_length = ReaderXML::getNodeAttr_Double( geometry, "Length" );
  m_radius = ReaderXML::getNodeAttr_Double( geometry, "Radius" );

  DOMNode* centerNode = ReaderXML::getNode( root, "Center" );
  m_center[X] = ReaderXML::getNodeAttr_Double( centerNode, "X" );
  m_center[Y] = ReaderXML::getNodeAttr_Double( centerNode, "Y" );
  m_center[Z] = ReaderXML::getNodeAttr_Double( centerNode, "Z" );

  DOMNode* panels = ReaderXML::getNode( root, "Panels" );
  m_thickness = ReaderXML::getNodeAttr_Double( panels, "Thickness" );
  m_nbPanels = size_t( ReaderXML::getNodeAttr_Int( panels, "Number" ));
  m_polygonAngle = 2.*PI/double(m_nbPanels);
  

  DOMNode* rayonInteraction = ReaderXML::getNode( root, "RayonInteraction" );
  m_overlapMax = ReaderXML::getNodeAttr_Double( rayonInteraction, "Delta" );

  DOMNode* materiau_ = ReaderXML::getNode( root, "Materiau" );
  m_nomMateriau = ReaderXML::getNodeValue_String( materiau_ );

  DOMNode* closeBox = NULL;
  closeBox = ReaderXML::getNode( root, "CloseBox" );
  if ( closeBox )
  {
    m_bottomPanel = ReaderXML::getNodeAttr_Int( closeBox, "Bottom" );
    m_topPanel = ReaderXML::getNodeAttr_Int( closeBox, "Top" );
  }
    
  // Obstacle seen by fluid ?
  DOMNode* statut = NULL;
  statut = ReaderXML::getNode( root, "Statut" );
  if ( statut )
    m_transfer2Fluid = ReaderXML::getNodeAttr_Int( statut, "ToFluid" ); 

  generate_sidePanels();
  
  if ( m_bottomPanel || m_topPanel )
    generate_topAndBotomPanels();

  EvalPosition();
  
}




// ----------------------------------------------------------------------------
// Destructeur
CylindricalBoxObstacle::~CylindricalBoxObstacle()
{}




// ----------------------------------------------------------------------------
// Generates side panels from informations read in the constructor
void CylindricalBoxObstacle::generate_sidePanels()
{
  double Lx=0., Ly=0., Lz=0., angle=0., polygonRadius=0., polygonApothem=0.;
  Transform panelPosition;
  Point panelCenter;
  FormeVdW *geoFormeVdw = NULL;
  Obstacle *obstacle = NULL;
  Convex *panel = NULL;
  
  // polygon edge-length for A_polygon = A_disk
  Lx = sqrt(4.*PI*pow(m_radius,2.)*tan(m_polygonAngle/2.)/
  	double(m_nbPanels));
  
  // We enlarge this length such that exterior corners collapse (estetic)
  Lx += 2.*m_thickness*tan(m_polygonAngle/2.);
  
  // Polygon Radius and Apothem for A_polygon = A_disk
  polygonRadius = sqrt( 2.*PI*pow(m_radius,2.)/
  	(double(m_nbPanels)*sin(2.*PI/double(m_nbPanels))) ) ;
  polygonApothem = polygonRadius * cos(m_polygonAngle/2.);
  
  Ly = m_length;
  Lz = m_thickness;

  for (size_t iNb=0; iNb!=m_nbPanels; iNb++)
  {
    angle = 2.*PI * double(iNb) / double(m_nbPanels);
    
    // Coordinates of each panel gravity center
    panelCenter[X] = m_center[X] + (polygonApothem+m_thickness/2.) * cos(angle);
    panelCenter[Y] = m_center[Y] + (polygonApothem+m_thickness/2.) * sin(angle);
    panelCenter[Z] = m_center[Z];
    panelPosition.setOrigin( panelCenter );
    
    Matrix panelRotation(
        -sin(angle), 0., -cos(angle),
        cos(angle), 0., -sin(angle),
        0., -1., 0. );
    panelPosition.setBasis( panelRotation );
    
    // We create a reference geometrical box
    panel = new Box(Lx,Ly,Lz);

    // We create a FormeVdW for each panel at the right position
    geoFormeVdw = new FormeVdW( panel, panelPosition );
    geoFormeVdw->setRayonInterAction(m_overlapMax);
    
    // Then we create an obstacle with those individual panels
    std::stringstream ss;
    ss << m_nom <<"_" << iNb;
    string name = ss.str();
    obstacle = new MonObstacle( geoFormeVdw, name, m_nomMateriau,
    	false );
    
    m_obstacles.push_back(obstacle);
  }
  
}




// ----------------------------------------------------------------------------
// Generates side panels from informations read in the constructor
void CylindricalBoxObstacle::generate_topAndBotomPanels()
{
  double Lx=0., Ly=0., Lz=0.;
  Transform panelPosition;
  Point panelCenter;
  FormeVdW *geoFormeVdw = NULL;
  Obstacle *obstacle = NULL;
  Convex *panel = NULL;

  Lx=2.*(m_radius + m_thickness);
  Ly=2.*(m_radius + m_thickness);
  Lz = m_thickness;

  bool cylinder = false;
  Matrix cylinderRotation(
      0., 0., -1.,
      1., 0., 0., 
      0., -1., 0. );
  if( cylinder ) panelPosition.setBasis( cylinderRotation );

  if ( m_bottomPanel )
  {
    // Coordinates of each panel gravity center
    panelCenter[X] = m_center[X] ;
    panelCenter[Y] = m_center[Y] ;
    panelCenter[Z] = m_center[Z]-(m_length+m_thickness)/2. ;
    panelPosition.setOrigin( panelCenter );

    // We create a convex
    if( cylinder ) 
      panel = new Cylinder( m_radius+m_thickness, m_thickness );
    else
      panel = new Box( Lx,Ly,Lz );

    // We create a FormeVdW for each panel at the right position
    geoFormeVdw = new FormeVdW( panel, panelPosition );
    geoFormeVdw->setRayonInterAction(m_overlapMax);
    
    // Then we create an obstacle with those individual panels
    std::stringstream ss;
    ss << m_nom <<"_Bottom";
    string name = ss.str();
    obstacle = new MonObstacle( geoFormeVdw, name, m_nomMateriau, 
    	false );
    
    m_obstacles.push_back(obstacle);
  }


  if ( m_topPanel )
  {
    // Coordinates of each panel gravity center
    panelCenter[X] = m_center[X] ;
    panelCenter[Y] = m_center[Y] ;
    panelCenter[Z] = m_center[Z]+(m_length+m_thickness)/2. ;
    panelPosition.setOrigin( panelCenter );

    // We create a reference geometrical box
    if( cylinder ) 
      panel = new Cylinder( m_radius+m_thickness, m_thickness );
    else
      panel = new Box( Lx,Ly,Lz );
    
    // We create a FormeVdW for each panel at the right position
    geoFormeVdw = new FormeVdW( panel, panelPosition );
    
    // Then we create an obstacle with those individual panels
    std::stringstream ss;
    ss << m_nom <<"_Top";
    string name = ss.str();
    obstacle = new MonObstacle( geoFormeVdw, name, m_nomMateriau, 
    	false );
    
    m_obstacles.push_back(obstacle);
  }
  
}




// ----------------------------------------------------------------------------
// Sauvegarde de l'obstacle compose pour Reload
void CylindricalBoxObstacle::write( ostream &fileSave, 
	Composant const* composant ) const
{
  fileSave << "<CylindricalBox>\t" << m_nom << '\n';
  fileSave << m_radius << '\n';
  fileSave << m_length << '\n';
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
  {
    (*obstacle)->write( fileSave );
  }
  fileSave << "</CylindricalBox>\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reload du groupe d'Obstacle & publication dans son referent
void CylindricalBoxObstacle::reload( Obstacle &obstacle, istream &file )
{
  string ttag;
  file >> ttag;
  
  while ( ttag != "</CylindricalBox>" ) 
  {
    Obstacle_BuilderFactory::reload( ttag, *this, file );
    file >> ttag;
  }
  EvalPosition();
  obstacle.append(this);
}    




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Liste des obstacles a transmettre au fluide
list<Obstacle*> CylindricalBoxObstacle::getObstaclesToFluid()
{
  list<Obstacle*> liste;
  if ( m_transfer2Fluid == true )
  {
    liste.push_back(this);
  }
  return liste;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le type de l'obstacle
// 
string CylindricalBoxObstacle::getObstacleType()
{
  return m_ObstacleType;
}  




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant le Convex pour visu dans une appli externe (Fluide)
vector<Point> CylindricalBoxObstacle::getEnveloppe() const
{
  Point point(0.,0.,0.);
  vector<Point> enveloppe(3,point);
  enveloppe[0][Y] = - m_length/2.;
  enveloppe[1][Y] = - m_length/2.;
  enveloppe[1][X] = m_radius;  
  enveloppe[2][Y] = m_length/2.; 
  return enveloppe;
}




