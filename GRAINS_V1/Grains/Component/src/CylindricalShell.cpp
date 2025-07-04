#include "CylindricalShell.hh"
#include "GrainsExec.hh"
#include "PointC.hh"
#include "ContactBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "Box.hh"
#include "Rectangle.hh"


// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
CylindricalShell::CylindricalShell( DOMNode* root ) :
  CompositeObstacle( "shell" )
{    
  assert( root != NULL );
  
  // Type
  m_type = "CylindricalShell";

  // Name
  m_name = ReaderXML::getNodeAttr_String( root, "name" );
  
  // Geometry
  DOMNode* nGeom = ReaderXML::getNode( root, "Geometry" );
  m_height = ReaderXML::getNodeAttr_Double( nGeom, "Height" );
  m_innerRadius = ReaderXML::getNodeAttr_Double( nGeom, "Radius" );
  m_shellWidth = ReaderXML::getNodeAttr_Double( nGeom, "Width" );
  m_nbBoxes = size_t(ReaderXML::getNodeAttr_Int( nGeom, "N" ));
  double crust_thickness = ReaderXML::getNodeAttr_Double( nGeom, 
  	"CrustThickness" );
  m_box = false;
  if ( ReaderXML::hasNodeAttr( nGeom, "ElementaryObstacle" ) )
  {
    string selem = ReaderXML::getNodeAttr_String( nGeom, "ElementaryObstacle" );
    if ( selem == "Box" ) m_box = true;
  }   

  // Center of mass and angular position of the cylindrical shell
  m_geoRBWC->getTransform()->load( root );
  
  // Crust thickness 
  m_geoRBWC->setCrustThickness( crust_thickness );  

  // Material
  DOMNode* materiau_ = ReaderXML::getNode( root, "Material" );
  m_materialName = ReaderXML::getNodeValue_String( materiau_ );
  ContactBuilderFactory::defineMaterial( m_materialName, true );

  // Obstacle to transfer to the fluid
  bool transferToFluid = false;
  DOMNode* status = ReaderXML::getNode( root, "Status" );
  if ( status )
    transferToFluid = ReaderXML::getNodeAttr_Int( status, "ToFluid" );

  // Create elementary box obstacles 
  double extRad = m_innerRadius; 
  if ( m_box ) extRad += + m_shellWidth / 2.;
  double delta = 2. * PI / double(m_nbBoxes), angle;
  double ly = 2. * ( m_innerRadius + m_shellWidth ) *
	( sin( delta ) - ( 1. - cos( delta ) ) / tan( delta ) );	
  string name;
  Matrix mrot, mm( 0., 0., -1., 0., 1., 0., 1., 0., 0. );
  Vector3 zerotranslation;  
  for (size_t i=0;i<m_nbBoxes;++i)
  {
    name = m_name + GrainsExec::intToString( int(i) );
    RigidBodyWithCrust* geoRBWC_box = NULL;
    if ( m_box )
    { 
      Box* box = new Box( m_shellWidth, ly, m_height );
      geoRBWC_box = new RigidBodyWithCrust( box, Transform(),
  	false, crust_thickness ); 
    }
    else
    {
      Rectangle* box = new Rectangle( m_height, ly, m_shellWidth, RVE_ZMINUS );
      geoRBWC_box = new RigidBodyWithCrust( box, Transform(),
  	false, crust_thickness );       
    }
    Obstacle* sbox = new SimpleObstacle( name, geoRBWC_box, 
	m_materialName, transferToFluid, true );
    angle = ( 0.5 + double(i) ) * delta;
    Point3 cg( extRad * cos( angle ), extRad * sin( angle ),  0. );
    sbox->setPosition( cg );
    mrot.setValue( cos( angle ), - sin( angle ), 0., 
    	sin( angle ), cos( angle ), 0., 
    	0., 0., 1.);
    if ( !m_box ) mrot *= mm;	
    sbox->getRigidBody()->getTransform()->setBasis( mrot );
    sbox->getRigidBody()->composeLeftByTransform( *m_geoRBWC->getTransform() );
    sbox->Translate( zerotranslation );
    m_obstacles.push_back( sbox );	    
  }
  
  computeVolumeCenterOfMass();
}




// ----------------------------------------------------------------------------
// Constructor with name as input parameter
CylindricalShell::CylindricalShell( string const& s )
  : CompositeObstacle( s )
{
  m_type = "CylindricalShell";  
} 




// ----------------------------------------------------------------------------
// Destructor
CylindricalShell::~CylindricalShell()
{}




// ----------------------------------------------------------------------------
// Outputs the cylindrical shell for reload
void CylindricalShell::write( ostream& fileSave ) const
{
  fileSave << "<Composite> " << m_name << " " << m_type << endl;
  fileSave << "*Properties " << m_height << " " << m_innerRadius << " " 
  	<< m_shellWidth << " " << m_nbBoxes << endl;  
  if ( m_CompositeObstacle_id ) m_torsor.write( fileSave );
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
  {
    (*obstacle)->write( fileSave );
    fileSave << endl;
  }
  fileSave << "</Composite>";
}




// ----------------------------------------------------------------------------
// Reloads the cylindrical shell and links it to the higher level 
// obstacle in the obstacle tree
void CylindricalShell::reload( Obstacle& mother, istream& file )
{
  string ttag, buffer;
  
  // Read extra properties 
  file >> buffer >> m_height >> m_innerRadius >> m_shellWidth >> m_nbBoxes;
  
  // Standard Composite Obstacle reload
  if ( m_CompositeObstacle_id ) m_torsor.read( file ); 
  file >> ttag;
  while ( ttag != "</Composite>" ) 
  {
    ObstacleBuilderFactory::reload( ttag, *this, file );
    file >> ttag;
  }
  computeVolumeCenterOfMass();
  mother.append( this );
} 
