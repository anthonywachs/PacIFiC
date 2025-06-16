#include "TruncatedConicalShell.hh"
#include "GrainsExec.hh"
#include "PointC.hh"
#include "ContactBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "TrapezoidalPrism.hh"


// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
TruncatedConicalShell::TruncatedConicalShell( DOMNode* root ) :
  CompositeObstacle( "shell" )
{    
  assert( root != NULL );
  
  // Type
  m_type = "TruncatedConicalShell";

  // Name
  m_name = ReaderXML::getNodeAttr_String( root, "name" );
  
  // Geometry
  DOMNode* nGeom = ReaderXML::getNode( root, "Geometry" );
  m_height = ReaderXML::getNodeAttr_Double( nGeom, "Height" );
  m_innerLargeRadius = ReaderXML::getNodeAttr_Double( nGeom, "LargeRadius" );
  m_innerSmallRadius = ReaderXML::getNodeAttr_Double( nGeom, "SmallRadius" );  
  m_shellWidth = ReaderXML::getNodeAttr_Double( nGeom, "Width" );
  m_nbBoxes = size_t(ReaderXML::getNodeAttr_Int( nGeom, "N" ));
  double crust_thickness = ReaderXML::getNodeAttr_Double( nGeom, 
  	"CrustThickness" ); 

  // Angular position of the truncated conical shell
  m_geoRBWC->getTransform()->load( root );
  
  // Bottom centre position of the truncated conical shell
  Point3 BottomCentre;
  DOMNode* npoint = ReaderXML::getNode( root, "BottomCentre" );
  BottomCentre[X] = ReaderXML::getNodeAttr_Double( npoint, "X" );
  BottomCentre[Y] = ReaderXML::getNodeAttr_Double( npoint, "Y" );
  BottomCentre[Z] = ReaderXML::getNodeAttr_Double( npoint, "Z" );
  
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
 
  double delta = 2. * PI / double(m_nbBoxes), angle, 
  	small_ly = 2. * ( m_innerSmallRadius + m_shellWidth ) *
		( sin( delta ) - ( 1. - cos( delta ) ) / tan( delta ) ),
	large_ly = 2. * ( m_innerLargeRadius + m_shellWidth ) *
		( sin( delta ) - ( 1. - cos( delta ) ) / tan( delta ) ),
	deltaR = m_innerLargeRadius - m_innerSmallRadius,
	lz = sqrt( m_height * m_height + deltaR * deltaR ),
	angleY = asin( deltaR / lz ),	
	extRad = m_innerSmallRadius 
  		+ ( small_ly + 2. * large_ly  ) * deltaR 
  		/ ( 3. * ( small_ly + large_ly ) ) + m_shellWidth / 2.,
	lt = ( small_ly + 2. * large_ly  ) * lz 
  		/ ( 3. * ( small_ly + large_ly ) ),	
	ht = sqrt( lt * lt - ( extRad - m_innerSmallRadius )
		* ( extRad - m_innerSmallRadius ) );
  string name;
  Matrix mrotZ, mrotY;
  Vector3 zerotranslation; 		 

  // We set the center of mass position such that the bottom is at the
  // specified position
  Point3 cg = BottomCentre;
  cg[Z] += ht;
  setPosition( cg );
  	
  // Create elememtary box obstacles  
  for (size_t i=0;i<m_nbBoxes;++i)
  {
    name = m_name + GrainsExec::intToString( int(i) );
    TrapezoidalPrism* box = new TrapezoidalPrism( m_shellWidth, small_ly, 
    	large_ly, lz );
    RigidBodyWithCrust* geoRBWC_box = new RigidBodyWithCrust( box, Transform(),
  	false, crust_thickness ); 
    Obstacle* sbox = new SimpleObstacle( name, geoRBWC_box, 
	m_materialName, transferToFluid, true );
    angle = ( 0.5 + double(i) ) * delta;
    cg.setValue( extRad * cos( angle ), extRad * sin( angle ),  0. );
    sbox->setPosition( cg );
    mrotY.setValue( cos( - angleY ), 0., - sin( - angleY ), 
    	0., 1., 0.,
	sin( - angleY ), 0., cos( - angleY ) );    
    mrotZ.setValue( cos( angle ), - sin( angle ), 0., 
    	sin( angle ), cos( angle ), 0., 
    	0., 0., 1. );
    sbox->getRigidBody()->getTransform()->setBasis( mrotZ * mrotY );
    sbox->getRigidBody()->composeLeftByTransform( *m_geoRBWC->getTransform() );
    sbox->Translate( zerotranslation );
    m_obstacles.push_back( sbox );	    
  }
  
  computeVolumeCenterOfMass();
}




// ----------------------------------------------------------------------------
// Constructor with name as input parameter
TruncatedConicalShell::TruncatedConicalShell( string const& s )
  : CompositeObstacle( s )
{
  m_type = "TruncatedConicalShell";  
} 




// ----------------------------------------------------------------------------
// Destructor
TruncatedConicalShell::~TruncatedConicalShell()
{}




// ----------------------------------------------------------------------------
// Outputs the truncated conical shell for reload
void TruncatedConicalShell::write( ostream& fileSave ) const
{
  fileSave << "<Composite> " << m_name << " " << m_type << endl;
  fileSave << "*Properties " << m_height << " " << m_innerSmallRadius << " " 
  	<< m_innerLargeRadius << " " << m_shellWidth << " " << m_nbBoxes 
	<< endl;  
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
// Reloads the truncated conical shell and links it to the higher level 
// obstacle in the obstacle tree
void TruncatedConicalShell::reload( Obstacle& mother, istream& file )
{
  string ttag, buffer;
  
  // Read extra properties 
  file >> buffer >> m_height >> m_innerSmallRadius >> m_innerLargeRadius >> 
  	m_shellWidth >> m_nbBoxes;
  
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
