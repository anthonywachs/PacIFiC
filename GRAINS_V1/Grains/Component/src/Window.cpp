#include "Window.hh"	
#include "GrainsExec.hh"
#include "RigidBody.hh"
#include "Box.hh"
#include "Cylinder.hh"
#include "Segment.hh"


// ----------------------------------------------------------------------------
// Default constructor
Window::Window()
  : m_ftype( WINDOW_NONE )
  , m_radius( 0. )
  , m_radius_int( 0. )
  , m_height( 0. )
  , m_axisdir( NONE )
  , m_name( "unknown" )
  , m_kinematics( NULL )
{}




// ----------------------------------------------------------------------------
// Copy constructor 
Window::Window( Window const& copy )
  : m_ftype( copy.m_ftype )
  , m_ptA( copy.m_ptA )
  , m_ptB( copy.m_ptB )
  , m_radius( copy.m_radius )
  , m_radius_int( copy.m_radius_int )
  , m_height( copy.m_height )
  , m_axisdir( copy.m_axisdir )
  , m_name( copy.m_name )
  , m_kinematics( copy.m_kinematics )
{}




// ----------------------------------------------------------------------------
// Destructor
Window::~Window()
{
  if ( m_kinematics ) delete m_kinematics;
}




// ----------------------------------------------------------------------------
// Reads a window from an XML node
bool Window::readWindow( DOMNode* nWindow, string const& oshift,
    	int const& rank )
{
  bool ok = true;

  // Insertion window type
  string iwindow_type = "Box";
  if ( ReaderXML::hasNodeAttr( nWindow, "Type" ) )
    iwindow_type = ReaderXML::getNodeAttr_String( nWindow,
		  	"Type" );

  if ( iwindow_type == "Cylinder" )
    m_ftype = WINDOW_CYLINDER;
  else if ( iwindow_type == "Annulus" )
    m_ftype = WINDOW_ANNULUS;
  else if ( iwindow_type == "Line" )
    m_ftype = WINDOW_LINE;
  else if ( iwindow_type == "Box" )
    m_ftype = WINDOW_BOX;
  else m_ftype = WINDOW_NONE;
  
  // Insertion window name
  if ( ReaderXML::hasNodeAttr( nWindow, "Name" ) )
    m_name = ReaderXML::getNodeAttr_String( nWindow,
		  	"Name" );    

  DOMNodeList* points = NULL;
  DOMNode* pointA = NULL;
  DOMNode* pointB = NULL;
  DOMNode* cylGeom = NULL;
  DOMNode* annGeom = NULL;
  string axisdir_str = "W";

  if ( rank == 0 )
  {
    cout << oshift << "Name = " << m_name << endl;
    cout << oshift << "Type = " << iwindow_type << endl;
  }

  // Read features depending on the type
  switch( m_ftype )
  {
    case WINDOW_BOX:
      points = ReaderXML::getNodes( nWindow );
      pointA = points->item( 0 );
      pointB = points->item( 1 );

      m_radius = m_radius_int = m_height = 0. ;
      m_axisdir = NONE ;
      m_ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      m_ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      m_ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      m_ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
      m_ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
      m_ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );

      if ( rank == 0 )
      {
	cout << oshift << GrainsExec::m_shift3 << "Point3 min = " <<
		m_ptA[X] << " " << m_ptA[Y] << " " <<
		m_ptA[Z] << endl;
        cout << oshift << GrainsExec::m_shift3 << "Point3 max = " <<
		m_ptB[X] << " " << m_ptB[Y] << " " <<
		m_ptB[Z] << endl;
      }
      break;

    case WINDOW_CYLINDER:
      pointA = ReaderXML::getNode( nWindow, "BottomCentre" );
      m_ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      m_ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      m_ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      m_ptB[X] = m_ptB[Y] = m_ptB[Z] = 0.;
      cylGeom = ReaderXML::getNode( nWindow, "Cylinder" );
      m_radius = ReaderXML::getNodeAttr_Double( cylGeom, "Radius" );
      m_radius_int = 0.;
      m_height = ReaderXML::getNodeAttr_Double( cylGeom, "Height" );
      axisdir_str = ReaderXML::getNodeAttr_String( cylGeom, "Direction" );
      if ( axisdir_str == "X" ) m_axisdir = X;
      else if ( axisdir_str == "Y" ) m_axisdir = Y;
      else if ( axisdir_str == "Z" ) m_axisdir = Z;
      else
      {
	if ( rank == 0 )
          cout << "Wrong axis direction in cylindrical "
		<< "insertion window; values: X, Y or Z" << endl;
        ok = false;
      }

      if ( rank == 0 )
      {
	cout << oshift << GrainsExec::m_shift3 << "Bottom centre = " <<
		m_ptA[X] << " " << m_ptA[Y] << " " <<
		m_ptA[Z] << endl;
	cout << oshift << GrainsExec::m_shift3 << "Radius = " <<
	  	m_radius << endl;
        cout << oshift << GrainsExec::m_shift3 << "Height = " <<
		m_height << endl;
        cout << oshift << GrainsExec::m_shift3 << "Direction = " <<
		axisdir_str << endl;
      }
      break;

    case WINDOW_ANNULUS:
      pointA = ReaderXML::getNode( nWindow, "BottomCentre" );
      m_ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      m_ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      m_ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      m_ptB[X] = m_ptB[Y] = m_ptB[Z] = 0.;
      annGeom = ReaderXML::getNode( nWindow, "Annulus" );
      m_radius = ReaderXML::getNodeAttr_Double( annGeom,
	  	"RadiusExt" );
      m_radius_int = ReaderXML::getNodeAttr_Double( annGeom,
		"RadiusInt" );
      m_height = ReaderXML::getNodeAttr_Double( annGeom,
		"Height" );
      axisdir_str = ReaderXML::getNodeAttr_String( annGeom,
		"Direction" );
      if ( axisdir_str == "X" ) m_axisdir = X;
      else if ( axisdir_str == "Y" ) m_axisdir = Y;
      else if ( axisdir_str == "Z" ) m_axisdir = Z;
      else
      {
        if ( rank == 0 )
          cout << "Wrong axis direction in cylindrical "
		<< "insertion window; values: X, Y or Z" << endl;
        ok = false;
      }

      if ( rank == 0 )
      {
	cout << oshift << GrainsExec::m_shift3 << "Bottom centre = " <<
		m_ptA[X] << " " << m_ptA[Y] << " " <<
		m_ptA[Z] << endl;
	cout << oshift << GrainsExec::m_shift3 << "External radius = " <<
		m_radius << endl;
        cout << oshift << GrainsExec::m_shift3 << "Internal radius = " <<
		m_radius_int << endl;
        cout << oshift << GrainsExec::m_shift3 << "Height = " <<
		m_height << endl;
        cout << oshift << GrainsExec::m_shift3 << "Direction = " <<
		axisdir_str << endl;
      }
      break;

    case WINDOW_LINE:
      points = ReaderXML::getNodes( nWindow );
      pointA = points->item( 0 );
      pointB = points->item( 1 );

      m_radius = m_radius_int = m_height = 0. ;
      m_axisdir = NONE ;
      m_ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      m_ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      m_ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      m_ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
      m_ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
      m_ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );

      if ( rank == 0 )
      {
	cout << oshift << GrainsExec::m_shift3 << "Point3 A = " <<
		m_ptA[X] << " " << m_ptA[Y] << " " <<
		m_ptA[Z] << endl;
        cout << oshift << GrainsExec::m_shift3 << "Point3 B = " <<
		m_ptB[X] << " " << m_ptB[Y] << " " <<
		m_ptB[Z] << endl;
      }
      break;

    default:
      if ( rank == 0 ) cout << "Unknown insertion window "
		"type" << endl;
      ok = false;
      break;
  }
  
  return( ok );
}




// ----------------------------------------------------------------------------
// Returns a point randomly selected within thw window
Point3 Window::getInsertionPoint() const
{
  Point3 P;
  double r = 0., theta = 0., axiscoor = 0.;
  
  switch( m_ftype )
  {
    case WINDOW_BOX:
      P[X] = m_ptA[X]
	+ ( double(random()) / double(INT_MAX) )
	* ( m_ptB[X] - m_ptA[X] );
      P[Y] = m_ptA[Y]
	+ ( double(random()) / double(INT_MAX) )
	* ( m_ptB[Y] - m_ptA[Y] );
      P[Z] = m_ptA[Z]
	+ ( double(random()) / double(INT_MAX) )
	* ( m_ptB[Z] - m_ptA[Z] );
        break;

    case WINDOW_CYLINDER:
      r = ( double(random()) / double(INT_MAX) ) * m_radius;
      theta = ( double(random()) / double(INT_MAX) ) * 2. * PI;
      axiscoor = ( double(random()) / double(INT_MAX) ) * m_height;
      switch ( m_axisdir )
      {
        case X:
	  P[X] = m_ptA[X] + axiscoor;
	  P[Y] = m_ptA[Y] + r * cos( theta );
	  P[Z] = m_ptA[Z] + r * sin( theta );
	  break;

        case Y:
	  P[X] = m_ptA[X] + r * sin( theta );
	  P[Y] = m_ptA[Y] + axiscoor;
	  P[Z] = m_ptA[Z] + r * cos( theta );
	  break;

        default:
	  P[X] = m_ptA[X] + r * cos( theta );
	  P[Y] = m_ptA[Y] + r * sin( theta );
	  P[Z] = m_ptA[Z] + axiscoor;
	  break;
      }
      break;

    case WINDOW_ANNULUS:
      r = m_radius_int + ( double(random()) / double(INT_MAX) )
      		* ( m_radius - m_radius_int );
      theta = ( double(random()) / double(INT_MAX) ) * 2. * PI;
      axiscoor = ( double(random()) / double(INT_MAX) ) * m_height;
      switch ( m_axisdir )
      {
        case X:
	  P[X] = m_ptA[X] + axiscoor;
	  P[Y] = m_ptA[Y] + r * cos( theta );
	  P[Z] = m_ptA[Z] + r * sin( theta );
	  break;

        case Y:
	  P[X] = m_ptA[X] + r * sin( theta );
	  P[Y] = m_ptA[Y] + axiscoor;
	  P[Z] = m_ptA[Z] + r * cos( theta );
	  break;

        default:
	  P[X] = m_ptA[X] + r * cos( theta );
	  P[Y] = m_ptA[Y] + r * sin( theta );
	  P[Z] = m_ptA[Z] + axiscoor;
	  break;
      }
      break;

    case WINDOW_LINE:
      P = m_ptA + ( double(random()) / double(INT_MAX) ) * ( m_ptB - m_ptA );
      break;

    default:
      break;
  }

  return ( P );
}




// ----------------------------------------------------------------------------
// Sets the window as a box
void Window::setAsBox( Point3 const& ptA, Point3 const& ptB )
{
  m_ftype = WINDOW_BOX;
  m_ptA = ptA;
  m_ptB = ptB;
  m_radius = m_radius_int = m_height = 0.;
  m_axisdir = NONE;
}




// ----------------------------------------------------------------------------
// Returns a pointer to point A
Point3 const* Window::getPointA() const
{
  return ( &m_ptA );
}




// ----------------------------------------------------------------------------
// Returns a pointer to point B
Point3 const* Window::getPointB() const
{
  return ( &m_ptB );
}




// ----------------------------------------------------------------------------
// Adds the window as a rigid body to a list of rigid bodies
void Window::addAsRigidBody( list<RigidBody*>& iwlist ) const
{
  Transform gcwindow; 
  Convex* ccw = NULL;
  RigidBody* ffw = NULL;
  Vector3 v_axis;
  Matrix mrot;
  double Le = 0., Lh = 0., Lt = 0., angle = 0., 
  	mean_radius = 0., thickness = 0. ;
  Point3 panelCenter, center;
  size_t nbPanels = 32;
  double polygonAngle = 2. * PI / double(nbPanels);
    
  switch( m_ftype )
  {
    case WINDOW_BOX:
      gcwindow.setOrigin( 0.5 * ( m_ptA + m_ptB ) );
      ccw = new Box( 0.5 * ( m_ptB - m_ptA ) );
      ffw = new RigidBody( ccw, gcwindow );
      iwlist.push_back( ffw );
      break;
      
    case WINDOW_CYLINDER:
      switch ( m_axisdir )
      {
        case X: 
          v_axis[X] = m_height; 
          mrot.setValue( cos( 0.5 * PI ), -sin( 0.5 * PI ), 0.,
                sin( 0.5 * PI ), cos( 0.5 * PI ), 0.,
                0., 0., 1. );
          break;

        case Y: 
          v_axis[Y] = m_height; 
          break;
  
        default: 
          v_axis[Z] = m_height; 
          mrot.setValue( 1., 0., 0.,
                0., cos( 0.5 * PI ), -sin( 0.5 * PI ),
                0., sin( 0.5 * PI ), cos( 0.5 * PI ) );
          break;
      }
      gcwindow.setOrigin( m_ptA + 0.5 * v_axis );
      ccw = new Cylinder( m_radius, m_height );
      ffw = new RigidBody( ccw, gcwindow );
      ffw->getTransform()->setBasis( mrot );
      iwlist.push_back( ffw );
      break;	

    case WINDOW_ANNULUS:
      mean_radius = m_radius_int + 0.5 * ( m_radius - m_radius_int ) ;
      thickness = m_radius - m_radius_int	;
      center = m_ptA;
      Le = thickness;       
      Lt = m_radius * tan( polygonAngle ); 
      Lh = m_height;
      switch ( m_axisdir )
      {
        case X: 
          center[X] += 0.5 * m_height;
          for (size_t iNb=0; iNb!=nbPanels; iNb++)
          {
            angle = 2. * PI * double(iNb) / double(nbPanels);
            panelCenter[X] = center[X];
            panelCenter[Y] = center[Y] + mean_radius * cos(angle);
            panelCenter[Z] = center[Z] + mean_radius * sin(angle);
            gcwindow.setOrigin( panelCenter );
 
            mrot.setValue( 1., 0., 0.,
                  0., cos(angle), -sin(angle),
                  0., sin(angle), cos(angle) );
            gcwindow.setBasis( mrot );
            ccw = new Box( Lh, Le, Lt );
            ffw = new RigidBody( ccw, gcwindow );
            iwlist.push_back( ffw );
          }
          break;

        case Y: 
          center[Y] += 0.5 * m_height;
          for (size_t iNb=0; iNb!=nbPanels; iNb++)
          {
            angle = 2. * PI * double(iNb) / double(nbPanels);
            panelCenter[X] = center[X] + mean_radius * sin(angle);
            panelCenter[Y] = center[Y];
            panelCenter[Z] = center[Z] + mean_radius * cos(angle);
            gcwindow.setOrigin( panelCenter );

            mrot.setValue( cos(angle), 0., sin(angle),
                  0., 1., 0.,
                  -sin(angle), 0., cos(angle) );
            gcwindow.setBasis( mrot );
            ccw = new Box( Lt, Lh, Le );
            ffw = new RigidBody( ccw, gcwindow );
            iwlist.push_back( ffw );
          }
          break;

        default:
          center[Z] += 0.5 * m_height;
          for (size_t iNb=0; iNb!=nbPanels; iNb++)
          {
            angle = 2. * PI * double(iNb) / double(nbPanels);
            panelCenter[X] = center[X] + mean_radius * cos(angle);
            panelCenter[Y] = center[Y] + mean_radius * sin(angle);
            panelCenter[Z] = center[Z];
            gcwindow.setOrigin( panelCenter );

            mrot.setValue( cos(angle), -sin(angle), 0.,
                  sin(angle), cos(angle), 0.,
                  0., 0., 1. );
            gcwindow.setBasis( mrot );
            ccw = new Box( Le, Lt, Lh );
            ffw = new RigidBody( ccw, gcwindow );
            iwlist.push_back( ffw );
          }
          break;
      }
      break;

    case WINDOW_LINE:
      gcwindow = Segment::computeTransform( 
            0.5 * ( m_ptB - m_ptA ), 0.5 * ( m_ptA + m_ptB ) );
      ccw = new Segment( Norm( m_ptB - m_ptA ) );
      ffw = new RigidBody( ccw, gcwindow );
      iwlist.push_back( ffw );
      break;

    default:
      break;
  }
}




// ----------------------------------------------------------------------------
// Returns the cylinder radius
double Window::getRadius() const
{
  return ( m_radius );
}




// ----------------------------------------------------------------------------
// Returns the cylinder height
double Window::getHeight() const
{
  return ( m_height );
}




// ----------------------------------------------------------------------------
// Returns the cylinder type
WindowType Window::getType() const
{
  return ( m_ftype );
}




// ----------------------------------------------------------------------------
// Returns the cylinder axis direction */
Direction Window::getAxisDirection() const
{
  return ( m_axisdir );
}




// ----------------------------------------------------------------------------
// Shifts the 2 points defined the window by a specified amount in a specified 
// direction
void Window::shiftWindow( double const& geoshift, Direction dir )
{
  m_ptA[dir] += geoshift;
  m_ptB[dir] += geoshift;      
}




// ----------------------------------------------------------------------------
// Associates the imposed velocity to the insertion window
bool Window::LinkImposedMotion( ObstacleImposedVelocity* impvel )
{
  bool status = false;
  if ( m_name == impvel->getObstacleName() )
  {
    if ( m_kinematics == NULL ) m_kinematics = new ObstacleKinematicsVelocity();
    m_kinematics->append( impvel );
    status = true;
  }

  return ( status );
}




// ----------------------------------------------------------------------------
// Moves the window if it has an imposed motion
void Window::Move( double time, double dt )
{
  if ( m_kinematics )
    if ( m_kinematics->ImposedMotion( time, dt, m_ptA ) )
    {
      Vector3 translation = *(m_kinematics->getTranslation());
      m_ptA += translation;
      m_ptB += translation;
    }    
} 




// ----------------------------------------------------------------------------
// Resets kinematics to 0 */
void Window::resetKinematics()
{
  if ( m_kinematics ) m_kinematics->reset();  
}
