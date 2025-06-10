#include "TruncatedCone.hh"
#include "GrainsExec.hh"

int TruncatedCone::m_visuNodeNbOnPer = 32;


// ----------------------------------------------------------------------------
// Constructor with radius and height as input parameters
TruncatedCone::TruncatedCone( double r_bot, double r_top, double h )
  : m_bottomRadius( r_bot )
  , m_topRadius( r_top )
  , m_halfHeight( h / 2. )
  , m_bottomHeight( ( h / 4. ) * ( r_bot * r_bot + 2. * r_bot * r_top 
  	+ 3. * r_top * r_top ) 
		/ ( r_bot * r_bot + r_bot * r_top + r_top * r_top ) )
  , m_topHeight( h - m_bottomHeight )
  , m_Hprim( h * m_bottomRadius / ( m_bottomRadius - m_topRadius ) )
  , m_sinAngle( m_bottomRadius / sqrt( m_bottomRadius * m_bottomRadius 
  	+ m_Hprim * m_Hprim ) )
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
TruncatedCone::TruncatedCone( istream &fileIn )
{
  readShape( fileIn );
}




// -----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
TruncatedCone::TruncatedCone( DOMNode* root )
{
  m_bottomRadius = ReaderXML::getNodeAttr_Double( root, "BottomRadius" );
  m_topRadius = ReaderXML::getNodeAttr_Double( root, "TopRadius" );
  m_halfHeight = ReaderXML::getNodeAttr_Double( root, "Height") / 2.;
  m_bottomHeight = ( 2. * m_halfHeight / 4. ) 
  	* ( m_bottomRadius * m_bottomRadius 
  	+ 2. * m_bottomRadius * m_topRadius + 3. * m_topRadius * m_topRadius )
	/ ( m_bottomRadius * m_bottomRadius + m_bottomRadius * m_topRadius 
		+ m_topRadius * m_topRadius );
  m_topHeight = ( 2. * m_halfHeight ) - m_bottomHeight;
  m_Hprim = 2. * m_halfHeight * m_bottomRadius 
  	/ ( m_bottomRadius - m_topRadius );
  m_sinAngle = m_bottomRadius / sqrt( m_bottomRadius * m_bottomRadius 
  	+ m_Hprim * m_Hprim );
}




// ----------------------------------------------------------------------------
// Destructor
TruncatedCone::~TruncatedCone()
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType TruncatedCone::getConvexType() const
{
  return ( TRUNCATEDCONE );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool TruncatedCone::BuildInertia( double* inertia, double* inertia_1 ) const
{
  double H = 2. * m_halfHeight;
  double V1 = ( m_Hprim * PI / 3.0 ) * ( m_bottomRadius * m_bottomRadius );
  double V2 = ( ( m_Hprim - H ) * PI / 3.0 ) * ( m_topRadius * m_topRadius );
  
  inertia[1] = inertia[2] = inertia[4] = 0.0;  
  inertia[0] = inertia[5] = ( 3. / 80. ) 
  	* V1 * ( 4. * m_bottomRadius * m_bottomRadius + m_Hprim * m_Hprim )
	+ V1 * ( m_Hprim / 4. - m_bottomHeight ) 
		* ( m_Hprim / 4. - m_bottomHeight )
	- ( 3. / 80. ) * V2 * ( 4. * m_topRadius * m_topRadius 
		+ ( m_Hprim - H ) * ( m_Hprim - H ) ) 
	- V2 * ( ( m_Hprim - H ) / 4. + m_topHeight )
		* ( ( m_Hprim - H ) /4. + m_topHeight );    
  inertia[3] = ( 3. / 10. ) * ( V1 * m_bottomRadius * m_bottomRadius 
  	- V2 * m_topRadius * m_topRadius );

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[5] = inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  return ( true );
  
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference sphere,
// i.e., without applying any transformation
double TruncatedCone::computeCircumscribedRadius() const
{
  return ( max( sqrt( m_bottomRadius * m_bottomRadius 
  		+ m_bottomHeight * m_bottomHeight ),
  	sqrt( m_topRadius * m_topRadius + m_topHeight * m_topHeight ) ) );
}




// ----------------------------------------------------------------------------
// Returns a clone of the truncated cone
Convex* TruncatedCone::clone() const
{
  return ( new TruncatedCone( m_bottomRadius, m_topRadius, 
  	2. * m_halfHeight ) );
}




// ----------------------------------------------------------------------------
// Returns the TruncatedCone volume
double TruncatedCone::getVolume() const
{
  return ( ( 2. * m_halfHeight * PI / 3.0 ) 
  	* ( m_bottomRadius * m_bottomRadius + m_topRadius * m_topRadius 
		+ m_bottomRadius * m_topRadius ) );
}




// ----------------------------------------------------------------------------
// Truncated cone support function, returns the support point P, i.e. the
// point on the surface of the truncated cone that satisfies max(P.v)
Point3 TruncatedCone::support( Vector3 const& v ) const
{
  double norm = Norm( v );
  if ( norm > EPSILON )
  {
    double s = sqrt( v[X] * v[X] + v[Z] * v[Z] );
    if ( s > EPSILON )
    {      
      if (v[Y] > norm * m_sinAngle)
      {
        double d = m_topRadius / s;
        return ( Point3( v[X] * d,   m_topHeight, v[Z] * d ) );
      }
      else
      {
        double d = m_bottomRadius / s;
        return ( Point3( v[X] * d, - m_bottomHeight, v[Z] * d ) );
      }
    }
    else
      return ( Point3( 0., v[Y] < 0. ? - m_bottomHeight : m_topHeight, 0. ) );
  }
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the surface of the
// truncated cone. Here simply returns 4 points as follows: center of bottom 
// circular face, an arbitrary point on the bottom perimeter, center of top 
// circular face and an arbitrary point on the top perimeter
vector<Point3> TruncatedCone::getEnvelope() const
{
  Point3 point(0.,0.,0.);
  vector<Point3> surface( 4, point );
  surface[0][Y] = - m_bottomHeight;
  surface[1][Y] = - m_bottomHeight;
  surface[1][X] = m_bottomRadius;
  surface[2][Y] = m_topHeight;
  surface[3][Y] = m_topHeight;
  surface[3][X] = m_topRadius;  
  return ( surface );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 8888
int TruncatedCone::getNbCorners() const
{
  return ( 8888 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* TruncatedCone::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Output operator
void TruncatedCone::writeShape( ostream& fileOut ) const
{
  fileOut << "*TruncatedCone " << m_bottomRadius << " " <<  m_topRadius << " " 
  	<< 2.0 * m_halfHeight << " *END";
}



// ----------------------------------------------------------------------------
// Input operator
void TruncatedCone::readShape( istream& fileIn )
{
  fileIn >> m_bottomRadius >> m_topRadius >> m_halfHeight;
  m_halfHeight /= 2.0;  
  m_bottomHeight = ( m_halfHeight / 2. ) * ( m_bottomRadius * m_bottomRadius 
  	+ 2. * m_bottomRadius * m_topRadius + 3. * m_topRadius * m_topRadius )
	/( m_bottomRadius * m_bottomRadius + m_bottomRadius * m_topRadius 
		+ m_topRadius * m_topRadius );
  m_topHeight = 2. * m_halfHeight - m_bottomHeight;
  m_Hprim = 2. * m_halfHeight * m_bottomRadius 
  	/ ( m_bottomRadius - m_topRadius );
  m_sinAngle = m_bottomRadius / sqrt( m_bottomRadius * m_bottomRadius 
  	+ m_Hprim * m_Hprim );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the truncated cone in a Paraview format
int TruncatedCone::numberOfPoints_PARAVIEW() const
{
  return ( 2 * m_visuNodeNbOnPer + 2 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the truncated cone in a
// Paraview format
int TruncatedCone::numberOfCells_PARAVIEW() const
{
  return ( m_visuNodeNbOnPer );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the truncated cone in a Paraview format
void TruncatedCone::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp, p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Lower disk rim
  p[Y] = - m_bottomHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_bottomRadius * cos ( i * dtheta );
    p[Z] = m_bottomRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Upper disk rim
  p[Y] = m_topHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_topRadius * cos ( i * dtheta );
    p[Z] = m_topRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Lower disk center
  p[X] = 0.;
  p[Y] = - m_bottomHeight;
  p[Z] = 0.;
  pp = transform( p );

  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

  // Upper disk center
  p[Y] = m_topHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the truncated cone in a Paraview format
list<Point3> TruncatedCone::get_polygonsPts_PARAVIEW( 
	Transform const& transform, Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp,p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Lower disk rim
  p[Y] = - m_bottomHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_bottomRadius * cos ( i * dtheta );
    p[Z] = m_bottomRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Upper disk rim
  p[Y] = m_topHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_topRadius * cos ( i * dtheta );
    p[Z] = m_topRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Lower disk center
  p[X] = 0.;
  p[Y] = - m_bottomHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  // Upper disk center
  p[Y] = m_topHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the truncated cone in a Paraview format
void TruncatedCone::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer);
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer
    	+ 1 );
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer
    	+ 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer + 1 );
  last_offset += 6;
  offsets.push_back( last_offset );
  cellstype.push_back( 13 );

  firstpoint_globalnumber += 2 * m_visuNodeNbOnPer + 2;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the truncated cone
bool TruncatedCone::isIn( Point3 const& pt ) const
{    
  return ( pt[Y] > - m_bottomHeight && pt[Y] < m_topHeight 
  	&& sqrt( pt[X] * pt[X] + pt[Z] * pt[Z] ) < 
		m_bottomRadius - ( m_bottomRadius - m_topRadius ) 
		* ( m_bottomHeight + pt[Y] ) / ( 2. * m_halfHeight ));
    
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to truncated cone
BVolume* TruncatedCone::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( Vector3( max( m_bottomRadius, m_topRadius ),
	max( m_bottomHeight, m_topHeight ), 
	max( m_bottomRadius, m_topRadius ) ), Matrix() );
  else if ( type == 2 ) // OBC
    bvol = new OBC( max( m_bottomRadius, m_topRadius ), 
    	2. * max( m_bottomHeight, m_topHeight ), Vector3(0., 1., 0.) );

  return( bvol );
}



// ----------------------------------------------------------------------------
// Sets the number of point over the truncated cone perimeter for Paraview 
// post-processing, i.e., controls the number of facets in the truncated cone 
// reconstruction in Paraview
void TruncatedCone::SetvisuNodeNbOverPer( int nbpts )
{
  m_visuNodeNbOnPer = nbpts;
}




// ----------------------------------------------------------------------------
// Performs advanced comparison of the two truncated cones and returns whether 
// they match
bool TruncatedCone::equalType_level2( Convex const* other ) const
{
  // We know that other points to a TruncatedCone, we dynamically cast it to 
  // actual type
  TruncatedCone const* other_ = dynamic_cast<TruncatedCone const*>(other);
  
  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() );    

  bool same = ( 
  	fabs( m_bottomRadius - other_->m_bottomRadius ) <  LOWEPS * lmin 
	&& fabs( m_topRadius - other_->m_topRadius ) <  LOWEPS * lmin 
	&& fabs( m_halfHeight - other_->m_halfHeight ) <  LOWEPS * lmin
	&& fabs( m_bottomHeight - other_->m_bottomHeight ) <  LOWEPS * lmin	
	&& fabs( m_topHeight - other_->m_topHeight ) <  LOWEPS * lmin	
	&& fabs( m_sinAngle - other_->m_sinAngle ) <  LOWEPS );
  
  return ( same );
}




// ----------------------------------------------------------------------------
// Writes the truncated cone in an OBJ format
void TruncatedCone::write_convex_OBJ( ostream& f, Transform  const& transform,
    	size_t& firstpoint_number ) const
{
  Point3 pp, p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Vertices  
  // Bottom disk rim
  p[Y] = - m_bottomHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_bottomRadius * cos ( i * dtheta );
    p[Z] = m_bottomRadius * sin ( i * dtheta );
    pp = transform( p );
    f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Z] ) << " " << endl;	
  }

  // Bottom disk center
  p[X] = 0.;
  p[Y] = - m_bottomHeight;
  p[Z] = 0.;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;

  // Top disk rim
  p[Y] = m_topHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_topRadius * cos ( i * dtheta );
    p[Z] = m_topRadius * sin ( i * dtheta );
    pp = transform( p );
    f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Z] ) << " " << endl;	
  }

  // Top disk center
  p[X] = 0.;
  p[Y] = m_topHeight;
  p[Z] = 0.;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;
  
  // Faces 
  // Triangular lateral faces 
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    f << "f "<< firstpoint_number + i << " "
    	<< firstpoint_number + i + 1 << " "
    	<< firstpoint_number + m_visuNodeNbOnPer + i + 1 << endl;
    f << "f "<< firstpoint_number + i + 1 << " "
    	<< firstpoint_number + m_visuNodeNbOnPer + i + 2 << " "
    	<< firstpoint_number + m_visuNodeNbOnPer + i + 1 << endl;	
  }
  f << "f " << firstpoint_number + m_visuNodeNbOnPer - 1 << " "
  	<< firstpoint_number << " "
  	<< firstpoint_number + 2 * m_visuNodeNbOnPer << endl;
  f << "f "<< firstpoint_number << " "
    	<< firstpoint_number + m_visuNodeNbOnPer + 1 << " "
    	<< firstpoint_number + 2 * m_visuNodeNbOnPer << endl;

  // Triangular bottom faces 
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    f << "f "<< firstpoint_number + i << " "
    	<< firstpoint_number + i + 1 << " "
    	<< firstpoint_number + m_visuNodeNbOnPer << endl;
  }
  f << "f " << firstpoint_number + m_visuNodeNbOnPer - 1 << " "
  	<< firstpoint_number << " "
  	<< firstpoint_number + m_visuNodeNbOnPer << endl;
	
  // Triangular top faces 
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    f << "f "<< firstpoint_number + i + m_visuNodeNbOnPer + 1 << " "
    	<< firstpoint_number + i + m_visuNodeNbOnPer + 2 << " "
    	<< firstpoint_number + 2 * m_visuNodeNbOnPer + 1 << endl;
  }
  f << "f "<< firstpoint_number + 2 * m_visuNodeNbOnPer << " "
    	<< firstpoint_number + m_visuNodeNbOnPer + 1 << " "
    	<< firstpoint_number + 2 * m_visuNodeNbOnPer + 1 << endl;	

  firstpoint_number += 2 * m_visuNodeNbOnPer + 2;    			  
}

