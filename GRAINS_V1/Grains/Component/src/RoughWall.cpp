#include "RoughWall.hh"
#include "GrainsExec.hh"
#include "PointC.hh"
#include "ContactBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "Box.hh"
#include "Sphere.hh"
#include "LinkedCell.hh"


// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
RoughWall::RoughWall( DOMNode* root ) :
  CompositeObstacle( "roughwall" )
{    
  assert( root != NULL );
  m_nb_spheres = new size_t[2];
  
  // Type
  m_type = "RoughWall";  

  // Name
  m_name = ReaderXML::getNodeAttr_String( root, "name" );

  // Position of the rough wall
  m_geoRBWC->getTransform()->load( root );

  // Material
  DOMNode* materiau_ = ReaderXML::getNode( root, "Material" );
  m_materialName = ReaderXML::getNodeValue_String( materiau_ );
  ContactBuilderFactory::defineMaterial( m_materialName, true );

  // Obstacle to transfer to the fluid
  bool transferToFluid = false;
  DOMNode* status = ReaderXML::getNode( root, "Status" );
  if ( status )
    transferToFluid = ReaderXML::getNodeAttr_Int( status, "ToFluid" );
 
  // Geometry
  DOMNode* nGeom = ReaderXML::getNode( root, "Geometry" );
  
  DOMNode* nRough = ReaderXML::getNode( nGeom, "Roughness" );  
  m_radius = ReaderXML::getNodeAttr_Double( nRough, "Radius" );
  m_shift = ReaderXML::getNodeAttr_Double( nRough, "Shift" );
  m_nb_spheres[0] = size_t(ReaderXML::getNodeAttr_Int( nRough, "NX" ));
  m_nb_spheres[1] = size_t(ReaderXML::getNodeAttr_Int( nRough, "NY" ));
  double crust_thickness = ReaderXML::getNodeAttr_Double( nRough, 
  	"CrustThickness" );
  m_random_mag = ReaderXML::getNodeAttr_Double( nRough, "RandMag" );
  
  m_periodic = false;
  m_restrict_box_geommotion = false;
  DOMNode* nPer = ReaderXML::getNode( nGeom, "Periodicity" ); 
  if ( nPer )
  {
    m_periodic = true; 
    m_periodic_ext = ReaderXML::getNodeAttr_Double( nPer, "Extension" );
    string sper =  ReaderXML::getNodeAttr_String( nPer, "Directions" );
    if ( sper == "X" ) m_periodic_directions.push_back( 0 );
    else if ( sper == "Y" ) m_periodic_directions.push_back( 1 );
    else if ( sper == "Z" ) m_periodic_directions.push_back( 2 );  
    else if ( sper == "XY" || sper == "YX" ) 
    {
      m_periodic_directions.push_back( 0 ); 
      m_periodic_directions.push_back( 1 );
    }       
    else if ( sper == "XZ" || sper == "ZX" ) 
    {
      m_periodic_directions.push_back( 0 ); 
      m_periodic_directions.push_back( 2 );
    }  
    else if ( sper == "YZ" || sper == "ZY" ) 
    {
      m_periodic_directions.push_back( 1 ); 
      m_periodic_directions.push_back( 2 );
    }
    else
      cout << "WARNING: unknown periodic directions in RoughWall " <<
     	 m_name << endl;
    if ( ReaderXML::hasNodeAttr( nPer, "RestrictBoxMotion" ) )
    {	 
      string sres = ReaderXML::getNodeAttr_String( nPer, "RestrictBoxMotion" );
      if ( sres == "True" ) m_restrict_box_geommotion = true;
    }
  }   
  
  // Crust thickness 
  m_geoRBWC->setCrustThickness( crust_thickness ); 
  
  // Create the box obstacle (smooth wall)
  DOMNode* nBox = ReaderXML::getNode( nGeom, "Box" );
  double lx = ReaderXML::getNodeAttr_Double( nBox, "LX" );
  double ly = ReaderXML::getNodeAttr_Double( nBox, "LY" );  
  double lz = ReaderXML::getNodeAttr_Double( nBox, "LZ" );    
  string name = m_name + "_Box";	
  Box* box = new Box( lx, ly, lz );
  RigidBodyWithCrust* geoRBWC_box = new RigidBodyWithCrust( box, Transform(),
	false, crust_thickness );
  Obstacle* sbox = new SimpleObstacle( name, geoRBWC_box, 
 	m_materialName, transferToFluid, true );
  Point3 cg, cgphys;
  sbox->setPosition( cg );
  if ( m_restrict_box_geommotion )
    sbox->setRestrictedGeomDirMotion( m_periodic_directions );
  m_obstacles.push_back( sbox );		 

  // Create elememtary sphere obstacles to create roughness
  int counter = 1;
  double dx = lx / double( m_nb_spheres[0] );
  double dy = ly / double( m_nb_spheres[1] ); 
  bool create = true; 
  for (size_t i=0;i<m_nb_spheres[0];++i)  
    for (size_t j=0;j<m_nb_spheres[1];++j)
    {
      // Position
      cg.setValue( - 0.5 * lx + ( 0.5 + double(i) ) * dx
      	+ ( 2. * (double)rand() / RAND_MAX - 1. ) * m_random_mag, 
	- 0.5 * ly + ( 0.5 + double(j) ) * dy
	+ ( 2. * (double)rand() / RAND_MAX - 1. ) * m_random_mag, 
	m_shift + ( 2. * (double)rand() / RAND_MAX - 1. ) * m_random_mag );
	
      // If periodic, check if in the global domain
      create = true;
      if ( GrainsExec::m_periodic )
      {
        cgphys = (*m_geoRBWC->getTransform())(cg);
	if ( !isInDomain( &cgphys ) ) create = false;
      }
    
      // Create sphere
      if ( create )
      {
        Sphere* sphere = new Sphere( m_radius );
        RigidBodyWithCrust* geoRBWC_sphere = new RigidBodyWithCrust( sphere, 
      		Transform(), false, crust_thickness ); 
        Obstacle* ssphere = new SimpleObstacle( "", geoRBWC_sphere, 
		m_materialName, transferToFluid, true );
        ssphere->setPosition( cg );
	name = m_name + "_Rough" + GrainsExec::intToString( ssphere->getID() );
	ssphere->setName( name );
        m_obstacles.push_back( ssphere );
	++counter;
      }	    
    }
    
  // Apply transform to all simple obstacles
  Vector3 zerotranslation;
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end();obstacle++)
    (*obstacle)->getRigidBody()->composeLeftByTransform( 
    	*m_geoRBWC->getTransform() ); 
  
  // Create periodic roughness  
  if ( GrainsExec::m_periodic && m_periodic )
  {
    App::get_origin( m_domain_origin[X], m_domain_origin[Y], 
    	m_domain_origin[Z] );
    App::get_size( m_domain_size[X], m_domain_size[Y], m_domain_size[Z] );
    list<Obstacle*> perspheres;
    for (size_t m=0;m<3;++m) 
      m_lper[m] = min( m_periodic_ext, m_domain_size[m] );
    size_t m, n;
    
    // Mono-periodicity
    list<size_t>::const_iterator id;
    for (id=m_periodic_directions.cbegin();id!=m_periodic_directions.cend();
    	id++)
    {
      obstacle = m_obstacles.begin();
      obstacle++;
      m = *id;
      for (; obstacle!=m_obstacles.end();obstacle++)
      {
	cg = *(*obstacle)->getPosition();
	if ( ( cg[m] > m_domain_origin[m] 
	  		&& cg[m] < m_domain_origin[m] + m_lper[m] )
		|| ( cg[m] > m_domain_origin[m] + m_domain_size[m] - m_lper[m]
			&& cg[m] < m_domain_origin[m] + m_domain_size[m] ) )
	{
          Sphere* sphere = new Sphere( m_radius );
          RigidBodyWithCrust* geoRBWC_sphere = new RigidBodyWithCrust( sphere,
      		Transform(), false, crust_thickness ); 
          Obstacle* ssphere = new SimpleObstacle( "", geoRBWC_sphere, 
		m_materialName, transferToFluid, false );
	  ssphere->setID( (*obstacle)->getID() );
	  cg[m] = cg[m] + ( cg[m] < m_domain_origin[m] + m_lper[m] ? 
	    	1. : -1. ) * m_domain_size[m];
          ssphere->setPosition( cg );
          name = m_name + "_Rough" + GrainsExec::intToString( ssphere->getID() );
	  ssphere->setName( name );	  
          perspheres.push_back( ssphere );
	  ++counter;	    
	}
      }    
    }
    
    // Bi-periodicity
    if ( m_periodic_directions.size() == 2 )
    {
      m = m_periodic_directions.front(); 
      n = m_periodic_directions.back();

      obstacle = m_obstacles.begin();
      obstacle++;
      for (; obstacle!=m_obstacles.end();obstacle++)
      {
	cg = *(*obstacle)->getPosition();
	create = false;
	if ( ( cg[m] > m_domain_origin[m] 
	  		&& cg[m] < m_domain_origin[m] + m_lper[m] )
	  	&& ( cg[n] > m_domain_origin[n] 
			&& cg[n] < m_domain_origin[n] + m_lper[n] ) )
	{
	  create = true;
	  cg[m] += m_domain_size[m];
	  cg[n] += m_domain_size[n];	      
	}	  

	if ( !create )
	  if ( ( cg[m] > m_domain_origin[m] 
	    		&& cg[m] < m_domain_origin[m] + m_lper[m] )
	  	&& ( cg[n] > m_domain_origin[n] + m_domain_size[n] - m_lper[n]
			&& cg[n] < m_domain_origin[n] + m_domain_size[n] ) )
	  {
	    create = true;
	    cg[m] += m_domain_size[m];
	    cg[n] -= m_domain_size[n];	      
	  }
	  
	if ( !create )
	  if ( ( cg[m] > m_domain_origin[m] + m_domain_size[m] - m_lper[m]
			&& cg[m] < m_domain_origin[m] + m_domain_size[m] )
	  	&& ( cg[n] > m_domain_origin[n] 
			&& cg[n] < m_domain_origin[n] + m_lper[n] ) )
	  {
	    create = true;
	    cg[m] -= m_domain_size[m];
	    cg[n] += m_domain_size[n];	      
	  }	  	  
	   
	if ( !create )
	  if ( ( cg[m] > m_domain_origin[m] + m_domain_size[m] - m_lper[m]
			&& cg[m] < m_domain_origin[m] + m_domain_size[m] )
	  	&& ( cg[n] > m_domain_origin[n] + m_domain_size[n] - m_lper[n]
			&& cg[n] < m_domain_origin[n] + m_domain_size[n] ) )
	  {
	    create = true;
	    cg[m] -= m_domain_size[m];
	    cg[n] -= m_domain_size[n];	      
	  }	   

        if ( create )
	{
          Sphere* sphere = new Sphere( m_radius );
          RigidBodyWithCrust* geoRBWC_sphere = new RigidBodyWithCrust( sphere,
      		Transform(), false, crust_thickness ); 
          Obstacle* ssphere = new SimpleObstacle( "", geoRBWC_sphere, 
		m_materialName, transferToFluid, false );
	  ssphere->setID( (*obstacle)->getID() );
          ssphere->setPosition( cg );
          name = m_name + "_Rough" + GrainsExec::intToString( ssphere->getID() );
	  ssphere->setName( name );	  
          perspheres.push_back( ssphere );
	  ++counter;	    
	}   
      } 
    }
                     
    // Add periodic spheres
    for (obstacle=perspheres.begin(); obstacle!=perspheres.end();obstacle++)
      m_obstacles.push_back( *obstacle );
      
    // Initialize position at previous time
    obstacle = m_obstacles.begin();
    obstacle++;
    for (; obstacle!=m_obstacles.end();obstacle++)
      if ( isInDomain( (*obstacle)->getPosition() ) ) 
        m_pos_nm1[(*obstacle)->getID()] = *(*obstacle)->getPosition();
  } 
      
  // Update simple obstacle bounding boxes  
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end();obstacle++)
    (*obstacle)->Translate( zerotranslation );      	
}




// ----------------------------------------------------------------------------
// Constructor with name as input parameter
RoughWall::RoughWall( string const& s )
  : CompositeObstacle( s )
{
  m_type = "RoughWall";
  m_nb_spheres = new size_t[2]; 
} 




// ----------------------------------------------------------------------------
// Destructor
RoughWall::~RoughWall()
{
  delete [] m_nb_spheres;
}




// ----------------------------------------------------------------------------
// Checks if there is anything special to do about periodicity and
// if there is applies periodicity
void RoughWall::periodicity( LinkedCell* LC )
{
  Point3 cg;
  list<Obstacle*>::iterator obstacle;
  string name;
  list<Obstacle*> perspheres;
  bool create = true;
  size_t m, n, ndel = 0; 
  
  if ( GrainsExec::m_periodic && m_periodic && m_ismoving )
  {
    // Delete spheres out of extended periodic domain
    list<size_t>::const_iterator id;
    for (id=m_periodic_directions.cbegin();id!=m_periodic_directions.cend();
    	id++)
    {
      obstacle = m_obstacles.begin();
      obstacle++;
      m = *id;
      for (; obstacle!=m_obstacles.end();)
      {
	cg = *(*obstacle)->getPosition();
	if ( cg[m] < m_domain_origin[m] - m_lper[m] 
		|| cg[m] > m_domain_origin[m] + m_domain_size[m] + m_lper[m] )
	{
	  LC->remove( (SimpleObstacle*)*obstacle );
	  delete *obstacle;
	  obstacle = m_obstacles.erase( obstacle );
	  ndel++;
	}
	else obstacle++;
      }
    }
    
    if ( ndel ) LC->resetListSimpleObstacles();
   
    // Create new periodic spheres
    // Mono-periodicity
    for (id=m_periodic_directions.cbegin();id!=m_periodic_directions.cend();
    	id++)
    {    
      obstacle = m_obstacles.begin();
      obstacle++;
      m = *id;
      for (; obstacle!=m_obstacles.end();obstacle++)
      {
	if ( isInDomain( (*obstacle)->getPosition() ) )
	{
	  create = false;
	  cg = *(*obstacle)->getPosition();

	  if ( cg[m] > m_domain_origin[m] 
		&& cg[m] < m_domain_origin[m] + m_lper[m] 
		&& m_pos_nm1[(*obstacle)->getID()][m] 
		> m_domain_origin[m] + m_lper[m]
		&& fabs( cg[m] - m_pos_nm1[(*obstacle)->getID()][m] ) 
			< m_radius )
          {
	    create = true;
	    cg[m] += m_domain_size[m];	      
	  }			
		
	  if ( !create )
	    if ( cg[m] > m_domain_origin[m] + m_domain_size[m] - m_lper[m]
		&& cg[m] < m_domain_origin[m] + m_domain_size[m] 
		&& m_pos_nm1[(*obstacle)->getID()][m] 
		< m_domain_origin[m] + m_domain_size[m] - m_lper[m] 
		&& fabs( cg[m] - m_pos_nm1[(*obstacle)->getID()][m] ) 
			< m_radius )
            {
	      create = true;
	      cg[m] -= m_domain_size[m];	      
	    }	    
	    
          if ( create )	    	    
	  {
            Sphere* sphere = new Sphere( m_radius );
            RigidBodyWithCrust* geoRBWC_sphere = new RigidBodyWithCrust( 
	      		sphere, Transform(), false, 
			(*obstacle)->getCrustThickness() ); 
            Obstacle* ssphere = new SimpleObstacle( "", geoRBWC_sphere, 
			(*obstacle)->getMaterial(), 
			(*obstacle)->getObstaclesToFluid().size() != 0, false );
	    ssphere->setID( (*obstacle)->getID() );	      
            ssphere->setPosition( cg );
	    ssphere->setName( (*obstacle)->getName() );
            perspheres.push_back( ssphere );
	  }
        }
      }    
    }
    
    // Bi-periodicity    
    if ( m_periodic_directions.size() == 2 )
    {
      m = m_periodic_directions.front(); 
      n = m_periodic_directions.back();    
      
      obstacle = m_obstacles.begin();
      obstacle++;
      for (; obstacle!=m_obstacles.end();obstacle++)
      {
	if ( isInDomain( (*obstacle)->getPosition() ) )
	{
	  cg = *(*obstacle)->getPosition();
	  create = false;
	  
	  if ( fabs( cg[m] - m_pos_nm1[(*obstacle)->getID()][m] ) 
			< m_radius
		&& fabs( cg[n] - m_pos_nm1[(*obstacle)->getID()][n] ) 
			< m_radius )
	  {
	    if ( ( cg[m] > m_domain_origin[m] 
			&& cg[m] < m_domain_origin[m] + m_lper[m] )
	  	&& ( cg[n] > m_domain_origin[n] 
			&& cg[n] < m_domain_origin[n] + m_lper[n] )
		&& ( m_pos_nm1[(*obstacle)->getID()][m] 
			> m_domain_origin[m] + m_lper[m] 
			|| m_pos_nm1[(*obstacle)->getID()][n] 
			> m_domain_origin[n] + m_lper[n] ) )			
	    {
	      create = true;
	      cg[m] += m_domain_size[m];
	      cg[n] += m_domain_size[n];	      
	    }	  

            if ( !create )
	      if ( ( cg[m] > m_domain_origin[m] 
	    		&& cg[m] < m_domain_origin[m] + m_lper[m] )
	  	&& ( cg[n] > m_domain_origin[n] + m_domain_size[n] - m_lper[n]
			&& cg[n] < m_domain_origin[n] + m_domain_size[n] ) 
		&& ( m_pos_nm1[(*obstacle)->getID()][m] 
			> m_domain_origin[m] + m_lper[m] 
			|| m_pos_nm1[(*obstacle)->getID()][n] 
			< m_domain_origin[n] + m_domain_size[n] - m_lper[n] ) )
	      {
	        create = true;
	        cg[m] += m_domain_size[m];
	        cg[n] -= m_domain_size[n];	      
	      }
	  
            if ( !create )
	      if ( ( cg[m] > m_domain_origin[m] + m_domain_size[m] - m_lper[m]
			&& cg[m] < m_domain_origin[m] + m_domain_size[m] )
	  	&& ( cg[n] > m_domain_origin[n] 
			&& cg[n] < m_domain_origin[n] + m_lper[n] ) 
		&& ( m_pos_nm1[(*obstacle)->getID()][m] 
			< m_domain_origin[m] + m_domain_size[m] - m_lper[m] 
			|| m_pos_nm1[(*obstacle)->getID()][n] 
			> m_domain_origin[n] + m_lper[n] ) )
	      {
	        create = true;
	        cg[m] -= m_domain_size[m];
	        cg[n] += m_domain_size[n];	      
	      }	  	  
	   
	    if ( !create )
	      if ( ( cg[m] > m_domain_origin[m] + m_domain_size[m] 
		  	- m_lper[m] && cg[m] < m_domain_origin[m] 
			+ m_domain_size[m] )
	  	&& ( cg[n] > m_domain_origin[n] + m_domain_size[n] - m_lper[n]
			&& cg[n] < m_domain_origin[n] + m_domain_size[n] ) 
		&& ( m_pos_nm1[(*obstacle)->getID()][m] 
			< m_domain_origin[m] + m_domain_size[m] - m_lper[m] 
			|| m_pos_nm1[(*obstacle)->getID()][n] 
			< m_domain_origin[n] + m_domain_size[n] - m_lper[n] ))
	      {
	        create = true;
	        cg[m] -= m_domain_size[m];
	        cg[n] -= m_domain_size[n];	      
	      }	   

            if ( create )
	    {
              Sphere* sphere = new Sphere( m_radius );
              RigidBodyWithCrust* geoRBWC_sphere = new RigidBodyWithCrust( 
		  	sphere, Transform(), false, 
			(*obstacle)->getCrustThickness() ); 
              Obstacle* ssphere = new SimpleObstacle( "", geoRBWC_sphere, 
			(*obstacle)->getMaterial(), 
			(*obstacle)->getObstaclesToFluid().size() != 0, false );
	      ssphere->setID( (*obstacle)->getID() );
              ssphere->setPosition( cg );
	      ssphere->setName( (*obstacle)->getName() );
              perspheres.push_back( ssphere );	    
	    }
	  }
        }
      }
    }
      
    // Add periodic spheres
    for (obstacle=perspheres.begin(); obstacle!=perspheres.end();obstacle++)
      m_obstacles.push_back( *obstacle );      
            
    // Update position at previous time
    obstacle = m_obstacles.begin();
    obstacle++;
    for (; obstacle!=m_obstacles.end();obstacle++)
      if ( isInDomain( (*obstacle)->getPosition() ) ) 
        m_pos_nm1[(*obstacle)->getID()] = *(*obstacle)->getPosition();    
  }
}    




// ----------------------------------------------------------------------------
// Returns whether a position is in the global domain in the directions 
// specified by m_periodic_directions
bool RoughWall::isInDomain( Point3 const* pos )
{
  bool isIn = false;
  
  switch ( m_periodic_directions.size() )
  {
    case 1:
      isIn = App::isInDomain( pos, m_periodic_directions.front() );
      break;
      
    default: // =2 here
    isIn = App::isInDomain( pos, m_periodic_directions.front() )
    	&& App::isInDomain( pos, m_periodic_directions.back() );
      break;
  }
  
  return ( isIn );
}




// ----------------------------------------------------------------------------
// Outputs the cylindrical shell for reload
void RoughWall::write( ostream& fileSave ) const
{
  fileSave << "<Composite> " << m_name << " " << m_type << endl;
  fileSave << "*Properties " << m_shift << " " << m_nb_spheres[0] << " " 
  	<< m_nb_spheres[1] << " " << m_radius << " " << m_random_mag << " "
	<< m_periodic_ext << " " << m_periodic << " "
	<< m_periodic_directions.size() << " ";
  if ( m_periodic_directions.size() )
    for (list<size_t>::const_iterator il=m_periodic_directions.begin();
    	il!=m_periodic_directions.end();il++) fileSave << *il << " ";
  fileSave << m_restrict_box_geommotion << " " << m_domain_origin << " " << 
  	m_domain_size << " " << m_lper << endl;
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
// Reloads the rough wall and links it to the higher level 
// obstacle in the obstacle tree
void RoughWall::reload( Obstacle& mother, istream& file )
{
  string ttag, buffer;
  size_t size, dir;
  
  // Read extra properties  
  file >> buffer >> m_shift >> m_nb_spheres[0] >> m_nb_spheres[1] >> m_radius >>
  	m_random_mag >> m_periodic_ext >> m_periodic >> size;
  if  ( size )
    for (size_t i=0;i<size;++i) 
    {
      file >> dir;
      m_periodic_directions.push_back( dir );
    }
  file >> m_restrict_box_geommotion >> m_domain_origin >> m_domain_size >> 
  	m_lper;
  
  // Standard Composite Obstacle reload  
  if ( m_CompositeObstacle_id ) m_torsor.read( file ); 
  file >> ttag;
  while ( ttag != "</Composite>" ) 
  {
    ObstacleBuilderFactory::reload( ttag, *this, file );
    file >> ttag;
  }
  computeCenterOfMass();
  mother.append( this );

  // Crust thickness 
  m_geoRBWC->setCrustThickness( m_obstacles.front()->getCrustThickness() ); 
  
  // Restricts the geometric directions of translational motion of the box 
  // and initialize position at previous time
  list<Obstacle*>::iterator obstacle = m_obstacles.begin();
  if ( m_restrict_box_geommotion )
    (*obstacle)->setRestrictedGeomDirMotion( m_periodic_directions );
  obstacle++;
  for (; obstacle!=m_obstacles.end();obstacle++)
  {
    if ( isInDomain( (*obstacle)->getPosition() ) ) 
      m_pos_nm1[(*obstacle)->getID()] = *(*obstacle)->getPosition();
  }  
} 
