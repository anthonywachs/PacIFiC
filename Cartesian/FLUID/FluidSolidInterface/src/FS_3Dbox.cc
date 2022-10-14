#include <FS_3Dbox.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
FS_3Dbox:: FS_3Dbox()
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_3Dbox:: FS_3Dbox" ) ;

  m_space_dimension = 3;
  m_shape_type = GEOM_3DBOX;

}




//---------------------------------------------------------------------------
FS_3Dbox:: FS_3Dbox( istream& in, size_t& id_ )
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_3Dbox:: FS_3Dbox" ) ;

  // Default parameter
  m_space_dimension = 3;
  m_Id = id_;
  m_shape_type = GEOM_3DBOX;

  // Resize parameters
  m_gravity_center.resize(3);
  m_translational_velocity.resize(3);
  m_angular_velocity.resize(3);
  m_hydro_force.resize(3);
  m_hydro_torque.resize(3);
  m_heat_flux.resize(3);

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
FS_3Dbox:: ~FS_3Dbox()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: ~FS_3Dbox" ) ;

}




//---------------------------------------------------------------------------
void FS_3Dbox:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: update" ) ;

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
void FS_3Dbox:: set( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: set" ) ;

  size_t ncorners, i, nfaces, nper ;
  geomVector node(3);

  // Read the input stream
  in >> m_type >> m_translational_velocity >>
  	m_angular_velocity >> m_density >> m_mass;
  in >> m_inertia[0][0];
  in >> m_inertia[0][1];
  in >> m_inertia[0][2];
  in >> m_inertia[1][1];
  in >> m_inertia[1][2];
  in >> m_inertia[2][2];
  m_inertia[1][0] = m_inertia[0][1];
  m_inertia[2][0] = m_inertia[0][2];
  m_inertia[2][1] = m_inertia[1][2];

  in >> m_rotation_matrix[0][0];
  in >> m_rotation_matrix[0][1];
  in >> m_rotation_matrix[0][2];
  in >> m_rotation_matrix[1][0];
  in >> m_rotation_matrix[1][1];
  in >> m_rotation_matrix[1][2];
  in >> m_rotation_matrix[2][0];
  in >> m_rotation_matrix[2][1];
  in >> m_rotation_matrix[2][2];

  in >> m_gravity_center;

  if ( m_periodic_directions )
  {
    delete m_periodic_directions;
    m_periodic_directions = NULL;
  }

  if ( m_type == "PP" )
  {
    geomVector pv( 3 );
    in >> nper;
    m_periodic_directions = new vector<geomVector>( nper, pv ) ;
    for (i=0;i<nper;++i)
      in >> (*m_periodic_directions)[i](0) >> (*m_periodic_directions)[i](1)
      	>> (*m_periodic_directions)[i](2);
  }
  in >> m_circumscribed_radius >> ncorners;

  // Read the corners and faces numbering of the 3D box

  // Build the polyhedron corners
  if ( m_agp_3dbox.corners.empty() ) {
    m_agp_3dbox.corners.reserve( ncorners );
    for (i = 0; i < ncorners; i++) {
      in >> node;
      m_agp_3dbox.corners.push_back(node);
    }
  } else {
    for (i = 0; i < ncorners; i++) in >> m_agp_3dbox.corners[i];
  }

  // Build the reference polyhedron corners
  if ( m_agp_3dbox.ref_corners.empty() ) {
    m_agp_3dbox.ref_corners.reserve( ncorners );
    for (i = 0; i < ncorners; i++) {
      m_agp_3dbox.ref_corners.push_back(node);
    }
  }

  compute_reverseTransformationOfCorners( );

  // build the polyhedron faces
  size_t nbFaceCorners, localCornerIdx;

  in >> nfaces;

  if ( m_agp_3dbox.facesVec.empty() ) {
    m_agp_3dbox.facesVec.reserve( nfaces );
    for ( i = 0; i <= nfaces-1; i++ ) {
      in >> nbFaceCorners;

      vector<size_t> localVect;
      localVect.reserve(nbFaceCorners);

      for ( size_t idx = 0; idx < nbFaceCorners; idx++ ) {
        in >> localCornerIdx;
        localVect.push_back( localCornerIdx );
      }
      m_agp_3dbox.facesVec.push_back( localVect );
    }
  } else {
    for ( i = 0; i <= nfaces-1; i++ ) {
      in >> nbFaceCorners;
      for ( size_t idx = 0; idx < nbFaceCorners; idx++ )
        in >> m_agp_3dbox.facesVec[i][idx];
    }
  }

  // Set volume
  m_volume = m_mass / m_density ;

  // Reset the rotation matrix
  double roll = (MAC::pi()/180.)*m_orientation(0);
  double pitch = (MAC::pi()/180.)*m_orientation(1);
  double yaw = (MAC::pi()/180.)*m_orientation(2);

  m_rotation_matrix[0][0] = MAC::cos(yaw)*MAC::cos(pitch);
  m_rotation_matrix[0][1] = MAC::cos(yaw)*MAC::sin(pitch)*MAC::sin(roll) - MAC::sin(yaw)*MAC::cos(roll);
  m_rotation_matrix[0][2] = MAC::cos(yaw)*MAC::sin(pitch)*MAC::cos(roll) + MAC::sin(yaw)*MAC::sin(roll);
  m_rotation_matrix[1][0] = MAC::sin(yaw)*MAC::cos(pitch);
  m_rotation_matrix[1][1] = MAC::sin(yaw)*MAC::sin(pitch)*MAC::sin(roll) + MAC::cos(yaw)*MAC::cos(roll);
  m_rotation_matrix[1][2] = MAC::sin(yaw)*MAC::sin(pitch)*MAC::cos(roll) - MAC::cos(yaw)*MAC::sin(roll);
  m_rotation_matrix[2][0] = -MAC::sin(pitch);
  m_rotation_matrix[2][1] = MAC::cos(pitch)*MAC::sin(roll);
  m_rotation_matrix[2][2] = MAC::cos(pitch)*MAC::cos(roll);
}




//---------------------------------------------------------------------------
void FS_3Dbox:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;
  out << space << "Shape type = " <<
  	FS_RigidBody::GEOMETRICSHAPE_name[m_shape_type] << endl;
  out << space << "Specific attributes" << endl;
  for (int i = 0; i < (int) m_agp_3dbox.corners.size(); i++)
     out << space << three << i+1 << "st Corner = " << m_agp_3dbox.corners[i](0)
                                           << three << m_agp_3dbox.corners[i](1)
                                           << three << m_agp_3dbox.corners[i](2)
                                           << endl;

  for (int i = 0; i < (int) m_agp_3dbox.facesVec.size(); i++) {
     out << space << three << i+1 << "st faceVec = ";
     for (int j = 0; j < (int) m_agp_3dbox.facesVec[i].size(); j++)
        out << three << m_agp_3dbox.facesVec[i][j];
     out << endl;
  }


  out << space << "General attributes" << endl;
  display_general( out, indent_width + 3 );

}




//---------------------------------------------------------------------------
bool FS_3Dbox:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: isIn(pt)" ) ;

  bool b_isIn = false;

  for (vector< vector<size_t> >::const_iterator
         faceIter = m_agp_3dbox.facesVec.begin();
         faceIter != m_agp_3dbox.facesVec.end() && !b_isIn; ++faceIter )
     for (size_t idxPts = 2; idxPts < 4 && !b_isIn; ++idxPts )
        if ( checkPointInTetrahedron( m_agp_3dbox.corners[ (*faceIter)[0] ],
                               m_agp_3dbox.corners[ (*faceIter)[idxPts-1] ],
                                 m_agp_3dbox.corners[ (*faceIter)[idxPts] ],
                                      m_gravity_center, pt ) )
           b_isIn = true;

  return ( b_isIn );

}




//---------------------------------------------------------------------------
bool FS_3Dbox:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: isIn(x,y,z)" ) ;

  bool b_isIn = false;

  geomVector pt( x, y, z );

  for (vector< vector<size_t> >::const_iterator
         faceIter = m_agp_3dbox.facesVec.begin();
         faceIter != m_agp_3dbox.facesVec.end() && !b_isIn; ++faceIter )
     for (size_t idxPts = 2; idxPts < 4 && !b_isIn; ++idxPts ) {
        if ( checkPointInTetrahedron( m_agp_3dbox.corners[ (*faceIter)[0] ],
                               m_agp_3dbox.corners[ (*faceIter)[idxPts-1] ],
                                 m_agp_3dbox.corners[ (*faceIter)[idxPts] ],
                                      m_gravity_center, pt ) )
           b_isIn = true;
        }

  return ( b_isIn );

}




//---------------------------------------------------------------------------
double FS_3Dbox:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: level_set_value(pt)" ) ;

  double value = 1.e10;

  for (vector< vector<size_t> >::const_iterator
           faceIter = m_agp_3dbox.facesVec.begin();
           faceIter != m_agp_3dbox.facesVec.end() && (value > 0.);
           ++faceIter )
     for (size_t idxPts = 2; idxPts < 4 && (value > 0.); ++idxPts ) {
        double temp = DistOfPointFromTetrahedron(
                                 m_agp_3dbox.corners[ (*faceIter)[0] ],
                                 m_agp_3dbox.corners[ (*faceIter)[idxPts-1] ],
                                 m_agp_3dbox.corners[ (*faceIter)[idxPts] ],
                                 m_gravity_center, pt );

        if (temp >= -1.e-7) {
           value = std::max(temp, value);
        } else {
           value = temp;
        }
     }

  return(value);

}




//---------------------------------------------------------------------------
double FS_3Dbox:: level_set_value( double const& x
                                 , double const& y
                                 , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: level_set_value(x,y,z)" ) ;

  double value = 1.e10;       // +ve outside, -ve inside

  geomVector pt( x, y, z );

  for (vector< vector<size_t> >::const_iterator
           faceIter = m_agp_3dbox.facesVec.begin();
           faceIter != m_agp_3dbox.facesVec.end() && (value > 0.);
           ++faceIter )
     for (size_t idxPts = 2; idxPts < 4 && (value > 0.); ++idxPts ) {
        double temp = DistOfPointFromTetrahedron(
                                 m_agp_3dbox.corners[ (*faceIter)[0] ],
                                 m_agp_3dbox.corners[ (*faceIter)[idxPts-1] ],
                                 m_agp_3dbox.corners[ (*faceIter)[idxPts] ],
                                 m_gravity_center, pt );

        if (temp >= -1.e-7) {
           value = std::min(temp, value);
        } else {
           value = temp;
        }
     }

  return(value);
}




//---------------------------------------------------------------------------
struct FS_3Dbox_Additional_Param const* FS_3Dbox::
	get_ptr_FS_3Dbox_Additional_Param() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: get_ptr_FS_3Dbox_Additional_Param()" ) ;

  return( &m_agp_3dbox );

}




//---------------------------------------------------------------------------
double FS_3Dbox::calcPointDeterm4by4( const geomVector &pointOne,
        const geomVector &pointTwo, const geomVector &pointThree,
        const geomVector &pointFour ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: calcPointDeterm4by4" ) ;

  double x2 = pointTwo(0);
  double y2 = pointTwo(1);
  double z2 = pointTwo(2);
  double x3 = pointThree(0);
  double y3 = pointThree(1);
  double z3 = pointThree(2);
  double x4 = pointFour(0);
  double y4 = pointFour(1);
  double z4 = pointFour(2);


  double retVal =
  pointOne(0) * ( y2*z3 + y3*z4 + y4*z2 - y4*z3 - y3*z2 - y2*z4 ) -
  pointOne(1) * ( x2*z3 + x3*z4 + x4*z2 - x4*z3 - x3*z2 - x2*z4 ) +
  pointOne(2) * ( x2*y3 + x3*y4 + x4*y2 - x4*y3 - x3*y2 - x2*y4 ) -
        ( x2*y3*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2 - x3*y2*z4 - x2*y4*z3 ) ;

  return retVal;

}




//---------------------------------------------------------------------------
bool FS_3Dbox::checkPointInTetrahedron( const geomVector &pointOne,
        const geomVector &pointTwo, const geomVector &pointThree,
        const geomVector &pointFour, const geomVector &pointToCheck ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: checkPointInTetrahedron" ) ;

// Link: http://steve.hollasch.net/cgindex/geometry/ptintet.html

  double detTot   = calcPointDeterm4by4( pointOne, pointTwo,
        pointThree, pointFour );
  double detOne   = calcPointDeterm4by4( pointToCheck, pointTwo,
        pointThree, pointFour );
  double detTwo   = calcPointDeterm4by4( pointOne, pointToCheck,
        pointThree, pointFour );
  double detThree = calcPointDeterm4by4( pointOne, pointTwo,
        pointToCheck, pointFour );
  double detFour  = calcPointDeterm4by4( pointOne, pointTwo,
        pointThree, pointToCheck );

  bool in=false;

  double sumSubElem = detOne + detTwo + detThree + detFour;

  if ( fabs( detTot - sumSubElem ) > 1.e-7  )
        std::cout << "ERROR: summation error in determinat 3D : "
        << detTot - sumSubElem << endl;

  if ( detTot == 0 ) {
        std::cout << "degenerated tetrahedron: det == 0 " << endl;
        abort();
  }

  if ( ( detTot*detOne >= -1.e-7 )  && ( detTot*detTwo >= -1.e-7 ) &&
       ( detTot*detThree >= -1.e-7 ) && ( detTot*detFour >= -1.e-7 ) )
      in=true;

  return in;

}




//---------------------------------------------------------------------------
double FS_3Dbox::DistOfPointFromTetrahedron( const geomVector &pointOne,
        const geomVector &pointTwo, const geomVector &pointThree,
        const geomVector &pointFour, const geomVector &pointToCheck ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: DistOfPointFromTetrahedron" ) ;

  double detTot   = calcPointDeterm4by4( pointOne, pointTwo,
        pointThree, pointFour );
  double detOne   = calcPointDeterm4by4( pointToCheck, pointTwo,
        pointThree, pointFour );
  double detTwo   = calcPointDeterm4by4( pointOne, pointToCheck,
        pointThree, pointFour );
  double detThree = calcPointDeterm4by4( pointOne, pointTwo,
        pointToCheck, pointFour );
  double detFour  = calcPointDeterm4by4( pointOne, pointTwo,
        pointThree, pointToCheck );

  double out_dist = 0.;

  double sumSubElem = detOne + detTwo + detThree + detFour;

  if ( fabs( detTot - sumSubElem ) > 1.e-7  )
     std::cout << "ERROR: summation error in determinat 3D : "
               << detTot - sumSubElem << endl;

  if ( detTot == 0 ) {
     std::cout << "degenerated tetrahedron: det == 0 " << endl;
     abort();
  }

  if ( ( detTot*detOne >= -1.e-7 )  && ( detTot*detTwo >= -1.e-7 ) &&
       ( detTot*detThree >= -1.e-7 ) && ( detTot*detFour >= -1.e-7 ) ) {
     out_dist = -1.*fabs(sumSubElem);
  } else {
     out_dist = fabs(sumSubElem);
  }

  return out_dist;

}




//---------------------------------------------------------------------------
void FS_3Dbox::compute_reverseTransformationOfCorners( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dbox:: compute_reverseTransformationOfCorners()" ) ;

  for (int i = 0; i < (int) m_agp_3dbox.ref_corners.size(); i++) {
     geomVector pt(m_agp_3dbox.corners[i]);

     m_agp_3dbox.ref_corners[i](0) = pt(0)*m_rotation_matrix[0][0]
                                   + pt(1)*m_rotation_matrix[1][0]
                                   + pt(2)*m_rotation_matrix[2][0];
     m_agp_3dbox.ref_corners[i](1) = pt(0)*m_rotation_matrix[0][1]
                                   + pt(1)*m_rotation_matrix[1][1]
                                   + pt(2)*m_rotation_matrix[2][1];
     m_agp_3dbox.ref_corners[i](2) = pt(0)*m_rotation_matrix[0][2]
                                   + pt(1)*m_rotation_matrix[1][2]
                                   + pt(2)*m_rotation_matrix[2][2];
  }

}
