#include <FS_3Dcylinder.hh>
#include <MAC.hh>
using std::endl;


//---------------------------------------------------------------------------
FS_3Dcylinder:: FS_3Dcylinder()
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_3Dcylinder:: FS_3Dcylinder" ) ;

  m_space_dimension = 3;
  m_shape_type = GEOM_3DCYLINDER;
  m_agp_3dcyl.BottomCenter.resize(3);
  m_agp_3dcyl.TopCenter.resize(3);
  m_agp_3dcyl.BottomToTopVec.resize(3);
  m_agp_3dcyl.RadialRefVec.resize(3);
  m_agp_3dcyl.cylinder_radius = 0.;
  m_agp_3dcyl.cylinder_height = 0.;
}




//---------------------------------------------------------------------------
FS_3Dcylinder:: FS_3Dcylinder( istream& in, size_t& id_ )
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_3Dcylinder:: FS_3Dcylinder" ) ;

  // Default parameter
  m_space_dimension = 3;
  m_Id = id_;
  m_shape_type = GEOM_3DCYLINDER;

  // Resize parameters
  m_gravity_center.resize(3);
  m_translational_velocity.resize(3);
  m_angular_velocity.resize(3);
  m_hydro_force.resize(3);
  m_hydro_torque.resize(3);
  m_heat_flux.resize(3);
  m_agp_3dcyl.BottomCenter.resize(3);
  m_agp_3dcyl.TopCenter.resize(3);
  m_agp_3dcyl.BottomToTopVec.resize(3);
  m_agp_3dcyl.RadialRefVec.resize(3);

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
FS_3Dcylinder:: ~FS_3Dcylinder()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: ~FS_3Dcylinder" ) ;

}




//---------------------------------------------------------------------------
void FS_3Dcylinder:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: update" ) ;

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
void FS_3Dcylinder:: set( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: set" ) ;

  size_t ncorners, i, nfaces, nper ;
  geomVector vv(3);

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

  // Read the bottom and top disk centers
  // Second point is used to compute the radial reference vector
  // and the cylinder radius
  // The radial reference vector is a vector pointing from the axis to the
  // cylinder surface in the plane of the disks
  in >> m_agp_3dcyl.BottomCenter;
  in >> vv;
  m_agp_3dcyl.RadialRefVec = vv - m_agp_3dcyl.BottomCenter ;
  m_agp_3dcyl.cylinder_radius = m_agp_3dcyl.RadialRefVec.calcNorm() ;
  in >> m_agp_3dcyl.TopCenter;
  m_agp_3dcyl.BottomToTopVec = m_agp_3dcyl.TopCenter - m_agp_3dcyl.BottomCenter ;
  m_agp_3dcyl.cylinder_height = m_agp_3dcyl.BottomToTopVec.calcNorm() ;
  in >> nfaces;

  // Set volume
  m_volume = m_mass / m_density ;
}




//---------------------------------------------------------------------------
void FS_3Dcylinder:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;
  out << space << "Shape type = " <<
  	FS_RigidBody::GEOMETRICSHAPE_name[m_shape_type] << endl;
  out << space << "Specific attributes" << endl;
  out << space << three << "Radius = " << m_agp_3dcyl.cylinder_radius << endl;
  out << space << three << "Height = " << m_agp_3dcyl.cylinder_height << endl;
  out << space << three << "Center of the bottom disk = "
  	                     << m_agp_3dcyl.BottomCenter(0)
               << three << m_agp_3dcyl.BottomCenter(1)
               << three << m_agp_3dcyl.BottomCenter(2) << endl;
  out << space << three << "Center of the top disk = "
  	                     << m_agp_3dcyl.TopCenter(0)
               << three << m_agp_3dcyl.TopCenter(1)
               << three << m_agp_3dcyl.TopCenter(2) << endl;
  out << space << three << "Axial vector from bottom to top disk center = "
  	                     << m_agp_3dcyl.BottomToTopVec(0)
               << three << m_agp_3dcyl.BottomToTopVec(1)
               << three << m_agp_3dcyl.BottomToTopVec(2) << endl;
  out << space << three << "Radial reference vector = "
  	                     << m_agp_3dcyl.RadialRefVec(0)
               << three << m_agp_3dcyl.RadialRefVec(1)
               << three << m_agp_3dcyl.RadialRefVec(2) << endl;
  out << space << "General attributes" << endl;
  display_general( out, indent_width + 3 );

}




//---------------------------------------------------------------------------
bool FS_3Dcylinder:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: isIn(pt)" ) ;

  bool b_isIn = false;

  geomVector BottomToPoint( pt - m_agp_3dcyl.BottomCenter );
  double dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec )
	                        /  m_agp_3dcyl.cylinder_height ;

  if ( dot <= m_agp_3dcyl.cylinder_height && dot >= 0. )
    if ( BottomToPoint.calcNormSquare() - dot * dot <=
         m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius )
      b_isIn = true ;

  return ( b_isIn );

}




//---------------------------------------------------------------------------
bool FS_3Dcylinder:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: isIn(x,y,z)" ) ;

  bool b_isIn = false;

  geomVector pt(x,y,z);

  geomVector BottomToPoint( pt - m_agp_3dcyl.BottomCenter );
  double dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec )
                           /  m_agp_3dcyl.cylinder_height ;

  if ( dot <= m_agp_3dcyl.cylinder_height && dot >= 0. )
    if ( BottomToPoint.calcNormSquare() - dot * dot <=
         m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius )
      b_isIn = true ;

  return ( b_isIn );

}




//---------------------------------------------------------------------------
double FS_3Dcylinder:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: level_set_value(pt)" ) ;

  geomVector BottomToPoint( pt - m_agp_3dcyl.BottomCenter );
  double dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec )
	                        /  m_agp_3dcyl.cylinder_height ;

  double rad_sq_dist = BottomToPoint.calcNormSquare() - dot * dot -
         m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius;

  double value = 0.;

  if ( dot < m_agp_3dcyl.cylinder_height && dot > 0. && rad_sq_dist < 0.) {
     value = MAC::abs(MAC::abs(dot) - m_agp_3dcyl.cylinder_height)*rad_sq_dist;
  } else {
     value = MAC::min(MAC::abs(MAC::abs(dot)-m_agp_3dcyl.cylinder_height)
                    , MAC::abs(rad_sq_dist));
  }

  return ( value );

}




//---------------------------------------------------------------------------
double FS_3Dcylinder:: level_set_value( double const& x
                                      , double const& y
                                      , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: level_set_value(x,y,z)" ) ;

  geomVector pt(x,y,z);

  geomVector BottomToPoint( pt - m_agp_3dcyl.BottomCenter );
  double dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec )
	                        /  m_agp_3dcyl.cylinder_height ;

  double rad_sq_dist = BottomToPoint.calcNormSquare() - dot * dot -
         m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius;

  double value = 0.;

  if ( dot < m_agp_3dcyl.cylinder_height && dot > 0. && rad_sq_dist < 0.) {
     value = MAC::abs(MAC::abs(dot) - m_agp_3dcyl.cylinder_height)*rad_sq_dist;
  } else {
     value = MAC::min(MAC::abs(MAC::abs(dot)-m_agp_3dcyl.cylinder_height)
                    , MAC::abs(rad_sq_dist));
  }

  return ( value );

}




//---------------------------------------------------------------------------
geomVector FS_3Dcylinder:: rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: rigid_body_velocity(pt)" ) ;

  return (m_translational_velocity
        + (m_angular_velocity^(pt-m_gravity_center)));

}




//---------------------------------------------------------------------------
struct FS_3Dcylinder_Additional_Param const* FS_3Dcylinder::
	get_ptr_FS_3Dcylinder_Additional_Param() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: isIn(x,y,z)" ) ;

  return( &m_agp_3dcyl );

}
