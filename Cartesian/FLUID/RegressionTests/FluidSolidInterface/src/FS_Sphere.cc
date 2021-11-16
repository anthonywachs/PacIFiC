#include <FS_Sphere.hh>
using std::endl;


//---------------------------------------------------------------------------
FS_Sphere:: FS_Sphere()
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_Sphere:: FS_Sphere" ) ;

  m_space_dimension = 3;
  m_shape_type = GEOM_SPHERE;
  m_agp_sphere.radius = 0.;
}




//---------------------------------------------------------------------------
FS_Sphere:: FS_Sphere( istream& in, size_t& id_ )
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_Sphere:: FS_Sphere" ) ;

  // Default parameter
  m_space_dimension = 3;
  m_Id = id_;
  m_shape_type = GEOM_SPHERE;

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
FS_Sphere:: ~FS_Sphere()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: ~FS_Sphere" ) ;

}




//---------------------------------------------------------------------------
void FS_Sphere:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: update" ) ;

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
void FS_Sphere:: set( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: set" ) ;

  size_t ncorners, i, nfaces, nper ;

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
  in >> m_gravity_center;
  in >> nfaces;

  // Set volume
  m_volume = m_mass / m_density ;

  // Set radius
  m_agp_sphere.radius = m_circumscribed_radius;

}




//---------------------------------------------------------------------------
void FS_Sphere:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;
  out << space << "Shape type = " <<
  	FS_RigidBody::GEOMETRICSHAPE_name[m_shape_type] << endl;
  out << space << "Specific attributes" << endl;
  out << space << three << "Radius = " << m_agp_sphere.radius << endl;
  out << space << "General attributes" << endl;

  display_general( out, indent_width + 3 );

}




//---------------------------------------------------------------------------
bool FS_Sphere:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: isIn(pt)" ) ;

  return ( m_gravity_center.calcDist( pt ) <= m_agp_sphere.radius );

}




//---------------------------------------------------------------------------
bool FS_Sphere:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: isIn(x,y,z)" ) ;

  return ( m_gravity_center.calcDist( x, y, z ) <= m_agp_sphere.radius );

}




//---------------------------------------------------------------------------
double FS_Sphere:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: level_set_value(pt)" ) ;

  return ( m_gravity_center.calcDist( pt ) - m_agp_sphere.radius );

}




//---------------------------------------------------------------------------
double FS_Sphere:: level_set_value( double const& x
                                , double const& y
                                , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: level_set_value(x,y,z)" ) ;

  return ( m_gravity_center.calcDist( x, y, z ) - m_agp_sphere.radius );

}




//---------------------------------------------------------------------------
tuple<bool,double,size_t> FS_Sphere:: distanceTo( geomVector const& pt,
      	size_t const& direction,
      	bool const& positive ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: distanceTo( pt, direction, positive )" ) ;

  tuple<bool,double,size_t> dist( false, 0., 0 );

  // TO DO
  // Example to assign values
  std::get<0>(dist) = true;
  std::get<1>(dist) = 0.1;
  std::get<2>(dist) = m_Id;

  return ( dist );

}




//---------------------------------------------------------------------------
struct FS_Sphere_Additional_Param const* FS_Sphere::
	get_ptr_FS_Sphere_Additional_Param() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Sphere:: isIn(x,y,z)" ) ;

  return ( &m_agp_sphere );

}
