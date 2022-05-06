#include <FS_2Dcylinder.hh>
#include <MAC.hh>
#include <math.h>
#define EPSILON 1.e-4         // 0.01% of particle size
using std::endl;


//---------------------------------------------------------------------------
FS_2Dcylinder:: FS_2Dcylinder()
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_2Dcylinder:: FS_2Dcylinder" ) ;

  m_space_dimension = 2;
  m_shape_type = GEOM_2DCYLINDER;
  m_agp_2Dcylinder.radius = 0.;
}




//---------------------------------------------------------------------------
FS_2Dcylinder:: FS_2Dcylinder( istream& in, size_t& id_ )
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_2Dcylinder:: FS_2Dcylinder" ) ;

  // Default parameter
  m_space_dimension = 2;
  m_Id = id_;
  m_shape_type = GEOM_2DCYLINDER;

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
FS_2Dcylinder:: ~FS_2Dcylinder()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: ~FS_2Dcylinder" ) ;

}




//---------------------------------------------------------------------------
void FS_2Dcylinder:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: update" ) ;

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
void FS_2Dcylinder:: set( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: set" ) ;

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

  // Force the mass equivalent to a disk
  m_mass = MAC::pi() * m_circumscribed_radius
                     * m_circumscribed_radius
                     * m_density;

  // Set volume
  m_volume = m_mass / m_density ;

  // Set radius
  m_agp_2Dcylinder.radius = m_circumscribed_radius;

}




//---------------------------------------------------------------------------
void FS_2Dcylinder:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;
  out << space << "Shape type = " <<
  	FS_RigidBody::GEOMETRICSHAPE_name[m_shape_type] << endl;
  out << space << "Specific attributes" << endl;
  out << space << three << "Radius = " << m_agp_2Dcylinder.radius << endl;
  out << space << "General attributes" << endl;

  display_general( out, indent_width + 3 );

}




//---------------------------------------------------------------------------
bool FS_2Dcylinder:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: isIn(pt)" ) ;

  bool status = ( m_gravity_center.calcDist( pt ) <= m_agp_2Dcylinder.radius
                                                   * (1. + EPSILON));

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        if (status) break;
        status = (m_gravity_center + (*m_periodic_directions)[i]).calcDist( pt )
                  <= m_agp_2Dcylinder.radius * ( 1. + EPSILON);
     }
  }

  return (status);

  // return ( m_gravity_center.calcDist( pt ) <= m_agp_2Dcylinder.radius );

}




//---------------------------------------------------------------------------
bool FS_2Dcylinder:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: isIn(x,y,z)" ) ;

  geomVector pt(x, y, z);

  bool status = ( m_gravity_center.calcDist( pt ) <= m_agp_2Dcylinder.radius
                                                   * (1. + EPSILON));

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        if (status) break;
        status = (m_gravity_center + (*m_periodic_directions)[i]).calcDist( pt )
                  <= m_agp_2Dcylinder.radius * ( 1. + EPSILON);
     }
  }

  return (status);

  // return ( m_gravity_center.calcDist( x, y, z ) <= m_agp_2Dcylinder.radius );

}




//---------------------------------------------------------------------------
double FS_2Dcylinder:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: level_set_value(pt)" ) ;

  double value = ( m_gravity_center.calcDist( pt ) - m_agp_2Dcylinder.radius
                                                   * (1. + EPSILON));

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        double temp = (m_gravity_center +
                     (*m_periodic_directions)[i]).calcDist( pt )
                     - m_agp_2Dcylinder.radius * ( 1. + EPSILON);
        value = MAC::min(temp,value);
     }
  }

  return (value);

  // return ( m_gravity_center.calcDist( pt ) - m_agp_2Dcylinder.radius );

}




//---------------------------------------------------------------------------
double FS_2Dcylinder:: level_set_value( double const& x
                                , double const& y
                                , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: level_set_value(x,y,z)" ) ;

  geomVector pt(x, y, z);

  double value = ( m_gravity_center.calcDist( pt ) - m_agp_2Dcylinder.radius
                                                   * (1. + EPSILON));

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        double temp = (m_gravity_center +
                     (*m_periodic_directions)[i]).calcDist( pt )
                     - m_agp_2Dcylinder.radius * ( 1. + EPSILON);
        value = MAC::min(temp,value);
     }
  }

  return (value);

  // return ( m_gravity_center.calcDist( x, y, z ) - m_agp_2Dcylinder.radius );

}




//---------------------------------------------------------------------------
struct FS_2Dcylinder_Additional_Param const* FS_2Dcylinder::
	get_ptr_FS_2Dcylinder_Additional_Param() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dcylinder:: get_ptr_FS_2Dcylinder_Additional_Param" ) ;

  return ( &m_agp_2Dcylinder );

}
