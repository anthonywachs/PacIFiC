#include <FS_RigidBody.hh>
using std::endl;


vector<string> FS_RigidBody::GEOMETRICSHAPE_name = { "Sphere", "2D cylinder",
	"3D cylinder", "General polygon", "General polyhedron", "Square",
	"Cube", "Equilateral triangle", "Regular tetrahedron" };


//---------------------------------------------------------------------------
FS_RigidBody:: FS_RigidBody()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: FS_RigidBody" ) ;

  m_periodic_directions = NULL;
  m_temperature = 0.;
  m_rotation_matrix[0][0] = m_rotation_matrix[0][1] = m_rotation_matrix[0][2]
  	= m_rotation_matrix[1][0] = m_rotation_matrix[1][1]
  	= m_rotation_matrix[1][2] = m_rotation_matrix[2][0]
	= m_rotation_matrix[2][1] = m_rotation_matrix[2][2] = 0.;

}




//---------------------------------------------------------------------------
FS_RigidBody:: ~FS_RigidBody()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: ~FS_RigidBody" ) ;

  if ( m_periodic_directions )
  {
    m_periodic_directions->clear();
    delete m_periodic_directions ;
  }

}




//---------------------------------------------------------------------------
GEOMETRICSHAPE FS_RigidBody:: get_shape_type() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: get_shape_type" ) ;

  return ( m_shape_type );

}




//---------------------------------------------------------------------------
string FS_RigidBody:: get_type() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: get_type" ) ;

  return ( m_type );

}




//---------------------------------------------------------------------------
void FS_RigidBody:: display_general( ostream& out,
	size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: display_general" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Space dimension = " << m_space_dimension << endl;
  out << space << "Identification number = " << m_Id << endl;
  out << space << "Rigid body type = " << m_type << endl;
  out << space << "Density = " << m_density << endl;
  out << space << "Volume = " << m_volume << endl;
  out << space << "Mass = " << m_mass << endl;
  out << space << "Circumscribed radius = " << m_circumscribed_radius << endl;
  out << space << "Moment of inertia matrix J" << endl;
  out << space << three << "Jxx = " << m_inertia[0][0] << endl;
  out << space << three << "Jxy = " << m_inertia[0][1] << endl;
  out << space << three << "Jxz = " << m_inertia[0][2] << endl;
  out << space << three << "Jyy = " << m_inertia[1][1] << endl;
  out << space << three << "Jyz = " << m_inertia[1][2] << endl;
  out << space << three << "Jzz = " << m_inertia[2][2] << endl;
  out << space << "Gravity center = " << m_gravity_center(0)
  									  << three << m_gravity_center(1)
									  << three << m_gravity_center(2) << endl;
  out << space << "Rotation matrix M" << endl;
  out << space << three << "Mxx = " << m_rotation_matrix[0][0] << endl;
  out << space << three << "Mxy = " << m_rotation_matrix[0][1] << endl;
  out << space << three << "Mxz = " << m_rotation_matrix[0][2] << endl;
  out << space << three << "Myx = " << m_rotation_matrix[1][0] << endl;
  out << space << three << "Myy = " << m_rotation_matrix[1][1] << endl;
  out << space << three << "Myz = " << m_rotation_matrix[1][2] << endl;
  out << space << three << "Mzx = " << m_rotation_matrix[2][0] << endl;
  out << space << three << "Mzy = " << m_rotation_matrix[2][1] << endl;
  out << space << three << "Mzz = " << m_rotation_matrix[2][2] << endl;
  out << space << "Translational velocity = " << m_translational_velocity(0)
  												 << three << m_translational_velocity(1)
												 << three << m_translational_velocity(2)
												 << endl;
  out << space << "Angular velocity = " << m_angular_velocity(0)
   									 << three << m_angular_velocity(1)
										 << three << m_angular_velocity(2) << endl;
  out << space << "Hydrodynamic force = " << m_hydro_force(0)
   										<< three << m_hydro_force(1)
											<< three << m_hydro_force(2) << endl;
  out << space << "Hydrodynamic torque = " << m_hydro_torque(0)
  											<< three << m_hydro_torque(1)
											<< three << m_hydro_torque(2) << endl;
  out << space << "Temperature = " << m_temperature << endl;
  out << space << "Heat flux = " << m_heat_flux(0)
   							<< three << m_heat_flux(1)
								<< three << m_heat_flux(2) << endl;
  if ( m_periodic_directions )
  {
    out << space << "Periodic directions:" << endl;
    for (size_t i=0;i<m_periodic_directions->size();++i)
      out << space << three << (*m_periodic_directions)[i] << endl;
  }

}




//---------------------------------------------------------------------------
void FS_RigidBody:: set_hydro_force( geomVector const& hf )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: set_hydro_force" ) ;

  m_hydro_force = hf ;

}




//---------------------------------------------------------------------------
void FS_RigidBody:: set_hydro_torque( geomVector const& ht )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: set_hydro_torque" ) ;

  m_hydro_torque = ht ;

}




//---------------------------------------------------------------------------
void FS_RigidBody:: nullify_velocity()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: nullify_velocity" ) ;

  m_translational_velocity.setVecZero();
  m_angular_velocity.setVecZero();

}




//---------------------------------------------------------------------------
void FS_RigidBody:: change_from_particle_to_obstacle()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: change_from_particle_to_obstacle" ) ;

  if ( m_type == "P" ) m_type = "O";
  else if ( m_type == "PP" ) m_type = "PO";

}
