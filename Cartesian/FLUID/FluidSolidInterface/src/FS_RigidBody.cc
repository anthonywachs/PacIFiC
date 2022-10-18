#include <FS_RigidBody.hh>
using std::endl;


vector<string> FS_RigidBody::GEOMETRICSHAPE_name = { "Sphere", "2D cylinder",
	"3D cylinder", "General polygon", "General polyhedron", "Square",
	"3D Box", "Equilateral triangle", "Regular tetrahedron" };


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

	m_orientation.resize(3);

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
      out << space << three << (*m_periodic_directions)[i](0)
		 				 << three << (*m_periodic_directions)[i](1)
						 << three << (*m_periodic_directions)[i](2) << endl;
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




//---------------------------------------------------------------------------
double FS_RigidBody:: distanceTo( geomVector const& source,
											 geomVector const& rayDir,
				        				    double const& delta )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: distanceTo" ) ;

  double t = 0.;
  // Consider threshold of intersection as 0.0001% change in the relative solution
  double threshold = 1.e-10;

  geomVector rayVec(3);

  // rayVec = source + t * rayDir;
  rayVec(0) = source(0) + t * rayDir(0);
  rayVec(1) = source(1) + t * rayDir(1);
  rayVec(2) = source(2) + t * rayDir(2);


  size_t max_iter = 100, iter = 0;

  // Find the point inside the rigid body
  while ((level_set_value (rayVec) > 0.) && (iter <= max_iter)) {
	  iter++;
	  // rayVec = source + t * rayDir;
	  rayVec(0) = source(0) + t * rayDir(0);
	  rayVec(1) = source(1) + t * rayDir(1);
	  rayVec(2) = source(2) + t * rayDir(2);

	  t += delta/max_iter;
  }

  geomVector a(3), b(3), c(3), c_old(3);

  // a = source;
  a(0) = source(0);
  a(1) = source(1);
  a(2) = source(2);

  // b = rayVec;
  b(0) = rayVec(0);
  b(1) = rayVec(1);
  b(2) = rayVec(2);

  // Initialize c
  // c = source + delta * rayDir ;
  c(0) = source(0) + delta * rayDir(0);
  c(1) = source(1) + delta * rayDir(1);
  c(2) = source(2) + delta * rayDir(2);

  c_old(0) = MAC::abs(a(0) - b(0))/2.;
  c_old(1) = MAC::abs(a(1) - b(1))/2.;
  c_old(2) = MAC::abs(a(2) - b(2))/2.;

  double eps = (b.calcNorm() > 1.e-12) ? a.calcDist(b)/b.calcNorm()
  													: a.calcDist(b);

  max_iter = 500; iter = 0;

  // Bisection method
  while ((eps >= threshold) && (iter <= max_iter)) {
	  // Find middle point
	 c = 0.5 * ( a + b );

	 if (c_old.calcNorm() > 1.e-12) {
		 eps = c.calcDist(c_old)/c_old.calcNorm();
	 } else {
		 eps = c.calcDist(c_old);
	 }

	 // Check if middle point is root
	 if (level_set_value (c) == 0.)
		 break;
	 // Decide the side to repeat the steps
	 else if (level_set_value (c) * level_set_value (a) < 0)
		  b = c;
	 else
		  a = c;

  	 iter += 1;

	 c_old = c;

	 if (iter == max_iter)
		 std::cout << "WARNING: Maxmimum iteration reached for intersection"
					  << " calculation. Proceed with Caution !!!" << endl;
  }

  double dist = 0.;
  for (size_t i = 0; i < 3; ++i)
	 dist += ( source(i) - c(i) ) * ( source(i) - c(i) );

  return (MAC::sqrt(dist));

  // return (source.calcDist(c));

}




//---------------------------------------------------------------------------
geomVector const* FS_RigidBody:: get_ptr_to_gravity_centre() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: get_ptr_to_gravity_centre" ) ;

  return ( &m_gravity_center );

}




//---------------------------------------------------------------------------
double FS_RigidBody:: get_circumscribed_radius() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: get_circumscribed_radius" ) ;

  return ( m_circumscribed_radius );

}




//---------------------------------------------------------------------------
std::tuple<double,double,double> FS_RigidBody:: get_mass_and_density_and_moi() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: get_mass_and_density()" ) ;

  return ( std::make_tuple(m_mass,m_density,m_inertia[0][0]) );

}




//---------------------------------------------------------------------------
void FS_RigidBody:: update_RB_position_and_velocity(geomVector const& pos,
												 					 geomVector const& vel,
																	 geomVector const& ang_vel,
											vector<geomVector> const& periodic_directions,
											double const& time_step)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: update_RB_position_and_velocity()" ) ;

  m_gravity_center(0) = pos(0);
  m_gravity_center(1) = pos(1);
  m_gravity_center(2) = pos(2);
  m_translational_velocity(0) = vel(0);
  m_translational_velocity(1) = vel(1);
  m_translational_velocity(2) = vel(2);
  m_angular_velocity(0) = ang_vel(0);
  m_angular_velocity(1) = ang_vel(1);
  m_angular_velocity(2) = ang_vel(2);

  m_orientation(0) += ang_vel(0)*time_step;
  m_orientation(1) += ang_vel(1)*time_step;
  m_orientation(2) += ang_vel(2)*time_step;

  // Reset the rotation matrix
  double roll = (MAC::pi()/180.)*m_orientation(0);
  double pitch = (MAC::pi()/180.)*m_orientation(1);
  double yaw = (MAC::pi()/180.)*m_orientation(2);

  m_rotation_matrix[0][0] = MAC::cos(yaw)*MAC::cos(pitch);
  m_rotation_matrix[0][1] = MAC::cos(yaw)*MAC::sin(pitch)*MAC::sin(roll)
  								  - MAC::sin(yaw)*MAC::cos(roll);
  m_rotation_matrix[0][2] = MAC::cos(yaw)*MAC::sin(pitch)*MAC::cos(roll)
  								  + MAC::sin(yaw)*MAC::sin(roll);
  m_rotation_matrix[1][0] = MAC::sin(yaw)*MAC::cos(pitch);
  m_rotation_matrix[1][1] = MAC::sin(yaw)*MAC::sin(pitch)*MAC::sin(roll)
  								  + MAC::cos(yaw)*MAC::cos(roll);
  m_rotation_matrix[1][2] = MAC::sin(yaw)*MAC::sin(pitch)*MAC::cos(roll)
  								  - MAC::cos(yaw)*MAC::sin(roll);
  m_rotation_matrix[2][0] = -MAC::sin(pitch);
  m_rotation_matrix[2][1] = MAC::cos(pitch)*MAC::sin(roll);
  m_rotation_matrix[2][2] = MAC::cos(pitch)*MAC::cos(roll);



  if ( m_periodic_directions )
  {
    delete m_periodic_directions;
    m_periodic_directions = NULL;
  }

  geomVector pv( 3 );
  m_periodic_directions = new vector<geomVector>( periodic_directions.size()
  																, pv ) ;
  for (size_t i = 0; i < periodic_directions.size(); ++i) {
	  (*m_periodic_directions)[i](0) = periodic_directions[i](0);
	  (*m_periodic_directions)[i](1) = periodic_directions[i](1);
	  (*m_periodic_directions)[i](2) = periodic_directions[i](2);
  }

}




//---------------------------------------------------------------------------
vector<geomVector> const* FS_RigidBody:: get_ptr_to_periodic_directions() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: get_ptr_to_periodic_directions" ) ;

  return ( m_periodic_directions );

}




//---------------------------------------------------------------------------
geomVector FS_RigidBody:: rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: rigid_body_velocity(pt)" ) ;

  return (m_translational_velocity + (m_angular_velocity^pt));

}




//---------------------------------------------------------------------------
geomVector FS_RigidBody:: rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: rigid_body_angular_velocity" ) ;

  return (m_angular_velocity);

}




//---------------------------------------------------------------------------
void FS_RigidBody:: rotate(geomVector* pt)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: rotate" ) ;

  double delta_x = pt->operator()(0)*m_rotation_matrix[0][0]
                 + pt->operator()(1)*m_rotation_matrix[0][1]
                 + pt->operator()(2)*m_rotation_matrix[0][2];
  double delta_y = pt->operator()(0)*m_rotation_matrix[1][0]
                 + pt->operator()(1)*m_rotation_matrix[1][1]
                 + pt->operator()(2)*m_rotation_matrix[1][2];
  double delta_z = pt->operator()(0)*m_rotation_matrix[2][0]
                 + pt->operator()(1)*m_rotation_matrix[2][1]
                 + pt->operator()(2)*m_rotation_matrix[2][2];

  pt->operator()(0) = delta_x;
  pt->operator()(1) = delta_y;
  pt->operator()(2) = delta_z;

}
