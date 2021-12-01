#include <DS_RigidBody.hh>
#include <FS_RigidBody.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_RigidBody:: DS_RigidBody()
//---------------------------------------------------------------------------
  : m_geometric_rigid_body( NULL )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_RigidBody:: DS_RigidBody( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : m_geometric_rigid_body( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

  m_halo_zone.reserve(2);
  m_halo_zone.push_back(new geomVector(3));
  m_halo_zone.push_back(new geomVector(3));

}




//---------------------------------------------------------------------------
DS_RigidBody:: ~DS_RigidBody()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: ~DS_RigidBody" ) ;

  if ( !m_surface_points.empty() ) m_surface_points.clear();

}




//---------------------------------------------------------------------------
bool DS_RigidBody:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: isIn(pt)" ) ;

  return ( m_geometric_rigid_body->isIn( pt ) );

}




//---------------------------------------------------------------------------
bool DS_RigidBody:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: isIn(x,y,z)" ) ;

  return ( m_geometric_rigid_body->isIn( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_RigidBody:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: level_set_value(pt)" ) ;

  return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_RigidBody:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: level_set_value(x,y,z)" ) ;

  return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

}




//---------------------------------------------------------------------------
void DS_RigidBody:: compute_surface_integrals_hydro_force_torque(
      	FV_DiscreteField const* PP,
	FV_DiscreteField const* UU )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: isIn(x,y,z)" ) ;

  geomVector hydro_force(3), hydro_torque(3);

  // Compute the surface intergals for the hydro force and torque
  // based on the distributed surface points
  MAC::out() << "DS_RigidBody:: compute_surface_integrals_hydro_force_torque - "
  	"Actual computation requires programming" << endl;

  // Store the values in the corresponding geometric sphere
  m_geometric_rigid_body->set_hydro_force( hydro_force );
  m_geometric_rigid_body->set_hydro_torque( hydro_torque );

}




//---------------------------------------------------------------------------
double DS_RigidBody:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_distanceTo" ) ;

  return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));

}




//---------------------------------------------------------------------------
geomVector DS_RigidBody:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
void DS_RigidBody:: initialize_surface_variables( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: initialize_surface_variables" ) ;

  m_surface_points.reserve( Ntot );
  m_surface_area.reserve( Ntot );
  m_surface_normal.reserve( Ntot );
  m_surface_Pforce.reserve( Ntot );
  m_surface_Vforce.reserve( Ntot );

  geomVector vvv(3);

   for (size_t i = 0; i < Ntot; ++i) {
      m_surface_points.push_back( new geomVector(3) );
      m_surface_area.push_back( new geomVector(1) );
      m_surface_normal.push_back( new geomVector(3) );
      m_surface_Pforce.push_back( vvv );
      m_surface_Vforce.push_back( vvv );
   }

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_RigidBody:: get_rigid_body_surface_points( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_rigid_body_surface_points" ) ;

  return (m_surface_points);

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_RigidBody:: get_rigid_body_surface_normals( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_rigid_body_surface_normals" ) ;

  return (m_surface_normal);

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_RigidBody:: get_rigid_body_surface_areas( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_rigid_body_surface_areas" ) ;

  return (m_surface_area);

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_RigidBody:: get_rigid_body_haloZone( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_rigid_body_haloZone" ) ;

  return (m_halo_zone);

}




//---------------------------------------------------------------------------
void DS_RigidBody:: write_surface_discretization( const std::string& file )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: write_surface_discretization()" ) ;

  std::ofstream out;
  out.open(file.c_str());
  out << "x ,y ,z ,nx ,ny ,nz ,area" << endl;

  for (size_t i = 0; i < m_surface_area.size(); i++) {
     out << m_surface_points[i]->operator()(0) << " ,"
         << m_surface_points[i]->operator()(1) << " ,"
         << m_surface_points[i]->operator()(2) << " ,"
         << m_surface_normal[i]->operator()(0) << " ,"
         << m_surface_normal[i]->operator()(1) << " ,"
         << m_surface_normal[i]->operator()(2) << " ,"
         << m_surface_area[i]->operator()(0) << endl;
  }

  out.close();

}
