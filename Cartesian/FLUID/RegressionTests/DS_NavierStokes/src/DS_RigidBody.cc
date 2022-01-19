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
geomVector const* DS_RigidBody:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_ptr_to_gravity_centre( )" ) ;

  return (dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre());

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
void DS_RigidBody:: update_Pforce_on_surface_point( size_t const& i
                                                  , geomVector const& value )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: update_Pforce_on_surface_point" ) ;

  m_surface_Pforce[i] = value;

}




//---------------------------------------------------------------------------
void DS_RigidBody:: update_Vforce_on_surface_point( size_t const& i
                                                  , geomVector const& value )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: update_Vforce_on_surface_point" ) ;

  m_surface_Vforce[i] = value;

}




//---------------------------------------------------------------------------
void DS_RigidBody:: write_surface_discretization( const std::string& file
                                      , MAC_Communicator const* m_macCOMM)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: write_surface_discretization" ) ;

  std::ofstream out;

  if (m_macCOMM->rank() == 0) {
     out.open(file.c_str());
     out << "x ,y ,z ,nx ,ny ,nz ,area ,Fpx ,Fpy ,Fpz ,Fvx ,Fvy , Fvz " << endl;
  }

  for (size_t i = 0; i < m_surface_area.size(); i++) {
     m_surface_Pforce[i](0) = m_macCOMM->sum(m_surface_Pforce[i](0));
     m_surface_Pforce[i](1) = m_macCOMM->sum(m_surface_Pforce[i](1));
     m_surface_Pforce[i](2) = m_macCOMM->sum(m_surface_Pforce[i](2));

     m_surface_Vforce[i](0) = m_macCOMM->sum(m_surface_Vforce[i](0));
     m_surface_Vforce[i](1) = m_macCOMM->sum(m_surface_Vforce[i](1));
     m_surface_Vforce[i](2) = m_macCOMM->sum(m_surface_Vforce[i](2));

     if (m_macCOMM->rank() == 0) {
        out << m_surface_points[i]->operator()(0) << " ,"
            << m_surface_points[i]->operator()(1) << " ,"
            << m_surface_points[i]->operator()(2) << " ,"
            << m_surface_normal[i]->operator()(0) << " ,"
            << m_surface_normal[i]->operator()(1) << " ,"
            << m_surface_normal[i]->operator()(2) << " ,"
            << m_surface_area[i]->operator()(0) << " ,"
            << m_surface_Pforce[i](0) << " ,"
            << m_surface_Pforce[i](1) << " ,"
            << m_surface_Pforce[i](2) << " ,"
            << m_surface_Vforce[i](0) << " ,"
            << m_surface_Vforce[i](1) << " ,"
            << m_surface_Vforce[i](2) << endl;
     }
  }

  if (m_macCOMM->rank() == 0) out.close();

}
