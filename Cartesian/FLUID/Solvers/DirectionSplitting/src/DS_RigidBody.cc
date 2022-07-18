#include <DS_RigidBody.hh>
#include <FS_RigidBody.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_RigidBody:: DS_RigidBody()
//---------------------------------------------------------------------------
  : m_geometric_rigid_body( NULL )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

  m_halo_zone.reserve(2);
  m_halo_zone.push_back(new geomVector(3));
  m_halo_zone.push_back(new geomVector(3));

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
void DS_RigidBody:: translate_surface_points( geomVector const& delta)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: translate_surface_points()" ) ;

  // Translate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
     m_surface_points[i]->translate(delta);
  }

}




//---------------------------------------------------------------------------
void DS_RigidBody:: initialize_surface_variables( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: initialize_surface_variables" ) ;

  if (m_surface_points.empty()) {
     m_surface_points.reserve( Ntot );
     m_surface_area.reserve( Ntot );
     m_surface_normal.reserve( Ntot );
     m_surface_Pforce.reserve( Ntot );
     m_surface_Vforce.reserve( Ntot );
     m_surface_Tgrad.reserve( Ntot );

     geomVector vvv(3);

     for (size_t i = 0; i < Ntot; ++i) {
        m_surface_points.push_back( new geomVector(3) );
        m_surface_area.push_back( new geomVector(1) );
        m_surface_normal.push_back( new geomVector(3) );
        m_surface_Pforce.push_back( vvv );
        m_surface_Vforce.push_back( vvv );
        m_surface_Tgrad.push_back( 0. );
     }
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
void DS_RigidBody:: update_Tgrad_on_surface_point( size_t const& i
                                                 , double const& value )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: update_Tgrad_on_surface_point" ) ;

  m_surface_Tgrad[i] = value;

}




//---------------------------------------------------------------------------
void DS_RigidBody:: correct_surface_discretization( FV_Mesh const* MESH )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: correct_surface_discretization( )" ) ;

  // vector<geomVector> const* p_pbc =
  //                 dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
  //                             ->get_ptr_to_periodic_directions();
  //
  // if (p_pbc) {
     boolVector const* periodic_comp = MESH->get_periodic_directions();
     size_t dim = MESH->nb_space_dimensions() ;

     for (size_t dir=0;dir < dim; dir++) {
        bool is_periodic = periodic_comp->operator()( dir );

        if (is_periodic) {
           double isize = MESH->get_main_domain_max_coordinate(dir)
                        - MESH->get_main_domain_min_coordinate(dir);
           double imin = MESH->get_main_domain_min_coordinate(dir);

           for (size_t i = 0; i < m_surface_area.size(); i++) {
              m_surface_points[i]->operator()(dir) =
                  m_surface_points[i]->operator()(dir)
                - MAC::floor((m_surface_points[i]->operator()(dir)-imin)/isize)
                  * isize;
           }
        }
     }
  // }
}




//---------------------------------------------------------------------------
void DS_RigidBody:: write_surface_discretization( const std::string& file)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: write_surface_discretization" ) ;

  std::ofstream out;

  out.open(file.c_str());
  out << "x ,y ,z ,nx ,ny ,nz ,area ,Fpx ,Fpy ,Fpz ,Fvx ,Fvy ,Fvz, Tgrad" << endl;

  for (size_t i = 0; i < m_surface_area.size(); i++) {
     if ((m_surface_Pforce[i].calcNorm() != 0) ||
         (m_surface_Vforce[i].calcNorm() != 0) ||
         (m_surface_Tgrad[i] != 0))
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
            << m_surface_Vforce[i](2) << " ,"
            << m_surface_Tgrad[i] << endl;
  }

  out.close();

}
