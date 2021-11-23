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
void DS_RigidBody:: compute_void_fraction_on_grid( FV_DiscreteField const* FF
                                                 , size_t_vector* void_fraction
                                                 , size_t_vector* rb_ID
                                                 , size_t const& parID )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBodies:: compute_void_fraction_on_grid" ) ;

  size_t nb_comps = FF->nb_components() ;
  size_t dim = FF->primary_grid()->nb_space_dimensions() ;

  boolVector const* periodic_comp =
                        FF->primary_grid()->get_periodic_directions();

  // Get local min and max indices;
  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  size_t i0_temp = 0;

  // Calculation on the indexes near the rigid body
  for (size_t comp = 0; comp < nb_comps; ++comp) {
     for (size_t dir = 0; dir < dim; ++dir) {
        // Calculations for solids on the total unknown on the proc
        min_unknown_index(dir) = FF->get_min_index_unknown_on_proc( comp, dir );
        max_unknown_index(dir) = FF->get_max_index_unknown_on_proc( comp, dir );

        bool is_periodic = periodic_comp->operator()( dir );
        double domain_min =
                FF->primary_grid()->get_main_domain_min_coordinate( dir );
        double domain_max =
                FF->primary_grid()->get_main_domain_max_coordinate( dir );

        bool found = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,dir)
                                    , m_halo_zone[0]->operator()(dir)
                                    , i0_temp) ;
        size_t index_min = (found) ? i0_temp : min_unknown_index(dir);


        found = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,dir)
                                , m_halo_zone[1]->operator()(dir)
                                , i0_temp) ;
        size_t index_max = (found) ? i0_temp : max_unknown_index(dir);

        if (is_periodic &&
           ((m_halo_zone[1]->operator()(dir) > domain_max)
         || (m_halo_zone[0]->operator()(dir) < domain_min))) {
           index_min = min_unknown_index(dir);
           index_max = max_unknown_index(dir);
        }

        min_unknown_index(dir) = MAC::max(min_unknown_index(dir),index_min);
        max_unknown_index(dir) = MAC::min(max_unknown_index(dir),index_max);

     }

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
        for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
           for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
              double zC = (dim == 2) ? 0.
                                  : FF->get_DOF_coordinate( k, comp, 2 ) ;
              size_t p = FF->DOF_local_number(i,j,k,comp);

              if (isIn(xC,yC,zC)) {
                 void_fraction->operator()(p) = 1 + parID;
                 rb_ID->operator()(p) = parID;
              }
           }
        }
     }
 }

}




//---------------------------------------------------------------------------
void
DS_RigidBody:: compute_grid_intersection_with_rigidbody (
                                          FV_DiscreteField const* FF
                                        , size_t_vector const* void_fraction
                                        , size_t_array2D* intersect_vector
                                        , doubleArray2D* intersect_distance
                                        , doubleArray2D* intersect_fieldValue
                                        , size_t const& parID)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: compute_grid_intersection_with_rigidbody" ) ;

  size_t nb_comps = FF->nb_components() ;
  size_t dim = FF->primary_grid()->nb_space_dimensions() ;

  boolVector const* periodic_comp =
                        FF->primary_grid()->get_periodic_directions();

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);
  size_t_vector min_nearest_index(3,0);
  size_t_vector max_nearest_index(3,0);
  size_t_vector ipos(3,0);
  size_t_array2D local_extents(3,2,0);
  vector<double> net_vel(3,0.);
  size_t i0_temp = 0;

  double delta = FF->primary_grid()->get_smallest_grid_size();

  for (size_t comp = 0; comp < nb_comps; comp++) {

     for (size_t dir = 0; dir < dim; dir++) {
        // To include knowns at dirichlet boundary in the intersection calculation
        // as well, modification to the looping extents are required
        min_unknown_index(dir) = FF->get_min_index_unknown_on_proc( comp, dir ) ;
        max_unknown_index(dir) = FF->get_max_index_unknown_on_proc( comp, dir ) ;
        local_extents(dir,0) = 0;
        local_extents(dir,1) = max_unknown_index(dir) - min_unknown_index(dir);

        bool is_periodic = periodic_comp->operator()( dir );

        double domain_min =
               FF->primary_grid()->get_main_domain_min_coordinate(dir);
        double domain_max =
               FF->primary_grid()->get_main_domain_max_coordinate(dir);
        bool found = FV_Mesh::between(FF->get_DOF_coordinates_vector( comp, dir )
                                    , m_halo_zone[0]->operator()(dir)
                                    , i0_temp) ;
        size_t index_min = (found) ? i0_temp : min_unknown_index(dir);
        found = FV_Mesh::between(FF->get_DOF_coordinates_vector( comp, dir )
                               , m_halo_zone[1]->operator()(dir)
                               , i0_temp) ;
        size_t index_max = (found) ? i0_temp : max_unknown_index(dir);

        if (is_periodic &&
           ((m_halo_zone[1]->operator()(dir) > domain_max)
        || (m_halo_zone[0]->operator()(dir) < domain_min))) {
           index_min = min_unknown_index(dir);
           index_max = max_unknown_index(dir);
        }

        min_nearest_index(dir) = MAC::max(min_unknown_index(dir),index_min);
        max_nearest_index(dir) = MAC::min(max_unknown_index(dir),index_max);

     }

     max_nearest_index(2) = (dim == 2) ? 1 : max_nearest_index(2);

     for (size_t i = min_nearest_index(0); i < max_nearest_index(0); ++i) {
        ipos(0) = i - min_unknown_index(0);
        double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
        for (size_t j = min_nearest_index(1); j < max_nearest_index(1); ++j) {
           ipos(1) = j - min_unknown_index(1);
           double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
           for (size_t k = min_nearest_index(2);k < max_nearest_index(2); ++k) {
              ipos(2) = k - min_unknown_index(2);
              double zC = (dim == 2) ? 0 : FF->get_DOF_coordinate( k, comp, 2 );

              size_t p = FF->DOF_local_number(i,j,k,comp);
              geomVector source(xC,yC,zC);

              if (void_fraction->operator()(p) == 0) {
                 for (size_t dir = 0; dir < dim; dir++) {
                    for (size_t off = 0; off < 2; off++) {
                       size_t col = 2*dir + off;

                       geomVector rayDir(0.,0.,0.);
                       rayDir(dir) = (off == 0) ? -1 : 1 ;

            		     // Checking if the nodes are on domain boundary or not,
            		     // if so, the check the intersection only on one side
            		     if (ipos(dir) != local_extents(dir,off)) {
                          geomVector ineigh((double)i,(double)j,(double)k);

                          ineigh += rayDir;

                          size_t neigh_num = FF->DOF_local_number(
                           (size_t)ineigh(0),(size_t)ineigh(1),(size_t)ineigh(2)
                                                                         ,comp);

                          if (void_fraction->operator()(neigh_num) == parID+1) {
                             double t = m_geometric_rigid_body
                                          ->distanceTo( source, rayDir, delta );
                             // Storing the direction with RB intersection
                             intersect_vector->operator()(p,col) = 1;
                             // Storing the intersection distance
                             intersect_distance->operator()(p,col) = t;
                             // Calculate the variable values on the
                             // intersection of grid and solid
                             geomVector netVel = rigid_body_velocity(source + t * rayDir);
                             // Value of variable at the surface of particle
                             if ( nb_comps == 1) { // i.e. PF
                                intersect_fieldValue->operator()(p,col)
                                                                 = net_vel[dir];
                             } else { // i.e. UF
                                intersect_fieldValue->operator()(p,col)
                                                                 = net_vel[comp];
                             }
                          }
                       }
                    }
   	           }
              }
           }
        }
     }
  }
}




//---------------------------------------------------------------------------
geomVector DS_RigidBody:: rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}
