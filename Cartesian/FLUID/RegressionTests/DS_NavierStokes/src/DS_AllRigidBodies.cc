#include <DS_AllRigidBodies.hh>
#include <FS_AllRigidBodies.hh>
#include <DS_RigidBody.hh>
#include <DS_RigidBody_BuilderFactory.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_AllRigidBodies:: DS_AllRigidBodies()
//---------------------------------------------------------------------------
  : m_npart( 0 )
  , m_nrb( 0 )
  , m_FSallrigidbodies( NULL )
{
  MAC_LABEL( "DS_AllRigidBodies:: DS_AllRigidBodies" ) ;

}




//---------------------------------------------------------------------------
DS_AllRigidBodies:: DS_AllRigidBodies( size_t& dimens
                                  , istream& in
                                  , bool const& b_particles_as_fixed_obstacles
                                  , FV_DiscreteField const* arb_UF
                                  , FV_DiscreteField const* arb_PF
                                  , double const& arb_scs)
//---------------------------------------------------------------------------
  : m_space_dimension( dimens )
  , UF ( arb_UF )
  , PF ( arb_PF )
  , surface_cell_scale ( arb_scs )
{
  MAC_LABEL( "DS_AllRigidBodies:: DS_AllRigidBodies(size_t&,istream&)" ) ;

  m_FSallrigidbodies = new FS_AllRigidBodies( m_space_dimension, in,
  	b_particles_as_fixed_obstacles );
  m_nrb = m_FSallrigidbodies->get_number_rigid_bodies();
  m_npart = m_FSallrigidbodies->get_number_particles();
  DS_RigidBody* dsrb = NULL;
  m_allDSrigidbodies.reserve( m_nrb );
  for (size_t i = 0; i < m_nrb; ++i)
  {
    m_allDSrigidbodies.push_back( dsrb );
    m_allDSrigidbodies[i] = DS_RigidBody_BuilderFactory::create(
    	m_FSallrigidbodies->get_ptr_rigid_body(i) );
  }

  build_solid_variables_on_grid();

  initialize_surface_variables_for_all_RB();

  compute_surface_variables_for_all_RB();

  write_surface_discretization_for_all_RB();

  compute_halo_zones_for_all_rigid_body();

  compute_void_fraction_on_grid(PF);
  compute_void_fraction_on_grid(UF);

  compute_grid_intersection_with_rigidbody(PF);
  compute_grid_intersection_with_rigidbody(UF);

}




//---------------------------------------------------------------------------
DS_AllRigidBodies:: ~DS_AllRigidBodies()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: ~DS_AllRigidBodies" ) ;

  for (size_t i = 0; i < m_nrb; ++i) delete m_allDSrigidbodies[i];
  m_allDSrigidbodies.clear();
  if ( m_FSallrigidbodies ) delete m_FSallrigidbodies;

}




//---------------------------------------------------------------------------
size_t DS_AllRigidBodies:: get_number_rigid_bodies() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: get_number_rigid_bodies" ) ;

  return ( m_nrb );

}




//---------------------------------------------------------------------------
size_t DS_AllRigidBodies:: get_number_particles() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: get_number_particles" ) ;

  return ( m_npart );

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: update" ) ;

  m_FSallrigidbodies->update( in );
  for (size_t i = 0; i < m_nrb; ++i) m_allDSrigidbodies[i]->update();

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: display_geometric( ostream& out,
	size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: display_geometric" ) ;

  m_FSallrigidbodies->display( out, indent_width );

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: display( ostream& out,
	size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Features of all Direction Splitting rigid bodies" << endl;
  out << space << three << "Space dimension = " << m_space_dimension << endl;
  out << space << three << "Total number of rigid bodies = " << m_nrb << endl;
  out << space << three << "Number of particles = " << m_npart << endl;
  out << space << three << "Number of obstacles = " << m_nrb - m_npart << endl;
  for (size_t i = 0; i < m_nrb; ++i)
  {
    out << endl;
    out << space << three << "Direction Splitting Rigid body " << i << endl;
    m_allDSrigidbodies[i]->display( out, indent_width + 6 );
  }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_hydro_force_torque( FV_DiscreteField const* PP,
	FV_DiscreteField const* UU )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_hydro_force_torque" ) ;

  for (size_t i = 0; i < m_nrb; ++i)
    m_allDSrigidbodies[i]->compute_hydro_force_torque( PP, UU );

}




//---------------------------------------------------------------------------
bool DS_AllRigidBodies:: isIn( size_t const& parID,
		                         geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn(pt)" ) ;

  return (m_allDSrigidbodies[parID]->isIn( pt ));

}




//---------------------------------------------------------------------------
bool DS_AllRigidBodies:: isIn( size_t const& parID,
		                         double const& x,
                               double const& y,
                               double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn(x,y,z)" ) ;

  return (m_allDSrigidbodies[parID]->isIn( x, y, z ));

}




//---------------------------------------------------------------------------
double DS_AllRigidBodies:: level_set_value( size_t const& parID,
		                         geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn(pt)" ) ;

  return (m_allDSrigidbodies[parID]->level_set_value( pt ));

}




//---------------------------------------------------------------------------
double DS_AllRigidBodies:: level_set_value( size_t const& parID,
		                         double const& x,
                               double const& y,
                               double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn(x,y,z)" ) ;

  return (m_allDSrigidbodies[parID]->level_set_value( x, y, z ));

}




//---------------------------------------------------------------------------
geomVector DS_AllRigidBodies:: rigid_body_velocity( size_t const& parID,
                                             geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: rigid_body_velocity(pt)" ) ;

  return (m_allDSrigidbodies[parID]->get_rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
FS_AllRigidBodies const* DS_AllRigidBodies:: get_ptr_FS_AllRigidBodies() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: get_ptr_FS_AllRigidBodies" ) ;

  return ( m_FSallrigidbodies );

}




//---------------------------------------------------------------------------
DS_RigidBody* DS_AllRigidBodies:: get_ptr_rigid_body( size_t i )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: get_ptr_rigid_body" ) ;

  return ( m_allDSrigidbodies[i] );

}




//---------------------------------------------------------------------------
DS_RigidBody const* DS_AllRigidBodies:: get_ptr_rigid_body( size_t i ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: get_ptr_rigid_body" ) ;

  return ( m_allDSrigidbodies[i] );

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_void_fraction_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_void_fraction_on_grid" ) ;

  size_t nb_comps = FF->nb_components() ;
  size_t field = (FF == PF) ? 0 : 1 ;
  size_t dim = FF->primary_grid()->nb_space_dimensions() ;

  boolVector const* periodic_comp =
                        FF->primary_grid()->get_periodic_directions();

  // Get local min and max indices;
  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  size_t i0_temp = 0;

  for (size_t parID = 0; parID < m_nrb; ++parID) {
     vector<geomVector*> haloZone = m_allDSrigidbodies[parID]
                                                   ->get_rigid_body_haloZone();
     // Calculation on the indexes near the rigid body
     for (size_t comp = 0; comp < nb_comps; ++comp) {
        for (size_t dir = 0; dir < dim; ++dir) {
           // Calculations for solids on the total unknown on the proc
           min_unknown_index(dir) =
                   FF->get_min_index_unknown_on_proc( comp, dir );
           max_unknown_index(dir) =
                   FF->get_max_index_unknown_on_proc( comp, dir );

           bool is_periodic = periodic_comp->operator()( dir );
           double domain_min =
                   FF->primary_grid()->get_main_domain_min_coordinate( dir );
           double domain_max =
                   FF->primary_grid()->get_main_domain_max_coordinate( dir );

           bool found =FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,dir)
                                       , haloZone[0]->operator()(dir)
                                       , i0_temp) ;
           size_t index_min = (found) ? i0_temp : min_unknown_index(dir);


           found = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,dir)
                                   , haloZone[1]->operator()(dir)
                                   , i0_temp) ;
           size_t index_max = (found) ? i0_temp : max_unknown_index(dir);

           if (is_periodic &&
              ((haloZone[1]->operator()(dir) > domain_max)
            || (haloZone[0]->operator()(dir) < domain_min))) {
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

                 if (m_allDSrigidbodies[parID]->isIn(xC,yC,zC))
                    void_fraction[field]->operator()(p) = 1 + parID;
             }
           }
        }
     }
  }
}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_halo_zones_for_all_rigid_body( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_halo_zones_for_all_rigid_body" ) ;

  for (size_t i = 0; i < m_nrb; ++i) {
     m_allDSrigidbodies[i]->compute_rigid_body_halozone();
  }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_grid_intersection_with_rigidbody(
                                                   FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_grid_intersection_with_rigidbody" ) ;

  size_t nb_comps = FF->nb_components() ;
  size_t dim = FF->primary_grid()->nb_space_dimensions() ;
  size_t field = (FF == PF) ? 0 : 1 ;

  boolVector const* periodic_comp =
                        FF->primary_grid()->get_periodic_directions();

  for (size_t parID = 0; parID < m_nrb; ++parID) {
     vector<geomVector*> haloZone = m_allDSrigidbodies[parID]
                                                   ->get_rigid_body_haloZone();
     size_t_vector min_unknown_index(3,0);
     size_t_vector max_unknown_index(3,0);
     size_t_vector min_nearest_index(3,0);
     size_t_vector max_nearest_index(3,0);
     size_t_vector ipos(3,0);
     size_t_array2D local_extents(3,2,0);
     size_t i0_temp = 0;

     double delta = FF->primary_grid()->get_smallest_grid_size();

     for (size_t comp = 0; comp < nb_comps; comp++) {

        for (size_t dir = 0; dir < dim; dir++) {
          // To include knowns at dirichlet boundary in the intersection
          // calculation as well, modification to the looping extents are required
          min_unknown_index(dir) =
                              FF->get_min_index_unknown_on_proc( comp, dir );
          max_unknown_index(dir) =
                              FF->get_max_index_unknown_on_proc( comp, dir );
          local_extents(dir,0) = 0;
          local_extents(dir,1) = max_unknown_index(dir)
                               - min_unknown_index(dir);

          bool is_periodic = periodic_comp->operator()( dir );

          double domain_min =
                  FF->primary_grid()->get_main_domain_min_coordinate(dir);
          double domain_max =
                  FF->primary_grid()->get_main_domain_max_coordinate(dir);
          bool found =
                  FV_Mesh::between(FF->get_DOF_coordinates_vector( comp, dir)
                                 , haloZone[0]->operator()(dir)
                                 , i0_temp) ;
          size_t index_min = (found) ? i0_temp : min_unknown_index(dir);
          found = FV_Mesh::between(FF->get_DOF_coordinates_vector( comp, dir )
                                  , haloZone[1]->operator()(dir)
                                  , i0_temp) ;
          size_t index_max = (found) ? i0_temp : max_unknown_index(dir);

          if (is_periodic &&
              ((haloZone[1]->operator()(dir) > domain_max)
          || (haloZone[0]->operator()(dir) < domain_min))) {
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
                 double zC = (dim == 2) ? 0
                                        : FF->get_DOF_coordinate( k, comp, 2 );

                 size_t p = FF->DOF_local_number(i,j,k,comp);
                 geomVector source(xC,yC,zC);

                 if (void_fraction[field]->operator()(p) == 0) {
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

                            if (void_fraction[field]->operator()(neigh_num)
                                                                  == parID+1) {
                               double t = m_allDSrigidbodies[parID]
                                    ->get_distanceTo( source, rayDir, delta );
                               // Storing the direction with RB intersection
                               intersect_vector[field]->operator()(p,col) = 1;
                               // Storing the intersection distance
                               intersect_distance[field]->operator()(p,col) = t;
                               // Calculate the variable values on the
                               // intersection of grid and solid
                               geomVector netVel =
                                rigid_body_velocity(parID, source + t * rayDir);
                               // Value of variable at the surface of particle
                               if ( nb_comps == 1) { // i.e. PF
                                  intersect_fieldValue[field]->operator()(p,col)
                                                                 = netVel(dir);
                               } else { // i.e. UF
                                  intersect_fieldValue[field]->operator()(p,col)
                                                                 = netVel(comp);
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

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: initialize_surface_variables_for_all_RB( )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_AllRigidBodies:: initialize_surface_variables_for_all_RB" ) ;

   double dx = UF->primary_grid()->get_smallest_grid_size();

   for (size_t i = 0; i < m_nrb; ++i) {
      m_allDSrigidbodies[i]->initialize_surface_variables(surface_cell_scale
                                                        , dx);
   }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_surface_variables_for_all_RB( )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_AllRigidBodies:: compute_surface_variables_for_all_RB" ) ;

   for (size_t i = 0; i < m_nrb; ++i)
      m_allDSrigidbodies[i]->compute_surface_points( );

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: write_surface_discretization_for_all_RB( )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_AllRigidBodies:: write_surface_discretization_for_all_RB" ) ;

   for (size_t i = 0; i < m_nrb; ++i) {
      std::ostringstream os2;
      os2 << "./DS_results/discretized_surface_parID_" << i << ".csv";
      std::string file = os2.str();
      m_allDSrigidbodies[i]->write_surface_discretization( file );
   }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: build_solid_variables_on_grid(  )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_AllRigidBodies:: build_solid_variables_on_grid" ) ;

   size_t PF_LOC_UNK = PF->nb_local_unknowns();
   size_t UF_LOC_UNK = UF->nb_local_unknowns();

   // void fraction on the computational grid
   // For PF and UF
   void_fraction.reserve(2);
   void_fraction.push_back(new size_t_vector(1,0));
   void_fraction.push_back(new size_t_vector(1,0));
   void_fraction[0]->re_initialize(PF_LOC_UNK);
   void_fraction[1]->re_initialize(UF_LOC_UNK);

   // Intersection parameters on the computational grid
   // For PF and UF
   intersect_vector.reserve(2);
   intersect_vector.push_back(new size_t_array2D(1,1,0));
   intersect_vector.push_back(new size_t_array2D(1,1,0));
   intersect_vector[0]->re_initialize(PF_LOC_UNK,6);
   intersect_vector[1]->re_initialize(UF_LOC_UNK,6);
   intersect_distance.reserve(2);
   intersect_distance.push_back(new doubleArray2D(1,1,0.));
   intersect_distance.push_back(new doubleArray2D(1,1,0.));
   intersect_distance[0]->re_initialize(PF_LOC_UNK,6);
   intersect_distance[1]->re_initialize(UF_LOC_UNK,6);
   intersect_fieldValue.reserve(2);
   intersect_fieldValue.push_back(new doubleArray2D(1,1,0.));
   intersect_fieldValue.push_back(new doubleArray2D(1,1,0.));
   intersect_fieldValue[0]->re_initialize(PF_LOC_UNK,6);
   intersect_fieldValue[1]->re_initialize(UF_LOC_UNK,6);

}




//---------------------------------------------------------------------------
size_t_vector* DS_AllRigidBodies:: get_void_fraction_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = (FF == PF) ? 0 : 1;

  return (void_fraction[field]);

}




//---------------------------------------------------------------------------
size_t_array2D* DS_AllRigidBodies:: get_intersect_vector_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = (FF == PF) ? 0 : 1;

  return (intersect_vector[field]);

}




//---------------------------------------------------------------------------
doubleArray2D* DS_AllRigidBodies:: get_intersect_distance_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = (FF == PF) ? 0 : 1;

  return (intersect_distance[field]);

}




//---------------------------------------------------------------------------
doubleArray2D* DS_AllRigidBodies:: get_intersect_fieldValue_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = (FF == PF) ? 0 : 1;

  return (intersect_fieldValue[field]);

}
