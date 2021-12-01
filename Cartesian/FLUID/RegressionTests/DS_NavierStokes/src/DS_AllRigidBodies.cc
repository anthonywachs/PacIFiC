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


  compute_halo_zones_for_all_rigid_body();

  compute_void_fraction_on_grid(PF);
  compute_void_fraction_on_grid(UF);

  compute_grid_intersection_with_rigidbody(PF);
  compute_grid_intersection_with_rigidbody(UF);

  compute_pressure_force_and_torque_for_allRB( );

  write_surface_discretization_for_all_RB();
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
void DS_AllRigidBodies:: compute_pressure_force_and_torque_for_allRB( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_pressure_force_and_torque" ) ;

  for (size_t i = 0; i < m_nrb; ++i)
    first_order_pressure_stress(i);

}




//---------------------------------------------------------------------------
int DS_AllRigidBodies:: isIn_any_RB( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn_any_RB(pt)" ) ;

  for (size_t i = 0; i < m_nrb; ++i)
     if (m_allDSrigidbodies[i]->isIn( pt )) return ((int)i);

  return (-1);

}




//---------------------------------------------------------------------------
int DS_AllRigidBodies:: isIn_any_RB( double const& x,
                                     double const& y,
                                     double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn_any_RB(x,y,z)" ) ;

  for (size_t i = 0; i < m_nrb; ++i)
     if (m_allDSrigidbodies[i]->isIn( x, y, z )) return ((int)i);

  return (-1);

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
  MAC_LABEL( "DS_AllRigidBodies:: level_set_value" ) ;

  return (m_allDSrigidbodies[parID]->level_set_value( pt ));

}




//---------------------------------------------------------------------------
double DS_AllRigidBodies:: level_set_value( size_t const& parID,
		                         double const& x,
                               double const& y,
                               double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: level_set_value" ) ;

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
        for (size_t dir = 0; dir < m_space_dimension; ++dir) {
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
                 double zC = (m_space_dimension == 2) ? 0.
                                     : FF->get_DOF_coordinate( k, comp, 2 ) ;
                 size_t p = FF->DOF_local_number(i,j,k,comp);

                 if (isIn(parID,xC,yC,zC))
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

        for (size_t dir = 0; dir < m_space_dimension; dir++) {
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

        max_nearest_index(2) = (m_space_dimension == 2) ? 1
                                                        : max_nearest_index(2);

        for (size_t i = min_nearest_index(0); i < max_nearest_index(0); ++i) {
          ipos(0) = i - min_unknown_index(0);
          double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
          for (size_t j = min_nearest_index(1); j < max_nearest_index(1); ++j) {
              ipos(1) = j - min_unknown_index(1);
              double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
              for (size_t k = min_nearest_index(2);k < max_nearest_index(2); ++k) {
                 ipos(2) = k - min_unknown_index(2);
                 double zC = (m_space_dimension == 2) ? 0
                                        : FF->get_DOF_coordinate( k, comp, 2 );

                 size_t p = FF->DOF_local_number(i,j,k,comp);
                 geomVector source(xC,yC,zC);

                 if (void_fraction[field]->operator()(p) == 0) {
                    for (size_t dir = 0; dir < m_space_dimension; dir++) {
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
      m_allDSrigidbodies[i]->compute_number_of_surface_variables(
                                          surface_cell_scale, dx);
      m_allDSrigidbodies[i]->initialize_surface_variables( );
   }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: first_order_pressure_stress( size_t const& parID )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllRigidBodies:: first_order_pressure_stress" ) ;

  boolVector const* periodic_comp =
                        PF->primary_grid()->get_periodic_directions();

  size_t comp = 0;

  size_t_vector min_unknown_index(m_space_dimension,0);
  size_t_vector max_unknown_index(m_space_dimension,0);
  // Domain length and minimum
  geomVector domain_length(3), domain_min(3);
  // Extents on the currect processor
  geomVector Dmin(3), Dmax(3);
  vector<geomVector*> surface_point = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_points();
  vector<geomVector*> surface_normal = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_normals();
  vector<geomVector*> surface_area = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_areas();


  // Get local min and max indices
  // One extra grid cell needs to considered, since ghost points can be
  // located in between the min/max index handled by the proc
  for (size_t l = 0; l < m_space_dimension; l++) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( comp, l );
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( comp, l );
     Dmin(l) = PF->primary_grid()->get_min_coordinate_on_current_processor( l );
     Dmax(l) = PF->primary_grid()->get_max_coordinate_on_current_processor( l );
     domain_length(l) = PF->primary_grid()->get_main_domain_max_coordinate( l )
                      - PF->primary_grid()->get_main_domain_min_coordinate( l );
     domain_min(l) = PF->primary_grid()->get_main_domain_min_coordinate( l );
  }

  for (size_t i = 0; i < surface_area.size(); i++) {
     double stress = 0.;
   	// Displacement correction in case of periodic
      // boundary condition in any direction
  //    for (size_t dir=0;dir<m_space_dimension;dir++) {
  //       bool is_periodic = periodic_comp->operator()( dir );
  //       if (is_periodic)
  //          surface_point[i]->operator()(dir) = surface_point[i]->operator()(dir)
  //                                - MAC::floor((surface_point[i]->operator()(dir)
  //                                                             - domain_min(dir))
  //                                   / domain_length(dir))*domain_length(dir);
  //    }
  //
     // Check it the point is in the current domain
     bool status = (surface_point[i]->operator()(0) > Dmin(0))
                && (surface_point[i]->operator()(0) <= Dmax(0))
                && (surface_point[i]->operator()(1) > Dmin(1))
                && (surface_point[i]->operator()(1) <= Dmax(1))
                && (surface_point[i]->operator()(2) > Dmin(2))
                && (surface_point[i]->operator()(2) <= Dmax(2));

     if (status) {
        // Finding the grid indexes next to ghost point
        size_t i0_temp;
        size_t_vector i0(3,0);
        bool found = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,0)
                                 , surface_point[i]->operator()(0), i0_temp);
        if (found == 1) i0(0) = i0_temp;

        found = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,1)
                                 , surface_point[i]->operator()(1), i0_temp);
        if (found == 1) i0(1) = i0_temp;

        if (m_space_dimension == 3) {
           found = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,2)
                                 , surface_point[i]->operator()(2), i0_temp);
           if (found == 1) i0(2) = i0_temp;
        }
        // Calculation of field variable on ghost point
        size_t_vector face_vector(3,0);
        face_vector(0) = 1; face_vector(1) = 1; face_vector(2) = 0;
        double press = (m_space_dimension == 2) ?
         Bilinear_interpolation(PF,comp,surface_point[i],i0,face_vector,{0,1}) :
         Trilinear_interpolation(PF,comp,surface_point[i],i0,{0,1}) ;
        stress = - press/2.;
     }

     // Ref: Keating thesis Pg-85
     // point_coord*(area) --> Component of area in particular direction
     geomVector value(3);
     double norm = surface_normal[i]->calcNorm();
     value(0) = stress*surface_normal[i]->operator()(0)/norm
                      *surface_area[i]->operator()(0);
     value(1) = stress*surface_normal[i]->operator()(1)/norm
                      *surface_area[i]->operator()(0);
     value(2) = stress*surface_normal[i]->operator()(2)/norm
                      *surface_area[i]->operator()(0);

     m_allDSrigidbodies[parID]->update_Pforce_on_surface_point(i,value);

     pressure_force[parID] += value;

     pressure_torque[parID] += surface_point[i]->operator^(value);
  }

}




//---------------------------------------------------------------------------
double DS_AllRigidBodies:: Trilinear_interpolation ( FV_DiscreteField const* FF
                                                   , size_t const& comp
                                                   , geomVector const* pt
                                                   , size_t_vector const& i0
                                                   , vector<size_t> const& list)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_AllRigidBodies:: Trilinear_interpolation" ) ;

   vector<size_t_vector> face_vec;
   face_vec.reserve(3);
   size_t_vector vvv(3,0);
   face_vec.push_back(vvv);
   face_vec.push_back(vvv);
   face_vec.push_back(vvv);
   face_vec[0](0) = 0; face_vec[0](1) = 1; face_vec[0](2) = 1;
   face_vec[1](0) = 1; face_vec[1](1) = 0; face_vec[1](2) = 1;
   face_vec[2](0) = 1; face_vec[2](1) = 1; face_vec[2](2) = 0;
   double dh = FF->primary_grid()->get_smallest_grid_size();
   double value = 0.;
   geomVector netVel(3);
   geomVector pt_at_face(3);
   int in_RB = -1;

   for (size_t face_norm = 0; face_norm < 3; face_norm++) {
      geomVector rayDir(3);
      double dl = 0., dr = 0.;
      double fl = 0., fr = 0.;

      // Left face
      size_t_vector i0_left = i0;
      pt_at_face = *pt;
      pt_at_face(face_norm) = FF->get_DOF_coordinate(i0_left(face_norm)
                                                   , comp
                                                   , face_norm);
      rayDir(face_norm) = -1.;

      in_RB = (FF == UF) ? isIn_any_RB(pt_at_face) : -1 ;

      if (in_RB != -1) {
         dl = m_allDSrigidbodies[in_RB]->get_distanceTo(*pt, rayDir, dh);
         netVel = rigid_body_velocity(in_RB, *pt + dl * rayDir);
         fl = netVel(comp);
      } else {
         fl = Bilinear_interpolation(FF, comp, &pt_at_face, i0_left
                                             , face_vec[face_norm], list);
         dl = MAC::abs(pt_at_face(face_norm) - pt->operator()(face_norm));
      }

      // Right face
      size_t_vector i0_right = i0;
      i0_right(face_norm) += 1;
      pt_at_face = *pt;
      pt_at_face(face_norm) = FF->get_DOF_coordinate(i0_right(face_norm)
                                                   , comp
                                                   , face_norm);
      rayDir(face_norm) = 1.;

      in_RB = (FF == UF) ? isIn_any_RB(pt_at_face) : -1 ;

      if (in_RB != -1) {
         dr = m_allDSrigidbodies[in_RB]->get_distanceTo(*pt, rayDir, dh);
         netVel = rigid_body_velocity(in_RB, *pt + dr * rayDir);
         fr = netVel(comp);
      } else {
         fr = Bilinear_interpolation(FF, comp, &pt_at_face, i0_right
                                             , face_vec[face_norm], list);
         dr = MAC::abs(pt_at_face(face_norm) - pt->operator()(face_norm));
      }

      value += (dr*fl + dl*fr)/(dl + dr);

   }

   value *= (1./3.);

   return(value);
}




//---------------------------------------------------------------------------
double DS_AllRigidBodies:: Bilinear_interpolation ( FV_DiscreteField const* FF
                                                , size_t const& comp
                                                , geomVector const* pt
                                                , size_t_vector const& i0
                                                , size_t_vector const& face_vec
                                                , vector<size_t> const& list)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_RigidBody:: Bilinear_interpolation" ) ;

   size_t field = (FF == PF) ? 0 : 1;
   double dh = FF->primary_grid()->get_smallest_grid_size();


   size_t_array2D p(2,2,0);
   // arrays of vertex indexes of the cube/square
   size_t_array2D ix(2,2,0);
   size_t_array2D iy(2,2,0);
   size_t_array2D iz(2,2,0);
   // Min and max of the cell containing ghost point
   doubleArray2D extents(m_space_dimension,2,0);
   // Field value at vertex of face
   doubleArray2D f(2,2,0);
   // Interpolated values at the walls
   doubleArray2D fwall(2,2,0);
   // Distance of grid/particle walls from the ghost point
   doubleArray2D del_wall(2,2,0);

   // Direction in the plane of face
   size_t dir1=0, dir2=0;

   if (face_vec(0) == 0) {
      ix(0,0) = i0(0),     iy(0,0) = i0(1),    iz(0,0) = i0(2);
      ix(1,0) = i0(0),     iy(1,0) = i0(1),    iz(1,0) = i0(2)+1;
      ix(0,1) = i0(0),     iy(0,1) = i0(1)+1,  iz(0,1) = i0(2);
      ix(1,1) = i0(0),     iy(1,1) = i0(1)+1,  iz(1,1) = i0(2)+1;
      dir1 = 2, dir2 = 1;
   } else if (face_vec(1) == 0) {
      ix(0,0) = i0(0),     iy(0,0) = i0(1),    iz(0,0) = i0(2);
      ix(1,0) = i0(0)+1,   iy(1,0) = i0(1),    iz(1,0) = i0(2);
      ix(0,1) = i0(0),     iy(0,1) = i0(1),    iz(0,1) = i0(2)+1;
      ix(1,1) = i0(0)+1,   iy(1,1) = i0(1),    iz(1,1) = i0(2)+1;
      dir1 = 0, dir2 = 2;
   } else if (face_vec(2) == 0) {
      ix(0,0) = i0(0),     iy(0,0) = i0(1),    iz(0,0) = i0(2);
      ix(1,0) = i0(0)+1,   iy(1,0) = i0(1),    iz(1,0) = i0(2);
      ix(0,1) = i0(0),     iy(0,1) = i0(1)+1,  iz(0,1) = i0(2);
      ix(1,1) = i0(0)+1,   iy(1,1) = i0(1)+1,  iz(1,1) = i0(2);
      dir1 = 0, dir2 = 1;
   }

   for (size_t i=0; i < 2; i++) {
      for (size_t j=0; j < 2; j++) {
         // Face vertex index
         p(i,j) = FF->DOF_local_number( ix(i,j), iy(i,j), iz(i,j), comp);
         // Vertex field values
         for (size_t level : list)
            f(i,j) += FF->DOF_value( ix(i,j), iy(i,j), iz(i,j), comp, level );
      }
   }

   // Min and max coordinate in the grid cell
   for (size_t dir = 0; dir < m_space_dimension; dir++) {
      extents(dir,0) = FF->get_DOF_coordinate(i0(dir),comp,dir) ;
      extents(dir,1) = FF->get_DOF_coordinate(i0(dir)+face_vec(dir),comp,dir) ;
   }

   // Contribution from left and right wall
   for (size_t i = 0; i < 2; i++) {     // 0 --> left; 1 --> right
      size_t col_top = 2*dir2 + 1;
      size_t col_bot = 2*dir2 + 0;
      if ((field == 0) || ((void_fraction[field]->operator()(p(i,0)) == 0)
                        && (void_fraction[field]->operator()(p(i,1)) == 0))) {
         fwall(0,i) = ((extents(dir2,1) - pt->operator()(dir2))*f(i,0)
                     + (pt->operator()(dir2) - extents(dir2,0))*f(i,1))
                     / (extents(dir2,1)-extents(dir2,0));
         del_wall(0,i) = MAC::abs(extents(dir1,i) - pt->operator()(dir1));
      // if bottom vertex is in fluid domain
      } else if ((void_fraction[field]->operator()(p(i,0)) == 0)
           && (intersect_vector[field]->operator()(p(i,0),col_top) == 1)) {
         double yint = intersect_distance[field]->operator()(p(i,0),col_top);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (yint >= (pt->operator()(dir2)-extents(dir2,0))) {
            fwall(0,i) = ((extents(dir2,0) + yint - pt->operator()(dir2))*f(i,0)
                     + (pt->operator()(dir2) - extents(dir2,0))
                     *intersect_fieldValue[field]->operator()(p(i,0),col_top))
                     / yint;
            del_wall(0,i) = MAC::abs(extents(dir1,i) - pt->operator()(dir1));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction[field]->operator()(p(i,1)) - 1;
            geomVector rayDir(3);
            rayDir(dir1) = (i == 0) ? -1 : 1 ;
            del_wall(0,i) = m_allDSrigidbodies[id]
                                       ->get_distanceTo( *pt, rayDir, dh );
            geomVector surface_point = *pt + del_wall(0,i)*rayDir;
            geomVector net_vel = rigid_body_velocity(id,surface_point);
            fwall(0,i) = net_vel(comp);
         }
      // if top vertex is in fluid domain
      } else if ((void_fraction[field]->operator()(p(i,1)) == 0)
           && (intersect_vector[field]->operator()(p(i,1),col_bot) == 1)) {
         double yint = intersect_distance[field]->operator()(p(i,1),col_bot);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (yint >= (extents(dir2,1)-pt->operator()(dir2))) {
            fwall(0,i) = ((pt->operator()(dir2) + yint - extents(dir2,1))*f(i,1)
                     + (extents(dir2,1) - pt->operator()(dir2))
                     *intersect_fieldValue[field]->operator()(p(i,1),col_bot))
                     / yint;
            del_wall(0,i) = MAC::abs(extents(dir1,i) - pt->operator()(dir1));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction[field]->operator()(p(i,0)) - 1;
            geomVector rayDir(3);
            rayDir(dir1) = (i == 0) ? -1 : 1 ;
            del_wall(0,i) = m_allDSrigidbodies[id]
                                       ->get_distanceTo( *pt, rayDir, dh );
            geomVector surface_point = *pt + del_wall(0,i)*rayDir;
            geomVector net_vel = rigid_body_velocity(id,surface_point);
            fwall(0,i) = net_vel(comp);
         }
      // if both vertex's are in solid domain
      } else if ((void_fraction[field]->operator()(p(i,0)) != 0)
              && (void_fraction[field]->operator()(p(i,1)) != 0)) {

         size_t id = void_fraction[field]->operator()(p(i,0)) - 1;
         geomVector rayDir(3);
         rayDir(dir1) = (i == 0) ? -1 : 1 ;
         del_wall(0,i) = m_allDSrigidbodies[id]
                                       ->get_distanceTo( *pt, rayDir, dh );
         geomVector surface_point = *pt + del_wall(0,i)*rayDir;
         geomVector net_vel = rigid_body_velocity(id,surface_point);
         fwall(0,i) = net_vel(comp);
      }
   }

   // Contribution from top and bottom wall
   for (size_t j = 0; j < 2; j++) {         // 0 --> bottom; 1 --> top
      size_t col_right = 2*dir1 + 1;
      size_t col_left = 2*dir1 + 0;
      if ((field == 0) || ((void_fraction[field]->operator()(p(0,j)) == 0)
                        && (void_fraction[field]->operator()(p(1,j)) == 0))) {
         fwall(1,j) = ((extents(dir1,1) - pt->operator()(dir1)) * f(0,j)
                     + (pt->operator()(dir1) - extents(dir1,0)) * f(1,j))
                     / (extents(dir1,1) - extents(dir1,0));
         del_wall(1,j) = MAC::abs(extents(dir2,j) - pt->operator()(dir2));
      // if left vertex is in fluid domain
      } else if ((void_fraction[field]->operator()(p(0,j)) == 0)
           && (intersect_vector[field]->operator()(p(0,j),col_right) == 1)) {
         double xint = intersect_distance[field]->operator()(p(0,j),col_right);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (xint >= (pt->operator()(dir1) - extents(dir1,0))) {
            fwall(1,j) = ((extents(dir1,0) + xint - pt->operator()(dir1)) * f(0,j)
                  + (pt->operator()(dir1) - extents(dir1,0))
                  * intersect_fieldValue[field]->operator()(p(0,j),col_right))
                  / xint;
            del_wall(1,j) = MAC::abs(extents(dir2,j) - pt->operator()(dir2));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction[field]->operator()(p(1,j)) - 1;
            geomVector rayDir(3);
            rayDir(dir2) = (j == 0) ? -1 : 1 ;
            del_wall(1,j) = m_allDSrigidbodies[id]
                                       ->get_distanceTo( *pt, rayDir, dh );
            geomVector surface_point = *pt + del_wall(1,j)*rayDir;
            geomVector net_vel = rigid_body_velocity(id,surface_point);
            fwall(1,j) = net_vel(comp);
         }
      // if right vertex is in fluid domain
      } else if ((void_fraction[field]->operator()(p(1,j)) == 0)
           && (intersect_vector[field]->operator()(p(1,j),col_left) == 1)) {

         double xint = intersect_distance[field]->operator()(p(1,j),col_left);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (xint >= (extents(dir1,1) - pt->operator()(dir1))) {
            fwall(1,j) = ((pt->operator()(dir1) + xint - extents(dir1,1)) * f(1,j)
                     + (extents(dir1,1) - pt->operator()(dir1))
                     * intersect_fieldValue[field]->operator()(p(1,j),col_left))
                     / xint;
            del_wall(1,j) = MAC::abs(extents(dir2,j) - pt->operator()(dir2));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction[field]->operator()(p(0,j)) - 1;
            geomVector rayDir(3);
            rayDir(dir2) = (j == 0) ? -1 : 1 ;
            del_wall(1,j) = m_allDSrigidbodies[id]
                                       ->get_distanceTo( *pt, rayDir, dh );
            geomVector surface_point = *pt + del_wall(1,j)*rayDir;
            geomVector net_vel = rigid_body_velocity(id,surface_point);
            fwall(1,j) = net_vel(comp);
         }
      // if both vertex's are in solid domain
   } else if ((void_fraction[field]->operator()(p(0,j)) != 0)
              && (void_fraction[field]->operator()(p(1,j)) != 0)) {
         size_t id = void_fraction[field]->operator()(p(0,j)) - 1 ;
         geomVector rayDir(3);
         rayDir(dir2) = (j == 0) ? -1 : 1 ;
         del_wall(1,j) = m_allDSrigidbodies[id]
                                       ->get_distanceTo( *pt, rayDir, dh );
         geomVector surface_point = *pt + del_wall(1,j)*rayDir;
         geomVector net_vel = rigid_body_velocity(id,surface_point);
         fwall(1,j) = net_vel(comp);
      }
   }

   double field_value = (1./2.)*((del_wall(0,1)*fwall(0,0)
                                + del_wall(0,0)*fwall(0,1))
                                / (del_wall(0,1)+del_wall(0,0))
                               + (del_wall(1,0)*fwall(1,1)
                                + del_wall(1,1)*fwall(1,0))
                                / (del_wall(1,0)+del_wall(1,1)));

   return (field_value);
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

   // force variable for all the rigid bodies
   viscous_force.reserve(m_nrb);
   pressure_force.reserve(m_nrb);
   viscous_torque.reserve(m_nrb);
   pressure_torque.reserve(m_nrb);

   geomVector vvv(3);

   for (size_t i = 0; i < m_nrb; ++i) {
      viscous_force.push_back( vvv );
      pressure_force.push_back( vvv );
      viscous_torque.push_back( vvv );
      pressure_torque.push_back( vvv );
   }

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
