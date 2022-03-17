#include <DS_AllRigidBodies.hh>
#include <FS_AllRigidBodies.hh>
#include <DS_RigidBody.hh>
#include <DS_RigidBody_BuilderFactory.hh>
#include <FV_TimeIterator.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <cmath>
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
                                  , double const& arb_rho
                                  , MAC_DoubleVector const* arb_gv
                                  , double const& arb_scs
                                  , MAC_Communicator const* arb_macCOMM
                                  , double const& arb_mu )
//---------------------------------------------------------------------------
  : m_space_dimension( dimens )
  , UF ( arb_UF )
  , PF ( arb_PF )
  , MESH ( UF->primary_grid() )
  , surface_cell_scale ( arb_scs )
  , m_macCOMM ( arb_macCOMM )
  , m_mu ( arb_mu )
  , m_rho ( arb_rho )
  , gravity_vector ( arb_gv )
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

  initialize_surface_variables_for_all_RB();

  compute_surface_variables_for_all_RB();

  compute_halo_zones_for_all_rigid_body();

  create_neighbour_list_for_AllRB();
}




//---------------------------------------------------------------------------
DS_AllRigidBodies:: DS_AllRigidBodies( size_t& dimens
                                  , istream& in
                                  , bool const& b_particles_as_fixed_obstacles
                                  , FV_DiscreteField const* arb_UF
                                  , FV_DiscreteField const* arb_PF
                                  , FV_DiscreteField const* arb_TF
                                  , double const& arb_scs
                                  , MAC_Communicator const* arb_macCOMM
                                  , double const& arb_mu )
//---------------------------------------------------------------------------
  : m_space_dimension( dimens )
  , UF ( arb_UF )
  , PF ( arb_PF )
  , TF ( arb_TF )
  , MESH ( UF->primary_grid() )
  , surface_cell_scale ( arb_scs )
  , m_macCOMM ( arb_macCOMM )
  , m_mu ( arb_mu )
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

  initialize_surface_variables_for_all_RB();

  compute_surface_variables_for_all_RB();

  compute_halo_zones_for_all_rigid_body();

  create_neighbour_list_for_AllRB();
}




//---------------------------------------------------------------------------
DS_AllRigidBodies:: DS_AllRigidBodies( size_t& dimens
                                  , istream& in
                                  , bool const& b_particles_as_fixed_obstacles
                                  , FV_DiscreteField const* arb_TF
                                  , double const& arb_scs
                                  , MAC_Communicator const* arb_macCOMM
                                  , double const& arb_mu )
//---------------------------------------------------------------------------
  : m_space_dimension( dimens )
  , TF ( arb_TF )
  , MESH ( TF->primary_grid() )
  , surface_cell_scale ( arb_scs )
  , m_macCOMM ( arb_macCOMM )
  , m_mu ( arb_mu )
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

  initialize_surface_variables_for_all_RB();

  compute_surface_variables_for_all_RB();

  compute_halo_zones_for_all_rigid_body();

  create_neighbour_list_for_AllRB();
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
void DS_AllRigidBodies:: compute_pressure_force_and_torque_for_allRB( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies::compute_pressure_force_and_torque_for_allRB") ;

  avg_pressure_force(0) = 0.;
  avg_pressure_force(1) = 0.;
  avg_pressure_force(2) = 0.;

  avg_pressure_torque(0) = 0.;
  avg_pressure_torque(1) = 0.;
  avg_pressure_torque(2) = 0.;

  pressure_force->set(0.);
  pressure_torque->set(0.);

  for (size_t i = 0; i < m_nrb; ++i) {
     first_order_pressure_stress(i);
  }

  m_macCOMM->sum_array(*pressure_force);
  m_macCOMM->sum_array(*pressure_torque);

  for (size_t i = 0; i < m_nrb; ++i) {
     avg_pressure_force(0) += pressure_force->operator()(i,0);
     avg_pressure_force(1) += pressure_force->operator()(i,1);
     avg_pressure_force(2) += pressure_force->operator()(i,2);

     avg_pressure_torque(0) += pressure_torque->operator()(i,0);
     avg_pressure_torque(1) += pressure_torque->operator()(i,1);
     avg_pressure_torque(2) += pressure_torque->operator()(i,2);
  }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_viscous_force_and_torque_for_allRB(
                                                string const& StressOrder )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_viscous_force_and_torque_for_allRB") ;

  avg_viscous_force(0) = 0.;
  avg_viscous_force(1) = 0.;
  avg_viscous_force(2) = 0.;

  avg_viscous_torque(0) = 0.;
  avg_viscous_torque(1) = 0.;
  avg_viscous_torque(2) = 0.;

  viscous_force->set(0.);
  viscous_torque->set(0.);

  for (size_t i = 0; i < m_nrb; ++i) {
     if (StressOrder == "first") {
        first_order_viscous_stress(i);
     } else if (StressOrder == "second") {
        second_order_viscous_stress(i);
     }
  }

  m_macCOMM->sum_array(*viscous_force);
  m_macCOMM->sum_array(*viscous_torque);

  string fileName = "./DS_results/particle_forces.csv" ;
  std::ofstream MyFile( fileName.c_str(), std::ios::app ) ;

  for (size_t i = 0; i < m_nrb; ++i) {
     avg_viscous_force(0) += viscous_force->operator()(i,0);
     avg_viscous_force(1) += viscous_force->operator()(i,1);
     avg_viscous_force(2) += viscous_force->operator()(i,2);

     avg_viscous_torque(0) += viscous_torque->operator()(i,0);
     avg_viscous_torque(1) += viscous_torque->operator()(i,1);
     avg_viscous_torque(2) += viscous_torque->operator()(i,2);

     geomVector const* pgc = m_allDSrigidbodies[i]->get_ptr_to_gravity_centre();
     geomVector pv = rigid_body_velocity(i,*pgc);

     if (m_macCOMM->rank() == 0) {
        MyFile << i << " , " << pgc->operator()(0)
                    << " , " << pgc->operator()(1)
                    << " , " << pgc->operator()(2)
                    << " , " << pv(0)
                    << " , " << pv(1)
                    << " , " << pv(2)
                    << " , " << pressure_force->operator()(i,0)
                    << " , " << pressure_force->operator()(i,1)
                    << " , " << pressure_force->operator()(i,2)
                    << " , " << viscous_force->operator()(i,0)
                    << " , " << viscous_force->operator()(i,1)
                    << " , " << viscous_force->operator()(i,2)
                    << " , " << pressure_torque->operator()(i,0)
                    << " , " << pressure_torque->operator()(i,1)
                    << " , " << pressure_torque->operator()(i,2)
                    << " , " << viscous_torque->operator()(i,0)
                    << " , " << viscous_torque->operator()(i,1)
                    << " , " << viscous_torque->operator()(i,2)
                    << endl;
     }
  }

  if (m_macCOMM->rank() == 0) MyFile.close( ) ;

  if (m_macCOMM->rank() == 0) {
     std::cout << "Average pressure force on RB: "
               << avg_pressure_force(0) / double(m_nrb) << " ,"
               << avg_pressure_force(1) / double(m_nrb) << " ,"
               << avg_pressure_force(2) / double(m_nrb) << endl;

     std::cout << "Average viscous force on RB: "
               << avg_viscous_force(0) / double(m_nrb) << " ,"
               << avg_viscous_force(1) / double(m_nrb) << " ,"
               << avg_viscous_force(2) / double(m_nrb) << endl;

     std::cout << "Average pressure torque on RB: "
               << avg_pressure_torque(0) / double(m_nrb) << " ,"
               << avg_pressure_torque(1) / double(m_nrb) << " ,"
               << avg_pressure_torque(2) / double(m_nrb) << endl;

     std::cout << "Average viscous torque on RB: "
               << avg_viscous_torque(0) / double(m_nrb) << " ,"
               << avg_viscous_torque(1) / double(m_nrb) << " ,"
               << avg_viscous_torque(2) / double(m_nrb) << endl;
  }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: solve_RB_equation_of_motion(
                     FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: solve_RB_equation_of_motion()" ) ;

  // Get gravity vector from NS solver
  doubleVector const& gg = gravity_vector->to_double_vector();


  for (size_t parID = 0; parID < m_nrb; ++parID) {
     // Get RB gravity center
     geomVector const* pgc = m_allDSrigidbodies[parID]
                                             ->get_ptr_to_gravity_centre();
     geomVector pv = rigid_body_velocity(parID,*pgc);

     // Get solid mass and density from FS class
     auto value = m_allDSrigidbodies[parID]->get_mass_and_density();
     double mass_p = std::get<0>(value);
     double rho_p = std::get<1>(value);

     geomVector pos(3,0);
     geomVector vel(3,0);
     geomVector acc(3,0);
     geomVector delta(3);

     // Solving equation of motion
     for (size_t dir=0;dir<m_space_dimension;dir++) {
        pos(dir) = pgc->operator()(dir);
        vel(dir) = pv(dir);

        acc(dir) = gg(dir)*(1-m_rho/rho_p) + (viscous_force->operator()(parID,dir)
                                         +  pressure_force->operator()(parID,dir))
                                         / mass_p ;
        vel(dir) = vel(dir) + acc(dir)*t_it->time_step();
        pos(dir) = periodic_transformation(pos(dir) + vel(dir)*t_it->time_step()
                                          , dir) ;
        delta(dir) = vel(dir)*t_it->time_step() ;
     }

     // Temporary added this part to create the clones if RB near boundary
     // Don't trust this for all the cases. Primary purpose for settling
     // sphere case
     boolVector const* periodic_comp = MESH->get_periodic_directions();
     double radius = m_allDSrigidbodies[parID]->get_circumscribed_radius();

     vector<geomVector> m_periodic_directions; /**< periodic clones whose
     position is gravity_center + (*periodic_directions)[i] */
     m_periodic_directions.reserve(3);

     for (size_t dir = 0; dir < m_space_dimension; dir++) {
        double Dmin = MESH->get_main_domain_min_coordinate( dir );
        double Dmax = MESH->get_main_domain_max_coordinate( dir );
        if (periodic_comp->operator()(dir)) {
           if (std::abs(pos(dir)-Dmin) <= radius) {
              geomVector delta_clone(3);
              delta_clone(dir) = (Dmax - Dmin);
              m_periodic_directions.push_back(delta_clone);
           } else if (std::abs(pos(dir)-Dmax) <= radius) {
              geomVector delta_clone(3);
              delta_clone(dir) = - (Dmax - Dmin);
              m_periodic_directions.push_back(delta_clone);
           }
        }
     }

     if (m_periodic_directions.size() == 2) {
        geomVector delta_clone(3);
        delta_clone(0) = m_periodic_directions[0](0) + m_periodic_directions[1](0);
        delta_clone(1) = m_periodic_directions[0](1) + m_periodic_directions[1](1);
        delta_clone(2) = m_periodic_directions[0](2) + m_periodic_directions[1](2);
        m_periodic_directions.push_back(delta_clone);
     }
     // // Writing istream to update the positions
     // write_new_parameters_to_stream();
     //
     // Updating the RB data
     m_allDSrigidbodies[parID]->update_RB_position_and_velocity(pos,vel
                                                ,m_periodic_directions);
     m_allDSrigidbodies[parID]->translate_surface_points(delta);
     m_allDSrigidbodies[parID]->correct_surface_discretization( MESH );

     if ( !m_periodic_directions.empty() ) m_periodic_directions.clear();
  }
}




//---------------------------------------------------------------------------
double DS_AllRigidBodies:: get_min_RB_coord( size_t const& dir )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: get_min_RB_coord()" ) ;

  double value = 1000000.;

  for (size_t parID = 0; parID < m_nrb; ++parID) {
     geomVector const* pgc = m_allDSrigidbodies[parID]
                                          ->get_ptr_to_gravity_centre();
     value = std::min(value,pgc->operator()(dir));
  }

  return(value);

}



//---------------------------------------------------------------------------
int DS_AllRigidBodies:: isIn_any_RB( size_t const& ownID
                                   , geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn_any_RB(ownID,pt)" ) ;

  for (size_t i = 0; i < neighbour_list[ownID].size(); ++i) {
     size_t neighID = neighbour_list[ownID][i];
     if (m_allDSrigidbodies[neighID]->isIn( pt )) return ((int)neighID);
  }

  return (-1);

}




//---------------------------------------------------------------------------
int DS_AllRigidBodies:: isIn_any_RB( size_t const& ownID,
                                     double const& x,
                                     double const& y,
                                     double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: isIn_any_RB(ownID,x,y,z)" ) ;

  for (size_t i = 0; i < neighbour_list[ownID].size(); ++i) {
     size_t neighID = neighbour_list[ownID][i];
     if (m_allDSrigidbodies[neighID]->isIn( x, y, z )) return ((int)neighID);
  }

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
  MAC_LABEL( "DS_AllRigidBodies:: rigid_body_velocity(pt)" ) ;

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
  size_t field = field_num(FF) ;

  void_fraction[field]->set(0);

  boolVector const* periodic_comp = MESH->get_periodic_directions();

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
           double domain_min = MESH->get_main_domain_min_coordinate( dir );
           double domain_max = MESH->get_main_domain_max_coordinate( dir );

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

           if (periodic_comp->operator()( 0 ))
              xC = periodic_transformation(xC,0);

           for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
             double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;

             if (periodic_comp->operator()( 1 ))
                yC = periodic_transformation(yC,1);

             for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
               double zC = (m_space_dimension == 2) ? 0.
                                  : FF->get_DOF_coordinate( k, comp, 2 ) ;

               if (periodic_comp->operator()( 2 ))
                  zC = periodic_transformation(zC,2);

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
double DS_AllRigidBodies:: periodic_transformation( double const& x
                                                , size_t const& dir)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: periodic_transformation" ) ;

  double isize = MESH->get_main_domain_max_coordinate(dir)
               - MESH->get_main_domain_min_coordinate(dir);
  double imin = MESH->get_main_domain_min_coordinate(dir);
  double value = x - MAC::floor((x-imin)/isize)*isize;

  return (value);
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
  size_t field = field_num(FF) ;

  intersect_vector[field]->set(0);
  intersect_distance[field]->set(0.);
  intersect_fieldValue[field]->set(0.);

  boolVector const* periodic_comp = MESH->get_periodic_directions();

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

     double delta = MESH->get_smallest_grid_size();

     for (size_t comp = 0; comp < nb_comps; comp++) {

        for (size_t dir = 0; dir < m_space_dimension; dir++) {
          // To include knowns at dirichlet boundary in the intersection
          // calculation as well, modification of the looping extents are required
          min_unknown_index(dir) =
                              FF->get_min_index_unknown_on_proc( comp, dir );
          max_unknown_index(dir) =
                              FF->get_max_index_unknown_on_proc( comp, dir );
          local_extents(dir,0) = 0;
          local_extents(dir,1) = max_unknown_index(dir)
                               - min_unknown_index(dir);

          bool is_periodic = periodic_comp->operator()( dir );

          double domain_min = MESH->get_main_domain_min_coordinate(dir);
          double domain_max = MESH->get_main_domain_max_coordinate(dir);
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

          if (periodic_comp->operator()( 0 ))
             xC = periodic_transformation(xC,0);

          for (size_t j = min_nearest_index(1); j < max_nearest_index(1); ++j) {
              ipos(1) = j - min_unknown_index(1);
              double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;

              if (periodic_comp->operator()( 1 ))
                 yC = periodic_transformation(yC,1);

              for (size_t k = min_nearest_index(2);k < max_nearest_index(2); ++k) {
                 ipos(2) = k - min_unknown_index(2);
                 double zC = (m_space_dimension == 2) ? 0
                                        : FF->get_DOF_coordinate( k, comp, 2 );

                 if (periodic_comp->operator()( 2 ))
                    zC = periodic_transformation(zC,2);

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
                               geomVector rayVec(3);
                               rayVec(0) = source(0) + t * rayDir(0);
                               rayVec(1) = source(1) + t * rayDir(1);
                               rayVec(2) = source(2) + t * rayDir(2);
                               geomVector netVel =
                                rigid_body_velocity(parID, rayVec);
                                // rigid_body_velocity(parID, source + t * rayDir);
                               // Value of variable at the surface of particle
                               if ( FF == PF ) { // i.e. PF
                                  intersect_fieldValue[field]->operator()(p,col)
                                                                 = netVel(dir);
                               } else if (FF == UF) { // i.e. UF
                                  intersect_fieldValue[field]->operator()(p,col)
                                                                 = netVel(comp);
                               } else if (FF == TF) {
                                  intersect_fieldValue[field]->operator()(p,col)
                                                                 = 0.;
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

   double dx = MESH->get_smallest_grid_size();

   for (size_t i = 0; i < m_nrb; ++i) {
      m_allDSrigidbodies[i]->compute_number_of_surface_variables(
                                          surface_cell_scale, dx);
      m_allDSrigidbodies[i]->initialize_surface_variables( );
   }

   // Reserving three location in memory for PF, UF, TF grid, respectively.
   void_fraction.reserve(3);
   void_fraction.push_back(new size_t_vector(1,0));
   void_fraction.push_back(new size_t_vector(1,0));
   void_fraction.push_back(new size_t_vector(1,0));
   intersect_vector.reserve(3);
   intersect_vector.push_back(new size_t_array2D(1,1,0));
   intersect_vector.push_back(new size_t_array2D(1,1,0));
   intersect_vector.push_back(new size_t_array2D(1,1,0));
   intersect_distance.reserve(3);
   intersect_distance.push_back(new doubleArray2D(1,1,0.));
   intersect_distance.push_back(new doubleArray2D(1,1,0.));
   intersect_distance.push_back(new doubleArray2D(1,1,0.));
   intersect_fieldValue.reserve(3);
   intersect_fieldValue.push_back(new doubleArray2D(1,1,0.));
   intersect_fieldValue.push_back(new doubleArray2D(1,1,0.));
   intersect_fieldValue.push_back(new doubleArray2D(1,1,0.));

   // Intialization of force and torque variables
   viscous_force = new doubleArray2D(1,1,0.);
   pressure_force = new doubleArray2D(1,1,0.);
   viscous_force->re_initialize(m_nrb,3);
   pressure_force->re_initialize(m_nrb,3);

   viscous_torque = new doubleArray2D(1,1,0.);
   pressure_torque = new doubleArray2D(1,1,0.);
   viscous_torque->re_initialize(m_nrb,3);
   pressure_torque->re_initialize(m_nrb,3);

   geomVector vvv(3);

   avg_pressure_force = vvv;
   avg_viscous_force = vvv;
   avg_pressure_torque = vvv;
   avg_viscous_torque = vvv;

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: first_order_pressure_stress( size_t const& parID )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllRigidBodies:: first_order_pressure_stress" ) ;

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
  geomVector const* pgc = m_allDSrigidbodies[parID]
                                          ->get_ptr_to_gravity_centre();


  // Get local min and max indices
  // One extra grid cell needs to considered, since ghost points can be
  // located in between the min/max index handled by the proc
  for (size_t l = 0; l < m_space_dimension; l++) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( comp, l );
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( comp, l );
     Dmin(l) = MESH->get_min_coordinate_on_current_processor( l );
     Dmax(l) = MESH->get_max_coordinate_on_current_processor( l );
     domain_length(l) = MESH->get_main_domain_max_coordinate( l )
                      - MESH->get_main_domain_min_coordinate( l );
     domain_min(l) = MESH->get_main_domain_min_coordinate( l );
  }

  size_t pfd = 0;
  double external_gradP = 0.;
  double isize = 0.;

  if ( MESH->is_periodic_flow() ) {
     pfd = MESH->get_periodic_flow_direction() ;

     external_gradP = MESH->get_periodic_pressure_drop() /
              ( MESH->get_main_domain_max_coordinate( pfd )
              - MESH->get_main_domain_min_coordinate( pfd ) ) ;

     isize = MESH->get_main_domain_max_coordinate(pfd)
           - MESH->get_main_domain_min_coordinate(pfd);
  }

  for (size_t i = 0; i < surface_area.size(); i++) {
     double stress = 0.;

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
         Trilinear_interpolation(PF,comp,surface_point[i],i0,parID,{0,1}) ;
        stress = - press/2.;

        if (MESH->is_periodic_flow() &&
            (external_gradP != 0.)) {
           double delta = surface_point[i]->operator()(pfd)
                        - pgc->operator()(pfd);
           delta = delta - round(delta/isize) * isize;
           stress += - (external_gradP * (delta + pgc->operator()(pfd)));
        }
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

     // value(0) = m_macCOMM->sum(value(0));
     // value(1) = m_macCOMM->sum(value(1));
     // value(2) = m_macCOMM->sum(value(2));

     m_allDSrigidbodies[parID]->update_Pforce_on_surface_point(i,value);

     pressure_force->operator()(parID,0) += value(0);
     pressure_force->operator()(parID,1) += value(1);
     pressure_force->operator()(parID,2) += value(2);

     pressure_torque->operator()(parID,0)
                               += value(2)*(surface_point[i]->operator()(1)
                                                       - pgc->operator()(1))
                                - value(1)*(surface_point[i]->operator()(2)
                                                       - pgc->operator()(2));
     pressure_torque->operator()(parID,1)
                               += value(0)*(surface_point[i]->operator()(2)
                                                       - pgc->operator()(2))
                                - value(2)*(surface_point[i]->operator()(0)
                                                       - pgc->operator()(0));
     pressure_torque->operator()(parID,2)
                               += value(1)*(surface_point[i]->operator()(0)
                                                       - pgc->operator()(0))
                                - value(0)*(surface_point[i]->operator()(1)
                                                       - pgc->operator()(1));

  }

}




//---------------------------------------------------------------------------
void
DS_AllRigidBodies:: second_order_viscous_stress(size_t const& parID)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllRigidBodies:: second_order_viscous_stress" ) ;

  size_t nb_comps = UF->nb_components() ;

  vector<geomVector*> surface_point = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_points();
  vector<geomVector*> surface_normal = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_normals();
  vector<geomVector*> surface_area = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_areas();
  geomVector const* pgc = m_allDSrigidbodies[parID]
                                          ->get_ptr_to_gravity_centre();

  // 6 ghost points and 1 surface point
  size_t_vector vvv(3,0);
  vector<geomVector> ghost_pt(7,0);
  vector<double> f(7,0.);
  vector<size_t_vector> i0_new(7,vvv);
  vector<int> in_parID(7,0);
  vector<int> sign(3,0);
  boolVector in_domain(7,true);
  // Required for the call to Bisection in case of 2D domain
  size_t_vector face_vector(3,0);
  face_vector(0) = 1; face_vector(1) = 1; face_vector(2) = 0;
  vector<double> net_vel(3,0.);

  size_t_vector min_unknown_index(m_space_dimension,0);
  size_t_vector max_unknown_index(m_space_dimension,0);
  // Domain length and minimum
  geomVector domain_length(3), domain_min(3);
  // Extents on the currect processor
  geomVector Dmin(3), Dmax(3);

  // Get local min and max indices
  // One extra grid cell needs to considered, since ghost points can be
  // located in between the min/max index handled by the proc
  for (size_t l = 0; l < m_space_dimension; l++) {
     Dmin(l) = MESH->get_min_coordinate_on_current_processor( l );
     Dmax(l) = MESH->get_max_coordinate_on_current_processor( l );
     domain_length(l) = MESH->get_main_domain_max_coordinate( l )
                      - MESH->get_main_domain_min_coordinate( l );
     domain_min(l) = MESH->get_main_domain_min_coordinate( l );
  }

  for (size_t i = 0; i < surface_area.size(); i++) {
     ghost_pt[0] = *surface_point[i];
     // Check it the point is in the current domain
     in_domain(0) = (ghost_pt[0](0) > Dmin(0)) && (ghost_pt[0](0) <= Dmax(0))
                 && (ghost_pt[0](1) > Dmin(1)) && (ghost_pt[0](1) <= Dmax(1))
                 && (ghost_pt[0](2) > Dmin(2)) && (ghost_pt[0](2) <= Dmax(2));

     doubleVector stress(6,0.);

     if (in_domain(0)) {

        for (size_t l = 0; l < m_space_dimension; l++)
           sign[l] = (surface_normal[i]->operator()(l) > 0.) ? 1 : -1 ;


        for (size_t comp = 0; comp < nb_comps; comp++) {
           // Finding the grid indexes next to ghost points
           for (size_t l = 0; l < m_space_dimension; l++) {
              size_t i0_temp;
              bool found = FV_Mesh::between(
                                 UF->get_DOF_coordinates_vector(comp,l),
                                 ghost_pt[0](l),
                                 i0_temp);
              if (found) i0_new[0](l) = i0_temp;
           }

           // Ghost points generation in each direction
           for (size_t dir = 0; dir < m_space_dimension; dir++) {
              // 1,2,3,4,5,6 are x1, x2, y1, y2, z1, z2 points, respectively.
              size_t col1 = 2*dir + 0 + 1;
              size_t col2 = 2*dir + 1 + 1;
              ghost_pt[col1] = *surface_point[i];
              ghost_pt[col2] = *surface_point[i];
              i0_new[col1] = i0_new[0];
              i0_new[col2] = i0_new[0];

              intVector i0_temp(2,0);

              // Ghost points in i for the calculation of i-derivative of field
              i0_temp(0) = (sign[dir] == 1) ? (int(i0_new[0](dir)) + 2*sign[dir])
                                            : (int(i0_new[0](dir)) + 1*sign[dir]);
              i0_temp(1) = (sign[dir] == 1) ? (int(i0_new[0](dir)) + 3*sign[dir])
                                            : (int(i0_new[0](dir)) + 2*sign[dir]);

              if ((i0_temp(0) >= 0) &&
                  (i0_temp(0) < (int)UF->get_local_nb_dof(comp,dir))) {
                 i0_new[col1](dir) = i0_temp(0);
                 ghost_pt[col1](dir) = UF->get_DOF_coordinate(i0_temp(0), comp, dir);
                 in_domain(col1) = 1;
              } else {
                 in_domain(col1) = 0;
              }

              if ((i0_temp(1) >= 0) &&
                  (i0_temp(1) < (int)UF->get_local_nb_dof(comp,dir))) {
                 i0_new[col2](dir) = i0_temp(1);
                 ghost_pt[col2](dir) = UF->get_DOF_coordinate(i0_temp(1), comp, dir);
                 in_domain(col2) = 1;
              } else {
                 in_domain(col2) = 0;
              }

              // Checking all the ghost points in the solid/fluid,
              // and storing the parID if present in solid
              in_parID[col1] = isIn_any_RB(parID,ghost_pt[col1]);
              in_parID[col2] = isIn_any_RB(parID,ghost_pt[col2]);
           }

           // Get local min and max indices
           for (size_t l = 0; l < m_space_dimension; l++) {
              min_unknown_index(l) =
                           UF->get_min_index_unknown_handled_by_proc( comp, l );
              max_unknown_index(l) =
                           UF->get_max_index_unknown_handled_by_proc( comp, l );
           }

           // Calculation of field variable on the surface point i
           geomVector netVel = rigid_body_velocity(parID, ghost_pt[0]);
           f[0] = netVel(comp);

           // Calculation of field variable on the ghost points
           for (size_t dir = 0; dir < m_space_dimension; dir++) {
              for (size_t ig = 0; ig < 2; ig++) {
                 // 1,2,3,4,5,6 are x1, x2, y1, y2, z1, z2 points, respectively.
                 size_t col = 2*dir + ig + 1;

                 if ((in_parID[col] == -1) && in_domain(col)) {
                    f[col] = (m_space_dimension == 2) ?
                                       Biquadratic_interpolation(UF
                                                           , comp
                                                           , &ghost_pt[col]
                                                           , i0_new[col]
                                                           , (dir == 0) ? 1 : 0
                                                           , sign[dir]
                                                           , {0})
                                     : Triquadratic_interpolation(UF
                                                            , comp
                                                            , &ghost_pt[col]
                                                            , i0_new[col]
                                                            , parID
                                                            , dir
                                                            , sign
                                                            , {0}) ;
                 } else if ((in_parID[col] != -1) && in_domain(col)) {
                    geomVector netVelg = rigid_body_velocity(in_parID[col]
                                                           , ghost_pt[col]);
                    f[col] = netVelg(comp);
                 }
              }
           }

           geomVector dfdi(3);

           // Calculation of derivatives in each direction
           for (size_t dir = 0; dir < m_space_dimension; dir++) {
              size_t col1 = 2*dir + 1;
              size_t col2 = 2*dir + 2;

              // Point 1 and 2 in computational domain
   	        if (in_domain(col1) && in_domain(col2)) {
                 if ((in_parID[col1] == -1) && (in_parID[col2] == -1)) {
                    double dx1 = ghost_pt[col1](dir) - ghost_pt[0](dir);
                    double dx2 = ghost_pt[col2](dir) - ghost_pt[0](dir);
                    dfdi(dir) = ((f[col1] - f[0])*dx2/dx1
                               - (f[col2] - f[0])*dx1/dx2)
                               / (dx2-dx1);
                 // Point 1 in fluid and 2 in the solid
                 } else if ((in_parID[col1] == -1) && (in_parID[col2] != -1)) {
                    double dx1 = ghost_pt[col1](dir) - ghost_pt[0](dir);
                    dfdi(dir) = (f[col1] - f[0])/dx1;
                 // Point 1 is present in solid
                 } else if (in_parID[col1] != -1) {
   	              double dx1 = ghost_pt[col1](dir) - ghost_pt[0](dir);
                    dfdi(dir) = (f[col1] - f[0])/dx1;
   	           }
              // Point 1 in computational domain
              } else if (in_domain(col1) && !in_domain(col2)) {
                 double dx1 = ghost_pt[col1](dir) - ghost_pt[0](dir);
                 dfdi(dir) = (f[col1] - f[0])/dx1;
              // Particle close to wall
              } else if (!in_domain(col1)) {
                 i0_new[col1](dir) = (sign[dir] == 1) ? (i0_new[0](dir) + 1*sign[dir])
                                                      : (i0_new[0](dir) + 0*sign[dir]);
                 ghost_pt[col1](dir) = UF->get_DOF_coordinate(i0_new[col1](dir), comp, dir);
                 f[col1] = (m_space_dimension == 2) ?
                                         Biquadratic_interpolation(UF
                                                             , comp
                                                             , &ghost_pt[col1]
                                                             , i0_new[col1]
                                                             , (dir == 0) ? 1 : 0
                                                             , sign[dir]
                                                             , {0})
                                       : Triquadratic_interpolation(UF
                                                              , comp
                                                              , &ghost_pt[col1]
                                                              , i0_new[col1]
                                                              , parID
                                                              , dir
                                                              , sign
                                                              , {0}) ;
                 double dx1 = ghost_pt[col1](dir) - ghost_pt[0](dir);
                 dfdi(dir) = (f[col1] - f[0])/dx1;
              }

              dfdi(dir) *= m_mu;//*sign[dir];
           }

           if (comp == 0) {
              stress(0) = 2.*dfdi(0);
              stress(3) = stress(3) + dfdi(1);
              stress(5) = stress(5) + dfdi(2);
           } else if (comp == 1) {
              stress(1) = 2.*dfdi(1);
              stress(3) = stress(3) + dfdi(0);
              stress(4) = stress(4) + dfdi(2);
           } else if (comp == 2) {
              stress(2) = 2.*dfdi(2);
              stress(4) = stress(4) + dfdi(1);
              stress(5) = stress(5) + dfdi(0);
           }
        }
     }

     // Ref: Keating thesis Pg-85
     // point_coord*(area) --> Component of area in particular direction
     geomVector value(3);
     double norm = surface_normal[i]->calcNorm();
     value(0) = (stress(0)*surface_normal[i]->operator()(0)
               + stress(3)*surface_normal[i]->operator()(1)
               + stress(5)*surface_normal[i]->operator()(2))/norm
               * surface_area[i]->operator()(0);

     value(1) = (stress(3)*surface_normal[i]->operator()(0)
               + stress(1)*surface_normal[i]->operator()(1)
               + stress(4)*surface_normal[i]->operator()(2))/norm
               * surface_area[i]->operator()(0);

     value(2) = (stress(5)*surface_normal[i]->operator()(0)
               + stress(4)*surface_normal[i]->operator()(1)
               + stress(2)*surface_normal[i]->operator()(2))/norm
               * surface_area[i]->operator()(0);

     m_allDSrigidbodies[parID]->update_Vforce_on_surface_point(i,value);

     viscous_force->operator()(parID,0) += value(0);
     viscous_force->operator()(parID,1) += value(1);
     viscous_force->operator()(parID,2) += value(2);

     viscous_torque->operator()(parID,0)
                              += value(2)*(surface_point[i]->operator()(1)
                                                      - pgc->operator()(1))
                               - value(1)*(surface_point[i]->operator()(2)
                                                      - pgc->operator()(2));
     viscous_torque->operator()(parID,1)
                              += value(0)*(surface_point[i]->operator()(2)
                                                      - pgc->operator()(2))
                               - value(2)*(surface_point[i]->operator()(0)
                                                      - pgc->operator()(0));
     viscous_torque->operator()(parID,2)
                              += value(1)*(surface_point[i]->operator()(0)
                                                      - pgc->operator()(0))
                               - value(0)*(surface_point[i]->operator()(1)
                                                      - pgc->operator()(1));

  }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: first_order_viscous_stress( size_t const& parID )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllRigidBodies:: first_order_viscous_stress" ) ;

  size_t nb_comps = UF->nb_components() ;
  double dh = MESH->get_smallest_grid_size();

  vector<geomVector*> surface_point = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_points();
  vector<geomVector*> surface_normal = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_normals();
  vector<geomVector*> surface_area = m_allDSrigidbodies[parID]
                                          ->get_rigid_body_surface_areas();
  geomVector const* pgc = m_allDSrigidbodies[parID]
                                          ->get_ptr_to_gravity_centre();

  // 6 ghost points and 1 surface point
  vector<geomVector> ghost_pt(7,0);
  vector<double> f(7,0.);
  boolArray2D found(7,3,0);
  size_t_vector vvv(3,0);
  vector<size_t_vector> i0_new(7,vvv);
  vector<int> in_parID(7,0);
  boolVector in_domain(7,true);
  // Required for the call to Bisection in case of 2D domain
  size_t_vector face_vector(3,0);
  face_vector(0) = 1; face_vector(1) = 1; face_vector(2) = 0;
  vector<double> net_vel(3,0.);

  size_t_vector min_unknown_index(m_space_dimension,0);
  size_t_vector max_unknown_index(m_space_dimension,0);
  // Domain length and minimum
  geomVector domain_length(3), domain_min(3);
  // Extents on the currect processor
  geomVector Dmin(3), Dmax(3);

  // Get local min and max indices
  // One extra grid cell needs to considered, since ghost points can be
  // located in between the min/max index handled by the proc
  for (size_t l = 0; l < m_space_dimension; l++) {
     Dmin(l) = MESH->get_min_coordinate_on_current_processor( l );
     Dmax(l) = MESH->get_max_coordinate_on_current_processor( l );
     domain_length(l) = MESH->get_main_domain_max_coordinate( l )
                      - MESH->get_main_domain_min_coordinate( l );
     domain_min(l) = MESH->get_main_domain_min_coordinate( l );
  }


  for (size_t i = 0; i < surface_area.size(); i++) {
     // Check it the point is in the current domain
     ghost_pt[0] = *surface_point[i];
     bool status = (ghost_pt[0](0) > Dmin(0))
                && (ghost_pt[0](0) <= Dmax(0))
                && (ghost_pt[0](1) > Dmin(1))
                && (ghost_pt[0](1) <= Dmax(1))
                && (ghost_pt[0](2) > Dmin(2))
                && (ghost_pt[0](2) <= Dmax(2));

     doubleVector stress(6,0.);

     if (status) {

        // Ghost points generation in each direction
        for (size_t dir = 0; dir < m_space_dimension; dir++) {
           for (size_t ig = 0; ig < 2; ig++) {
              // 1,2,3,4,5,6 are x1, x2, y1, y2, z1, z2 points, respectively.
              size_t col = 2*dir + ig + 1;
              ghost_pt[col] = *surface_point[i];
              double sign = (surface_normal[i]->operator()(dir) > 0.) ? 1.:-1.;
              ghost_pt[col](dir) += double(ig + 1) * sign * dh;
              ghost_pt[col](dir) += - MAC::floor((ghost_pt[col](dir)
                                                 -domain_min(dir))
                                                /domain_length(dir))
                                    *domain_length(dir);

              // Checking all the ghost points in the solid/fluid,
              // and storing the parID if present in solid
              in_parID[col] = isIn_any_RB(parID,ghost_pt[col]);
           }
        }

        for (size_t comp = 0; comp < nb_comps; comp++) {
           // Get local min and max indices
           for (size_t l = 0; l < m_space_dimension; l++) {
              min_unknown_index(l) =
                           UF->get_min_index_unknown_handled_by_proc( comp, l );
              max_unknown_index(l) =
                           UF->get_max_index_unknown_handled_by_proc( comp, l );
           }

           // Checking the grid indexes for each ghost point
           for (size_t col = 0; col < 7; col++) {
              for (size_t dir = 0; dir < m_space_dimension; dir++) {
                 size_t i0_temp;
                 found(col,dir) = FV_Mesh::between(
                                       UF->get_DOF_coordinates_vector(comp,dir)
                                     , ghost_pt[col](dir)
                                     , i0_temp);
                 if (found(col,dir) == 1) i0_new[col](dir) = i0_temp;

              }
           }

           // Calculation of field variable on the surface point i
           geomVector netVel = rigid_body_velocity(parID, *surface_point[i]);
           f[0] = netVel(comp);

           // Calculation of field variable on the ghost points
           for (size_t dir = 0; dir < m_space_dimension; dir++) {
              for (size_t ig = 0; ig < 2; ig++) {
                 // 1,2,3,4,5,6 are x1, x2, y1, y2, z1, z2 points, respectively.
                 size_t col = 2*dir + ig + 1;
                 in_domain(col) = (m_space_dimension == 2) ? found(col,0)
                                                          && found(col,1)
                                                           : found(col,0)
                                                          && found(col,1)
                                                          && found(col,2);

                 if ((in_parID[col] == -1) && in_domain(col)) {
                    f[col] = (m_space_dimension == 2) ?
                                       Bilinear_interpolation(UF
                                                           , comp
                                                           , &ghost_pt[col]
                                                           , i0_new[col]
                                                           , face_vector
                                                           , {0})
                                     : Trilinear_interpolation(UF
                                                            , comp
                                                            , &ghost_pt[col]
                                                            , i0_new[col]
                                                            , parID
                                                            , {0}) ;
                 } else if ((in_parID[col] != -1) && in_domain(col)) {
                    geomVector netVelg = rigid_body_velocity(in_parID[col]
                                                           , ghost_pt[col]);
                    f[col] = netVelg(comp);
                 }
              }
           }

           geomVector dfdi(3);

           // Calculation of derivatives in each direction
           for (size_t dir = 0; dir < m_space_dimension; dir++) {
              double sign = (surface_normal[i]->operator()(dir) > 0.) ? 1.:-1.;
              size_t col1 = 2*dir + 1;
              size_t col2 = 2*dir + 2;

              if ((in_parID[col1] == -1) && (in_parID[col2] == -1)
               && (in_domain(col1) && in_domain(col2))) {
                 dfdi(dir) = (-f[col2] + 4.*f[col1] - 3.*f[0])/2./dh;
              // Point 1 in fluid and 2 is either in the solid
              // or out of the computational domain
              } else if ((in_parID[col1] == -1) && ((in_parID[col2] != -1)
                     || ((in_domain(col2) == 0) && (in_domain(col1) == 1)))) {
                 dfdi(dir) = (f[col1] - f[0])/dh;
              // Point 1 is present in solid
              } else if (in_parID[col1] != -1) {
                 geomVector netVelg = rigid_body_velocity(in_parID[col1]
                                                        , ghost_pt[col1]);
                 dfdi(dir) = (netVelg(comp) - f[0])/dh;
              // Point 1 is out of the computational domain
              } else if (in_domain(col1) == 0) {
                 double dh_wall = (sign > 0.) ?
                     MAC::abs(surface_point[i]->operator()(dir)
                   - MESH->get_main_domain_max_coordinate(dir)) :
                     MAC::abs(surface_point[i]->operator()(dir)
                   - MESH->get_main_domain_min_coordinate(dir)) ;

                 size_t_vector i0(3,0);
                 // Point on the domain boundary
                 geomVector pt(3);
                 pt = ghost_pt[0];
                 pt(dir) += sign*dh_wall;

                 for (size_t l = 0; l < m_space_dimension; l++) {
                    size_t i0_temp;
                    bool found_temp =
                     FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,l)
                                    , pt(l)
                                    , i0_temp);
                    if (found_temp == 1) i0(l) = i0_temp;
                 }
                 f[col1] = (m_space_dimension == 2) ?
                                          Bilinear_interpolation(UF
                                                               , comp
                                                               , &pt
                                                               , i0
                                                               , face_vector
                                                               , {0})
                                        : Trilinear_interpolation(UF
                                                               , comp
                                                               , &pt
                                                               , i0
                                                               , parID
                                                               , {0}) ;
                 dfdi(dir) = (f[col1] - f[0])/dh_wall;
              }
              dfdi(dir) *= m_mu*sign;
           }

           if (comp == 0) {
              stress(0) = 2.*dfdi(0);
              stress(3) = stress(3) + dfdi(1);
              stress(5) = stress(5) + dfdi(2);
           } else if (comp == 1) {
              stress(1) = 2.*dfdi(1);
              stress(3) = stress(3) + dfdi(0);
              stress(4) = stress(4) + dfdi(2);
           } else if (comp == 2) {
              stress(2) = 2.*dfdi(2);
              stress(4) = stress(4) + dfdi(1);
              stress(5) = stress(5) + dfdi(0);
           }
        }
     }

     // Ref: Keating thesis Pg-85
     // point_coord*(area) --> Component of area in particular direction
     geomVector value(3);
     double norm = surface_normal[i]->calcNorm();

     value(0) = (stress(0)*surface_normal[i]->operator()(0)
               + stress(3)*surface_normal[i]->operator()(1)
               + stress(5)*surface_normal[i]->operator()(2))/norm
               * surface_area[i]->operator()(0);

     value(1) = (stress(3)*surface_normal[i]->operator()(0)
               + stress(1)*surface_normal[i]->operator()(1)
               + stress(4)*surface_normal[i]->operator()(2))/norm
               * surface_area[i]->operator()(0);

     value(2) = (stress(5)*surface_normal[i]->operator()(0)
               + stress(4)*surface_normal[i]->operator()(1)
               + stress(2)*surface_normal[i]->operator()(2))/norm
               * surface_area[i]->operator()(0);

     // value(0) = m_macCOMM->sum(value(0));
     // value(1) = m_macCOMM->sum(value(1));
     // value(2) = m_macCOMM->sum(value(2));

     m_allDSrigidbodies[parID]->update_Vforce_on_surface_point(i,value);

     viscous_force->operator()(parID,0) += value(0);
     viscous_force->operator()(parID,1) += value(1);
     viscous_force->operator()(parID,2) += value(2);

     viscous_torque->operator()(parID,0)
                              += value(2)*(surface_point[i]->operator()(1)
                                                      - pgc->operator()(1))
                               - value(1)*(surface_point[i]->operator()(2)
                                                      - pgc->operator()(2));
     viscous_torque->operator()(parID,1)
                              += value(0)*(surface_point[i]->operator()(2)
                                                      - pgc->operator()(2))
                               - value(2)*(surface_point[i]->operator()(0)
                                                      - pgc->operator()(0));
     viscous_torque->operator()(parID,2)
                              += value(1)*(surface_point[i]->operator()(0)
                                                      - pgc->operator()(0))
                               - value(0)*(surface_point[i]->operator()(1)
                                                      - pgc->operator()(1));

  }

}




//---------------------------------------------------------------------------
double DS_AllRigidBodies:: Trilinear_interpolation ( FV_DiscreteField const* FF
                                                   , size_t const& comp
                                                   , geomVector const* pt
                                                   , size_t_vector const& i0
                                                   , size_t const& parID
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
   double dh = MESH->get_smallest_grid_size();
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

      in_RB = (FF == UF) ? isIn_any_RB(parID, pt_at_face) : -1 ;

      if (in_RB != -1) {
         dl = m_allDSrigidbodies[in_RB]->get_distanceTo(*pt, rayDir, dh);
         geomVector rayVec(3);
         rayVec(0) = pt->operator()(0) + dl * rayDir(0);
         rayVec(1) = pt->operator()(1) + dl * rayDir(1);
         rayVec(2) = pt->operator()(2) + dl * rayDir(2);
         netVel = rigid_body_velocity(in_RB, rayVec);
         // rigid_body_velocity(in_RB, *pt + dl * rayDir);
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

      in_RB = (FF == UF) ? isIn_any_RB(parID, pt_at_face) : -1 ;

      if (in_RB != -1) {
         dr = m_allDSrigidbodies[in_RB]->get_distanceTo(*pt, rayDir, dh);
         geomVector rayVec(3);
         rayVec(0) = pt->operator()(0) + dr * rayDir(0);
         rayVec(1) = pt->operator()(1) + dr * rayDir(1);
         rayVec(2) = pt->operator()(2) + dr * rayDir(2);
         netVel = rigid_body_velocity(in_RB, rayVec);
         // rigid_body_velocity(in_RB, *pt + dr * rayDir);
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
   MAC_LABEL("DS_AllRigidBodies:: Bilinear_interpolation" ) ;

   size_t field = field_num(FF);
   double dh = MESH->get_smallest_grid_size();


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
            geomVector surface_point(3);// = *pt + del_wall(0,i)*rayDir;
            surface_point(0) = pt->operator()(0) + del_wall(0,i)*rayDir(0);
            surface_point(1) = pt->operator()(1) + del_wall(0,i)*rayDir(1);
            surface_point(2) = pt->operator()(2) + del_wall(0,i)*rayDir(2);
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
            geomVector surface_point(3);// = *pt + del_wall(0,i)*rayDir;
            surface_point(0) = pt->operator()(0) + del_wall(0,i)*rayDir(0);
            surface_point(1) = pt->operator()(1) + del_wall(0,i)*rayDir(1);
            surface_point(2) = pt->operator()(2) + del_wall(0,i)*rayDir(2);
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
         geomVector surface_point(3);// = *pt + del_wall(0,i)*rayDir;
         surface_point(0) = pt->operator()(0) + del_wall(0,i)*rayDir(0);
         surface_point(1) = pt->operator()(1) + del_wall(0,i)*rayDir(1);
         surface_point(2) = pt->operator()(2) + del_wall(0,i)*rayDir(2);
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
            geomVector surface_point(3);// = *pt + del_wall(1,j)*rayDir;
            surface_point(0) = pt->operator()(0) + del_wall(1,j)*rayDir(0);
            surface_point(1) = pt->operator()(1) + del_wall(1,j)*rayDir(1);
            surface_point(2) = pt->operator()(2) + del_wall(1,j)*rayDir(2);
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
            geomVector surface_point(3);// = *pt + del_wall(1,j)*rayDir;
            surface_point(0) = pt->operator()(0) + del_wall(1,j)*rayDir(0);
            surface_point(1) = pt->operator()(1) + del_wall(1,j)*rayDir(1);
            surface_point(2) = pt->operator()(2) + del_wall(1,j)*rayDir(2);
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
         geomVector surface_point(3);// = *pt + del_wall(1,j)*rayDir;
         surface_point(0) = pt->operator()(0) + del_wall(1,j)*rayDir(0);
         surface_point(1) = pt->operator()(1) + del_wall(1,j)*rayDir(1);
         surface_point(2) = pt->operator()(2) + del_wall(1,j)*rayDir(2);
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
double
DS_AllRigidBodies:: Biquadratic_interpolation ( FV_DiscreteField const* FF
                                             , size_t const& comp
                                             , geomVector const* pt
                                             , size_t_vector const& i0
                                             , size_t const& interpol_dir
                                             , int const& sign
                                             , vector<size_t> const& list)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_AllRigidBodies:: Biquadratic_interpolation" ) ;

// Calculates the field value at the ghost points
// near the particle boundary using the quadratic interpolation
// inspired from Johansen 1998;
// xp,yp,zp are the ghost point coordinated; interpol_dir is the direction
// in which the additional points will be used for quadratic interpolation

   size_t field = field_num(FF) ;

   // Directional index of point
   boolVector in_solid(3,0);
   boolVector in_domain(3,1);
   // Directional indexes of ghost points
   vector<size_t_vector> i0_ghost(3,i0);
   // Local node index of ghost points
   size_t_vector node_index(3,0);
   // Decide which scheme to use
   scheme_list scheme = quadratic;

   geomVector xi(3), fi(3);

   // Creating ghost points for quadratic interpolation
   intVector i0_temp(3,0);

   if (sign > 0) {
      i0_temp(0) = int(i0(interpol_dir));
      i0_temp(1) = int(i0(interpol_dir)) + 1;
      i0_temp(2) = int(i0(interpol_dir)) + 2;
   } else if (sign <= 0) {
      i0_temp(0) = int(i0(interpol_dir)) - 1;
      i0_temp(1) = int(i0(interpol_dir));
      i0_temp(2) = int(i0(interpol_dir)) + 1;
   }

   // Checking the ghost points in domain or not
   for (size_t l = 0; l < 3; l++) {
      if ((i0_temp(l) < 0) ||
          (i0_temp(l) >= (int)FF->get_local_nb_dof(comp,interpol_dir))) {
         in_domain(l) = 0;
      } else {
         in_domain(l) = 1;
      }
      i0_ghost[l](interpol_dir) = i0_temp(l);
   }

   // Assume all the ghost points in fluid
   // Storing the field values assuming all ghost points in fluid and domain
   // Check weather the ghost points are in solid or not; TRUE if they are
   for (size_t l = 0; l < 3; l++) {
      if (in_domain(l)) {
         xi(l) = FF->get_DOF_coordinate( i0_ghost[l](interpol_dir)
                                       , comp
                                       , interpol_dir);
         fi(l) = 0.;
         for (size_t level : list)
            fi(l) += FF->DOF_value( i0_ghost[l](0)
                                  , i0_ghost[l](1)
                                  , i0_ghost[l](2), comp, level );
         fi(l) /= (double)list.size();

         node_index(l) = FF->DOF_local_number( i0_ghost[l](0)
                                             , i0_ghost[l](1)
                                             , i0_ghost[l](2),comp);
         in_solid(l) = void_fraction[field]->operator()(node_index(l));
      }
   }

   // Ghost points corrections
   // All points in domain
   if (in_domain(0) && in_domain(1) && in_domain(2)) {
      // 0 in solid, rest in fluid
      if ((in_solid(0) != 0) &&
          (in_solid(1) == 0) &&
          (in_solid(2) == 0)) {
	      if (FF == UF) {
            xi(0) = FF->get_DOF_coordinate(i0_ghost[1](interpol_dir), comp, interpol_dir)
               - intersect_distance[field]->operator()(node_index(1),2*interpol_dir + 0);
            fi(0) = intersect_fieldValue[field]->operator()(node_index(1),2*interpol_dir + 0);
         } else {
            scheme = linear12;
         }
      // 2 in solid, rest in fluid
      } else if ((in_solid(0) == 0) &&
                 (in_solid(1) == 0) &&
                 (in_solid(2) != 0)) {
         if (FF == UF) {
            xi(2) = FF->get_DOF_coordinate(i0_ghost[1](interpol_dir), comp, interpol_dir)
               + intersect_distance[field]->operator()(node_index(1),2*interpol_dir + 1);
            fi(2) = intersect_fieldValue[field]->operator()(node_index(1),2*interpol_dir + 1);
         } else {
            scheme = linear01;
         }
      // 0, 2 in solid; 1 in fluid
      } else if ((in_solid(0) != 0) &&
                 (in_solid(1) == 0) &&
                 (in_solid(2) != 0)) {
         if (FF == UF) {
            xi(0) = FF->get_DOF_coordinate(i0_ghost[1](interpol_dir), comp, interpol_dir)
               - intersect_distance[field]->operator()(node_index(1),2*interpol_dir + 0);
            fi(0) = intersect_fieldValue[field]->operator()(node_index(1),2*interpol_dir + 0);
            xi(2) = FF->get_DOF_coordinate(i0_ghost[1](interpol_dir), comp, interpol_dir)
               + intersect_distance[field]->operator()(node_index(1),2*interpol_dir + 1);
            fi(2) = intersect_fieldValue[field]->operator()(node_index(1),2*interpol_dir + 1);
         } else {
            scheme = linear1;
         }
      // 0, 1 in solid; 2 in fluid
      } else if ((in_solid(0) != 0) &&
                 (in_solid(1) != 0) &&
                 (in_solid(2) == 0)) {
         if (FF == UF) {
            xi(1) = FF->get_DOF_coordinate(i0_ghost[2](interpol_dir), comp, interpol_dir)
               - intersect_distance[field]->operator()(node_index(2),2*interpol_dir + 0);
            fi(1) = intersect_fieldValue[field]->operator()(node_index(2),2*interpol_dir + 0);
            scheme = linear12;
         } else {
            scheme = linear2;
         }
      // 1, 2 in solid; 0 in fluid
      } else if ((in_solid(0) == 0) &&
                 (in_solid(1) != 0) &&
                 (in_solid(2) != 0)) {
         if (FF == UF) {
            xi(1) = FF->get_DOF_coordinate(i0_ghost[0](interpol_dir), comp, interpol_dir)
            + intersect_distance[field]->operator()(node_index(0),2*interpol_dir + 1);
            fi(1) = intersect_fieldValue[field]->operator()(node_index(0),2*interpol_dir + 1);
            scheme = linear01;
         } else {
            scheme = linear0;
         }
      }
   // Point 0 and 1 are in domain, 2 not in domain
   } else if (in_domain(0) && in_domain(1) && !in_domain(2)) {
      scheme = linear01;
      // 0 in fluid; 1 in solid
      if ((in_solid(0) == 0) &&
          (in_solid(1) != 0)) {
         if (FF == UF) {
            xi(1) = FF->get_DOF_coordinate(i0_ghost[0](interpol_dir), comp, interpol_dir)
               + intersect_distance[field]->operator()(node_index(0),2*interpol_dir + 1);
            fi(1) = intersect_fieldValue[field]->operator()(node_index(0),2*interpol_dir + 1);
         } else {
            scheme = linear0;
         }
      // 0 in solid, 1 in fluid
      } else if ((in_solid(0) != 0) &&
                 (in_solid(1) == 0)) {
        if (FF == UF) {
            xi(0) = FF->get_DOF_coordinate(i0_ghost[1](interpol_dir), comp, interpol_dir)
               - intersect_distance[field]->operator()(node_index(1),2*interpol_dir + 0);
            fi(0) = intersect_fieldValue[field]->operator()(node_index(1),2*interpol_dir + 0);
        } else {
            scheme = linear1;
        }
      }
   // Point 1 and 2 are in domain, 0 not in domain
   } else if (!in_domain(0) && in_domain(1) && in_domain(2)) {
      scheme = linear12;
      // 1 in fluid; 2 in solid
      if ((in_solid(1) == 0) &&
          (in_solid(2) != 0)) {
         if (FF == UF) {
            xi(2) = FF->get_DOF_coordinate(i0_ghost[1](interpol_dir), comp, interpol_dir)
               + intersect_distance[field]->operator()(node_index(1),2*interpol_dir + 1);
            fi(2) = intersect_fieldValue[field]->operator()(node_index(1),2*interpol_dir + 1);
         } else {
            scheme = linear1;
         }
      // 1 in solid, 2 in fluid
      } else if ((in_solid(1) != 0) &&
                 (in_solid(2) == 0)) {
         if (FF == UF) {
            xi(1) = FF->get_DOF_coordinate(i0_ghost[2](interpol_dir), comp, interpol_dir)
               - intersect_distance[field]->operator()(node_index(2),2*interpol_dir + 0);
            fi(1) = intersect_fieldValue[field]->operator()(node_index(2),2*interpol_dir + 0);
         } else {
            scheme = linear2;
         }
      }
   }

   double value = 0.;
   double l0=0.,l1=0.,l2=0.;

   switch(scheme) {
      case quadratic:
         l0 = (pt->operator()(interpol_dir) - xi(1))
             *(pt->operator()(interpol_dir) - xi(2))
             /(xi(0) - xi(1))/(xi(0) - xi(2));
         l1 = (pt->operator()(interpol_dir) - xi(0))
             *(pt->operator()(interpol_dir) - xi(2))
             /(xi(1) - xi(0))/(xi(1) - xi(2));
         l2 = (pt->operator()(interpol_dir) - xi(0))
             *(pt->operator()(interpol_dir) - xi(1))
             /(xi(2) - xi(0))/(xi(2) - xi(1));
         value = fi(0)*l0 + fi(1)*l1 + fi(2)*l2;
         break;
      case linear01:
         l0 = (pt->operator()(interpol_dir) - xi(1)) / (xi(0) - xi(1));
         l1 = (pt->operator()(interpol_dir) - xi(0)) / (xi(1) - xi(0));
         value = fi(0)*l0 + fi(1)*l1;
         break;
      case linear12:
         l1 = (pt->operator()(interpol_dir) - xi(2)) / (xi(1) - xi(2));
         l2 = (pt->operator()(interpol_dir) - xi(1)) / (xi(2) - xi(1));
         value = fi(1)*l1 + fi(2)*l2;
         break;
      case linear0:
         value = fi(0);
         break;
      case linear1:
         value = fi(1);
         break;
      case linear2:
         value = fi(2);
         break;
   }

   return(value);
}




//---------------------------------------------------------------------------
double
DS_AllRigidBodies:: Triquadratic_interpolation ( FV_DiscreteField const* FF
                                               , size_t const& comp
                                               , geomVector const* pt
                                               , size_t_vector const& i0
                                               , size_t const& parID
                                               , size_t const& ghost_points_dir
                                               , vector<int> const& sign
                                               , vector<size_t> const& list)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllRigidBodies:: Triquadratic_interpolation" ) ;

  geomVector point(3);
  // Directional indexes of ghost points
  size_t_vector index(i0);
  // Directional indexes of ghost points
  vector<size_t_vector> i0_ghost(3,i0);
  // Coordinates of secondary ghost points
  vector<geomVector> coord_g(3,0.);
  // Store particle ID if level_set becomes negative
  vector<int> in_parID(3,0);
  // Presence in domain or not
  boolVector in_domain(3,1);
  vector<double> net_vel(3,0.);
  // Decide which scheme to use
  scheme_list scheme = quadratic;

  size_t sec_ghost_dir = 0;
  size_t sec_interpol_dir = 0;

  point(0) = pt->operator()(0);
  point(1) = pt->operator()(1);
  point(2) = pt->operator()(2);

  // Ghost points generated in y and then quadratic interpolation
  // in z will generate the same stencil if ghost points are
  // generated in z and the quadratic interpolation done in y
  if (ghost_points_dir == 0) {
     sec_ghost_dir = 1;
     sec_interpol_dir = 2;
  } else if (ghost_points_dir == 1) {
     sec_ghost_dir = 0;
     sec_interpol_dir = 2;
  } else if (ghost_points_dir == 2) {
     sec_ghost_dir = 0;
     sec_interpol_dir = 1;
  }

  // Creating ghost points for quadratic interpolation
  intVector i0_temp(3,0);

  if (sign[sec_ghost_dir] > 0.) {
     i0_temp(0) = int(i0(sec_ghost_dir));
     i0_temp(1) = int(i0(sec_ghost_dir)) + 1;
     i0_temp(2) = int(i0(sec_ghost_dir)) + 2;
  } else if (sign[sec_ghost_dir] <= 0.) {
     i0_temp(0) = int(i0(sec_ghost_dir)) - 1;
     i0_temp(1) = int(i0(sec_ghost_dir));
     i0_temp(2) = int(i0(sec_ghost_dir)) + 1;
  }

  // Checking the ghost points in domain or not; loop on the ghost points
  for (size_t l = 0; l < 3; l++) {
     if ((i0_temp(l) < 0) ||
         (i0_temp(l) >= (int)FF->get_local_nb_dof(comp,sec_ghost_dir))) {
        in_domain(l) = false;
     } else {
        in_domain(l) = true;
     }
     i0_ghost[l](sec_ghost_dir) = i0_temp(l);
  }

  // Assume all secondary ghost points in fluid
  double x0 = in_domain(0) ?
      FF->get_DOF_coordinate(i0_ghost[0](sec_ghost_dir), comp, sec_ghost_dir)
                           : 0. ;
  double x1 = in_domain(1) ?
      FF->get_DOF_coordinate(i0_ghost[1](sec_ghost_dir), comp, sec_ghost_dir)
                           : 0. ;
  double x2 = in_domain(2) ?
      FF->get_DOF_coordinate(i0_ghost[2](sec_ghost_dir), comp, sec_ghost_dir)
                           : 0. ;

  if (sec_ghost_dir == 0) {
     coord_g[0](0) = x0; coord_g[0](1) = point(1); coord_g[0](2) = point(2);
     coord_g[1](0) = x1; coord_g[1](1) = point(1); coord_g[1](2) = point(2);
     coord_g[2](0) = x2; coord_g[2](1) = point(1); coord_g[2](2) = point(2);
  } else if (sec_ghost_dir == 1) {
     coord_g[0](0) = point(0); coord_g[0](1) = x0; coord_g[0](2) = point(2);
     coord_g[1](0) = point(0); coord_g[1](1) = x1; coord_g[1](2) = point(2);
     coord_g[2](0) = point(0); coord_g[2](1) = x2; coord_g[2](2) = point(2);
  } else if (sec_ghost_dir == 2) {
     coord_g[0](0) = point(0); coord_g[0](1) = point(1); coord_g[0](2) = x0;
     coord_g[1](0) = point(0); coord_g[1](1) = point(1); coord_g[1](2) = x1;
     coord_g[2](0) = point(0); coord_g[2](1) = point(1); coord_g[2](2) = x2;
  }

  double dh = MESH->get_smallest_grid_size();



  // Stores rigid body ID if inside solid; -1 otherwise
  in_parID[0] = isIn_any_RB(parID, coord_g[0]);
  in_parID[1] = isIn_any_RB(parID, coord_g[1]);
  in_parID[2] = isIn_any_RB(parID, coord_g[2]);

  double del = 0.;

  // Estimate the field values at the secondary ghost points
  double f0 = (in_domain(0)) ?
               Biquadratic_interpolation(FF
                                       , comp
                                       , &coord_g[0]
                                       , i0_ghost[0]
                                       , sec_interpol_dir
                                       , sign[sec_interpol_dir],list)
               : 0.;
  double f1 = (in_domain(1)) ?
               Biquadratic_interpolation(FF
                                       , comp
                                       , &coord_g[1]
                                       , i0_ghost[1]
                                       , sec_interpol_dir
                                       , sign[sec_interpol_dir],list)
               : 0.;
  double f2 = (in_domain(2)) ?
               Biquadratic_interpolation(FF
                                       , comp
                                       , &coord_g[2]
                                       , i0_ghost[2]
                                       , sec_interpol_dir
                                       , sign[sec_interpol_dir],list)
               : 0.;

  // Ghost points corrections
  if (in_domain(0) && in_domain(1) && in_domain(2)) {
     // 0 in solid, rest in fluid
     if ((in_parID[0] != -1) &&
         (in_parID[1] == -1) &&
         (in_parID[2] == -1)) {
        geomVector rayDir(3);
        rayDir(sec_ghost_dir) = -1. ;
        if (FF == UF) {
           del = m_allDSrigidbodies[in_parID[0]]->
                                       get_distanceTo(coord_g[1], rayDir, dh);
           geomVector rayVec(3);
           rayVec(0) = coord_g[1](0) + del * rayDir(0);
           rayVec(1) = coord_g[1](1) + del * rayDir(1);
           rayVec(2) = coord_g[1](2) + del * rayDir(2);
           geomVector netVel =
                  rigid_body_velocity(in_parID[0], rayVec);
                  // rigid_body_velocity(in_parID[0], coord_g[1] + del * rayDir);

           x0 = x1 - del;
           f0 = netVel(comp);
	     } else {
           scheme = linear12;
	     }
     // 2 in solid, rest in fluid
     } else if ((in_parID[0] == -1) &&
                (in_parID[1] == -1) &&
                (in_parID[2] != -1)) {
         geomVector rayDir(3);
         rayDir(sec_ghost_dir) = 1. ;
         if (FF == UF) {
            del = m_allDSrigidbodies[in_parID[2]]->
                                        get_distanceTo(coord_g[1], rayDir, dh);

            geomVector rayVec(3);
            rayVec(0) = coord_g[1](0) + del * rayDir(0);
            rayVec(1) = coord_g[1](1) + del * rayDir(1);
            rayVec(2) = coord_g[1](2) + del * rayDir(2);
            geomVector netVel =
                   rigid_body_velocity(in_parID[2], rayVec);
                   // rigid_body_velocity(in_parID[2], coord_g[1] + del * rayDir);

            x2 = x1 + del;
            f2 = netVel(comp);
   	   } else {
            scheme = linear01;
	      }
      // 0, 2 in solid; 1 in fluid
      } else if ((in_parID[0] != -1) &&
                 (in_parID[1] == -1) &&
                 (in_parID[2] != -1)) {
         geomVector rayDir(3);
         rayDir(sec_ghost_dir) = -1. ;
         if (FF == UF) {
            del = m_allDSrigidbodies[in_parID[0]]->
                                        get_distanceTo(coord_g[1], rayDir, dh);

            geomVector rayVec(3);
            rayVec(0) = coord_g[1](0) + del * rayDir(0);
            rayVec(1) = coord_g[1](1) + del * rayDir(1);
            rayVec(2) = coord_g[1](2) + del * rayDir(2);
            geomVector netVel =
                   rigid_body_velocity(in_parID[0], rayVec);
                   // rigid_body_velocity(in_parID[0], coord_g[1] + del * rayDir);

            x0 = x1 - del;
            f0 = netVel(comp);

            rayDir(sec_ghost_dir) = 1. ;

            del = m_allDSrigidbodies[in_parID[2]]->
                                        get_distanceTo(coord_g[1], rayDir, dh);

            rayVec(0) = coord_g[1](0) + del * rayDir(0);
            rayVec(1) = coord_g[1](1) + del * rayDir(1);
            rayVec(2) = coord_g[1](2) + del * rayDir(2);
            netVel =
                   rigid_body_velocity(in_parID[2], rayVec);
                   // rigid_body_velocity(in_parID[2], coord_g[1] + del * rayDir);

            x2 = x1 + del;
            f2 = netVel(comp);
	      } else {
	         scheme = linear1;
         }
      // 0, 1 in solid; 2 in fluid
      } else if ((in_parID[0] != -1) &&
                 (in_parID[1] != -1) &&
                 (in_parID[2] == -1)) {
         geomVector rayDir(3);
         rayDir(sec_ghost_dir) = -1. ;
         if (FF == UF) {
            del = m_allDSrigidbodies[in_parID[1]]->
                                        get_distanceTo(coord_g[2], rayDir, dh);

            geomVector rayVec(3);
            rayVec(0) = coord_g[2](0) + del * rayDir(0);
            rayVec(1) = coord_g[2](1) + del * rayDir(1);
            rayVec(2) = coord_g[2](2) + del * rayDir(2);
            geomVector netVel =
                   rigid_body_velocity(in_parID[1], rayVec);
                   // rigid_body_velocity(in_parID[1], coord_g[2] + del * rayDir);

            x1 = x2 - del;
            f1 = netVel(comp);
            scheme = linear12;
         } else {
            scheme = linear2;
         }
      // 1, 2 in solid; 0 in fluid
      } else if ((in_parID[0] == -1) &&
                 (in_parID[1] != -1) &&
                 (in_parID[2] != -1)) {
         geomVector rayDir(3);
         rayDir(sec_ghost_dir) = 1. ;
         if (FF == UF) {
            del = m_allDSrigidbodies[in_parID[1]]->
                                        get_distanceTo(coord_g[0], rayDir, dh);

            geomVector rayVec(3);
            rayVec(0) = coord_g[0](0) + del * rayDir(0);
            rayVec(1) = coord_g[0](1) + del * rayDir(1);
            rayVec(2) = coord_g[0](2) + del * rayDir(2);
            geomVector netVel =
                   rigid_body_velocity(in_parID[1], rayVec);
                   // rigid_body_velocity(in_parID[1], coord_g[0] + del * rayDir);

            x1 = x0 + del;
            f1 = netVel(comp);
            scheme = linear01;
         } else {
            scheme = linear0;
         }
      }
   // Point 0 and 1 are in domain, 2 not in domain
   } else if (in_domain(0) && in_domain(1) && !in_domain(2)) {
      scheme = linear01;
      // 0 in fluid; 1 in solid
      if ((in_parID[0] == -1) &&
          (in_parID[1] != -1)) {
          geomVector rayDir(3);
          rayDir(sec_ghost_dir) = 1. ;
          if (FF == UF) {
             del = m_allDSrigidbodies[in_parID[1]]->
                                        get_distanceTo(coord_g[0], rayDir, dh);

             geomVector rayVec(3);
             rayVec(0) = coord_g[0](0) + del * rayDir(0);
             rayVec(1) = coord_g[0](1) + del * rayDir(1);
             rayVec(2) = coord_g[0](2) + del * rayDir(2);
             geomVector netVel =
                   rigid_body_velocity(in_parID[1], rayVec);
                   // rigid_body_velocity(in_parID[1], coord_g[0] + del * rayDir);

             x1 = x0 + del;
             f1 = netVel(comp);
	       } else {
             scheme = linear0;
          }
       // 0 in solid, 1 in fluid
       } else if ((in_parID[0] != -1) &&
                  (in_parID[1] == -1)) {
          geomVector rayDir(3);
          rayDir(sec_ghost_dir) = -1. ;
          if (FF == UF) {
             del = m_allDSrigidbodies[in_parID[0]]->
                                        get_distanceTo(coord_g[1], rayDir, dh);

             geomVector rayVec(3);
             rayVec(0) = coord_g[1](0) + del * rayDir(0);
             rayVec(1) = coord_g[1](1) + del * rayDir(1);
             rayVec(2) = coord_g[1](2) + del * rayDir(2);
             geomVector netVel =
                   rigid_body_velocity(in_parID[0], rayVec);
                   // rigid_body_velocity(in_parID[0], coord_g[1] + del * rayDir);

             x0 = x1 - del;
             f0 = netVel(comp);
      	 } else {
      	    scheme = linear1;
      	 }
       }
   // Point 1 and 2 are in domain, 0 not in domain
   } else if (!in_domain(0) && in_domain(1) && in_domain(2)) {
       scheme = linear12;
       // 1 in fluid; 2 in solid
       if ((in_parID[1] == -1) &&
           (in_parID[2] != -1)) {
           geomVector rayDir(3);
           rayDir(sec_ghost_dir) = 1. ;
           if (FF == UF) {
              del = m_allDSrigidbodies[in_parID[2]]->
                                         get_distanceTo(coord_g[1], rayDir, dh);

              geomVector rayVec(3);
              rayVec(0) = coord_g[1](0) + del * rayDir(0);
              rayVec(1) = coord_g[1](1) + del * rayDir(1);
              rayVec(2) = coord_g[1](2) + del * rayDir(2);
              geomVector netVel =
                    rigid_body_velocity(in_parID[2], rayVec);
                    // rigid_body_velocity(in_parID[2], coord_g[1] + del * rayDir);

              x2 = x1 + del;
              f2 = netVel(comp);
           } else {
              scheme = linear1;
           }
       // 1 in solid, 2 in fluid
       } else if ((in_parID[1] != -1) &&
                  (in_parID[2] == -1)) {
           geomVector rayDir(3);
           rayDir(sec_ghost_dir) = -1. ;
           if (FF == UF) {
              del = m_allDSrigidbodies[in_parID[1]]->
                                        get_distanceTo(coord_g[2], rayDir, dh);

              geomVector rayVec(3);
              rayVec(0) = coord_g[2](0) + del * rayDir(0);
              rayVec(1) = coord_g[2](1) + del * rayDir(1);
              rayVec(2) = coord_g[2](2) + del * rayDir(2);
              geomVector netVel =
                   rigid_body_velocity(in_parID[1], rayVec);
                   // rigid_body_velocity(in_parID[1], coord_g[2] + del * rayDir);

              x1 = x2 - del;
              f1 = netVel(comp);
           } else {
              scheme = linear2;
           }
       }
   }

  double l0 = 0., l1 = 0., l2 = 0.;
  double result = 0.;

  switch (scheme) {
     case quadratic:
        l0 = (point(sec_ghost_dir) - x1)
           * (point(sec_ghost_dir) - x2)
           / (x0 - x1) / (x0 - x2);
        l1 = (point(sec_ghost_dir) - x0)
           * (point(sec_ghost_dir) - x2)
           / (x1 - x0) / (x1 - x2);
        l2 = (point(sec_ghost_dir) - x0)
           * (point(sec_ghost_dir) - x1)
           / (x2 - x0) / (x2 - x1);
        result = f0*l0 + f1*l1 + f2*l2;
        break;
     case linear01:
        l0 = (point(sec_ghost_dir) - x1)/(x0 - x1);
        l1 = (point(sec_ghost_dir) - x0)/(x1 - x0);
        result = f0*l0 + f1*l1;
        break;
     case linear12:
        l1 = (point(sec_ghost_dir) - x2)/(x1 - x2);
        l2 = (point(sec_ghost_dir) - x1)/(x2 - x1);
        result = f1*l1 + f2*l2;
        break;
     case linear0:
        result = f0;
        break;
     case linear1:
        result = f1;
        break;
     case linear2:
        result = f2;
        break;
  }

  return(result);
}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_surface_variables_for_all_RB( )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_AllRigidBodies:: compute_surface_variables_for_all_RB" ) ;

   for (size_t i = 0; i < m_nrb; ++i) {
      m_allDSrigidbodies[i]->compute_surface_points( );
      m_allDSrigidbodies[i]->correct_surface_discretization( MESH );
   }

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
      m_allDSrigidbodies[i]->write_surface_discretization( file, m_macCOMM );
   }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: build_solid_variables_on_fluid_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_AllRigidBodies:: build_solid_variables_on_fluid_grid" ) ;

   int field = field_num(FF);

   size_t FF_LOC_UNK = FF->nb_local_unknowns();

   // void fraction on the computational grid
   void_fraction[field]->re_initialize(FF_LOC_UNK);

   // Intersection parameters on the computational grid
   // For PF and UF
   intersect_vector[field]->re_initialize(FF_LOC_UNK,6);
   intersect_distance[field]->re_initialize(FF_LOC_UNK,6);
   intersect_fieldValue[field]->re_initialize(FF_LOC_UNK,6);

}




//---------------------------------------------------------------------------
int DS_AllRigidBodies:: field_num( FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  int number = 0;
  if (FF == PF) {
     number = 0;
  } else if (FF == UF) {
     number = 1;
  } else if (FF == TF) {
     number = 2;
  }

  return (number);

}




//---------------------------------------------------------------------------
size_t_vector* DS_AllRigidBodies:: get_void_fraction_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = field_num(FF);

  return (void_fraction[field]);

}




//---------------------------------------------------------------------------
size_t_array2D* DS_AllRigidBodies:: get_intersect_vector_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = field_num(FF);

  return (intersect_vector[field]);

}




//---------------------------------------------------------------------------
doubleArray2D* DS_AllRigidBodies:: get_intersect_distance_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = field_num(FF);

  return (intersect_distance[field]);

}




//---------------------------------------------------------------------------
doubleArray2D* DS_AllRigidBodies:: get_intersect_fieldValue_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = field_num(FF);

  return (intersect_fieldValue[field]);

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: create_neighbour_list_for_AllRB( )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_AllRigidBodies:: create_neighbour_list_for_AllRB()" ) ;

   boolVector const* periodic_comp = MESH->get_periodic_directions();

   double dh = MESH->get_smallest_grid_size();

   for (size_t i = 0; i < m_nrb; ++i) {
      vector<size_t> v1;
      v1.push_back(i);
      double r_i = m_allDSrigidbodies[i]->get_circumscribed_radius();
      geomVector const* pgc_i = m_allDSrigidbodies[i]
                                              ->get_ptr_to_gravity_centre();
      for (size_t j = 0; j < m_nrb; ++j) {
         double r_j = m_allDSrigidbodies[j]->get_circumscribed_radius();
         geomVector const* pgc_j = m_allDSrigidbodies[j]
                                                 ->get_ptr_to_gravity_centre();

         size_t dim = MESH->nb_space_dimensions() ;

         geomVector delta(3);

         for (size_t dir = 0; dir < dim; dir++) {
            delta(dir) = pgc_j->operator()(dir) - pgc_i->operator()(dir);

            bool is_periodic = periodic_comp->operator()( dir );

            if (is_periodic) {
               double isize = MESH->get_main_domain_max_coordinate(dir)
                            - MESH->get_main_domain_min_coordinate(dir);
               delta(dir) = delta(dir) - round(delta(dir)/isize) * isize;
            }
         }
         double dist = delta.calcNorm();

         // Store if two surfaces are less than 10 grid cells away
         if (dist <= (r_i + r_j + 10*dh))
            v1.push_back(j);
      }
      if (!v1.empty()) {
         neighbour_list.push_back(v1);
      } else {
         neighbour_list.push_back({0});
      }
   }
}
