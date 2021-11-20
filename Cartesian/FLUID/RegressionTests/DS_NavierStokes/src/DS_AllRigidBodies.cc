#include <DS_AllRigidBodies.hh>
#include <FS_AllRigidBodies.hh>
#include <DS_RigidBody.hh>
#include <DS_RigidBody_BuilderFactory.hh>
#include <FV_DiscreteField.hh>
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
                                  , FV_DiscreteField const* arb_PF )
//---------------------------------------------------------------------------
  : m_space_dimension( dimens )
  , UF ( arb_UF )
  , PF ( arb_PF )
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

  compute_void_fraction_on_grid();

  compute_grid_intersection_with_rigidbody();

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
void DS_AllRigidBodies:: compute_void_fraction_on_grid( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_void_fraction_on_grid" ) ;

  for (size_t i = 0; i < m_nrb; ++i) {
     m_allDSrigidbodies[i]->
            compute_void_fraction_on_grid(PF,void_fraction[0],rb_ID[0],i);
     m_allDSrigidbodies[i]->
            compute_void_fraction_on_grid(UF,void_fraction[1],rb_ID[1],i);
  }

}




//---------------------------------------------------------------------------
void DS_AllRigidBodies:: compute_grid_intersection_with_rigidbody( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllRigidBodies:: compute_grid_intersection_with_rigidbody" ) ;

  for (size_t i = 0; i < m_nrb; ++i) {
     m_allDSrigidbodies[i]->
      compute_grid_intersection_with_rigidbody(PF
                                              ,void_fraction[0]
                                              ,intersect_vector[0]
                                              ,intersect_distance[0]
                                              ,intersect_fieldValue[0]);
     m_allDSrigidbodies[i]->
      compute_grid_intersection_with_rigidbody(UF
                                              ,void_fraction[1]
                                              ,intersect_vector[1]
                                              ,intersect_distance[1]
                                              ,intersect_fieldValue[1]);
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

   // ID of rigid body on the computational grid, if present
   // For PF and UF
   rb_ID.reserve(2);
   rb_ID.push_back(new size_t_vector(1,0));
   rb_ID.push_back(new size_t_vector(1,0));
   rb_ID[0]->re_initialize(PF_LOC_UNK);
   rb_ID[1]->re_initialize(UF_LOC_UNK);

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
size_t_vector* DS_AllRigidBodies:: get_rigidbodyIDs_on_grid(
                                                FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
  size_t field = (FF == PF) ? 0 : 1;

  return (rb_ID[field]);

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
