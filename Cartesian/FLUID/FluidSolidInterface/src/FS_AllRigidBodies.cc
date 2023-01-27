#include <FS_AllRigidBodies.hh>
#include <FS_RigidBody.hh>
#include <FS_RigidBody_BuilderFactory.hh>
using std::endl;


//---------------------------------------------------------------------------
FS_AllRigidBodies:: FS_AllRigidBodies()
//---------------------------------------------------------------------------
  : m_space_dimension( 3 )
  , m_npart( 0 )
  , m_nrb( 0 )
{
  MAC_LABEL( "FS_AllRigidBodies:: FS_AllRigidBodies" ) ;

}




//---------------------------------------------------------------------------
FS_AllRigidBodies:: FS_AllRigidBodies( size_t& dimens, istream& in,
	bool const& b_particles_as_fixed_obstacles )
//---------------------------------------------------------------------------
  : m_space_dimension( dimens )
{
  MAC_LABEL( "FS_AllRigidBodies:: FS_AllRigidBodies(size_t&,istream&)" ) ;

  FS_RigidBody* prb = NULL;
  size_t ncorners = 0, pIdsolidsolver = 0;

  // Read the total number of rigid bodies
  in >> m_nrb;

  // Allocate the vector of all rigid bodies
  m_allrigidbodies.reserve( m_nrb );

  // Read the input stream and create the rigid bodies
  m_npart = 0;
  for (size_t i = 0; i < m_nrb; i++)
  {
    in >> pIdsolidsolver >> ncorners;
    m_allrigidbodies.push_back( prb );
    m_allrigidbodies[i] = FS_RigidBody_BuilderFactory::create(
    	m_space_dimension, in, ncorners, i );

    if ( b_particles_as_fixed_obstacles )
    {
      m_allrigidbodies[i]->nullify_velocity();
      m_allrigidbodies[i]->change_from_particle_to_obstacle();
    }

    if ( m_allrigidbodies[i]->get_type() != "O" &&
    	m_allrigidbodies[i]->get_type() != "PO" ) ++m_npart;
  }

}




//---------------------------------------------------------------------------
FS_AllRigidBodies:: ~FS_AllRigidBodies()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_AllRigidBodies:: ~FS_AllRigidBodies" ) ;

  for (size_t i=0;i<m_nrb;i++) delete m_allrigidbodies[i];
  m_allrigidbodies.clear();

}




//---------------------------------------------------------------------------
size_t FS_AllRigidBodies:: get_number_rigid_bodies() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_AllRigidBodies:: get_number_rigid_bodies" ) ;

  return ( m_nrb );

}




//---------------------------------------------------------------------------
size_t FS_AllRigidBodies:: get_number_particles() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_AllRigidBodies:: get_number_particles" ) ;

  return ( m_npart );

}




//---------------------------------------------------------------------------
FS_RigidBody* FS_AllRigidBodies:: get_ptr_rigid_body( size_t i ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_AllRigidBodies:: get_ptr_rigid_body" ) ;

  return ( m_allrigidbodies[i] );

}




//---------------------------------------------------------------------------
void FS_AllRigidBodies:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_AllRigidBodies:: update" ) ;

  // IMPORTANT: we assume that the number of rigid bodies, the numbering
  // of rigid bodies and the type of rigid bodies is unchanged over time

  size_t ncorners = 0, pIdsolidsolver = 0;

  // Read the total number of rigid bodies
  in >> m_nrb;

  // Read the input stream and update the rigid bodies
  m_npart = 0;
  for (size_t i = 0; i < m_nrb; i++)
  {
    in >> pIdsolidsolver >> ncorners;
    m_allrigidbodies[i]->update( in );

    if ( m_allrigidbodies[i]->get_type() != "O" ) ++m_npart;
  }
}




//---------------------------------------------------------------------------
void FS_AllRigidBodies:: display( ostream& out,
	size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_AllRigidBodies:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Features of all geometric rigid bodies" << endl;
  out << space << three << "Space dimension = " << m_space_dimension << endl;
  out << space << three << "Total number of rigid bodies = " << m_nrb << endl;
  out << space << three << "Number of particles = " << m_npart << endl;
  out << space << three << "Number of obstacles = " << m_nrb - m_npart << endl;
  for (size_t i = 0; i < m_nrb; i++)
  {
    out << endl;
    out << space << three << "Rigid body " << i << endl;
    m_allrigidbodies[i]->display( out, indent_width + 6 );
  }

}




//---------------------------------------------------------------------------
void FS_AllRigidBodies:: nullify_velocity()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_AllRigidBodies:: nullify_velocity" ) ;

  for (size_t i = 0; i < m_nrb; i++)
    m_allrigidbodies[i]->nullify_velocity();

}
