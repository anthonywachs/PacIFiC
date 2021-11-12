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
DS_AllRigidBodies:: DS_AllRigidBodies( size_t& dimens, istream& in,
	bool const& b_particles_as_fixed_obstacles )
//--------------------------------------------------------------------------- 
  : m_space_dimension( dimens )
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
