#include <DS_RigidBody.hh>
#include <FS_RigidBody.hh>
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
tuple<bool,double,size_t> DS_RigidBody:: distanceTo( geomVector const& pt, 
      	size_t const& direction,
      	bool const& positive ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: distanceTo( pt, direction, positive )" ) ;

  return ( m_geometric_rigid_body->distanceTo( pt, direction, positive ) );
  
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
