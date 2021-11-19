#include <DS_Sphere.hh>
#include <FS_RigidBody.hh>
#include <FS_Sphere.hh>
#include <FV_DiscreteField.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_Sphere:: DS_Sphere()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_Sphere:: DS_Sphere" ) ;

}




//---------------------------------------------------------------------------
DS_Sphere:: DS_Sphere( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

  m_halo_zone = new doubleArray2D (3,2,0);

}




//---------------------------------------------------------------------------
DS_Sphere:: ~DS_Sphere()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: ~DS_Sphere" ) ;

}




//---------------------------------------------------------------------------
void DS_Sphere:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_Sphere:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
doubleArray2D* DS_Sphere:: compute_rigid_body_halozone( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: compute_rigid_body_halozone" ) ;

  struct FS_Sphere_Additional_Param const* pagp =
   dynamic_cast<FS_Sphere*>(m_geometric_rigid_body)
      ->get_ptr_FS_Sphere_Additional_Param();

  geomVector const* pgs = dynamic_cast<FS_Sphere*>(m_geometric_rigid_body)
                              ->get_ptr_FS_Sphere_gravity_centre();

  m_halo_zone->operator()(0,0) = pgs->operator()(0) - 1.5*pagp->radius;
  m_halo_zone->operator()(0,1) = pgs->operator()(0) + 1.5*pagp->radius;

  m_halo_zone->operator()(1,0) = pgs->operator()(1) - 1.5*pagp->radius;
  m_halo_zone->operator()(1,1) = pgs->operator()(1) + 1.5*pagp->radius;

  m_halo_zone->operator()(2,0) = pgs->operator()(2) - 1.5*pagp->radius;
  m_halo_zone->operator()(2,1) = pgs->operator()(2) + 1.5*pagp->radius;

  return(m_halo_zone);

}




//---------------------------------------------------------------------------
void DS_Sphere:: compute_hydro_force_torque( FV_DiscreteField const* PP,
	FV_DiscreteField const* UU )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: compute_hydro_force_torque" ) ;

  // Access additional geometric parameters of the sphere
  struct FS_Sphere_Additional_Param const* pagp =
  	dynamic_cast<FS_Sphere*>(m_geometric_rigid_body)
		->get_ptr_FS_Sphere_Additional_Param();
  MAC::out() << pagp->radius << endl; // example, delete later

  // Determine the number of surface points and allocate the vector if empty
  if ( m_surface_points.empty() )
  {
    size_t npts = 0, i;
    geomVector vvv(3);

    npts = 1; // NEEDS TO BE SET TO THE PROPER VALUE
    m_surface_points.reserve( npts );
    for (i = 0; i < npts; ++i) m_surface_points.push_back( vvv );
  }

  // Determine the surface point coordinates
  MAC::out() << "DS_Sphere:: compute_hydro_force_torque - "
  	"Determining the set of surface points requires programming" << endl;

  // Compute the surface integrals and stores the force & torque value in
  // the corresponding geometric rigid body
  compute_surface_integrals_hydro_force_torque( PP, UU );

}
