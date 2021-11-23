#include <DS_3Dcylinder.hh>
#include <FS_RigidBody.hh>
#include <FS_3Dcylinder.hh>
#include <FV_DiscreteField.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_3Dcylinder:: DS_3Dcylinder()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_3Dcylinder:: DS_3Dcylinder" ) ;

}




//---------------------------------------------------------------------------
DS_3Dcylinder:: DS_3Dcylinder( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_3Dcylinder:: ~DS_3Dcylinder()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: ~DS_3Dcylinder" ) ;

}




//---------------------------------------------------------------------------
void DS_3Dcylinder:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_3Dcylinder:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
void DS_3Dcylinder:: compute_rigid_body_halozone( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: compute_rigid_body_halozone" ) ;

  struct FS_3Dcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_3Dcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_3Dcylinder_Additional_Param();

  geomVector const* pgs = dynamic_cast<FS_3Dcylinder*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  double radius_equivalent = 3.0 *
                     MAC::sqrt(pagp->cylinder_radius * pagp->cylinder_radius
                             + pagp->cylinder_height * pagp->cylinder_height);

  m_halo_zone[0]->operator()(0) = pgs->operator()(0) - radius_equivalent;
  m_halo_zone[0]->operator()(1) = pgs->operator()(1) - radius_equivalent;
  m_halo_zone[0]->operator()(2) = pgs->operator()(2) - radius_equivalent;

  m_halo_zone[1]->operator()(0) = pgs->operator()(0) + radius_equivalent;
  m_halo_zone[1]->operator()(1) = pgs->operator()(1) + radius_equivalent;
  m_halo_zone[1]->operator()(2) = pgs->operator()(2) + radius_equivalent;


}




//---------------------------------------------------------------------------
void DS_3Dcylinder:: compute_surface_points( size_t const& Np )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: compute_surface_points" ) ;



}





//---------------------------------------------------------------------------
void DS_3Dcylinder:: initialize_variable_for_each_rigidBody( size_t const& Np )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: initialize_variable_for_each_rigidBody" ) ;

  size_t Ntot = 2*Np;

  m_surface_points.reserve( Ntot );
  m_surface_area.reserve( Ntot );
  m_surface_normal.reserve( Ntot );
   for (size_t i = 0; i < Ntot; ++i) {
      m_surface_points.push_back( new geomVector(3) );
      m_surface_area.push_back( new geomVector(3) );
      m_surface_normal.push_back( new geomVector(3) );
   }

}




//---------------------------------------------------------------------------
void DS_3Dcylinder:: compute_hydro_force_torque( FV_DiscreteField const* PP,
	FV_DiscreteField const* UU )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: compute_hydro_force_torque" ) ;

  // Access additional geometric parameters of the 3D cylinder
  struct FS_3Dcylinder_Additional_Param const* pagp =
  	dynamic_cast<FS_3Dcylinder*>(m_geometric_rigid_body)
		->get_ptr_FS_3Dcylinder_Additional_Param();
  MAC::out() << pagp->cylinder_height << endl; // example, delete later

  // Determine the number of surface points and allocate the vector if empty
  // if ( m_surface_points.empty() )
  // {
  //   size_t npts = 0, i;
  //   geomVector vvv(3);
  //
  //   npts = 1; // NEEDS TO BE SET TO THE PROPER VALUE
  //   m_surface_points.reserve( npts );
  //   for (i = 0; i < npts; ++i) m_surface_points.push_back( vvv );
  // }

  // Determine the surface point coordinates
  MAC::out() << "DS_3Dcylinder:: compute_hydro_force_torque - "
  	"Determining the set of surface points requires programming" << endl;

  // Compute the surface integrals and stores the force & torque value in
  // the corresponding geometric rigid body
  compute_surface_integrals_hydro_force_torque( PP, UU );

}
