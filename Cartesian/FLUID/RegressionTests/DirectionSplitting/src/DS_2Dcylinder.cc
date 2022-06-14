#include <DS_2Dcylinder.hh>
#include <FS_RigidBody.hh>
#include <FS_2Dcylinder.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_2Dcylinder:: DS_2Dcylinder()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_2Dcylinder:: DS_2Dcylinder" ) ;

}




//---------------------------------------------------------------------------
DS_2Dcylinder:: DS_2Dcylinder( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_2Dcylinder:: ~DS_2Dcylinder()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dcylinder:: ~DS_2Dcylinder" ) ;

}




//---------------------------------------------------------------------------
void DS_2Dcylinder:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dcylinder:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_2Dcylinder:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dcylinder:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
void DS_2Dcylinder:: compute_rigid_body_halozone( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dcylinder:: compute_rigid_body_halozone" ) ;

  struct FS_2Dcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_2Dcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_2Dcylinder_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  double r_equi = 3.0*pagp->radius;

  geomVector delta(r_equi, r_equi, r_equi);

  m_halo_zone[0]->operator=(*pgc);
  m_halo_zone[1]->operator=(*pgc);

  m_halo_zone[0]->operator-=(delta);
  m_halo_zone[1]->operator+=(delta);

}




//---------------------------------------------------------------------------
void DS_2Dcylinder:: compute_surface_points( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dcylinder:: compute_surface_points" ) ;

  // Pointers to location and additional parameters
  struct FS_2Dcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_2Dcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_2Dcylinder_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                            ->get_ptr_to_gravity_centre();

  size_t Npoints = m_surface_area.size();
  double d_theta = 2.*MAC::pi()/((double)Npoints);
  double theta = 0.5*d_theta;

  for (size_t i = 0; i < Npoints; i++) {
     theta = theta + d_theta;

     geomVector point( pagp->radius*MAC::cos(theta)
                     , pagp->radius*MAC::sin(theta)
                     , 0. );

     m_surface_points[i]->operator=(point);
     m_surface_area[i]->operator()(0) = pagp->radius*d_theta;

     // Create surface normal vectors
     m_surface_normal[i]->operator=(point);
  }

  // Translate and rotate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
//     m_geometric_rigid_body->rotate(m_surface_points[i]);
//     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_surface_points[i]->translate(*pgc);
  }


}




//---------------------------------------------------------------------------
void DS_2Dcylinder:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dcylinder:: compute_number_of_surface_variables" ) ;

  struct FS_2Dcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_2Dcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_2Dcylinder_Additional_Param();

  size_t temp = (size_t) ((1./surface_cell_scale)
               *(2.*MAC::pi()*pagp->radius)
               /(dx));

  // Getting the nearest even number
  Ntot = (size_t) (round((double)temp * 0.5) * 2.);

}
