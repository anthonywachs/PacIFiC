#include <DS_3Dbox.hh>
#include <FS_RigidBody.hh>
#include <FS_3Dbox.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_3Dbox:: DS_3Dbox()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_3Dbox:: DS_3Dbox" ) ;

}




//---------------------------------------------------------------------------
DS_3Dbox:: DS_3Dbox( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_3Dbox:: ~DS_3Dbox()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: ~DS_3Dbox" ) ;

}




//---------------------------------------------------------------------------
void DS_3Dbox:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_3Dbox:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
void DS_3Dbox:: compute_rigid_body_halozone( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: compute_rigid_body_halozone" ) ;

  struct FS_3Dbox_Additional_Param const* pagp =
   dynamic_cast<FS_3Dbox*>(m_geometric_rigid_body)
      ->get_ptr_FS_3Dbox_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_3Dbox*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  geomVector delta = pagp->corners[0] - *pgc;

  double r_equi = 3.0 * delta.calcNorm();

  delta(0) = r_equi;
  delta(1) = r_equi;
  delta(2) = r_equi;

  m_halo_zone[0]->operator=(*pgc);
  m_halo_zone[1]->operator=(*pgc);

  m_halo_zone[0]->operator-=(delta);
  m_halo_zone[1]->operator+=(delta);

}




//---------------------------------------------------------------------------
bool DS_3Dbox:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: isIn(pt)" ) ;

  return ( m_geometric_rigid_body->isIn( pt ) );

}




//---------------------------------------------------------------------------
bool DS_3Dbox:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: isIn(x,y,z)" ) ;

  return ( m_geometric_rigid_body->isIn( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_3Dbox:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: level_set_value(pt)" ) ;

  return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_3Dbox:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: level_set_value(x,y,z)" ) ;

  return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_3Dbox:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: get_distanceTo" ) ;

  return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));

}




//---------------------------------------------------------------------------
geomVector DS_3Dbox:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
geomVector DS_3Dbox:: get_rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: rigid_body_angular_velocity()" ) ;

  return (m_geometric_rigid_body->rigid_body_angular_velocity());

}




//---------------------------------------------------------------------------
std::tuple<double,double> DS_3Dbox:: get_mass_and_density() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: get_mass_and_density()" ) ;

  return ( m_geometric_rigid_body->get_mass_and_density() );

}




//---------------------------------------------------------------------------
double DS_3Dbox:: get_circumscribed_radius( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: get_circumscribed_radius()" ) ;

  return (m_geometric_rigid_body->get_circumscribed_radius());

}




//---------------------------------------------------------------------------
geomVector const* DS_3Dbox:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: get_ptr_to_gravity_centre( )" ) ;

  return (dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre());

}




//---------------------------------------------------------------------------
void DS_3Dbox:: update_RB_position_and_velocity(geomVector const& pos,
                                                    geomVector const& vel,
                                                    geomVector const& ang_vel,
                                   vector<geomVector> const& periodic_directions)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: update_RB_position_and_velocity" ) ;

  return (m_geometric_rigid_body->update_RB_position_and_velocity(pos,vel
                                                                  ,ang_vel
                                                         ,periodic_directions));

}




//---------------------------------------------------------------------------
void DS_3Dbox:: compute_surface_points(  )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: compute_surface_points" ) ;

  doubleVector di(3,0.);

  di(0) = (box_max[0] - box_min[0])/(double)Npoints[0];
  di(1) = (box_max[1] - box_min[1])/(double)Npoints[1];
  di(2) = (box_max[2] - box_min[2])/(double)Npoints[2];

  size_t cntr = 0;

  // x face
  for (size_t i = 0; i < Npoints[1]; i++) {
     for (size_t j = 0; j < Npoints[2]; j++) {
        geomVector point_min (box_min[0]
                            , box_min[1] + (0.5+(double)i)*di(1)
                            , box_min[2] + (0.5+(double)j)*di(2) );
        geomVector point_max (box_max[0]
                            , box_min[1] + (0.5+(double)i)*di(1)
                            , box_min[2] + (0.5+(double)j)*di(2) );
        m_surface_points[cntr]->operator=(point_min);
        m_surface_points[Npoints[1]*Npoints[2] + cntr]->operator=(point_max);

        m_surface_area[cntr]->operator()(0) = di(1)*di(2);
        m_surface_area[Npoints[1]*Npoints[2]+cntr]->operator()(0) = di(1)*di(2);

        geomVector normal(-1., 0., 0.);
        m_surface_normal[cntr]->operator=(normal);
        m_surface_normal[Npoints[1]*Npoints[2] + cntr]->operator=(-1.*normal);

        cntr++;
     }
  }

  cntr = cntr + Npoints[1]*Npoints[2];

  // y face
  for (size_t i = 0; i < Npoints[0]; i++) {
     for (size_t j = 0; j < Npoints[2]; j++) {
        geomVector point_min (box_min[0] + (0.5+(double)i)*di(0)
                            , box_min[1]
                            , box_min[2] + (0.5+(double)j)*di(2) );
        geomVector point_max (box_min[0] + (0.5+(double)i)*di(0)
                            , box_max[1]
                            , box_min[2] + (0.5+(double)j)*di(2) );
        m_surface_points[cntr]->operator=(point_min);
        m_surface_points[Npoints[0]*Npoints[2] + cntr]->operator=(point_max);

        m_surface_area[cntr]->operator()(0) = di(0)*di(2);
        m_surface_area[Npoints[0]*Npoints[2]+cntr]->operator()(0) = di(0)*di(2);

        geomVector normal(0., -1., 0.);
        m_surface_normal[cntr]->operator=(normal);
        m_surface_normal[Npoints[0]*Npoints[2] + cntr]->operator=(-1.*normal);

        cntr++;
     }
  }

  cntr = cntr + Npoints[0]*Npoints[2];

  // z face
  for (size_t i = 0; i < Npoints[0]; i++) {
     for (size_t j = 0; j < Npoints[1]; j++) {
        geomVector point_min (box_min[0] + (0.5+(double)i)*di(0)
                            , box_min[1] + (0.5+(double)j)*di(1)
                            , box_min[2] );
        geomVector point_max (box_min[0] + (0.5+(double)i)*di(0)
                            , box_min[1] + (0.5+(double)j)*di(1)
                            , box_max[2] );
        m_surface_points[cntr]->operator=(point_min);
        m_surface_points[Npoints[0]*Npoints[1] + cntr]->operator=(point_max);

        m_surface_area[cntr]->operator()(0) = di(0)*di(1);
        m_surface_area[Npoints[0]*Npoints[1]+cntr]->operator()(0) = di(0)*di(1);

        geomVector normal(0., 0., -1.);
        m_surface_normal[cntr]->operator=(normal);
        m_surface_normal[Npoints[0]*Npoints[1] + cntr]->operator=(-1.*normal);

        cntr++;
     }
  }

  // Translate and rotate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
     m_geometric_rigid_body->rotate(m_surface_points[i]);
     m_geometric_rigid_body->rotate(m_surface_normal[i]);
  }

}




//---------------------------------------------------------------------------
void DS_3Dbox:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dbox:: compute_number_of_surface_variables" ) ;

  struct FS_3Dbox_Additional_Param const* pagp =
                        dynamic_cast<FS_3Dbox*>(m_geometric_rigid_body)
                           ->get_ptr_FS_3Dbox_Additional_Param();

  if (box_min.empty()) box_min.assign(3,1.e14);
  if (box_max.empty()) box_max.assign(3,-1.e14);

  for (int i = 0; i < (int) pagp->ref_corners.size(); i++) {
     for (int dir = 0; dir < 3; dir++) {
        if (pagp->ref_corners[i](dir) < box_min[dir])
           box_min[dir] = pagp->ref_corners[i](dir);

        if (pagp->ref_corners[i](dir) > box_max[dir])
           box_max[dir] = pagp->ref_corners[i](dir);

     }
  }

  doubleVector delta(3,0.);

  delta(0) = (box_max[0]-box_min[0]);
  delta(1) = (box_max[1]-box_min[1]);
  delta(2) = (box_max[2]-box_min[2]);

  double scale = 1. / sqrt(surface_cell_scale) / dx;

  if (Npoints.empty()) Npoints.reserve(3);

  Npoints[0] = (size_t) round(scale * delta(0) );
  Npoints[1] = (size_t) round(scale * delta(1) );
  Npoints[2] = (size_t) round(scale * delta(2) );

  Ntot = 2 * (Npoints[0]*Npoints[1]
            + Npoints[1]*Npoints[2]
            + Npoints[2]*Npoints[0]);

}
