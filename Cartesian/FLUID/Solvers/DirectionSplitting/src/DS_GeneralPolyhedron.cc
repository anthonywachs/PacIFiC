#include <DS_GeneralPolyhedron.hh>
#include <FS_RigidBody.hh>
#include <FS_GeneralPolyhedron.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_GeneralPolyhedron:: DS_GeneralPolyhedron()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_GeneralPolyhedron:: DS_GeneralPolyhedron" ) ;

}




//---------------------------------------------------------------------------
DS_GeneralPolyhedron:: DS_GeneralPolyhedron( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_GeneralPolyhedron:: ~DS_GeneralPolyhedron()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: ~DS_GeneralPolyhedron" ) ;

}




//---------------------------------------------------------------------------
void DS_GeneralPolyhedron:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_GeneralPolyhedron:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
void DS_GeneralPolyhedron:: compute_rigid_body_halozone( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: compute_rigid_body_halozone" ) ;

  struct FS_GeneralPolyhedron_Additional_Param const* pagp =
   dynamic_cast<FS_GeneralPolyhedron*>(m_geometric_rigid_body)
      ->get_ptr_FS_GeneralPolyhedron_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_GeneralPolyhedron*>(m_geometric_rigid_body)
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
bool DS_GeneralPolyhedron:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: isIn(pt)" ) ;

  return ( m_geometric_rigid_body->isIn( pt ) );

}




//---------------------------------------------------------------------------
bool DS_GeneralPolyhedron:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: isIn(x,y,z)" ) ;

  return ( m_geometric_rigid_body->isIn( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_GeneralPolyhedron:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: level_set_value(pt)" ) ;

  return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_GeneralPolyhedron:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: level_set_value(x,y,z)" ) ;

  return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_GeneralPolyhedron:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: get_distanceTo" ) ;

  return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));

}




//---------------------------------------------------------------------------
geomVector DS_GeneralPolyhedron:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
geomVector DS_GeneralPolyhedron:: get_rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: rigid_body_angular_velocity()" ) ;

  return (m_geometric_rigid_body->rigid_body_angular_velocity());

}




//---------------------------------------------------------------------------
std::tuple<double,double,double> DS_GeneralPolyhedron:: get_mass_and_density_and_moi() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: get_mass_and_density()" ) ;

  return ( m_geometric_rigid_body->get_mass_and_density_and_moi() );

}




//---------------------------------------------------------------------------
double DS_GeneralPolyhedron:: get_circumscribed_radius( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: get_circumscribed_radius()" ) ;

  return (m_geometric_rigid_body->get_circumscribed_radius());

}




//---------------------------------------------------------------------------
geomVector const* DS_GeneralPolyhedron:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: get_ptr_to_gravity_centre( )" ) ;

  return (dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre());

}




//---------------------------------------------------------------------------
void DS_GeneralPolyhedron:: update_RB_position_and_velocity(geomVector const& pos,
                                                    geomVector const& vel,
                                                    geomVector const& ang_vel,
                                vector<geomVector> const& periodic_directions,
                                   double const& time_step)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: update_RB_position_and_velocity" ) ;

  return (m_geometric_rigid_body->update_RB_position_and_velocity(pos,vel
                                                                  ,ang_vel
                                                         ,periodic_directions
                                                         , time_step));

}




//---------------------------------------------------------------------------
void DS_GeneralPolyhedron:: update_additional_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: update_additional_parameters" ) ;

  m_geometric_rigid_body->update_additional_parameters();

}




//---------------------------------------------------------------------------
void DS_GeneralPolyhedron:: compute_surface_points(  )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: compute_surface_points" ) ;

  size_t cntr = 0;

  struct FS_GeneralPolyhedron_Additional_Param const* pagp =
   dynamic_cast<FS_GeneralPolyhedron*>(m_geometric_rigid_body)
      ->get_ptr_FS_GeneralPolyhedron_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_GeneralPolyhedron*>(m_geometric_rigid_body)
                          ->get_ptr_to_gravity_centre();

  geomVector gc(3);
  gc(0) = pgc->operator()(0);
  gc(1) = pgc->operator()(1);
  gc(2) = pgc->operator()(2);

  size_t nfaces = pagp->faceCen.size();

  for (size_t iface = 0; iface < nfaces; iface++) {
     size_t n_vertex = pagp->facesVec[iface].size();
     if (n_vertex == 3) {
        // Vertex ID
        size_t v1 = pagp->facesVec[iface][0];
        size_t v2 = pagp->facesVec[iface][1];
        size_t v3 = pagp->facesVec[iface][2];

        geomVector p1 = pagp->corners[v1];
        geomVector p2 = pagp->corners[v2];
        geomVector p3 = pagp->corners[v3];

        geomVector vec1 = (1./(dis_level+1.)) * (p2 - p1);
        geomVector vec2 = (1./(dis_level+1.)) * (p3 - p1);
        geomVector cross = (vec1^vec2);
        double area = 0.5*cross.calcNorm();

        // Test the direction of normal vector to be away from RB center
        geomVector delta = pagp->ref_corners[v1];

        if ((delta(0)*cross(0)
           + delta(1)*cross(1)
           + delta(2)*cross(2)) < 0.) {
           cross = -1.*cross;
        }
        cross *= (1./cross.calcNorm());
        //-----------------------------------------------------------------

        for (size_t i = 0; i <= dis_level; i++) {
           for (size_t j = 0; j <= dis_level - i; j++) {
             // four points of parallelogram
             geomVector par1 = p1 + i * vec1 + j * vec2;
             geomVector par2 = par1 + vec1;
             geomVector par3 = par1 + vec2;
             geomVector par4 = par1 + vec1 + vec2;
             m_surface_points[cntr]->operator=((1./3.) * (par1 + par2 + par3));
             m_surface_normal[cntr]->operator=(cross);
             m_surface_area[cntr]->operator()(0) = area;
             cntr++;
             if (j < dis_level - i) {
                m_surface_points[cntr]->operator=((1./3.) * (par2 + par3 + par4));
                m_surface_normal[cntr]->operator=(cross);
                m_surface_area[cntr]->operator()(0) = area;
                cntr++;
             }
           }
        }
     } else {
        for (size_t idxPts = 0; idxPts < n_vertex; ++idxPts ) {
           // Vertex index
           size_t i0 = idxPts;
           size_t i1 = (idxPts + 1 == n_vertex) ? 0 : idxPts + 1;
           // Vertex ID
           size_t v1 = pagp->facesVec[iface][i0];
           size_t v2 = pagp->facesVec[iface][i1];

           geomVector p1 = pagp->faceCen[iface];
           geomVector p2 = pagp->corners[v1];
           geomVector p3 = pagp->corners[v2];

           geomVector vec1 = (1./(dis_level+1.)) * (p2 - p1);
           geomVector vec2 = (1./(dis_level+1.)) * (p3 - p1);
           geomVector cross = (vec1^vec2);
           double area = 0.5*cross.calcNorm();

           // Test the direction of normal vector to be away from RB center
           geomVector delta = pagp->ref_corners[v1];

           if ((delta(0)*cross(0)
              + delta(1)*cross(1)
              + delta(2)*cross(2)) < 0.) {
              cross = -1.*cross;
           }
           cross *= (1./cross.calcNorm());
           //-----------------------------------------------------------------

           for (size_t i = 0; i <= dis_level; i++) {
              for (size_t j = 0; j <= dis_level - i; j++) {
                // four points of parallelogram
                geomVector par1 = p1 + i * vec1 + j * vec2;
                geomVector par2 = par1 + vec1;
                geomVector par3 = par1 + vec2;
                geomVector par4 = par1 + vec1 + vec2;
                m_surface_points[cntr]->operator=((1./3.) * (par1 + par2 + par3));
                m_surface_normal[cntr]->operator=(cross);
                m_surface_area[cntr]->operator()(0) = area;
                cntr++;
                if (j < dis_level - i) {
                   m_surface_points[cntr]->operator=((1./3.) * (par2 + par3 + par4));
                   m_surface_normal[cntr]->operator=(cross);
                   m_surface_area[cntr]->operator()(0) = area;
                   cntr++;
                }
             }
           }
        }
     }
  }
}




//---------------------------------------------------------------------------
void DS_GeneralPolyhedron:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_GeneralPolyhedron:: compute_number_of_surface_variables" ) ;

  struct FS_GeneralPolyhedron_Additional_Param const* pagp =
                        dynamic_cast<FS_GeneralPolyhedron*>(m_geometric_rigid_body)
                           ->get_ptr_FS_GeneralPolyhedron_Additional_Param();

  double size = get_circumscribed_radius();

  double scale = 1. / sqrt(surface_cell_scale) / dx;

  dis_level = MAC::floor(size*scale);

  size_t nfaces = pagp->faceCen.size();

  size_t faceEdges = pagp->facesVec[0].size();

  Ntot = (faceEdges == 3) ?
         (size_t) nfaces * (dis_level + 1.) * (dis_level + 1.)
       : (size_t) nfaces * (dis_level + 1.) * (dis_level + 1.) * faceEdges;

}
