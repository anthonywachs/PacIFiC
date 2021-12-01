#include <DS_RigidBody.hh>
#include <FS_RigidBody.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
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

  m_halo_zone.reserve(2);
  m_halo_zone.push_back(new geomVector(3));
  m_halo_zone.push_back(new geomVector(3));

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
double DS_RigidBody:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: level_set_value(pt)" ) ;

  return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_RigidBody:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: level_set_value(x,y,z)" ) ;

  return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

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




//---------------------------------------------------------------------------
double DS_RigidBody:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_distanceTo" ) ;

  return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));

}




//---------------------------------------------------------------------------
geomVector DS_RigidBody:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
void DS_RigidBody:: initialize_surface_variables( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: initialize_surface_variables" ) ;

  m_surface_points.reserve( Ntot );
  m_surface_area.reserve( Ntot );
  m_surface_normal.reserve( Ntot );
  m_surface_Pforce.reserve( Ntot );
  m_surface_Vforce.reserve( Ntot );

  geomVector vvv(3);

   for (size_t i = 0; i < Ntot; ++i) {
      m_surface_points.push_back( new geomVector(3) );
      m_surface_area.push_back( 0. );
      m_surface_normal.push_back( new geomVector(3) );
      m_surface_Pforce.push_back( vvv );
      m_surface_Vforce.push_back( vvv );
   }

}




// //---------------------------------------------------------------------------
// vector<geomVector*> DS_RigidBody:: get_rigid_body_surface_points( ) const
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL( "DS_RigidBody:: get_rigid_body_surface_points" ) ;
//
//   return (&m_surface_points);
//
// }
//
//
//
//
// //---------------------------------------------------------------------------
// vector<geomVector*> DS_RigidBody:: get_rigid_body_surface_normals( ) const
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL( "DS_RigidBody:: get_rigid_body_surface_normals" ) ;
//
//   return (&m_surface_normal);
//
// }
//
//
//
//
// //---------------------------------------------------------------------------
// vector<double*> DS_RigidBody:: get_rigid_body_surface_areas( ) const
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL( "DS_RigidBody:: get_rigid_body_surface_areas" ) ;
//
//   return (&m_surface_area);
//
// }




//---------------------------------------------------------------------------
vector<geomVector*> DS_RigidBody:: get_rigid_body_haloZone( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody:: get_rigid_body_haloZone" ) ;

  return (m_halo_zone);

}




//---------------------------------------------------------------------------
double DS_RigidBody:: Bilinear_interpolation ( FV_DiscreteField const* FF
                                             , size_t const& comp
                                             , geomVector const& point
                                             , geomVector const& face_vec
                                             , size_t_vector const* void_fraction
                                             , size_t_array2D* intersect_vector
                                             , doubleArray2D* intersect_distance
                                             , doubleArray2D* intersect_fieldValue
                                             , size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_RigidBody:: Bilinear_interpolation" ) ;

   size_t nb_comps = FF->nb_components() ;
   size_t field = (nb_comps == 1) ? 0 : 1;
   size_t i0_temp;
   size_t_vector i0(3,0);
   size_t dim = FF->primary_grid()->nb_space_dimensions() ;
   double dh = FF->primary_grid()->get_smallest_grid_size();

   // Finding the grid indexes next to ghost points
   bool found = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,0)
                                     , point(0), i0_temp);
   if (found == 1) i0(0) = i0_temp;

   found = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,1)
                                     , point(1), i0_temp);
   if (found == 1) i0(1) = i0_temp;

   if (dim == 3) {
      found = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,2)
                                     , point(2), i0_temp);
      if (found == 1) i0(2) = i0_temp;
   }

   size_t_array2D p(2,2,0);
   // arrays of vertex indexes of the cube/square
   size_t_array2D ix(2,2,0);
   size_t_array2D iy(2,2,0);
   size_t_array2D iz(2,2,0);
   // Min and max of the cell containing ghost point
   doubleArray2D extents(dim,2,0);
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
         f(i,j) = FF->DOF_value( ix(i,j), iy(i,j), iz(i,j), comp, level );
      }
   }

   // Min and max coordinate in the grid cell
   for (size_t dir = 0; dir < dim; dir++) {
      extents(dir,0) = FF->get_DOF_coordinate(i0(dir),comp,dir) ;
      extents(dir,1) = FF->get_DOF_coordinate(i0(dir)+face_vec(dir),comp,dir) ;
   }

   // Contribution from left and right wall
   for (size_t i = 0; i < 2; i++) {     // 0 --> left; 1 --> right
      size_t col_top = 2*dir2 + 1;
      size_t col_bot = 2*dir2 + 0;
      if ((field == 0) || ((void_fraction->operator()(p(i,0)) == 0)
                        && (void_fraction->operator()(p(i,1)) == 0))) {
         fwall(0,i) = ((extents(dir2,1) - point(dir2))*f(i,0)
                     + (point(dir2) - extents(dir2,0))*f(i,1))
                     / (extents(dir2,1)-extents(dir2,0));
         del_wall(0,i) = MAC::abs(extents(dir1,i) - point(dir1));
      // if bottom vertex is in fluid domain
      } else if ((void_fraction->operator()(p(i,0)) == 0)
              && (intersect_vector->operator()(p(i,0),col_top) == 1)) {
         double yint = intersect_distance->operator()(p(i,0),col_top);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (yint >= (point(dir2)-extents(dir2,0))) {
            fwall(0,i) = ((extents(dir2,0) + yint - point(dir2))*f(i,0)
                        + (point(dir2) - extents(dir2,0))
                              *intersect_fieldValue->operator()(p(i,0),col_top))
                        / yint;
            del_wall(0,i) = MAC::abs(extents(dir1,i) - point(dir1));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction->operator()(p(i,1)) - 1;
            geomVector rayDir(3);
            rayDir(dir1) = (i == 0) ? -1 : 1 ;
            del_wall(0,i) = m_geometric_rigid_body
                                             ->distanceTo( point, rayDir, dh );
            geomVector surface_point = point + del_wall(0,i)*rayDir;
            geomVector net_vel = get_rigid_body_velocity(surface_point);
            fwall(0,i) = net_vel(comp);
         }
      // if top vertex is in fluid domain
      } else if ((void_fraction->operator()(p(i,1)) == 0)
              && (intersect_vector->operator()(p(i,1),col_bot) == 1)) {
         double yint = intersect_distance->operator()(p(i,1),col_bot);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (yint >= (extents(dir2,1)-point(dir2))) {
            fwall(0,i) = ((point(dir2) + yint - extents(dir2,1))*f(i,1)
                        + (extents(dir2,1) - point(dir2))
                              *intersect_fieldValue->operator()(p(i,1),col_bot))
                        / yint;
            del_wall(0,i) = MAC::abs(extents(dir1,i) - point(dir1));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction->operator()(p(i,0)) - 1;
            geomVector rayDir(3);
            rayDir(dir1) = (i == 0) ? -1 : 1 ;
            del_wall(0,i) = m_geometric_rigid_body
                                             ->distanceTo( point, rayDir, dh );
            geomVector surface_point = point + del_wall(0,i)*rayDir;
            geomVector net_vel = get_rigid_body_velocity(surface_point);
            fwall(0,i) = net_vel(comp);
         }
      // if both vertex's are in solid domain
      } else if ((void_fraction->operator()(p(i,0)) != 0)
              && (void_fraction->operator()(p(i,1)) != 0)) {

         size_t id = void_fraction->operator()(p(i,0)) - 1;
         geomVector rayDir(3);
         rayDir(dir1) = (i == 0) ? -1 : 1 ;
         del_wall(0,i) = m_geometric_rigid_body
                                          ->distanceTo( point, rayDir, dh );
         geomVector surface_point = point + del_wall(0,i)*rayDir;
         geomVector net_vel = get_rigid_body_velocity(surface_point);
         fwall(0,i) = net_vel(comp);
      }
   }

   // Contribution from top and bottom wall
   for (size_t j = 0; j < 2; j++) {         // 0 --> bottom; 1 --> top
      size_t col_right = 2*dir1 + 1;
      size_t col_left = 2*dir1 + 0;
      if ((field == 0) || ((void_fraction->operator()(p(0,j)) == 0)
                        && (void_fraction->operator()(p(1,j)) == 0))) {
         fwall(1,j) = ((extents(dir1,1) - point(dir1)) * f(0,j)
                     + (point(dir1) - extents(dir1,0)) * f(1,j))
                     / (extents(dir1,1) - extents(dir1,0));
         del_wall(1,j) = MAC::abs(extents(dir2,j) - point(dir2));
      // if left vertex is in fluid domain
      } else if ((void_fraction->operator()(p(0,j)) == 0)
              && (intersect_vector->operator()(p(0,j),col_right) == 1)) {
         double xint = intersect_distance->operator()(p(0,j),col_right);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (xint >= (point(dir1)-extents(dir1,0))) {
            fwall(1,j) = ((extents(dir1,0) + xint - point(dir1)) * f(0,j)
                        + (point(dir1) - extents(dir1,0))
                           * intersect_fieldValue->operator()(p(0,j),col_right))
                        / xint;
            del_wall(1,j) = MAC::abs(extents(dir2,j) - point(dir2));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction->operator()(p(1,j)) - 1;
            geomVector rayDir(3);
            rayDir(dir2) = (j == 0) ? -1 : 1 ;
            del_wall(1,j) = m_geometric_rigid_body
                                             ->distanceTo( point, rayDir, dh );
            geomVector surface_point = point + del_wall(1,j)*rayDir;
            geomVector net_vel = get_rigid_body_velocity(surface_point);
            fwall(1,j) = net_vel(comp);
         }
      // if right vertex is in fluid domain
      } else if ((void_fraction->operator()(p(1,j)) == 0)
              && (intersect_vector->operator()(p(1,j),col_left) == 1)) {

         double xint = intersect_distance->operator()(p(1,j),col_left);
         // Condition where intersection distance is more than
         // ghost point distance, it means that the ghost point
         // can be projected on the cell wall
         if (xint >= (extents(dir1,1)-point(dir1))) {
            fwall(1,j) = ((point(dir1) + xint - extents(dir1,1)) * f(1,j)
                        + (extents(dir1,1) - point(dir1))
                           * intersect_fieldValue->operator()(p(1,j),col_left))
                        / xint;
            del_wall(1,j) = MAC::abs(extents(dir2,j) - point(dir2));
         // Ghost point cannot be projected on the cell wall
         // as the solid surface come first
         } else {
            size_t id = void_fraction->operator()(p(0,j)) - 1;
            geomVector rayDir(3);
            rayDir(dir2) = (j == 0) ? -1 : 1 ;
            del_wall(1,j) = m_geometric_rigid_body
                                             ->distanceTo( point, rayDir, dh );
            geomVector surface_point = point + del_wall(1,j)*rayDir;
            geomVector net_vel = get_rigid_body_velocity(surface_point);
            fwall(1,j) = net_vel(comp);
         }
      // if both vertex's are in solid domain
      } else if ((void_fraction->operator()(p(0,j)) != 0)
              && (void_fraction->operator()(p(1,j)) != 0)) {
         size_t id = void_fraction->operator()(p(0,j)) - 1 ;
         geomVector rayDir(3);
         rayDir(dir2) = (j == 0) ? -1 : 1 ;
         del_wall(1,j) = m_geometric_rigid_body
                                          ->distanceTo( point, rayDir, dh );
         geomVector surface_point = point + del_wall(1,j)*rayDir;
         geomVector net_vel = get_rigid_body_velocity(surface_point);
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
