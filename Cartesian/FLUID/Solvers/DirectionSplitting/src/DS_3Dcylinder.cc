#include <DS_3Dcylinder.hh>
#include <FS_RigidBody.hh>
#include <FS_3Dcylinder.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
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

  geomVector const* pgc = dynamic_cast<FS_3Dcylinder*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  double r_equi = 3.0 *
                  MAC::sqrt(pagp->cylinder_radius * pagp->cylinder_radius
                          + pagp->cylinder_height * pagp->cylinder_height);

  geomVector delta(r_equi, r_equi, r_equi);

  m_halo_zone[0]->operator=(*pgc);
  m_halo_zone[1]->operator=(*pgc);

  m_halo_zone[0]->operator-=(delta);
  m_halo_zone[1]->operator+=(delta);

}




//---------------------------------------------------------------------------
void DS_3Dcylinder:: compute_surface_points(  )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: compute_surface_points" ) ;

  // Pointers to location and additional parameters
  struct FS_3Dcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_3Dcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_3Dcylinder_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_3Dcylinder*>(m_geometric_rigid_body)
                            ->get_ptr_to_gravity_centre();

  // Reference: Becker and Becker, Computational Geometry 45 (2012) 275-283
  // Estimating the number of rings on the hemisphere, assuming Pmin=3
  // and aspect ratio(ar) as 1
  double p = MAC::pi()/ar;
  double k_temp = Ndisk;
  size_t cntr = 0;

  while (k_temp > double(Pmin+2)) {
     k_temp = round(pow(MAC::sqrt(k_temp) - MAC::sqrt(p),2.));
     cntr++;
  }

  size_t Nrings = cntr+1;

  // Summation of total discretized points with
  // increase in number of rings radially
  doubleVector k(Nrings,0.);
  // Assigning the maximum number of discretized
  // points to the last element of the array
  k(Nrings-1) = (double) Ndisk;

  for (int i=(int)Nrings-2; i>=0; --i) {
     k(i) = round(pow(MAC::sqrt(k(i+1)) - MAC::sqrt(p),2.));
     if (i==0) k(0) = (double) Pmin;
  }

  // Radius of the rings in lamber projection plane
  doubleVector Rring(Nrings,0.);

  Rring(Nrings-1) = 1.;

  size_t maxby2 = (size_t) k(Nrings-1);

  // Calculation for all rings except at the pole
  for (int i=(int)Nrings-1; i>0; --i) {
     double Ri = Rring(i);
     Rring(i-1) = MAC::sqrt(k(i-1)/k(i))*Rring(i);
     Rring(i) = (Rring(i) + Rring(i-1))/2.;
     double d_theta = 2.*MAC::pi()/(k(i)-k(i-1));
     // Initialize theta initialize as 1% to avoid
     // point overlap with mesh gridlines
     double theta = 0.01*d_theta;
     for (int j=(int)k(i-1); j<k(i); j++) {
        theta = theta + d_theta;

        geomVector point (pagp->cylinder_radius*Rring(i)*MAC::cos(theta)
                        , pagp->cylinder_height*1./2.
                        , pagp->cylinder_radius*Rring(i)*MAC::sin(theta) );
        geomVector point_mirror (pagp->cylinder_radius*Rring(i)*MAC::cos(theta)
                              , pagp->cylinder_height*-1./2.
                              , pagp->cylinder_radius*Rring(i)*MAC::sin(theta));
        // For top disk
        m_surface_points[j]->operator=(point);
        m_surface_area[j]->operator()(0) = 0.5*pagp->cylinder_radius
                                        *pagp->cylinder_radius
                                        *d_theta*(pow(Ri,2)-pow(Rring(i-1),2));

     	  // For bottom disk
        m_surface_points[maxby2+j]->operator=(point_mirror);
        m_surface_area[maxby2+j] = m_surface_area[j];

        // Create surface normal vectors
        geomVector normal(0., pagp->cylinder_height*1./2., 0.);
        m_surface_normal[j]->operator=(normal);
        m_surface_normal[maxby2+j]->operator=(-1.*normal);
     }
  }

  // Calculation at the ring on pole (i=0)
  double Ri = Rring(0);
  Rring(0) = Rring(0)/2.;
  double d_theta = 2.*MAC::pi()/(k(0));
  // Initialize theta as 1% of the d_theta to avoid
  // point overlap with mesh gridlines
  double theta = 0.01*d_theta;
  if (k(0)>1) {
     for (int j=0; j < k(0); j++) {
        theta = theta + d_theta;

        geomVector point (pagp->cylinder_radius*Rring(0)*MAC::cos(theta)
                        , pagp->cylinder_height*1./2.
                        , pagp->cylinder_radius*Rring(0)*MAC::sin(theta) );
        geomVector point_mirror (pagp->cylinder_radius*Rring(0)*MAC::cos(theta)
                              , pagp->cylinder_height*-1./2.
                              , pagp->cylinder_radius*Rring(0)*MAC::sin(theta));
        // For top disk
        m_surface_points[j]->operator=(point);
        m_surface_area[j]->operator()(0) = 0.5*pagp->cylinder_radius
                                          *pagp->cylinder_radius
                                          *d_theta*pow(Ri,2);
        // For bottom disk
        m_surface_points[maxby2+j]->operator=(point_mirror);
        m_surface_area[maxby2+j] = m_surface_area[j];
        // Create surface normal vectors
        geomVector normal(0., pagp->cylinder_height*1./2., 0.);
        m_surface_normal[j]->operator=(normal);
        m_surface_normal[maxby2+j]->operator=(-1.*normal);
     }
  } else {
     geomVector normal(0., pagp->cylinder_height*1./2., 0.);
     // For top disk
     m_surface_points[0]->operator=(normal);
     m_surface_area[0]->operator()(0) = 0.5*pagp->cylinder_radius
                                       *pagp->cylinder_radius
                                       *d_theta*pow(Ri,2.);
     // For bottom disk
     m_surface_points[maxby2]->operator=(-1.*normal);
     m_surface_area[maxby2] = m_surface_area[0];
     // Create surface normal vectors
     m_surface_normal[0]->operator=(normal);
     m_surface_normal[maxby2]->operator=(-1.*normal);
  }

  // Generating one ring of points on cylindrical surface
  // Can be used to calculate stress on whole surface by
  // a constant shift of points

  int pts_1_ring = (int)(k(Nrings-1) - k(Nrings-2));
  int cyl_rings = (int)round((pagp->cylinder_height/pagp->cylinder_radius)
                               /(1-MAC::sqrt(k(Nrings-2)/k(Nrings-1))));
  double cell_area = 2.*pagp->cylinder_radius*pagp->cylinder_height*MAC::pi()
                   / ((double) pts_1_ring*cyl_rings);

  d_theta = 2.*MAC::pi()/pts_1_ring;
  for (int j=0; j<cyl_rings; j++) {
     theta = 0.01*d_theta;
     for (int ij=0; ij<pts_1_ring; ij++) {
        theta = theta + d_theta;
        int n = 2*(int)k(Nrings-1) + j*pts_1_ring + ij;
        geomVector point (pagp->cylinder_radius*MAC::cos(theta)
                        , pagp->cylinder_height
                              *(-1./2. + (j+0.5)/(double(cyl_rings)))
                        , pagp->cylinder_radius*MAC::sin(theta) );
        geomVector normal (pagp->cylinder_radius*MAC::cos(theta)
                         , 0.
                         , pagp->cylinder_radius*MAC::sin(theta));
        m_surface_points[n]->operator=(point);
        m_surface_area[n]->operator()(0) = cell_area;
        m_surface_normal[n]->operator=(normal);
     }
  }

  // Translate and rotate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
     m_geometric_rigid_body->rotate(m_surface_points[i]);
     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_surface_points[i]->translate(*pgc);
  }

}




//---------------------------------------------------------------------------
void DS_3Dcylinder:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dcylinder:: compute_number_of_surface_variables" ) ;

  struct FS_3Dcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_3Dcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_3Dcylinder_Additional_Param();

  Ndisk = round((1./surface_cell_scale)
              * ( MAC::pi()*pagp->cylinder_radius*pagp->cylinder_radius )
              / ( dx*dx ));

  double Npm1 = round(pow(MAC::sqrt(Ndisk) - MAC::sqrt(MAC::pi()/ar),2.));
  double dh = 1. - MAC::sqrt(Npm1/Ndisk);
  double Nr = round((pagp->cylinder_height/pagp->cylinder_radius)/dh);
  Ntot = (size_t) (2*Ndisk + Nr*(Ndisk - Npm1));

}
