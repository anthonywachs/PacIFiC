#include <DS_Sphere.hh>
#include <FS_RigidBody.hh>
#include <FS_Sphere.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
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
void DS_Sphere:: compute_rigid_body_halozone( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: compute_rigid_body_halozone" ) ;

  struct FS_Sphere_Additional_Param const* pagp =
   dynamic_cast<FS_Sphere*>(m_geometric_rigid_body)
      ->get_ptr_FS_Sphere_Additional_Param();

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
void DS_Sphere:: compute_surface_points( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: compute_surface_points" ) ;

  // Pointers to location and additional parameters
  struct FS_Sphere_Additional_Param const* pagp =
   dynamic_cast<FS_Sphere*>(m_geometric_rigid_body)
      ->get_ptr_FS_Sphere_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                            ->get_ptr_to_gravity_centre();

  size_t kmax = m_surface_area.size()/2;

  // Reference: Becker and Becker, Computational Geometry 45 (2012) 275-283
  double eta_temp = MAC::pi()/2.;
  double k_temp = (double) kmax;
  double Ro_temp = MAC::sqrt(2);
  double Rn_temp = MAC::sqrt(2);
  size_t cntr = 0;

  // Estimating the number of rings on the hemisphere, assuming Pmin=3
  // and aspect ratio(ar) as 1
  size_t Pmin = 3;
  double ar = 1.;
  while (k_temp > double(Pmin+2)) {
     eta_temp = eta_temp - (2./ar)
                         * MAC::sqrt(MAC::pi()/k_temp)
                         * MAC::sin(eta_temp/2.);
     Rn_temp = 2.*MAC::sin(eta_temp/2.);
     k_temp = round(k_temp*MAC::pow(Rn_temp/Ro_temp,2.));
     Ro_temp = Rn_temp;
     cntr++;
  }

  size_t Nrings = cntr+1;

  // Summation of total discretized points
  // with increase in number of rings radially
  doubleVector k(Nrings,0.);
  // Zenithal angle for the sphere
  doubleVector eta(Nrings,0.);
  // Radius of the rings in lamber projection plane
  doubleVector Rring(Nrings,0.);

  // Assigning the maximum number of discretized
  // points to the last element of the array
  k(Nrings-1) = (double) kmax;
  // Zenithal angle for the last must be pi/2.
  eta(Nrings-1) = MAC::pi()/2.;
  // Radius of last ring in lamber projection plane
  Rring(Nrings-1) = MAC::sqrt(2.);

  for (int i=int(Nrings)-2; i>=0; --i) {
     eta(i) = eta(i+1) - (2./ar)
                       * MAC::sqrt(MAC::pi()/k(i+1))
                       * MAC::sin(eta(i+1)/2.);
     Rring(i) = 2.*MAC::sin(eta(i)/2.);
     k(i) = round(k(i+1)*MAC::pow(Rring(i)/Rring(i+1),2.));
     if (i==0) k(0) = (double) Pmin;
  }

  size_t maxby2 = (size_t) k(Nrings-1);
  // Calculation for all rings except at the pole
  for (int i=(int)Nrings-1; i>0; --i) {
     double Ri = Rring(i);
     Rring(i) = (Rring(i) + Rring(i-1))/2.;
     eta(i) = (eta(i) + eta(i-1))/2.;
     double d_theta = 2.*MAC::pi()/(k(i)-k(i-1));
     // Initialize theta as 1% of the d_theta to avoid
     // point overlap with mesh gridlines
     double theta = 0.01*d_theta;

     for (int j=(int)k(i-1); j<k(i); j++) {
        theta = theta + d_theta;

        geomVector point( pagp->radius*MAC::cos(eta(i))
                        , pagp->radius*MAC::cos(theta)*MAC::sin(eta(i))
                        , pagp->radius*MAC::sin(theta)*MAC::sin(eta(i)) );

        geomVector mirror_point(-pagp->radius*MAC::cos(eta(i))
                               , pagp->radius*MAC::cos(theta)*MAC::sin(eta(i))
                               , pagp->radius*MAC::sin(theta)*MAC::sin(eta(i)));

        m_surface_points[j] = point;
        m_surface_area[j] = pagp->radius*pagp->radius*
                            (0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
        // For second half of sphere
        m_surface_points[maxby2+j] = mirror_point;
        m_surface_area[maxby2+j] = m_surface_area[j];

  	     // Create surface normal vectors
        m_surface_normal[j] = m_surface_points[j];
  	     m_surface_normal[maxby2+j] = m_surface_points[maxby2+j];
     }
  }

  // Calculation at the ring on pole (i=0)
  double Ri = Rring(0);
  Rring(0) = Rring(0)/2.;
  eta(0) = eta(0)/2.;
  double d_theta = 2.*MAC::pi()/(k(0));

  // Initialize theta as 1% of the d_theta to avoid
  // point overlap with mesh gridlines
  double theta = 0.01*d_theta;
  if (k(0)>1) {
     for (int j=0; j < k(0); j++) {
        theta = theta + d_theta;

        geomVector point( pagp->radius*MAC::cos(eta(0))
                        , pagp->radius*MAC::cos(theta)*MAC::sin(eta(0))
                        , pagp->radius*MAC::sin(theta)*MAC::sin(eta(0)) );

        geomVector mirror_point(-pagp->radius*MAC::cos(eta(0))
                               , pagp->radius*MAC::cos(theta)*MAC::sin(eta(0))
                               , pagp->radius*MAC::sin(theta)*MAC::sin(eta(0)));

        m_surface_points[j] = point;
        m_surface_area[j] = pagp->radius*pagp->radius*0.5*d_theta*pow(Ri,2.);

        // For second half of sphere
        m_surface_points[maxby2+j] = mirror_point;
        m_surface_area[maxby2+j] = m_surface_area[j];

        // Create surface normal vectors
        m_surface_normal[j] = m_surface_points[j];
        m_surface_normal[maxby2+j] = m_surface_points[maxby2+j];
     }
  } else {
     geomVector point( pagp->radius*1., 0., 0.);
     geomVector mirror_point( -pagp->radius*1., 0., 0.);

     m_surface_points[0] = point;
     m_surface_area[0] = pagp->radius*pagp->radius*0.5*d_theta*pow(Ri,2.);

     //  For second half of sphere
     m_surface_points[maxby2] = mirror_point;
     m_surface_area[maxby2] = m_surface_area[0];

     // Create surface normal vectors
     m_surface_normal[0] = m_surface_points[0];
     m_surface_normal[maxby2] = m_surface_points[maxby2];
  }

  // Translate and rotate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
     m_geometric_rigid_body->rotate(m_surface_points[i]);
     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_surface_points[i].translate(*pgc);
  }

}




//---------------------------------------------------------------------------
void DS_Sphere:: initialize_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: initialize_variable_for_each_rigidBody" ) ;

  struct FS_Sphere_Additional_Param const* pagp =
   dynamic_cast<FS_Sphere*>(m_geometric_rigid_body)
      ->get_ptr_FS_Sphere_Additional_Param();

  size_t Ntot = (size_t) (round(1./surface_cell_scale)
               *(4.*MAC::pi()*pagp->radius*pagp->radius)
               /(dx*dx));

  // Getting the nearest even number
  Ntot = (size_t) (round((double)Ntot * 0.5) * 2.);

  m_surface_points.reserve( Ntot );
  m_surface_area.reserve( Ntot );
  m_surface_normal.reserve( Ntot );

  geomVector vvv(3);

   for (size_t i = 0; i < Ntot; ++i) {
      m_surface_points.push_back( vvv );
      m_surface_area.push_back( 0. );
      m_surface_normal.push_back( vvv );
   }

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

  // Determine the surface point coordinates
  MAC::out() << "DS_Sphere:: compute_hydro_force_torque - "
  	"Determining the set of surface points requires programming" << endl;

  // Compute the surface integrals and stores the force & torque value in
  // the corresponding geometric rigid body
  compute_surface_integrals_hydro_force_torque( PP, UU );

}
