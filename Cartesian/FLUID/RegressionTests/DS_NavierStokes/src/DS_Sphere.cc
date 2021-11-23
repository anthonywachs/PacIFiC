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

  geomVector const* pgs = dynamic_cast<FS_Sphere*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  double radius_equivalent = 3.0*pagp->radius;

  m_halo_zone[0]->operator()(0) = pgs->operator()(0) - radius_equivalent;
  m_halo_zone[0]->operator()(1) = pgs->operator()(1) - radius_equivalent;
  m_halo_zone[0]->operator()(2) = pgs->operator()(2) - radius_equivalent;

  m_halo_zone[1]->operator()(0) = pgs->operator()(0) + radius_equivalent;
  m_halo_zone[1]->operator()(1) = pgs->operator()(1) + radius_equivalent;
  m_halo_zone[1]->operator()(2) = pgs->operator()(2) + radius_equivalent;

}




//---------------------------------------------------------------------------
void DS_Sphere:: compute_surface_points( size_t const& Np )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Sphere:: compute_surface_points" ) ;

  size_t kmax = Np;

  // Reference: Becker and Becker, Computational Geometry 45 (2012) 275-283
  double eta_temp = MAC::pi()/2.;
  double k_temp = (double) kmax;
  double Ro_temp = MAC::sqrt(2);
  double Rn_temp = MAC::sqrt(2);
  size_t cntr = 0;
  size_t Pmin = 3;
  double ar = 1.;

  // Estimating the number of rings on the hemisphere, assuming Pmin=3
  // and aspect ratio(ar) as 1
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
