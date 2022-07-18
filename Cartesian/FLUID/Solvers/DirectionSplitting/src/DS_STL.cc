#include <DS_STL.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_STL:: DS_STL(FV_Mesh const* MESH, istream& STL_input)
//---------------------------------------------------------------------------
: m_MESH( MESH )
{
  MAC_LABEL( "DS_STL:: DS_STL" ) ;

  std::getline(STL_input,filename);
  STL_input >> Nopx >> Nopy >> Nopz;
  STL_input >> cenh;
  STL_input >> invertSTL;

  // Dynamic allocation of the variables of the STL

  tridx_xz = new doubleArray2D(1,1,0.);
  tridx_xz->re_initialize(Nopx,Nopz);

  tridx_xy = new doubleArray2D(1,1,0.);
  tridx_xy->re_initialize(Nopx,Nopy);

  tridx_yz = new doubleArray2D(1,1,0.);
  tridx_yz->re_initialize(Nopy,Nopz);

  vz = vector<vector<tuple<double,double,double>>>(Nopz,v);
  vy = vector<vector<tuple<double,double,double>>>(Nopy,v);

  ttrbox_xz = vector<vector<vector<tuple<double,double,double>>>>(Nopx,vz);
  ttrbox_xy = vector<vector<vector<tuple<double,double,double>>>>(Nopx,vy);
  ttrbox_yz = vector<vector<vector<tuple<double,double,double>>>>(Nopy,vz);

  readSTL();

  // Declating the variable to store gravity centre;
  ptr_gravity_centre = new geomVector(3);

  std::cout << "Construction of STL object completed" << endl;
}




//---------------------------------------------------------------------------
DS_STL:: ~DS_STL()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: ~DS_STL" ) ;

  if ( !m_surface_points.empty() ) m_surface_points.clear();

}




//---------------------------------------------------------------------------
void DS_STL:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: update" ) ;

}





//---------------------------------------------------------------------------
void DS_STL:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  std::cout << "- Displaying Matrix of filtered triangles xz" << endl;
  for (int i=0; i<Nopx; i++)
  {
     for (int k=0; k<Nopz; k++)
     {
        std::cout << tridx_xz->operator()(i,k)/3 << "\t";
     }
     std::cout << endl;
  }
  std::cout << "- Displaying Matrix of filtered triangles xy" << endl;
  for (int i=0; i<Nopx; i++)
  {
     for (int k=0; k<Nopy; k++)
     {
        std::cout << tridx_xy->operator()(i,k)/3 << "\t";
     }
     std::cout << endl;
  }
  std::cout << "- Displaying Matrix of filtered triangles yz" << endl;
  for (int i=0; i<Nopy; i++)
  {
     for (int k=0; k<Nopz; k++)
     {
        std::cout << tridx_yz->operator()(i,k)/3 << "\t";
     }
     std::cout << endl;
  }

}




//---------------------------------------------------------------------------
void DS_STL:: compute_rigid_body_halozone( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: compute_rigid_body_halozone" ) ;

  double radius = get_circumscribed_radius();

  geomVector const* pgc = get_ptr_to_gravity_centre();

  double r_equi = 3.0*radius;

  geomVector delta(r_equi, r_equi, r_equi);

  m_halo_zone[0]->operator=(*pgc);
  m_halo_zone[1]->operator=(*pgc);

  m_halo_zone[0]->operator-=(delta);
  m_halo_zone[1]->operator+=(delta);

}




//---------------------------------------------------------------------------
void DS_STL:: compute_surface_points( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: compute_surface_points" ) ;

  geomVector const* pgc = get_ptr_to_gravity_centre();
  geomVector normal(1.,1.,1.);

  m_surface_points[0]->operator=(*pgc);
  m_surface_area[0]->operator()(0) = 0.;
  m_surface_normal[0]->operator=(normal);



}




//---------------------------------------------------------------------------
void DS_STL:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: compute_number_of_surface_variables" ) ;

  // struct FS_STL_Additional_Param const* pagp =
  //  dynamic_cast<FS_STL*>(m_geometric_rigid_body)
  //     ->get_ptr_FS_STL_Additional_Param();
  //
  // size_t temp = (size_t) ((1./surface_cell_scale)
  //              *(4.*MAC::pi()*pagp->radius*pagp->radius)
  //              /(dx*dx));
  //
  // // Getting the nearest even number
  // Ntot = (size_t) (round((double)temp * 0.5) * 2.);
  Ntot = 1;

}




//---------------------------------------------------------------------------
bool DS_STL:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: isIn(pt)" ) ;

  return(isIn(pt(0), pt(1), pt(2)));

}




//---------------------------------------------------------------------------
bool DS_STL:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: isIn(x,y,z)" ) ;

  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double xC, yC, zC;

  double Xmin = m_MESH->get_main_domain_min_coordinate(0);
  double Xmax = m_MESH->get_main_domain_max_coordinate(0);
  double Ymin = m_MESH->get_main_domain_min_coordinate(1);
  double Zmin = m_MESH->get_main_domain_min_coordinate(2);
  double Zmax = m_MESH->get_main_domain_max_coordinate(2);

  double trdx = (Xmax - Xmin) / Nopx;
  double trdz = (Zmax - Zmin) / Nopz;

  double eps = 0.01*trdx;

  xC = x; yC = y; zC = z;

  int ii=-1; int kk=-1;

  // Searching for the right 2D Halo
  for (int i=0; i<Nopx; i++)
  {
     double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;

     for (int k=0; k<Nopz; k++)
     {
         double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;

         if ((xC > xmin-eps) && (xC < xmax+eps) && (zC > zmin-eps) && (zC < zmax+eps))
	 {
	     ii=i;
	     kk=k;
         }
     }
  }

  int itrscount = 0;
  double tri1[3], tri2[3], tri3[3], q1[3], q2[3];

  if (ii > -1 && kk > -1)
  {
     for (int m=0; m < tridx_xz->operator()(ii,kk)/3; m++)
     {
        // first point
	x1 = std::get<0>(ttrbox_xz[ii][kk][3*m  ]);
	y1 = std::get<1>(ttrbox_xz[ii][kk][3*m  ]);
	z1 = std::get<2>(ttrbox_xz[ii][kk][3*m  ]);

	tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;

	// second point
	x2 = std::get<0>(ttrbox_xz[ii][kk][3*m + 1]);
	y2 = std::get<1>(ttrbox_xz[ii][kk][3*m + 1]);
	z2 = std::get<2>(ttrbox_xz[ii][kk][3*m + 1]);
	tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;

	// third point
        x3 = std::get<0>(ttrbox_xz[ii][kk][3*m + 2]);
        y3 = std::get<1>(ttrbox_xz[ii][kk][3*m + 2]);
	z3 = std::get<2>(ttrbox_xz[ii][kk][3*m + 2]);
	tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

	q1[0] = xC; q1[1] = yC;        q1[2] = zC;
        q2[0] = xC; q2[1] = Ymin-1.0;  q2[2] = zC;

	if (intersect3d(q1,q2,tri1,tri2,tri3))
	{
	    itrscount = itrscount + 1;
	}
     }
  }

  bool level_set = false;

  if ( itrscount % 2 == 0 ) // even number -> outside
     level_set = true;

  if (invertSTL) level_set = !level_set;

  return(level_set);

}




//---------------------------------------------------------------------------
double DS_STL:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: level_set_value(pt)" ) ;

  bool flag = isIn(pt(0), pt(1), pt(2));

  double level_set = -1;
  if (flag)
     level_set = 1;

  return ( level_set );

}




//---------------------------------------------------------------------------
double DS_STL:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: level_set_value(x,y,z)" ) ;

  bool flag = isIn(x, y, z);

  double level_set = -1;
  if (flag)
     level_set = 1;

  return ( level_set );

}




//---------------------------------------------------------------------------
double DS_STL:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_distanceTo" ) ;

  double tri1[3], tri2[3], tri3[3], q1[3], q2[3];
  double x1,y1,z1,x2,y2,z2,x3,y3,z3;

  double Xmin = m_MESH->get_main_domain_min_coordinate(0);
  double Xmax = m_MESH->get_main_domain_max_coordinate(0);
  double Ymin = m_MESH->get_main_domain_min_coordinate(1);
  double Ymax = m_MESH->get_main_domain_max_coordinate(1);
  double Zmin = m_MESH->get_main_domain_min_coordinate(2);
  double Zmax = m_MESH->get_main_domain_max_coordinate(2);

  double trdx = (Xmax - Xmin) / Nopx;
  double trdy = (Ymax - Ymin) / Nopy;
  double trdz = (Zmax - Zmin) / Nopz;

  int ii = 0, kk = 0;

  double eps = 0.01*trdx;

  geomVector p1(3), p2(3);
  p1(0) = source(0);
  p1(1) = source(1);
  p1(2) = source(2);
  p2(0) = source(0) + delta*rayDir(0);
  p2(1) = source(1) + delta*rayDir(1);
  p2(2) = source(2) + delta*rayDir(2);

  double xcenter2 = 0.;

  // yz
  if ( rayDir(0) )
  {
 	  for (int i=0; i<Nopy; i++)
  	  {
	     double ymin = Ymin + i * trdy;
        double ymax = Ymin + (i+1) * trdy;

        for (int k=0; k<Nopz; k++)
	     {
           double zmin = Zmin + k * trdz;
           double zmax = Zmin + (k+1) * trdz;

           if ((p1(1) > ymin-eps) &&
               (p1(1) < ymax+eps) &&
               (p1(2) > zmin-eps) &&
               (p1(2) < zmax+eps))
           {
		       ii=i;
		       kk=k;
           }
	     }
     }

     for (int m=0; m<tridx_yz->operator()(ii,kk)/3; m++)
     {
	     // first point
        x1 = std::get<0>(ttrbox_yz[ii][kk][3*m  ]);
        y1 = std::get<1>(ttrbox_yz[ii][kk][3*m  ]);
		  z1 = std::get<2>(ttrbox_yz[ii][kk][3*m  ]);
		  tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;

		  // second point
	     x2 = std::get<0>(ttrbox_yz[ii][kk][3*m + 1]);
		  y2 = std::get<1>(ttrbox_yz[ii][kk][3*m + 1]);
		  z2 = std::get<2>(ttrbox_yz[ii][kk][3*m + 1]);
		  tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;

	     // third point
        x3 = std::get<0>(ttrbox_yz[ii][kk][3*m + 2]);
		  y3 = std::get<1>(ttrbox_yz[ii][kk][3*m + 2]);
		  z3 = std::get<2>(ttrbox_yz[ii][kk][3*m + 2]);
		  tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

        q1[0] = p1(0);  q1[1] = p1(1);  q1[2] = p1(2);
        q2[0] = p2(0); q2[1] = p2(1);  q2[2] = p2(2);

		  if (intersect3d(q1,q2,tri1,tri2,tri3))
		  {
	        // 974
		     double tmp1[3], tmp2[3], tmp3[3], tmp4[3];

           diffProduct(tri2,tri1,tmp1);
		     diffProduct(tri3,tri1,tmp2);
           double N[3];
		     crossProduct(tmp1,tmp2,N);
		     diffProduct(q1,tri1,tmp3);
		     diffProduct(q2,q1,tmp4);
		     double t = -dotProduct(tmp3,N)/dotProduct(tmp4,N);
		     xcenter2 = q1[0] + t * tmp4[0];
		  }
     }
  }
  // xz
  if ( rayDir(1) )
  {
 	for (int i=0; i<Nopx; i++)
  	{
	   double xmin = Xmin + i * trdx;
      double xmax = Xmin + (i+1) * trdx;

      for (int k=0; k<Nopz; k++)
	   {
         double zmin = Zmin + k * trdz;
         double zmax = Zmin + (k+1) * trdz;

         if ((p1(0) > xmin-eps) &&
             (p1(0) < xmax+eps) &&
             (p1(2) > zmin-eps) &&
             (p1(2) < zmax+eps))
	      {
		       ii=i;
		       kk=k;
         }
	   }
   }

   for (int m=0; m<tridx_xz->operator()(ii,kk)/3; m++)
	{
		// first point
      x1 = std::get<0>(ttrbox_xz[ii][kk][3*m  ]);
		y1 = std::get<1>(ttrbox_xz[ii][kk][3*m  ]);
		z1 = std::get<2>(ttrbox_xz[ii][kk][3*m  ]);
		tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;

		// second point
      x2 = std::get<0>(ttrbox_xz[ii][kk][3*m + 1]);
		y2 = std::get<1>(ttrbox_xz[ii][kk][3*m + 1]);
		z2 = std::get<2>(ttrbox_xz[ii][kk][3*m + 1]);
		tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;

    	// third point
      x3 = std::get<0>(ttrbox_xz[ii][kk][3*m + 2]);
		y3 = std::get<1>(ttrbox_xz[ii][kk][3*m + 2]);
		z3 = std::get<2>(ttrbox_xz[ii][kk][3*m + 2]);
		tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

      q1[0] = p1(0);  q1[1] = p1(1);  q1[2] = p1(2);
      q2[0] = p2(0); q2[1] = p2(1);  q2[2] = p2(2);

		if (intersect3d(q1,q2,tri1,tri2,tri3))
		{
		   double tmp1[3], tmp2[3], tmp3[3], tmp4[3];
	      diffProduct(tri2,tri1,tmp1);
		   diffProduct(tri3,tri1,tmp2);
         double N[3];
		   crossProduct(tmp1,tmp2,N);
		   diffProduct(q1,tri1,tmp3);
		   diffProduct(q2,q1,tmp4);
		   double t = -dotProduct(tmp3,N)/dotProduct(tmp4,N);
		   xcenter2 = q1[1] + t * tmp4[1];
		}
	}
  }

  // xy
  if ( rayDir(2) )
  {
 	for (int i=0; i<Nopx; i++)
  	{
	    double xmin = Xmin + i * trdx;
       double xmax = Xmin + (i+1) * trdx;

       for (int k=0; k<Nopy; k++)
	    {
          double ymin = Ymin + k * trdy;
          double ymax = Ymin + (k+1) * trdy;

          if ((p1(0) > xmin-eps) &&
              (p1(0) < xmax+eps) &&
              (p1(1) > ymin-eps) &&
              (p1(1) < ymax+eps))
          {
		       ii=i;
		       kk=k;
	       }
	    }
    }

    for (int m=0; m<tridx_xy->operator()(ii,kk)/3; m++)
	 {
	    // first point
       x1 = std::get<0>(ttrbox_xy[ii][kk][3*m  ]);
	    y1 = std::get<1>(ttrbox_xy[ii][kk][3*m  ]);
	    z1 = std::get<2>(ttrbox_xy[ii][kk][3*m  ]);
		 tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;

		 // second point
	    x2 = std::get<0>(ttrbox_xy[ii][kk][3*m + 1]);
	    y2 = std::get<1>(ttrbox_xy[ii][kk][3*m + 1]);
	    z2 = std::get<2>(ttrbox_xy[ii][kk][3*m + 1]);
		 tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;

	    // third point
       x3 = std::get<0>(ttrbox_xy[ii][kk][3*m + 2]);
	    y3 = std::get<1>(ttrbox_xy[ii][kk][3*m + 2]);
	    z3 = std::get<2>(ttrbox_xy[ii][kk][3*m + 2]);
		 tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

       q1[0] = p1(0);  q1[1] = p1(1);  q1[2] = p1(2);
       q2[0] = p2(0); q2[1] = p2(1);  q2[2] = p2(2);

		if (intersect3d(q1,q2,tri1,tri2,tri3))
		{
		   double tmp1[3], tmp2[3], tmp3[3], tmp4[3];
	      diffProduct(tri2,tri1,tmp1);
		   diffProduct(tri3,tri1,tmp2);
         double N[3];
		   crossProduct(tmp1,tmp2,N);
		   diffProduct(q1,tri1,tmp3);
		   diffProduct(q2,q1,tmp4);
		   double t = -dotProduct(tmp3,N)/dotProduct(tmp4,N);
		   xcenter2 = q1[2] + t * tmp4[2];
		}
	}
  }

  double dx = 0.;

  for (size_t i = 0; i < 3; i++) {
     if (rayDir(i) != 0)
        dx = MAC::abs(p1(i) - xcenter2);
  }

  return (dx);

}




//---------------------------------------------------------------------------
geomVector DS_STL:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: rigid_body_velocity(pt)" ) ;

  geomVector value(0.,0.,0.);

  return (value);

}




//---------------------------------------------------------------------------
geomVector DS_STL:: get_rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: rigid_body_angular_velocity()" ) ;

  geomVector value(0.,0.,0.);

  return (value);

}




//---------------------------------------------------------------------------
std::tuple<double,double> DS_STL:: get_mass_and_density() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_mass_and_density()" ) ;

  // return ( m_geometric_rigid_body->get_mass_and_density() );

}




//---------------------------------------------------------------------------
double DS_STL:: get_circumscribed_radius( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_circumscribed_radius()" ) ;

  double Xmin = m_MESH->get_main_domain_min_coordinate(0);
  double Xmax = m_MESH->get_main_domain_max_coordinate(0);
  double Ymin = m_MESH->get_main_domain_min_coordinate(1);
  double Ymax = m_MESH->get_main_domain_max_coordinate(1);
  double Zmin = m_MESH->get_main_domain_min_coordinate(2);
  double Zmax = m_MESH->get_main_domain_max_coordinate(2);
  
  double value = MAC::max(Xmax - Xmin, MAC::max(Ymax - Ymin, Zmax - Zmin));

  return(value);

}




//---------------------------------------------------------------------------
geomVector const* DS_STL:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_ptr_to_gravity_centre( )" ) ;

  double Xmin = m_MESH->get_main_domain_min_coordinate(0);
  double Xmax = m_MESH->get_main_domain_max_coordinate(0);
  double Ymin = m_MESH->get_main_domain_min_coordinate(1);
  double Ymax = m_MESH->get_main_domain_max_coordinate(1);
  double Zmin = m_MESH->get_main_domain_min_coordinate(2);
  double Zmax = m_MESH->get_main_domain_max_coordinate(2);

  ptr_gravity_centre->operator()(0) = 0.5*(Xmin + Xmax);
  ptr_gravity_centre->operator()(1) = 0.5*(Ymin + Ymax);
  ptr_gravity_centre->operator()(2) = 0.5*(Zmin + Zmax);

  return(ptr_gravity_centre);

}




//---------------------------------------------------------------------------
void DS_STL:: update_RB_position_and_velocity(geomVector const& pos,
                                                    geomVector const& vel,
                                                    geomVector const& ang_vel,
                                   vector<geomVector> const& periodic_directions)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: update_RB_position_and_velocity" ) ;


}




//---------------------------------------------------------------------------
double DS_STL:: dotProduct(double vect_A[], double vect_B[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: dotProduct" ) ;

  double product = 0;

  for (int i = 0; i < 3; i++)
     product = product + vect_A[i] * vect_B[i];

  return product;
}

//---------------------------------------------------------------------------
void DS_STL:: crossProduct(double vect_A[], double vect_B[]
		                          , double cross_P[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: crossProduct" ) ;

  cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
  cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
  cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

//---------------------------------------------------------------------------
void DS_STL:: diffProduct(double vect_A[], double vect_B[]
		                         , double diff_P[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: diffProduct" ) ;

  diff_P[0] = vect_A[0] - vect_B[0];
  diff_P[1] = vect_A[1] - vect_B[1];
  diff_P[2] = vect_A[2] - vect_B[2];
}

//---------------------------------------------------------------------------
int DS_STL:: orient3d(double vect_A[], double vect_B[], double vect_C[]
		                                 , double vect_D[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: orient3d" ) ;

  double tmp1[3], tmp2[3], tmp3[3], tmp4[3];

  diffProduct(vect_B,vect_A,tmp1);
  diffProduct(vect_C,vect_A,tmp2);
  diffProduct(vect_D,vect_A,tmp3);
  crossProduct(tmp1,tmp2,tmp4);
  if (dotProduct(tmp4,tmp3) > 0)
	  return 0;
  else
	  return 1;
}

//---------------------------------------------------------------------------
int DS_STL:: intersect3d(double q1[], double q2[], double tri1[]
		                    , double tri2[], double tri3[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: intersect3d" ) ;

  int s1 = orient3d(q1,tri1,tri2,tri3);
  int s2 = orient3d(q2,tri1,tri2,tri3);
  // Test whether the two extermities of the segment
  // are on the same side of the supporting plane of
  // the triangle
  if (s1 == s2)
     return 0;

  // Now we know that the segment 'straddles' the supporing
  // plane. We need to test whether the three tetrahedra formed
  // by the segment and the three edges of the triangle have
  // the same orientation
  int s3 = orient3d(q1,q2,tri1,tri2);
  int s4 = orient3d(q1,q2,tri2,tri3);
  int s5 = orient3d(q1,q2,tri3,tri1);
  return (s3 == s4 && s4 == s5);
}

//---------------------------------------------------------------------------
void DS_STL:: readSTL()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: readSTL" ) ;

  double xn, yn ,zn, x1,y1,z1,x2,y2,z2,x3,y3,z3;
  int kk = 0;
  string notsodummy,dummy;
  clock_t start;
  double duration;
  int dim = 3; // to be removed

  start = clock(); // start time for reading the set of triangles

  // reading STL file format
  // filename
  std::ostringstream os2;
  os2 << "./InputFiles/" << filename;
  std::string str_file = os2.str();

  // check if the file is ASCII or binary
  std::ifstream inFilep(str_file);

  int c;
  //if ( my_rank == is_master && field == 0 )
  {
     std::cout << endl;
     std::cout << "***" << endl;
     std::cout << "================= STL FILE =================" << endl;
	  std::cout << "   Read STL file: *" << str_file << "*" << endl;
  }

  while( (c = inFilep.get()) != EOF && c <= 127);

  //if ( my_rank == is_master && field == 0 )
  {
     if ( c == EOF )
	  {
        std::cout << "   File is in ASCII format" << endl;
     }
     else
	     std::cout << "   File is in Binary format" << endl;
  }

  inFilep.close();
  // end reading STL file format

  // reading triangulation
  if ( c == EOF ) // ASCII
  {
 	  std::ifstream inFile(str_file);
     //inFile.open(filename.c_str());
	  string line;

     if ( dim == 3 )
	     getline(inFile,line);

  	  while(getline(inFile,line))
  	  {
        std::istringstream iss1(line);
  	     iss1 >> notsodummy;
        if ( dim == 3  && notsodummy.compare("endsolid") == 0 )
           break;
	        iss1 >> dummy;
  	     iss1 >> xn;
  	     iss1 >> yn;
  	     iss1 >> zn;
        getline(inFile,line);
        getline(inFile,line);
        std::istringstream iss2(line);
        iss2 >> dummy;
  	     iss2 >> x1;
  	     iss2 >> y1;
  	     iss2 >> z1;
        getline(inFile,line);
        std::istringstream iss3(line);
        iss3 >> dummy;
  	     iss3 >> x2;
  	     iss3 >> y2;
   	  iss3 >> z2;
        if ( dim == 3 )
        {
           getline(inFile,line);
           std::istringstream iss4(line);
           iss4 >> dummy;
  	        iss4 >> x3;
  	        iss4 >> y3;
   	     iss4 >> z3;
        }
        getline(inFile,line);
        getline(inFile,line);

        // Vertices
        Llvls.push_back(std::make_tuple(x1, y1, z1));
        Llvls.push_back(std::make_tuple(x2, y2, z2));
        if ( dim == 3 )
           Llvls.push_back(std::make_tuple(x3, y3, z3));

        // Normals
        Llvns.push_back(std::make_tuple(xn, yn, zn));

        if ( dim == 2 )
           kk=kk+2;
        if ( dim == 3 )
           kk=kk+3;

     }

     inFile.close();
     Npls=kk;

  } // end ASCII
  else // binary
  {
     std::ifstream inFileb(str_file, std::ifstream::binary);

     // rdbuf returns a streambuf object associated with the
     // input fstream object ifs.

     std::filebuf* pbuf = inFileb.rdbuf();

     // Calculate the file's size.

     auto size = pbuf->pubseekoff(0, inFileb.end);

     // Set the position pointer to the beginning of the file.

     pbuf->pubseekpos(0);

     // Allocate memory to contain file data.

     char* buffer = new char[(size_t)size];

	  // Get file data. sgetn grabs all the characters from the streambuf
	  // object 'pbuf'. The return value of sgetn is the number of characters
	  // obtained - ordinarily, this value should be checked for equality
	  // against the number of characters requested.

	  pbuf->sgetn(buffer, size);

	  char * bufptr = buffer;

	  bufptr += 80;  // Skip past the header.
	  bufptr += 4;   // Skip past the number of triangles.

     while (bufptr < buffer + size)
     {

        xn = *(float *)(bufptr);
        yn = *(float *)(bufptr + 4);
        zn = *(float *)(bufptr + 8);
        bufptr += 12;

        Llvns.push_back(std::make_tuple(xn, yn, zn));

        x1 = *(float *)(bufptr);
   	  y1 = *(float *)(bufptr + 4);
   	  z1 = *(float *)(bufptr + 8);
   	  bufptr += 12;

   	  x2 = *(float *)(bufptr);
   	  y2 = *(float *)(bufptr + 4);
   	  z2 = *(float *)(bufptr + 8);
   	  bufptr += 12;

   	  x3 = *(float *)(bufptr);
   	  y3 = *(float *)(bufptr + 4);
   	  z3 = *(float *)(bufptr + 8);
   	  bufptr += 12;

        Llvls.push_back(std::make_tuple(x1, y1, z1));
        Llvls.push_back(std::make_tuple(x2, y2, z2));
        Llvls.push_back(std::make_tuple(x3, y3, z3));

        kk=kk+3;
        bufptr += 2;
     }

     inFileb.close();
     Npls=kk;
     delete[] buffer;
  }
  // end reading triangulation

  //if ( my_rank == is_master && field == 0 )
  {
      std::cout << "   Delaunay triangles: " << Npls/3 << endl;
      duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
      std::cout << "   File read in " << duration << "s" << endl;
  }

  // Filtering operations

  int filter_trbox = 1; // 1D distances approach

  if ( filter_trbox )
  {

     double Xmin = m_MESH->get_main_domain_min_coordinate(0);
     double Xmax = m_MESH->get_main_domain_max_coordinate(0);
     double Ymin = m_MESH->get_main_domain_min_coordinate(1);
     double Ymax = m_MESH->get_main_domain_max_coordinate(1);
     double Zmin = m_MESH->get_main_domain_min_coordinate(2);
     double Zmax = m_MESH->get_main_domain_max_coordinate(2);

     //if ( field == 0 && my_rank == is_master )
     {
    	   std::cout << "   - Box dimensions" << endl;
    	   std::cout << "   x [" << Xmin  << " " << Xmax << "]" <<  endl;
         std::cout << "   y [" << Ymin  << " " << Ymax << "]" <<  endl;
         std::cout << "   z [" << Zmin  << " " << Zmax << "]" <<  endl << endl;
     }

     double trdx = (Xmax - Xmin) / Nopx;
     double trdy = (Ymax - Ymin) / Nopy;
     double trdz = (Zmax - Zmin) / Nopz;

     double enhx = trdx * cenh;
     double enhy = trdy * cenh;
     double enhz = trdz * cenh;

        // xz
     for (int i=0; i<Nopx; i++)
     {
        double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;
        double xxmin = xmin - enhx;        double xxmax = xmax + enhx;

        for (int k=0; k<Nopz; k++)
        {
           double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;
           double zzmin = zmin - enhz;        double zzmax = zmax + enhz;

	        int m=0;

   	     for (int l=0; l<Npls/3; l++)
           {
              if  ( ((std::get<0>(Llvls[3*l])   > xxmin && std::get<0>(Llvls[3*l])   < xxmax)   ||
                     (std::get<0>(Llvls[3*l+1]) > xxmin && std::get<0>(Llvls[3*l+1]) < xxmax)   ||
   	               (std::get<0>(Llvls[3*l+2]) > xxmin && std::get<0>(Llvls[3*l+2]) < xxmax))  &&
   	              ((std::get<2>(Llvls[3*l])   > zzmin && std::get<2>(Llvls[3*l])   < zzmax)   ||
   	               (std::get<2>(Llvls[3*l+1]) > zzmin && std::get<2>(Llvls[3*l+1]) < zzmax)   ||
   	               (std::get<2>(Llvls[3*l+2]) > zzmin && std::get<2>(Llvls[3*l+2]) < zzmax)) )
   	        {
                  x1 = std::get<0>(Llvls[3*l  ]);
   	            x2 = std::get<0>(Llvls[3*l+1]);
                  x3 = std::get<0>(Llvls[3*l+2]);
   	            y1 = std::get<1>(Llvls[3*l  ]);
   	            y2 = std::get<1>(Llvls[3*l+1]);
                  y3 = std::get<1>(Llvls[3*l+2]);
                  z1 = std::get<2>(Llvls[3*l  ]);
   	            z2 = std::get<2>(Llvls[3*l+1]);
                  z3 = std::get<2>(Llvls[3*l+2]);

         		   ttrbox_xz[i][k].push_back(std::make_tuple(x1, y1, z1));
         		   ttrbox_xz[i][k].push_back(std::make_tuple(x2, y2, z2));
         		   ttrbox_xz[i][k].push_back(std::make_tuple(x3, y3, z3));

                  m=m+3;
   	        }
   	     }
           tridx_xz->operator()(i,k) = m;
       }
    }

    // xy
    for (int i=0; i<Nopx; i++)
    {
	    double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;
       double xxmin = xmin - enhx;        double xxmax = xmax + enhx;

       for (int k=0; k<Nopy; k++)
       {
          double ymin = Ymin + k * trdy;     double ymax = Ymin + (k+1) * trdy;
          double yymin = ymin - enhy;        double yymax = ymax + enhy;

	       int m=0;
          for (int l=0; l<Npls/3; l++)
          {
              if  ( ((std::get<0>(Llvls[3*l])   > xxmin && std::get<0>(Llvls[3*l])   < xxmax)   ||
	                  (std::get<0>(Llvls[3*l+1]) > xxmin && std::get<0>(Llvls[3*l+1]) < xxmax)   ||
	                  (std::get<0>(Llvls[3*l+2]) > xxmin && std::get<0>(Llvls[3*l+2]) < xxmax))  &&
		              ((std::get<1>(Llvls[3*l])   > yymin && std::get<1>(Llvls[3*l])   < yymax)   ||
	                  (std::get<1>(Llvls[3*l+1]) > yymin && std::get<1>(Llvls[3*l+1]) < yymax)   ||
	                  (std::get<1>(Llvls[3*l+2]) > yymin && std::get<1>(Llvls[3*l+2]) < yymax)) )
              {

                 x1 = std::get<0>(Llvls[3*l  ]);
                 x2 = std::get<0>(Llvls[3*l+1]);
                 x3 = std::get<0>(Llvls[3*l+2]);
                 y1 = std::get<1>(Llvls[3*l  ]);
                 y2 = std::get<1>(Llvls[3*l+1]);
                 y3 = std::get<1>(Llvls[3*l+2]);
                 z1 = std::get<2>(Llvls[3*l  ]);
                 z2 = std::get<2>(Llvls[3*l+1]);
                 z3 = std::get<2>(Llvls[3*l+2]);

      		     ttrbox_xy[i][k].push_back(std::make_tuple(x1, y1, z1));
      		     ttrbox_xy[i][k].push_back(std::make_tuple(x2, y2, z2));
      		     ttrbox_xy[i][k].push_back(std::make_tuple(x3, y3, z3));

                 m=m+3;
              }
           }

           tridx_xy->operator()(i,k) = m;
        }
     }

     // yz
     for (int i=0; i<Nopy; i++)
     {
        double ymin = Ymin + i * trdy;     double ymax = Ymin + (i+1) * trdy;
        double yymin = ymin - enhy;        double yymax = ymax + enhy;

        for (int k=0; k<Nopz; k++)
	     {
           double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;
           double zzmin = zmin - enhz;        double zzmax = zmax + enhz;

	        int m=0;
           for (int l=0; l<Npls/3; l++)
           {
              if  ( ((std::get<1>(Llvls[3*l])   > yymin && std::get<1>(Llvls[3*l])   < yymax)   ||
	                  (std::get<1>(Llvls[3*l+1]) > yymin && std::get<1>(Llvls[3*l+1]) < yymax)   ||
	                  (std::get<1>(Llvls[3*l+2]) > yymin && std::get<1>(Llvls[3*l+2]) < yymax))  &&
		              ((std::get<2>(Llvls[3*l])   > zzmin && std::get<2>(Llvls[3*l])   < zzmax)   ||
	                  (std::get<2>(Llvls[3*l+1]) > zzmin && std::get<2>(Llvls[3*l+1]) < zzmax)   ||
	                  (std::get<2>(Llvls[3*l+2]) > zzmin && std::get<2>(Llvls[3*l+2]) < zzmax)) )
	           {
   	           x1 = std::get<0>(Llvls[3*l  ]);
   	           x2 = std::get<0>(Llvls[3*l+1]);
                 x3 = std::get<0>(Llvls[3*l+2]);
                 y1 = std::get<1>(Llvls[3*l  ]);
	              y2 = std::get<1>(Llvls[3*l+1]);
                 y3 = std::get<1>(Llvls[3*l+2]);
	              z1 = std::get<2>(Llvls[3*l  ]);
	              z2 = std::get<2>(Llvls[3*l+1]);
                 z3 = std::get<2>(Llvls[3*l+2]);

      		     ttrbox_yz[i][k].push_back(std::make_tuple(x1, y1, z1));
      		     ttrbox_yz[i][k].push_back(std::make_tuple(x2, y2, z2));
      		     ttrbox_yz[i][k].push_back(std::make_tuple(x3, y3, z3));

                 m=m+3;
              }
	        }
           tridx_yz->operator()(i,k) = m;
        }
     }
   } // end filtertrbox
}
