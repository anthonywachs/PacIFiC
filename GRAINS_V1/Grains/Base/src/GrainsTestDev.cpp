#include "MPINeighbors.hh"
#include "GrainsTestDev.hh"
#include "GrainsBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "ContactBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "STLObstacle.hh"
#include "ConvexBuilderFactory.hh"
#include "PostProcessingWriterBuilderFactory.hh"
#include "RigidBody.hh"
#include "RigidBodyWithCrust.hh"
#include "Sphere.hh"
#include "Rectangle.hh"
#include <stdlib.h>
#include <time.h>
#include <random>
#include <thread>

#include <chrono>
#include "Box.hh"
#include "Cylinder.hh"
#include "Superquadric.hh"
#include "PointContact.hh"
#include "Particle.hh"
#include "KinematicsBuilderFactory.hh"


using namespace std;
using namespace std::chrono;


// ----------------------------------------------------------------------------
// Returns the solutions to to the quartic equation x^4 + bx^3 + cx^2 + dx + e
inline void solveQuartic( double const b, double const c, double const d,
                          double const e, double sol[4], int& nbRoots )
{
  // reseting the number of roots
  nbRoots = 0;
  // Deprressed quartic: y^4 + p*y^2 + q*y + r = 0
  double const b2 = b*b;
  double const p = c - 3.*b2/8.;
  double const q = b2*b/8. - b*c/2. + d;
  double const r = -3.*b2*b2/256. + e - b*d/4. + b2*c/16.;
  double const p2 = p*p;

  // Solve
  if ( fabs( q ) < EPSILON )
  {
    // finding solutions to the quadratic equation x^2 + px + r = 0.
    double const del = p2 / 4. - r; // this is actually del/4.!
    if ( del < 0. )
      return;
    else
    {
      double const m1 = - p / 2. + sqrt( del );
      double const m2 = - p / 2. - sqrt( del );
      if ( m1 > 0. )
      {
        sol[ nbRoots++ ] = sqrt( m1 ) - b / 4.;
        sol[ nbRoots++ ] = - sqrt( m1 ) - b / 4.;
      }
      if ( m2 > 0. )
      {
        sol[ nbRoots++ ] = sqrt( m2 ) - b / 4.;
        sol[ nbRoots++ ] = - sqrt( m2 ) - b / 4.;
      }
    }
  }
  else
  {
    // finding a real root to cubic equation x^3 + px^2 + (p*p/4. - r)x - q*q/8.
    double const u = -p2/36. - r/3.; // this is actually p/3.!
    double const v = -p2*p/216. + r*p/6. - q*q/16.; // this is actually v/2.!

    double const del = u*u*u + v*v;
    double m = 0.;
    if ( del < 0 )
      m = 2. * sqrt( -u ) * cos( acos( v / sqrt( -u ) / u ) / 3. ) - p / 3.;
    else
    {
      m = cbrt( -v + sqrt( del ) );
      m = m - u / m - p / 3.;
    }

    // roots
    if ( m < 0. )
      return;
    else
    {
      double const sqrt_mhalf = sqrt( m / 2. );
      double const first_var = - p / 2. - m / 2. - q / sqrt_mhalf / 4.;
      double const second_var = first_var + q / sqrt_mhalf / 2.;

      if ( first_var > 0. )
      {
        sol[ nbRoots++ ] = sqrt_mhalf + sqrt( first_var ) - b / 4.;
        sol[ nbRoots++ ] = sqrt_mhalf - sqrt( first_var ) - b / 4.;
      }
      if ( second_var > 0. )
      {
        sol[ nbRoots++ ] = - sqrt_mhalf + sqrt( second_var ) - b / 4.;
        sol[ nbRoots++ ] = - sqrt_mhalf - sqrt( second_var ) - b / 4.;
      }
    }
  }
}




// ----------------------------------------------------------------------------
// Sign function
template < typename T >
inline int sgn( T const val )
{
    return ( ( T(0) < val ) - ( val < T(0) ) );
}





// ----------------------------------------------------------------------------
// Default constructor
GrainsTestDev::GrainsTestDev() 
  : Grains()
{}




// ----------------------------------------------------------------------------
// Destructor
GrainsTestDev::~GrainsTestDev()
{}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsTestDev::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "Grains3D Test - Development - For developers only" << endl;
}



// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsTestDev::do_before_time_stepping( DOMElement* rootElement )
{
  Grains::do_before_time_stepping( rootElement );       
}






// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsTestDev::do_after_time_stepping()
{
  Grains::do_after_time_stepping();
}







// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition
void GrainsTestDev::Construction( DOMElement* rootElement )
{
  Grains::Construction( rootElement );
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsTestDev::Forces( DOMElement* rootElement )
{
  Grains::Forces( rootElement );
}





// ----------------------------------------------------------------------------
// Additional features of the simulation: time features, insertion, 
// post-processing
void GrainsTestDev::AdditionalFeatures( DOMElement* rootElement )
{
  Grains::AdditionalFeatures( rootElement );
}



// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsTestDev::Simulation( double time_interval )
{
// 
// //   int rankproc = 0, nprocs = 0;
// //   MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
// //   MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
// // 
// //   // Input file
// //   string filename = "Grains/Init/insert.xml";
// //   size_t error = 0;
// //   size_t pos = filename.find(".xml");
// //   if ( pos == string::npos )
// //   {
// //     cout << "ERROR : input file need the .xml extension" << endl;
// //     error = 1;
// //   }
// // 
// //   // Creating STL file and reading it
// //  // string str1 = "wall_double_ramp.stl";
// //   //string str1 = "wall_unit.stl";
// //   string str1 = "long_wall.stl";
// //     
// //   Obstacle *stlob = new STLObstacle( "testSTL", str1 );  
// // 
// //   cout << "Creating STLObstacle..." << endl;
// // 
// // 
// //   // ******************************
// //   //           TESTING 
// //   // ******************************
// // 
// //   int test1 = 0, test2 = 1;
// // 
// //   if (test1)
// //   {
// // 
// //   double xp = 0.5,   yp = 0.5,   zp = 0.5;
// //   double R = 0.1;
// // 
// //   // Test 1: projection belongs to triangle
// //   
// //   Point3 Pa(0., 0.5, 0.);
// //   Point3 Pb(-1000., 0.5, -1000.);
// //   Point3 Pc(-1., 0., 1.);
// //   Point3 Pd(-1., 0.1, 1.);
// //   Point3 P1(-1., 0., 1.);
// //   Point3 P2(1., 0., 1.);
// //   Point3 P3(0., 0., -1.);
// //  
// //   
// //   cout << "IsInter(Pa): " << STLObstacle::intersect(Pa, P1, P2, P3) << endl;
// //   cout << "IsInter(Pb): " << STLObstacle::intersect(Pb, P1, P2, P3) << endl;
// //   cout << "IsInter(Pc): " << STLObstacle::intersect(Pc, P1, P2, P3) << endl;
// //   cout << "IsInter(Pd): " << STLObstacle::intersect(Pd, P1, P2, P3) << endl;
// // 
// //   }
// // 
// //   cout << "Testing again..." << endl;
// // 
// //   // Test 2: triangle area
// //  
// //   if (test2)
// //   {
// // 
// //   double x1, y1, z1, x2, y2, z2, x3, y3, z3;
// // 
// //   x1 = 0.0;
// //   y1 = 0.0;
// //   z1 = 0.0;
// // 
// //   x2 = 1.0;
// //   y2 = 0.0;
// //   z2 = 0.0;
// // 
// //   x3 = 0.0;
// //   y3 = 0.0;
// //   z3 = 1.0;
// // 
// //   Vector3 n;
// //   n[0] = 0.; n[1] = 0.; n[2] = 0.;
// // 
// //   size_t vid = 0;
// //   STLVertex *v1 = new STLVertex(x1,y1,z1,n,vid); 
// //   STLVertex *v2 = new STLVertex(x2,y2,z2,n,vid);  
// //   STLVertex *v3 = new STLVertex(x3,y3,z3,n,vid);
// // 
// //   size_t tid = 0;
// //   tuple<STLVertex*,STLVertex*,STLVertex*> vetest;
// //   vetest = std::make_tuple(v1, v2, v3);
// //   STLTriangle trtest(vetest,n,tid); 
// // 
// //   cout << "Testing area computation: " << trtest.getSurfaceArea() << endl;
// // 
// //   }
// // 
// //   // Test 3: RBWC
// // //   
// // //   Rectangle *Rect = new Rectangle( 1.0, 1.0 );
// // //   Rect->ClosestPoint( RigidBodyWithCrust &Rect );


// //   // Input file
// //   string filename = "Grains/Init/insert.xml";
// //   size_t error = 0;
// //   size_t pos = filename.find(".xml");
// //   if ( pos == string::npos )
// //   {
// //     cout << "ERROR : input file need the .xml extension" << endl;
// //     error = 1;
// //   }
// // 

  
//   // Test: time integration of the angular motion of a cylinder submitted
//   // to q constant torque in the body-fixed space
//   // NOTE: requires to set the torque in the body-fixed space in 
//   // ParticleKinematics.cpp to the required value    
//   double vmax = 0., vmean = 0. ;
//   bool forcestats = false;
//   Particle const* pp = m_allcomponents.getParticle( 1 );
  
//   Vector3 omega, omb;
//   Quaternion qq, qqConj;
//   string savename = "dt" + to_string( (int) -log10( m_dt ) ) + ".txt";  
//   ofstream MyFile( savename, ios::out );
//   MyFile << std::fixed << setprecision(15);
  
//   qq = *( pp->getQuaternionRotation() );
//   qqConj = qq.Conjugate();
//   omega = *( pp->getAngularVelocity() );
//   omb = qqConj.multToVector3( ( omega , qq ) );
  
//   MyFile << m_time << " " << omb << " " << qq << endl;  

//   // Timers
//   SCT_insert_app( "ParticlesInsertion" );
//   SCT_insert_app( "ComputeForces" );
//   SCT_insert_app( "Move" );
//   SCT_insert_app( "UpdateParticleActivity" );
//   SCT_insert_app( "LinkUpdate" );
//   SCT_insert_app( "OutputResults" );

//   // Simulation: time marching algorithm
//   cout << "Time \t TO \tend \tParticles \tIn \tOut" << endl;
//   while ( m_tend - m_time > 0.01 * m_dt )
//   {
//     try
//     {
//       m_time += m_dt;
//       GrainsExec::m_time_counter++;
      

//       // Check whether data are output at this time
//       m_lastTime_save = false;
//       GrainsExec::m_output_data_at_this_time = false;
//       GrainsExec::m_postprocess_forces_at_this_time = false;
//       if ( m_timeSave != m_save.end() )
//         if ( *m_timeSave - m_time < 0.01 * m_dt )
// 	{
// 	  // Set the global data output boolean to true
// 	  GrainsExec::m_output_data_at_this_time = true;

// 	  // Next time of writing files
// 	  m_timeSave++;

// 	  m_lastTime_save = true;
// 	}
//       forcestats = m_collision->outputForceStatsAtThisTime( false, false );
//       if ( GrainsExec::m_output_data_at_this_time || forcestats )
//       {
//         GrainsExec::m_postprocess_forces_at_this_time = true;
// 	m_collision->resetPPForceIndex();
//       }
	

//       // Insertion of particles
//       SCT_set_start( "ParticlesInsertion" );      
//       m_npwait_nm1 = m_allcomponents.getNumberPhysicalParticlesToInsert();      
//       if ( m_insertion_mode == IM_OVERTIME )
//         insertParticle( m_insertion_order );      
//       m_allcomponents.computeNumberParticles( m_wrapper );
//       if ( m_npwait_nm1 
//       	!= m_allcomponents.getNumberPhysicalParticlesToInsert() )
//           cout << "\r                                              "
//                << "                 " << flush;
//       ostringstream oss;
//       oss.width(10);
//       oss << left << m_time;
//       cout << '\r' << oss.str() << "  \t" << m_tend << "\t\t\t"
//            << m_allcomponents.getNumberActiveParticlesOnProc() << '\t'
//            << m_allcomponents.getNumberPhysicalParticlesToInsert() 
// 	   << flush;
//       SCT_get_elapsed_time( "ParticlesInsertion" );


//       if ( GrainsExec::m_TIScheme == "SecondOrderLeapFrog" )
//       {
//         // We implement the kick-drift-kick version of the LeapFrog scheme

// 	// Move particles and obstacles
// 	// Update particle velocity over dt/2 and particle position over dt,
// 	// obstacle velocity and position over dt
// 	// v_i+1/2 = v_i + a_i * dt / 2
// 	// x_i+1 = x_i + v_i+1/2 * dt
//         moveParticlesAndObstacles( 0.5 * m_dt, m_dt, m_dt );
	
//         // Compute particle forces and torque
// 	// Compute f_i+1 and a_i+1 as a function of (x_i+1,v_i+1/2)
// 	SCT_set_start( "ComputeForces" );
//         computeParticlesForceAndTorque();
//         SCT_get_elapsed_time( "ComputeForces" );	
	
// 	// Update particle velocity over dt/2
// 	// v_i+1 = v_i+1/2 + a_i+1 * dt / 2 
// 	m_allcomponents.advanceParticlesVelocity( m_time, 0.5 * m_dt );
//       }
//       else
//       {
//         // Compute particle forces and torque
// 	// Compute f_i and a_i as a function of (x_i,v_i) 
// 	SCT_set_start( "ComputeForces" );
//         computeParticlesForceAndTorque();
//         SCT_get_elapsed_time( "ComputeForces" );
      
//         // Move particles and obstacles
// 	// x_i+1 = x_i + g(v_i,v_i-1,a_i,dt)
// 	// v_i+1 = v_i + g(a_i,a_i-1,dt)
//         moveParticlesAndObstacles( m_dt, m_dt, m_dt );
//       }


//       // Compute and write force & torque exerted on obstacles
//       m_allcomponents.computeObstaclesLoad( m_time, m_dt ); 
//       m_allcomponents.outputObstaclesLoad( m_time, m_dt, false, false, m_rank );
      
      
//       // Compute and write force statistics
//       if ( forcestats ) m_collision->outputForceStats( m_time, m_dt, m_rank, 
//       	m_wrapper );
        

//       // Write postprocessing and reload files
//       if ( GrainsExec::m_output_data_at_this_time )
//       {
// 	SCT_set_start( "OutputResults" );

// 	// Track component max and mean velocity
// 	m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean );
//         cout << endl << "Component velocity : max = " << vmax
//                << " average = " << vmean << endl;
//         fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
//                << "\t" << GrainsExec::doubleToString( ios::scientific, 6,
//                vmax ) << "\t" << GrainsExec::doubleToString( ios::scientific,
//                6, vmean ) << endl;

// 	// Display memory used by Grains
// 	display_used_memory();

// 	// Write reload files
// 	saveReload( m_time );

// 	// Write data
// 	qq = *( pp->getQuaternionRotation() );
//         qqConj = qq.Conjugate();
//         omega = *( pp->getAngularVelocity() );
//         omb = qqConj.multToVector3( ( omega , qq ) );	
//         MyFile << m_time << " " << omb << " " << qq << endl;  

// 	SCT_get_elapsed_time( "OutputResults" );
//       }
//     }
//     catch ( ContactError& errContact )
//     {
//       // Max overlap exceeded
//       cout << endl;
//       m_allcomponents.PostProcessingErreurComponents( "ContactError",
//             errContact.getComponents() );
//       errContact.Message( cout );
//       m_error_occured = true;
//       break;
//     }
//     catch ( MotionError& errMotion )
//     {
//       // Particle motion over dt is too large
//       cout << endl;
//       m_allcomponents.PostProcessingErreurComponents( "MotionError",
//             errMotion.getComponent() );
//       errMotion.Message(cout);
//       m_error_occured = true;
//       break;
//     }
//     catch ( SimulationError& errSimulation )
//     {
//       // Simulation error
//       cout << endl;
//       errSimulation.Message(cout);
//       m_error_occured = true;
//       break;
//     }
//   }

//   MyFile.close();

  GrainsExec::m_colDetBoundingVolume = 2;
  // ------
  double aX = 0., aY = M_PI / 2., aZ = 0.;
  // double arr1[12] = { cosl(aZ)*cosl(aY),
  //                     cosl(aZ)*sinl(aY)*sinl(aX) - sinl(aZ)*cosl(aX),
  //                     cosl(aZ)*sinl(aY)*cosl(aX) + sinl(aZ)*sinl(aX),
  //                     sinl(aZ)*cosl(aY),
  //                     sinl(aZ)*sinl(aY)*sinl(aX) + cosl(aZ)*cosl(aX),
  //                     sinl(aZ)*sinl(aY)*cosl(aX) - cosl(aZ)*sinl(aX),
  //                     -sinl(aY),
  //                     cosl(aY)*sinl(aX),
  //                     cosl(aY)*cosl(aX),
  //                     0., 0., 0. };
  double arr1[12] = { 0.90777818599267257, 
        -0.34394622422652388, -0.2400828187858734, -0.16991187547696518, 
        0.22177944270102096, -0.96017906317890767, 0.38349539711416725, 
        0.91242253021039654, 0.14288599216404918, 0.38349539711416725, 
        0.91242253021039654, 0.14288599216404918 };


        //               {0.68725783191399126, -0.72152872679487934, -0.0841009446001688, 0.52167773850176191, 0.4096726220059298, 
        // 0.74834529458698007, -0.50549877311875591, -0.55817975526576302, 0.65795619245272363, -0.071112092542979913, 0.017270385854537773, 
        // -1.1023320524838516}
  aY = 1. * M_PI / 180., aX = 0., aZ = 0.;
  // double arr2[12] = { cosl(aZ)*cosl(aY),
  //                     cosl(aZ)*sinl(aY)*sinl(aX) - sinl(aZ)*cosl(aX),
  //                     cosl(aZ)*sinl(aY)*cosl(aX) + sinl(aZ)*sinl(aX),
  //                     sinl(aZ)*cosl(aY),
  //                     sinl(aZ)*sinl(aY)*sinl(aX) + cosl(aZ)*cosl(aX),
  //                     sinl(aZ)*sinl(aY)*cosl(aX) - cosl(aZ)*sinl(aX),
  //                     -sinl(aY),
  //                     cosl(aY)*sinl(aX),
  //                     cosl(aY)*cosl(aX),
  //                     0., 0., 1. };
  double arr2[12] = { 0.96441649654648709, 
        -0.15269868060464081, 0.2158331163912981, -0.045978180344002272, 
        0.70704633989305199, 0.70567094327089375, -0.26035903695059331, 
        -0.69048431277341682, 0.67486634654036903, 0.0079795393380890246, -0.033813546577751048, 
        -1.1144920343735891 };

        //               {0.94129640473826892, -0.22655597386054568, 0.25026679591055367, 0.00052797640895790713, 0.74233994462123287, 
        // 0.67002322934453262, -0.33758080465001566, -0.63055832190881977, 0.69888153717492618, 
        // 0.0060606461208750533, -0.032212746771225431, -1.1155444231121658 };
  Transform const& a2w( arr1 );
  Transform const& b2w( arr2 );
  Convex* cvxA = new Cylinder( 0.035225199999999998, 0.070450399999999996 );
  Convex* cvxB = new Box( 1.56, 1.56, 0.005 );
  RigidBody rbA( cvxA, a2w );
  RigidBody rbB( cvxA, b2w );
  OBC const& obcA = (OBC const&) rbA.getBVolume();
  OBC const& obcB = (OBC const&) rbB.getBVolume();
  // cout << a2w << endl << b2w << endl;
  cout << obcA << endl << obcB << endl;

  // -----------------------------------------------------------
  // Variables
  double const r1 = obcA.getRadius();
  double const r2 = obcB.getRadius();
  double const h1 = obcA.getHeight() / 2.; // half-height
  double const h2 = obcB.getHeight() / 2.; // half-height
  Vector3 const eZ = a2w.getBasis() * obcA.getInitOrientation(); // e1 = eZ
  Vector3 e2 = b2w.getBasis() * obcB.getInitOrientation();

  // Variables to represent the 2nd cyl in the 1st cyl local coordinates
  Point3 xPt( *( b2w.getOrigin() ) - *( a2w.getOrigin() ) );
  Vector3 eX( ( e2 ^ eZ ).normalized() );
  Vector3 eY( ( eZ ^ eX ).normalized() );
  // Vector3 eZ( e1 );
  const double x = eX * xPt;
  const double y = eY * xPt;
  const double z = eZ * xPt;
  const double ey = eY * e2;
  const double ez = eZ * e2;


  /* Step one: Shortest distance -- Projection onto XY */
  if ( fabs( x ) >= r1 + r2 )
    cout << "check -1: " << r1 << " " << r2 << endl;
  else 
  {
    double s2 = y / ey;
    double s1 = z + s2 * ez;
    if ( fabs( s1 ) < h1 && fabs( s2 ) < h2 )
      cout << "check 0: " << s1 << " " << s2 << endl;
  }
  
  /* Step two: Projection onto YZ and check rectangles intersection */
  // We project the rectangles onto four axes; local coords of the
  // first and second rectangles.
  const double fey = fabs( ey ) + 1.e-5;
  const double fez = fabs( ez ) + 1.e-5;
  // onto y-axis
  if ( fabs( y ) > r1 + h2 * fey + r2 * fez )
    cout << "check 1: " << fabs( y ) << " " << r1 + h2 * fey + r2 * fez << endl;
  // onto z-axis
  if ( fabs( z ) > h1 + h2 * fez + r2 * fey )
    cout << "check 2: " << fabs( z ) << " " << h1 + h2 * fez + r2 * fey << endl;
  // onto e_normal
  if ( fabs( y * ez - z * ey ) > h2 + h1 * fey + r1 * fez )
    cout << "check 3: " << fabs( y * ez - z + ey ) << " " << h2 + h1 * fey + r1 * fez << endl;
  // onto e
  if ( fabs( y * ey + z * ez ) > r2 + h1 * fez + r1 * fey )
    cout << "check 4: " << fabs( y * ey + z + ez ) << " " << r2 + h1 * fez + r1 * fey << endl;



  /* Step three: Projection onto XY */
  // We check the overlap of an ellipse w/ a circle
  // First, secondary is an ellipse
  {
    const double topy = y + h2 * ey;
    const double bottomy = y - h2 * ey;
    if ( std::signbit( topy ) == std::signbit( bottomy ) )
    {
      const double x0 = x;
      const double y0 = abs( topy ) > abs( bottomy ) ? bottomy : topy;
      // check if origin is not in the ellipse
      if ( x0 * x0 + y0 * y0 / ez / ez > r2 * r2 )
      {
        const double ey_sq = ey * ey;
        // const double A = r2 * r2 * ey_sq * ey_sq;
        // const double B = 2. * y0 * ez / r2 / ey_sq;
        // const double C = x0 / A;
        // const double D = 0.25 * B * B;
        const double A = r2 * r2 * ey_sq * ey_sq;
        const double B = 2. * y0 * ez / r2 / ey_sq;
        const double C = x0 * x0 / A;
        const double D = y0 * y0 * ez * ez / A;

        double sint[4] = { 0. };
        int nbRoots;
        solveQuartic( -B, C + D - 1., B, -D, sint, nbRoots );
        // solveQuartic( B, C + D - 1., -B, -D, sint, nbRoots );


        cout << -B << " " << C + D - 1. << " " << B << " " << -D << endl;
        cout << -B << " " << C + D - 1. << " " << B << " " << -D << endl;
        for ( int i = 0; i < nbRoots; i++ )
        {
          if ( std::signbit( sint[i] ) != std::signbit( y0 ) )
          {
            const double cost = - sgn( x0 ) * sqrt( 1. - sint[i] * sint[i] );
            const double ptX = x0 + r2 * cost;
            const double ptY = y0 + r2 * ez * sint[i];
            if ( ptX * ptX + ptY * ptY > r1 * r1 )
              cout << "check 5:" << endl;
          }
        }
      }
    }
  }

  cout << "Relative Orientaion: " << eX * e2 << " " << ey << " " << ez << endl;
  cout << "Relative Center: " << x << " " << y << " " << z << endl;


  cout << "They are in OBC-Contact: " << isContactBVolume( (OBC const&) rbA.getBVolume(), 
                    (OBC const&) rbA.getBVolume(),
                              a2w, 
                              b2w ) << endl;

  Vector3 xVec = Vector3( xPt );
  cout << "They are in GJK-Contact: " << intersect( *cvxA, *cvxA, a2w, b2w, xVec ) << endl;
}
