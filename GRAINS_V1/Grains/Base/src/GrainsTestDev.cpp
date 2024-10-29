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
//   int rankproc = 0, nprocs = 0;
//   MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
//   MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
// 
//   // Input file
//   string filename = "Grains/Init/insert.xml";
//   size_t error = 0;
//   size_t pos = filename.find(".xml");
//   if ( pos == string::npos )
//   {
//     cout << "ERROR : input file need the .xml extension" << endl;
//     error = 1;
//   }
// 
//   // Creating STL file and reading it
//  // string str1 = "wall_double_ramp.stl";
//   //string str1 = "wall_unit.stl";
//   string str1 = "long_wall.stl";
//     
//   Obstacle *stlob = new STLObstacle( "testSTL", str1 );  
// 
//   cout << "Creating STLObstacle..." << endl;
// 
// 
//   // ******************************
//   //           TESTING 
//   // ******************************
// 
//   int test1 = 0, test2 = 1;
// 
//   if (test1)
//   {
// 
//   double xp = 0.5,   yp = 0.5,   zp = 0.5;
//   double R = 0.1;
// 
//   // Test 1: projection belongs to triangle
//   
//   Point3 Pa(0., 0.5, 0.);
//   Point3 Pb(-1000., 0.5, -1000.);
//   Point3 Pc(-1., 0., 1.);
//   Point3 Pd(-1., 0.1, 1.);
//   Point3 P1(-1., 0., 1.);
//   Point3 P2(1., 0., 1.);
//   Point3 P3(0., 0., -1.);
//  
//   
//   cout << "IsInter(Pa): " << STLObstacle::intersect(Pa, P1, P2, P3) << endl;
//   cout << "IsInter(Pb): " << STLObstacle::intersect(Pb, P1, P2, P3) << endl;
//   cout << "IsInter(Pc): " << STLObstacle::intersect(Pc, P1, P2, P3) << endl;
//   cout << "IsInter(Pd): " << STLObstacle::intersect(Pd, P1, P2, P3) << endl;
// 
//   }
// 
//   cout << "Testing again..." << endl;
// 
//   // Test 2: triangle area
//  
//   if (test2)
//   {
// 
//   double x1, y1, z1, x2, y2, z2, x3, y3, z3;
// 
//   x1 = 0.0;
//   y1 = 0.0;
//   z1 = 0.0;
// 
//   x2 = 1.0;
//   y2 = 0.0;
//   z2 = 0.0;
// 
//   x3 = 0.0;
//   y3 = 0.0;
//   z3 = 1.0;
// 
//   Vector3 n;
//   n[0] = 0.; n[1] = 0.; n[2] = 0.;
// 
//   size_t vid = 0;
//   STLVertex *v1 = new STLVertex(x1,y1,z1,n,vid); 
//   STLVertex *v2 = new STLVertex(x2,y2,z2,n,vid);  
//   STLVertex *v3 = new STLVertex(x3,y3,z3,n,vid);
// 
//   size_t tid = 0;
//   tuple<STLVertex*,STLVertex*,STLVertex*> vetest;
//   vetest = std::make_tuple(v1, v2, v3);
//   STLTriangle trtest(vetest,n,tid); 
// 
//   cout << "Testing area computation: " << trtest.getSurfaceArea() << endl;
// 
//   }
// 
//   // Test 3: RBWC
// //   
// //   Rectangle *Rect = new Rectangle( 1.0, 1.0 );
// //   Rect->ClosestPoint( RigidBodyWithCrust &Rect );


//   // Input file
//   string filename = "Grains/Init/insert.xml";
//   size_t error = 0;
//   size_t pos = filename.find(".xml");
//   if ( pos == string::npos )
//   {
//     cout << "ERROR : input file need the .xml extension" << endl;
//     error = 1;
//   }
// 

  
  // Test: time integration of the angular motion of a cylinder submitted
  // to q constant torque in the body-fixed space
  // NOTE: requires to set the torque in the body-fixed space in 
  // ParticleKinematics.cpp to the required value    
  double vmax = 0., vmean = 0. ;
  bool forcestats = false;
  Particle const* pp = m_allcomponents.getParticle( 1 );
  
  Vector3 omega, omb;
  Quaternion qq, qqConj;
  string savename = "dt" + to_string( (int) -log10( m_dt ) ) + ".txt";  
  ofstream MyFile( savename, ios::out );
  MyFile << std::fixed << setprecision(15);
  
  qq = *( pp->getQuaternionRotation() );
  qqConj = qq.Conjugate();
  omega = *( pp->getAngularVelocity() );
  omb = qqConj.multToVector3( ( omega , qq ) );
  
  MyFile << m_time << " " << omb << " " << qq << endl;  

  // Timers
  SCT_insert_app( "ParticlesInsertion" );
  SCT_insert_app( "ComputeForces" );
  SCT_insert_app( "Move" );
  SCT_insert_app( "UpdateParticleActivity" );
  SCT_insert_app( "LinkUpdate" );
  SCT_insert_app( "OutputResults" );

  // Simulation: time marching algorithm
  cout << "Time \t TO \tend \tParticles \tIn \tOut" << endl;
  while ( m_tend - m_time > 0.01 * m_dt )
  {
    try
    {
      m_time += m_dt;
      GrainsExec::m_time_counter++;
      

      // Check whether data are output at this time
      m_lastTime_save = false;
      GrainsExec::m_output_data_at_this_time = false;
      GrainsExec::m_postprocess_forces_at_this_time = false;
      if ( m_timeSave != m_save.end() )
        if ( *m_timeSave - m_time < 0.01 * m_dt )
	{
	  // Set the global data output boolean to true
	  GrainsExec::m_output_data_at_this_time = true;

	  // Next time of writing files
	  m_timeSave++;

	  m_lastTime_save = true;
	}
      forcestats = m_collision->outputForceStatsAtThisTime( false, false );
      if ( GrainsExec::m_output_data_at_this_time || forcestats )
      {
        GrainsExec::m_postprocess_forces_at_this_time = true;
	m_collision->resetPPForceIndex();
      }
	

      // Insertion of particles
      SCT_set_start( "ParticlesInsertion" );      
      m_npwait_nm1 = m_allcomponents.getNumberPhysicalParticlesToInsert();      
      if ( m_insertion_mode == IM_OVERTIME )
        insertParticle( m_insertion_order );      
      m_allcomponents.computeNumberParticles( m_wrapper );
      if ( m_npwait_nm1 
      	!= m_allcomponents.getNumberPhysicalParticlesToInsert() )
          cout << "\r                                              "
               << "                 " << flush;
      ostringstream oss;
      oss.width(10);
      oss << left << m_time;
      cout << '\r' << oss.str() << "  \t" << m_tend << "\t\t\t"
           << m_allcomponents.getNumberActiveParticlesOnProc() << '\t'
           << m_allcomponents.getNumberPhysicalParticlesToInsert() 
	   << flush;
      SCT_get_elapsed_time( "ParticlesInsertion" );


      if ( GrainsExec::m_TIScheme == "SecondOrderLeapFrog" )
      {
        // We implement the kick-drift-kick version of the LeapFrog scheme

	// Move particles and obstacles
	// Update particle velocity over dt/2 and particle position over dt,
	// obstacle velocity and position over dt
	// v_i+1/2 = v_i + a_i * dt / 2
	// x_i+1 = x_i + v_i+1/2 * dt
        moveParticlesAndObstacles( 0.5 * m_dt, m_dt, m_dt );
	
        // Compute particle forces and torque
	// Compute f_i+1 and a_i+1 as a function of (x_i+1,v_i+1/2)
	SCT_set_start( "ComputeForces" );
        computeParticlesForceAndTorque();
        SCT_get_elapsed_time( "ComputeForces" );	
	
	// Update particle velocity over dt/2
	// v_i+1 = v_i+1/2 + a_i+1 * dt / 2 
	m_allcomponents.advanceParticlesVelocity( m_time, 0.5 * m_dt );
      }
      else
      {
        // Compute particle forces and torque
	// Compute f_i and a_i as a function of (x_i,v_i) 
	SCT_set_start( "ComputeForces" );
        computeParticlesForceAndTorque();
        SCT_get_elapsed_time( "ComputeForces" );
      
        // Move particles and obstacles
	// x_i+1 = x_i + g(v_i,v_i-1,a_i,dt)
	// v_i+1 = v_i + g(a_i,a_i-1,dt)
        moveParticlesAndObstacles( m_dt, m_dt, m_dt );
      }


      // Compute and write force & torque exerted on obstacles
      m_allcomponents.computeObstaclesLoad( m_time, m_dt ); 
      m_allcomponents.outputObstaclesLoad( m_time, m_dt, false, false, m_rank );
      
      
      // Compute and write force statistics
      if ( forcestats ) m_collision->outputForceStats( m_time, m_dt, m_rank, 
      	m_wrapper );
        

      // Write postprocessing and reload files
      if ( GrainsExec::m_output_data_at_this_time )
      {
	SCT_set_start( "OutputResults" );

	// Track component max and mean velocity
	m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean );
        cout << endl << "Component velocity : max = " << vmax
               << " average = " << vmean << endl;
        fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
               << "\t" << GrainsExec::doubleToString( ios::scientific, 6,
               vmax ) << "\t" << GrainsExec::doubleToString( ios::scientific,
               6, vmean ) << endl;

	// Display memory used by Grains
	display_used_memory();

	// Write reload files
	saveReload( m_time );

	// Write data
	qq = *( pp->getQuaternionRotation() );
        qqConj = qq.Conjugate();
        omega = *( pp->getAngularVelocity() );
        omb = qqConj.multToVector3( ( omega , qq ) );	
        MyFile << m_time << " " << omb << " " << qq << endl;  

	SCT_get_elapsed_time( "OutputResults" );
      }
    }
    catch ( ContactError& errContact )
    {
      // Max overlap exceeded
      cout << endl;
      m_allcomponents.PostProcessingErreurComponents( "ContactError",
            errContact.getComponents() );
      errContact.Message( cout );
      m_error_occured = true;
      break;
    }
    catch ( MotionError& errMotion )
    {
      // Particle motion over dt is too large
      cout << endl;
      m_allcomponents.PostProcessingErreurComponents( "MotionError",
            errMotion.getComponent() );
      errMotion.Message(cout);
      m_error_occured = true;
      break;
    }
    catch ( SimulationError& errSimulation )
    {
      // Simulation error
      cout << endl;
      errSimulation.Message(cout);
      m_error_occured = true;
      break;
    }
  }

  MyFile.close();
}
