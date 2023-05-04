#include "MPINeighbors.hh"
#include "GrainsTestDev.hh"
#include "GrainsBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include <stdlib.h>
#include <time.h>
#include <random>

#include <chrono>
#include "Box.hh"
#include "Cylinder.hh"
#include "PointContact.hh"

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
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement );
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping
void GrainsTestDev::do_after_time_stepping()
{}




// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain
// decomposition
void GrainsTestDev::Construction( DOMElement* rootElement )
{}




// ----------------------------------------------------------------------------
// External force definition
void GrainsTestDev::Forces( DOMElement* rootElement )
{}




// ----------------------------------------------------------------------------
// Additional features of the simulation: time features, insertion,
// post-processing
void GrainsTestDev::AdditionalFeatures( DOMElement* rootElement )
{}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsTestDev::Simulation( double time_interval )
{
  if ( m_processorIsActive )
  {
    // Use current time as seed for random generator
    srand(1444);

    // Comparing GJK and BCylinder
  //   int nIter = 1e6;

  //   // bodies
  //   Convex* convexA = new Box( 1., 1., 16. );
  //   // Convex* convexA  = new Cylinder( 1., 4. );
  //   // RigidBodyWithCrust RBWCa( convexA, trA );
  //   // RBWCa.setCrustThickness( 0.1 );

  //   Convex* convexB = new Box( 16., 1., 1. );
  //   // Convex* convexB = new Cylinder( 1., 4. );
  //   // RigidBodyWithCrust RBWCb( convexB, trB );
  //   // RBWCb.setCrustThickness( 0.1 );


  //   // Transformation
  //   double aX, aY, aZ;
  //   aX = 0.; aY = aZ = 0.;
  //   double const AAA[12] =
  //   { cos(aZ)*cos(aY), cos(aZ)*sin(aY)*sin(aX) - sin(aZ)*cos(aX), cos(aZ)*sin(aY)*cos(aX) + sin(aZ)*sin(aX),
  //     sin(aZ)*cos(aY), sin(aZ)*sin(aY)*sin(aX) + cos(aZ)*cos(aX), sin(aZ)*sin(aY)*cos(aX) - cos(aZ)*sin(aX),
  //     -sin(aY), cos(aY)*sin(aX), cos(aY)*cos(aX),
  //     0., 0., 0.};
  //   Transform const* trA = new Transform( AAA );

  //   aX = 0.; aY = 0.; aZ = M_PI/6.;
  //   double const BBB[12] =
  //   { cos(aZ)*cos(aY), cos(aZ)*sin(aY)*sin(aX) - sin(aZ)*cos(aX), cos(aZ)*sin(aY)*cos(aX) + sin(aZ)*sin(aX),
  //     sin(aZ)*cos(aY), sin(aZ)*sin(aY)*sin(aX) + cos(aZ)*cos(aX), sin(aZ)*sin(aY)*cos(aX) - cos(aZ)*sin(aX),
  //     -sin(aY), cos(aY)*sin(aX), cos(aY)*cos(aX),
  //     0., 0., 2.5};
  //   Transform const* trB = new Transform( BBB );


  //   int ctr = 1;
  //   bool cntct = false;

  //   // Call GJK
  //   {
  //     Point3 pointA, pointB;
  //     int nbIterGJK = 0;
  //     auto start = high_resolution_clock::now();
  //     for ( ctr = 1; ctr <= nIter && cntct == false; ctr++ )
  //     {
  //       double distance = closest_points( *convexA, *convexB, *trA, *trB,
  //                                         pointA, pointB, nbIterGJK );
  //       if ( distance <= 0 )
  //         cntct = true;
  //     }
  //     auto stop = high_resolution_clock::now();
  //     auto duration = duration_cast<microseconds>( stop - start );
  //     cout << "Contact " << cntct << ", Iter " << ctr << ", No. GJK Iter. "
  //          << nbIterGJK  << ", Time taken by GJK: " << duration.count() << " us" << endl;
  //   }

  //   // Call analCyl
  //   {
  //     PointContact pt;
  //     BCylinder bCylA = convexA->bcylinder();
  //     BCylinder bCylB = convexB->bcylinder();
  //     // Transform const* a2wNoCrust = rbA.getTransform();
  //     // Transform const* b2wNoCrust = rbB.getTransform();
  //     auto start = high_resolution_clock::now();
  //     for ( ctr = 1; ctr <= nIter && cntct == false; ctr++ )
  //     {
  //       // pt = intersect( bCylA, bCylB, *trA, *trB );
  //       cntct = isContact( bCylA, bCylB, *trA, *trB );
  //     }
  //     auto stop = high_resolution_clock::now();
  //     auto duration = duration_cast<microseconds>( stop - start );
  //     cout << "Contact " << cntct << ", Iter " << ctr <<
  //             ", Time taken by Cyl: " << duration.count() << " us" << endl;
  //   }

  }
}
