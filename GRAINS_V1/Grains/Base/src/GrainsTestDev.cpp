#include "MPINeighbors.hh"
#include "GrainsTestDev.hh"
#include "GrainsBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include <stdlib.h>
#include <time.h>


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

    // Test data
    Vector3 v(1.,1.,1.);
    double alpha = 0., err = 0., cumerr = 0.;
    Vector3 unitvec, vrotm, vrotq;
    Matrix mrot;
    Quaternion qrot, qrot_conj;
    size_t n = 10000000;
    cout << "Number of random configurations = " << n << endl;

    // Quaternion to rotation matrix    
    cout << "Quaternion to rotation matrix" << endl;
    for (size_t i=0;i<n;++i)
    {
      alpha = 2. * PI * (double)rand() / RAND_MAX;
      unitvec = GrainsExec::RandomUnitVector( 3 );
      qrot.setQuaternion( sin( alpha / 2. ) * unitvec, cos( alpha / 2. ) );

      qrot_conj = qrot.Conjugate();
      vrotq = qrot.multToVector3( ( v , qrot_conj ) ); 
      //cout << "Vector rotation by quaternion" << endl;
      //cout << vrotq << endl;

      mrot.setRotation( qrot );
      vrotm = mrot * v;    
      //cout << "Vector rotation by matrix" << endl;
      //cout << vrotm << endl;
      
      err = Norm( vrotq - vrotm );
      cumerr += err;
      if ( err > EPSILON ) cout << "Error = " << err << endl;
    }
    cout << "Average error = " << cumerr / double(n) << endl;

    // Rotation matrix to quaternion   
    cout << "Rotation matrix to quaternion" << endl;    
    cumerr = 0.;
    for (size_t i=0;i<n;++i)
    {    
      mrot = GrainsExec::RandomRotationMatrix( 3 );
      //cout << mrot << endl;

      qrot.setQuaternion( mrot );    
      //cout << qrot << endl; 

      vrotm = mrot * v;    
      //cout << "Vector rotation by matrix" << endl;
      //cout << vrotm << endl;
    
      qrot_conj = qrot.Conjugate();
      vrotq = qrot.multToVector3( ( v , qrot_conj ) ); 
      //cout << "Vector rotation by quaternion" << endl;
      //cout << vrotq << endl;

      err = Norm( vrotq - vrotm );
      cumerr += err;
      if ( err > EPSILON ) cout << "Error = " << err << endl;
    } 
    cout << "Average error = " << cumerr / double(n) << endl;                   
  }
}
