/**
# Coupling interface for Grains3D and Basilisk 
*/

#include "GrainsBuilderFactory.hh"
#include "InterfaceGrains3DBasilisk.h"
#include <iostream>
#include <cstring>
#include <string>

#ifdef __cplusplus
extern "C" {
#endif
  
  static GrainsCoupledWithFluid* grains = NULL;

  
  void Init_Grains ( char const* inputfile, 
  	double fluid_density, const bool b_restart ) 
  {
    string simulation_file( inputfile );
    
    ReaderXML::initialize();
        
    string simulation_file_exe = GrainsBuilderFactory::init( simulation_file, 
    	0, 1 );
  
    DOMElement* rootNode = ReaderXML::getRoot( simulation_file_exe );
    
    grains = GrainsBuilderFactory::createCoupledWithFluid( rootNode, 
    	fluid_density );
    if ( b_restart ) grains->setReloadSame();
    grains->do_before_time_stepping( rootNode );
    ReaderXML::terminate();
    
    string cmd = "/bin/rm " + simulation_file_exe;
    system( cmd.c_str() ); 
         
    cout << "Construction of Grains completed" << endl;
  }




  void Simu_Grains( const double dt_fluid ) 
  {
    grains->Simulation( dt_fluid );	
  }

 
 
  
  char* GrainsToBasilisk( int* pstrsize )
  {
    // We use the interface function of PeliGRIFF
    istringstream iss;
    grains->GrainsToFluid( iss );

    // We remove the formatting and separate each entry by " "
    string buff;
    *pstrsize = int(iss.str().size()) + 1;    
    char* pstr = new char [*pstrsize];
    int pos = 0;
    while ( iss >> buff )
    {
      std::strcpy( &pstr[pos], buff.c_str() );
      pos += int(buff.size());
      std::strcpy( &pstr[pos], " " ); 
      pos += 1;     
    }    
    
    // Return the pointer to the char
    return pstr;
  } 




  void SaveResults_Grains() 
  {
    static unsigned int ppcounter = 0;
    if ( !ppcounter ) 
      grains->InitialPostProcessing( 6 );
    else
      grains->doPostProcessing( 6 );
    ++ppcounter; 
  }  




  void checkParaviewPostProcessing_Grains( char* solid_resDir )
  {    
    grains->checkParaviewPostProcessing( "grains", solid_resDir, true );
  }



  
  void UpdateVelocityGrains( double arrayv[][6], const int m, 
  	bool bsplit_explicit_acceleration ) 
  {
    // Transfer into a vector< vector<double> >
    vector<double> buf( 6, 0.);
    vector< vector<double> > vecv( m, buf );
    for (size_t i=0;i<size_t(m);++i)
      for (size_t j=0;j<6;++j)
        vecv[i][j] = arrayv[i][j];

    // We use the interface function of PeliGRIFF
    grains->updateParticlesVelocity( vecv, bsplit_explicit_acceleration );
  }  



  
  void SetInitialCycleNumber( int cycle0 ) 
  {
    grains->setInitialCycleNumber( cycle0 );
  }



  
  void SetInitialTime( double tinit ) 
  {
    grains->setInitialTime( tinit );
  }    
  
#ifdef __cplusplus
}
#endif
