#include <FS_Grains3DPlugIn.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include "GrainsBuilderFactory.hh"
#include "GrainsCoupledWithFluid.hh"
#include "ReaderXML.hh"


//---------------------------------------------------------------------------
FS_Grains3DPlugIn::FS_Grains3DPlugIn( string const& insertion_file_,
        string const& simulation_file_,
        double const& fluid_density,
	bool const& correct_particle_acceleration,
        bool const& b_restart,
        double const& grid_size,
        bool const& is_solidsolver_parallel,	
        int& error )
//---------------------------------------------------------------------------
  : FS_SolidPlugIn( insertion_file_, simulation_file_ )
  , m_Grains3D( NULL )
  , m_is_master( 0 )
  , m_Grains3D_parallel_mode( is_solidsolver_parallel )
{
  MAC_LABEL( "FS_Grains3DPlugIn:: FS_Grains3DPlugIn" ) ;

  // Set the MPI features
  m_macCOMM = MAC_Exec::communicator();
  m_my_rank = m_macCOMM->rank();
  m_nb_ranks = m_macCOMM->nb_ranks();
  m_Grains3D_active_on_this_rank = m_Grains3D_parallel_mode 
  	|| ( m_my_rank == m_is_master ); 
  
  // Output message
  if ( m_my_rank == m_is_master )
    MAC::out() << "Construction of Grains3D with " << m_simulation_file << endl;
  
  if ( m_Grains3D_active_on_this_rank )
  { 
    // Initialize XML reader
    ReaderXML::initialize();

    // Create an input file with the XML requirements from simulation_file
    // Only done by the master process 0
    string simulation_file_exe = 
      GrainsBuilderFactory::init( m_simulation_file, m_my_rank, 
        m_Grains3D_parallel_mode ? m_nb_ranks : 1 );

    // Create the Grains3D coupled to the fluid application
    DOMElement* rootNode = ReaderXML::getRoot( simulation_file_exe );
    m_Grains3D = GrainsBuilderFactory::createCoupledWithFluid( rootNode,
        fluid_density );
    m_Grains3D->setFluidCorrectedAcceleration( correct_particle_acceleration );
    if ( b_restart ) m_Grains3D->setReloadSame() ;
    m_Grains3D->do_before_time_stepping( rootNode );

    // Finalize XML reader
    ReaderXML::terminate(); 

    // Delete input file with XML requirements, done by master process only
    if ( m_my_rank == m_is_master ) 
    { 
      string cmd = "/bin/rm " + simulation_file_exe;
      system( cmd.c_str() );
    }
  }
    
  // Output message  
  if ( m_my_rank == m_is_master )  
    MAC::out() << "Construction of Grains3D completed" << endl;  

}




//---------------------------------------------------------------------------
FS_Grains3DPlugIn:: ~FS_Grains3DPlugIn()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: ~FS_Grains3DPlugIn" ) ;

  if ( m_Grains3D ) delete m_Grains3D;
  
}




//---------------------------------------------------------------------------
void FS_Grains3DPlugIn:: Simulation( double const& time_interval,
	bool const& predictor,
        bool const& isPredictorCorrector,
        double const& contact_force_coef,
        bool const& explicit_added_mass )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: Simulation" ) ;

  if ( m_Grains3D_active_on_this_rank )
    m_Grains3D->Simulation( time_interval ); 

}




// //---------------------------------------------------------------------------
// void FS_Grains3DPlugIn:: getSolidBodyFeatures( 
// 	vector<string>& solidbodyfeatures )
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL( "FS_Grains3DPlugIn:: getSolidBodyFeatures" ) ;
// 
//   if ( m_Grains3D_active_on_this_rank )
//     m_Grains3D->WriteParticulesInFluid( solidbodyfeatures );  
// 
// }




//---------------------------------------------------------------------------
void FS_Grains3DPlugIn:: getSolidBodyFeatures( istringstream* & is )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: getSolidBodyFeatures" ) ;

  // First delete the content of is and recreate a new istringstream 
  if ( !is ) delete is;
  is = new istringstream;
  
  // Get solid body features from Grains3D
  if ( m_Grains3D_active_on_this_rank )
    m_Grains3D->GrainsToFluid( *is ); 
    
//  cout << is->str() << endl;  
    
  // If Grains3D runs in serial and the fluid runs in parallel, we need
  // to broadcast the solid body features stream to all processes
  if ( !m_Grains3D_parallel_mode )
  {
    string transferString;
    if ( m_my_rank == m_is_master ) transferString = is->str();
    m_macCOMM->broadcast( transferString, m_is_master );
    is->str( transferString );   
  }   

}




//---------------------------------------------------------------------------
void FS_Grains3DPlugIn:: saveResults( string const& filename, 
	double const& time,
        int const& cycleNumber )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: saveResults" ) ;

  static int counter = 0;

  if ( m_Grains3D_active_on_this_rank )
  {
    if ( !counter ) 
    {
      m_Grains3D->setInitialCycleNumber( cycleNumber );
      m_Grains3D->InitialPostProcessing();
    }
    else m_Grains3D->doPostProcessing();
  }

  ++counter;

}




//---------------------------------------------------------------------------
void FS_Grains3DPlugIn:: transferHydroFTtoSolid( 
      	vector< vector<double> > const* hydroFT ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: transferHydroFTtoSolid" ) ;
  
  if ( m_Grains3D_active_on_this_rank )
    m_Grains3D->updateParticlesHydroFT( hydroFT );
  
}  




//---------------------------------------------------------------------------
void FS_Grains3DPlugIn:: checkParaviewPostProcessing( 
	string const& solid_resDir ) 
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: checkParaviewPostProcessing" ) ;
  
  if ( m_Grains3D_active_on_this_rank )
    m_Grains3D->checkParaviewPostProcessing( 
             "grains", solid_resDir, true ) ;
  
} 




//---------------------------------------------------------------------------
void FS_Grains3DPlugIn:: setParaviewPostProcessingTranslationVector( 
      	double const& tvx, double const& tvy, double const& tvz )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: setParaviewPostProcessingTranslationVector" );
  
  if ( m_Grains3D_active_on_this_rank )
    m_Grains3D->setParaviewPostProcessingTranslationVector( 
      	tvx, tvy, tvz );
  
} 




//---------------------------------------------------------------------------
void FS_Grains3DPlugIn:: setInitialTime( double const& time0 )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Grains3DPlugIn:: setInitialTime" );
  
  if ( m_Grains3D_active_on_this_rank )
    m_Grains3D->setInitialTime( time0 );
  
} 
