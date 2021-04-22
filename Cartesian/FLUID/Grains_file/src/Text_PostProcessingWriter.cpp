#include "Grains_Exec.hh"
#include "Text_PostProcessingWriter.hh"
#include "Particule.H"
#include "Obstacle.H"
#include "MPIWrapperGrains.hh"


/* ---- M.BERNARD - Janv 2013 - Creation --------------------------------------
Each processor (even the master one) sends particle kinematics.
The master proc receives messages from every proc and gathers them
 into "cinematique_Global" vector.
Then the master proc prints cinematique_Global values in corresponding
files.
We do not re-build all particles anymore on the master proc as it was done
with the previous version of PositionVitesse
-----------------------------------------------------------------------------*/

bool Text_PostProcessingWriter::b_totalForce = false;
bool Text_PostProcessingWriter::b_cumulatedContactForce = false;
bool Text_PostProcessingWriter::b_instantaneousContactForce = false;
bool Text_PostProcessingWriter::b_cumulatedLubriForce= false;
bool Text_PostProcessingWriter::b_hydroForce = false;
bool Text_PostProcessingWriter::b_slipVelocity = false;
bool Text_PostProcessingWriter::b_temperature = false;
bool Text_PostProcessingWriter::b_stressTensor = false;
bool Text_PostProcessingWriter::b_particleStressTensor = false;

/* Constructeur par defaut
--------------------------*/
Text_PostProcessingWriter::Text_PostProcessingWriter( DOMNode* dn,
    int const& rank_,int const& nbranks_ ):
  PostProcessingWriter(dn,rank_,nbranks_ )
{
  simul = ReaderXML::getNodeAttr_String( dn, "Name" );

  if( m_rank == 0 )
  {
    cout << "Text PostProcessingWriter(s) :" << endl;
    cout << "  * Position/Vitesse/Coordination (default ouputs)" << endl;
  }
  // For the lubrication force, if it is activated we make an output by default
  if (Grains_Exec::m_withlubrication)
  {
    if( m_rank == 0 ) cout << "  * cumulatedLubricationForce" << endl;
    b_cumulatedLubriForce = true;
  }

  if ( Grains_Exec::m_stressTensor )
    if ( m_rank == 0 )
    {
      b_stressTensor = true;
      cout << "  * Stress tensor" << endl;
    }

  DOMNodeList* allOutputs = ReaderXML::getNodes( dn );
  for( XMLSize_t i=0; i<allOutputs->getLength(); i++ )
  {
    DOMNode* nOutput = allOutputs->item( i );
    if( ReaderXML::getNodeName( nOutput ) == "totalForce" ||
             Grains_Exec::m_withCohesion )
    {
      b_totalForce = true;
      if( m_rank == 0 ) cout << "  * totalForce" << endl;
    }

    else if( ReaderXML::getNodeName( nOutput ) == "cumulatedContactForce" )
    {
      Grains_Exec::m_ContactforceOutput = true;
      b_cumulatedContactForce = true;
      if( m_rank == 0 ) cout << "  * cumulatedContactForce" << endl;
    }

    else if( ReaderXML::getNodeName( nOutput ) == "instantaneousContactForce" )
    {
      Grains_Exec::m_ContactforceOutput_instantaneous = true;
      b_instantaneousContactForce = true;
      if( m_rank == 0 ) cout << "  * instantaneousContactForce" << endl;
    }
    else if( ReaderXML::getNodeName( nOutput ) == "ParticleST" )
    {
      Grains_Exec::m_particleStressTensor = b_particleStressTensor = true;
      cout << "  * Individual Particle Stress tensor" << endl;
    }

    else if( ReaderXML::getNodeName( nOutput ) == "hydroForce" )
    {
      if( Grains_Exec::m_withHydroForce )
      {
        b_hydroForce = true;
        if( m_rank == 0 ) cout << "  * hydroForce" << endl;
      }
      else if( m_rank == 0 )
        cout << "WARNING : hydroForce Text PostProcessingWriter requested"
             << endl
             << "          whereas DRAG module is not activated." << endl
             << "    ==>   Output request will be neglected." << endl;
    }

    else if( ReaderXML::getNodeName( nOutput ) == "slipVelocity" )
    {
      if( Grains_Exec::m_withHydroForce )
      {
        b_slipVelocity = true;
        if( m_rank == 0 ) cout << "  * slipVelocity" << endl;
      }
      else if( m_rank == 0 )
        cout << "WARNING : slipVelocity Text PostProcessingWriter requested"
             << endl
             << "          whereas DRAG module is not activated." << endl
             << "    ==>   Output request will be neglected." << endl;
    }

    else if( ReaderXML::getNodeName( nOutput ) == "temperature" )
    {
      if( Grains_Exec::m_withFluidTemperature ||
          Grains_Exec::m_withSolidTemperature )
      {
        b_temperature = true;
        if( m_rank == 0 ) cout << "  * temperature" << endl;
      }
      else if( m_rank == 0 )
        cout << "WARNING : temperature Text PostProcessingWriter requested"
             << endl
             << "          whereas solid nor fluid temperature" << endl
             << "          module is activated." << endl
             << "    ==>   Output request will be neglected." << endl;
    }
  }

}




/* Destructeur
--------------*/
Text_PostProcessingWriter::~Text_PostProcessingWriter()
{}




/* Initialisation du post processeur
------------------------------------*/
void Text_PostProcessingWriter::PostProcessing_start(
    Scalar const& temps,
    Scalar const& dt,
    list<Particule*> const* particules,
    list<Particule*> const* pwait,
    list<Particule*> const* pperiodiques,
    vector<Particule*> const* ParticuleClassesReference,
    Obstacle *obstacle,
    LinkedCell const* LC,
    vector<Fenetre> const& insert_windows )
{
  MPIWrapperGrains const* wrapper = Grains_Exec::getComm() ;
  size_t nb_total_part = Grains_Exec::nbreParticulesOnAllProc() ;

  if ( m_rank == 0 )
  {
    ios_base::openmode mode = ios::app;
    if ( Grains_Exec::m_ReloadType == "new" )
    {
      mode = ios::out;
      clearResultFiles();
    }
    prepareResultFiles(mode);
  }

  // IF PARALLEL
  if( wrapper )
  {
    vector< vector<double> >* cinematique_Global;

    // R�unir la cinematique des particules de tous les procs sur le master
    cinematique_Global = wrapper->GatherPositionVitesse_PostProcessing(
    	  *particules, nb_total_part );

    // Ecrire les resultats contenus dans cinematique_Global
    if ( m_rank == 0 )
      if ( Grains_Exec::m_ReloadType == "new" )
        one_output_MPI( temps, nb_total_part, cinematique_Global ) ;

    // Gather particles class from every proc on master proc
    vector< vector<double> >* classes_Global =
        wrapper->GatherParticlesClass_PostProcessing( *particules,
            nb_total_part );

    // Write down particles class only once at the begining
    if ( m_rank == 0 )
    {
      for (size_t i=0; i<nb_total_part; i++)
        particle_class << (*classes_Global)[0][i] << " " ;
      particle_class << endl ;
    }

    if( cinematique_Global )
      delete cinematique_Global ;

    if( classes_Global )
      delete classes_Global ;
  }
  else
    one_output_Standard( temps, particules, pwait );
}




/* Ecriture d'evolution
-----------------------*/
void Text_PostProcessingWriter::PostProcessing( Scalar const& temps,
    Scalar const& dt,
    list<Particule*> const* particules,
    list<Particule*> const* pwait,
    list<Particule*> const* pperiodiques,
    vector<Particule*> const* ParticuleClassesReference,
    Obstacle *obstacle,
    LinkedCell const* LC )
{

  MPIWrapperGrains const* wrapper = Grains_Exec::getComm() ;
  size_t nb_total_part = Grains_Exec::nbreParticulesOnAllProc() ;
  // Following post-processings need a call to
  // ERContact::computeForcesPostProcessing, but we do not need
  // a list of PointForcePostProcessing in return (this method is
  // originally used in Paraview_PostProcessingWriter where the list of
  // contacts is needed). Since the method returns a list of
  // PointForcePostProcessing, we call it and then delete the pointer
  // right after. This can be improved in future (Amir)
  if ( b_instantaneousContactForce || b_stressTensor ||
          b_particleStressTensor )
  {
    list<struct PointForcePostProcessing>* pallContacts =
        LC->CalculerForcesPostProcessing( particules, dt );
    delete pallContacts;
  }

  if( wrapper )
  {
    vector< vector<double> >* cinematique_Global =
        wrapper->GatherPositionVitesse_PostProcessing( *particules,
        nb_total_part );

    // Ecrire les r�sultats contenus dans cinematique_Global
    if ( m_rank == 0 )
      one_output_MPI( temps, nb_total_part, cinematique_Global ) ;

    if( cinematique_Global )
      delete cinematique_Global ;
  }
  else
    one_output_Standard( temps, particules, pwait );
}




/* Clot les ecritures
---------------------*/
void Text_PostProcessingWriter::PostProcessing_end()
{
  if ( m_rank == 0 )
  {
    gc_coordinates_x.close();
    gc_coordinates_y.close();
    gc_coordinates_z.close();
    gc_velocity_x.close();
    gc_velocity_y.close();
    gc_velocity_z.close();
    gc_rotation_x.close();
    gc_rotation_y.close();
    gc_rotation_z.close();
    coordination_number.close();
    particle_class.close();
    if( b_totalForce )
    {
      gc_force_x.close();
      gc_force_y.close();
      gc_force_z.close();
    }
    if ( b_cumulatedContactForce )
    {
      contact_force_x.close();
      contact_force_y.close();
      contact_force_z.close();
    }
    if ( b_instantaneousContactForce )
    {
      contact_force_inst_x.close();
      contact_force_inst_y.close();
      contact_force_inst_z.close();
    }
    if ( b_cumulatedLubriForce )
    {
      lubri_force_x.close();
      lubri_force_y.close();
      lubri_force_z.close();
    }
    if ( b_hydroForce )
    {
      demcfd_HydroForce_x.close();
      demcfd_HydroForce_y.close();
      demcfd_HydroForce_z.close();
    }
    if ( b_slipVelocity )
    {
      demcfd_SlipVel_x.close();
      demcfd_SlipVel_y.close();
      demcfd_SlipVel_z.close();
    }
    if( b_temperature )
    {
      demcfd_ParticleTemperature.close();
      demcfd_ParticleHeatflux.close();
      demcfd_ParticleNusselt.close();
      demcfd_FluidTemperature.close();
    }
    if ( b_stressTensor )
      fIntMmt.close();
    if ( b_particleStressTensor )
    {
      fIntMmt0.close();
      fIntMmt1.close();
      fIntMmt2.close();
      fIntMmt3.close();
      fIntMmt4.close();
      fIntMmt5.close();
      fIntMmt6.close();
      fIntMmt7.close();
      fIntMmt8.close();
    }
  }
}





/* Initialisation du post processeur
------------------------------------*/
void Text_PostProcessingWriter::one_output_MPI(Scalar const& temps,
    size_t &nb_total_part,
    vector< vector<double> > const* cinematique_Global)
{
  int nTF=0, nCCF=0, nICF=0, nHF=0, nSV=0, nT=0, nCLF=0;
  vector<Scalar> InternalFeatures(9,0.);

  gc_coordinates_x << temps;
  gc_coordinates_y << temps;
  gc_coordinates_z << temps;
  gc_velocity_x << temps;
  gc_velocity_y << temps;
  gc_velocity_z << temps;
  gc_rotation_x << temps;
  gc_rotation_y << temps;
  gc_rotation_z << temps;

  coordination_number << temps;
  if( b_totalForce )
  {
    gc_force_x << temps;
    gc_force_y << temps;
    gc_force_z << temps;
    nTF = 3;
  }
  if ( b_cumulatedContactForce )
  {
    contact_force_x << temps;
    contact_force_y << temps;
    contact_force_z << temps;
    nCCF = 3;
  }
  if ( b_instantaneousContactForce )
  {
    contact_force_inst_x << temps;
    contact_force_inst_y << temps;
    contact_force_inst_z << temps;
    nICF = 3;
  }
  if( b_hydroForce )
  {
    demcfd_HydroForce_x << temps;
    demcfd_HydroForce_y << temps;
    demcfd_HydroForce_z << temps;
    nHF = 3;
  }
  if ( b_slipVelocity )
  {
    demcfd_SlipVel_x << temps;
    demcfd_SlipVel_y << temps;
    demcfd_SlipVel_z << temps;
    nSV = 3;
  }
  if( b_temperature )
  {
    demcfd_ParticleTemperature << temps;
    demcfd_ParticleHeatflux << temps;
    demcfd_ParticleNusselt << temps;
    demcfd_FluidTemperature << temps;
    nT = 4;
  }
  if ( b_cumulatedLubriForce )
  {
    lubri_force_x << temps;
    lubri_force_y << temps;
    lubri_force_z << temps;
    //nCLF = 3;
  }
  if ( b_stressTensor || b_particleStressTensor )
  {
    if ( b_stressTensor )
      fIntMmt << temps;
    else
    {
      fIntMmt0 << temps;
      fIntMmt1 << temps;
      fIntMmt2 << temps;
      fIntMmt3 << temps;
      fIntMmt4 << temps;
      fIntMmt5 << temps;
      fIntMmt6 << temps;
      fIntMmt7 << temps;
      fIntMmt8 << temps;
    }
  }

  // Dans le cas d'une insertion, tant que la particule n'est pas ins�r�e
  // la vitesse & la position sont nulles
  for (size_t i=0; i<nb_total_part; i++)
  {
    // Position du centre de gravit�
    gc_coordinates_x << " " << (*cinematique_Global)[0][i] ;
    gc_coordinates_y << " " << (*cinematique_Global)[1][i] ;
    gc_coordinates_z << " " << (*cinematique_Global)[2][i] ;

    // Vitesse translationnelle du centre de gravit�
    gc_velocity_x << " " << (*cinematique_Global)[3][i] ;
    gc_velocity_y << " " << (*cinematique_Global)[4][i] ;
    gc_velocity_z << " " << (*cinematique_Global)[5][i] ;

    // Vitesse de rotation du centre de gravit�
    gc_rotation_x << " " << (*cinematique_Global)[6][i] ;
    gc_rotation_y << " " << (*cinematique_Global)[7][i] ;
    gc_rotation_z << " " << (*cinematique_Global)[8][i] ;

    // Coordination number
    coordination_number << " " << (*cinematique_Global)[9][i] ;

    if( b_totalForce )
    {
      gc_force_x << " " << (*cinematique_Global)[10][i] ;
      gc_force_y << " " << (*cinematique_Global)[11][i] ;
      gc_force_z << " " << (*cinematique_Global)[12][i] ;
    }
    if ( b_cumulatedContactForce )
    {
      contact_force_x << " " << (*cinematique_Global)[10+nTF][i];
      contact_force_y << " " << (*cinematique_Global)[11+nTF][i];
      contact_force_z << " " << (*cinematique_Global)[12+nTF][i];
    }
    if ( b_instantaneousContactForce )
    {
      contact_force_inst_x << " " << (*cinematique_Global)[10+nTF+nCCF][i];
      contact_force_inst_y << " " << (*cinematique_Global)[11+nTF+nCCF][i];
      contact_force_inst_z << " " << (*cinematique_Global)[12+nTF+nCCF][i];
    }
    if ( b_hydroForce )
    {
      // Drag Force + eventual lift force
      demcfd_HydroForce_x << " " << (*cinematique_Global)[10+nTF+nCCF+nICF][i] ;
      demcfd_HydroForce_y << " " << (*cinematique_Global)[11+nTF+nCCF+nICF][i] ;
      demcfd_HydroForce_z << " " << (*cinematique_Global)[12+nTF+nCCF+nICF][i] ;
    }
    if ( b_slipVelocity )
    {
      demcfd_SlipVel_x << " " << (*cinematique_Global)[10+nTF+nCCF+nICF+nHF][i] ;
      demcfd_SlipVel_y << " " << (*cinematique_Global)[11+nTF+nCCF+nICF+nHF][i] ;
      demcfd_SlipVel_z << " " << (*cinematique_Global)[12+nTF+nCCF+nICF+nHF][i] ;
    }
    if ( b_temperature )
    {
      demcfd_ParticleTemperature << " "
          << (*cinematique_Global)[10+nTF+nCCF+nICF+nHF+nSV][i] ;
      demcfd_ParticleHeatflux << " "
          << (*cinematique_Global)[11+nTF+nCCF+nHF+nSV][i] ;
      demcfd_ParticleNusselt << " "
          << (*cinematique_Global)[12+nTF+nCCF+nHF+nSV][i] ;
      demcfd_FluidTemperature << " "
          << (*cinematique_Global)[13+nTF+nCCF+nICF+nHF+nSV][i] ;
    }
    if ( b_cumulatedLubriForce )
    {
      lubri_force_x << " " << (*cinematique_Global)[10+nTF+nCCF+nICF+nHF+nSV+nT][i];
      lubri_force_y << " " << (*cinematique_Global)[11+nTF+nCCF+nICF+nHF+nSV+nT][i];
      lubri_force_z << " " << (*cinematique_Global)[12+nTF+nCCF+nICF+nHF+nSV+nT][i];
    }
    if ( b_stressTensor || b_particleStressTensor )
    {
      if ( b_stressTensor )
      {
	InternalFeatures[0] +=
	    (*cinematique_Global)[10+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[1] +=
	    (*cinematique_Global)[11+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[2] +=
	    (*cinematique_Global)[12+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[3] +=
	    (*cinematique_Global)[13+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[4] +=
	    (*cinematique_Global)[14+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[5] +=
	    (*cinematique_Global)[15+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[6] +=
	    (*cinematique_Global)[16+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[7] +=
	    (*cinematique_Global)[17+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[8] +=
	    (*cinematique_Global)[18+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
      }
      else
      {
	InternalFeatures[0] =
	    (*cinematique_Global)[10+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[1] =
	    (*cinematique_Global)[11+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[2] =
	    (*cinematique_Global)[12+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[3] =
	    (*cinematique_Global)[13+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[4] =
	    (*cinematique_Global)[14+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[5] =
	    (*cinematique_Global)[15+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[6] =
	    (*cinematique_Global)[16+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[7] =
	    (*cinematique_Global)[17+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];
	InternalFeatures[8] =
	    (*cinematique_Global)[18+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][i];

	fIntMmt0 << " " << InternalFeatures[0];
	fIntMmt1 << " " << InternalFeatures[1];
	fIntMmt2 << " " << InternalFeatures[2];
	fIntMmt3 << " " << InternalFeatures[3];
	fIntMmt4 << " " << InternalFeatures[4];
	fIntMmt5 << " " << InternalFeatures[5];
	fIntMmt6 << " " << InternalFeatures[6];
	fIntMmt7 << " " << InternalFeatures[7];
	fIntMmt8 << " " << InternalFeatures[8];
      }
    }

  }
  gc_coordinates_x << endl ;
  gc_coordinates_y << endl ;
  gc_coordinates_z << endl ;
  gc_velocity_x << endl ;
  gc_velocity_y << endl ;
  gc_velocity_z << endl ;
  gc_rotation_x << endl ;
  gc_rotation_y << endl ;
  gc_rotation_z << endl ;
  coordination_number << endl ;
  if( b_totalForce )
  {
    gc_force_x << endl ;
    gc_force_y << endl ;
    gc_force_z << endl ;
  }
  if ( b_cumulatedContactForce )
  {
    contact_force_x << endl ;
    contact_force_y << endl ;
    contact_force_z << endl ;
  }
  if ( b_instantaneousContactForce )
  {
    contact_force_inst_x << endl ;
    contact_force_inst_y << endl ;
    contact_force_inst_z << endl ;
  }
  if( b_hydroForce )
  {
    demcfd_HydroForce_x << endl ;
    demcfd_HydroForce_y << endl ;
    demcfd_HydroForce_z << endl ;
  }
  if( b_slipVelocity )
  {
    demcfd_SlipVel_x << endl ;
    demcfd_SlipVel_y << endl ;
    demcfd_SlipVel_z << endl ;
  }
  if( b_temperature )
  {
    demcfd_ParticleTemperature << endl ;
    demcfd_ParticleHeatflux << endl ;
    demcfd_ParticleNusselt << endl ;
    demcfd_FluidTemperature << endl ;
  }
  if ( b_cumulatedLubriForce )
  {
    lubri_force_x << endl ;
    lubri_force_y << endl ;
    lubri_force_z << endl ;
  }
  if ( b_stressTensor )
  {
    fIntMmt << " " << InternalFeatures[0] << " "
            << " " << InternalFeatures[1] << " "
            << " " << InternalFeatures[2] << " "
            << " " << InternalFeatures[3] << " "
            << " " << InternalFeatures[4] << " "
            << " " << InternalFeatures[5] << " "
            << " " << InternalFeatures[6] << " "
            << " " << InternalFeatures[7] << " "
            << " " << InternalFeatures[8] << endl;
  }
  if ( b_particleStressTensor )
  {
    fIntMmt0 << endl;
    fIntMmt1 << endl;
    fIntMmt2 << endl;
    fIntMmt3 << endl;
    fIntMmt4 << endl;
    fIntMmt5 << endl;
    fIntMmt6 << endl;
    fIntMmt7 << endl;
    fIntMmt8 << endl;
  }
}




/* Efface les fichiers resultats
--------------------------------*/
void Text_PostProcessingWriter::clearResultFiles() const
{
  if ( m_rank == 0 )
  {
    string cmd = "bash " + Grains_Exec::m_GRAINS_HOME
        + "/ExecScripts/Text_clear.exec " + simul;
    system( cmd.c_str() );
  }
}




/* Creates output files and open streams
----------------------------------------*/
void Text_PostProcessingWriter::prepareResultFiles( ios_base::openmode mode )
{
  string file;
  file = simul+"_position_x.dat";
  gc_coordinates_x.open(file.c_str(),mode);
  file = simul+"_position_y.dat";
  gc_coordinates_y.open(file.c_str(),mode);
  file = simul+"_position_z.dat";
  gc_coordinates_z.open(file.c_str(),mode);

  file = simul+"_velocity_x.dat";
  gc_velocity_x.open(file.c_str(),mode);
  file = simul+"_velocity_y.dat";
  gc_velocity_y.open(file.c_str(),mode);
  file = simul+"_velocity_z.dat";
  gc_velocity_z.open(file.c_str(),mode);

  file = simul+"_rotation_x.dat";
  gc_rotation_x.open(file.c_str(),mode);
  file = simul+"_rotation_y.dat";
  gc_rotation_y.open(file.c_str(),mode);
  file = simul+"_rotation_z.dat";
  gc_rotation_z.open(file.c_str(),mode);

  file = simul+"_coordinationNumber.dat";
  coordination_number.open( file.c_str(), mode );

  if( b_totalForce )
  {
    file = simul+"_totalForce_x.dat";
    gc_force_x.open( file.c_str(), mode );
    file = simul+"_totalForce_y.dat";
    gc_force_y.open( file.c_str(), mode );
    file = simul+"_totalForce_z.dat";
    gc_force_z.open( file.c_str(), mode );
  }
  if ( b_cumulatedContactForce )
  {
    file = simul+"_cumulatedContactForce_x.dat";
    contact_force_x.open(file.c_str(),mode);
    file = simul+"_cumulatedContactForce_y.dat";
    contact_force_y.open(file.c_str(),mode);
    file = simul+"_cumulatedContactForce_z.dat";
    contact_force_z.open(file.c_str(),mode);
  }
  if ( b_instantaneousContactForce )
  {
    file = simul+"_instantaneousContactForce_x.dat";
    contact_force_inst_x.open(file.c_str(),mode);
    file = simul+"_instantaneousContactForce_y.dat";
    contact_force_inst_y.open(file.c_str(),mode);
    file = simul+"_instantaneousContactForce_z.dat";
    contact_force_inst_z.open(file.c_str(),mode);
  }
  if ( b_hydroForce )
  {
    file = simul+"_hydroForce_x.dat";
    demcfd_HydroForce_x.open( file.c_str(), mode );
    file = simul+"_hydroForce_y.dat";
    demcfd_HydroForce_y.open( file.c_str(), mode );
    file = simul+"_hydroForce_z.dat";
    demcfd_HydroForce_z.open( file.c_str(), mode );
  }
  if ( b_slipVelocity )
  {
    file = simul+"_slipVelocity_x.dat";
    demcfd_SlipVel_x.open( file.c_str(), mode );
    file = simul+"_slipVelocity_y.dat";
    demcfd_SlipVel_y.open( file.c_str(), mode );
    file = simul+"_slipVelocity_z.dat";
    demcfd_SlipVel_z.open( file.c_str(), mode );
  }
  if ( b_temperature )
  {
    file = simul+"_particleTemperature.dat";
    demcfd_ParticleTemperature.open( file.c_str(), mode );
    file = simul+"_particleHeatflux.dat";
    demcfd_ParticleHeatflux.open( file.c_str(), mode );
    file = simul+"_particleNusselt.dat";
    demcfd_ParticleNusselt.open( file.c_str(), mode );
    file = simul+"_fluidTemperature.dat";
    demcfd_FluidTemperature.open( file.c_str(), mode );
  }
  if ( b_cumulatedLubriForce )
  {
    file = simul+"_cumulatedLubricationForce_x.dat";
    lubri_force_x.open(file.c_str(),mode);
    file = simul+"_cumulatedLubricationForce_y.dat";
    lubri_force_y.open(file.c_str(),mode);
    file = simul+"_cumulatedLubricationForce_z.dat";
    lubri_force_z.open(file.c_str(),mode);
  }
  if ( b_stressTensor || b_particleStressTensor )
  {
    if ( b_stressTensor )
    {
      file = simul+"_SumInternalMoments.dat";
      fIntMmt.open( file.c_str(), mode );
    }
    else
    {
      file = simul+"_InternalMoment0.dat";
      fIntMmt0.open( file.c_str(), mode );
      file = simul+"_InternalMoment1.dat";
      fIntMmt1.open( file.c_str(), mode );
      file = simul+"_InternalMoment2.dat";
      fIntMmt2.open( file.c_str(), mode );
      file = simul+"_InternalMoment3.dat";
      fIntMmt3.open( file.c_str(), mode );
      file = simul+"_InternalMoment4.dat";
      fIntMmt4.open( file.c_str(), mode );
      file = simul+"_InternalMoment5.dat";
      fIntMmt5.open( file.c_str(), mode );
      file = simul+"_InternalMoment6.dat";
      fIntMmt6.open( file.c_str(), mode );
      file = simul+"_InternalMoment7.dat";
      fIntMmt7.open( file.c_str(), mode );
      file = simul+"_InternalMoment8.dat";
      fIntMmt8.open( file.c_str(), mode );
    }
  }

  file = simul+"_particleClass.dat";
  particle_class.open( file.c_str(), mode );

}




/* Write down one output
------------------------*/
void Text_PostProcessingWriter::one_output_Standard(Scalar const& temps,
    list<Particule*> const* particules,
    list<Particule*> const* pwait)
{
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;
  Vecteur const* demcfdHydroForce = NULL;
  Vecteur const* demcfdSlipVelocity = NULL;
  Vecteur const* force =  NULL;
  Vecteur const* contactF = NULL ;
  Vecteur const* contactF_inst = NULL ;
  Vecteur const* lubriF = NULL ;
  vector<Scalar> const* intMmt = NULL ;
  vector<Scalar> InternalFeatures(9,0.);
  Vecteur lubriFwait(0.,0.,0.);

  // Dans le cas d'une insertion, tant que la particule n'est pas ins�r�e
  // la vitesse & la position sont nulles
  // En particulier � t=0
  gc_coordinates_x << temps;
  gc_coordinates_y << temps;
  gc_coordinates_z << temps;
  gc_velocity_x << temps;
  gc_velocity_y << temps;
  gc_velocity_z << temps;
  gc_rotation_x << temps;
  gc_rotation_y << temps;
  gc_rotation_z << temps;
  coordination_number << temps;
  if( b_totalForce )
  {
    gc_force_x << temps;
    gc_force_y << temps;
    gc_force_z << temps;
  }
  if( b_cumulatedContactForce )
  {
    contact_force_x << temps;
    contact_force_y << temps;
    contact_force_z << temps;
  }
  if( b_instantaneousContactForce )
  {
    contact_force_inst_x << temps;
    contact_force_inst_y << temps;
    contact_force_inst_z << temps;
  }
  if( b_cumulatedLubriForce )
  {
    lubri_force_x << temps;
    lubri_force_y << temps;
    lubri_force_z << temps;
  }
  if ( b_hydroForce )
  {
    demcfd_HydroForce_x << temps;
    demcfd_HydroForce_y << temps;
    demcfd_HydroForce_z << temps;
  }
  if ( b_slipVelocity )
  {
    demcfd_SlipVel_x << temps;
    demcfd_SlipVel_y << temps;
    demcfd_SlipVel_z << temps;
  }
  if( b_temperature )
  {
    demcfd_ParticleTemperature << temps;
    demcfd_ParticleHeatflux << temps;
    demcfd_ParticleNusselt << temps;
    demcfd_FluidTemperature << temps;
  }
  if ( b_stressTensor )
    fIntMmt << temps;
  if ( b_particleStressTensor )
  {
    fIntMmt0 << temps;
    fIntMmt1 << temps;
    fIntMmt2 << temps;
    fIntMmt3 << temps;
    fIntMmt4 << temps;
    fIntMmt5 << temps;
    fIntMmt6 << temps;
    fIntMmt7 << temps;
    fIntMmt8 << temps;
  }
  // Dans le cas d'une insertion, tant que la particule n'est pas ins�r�e
  // la vitesse & la position sont nulles
  list<Particule*>::const_iterator particule;
  for (particule=particules->begin(); particule!=particules->end();particule++)
  {
    // Position du centre de gravit�
    centre = (*particule)->getPosition();
    gc_coordinates_x << " " << (*centre)[X];
    gc_coordinates_y << " " << (*centre)[Y];
    gc_coordinates_z << " " << scientific << (*centre)[Z];

    // Vitesse translationnelle du centre de gravit�
    vitesseT = (*particule)->getVitesseTranslation();
    gc_velocity_x << " " << (*vitesseT)[X];
    gc_velocity_y << " " << (*vitesseT)[Y];
    gc_velocity_z << " " << (*vitesseT)[Z];

    // Vitesse de rotation du centre de gravit�
    vitesseR = (*particule)->getVitesseRotation();
    gc_rotation_x << " " << (*vitesseR)[X];
    gc_rotation_y << " " << (*vitesseR)[Y];
    gc_rotation_z << " " << (*vitesseR)[Z];

    // Nombre de contacts
    coordination_number << " " << (*particule)->getCoordinationNumber();

    if( b_totalForce )
    {
      // Force
      force = (*particule)->getForce();
      gc_force_x << " " << (*force)[X];
      gc_force_y << " " << (*force)[Y] ;
      gc_force_z << " " << (*force)[Z] ;
    }
    if ( b_cumulatedContactForce )
    {
      contactF = (*particule)->getForceContactPP();
      contact_force_x << " " << (*contactF)[X];
      contact_force_y << " " << (*contactF)[Y];
      contact_force_z << " " << (*contactF)[Z];
    }
    if ( b_instantaneousContactForce )
    {
      contactF_inst = (*particule)->getForceContactPP_instantaneous();
      contact_force_inst_x << " " << (*contactF_inst)[X];
      contact_force_inst_y << " " << (*contactF_inst)[Y];
      contact_force_inst_z << " " << (*contactF_inst)[Z];
    }
    if( b_hydroForce )
    {
      // Drag Force + eventual lift force
      demcfdHydroForce = (*particule)->getParticleHydroForce();
      demcfd_HydroForce_x << " " << (*demcfdHydroForce)[X];
      demcfd_HydroForce_y << " " << (*demcfdHydroForce)[Y];
      demcfd_HydroForce_z << " " << (*demcfdHydroForce)[Z];
    }
    if ( b_slipVelocity )
    {
      demcfdSlipVelocity = (*particule)->getParticleSlipVel();
      demcfd_SlipVel_x << " " << (*demcfdSlipVelocity)[X];
      demcfd_SlipVel_y << " " << (*demcfdSlipVelocity)[Y];
      demcfd_SlipVel_z << " " << (*demcfdSlipVelocity)[Z];
    }
    if( b_temperature )
    {
      demcfd_ParticleTemperature << " "
                                 << (*particule)->get_solidTemperature();
      demcfd_ParticleHeatflux << " "
                                 << (*particule)->get_fluidSolidHeatFlux();
      demcfd_ParticleNusselt << " "
                                 << (*particule)->get_solidNusselt();
      demcfd_FluidTemperature << " "
                                 << (*particule)->get_DEMCFD_fluidTemperature();
    }
    if ( b_cumulatedLubriForce )
    {
      lubriF = (*particule)->getForceLubriPP();
      lubri_force_x << " " << (*lubriF)[X];
      lubri_force_y << " " << (*lubriF)[Y];
      lubri_force_z << " " << (*lubriF)[Z];
    }
    if ( b_stressTensor || b_particleStressTensor )
    {
      if ( b_stressTensor )
      {
	intMmt = (*particule)->getInternalMoment();
	// Write sum of internal moments to compute the macro. stress tensor
	InternalFeatures[0] += (*intMmt)[0];
	InternalFeatures[1] += (*intMmt)[1];
	InternalFeatures[2] += (*intMmt)[2];
	InternalFeatures[3] += (*intMmt)[3];
	InternalFeatures[4] += (*intMmt)[4];
	InternalFeatures[5] += (*intMmt)[5];
	InternalFeatures[6] += (*intMmt)[6];
	InternalFeatures[7] += (*intMmt)[7];
	InternalFeatures[8] += (*intMmt)[8];
      }
      else
      {
	intMmt = (*particule)->getStressTensor();
	fIntMmt0 << " " << (*intMmt)[0];
	fIntMmt1 << " " << (*intMmt)[1];
	fIntMmt2 << " " << (*intMmt)[2];
	fIntMmt3 << " " << (*intMmt)[3];
	fIntMmt4 << " " << (*intMmt)[4];
	fIntMmt5 << " " << (*intMmt)[5];
	fIntMmt6 << " " << (*intMmt)[6];
	fIntMmt7 << " " << (*intMmt)[7];
	fIntMmt8 << " " << (*intMmt)[8];
      }
    }
  }

  if ( b_stressTensor )
  {
    fIntMmt << "\t" << InternalFeatures[0] << "\t"
            << "\t" << InternalFeatures[1] << "\t"
            << "\t" << InternalFeatures[2] << "\t"
            << "\t" << InternalFeatures[3] << "\t"
            << "\t" << InternalFeatures[4] << "\t"
            << "\t" << InternalFeatures[5] << "\t"
            << "\t" << InternalFeatures[6] << "\t"
            << "\t" << InternalFeatures[7] << "\t"
            << "\t" << InternalFeatures[8] << endl;
  }

  if( pwait )
    for (particule=pwait->begin(); particule!=pwait->end();particule++)
    {
      // Position du centre de gravit�
      centre = (*particule)->getPosition();
      gc_coordinates_x << " " << (*centre)[X];
      gc_coordinates_y << " " << (*centre)[Y];
      gc_coordinates_z << " " << (*centre)[Z];

      // Vitesse translationnelle du centre de gravit�
      vitesseT = (*particule)->getVitesseTranslation();
      gc_velocity_x << " " << (*vitesseT)[X];
      gc_velocity_y << " " << (*vitesseT)[Y];
      gc_velocity_z << " " << (*vitesseT)[Z];

      // Vitesse de rotation du centre de gravit�
      vitesseR = (*particule)->getVitesseRotation();
      gc_rotation_x << " " << (*vitesseR)[X];
      gc_rotation_y << " " << (*vitesseR)[Y];
      gc_rotation_z << " " << (*vitesseR)[Z];

      // Nombre de contacts
      coordination_number << " 0";

      if( b_totalForce )
      {
        force = (*particule)->getForce();
        gc_force_x << " " << (*force)[X];
        gc_force_y << " " << (*force)[Y] ;
        gc_force_z << " " << (*force)[Z] ;
      }
      if ( b_cumulatedContactForce )
      {
        contactF = (*particule)->getForce();
        contact_force_x << " " << (*contactF)[X];
        contact_force_y << " " << (*contactF)[Y];
        contact_force_z << " " << (*contactF)[Z];
      }
      if ( b_instantaneousContactForce )
      {
        contactF_inst = (*particule)->getForce();
        contact_force_inst_x << " " << (*contactF)[X];
        contact_force_inst_y << " " << (*contactF)[Y];
        contact_force_inst_z << " " << (*contactF)[Z];
      }

      if( b_hydroForce )
      {
        demcfdHydroForce = (*particule)->getParticleHydroForce();
        demcfd_HydroForce_x << " " << (*demcfdHydroForce)[X];
        demcfd_HydroForce_y << " " << (*demcfdHydroForce)[Y];
        demcfd_HydroForce_z << " " << (*demcfdHydroForce)[Z];
      }
      if ( b_slipVelocity )
      {
        demcfdSlipVelocity = (*particule)->getParticleSlipVel();
        demcfd_SlipVel_x << " " << (*demcfdSlipVelocity)[X];
        demcfd_SlipVel_y << " " << (*demcfdSlipVelocity)[Y];
        demcfd_SlipVel_z << " " << (*demcfdSlipVelocity)[Z];
      }
      if( b_temperature )
      {
        demcfd_ParticleTemperature << " "
                                   << (*particule)->get_solidTemperature();
        demcfd_ParticleHeatflux << " "
                                   << (*particule)->get_fluidSolidHeatFlux();
	    demcfd_ParticleNusselt << " "
                                   << (*particule)->get_solidNusselt();
        demcfd_FluidTemperature << " "
                                 << (*particule)->get_DEMCFD_fluidTemperature();
      }
      if ( b_cumulatedLubriForce )
      {
        lubri_force_x << " " << lubriFwait[X];
        lubri_force_y << " " << lubriFwait[Y];
        lubri_force_z << " " << lubriFwait[Z];
      }
    }

  gc_coordinates_x << endl;
  gc_coordinates_y << endl;
  gc_coordinates_z << endl;
  gc_velocity_x << endl;
  gc_velocity_y << endl;
  gc_velocity_z << endl;
  gc_rotation_x << endl;
  gc_rotation_y << endl;
  gc_rotation_z << endl;
  coordination_number << endl;
  if( b_totalForce )
  {
    gc_force_x << endl;
    gc_force_y << endl;
    gc_force_z << endl;
  }
  if ( b_cumulatedContactForce )
  {
    contact_force_x << endl ;
    contact_force_y << endl ;
    contact_force_z << endl ;
  }
  if ( b_instantaneousContactForce )
  {
    contact_force_inst_x << endl ;
    contact_force_inst_y << endl ;
    contact_force_inst_z << endl ;
  }
  if ( b_cumulatedLubriForce )
  {
    lubri_force_x << endl ;
    lubri_force_y << endl ;
    lubri_force_z << endl ;
  }
  if( b_hydroForce )
  {
    demcfd_HydroForce_x << endl;
    demcfd_HydroForce_y << endl;
    demcfd_HydroForce_z << endl;
  }
  if( b_slipVelocity )
  {
    demcfd_SlipVel_x << endl;
    demcfd_SlipVel_y << endl;
    demcfd_SlipVel_z << endl;
  }
  if( b_temperature )
  {
    demcfd_ParticleTemperature << endl;
    demcfd_ParticleHeatflux << endl;
    demcfd_ParticleNusselt << endl;
    demcfd_FluidTemperature << endl;
  }
  if ( b_stressTensor )
    fIntMmt << endl;
  if ( b_particleStressTensor )
  {
    fIntMmt0 << endl;
    fIntMmt1 << endl;
    fIntMmt2 << endl;
    fIntMmt3 << endl;
    fIntMmt4 << endl;
    fIntMmt5 << endl;
    fIntMmt6 << endl;
    fIntMmt7 << endl;
    fIntMmt8 << endl;
  }
}
