#include "Grains_Exec.hh"
#include "DLMFD_PostProcessingWriter.hh"
#include "Particule.H"
#include "Obstacle.H"


/* ---- A.ESTEGHAMATIAN - Aug 2015 - Creation --------------------------------------
PostProcessing writer for DLMFD solver (Therefore it is by default in serial).
For the moment, it is created in order to make a text output of time-averaged
contact force on each particle. (The equivalent method can be found in 
Text_PostProcessing) In general, it can be used
to produce any grains text output when DLMFD method is used. 
-----------------------------------------------------------------------------*/

bool DLMFD_PostProcessingWriter::b_cumulatedContactForce = false;
bool DLMFD_PostProcessingWriter::b_cumulatedLubriForce= false;
/* Constructeur par defaut 
--------------------------*/
DLMFD_PostProcessingWriter::DLMFD_PostProcessingWriter(DOMNode* dn,
	int const& rank_,int const& nbranks_):
  PostProcessingWriter(dn,rank_,nbranks_)
{ 
  string is_contact_output;
//  if ( ReaderXML::hasNodeAttr_String(dn,"SlipVelocityOutput") )

  if ( ReaderXML::hasNodeAttr_String(dn,"ContactForceOutput") )
    is_contact_output = ReaderXML::getNodeAttr_String(dn, "ContactForceOutput");
  
  if( is_contact_output == "yes" )
  {
     Grains_Exec::m_ContactforceOutput = true;
     b_cumulatedContactForce = true;
     if ( m_rank == 0 ) cout << "  With contact force output in text format \n";
  }  
  simul = ReaderXML::getNodeAttr_String(dn,"Name");  
  // For the lubrication force, if it is activated we make an output by default
  if (Grains_Exec::m_withlubrication)
  {
    if( m_rank == 0 ) cout << "  * cumulatedLubricationForce" << endl;
    b_cumulatedLubriForce = true;
  }
}




/* Destructeur 
--------------*/
DLMFD_PostProcessingWriter::~DLMFD_PostProcessingWriter()
{}




/* Initialisation du post processeur 
------------------------------------*/
void DLMFD_PostProcessingWriter::PostProcessing_start(
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
  if ( Grains_Exec::m_ReloadType == "new" ) 
      one_output_Standard( temps, particules, pwait );  
}




/* Ecriture d'evolution 
-----------------------*/
void DLMFD_PostProcessingWriter::PostProcessing( Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC )
{
  one_output_Standard( temps, particules, pwait );
}




/* Clot les ecritures 
---------------------*/
void DLMFD_PostProcessingWriter::PostProcessing_end()
{
  if ( m_rank == 0 )
  {
    if ( b_cumulatedContactForce )
    {
      contact_force_x.close();
      contact_force_y.close();   
      contact_force_z.close();
    }  
    coordination_number.close();
    if ( b_cumulatedLubriForce )
    {
      lubri_force_x.close();
      lubri_force_y.close();   
      lubri_force_z.close();
    }
  }
}




/* Efface les fichiers resultats
--------------------------------*/
void DLMFD_PostProcessingWriter::clearResultFiles() const
{
  if ( m_rank == 0 ) 
  {
    string cmd = "bash " + Grains_Exec::m_GRAINS_HOME 
     	+ "/ExecScripts/DLMFD_clear.exec " + simul;
    system( cmd.c_str() );
  }     
}




/* Creates output files and open streams
----------------------------------------*/
void DLMFD_PostProcessingWriter::prepareResultFiles( ios_base::openmode mode )
{  
  string file;  
  file = simul+"_coordinationNumber.res";
  coordination_number.open( file.c_str(), mode );  
  if ( b_cumulatedContactForce )
  {
    file = simul+"_contact_force_x.res";
    contact_force_x.open(file.c_str(),mode);
    file = simul+"_contact_force_y.res";
    contact_force_y.open(file.c_str(),mode);  
    file = simul+"_contact_force_z.res";
    contact_force_z.open(file.c_str(),mode);
  }        
  if ( b_cumulatedLubriForce)
  {
    file = simul+"_lubri_force_x.res";
    lubri_force_x.open(file.c_str(),mode);
    file = simul+"_lubri_force_y.res";
    lubri_force_y.open(file.c_str(),mode);  
    file = simul+"_lubri_force_z.res";
    lubri_force_z.open(file.c_str(),mode);
  }
}




/* Write down one output
------------------------*/
void DLMFD_PostProcessingWriter::one_output_Standard(Scalar const& temps,
	list<Particule*> const* particules,
  	list<Particule*> const* pwait)
{ 
  Vecteur const* contactF = NULL ;
  Vecteur const* lubriF = NULL ;
  Vecteur lubriFpwait(0.,0.,0.);

  if ( b_cumulatedContactForce )
  { 
    contact_force_x << temps;
    contact_force_y << temps;
    contact_force_z << temps;
  }  
  coordination_number << temps;
  if ( b_cumulatedLubriForce )
  { 
    lubri_force_x << temps;
    lubri_force_y << temps;
    lubri_force_z << temps;
  }  
   
  list<Particule*>::const_iterator particule;
  for (particule=particules->begin(); particule!=particules->end();particule++)
  {   
    // Nombre de contacts
    coordination_number << " " << (*particule)->getCoordinationNumber();
    
    // Contact Force 
    if ( b_cumulatedContactForce )
    {
      contactF = (*particule)->getForceContactPP();
      contact_force_x << " " << (*contactF)[X];
      contact_force_y << " " << (*contactF)[Y];
      contact_force_z << " " << (*contactF)[Z];
    }  
    // Lubrication Force 
    if ( b_cumulatedLubriForce )
    {
      lubriF = (*particule)->getForceLubriPP();
      lubri_force_x << " " << (*lubriF)[X];
      lubri_force_y << " " << (*lubriF)[Y];
      lubri_force_z << " " << (*lubriF)[Z];
    }
  }

  if ( pwait )  
    for (particule=pwait->begin(); particule!=pwait->end();particule++)
    {   
      // Nombre de contacts
      coordination_number << " 0";     
      
      // Contact Force 
      if ( b_cumulatedContactForce )
      {
        contactF = (*particule)->getForce();
        contact_force_x << " " << (*contactF)[X];
        contact_force_y << " " << (*contactF)[Y];
        contact_force_z << " " << (*contactF)[Z];             
      }
      if ( b_cumulatedLubriForce )
      {
        lubri_force_x << " " << lubriFpwait[X];
        lubri_force_y << " " << lubriFpwait[Y];
        lubri_force_z << " " << lubriFpwait[Z];             
      }
    }
  coordination_number << endl;
  if ( b_cumulatedContactForce )
  {
    contact_force_x << endl ;
    contact_force_y << endl ;
    contact_force_z << endl ;
  }
  if ( b_cumulatedLubriForce )
  {
    lubri_force_x << endl ;
    lubri_force_y << endl ;
    lubri_force_z << endl ;
  }
}  
