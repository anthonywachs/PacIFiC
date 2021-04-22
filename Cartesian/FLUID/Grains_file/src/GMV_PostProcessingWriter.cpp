#include "GMV_PostProcessingWriter.hh"
#include "Particule.H"
#include "Obstacle.H"


/* Constructeur par defaut 
--------------------------*/
GMV_PostProcessingWriter::GMV_PostProcessingWriter(DOMNode* dn,
	int const& rank_,int const& nbranks_):
  PostProcessingWriter(dn,rank_,nbranks_),
  counter( 0 )
{
  string simul = ReaderXML::getNodeAttr_String(dn,"Name");
  GMVFilename = simul+"_GMV.";
  counter = ReaderXML::getNodeAttr_Int(dn,"InitialCycleNumber");  
}




/* Destructeur 
--------------*/
GMV_PostProcessingWriter::~GMV_PostProcessingWriter()
{}




/* Initialisation du post processeur 
------------------------------------*/
void GMV_PostProcessingWriter::PostProcessing_start(Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC,
	vector<Fenetre> const& insert_windows )
{
  if (m_rank == 0)
    one_output(temps,particules,obstacle);
}




/* Ecriture d'evolution 
-----------------------*/
void GMV_PostProcessingWriter::PostProcessing(Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC)
{
  if (m_rank == 0)
    one_output(temps,particules,obstacle);
}




/* Clot les ecritures 
---------------------*/
void GMV_PostProcessingWriter::PostProcessing_end()
{}




/* Ecriture
-----------*/
void GMV_PostProcessingWriter::one_output(Scalar const& temps,
  	list<Particule*> const* particules,
  	Obstacle *obstacle)
{
  // Nom de fichier de sortie GMV courant
  ostringstream oss;
  oss.width(6);
  oss.fill('0');
  oss << counter;    
  string currentGMVfile = GMVFilename+"dat."+oss.str();

  // GMV Obstacles
  ofstream GMVpost(currentGMVfile.c_str(),ios::out);
  GMVpost << "gmvinput ascii" << endl;
  GMVpost << "probtime " << temps << endl;
  GMVpost << "polygons" << endl;    
  obstacle->GMVoutput(GMVpost);
 
  // GMV Particules
  for (list<Particule*>::const_iterator particule=particules->begin(); 
	particule!=particules->end();particule++)
    if ((*particule)->getActivity() == COMPUTE)
      (*particule)->GMVoutput(GMVpost);
  GMVpost << "endpoly" << endl;
  GMVpost << "endgmv" << endl;    
  GMVpost.close(); 
  
  counter++;
}	 
