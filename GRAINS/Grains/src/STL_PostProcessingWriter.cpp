#include "STL_PostProcessingWriter.hh"
#include "Particule.H"
#include "Obstacle.H"


/* Constructeur par defaut 
--------------------------*/
STL_PostProcessingWriter::STL_PostProcessingWriter(DOMNode* dn,
	int const& rank_,int const& nbranks_):
  PostProcessingWriter(dn,rank_,nbranks_)
{
  simul = ReaderXML::getNodeAttr_String(dn,"Name");
  simul += ".stl";
}




/* Destructeur 
--------------*/
STL_PostProcessingWriter::~STL_PostProcessingWriter()
{}




/* Initialisation du post processeur 
------------------------------------*/
void STL_PostProcessingWriter::PostProcessing_start(Scalar const& temps, 
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
    one_output(temps,particules); 
}




/* Ecriture d'evolution 
-----------------------*/
void STL_PostProcessingWriter::PostProcessing(Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC)
{
  if (m_rank == 0)
    one_output(temps,particules); 
}




/* Clot les ecritures 
---------------------*/
void STL_PostProcessingWriter::PostProcessing_end()
{}





/* Initialisation du post processeur
------------------------------------*/
void STL_PostProcessingWriter::one_output(Scalar const& temps,
	list<Particule*> const* particules)
{
  list<Particule*>::const_iterator particule;
  
  ofstream fileOUT(simul.c_str(),ios::out);
  
  for (particule=particules->begin();particule!=particules->end();particule++)
  {
    fileOUT << "solid Particule" << (*particule)->getID() << endl;
    (*particule)->getForme()->write_convex_STL(fileOUT);
    fileOUT << "endsolid" << endl;    
  }
  
  fileOUT.close();
    
}  
