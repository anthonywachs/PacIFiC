#include "MyContact.H"
#include "Composant.H"
#include "Memento.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec un lot de parametres
MyContact::MyContact( map<string,double>& parameters ) :
  ContactLaw(),
  m_forceRugosite(0.)
{
  stiff = parameters["stiff"];
  muen  = parameters["mun"  ];
  muet  = parameters["mut"  ];
  muec  = parameters["muc"  ];
  k_s   = parameters["ks"   ];
  k_m_s = parameters["kms"  ];
}




// ----------------------------------------------------------------------------
// @brief Destructeur
MyContact::~MyContact()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces & moments de contact
bool MyContact::computeForces( Composant* p0_,
	Composant* p1_,
	PointContact &contactInfos,
	LinkedCell *LC,
	Scalar dt, int nbContact )
{
  // A FAIRE

  return false ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces & moments de contact pour le post processing
void MyContact::computeForcesPostProcessing( Composant* p0_,
	Composant* p1_, Scalar dt,
	PointContact &contactInfos,
	list<PointForcePostProcessing>* listOfContacts )
{
  // A FAIRE
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Define the parameters values for MyContact
map<string,double> MyContact::defineParameters( DOMNode* &root )
{
  map<string,double> parameters;

  DOMNode* parameter;
  double   value;
  parameter = ReaderXML::getNode(root, "stiff");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["stiff"]  = value;
  parameter = ReaderXML::getNode(root, "muc");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["muc"]    = value;
  parameter = ReaderXML::getNode(root, "mun");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["mun"]    = value;
  parameter = ReaderXML::getNode(root, "mut");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["mut"]    = value;
  parameter = ReaderXML::getNode(root, "ks");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["ks"]     = value;
  parameter = ReaderXML::getNode(root, "kms");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["kms"]    = value;
  parameter = ReaderXML::getNode(root, "color");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["color"]  = value;

  return parameters;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul une estimation du temps de contact et de la penetration
// maximale d'un choc elastique binaire avec dissipation d'energie, puis ecrit
// le resultat dans un flux de sortie
void MyContact::computeAndWriteEstimates( Composant* p0_,
	Composant* p1_,
  	const double &v0,
	ostream &OUT ) const
{
  // A FAIRE
}
