#include "Grains_Exec.hh"
#include "Contact_BuilderFactory.hh"
#include "ContactLaw.hh"
#include "MyContact.H"
#include "ERContact.H"
#include "ERHContact.H"
#include "CohContact.H"
#include "Particule.H"
#include "WriterXML.hh"

#include <assert.h>


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
map<string,int> Contact_BuilderFactory::m_materials;
map<string,bool> Contact_BuilderFactory::m_materialsForObstaclesOnly;
map<int,Contact_BuilderFactory::ContactFeatures>
	Contact_BuilderFactory::m_contactParametres;
map<int,ContactLaw*> Contact_BuilderFactory::m_contactLaws;

int Contact_BuilderFactory::m_value = 0x01;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie la loi de contact associee aux materiaux
ContactLaw* Contact_BuilderFactory::contactForceModel( const string &matA,
	const string &matB )
{
  return m_contactLaws[m_materials[matA] | m_materials[matB]];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Verifie que toutes les lois de contact existent pour l'ensemble
// des materiaux definis
bool Contact_BuilderFactory::checkContactLawsExist( string &matA, string &matB )
{
  map<string,int>::const_iterator imatA, imatB;
  map<string,bool>::const_iterator iobsA, iobsB;
  map<int,ContactFeatures>::const_iterator icp;
  bool ok = true, bA = true, bB = true ;

  for (imatA=m_materials.begin();imatA!=m_materials.end() && ok;imatA++,iobsA++)
  {
    matA = imatA->first;
    iobsA = m_materialsForObstaclesOnly.find(matA);
    if ( iobsA == m_materialsForObstaclesOnly.end() ) bA = true;
    else bA = iobsA->second;

    for (imatB=m_materials.begin();imatB!=m_materials.end() && ok;
    	imatB++,iobsB++)
    {
      matB = imatB->first;
      iobsB = m_materialsForObstaclesOnly.find(matB);
      if ( iobsB == m_materialsForObstaclesOnly.end() ) bB = true;
      else bB = iobsB->second;

      // Seuls les couples de materiaux ou au moins un des 2 est affecte a un
      // composant mobile sont testes
      // Si les 2 materiaux ne sont affectes qu'a des obstacles, aucune
      // interaction de contact n'est possible pour ce couple et il n'est pas
      // necessaire de definir une loi de contact
      if ( bA && bB ) ok = true ;
      else
      {
        int contactValue  = imatA->second | imatB->second;
        icp = m_contactParametres.find(contactValue);
        if ( icp == m_contactParametres.end() ) ok = false;
      }
    }
  }

  return ok;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition du type de contact utilise dans la simulation
void Contact_BuilderFactory::define( DOMNode* root )
{
  assert(root != NULL);

  DOMNodeList* allContacts = ReaderXML::getNodes(root);
  for (XMLSize_t i=0; i<allContacts->getLength(); i++) {
    DOMNode* contact = allContacts->item(i);
    string   name    = ReaderXML::getNodeName(contact);

    int contactValue;
    if (name == "Default")
    {
      // Contact Default
      contactValue = 0x00;

    }
    else
    {
      // Contact Pair
      DOMNode* material = ReaderXML::getNode(contact,"Material");
      string matA = ReaderXML::getNodeAttr_String(material, "materiauA");
      string matB = ReaderXML::getNodeAttr_String(material, "materiauB");
      contactValue = m_materials[matA] | m_materials[matB];
    }

    pair<Contact_BuilderFactory::ContactFeatures,ContactLaw*> forceLaw =
    	defineParameters(contact);
    m_contactParametres[contactValue] = forceLaw.first;
    m_contactLaws[contactValue] = forceLaw.second;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition des parametres associes au Contact
pair<Contact_BuilderFactory::ContactFeatures,ContactLaw*>
	Contact_BuilderFactory::defineParameters( DOMNode *root )
{
  pair<Contact_BuilderFactory::ContactFeatures,ContactLaw*>
  	forceLaw;
  DOMNodeList* allTags = ReaderXML::getNodes(root);
  // allTags->getLength() = 2 car
  // item(O) = Material
  // item(1) = MyContact, ERContact, CohContact or ERHContact
  for (XMLSize_t i=0; i<allTags->getLength(); i++)
  {
    DOMNode* contact = allTags->item(i);
    string   type    = ReaderXML::getNodeName(contact);
    if( type == "MyContact" )
    {
      forceLaw.first.name   = MYCONTACT;
      forceLaw.first.values = MyContact::defineParameters(contact);
      forceLaw.second = new MyContact(forceLaw.first.values);
    }
    else if( type == "ERContact" )
    {
      // cout << " TEMP ERContact" << endl;
      forceLaw.first.name   = ERCONTACT;
      forceLaw.first.values = ERContact::defineParameters(contact);
      forceLaw.second = new ERContact(forceLaw.first.values);
    }
    else if( type == "CohContact" )
    {
      // cout << " TEMP CohContact" << endl;
      forceLaw.first.name   = COHCONTACT;
      forceLaw.first.values = CohContact::defineParameters(contact);
      forceLaw.second = new CohContact(forceLaw.first.values);
    }
    else if( type == "ERHContact" )
    {
      //cout << " TEMP ERHContact" << endl;
      forceLaw.first.name   = ERHCONTACT;
      forceLaw.first.values = ERHContact::defineParameters(contact);
      forceLaw.second = new ERHContact(forceLaw.first.values);
    }
  }
  return forceLaw;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition des materiaux presents dans la simulation
void Contact_BuilderFactory::defineMaterial( const string &materiau,
  	const bool& matForObstacle )
{
  map<string,int>::const_iterator imat = m_materials.find(materiau);
  if ( imat == m_materials.end() )
  {
    m_materials[materiau] = Contact_BuilderFactory::m_value;
    m_materialsForObstaclesOnly[materiau] = matForObstacle;
    Contact_BuilderFactory::m_value <<= 1;
  }
  else
  {
    map<string,bool>::iterator iobs =
    	m_materialsForObstaclesOnly.find(materiau);
    iobs->second = !iobs->second ? false : matForObstacle;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture du fichier de persistance des contacts definis
// WRN : Le reload doit etre realise avant lecture de nouveau composant
void Contact_BuilderFactory::reload( istream &file )
{
  assert( Contact_BuilderFactory::m_value == 0x01 );

  // Lecture du nom du fichier de contacts
  string buffer, xmlFile;
  file >> buffer >> xmlFile >> buffer;

  // Suppression du directory (si il existe, ancien format de reload) et
  // ajout du bon directory de reload
  xmlFile = Grains_Exec::m_ReloadDirectory + "/"
  	+ Grains_Exec::extractFileName( xmlFile ) ;

  // Lecture du fichier de contacts
  DOMNode* root = ReaderXML::getRoot(xmlFile);

  if ( root )
  {
    DOMNode*     materiaux    = ReaderXML::getNode( root, "Materiaux" );
    DOMNodeList* allMateriaux = ReaderXML::getNodes( materiaux );
    for (XMLSize_t i=0; i<allMateriaux->getLength(); i++)
    {
      DOMNode* materiau = allMateriaux->item(i);
      int    value = ReaderXML::getNodeAttr_Int( materiau, "value" );
      string mat   = ReaderXML::getNodeValue_String( materiau );
      m_materials[mat] = value;
      Contact_BuilderFactory::m_value = value;
      Contact_BuilderFactory::m_value <<= 1;
    }

    DOMNode*     contacts    = ReaderXML::getNode( root, "Contacts" );
    DOMNodeList* allContacts = ReaderXML::getNodes( contacts );
    for (XMLSize_t i=0; i<allContacts->getLength(); i++)
    {
      DOMNode* contact = allContacts->item(i);
      int type  = ReaderXML::getNodeAttr_Int( contact, "type" );
      int value = ReaderXML::getNodeAttr_Int( contact, "value" );

      Contact_BuilderFactory::ContactFeatures parameters;
      parameters.name = (ContactType)type;

      DOMNode*     parametersNode = ReaderXML::getNode( contact, "Parameters" );
      DOMNodeList* allParameters  = ReaderXML::getNodes( parametersNode );
      for (XMLSize_t j=0; j<allParameters->getLength(); j++)
      {
        DOMNode* parameter = allParameters->item(j);
        string name  = ReaderXML::getNodeAttr_String ( parameter, "name" );
        double value_ = ReaderXML::getNodeValue_Double( parameter );
        parameters.values[name] = value_;
      }
      m_contactParametres[value] = parameters;
      if( parameters.name == MYCONTACT )
        m_contactLaws[value] = new MyContact( parameters.values );
      else if( parameters.name == ERCONTACT )
        m_contactLaws[value] = new ERContact( parameters.values );
      else if( parameters.name == COHCONTACT )
        m_contactLaws[value] = new CohContact( parameters.values );
      else if( parameters.name == ERHCONTACT )
        m_contactLaws[value] = new ERHContact( parameters.values );
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture du fichier de persistance des contacts definis
void Contact_BuilderFactory::save( ostream& file, const string& contactFile,
	const int& rank )
{
  string xmlFile( contactFile + "_MatContact.xml" );

  // Dans le fichier de reload des particules & obstacles, on indique le nom
  // du fichier de materiaux & contacts sans le chemin vers le directory
  string xmlFileNameInDir = Grains_Exec::extractFileName( xmlFile ) ;
  file << "<Materiaux> " << xmlFileNameInDir << " </Materiaux>\n";

  static int counter = 0 ;

  if ( rank == 0 && !counter )
  {
    DOMElement* root      = WriterXML::initialize("GRAINS");

    DOMElement* materiaux = WriterXML::createNode(root, "Materiaux");
    map<string,int>::const_iterator material;
    for (material=m_materials.begin(); material!=m_materials.end(); material++)
    {
      DOMElement* materialNode = WriterXML::createNode(materiaux, "Materiau");
      WriterXML::createNodeValue(materialNode, (*material).first);
      WriterXML::createNodeAttr (materialNode, "value", (*material).second);
    }

    DOMElement* contacts = WriterXML::createNode(root, "Contacts");
    map<int,Contact_BuilderFactory::ContactFeatures>::const_iterator contact;
    for (contact=m_contactParametres.begin();
    	contact!=m_contactParametres.end(); contact++)
    {
      DOMElement* contactNode = WriterXML::createNode(contacts, "Contact");
      WriterXML::createNodeAttr (contactNode, "type", (*contact).second.name);
      WriterXML::createNodeAttr (contactNode, "value", (*contact).first);

      DOMElement* valuesNode = WriterXML::createNode(contactNode, "Parameters");
      const map<string,double>& values = (*contact).second.values;
      map<string,double>::const_iterator value;
      for (value=values.begin(); value!=values.end(); value++)
      {
        DOMElement* valueNode = WriterXML::createNode(valuesNode, "Parameter");
        WriterXML::createNodeValue(valueNode, (*value).second);
        WriterXML::createNodeAttr (valueNode, "name", (*value).first);
      }
    }

    WriterXML::terminate(xmlFile);
    ++counter;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nom du mat�riau
string Contact_BuilderFactory::nomMaterial( int num )
{
  string res;
  map<string,int>::const_iterator im;

  for (im=m_materials.begin();im!=m_materials.end();im++)
    if ( im->second == num ) res=im->first;

  return res;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le num�ro du mat�riau
int Contact_BuilderFactory::numeroMaterial( const string &mat )
{
  return m_materials[mat];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destruit les lois de contact
void Contact_BuilderFactory::eraseAllContactLaws()
{
  map<int,ContactLaw*>::iterator icl;
  for (icl=m_contactLaws.begin();icl!=m_contactLaws.end();icl++)
    delete icl->second;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Verifie si un materiau ne s'applique qu'a un obstacle en cas de
// reload car cette info n'est pas conservee dans le fichier de reload
void Contact_BuilderFactory::set_materialsForObstaclesOnly_reload(
  	vector<Particule*> const* ParticuleClassesReference )
{
  map<string,int>::const_iterator imatA;
  vector<Particule*>::const_iterator particule;

  for (imatA=m_materials.begin();imatA!=m_materials.end();imatA++)
    m_materialsForObstaclesOnly[imatA->first] = true ;

  for (particule=ParticuleClassesReference->begin();
  	particule!=ParticuleClassesReference->end();particule++)
    m_materialsForObstaclesOnly[(*particule)->materiau()] = false ;
}
