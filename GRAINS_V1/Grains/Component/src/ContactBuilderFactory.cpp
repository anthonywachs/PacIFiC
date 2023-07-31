#include "GrainsExec.hh"
#include "ContactBuilderFactory.hh"
#include "ContactForceModel.hh"
#include "HODCContactForceModel.hh"
#include "MemoryContactForceModel.hh"
#include "Particle.hh"
#include "WriterXML.hh"
#include <assert.h>

map<string,int> ContactBuilderFactory::m_materials;
map<string,bool> ContactBuilderFactory::m_materialsForObstaclesOnly;
map<int,ContactBuilderFactory::ContactFeatures>
	ContactBuilderFactory::m_contactParametres;
map<int,ContactForceModel*> ContactBuilderFactory::m_contactForceModels;
int ContactBuilderFactory::m_value = 0x01;


// ----------------------------------------------------------------------------
// Returns the contact associated to a couple of materials
ContactForceModel* ContactBuilderFactory::contactForceModel( const string &matA,
	const string &matB )
{
  return ( m_contactForceModels[ m_materials[matA] | m_materials[matB] ] );
}




// ----------------------------------------------------------------------------
// Checks that all couples of materials in which one of the
/// materials is assigned to a particle has an associated contact force model.
// If not, the two material names are copied in matA and matB
bool ContactBuilderFactory::checkContactForceModelsExist( string& matA,
	string& matB )
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

      // Only pairs of materials in which one of the materials is assigned
      // to a freely moving component are checked
      // If both materials are assigned to obstacles (stationary or with a
      // prescribed motion), there is no need to define a contact force model
      if ( bA && bB ) ok = true ;
      else
      {
        int contactValue  = imatA->second | imatB->second;
        icp = m_contactParametres.find(contactValue);
        if ( icp == m_contactParametres.end() ) ok = false;
      }
    }
  }

  return ( ok );
}




// ----------------------------------------------------------------------------
// Creates the table (couple of materials)-(contact force model) from an XML
// node
void ContactBuilderFactory::define( DOMNode* root )
{
  assert(root != NULL);

  DOMNodeList* allContacts = ReaderXML::getNodes(root);
  for (XMLSize_t i=0; i<allContacts->getLength(); i++) {
    DOMNode* contact = allContacts->item(i);
    string   name    = ReaderXML::getNodeName(contact);

    int contactValue;
    if ( name == "Default" )
    {
      // Contact Default
      contactValue = 0x00;
    }
    else
    {
      // Contact Pair
      DOMNode* material = ReaderXML::getNode(contact,"Material");
      string matA = ReaderXML::getNodeAttr_String(material, "materialA");
      string matB = ReaderXML::getNodeAttr_String(material, "materialB");
      contactValue = m_materials[matA] | m_materials[matB];
    }

    pair<ContactBuilderFactory::ContactFeatures,ContactForceModel*> forceLaw =
    	defineParameters(contact);
    m_contactParametres[contactValue] = forceLaw.first;
    m_contactForceModels[contactValue] = forceLaw.second;
  }
}




// ----------------------------------------------------------------------------
// Reads and returns the contact parameters of the contact force model from
// an XML node
pair<ContactBuilderFactory::ContactFeatures,ContactForceModel*>
	ContactBuilderFactory::defineParameters( DOMNode* root )
{
  pair<ContactBuilderFactory::ContactFeatures,ContactForceModel*>
  	forceLaw;
  DOMNodeList* allTags = ReaderXML::getNodes(root);
  // Remark: allTags->getLength() = 2 because
  // item(O) = Material
  // item(1) = HODCContactForceModel or other contact force model
  for (XMLSize_t i=0; i<allTags->getLength(); i++)
  {
    DOMNode* contact = allTags->item(i);
    string   type    = ReaderXML::getNodeName(contact);
    if ( type == "HODC" )
    {
      forceLaw.first.name   = HODC;
      forceLaw.first.values = HODCContactForceModel::defineParameters(contact);
      forceLaw.second = new HODCContactForceModel(forceLaw.first.values);
    }
    else if ( type == "Memory" )
    {
			forceLaw.first.name   = Memory;
      forceLaw.first.values =
				MemoryContactForceModel::defineParameters(contact);
      forceLaw.second = new MemoryContactForceModel(forceLaw.first.values);
    }
  }

  return ( forceLaw );
}




// ----------------------------------------------------------------------------
// Adds a material to the simulation
void ContactBuilderFactory::defineMaterial( string const& material,
	bool const& matForObstacle )
{
  map<string,int>::const_iterator imat = m_materials.find(material);
  if ( imat == m_materials.end() )
  {
    m_materials[material] = ContactBuilderFactory::m_value;
    m_materialsForObstaclesOnly[material] = matForObstacle;
    ContactBuilderFactory::m_value <<= 1;
  }
  else
  {
    map<string,bool>::iterator iobs =
    	m_materialsForObstaclesOnly.find(material);
    iobs->second = !iobs->second ? false : matForObstacle;
  }
}




// ----------------------------------------------------------------------------
// Reloads and recreates contact force models and tables
void ContactBuilderFactory::reload( istream& file )
{
  assert( ContactBuilderFactory::m_value == 0x01 );

  // Read the contact file name
  string buffer, xmlFile;
  file >> buffer >> xmlFile >> buffer;

  // Delete the directory path if it exists and adds the correct directory path
  xmlFile = GrainsExec::m_ReloadDirectory + "/"
  	+ GrainsExec::extractFileName( xmlFile ) ;

  // Read the contact file
  DOMNode* root = ReaderXML::getRoot( xmlFile );

  if ( root )
  {
    DOMNode* materiaux = ReaderXML::getNode( root, "Materials" );
    DOMNodeList* allMateriaux = ReaderXML::getNodes( materiaux );
    for (XMLSize_t i=0; i<allMateriaux->getLength(); i++)
    {
      DOMNode* materiau = allMateriaux->item(i);
      int value = ReaderXML::getNodeAttr_Int( materiau, "value" );
      string mat = ReaderXML::getNodeValue_String( materiau );
      m_materials[mat] = value;
      ContactBuilderFactory::m_value = value;
      ContactBuilderFactory::m_value <<= 1;
    }

    DOMNode* contacts = ReaderXML::getNode( root, "ContactForceModels" );
    DOMNodeList* allContacts = ReaderXML::getNodes( contacts );
    for (XMLSize_t i=0; i<allContacts->getLength(); i++)
    {
      DOMNode* contact = allContacts->item(i);
      int type = ReaderXML::getNodeAttr_Int( contact, "type" );
      int value = ReaderXML::getNodeAttr_Int( contact, "value" );

      ContactBuilderFactory::ContactFeatures parameters;
      parameters.name = (ContactType)type;

      DOMNode* parametersNode = ReaderXML::getNode( contact, "Parameters" );
      DOMNodeList* allParameters = ReaderXML::getNodes( parametersNode );
      for (XMLSize_t j=0; j<allParameters->getLength(); j++)
      {
        DOMNode* parameter = allParameters->item(j);
        string name = ReaderXML::getNodeAttr_String ( parameter, "name" );
        double value_ = ReaderXML::getNodeValue_Double( parameter );
        parameters.values[name] = value_;
      }
      m_contactParametres[value] = parameters;
      if ( parameters.name == HODC )
        m_contactForceModels[value] = new HODCContactForceModel(
				parameters.values );
      else if( parameters.name == Memory )
        m_contactForceModels[value] = new MemoryContactForceModel( parameters.values );
    }
  }
}




// ----------------------------------------------------------------------------
// Writes the whole contact force model data in a stream
void ContactBuilderFactory::save( ostream& file, const string& contactFile,
	const int& rank )
{
  string xmlFile( contactFile + "_MatContact.xml" );

  // In the particle & obstacle reload file, we write the name of the
  // contact file without the directory path
  string xmlFileNameInDir = GrainsExec::extractFileName( xmlFile ) ;
  file << "<ContactForceModels> " << xmlFileNameInDir
  	<< " </ContactForceModels>" << endl;

  static int counter = 0 ;

  if ( rank == 0 && !counter )
  {
    DOMElement* root = WriterXML::initialize( "GRAINS" );

    DOMElement* materiaux = WriterXML::createNode( root, "Materials" );
    map<string,int>::const_iterator material;
    for (material=m_materials.begin(); material!=m_materials.end(); material++)
    {
      DOMElement* materialNode = WriterXML::createNode( materiaux, "Material" );
      WriterXML::createNodeValue( materialNode, (*material).first );
      WriterXML::createNodeAttr( materialNode, "value", (*material).second );
    }

    DOMElement* contacts = WriterXML::createNode( root, "ContactForceModels" );
    map<int,ContactBuilderFactory::ContactFeatures>::const_iterator contact;
    for (contact=m_contactParametres.begin();
    	contact!=m_contactParametres.end(); contact++)
    {
      DOMElement* contactNode = WriterXML::createNode( contacts,
      	"ContactForceModel" );
      WriterXML::createNodeAttr( contactNode, "type", (*contact).second.name );
      WriterXML::createNodeAttr( contactNode, "value", (*contact).first );

      DOMElement* valuesNode = WriterXML::createNode( contactNode,
      	"Parameters" );
      const map<string,double>& values = (*contact).second.values;
      map<string,double>::const_iterator value;
      for (value=values.begin(); value!=values.end(); value++)
      {
        DOMElement* valueNode = WriterXML::createNode( valuesNode,
		"Parameter" );
        WriterXML::createNodeValue( valueNode, (*value).second );
        WriterXML::createNodeAttr( valueNode, "name", (*value).first );
      }
    }

    WriterXML::terminate( xmlFile );
    ++counter;
  }
}




// ----------------------------------------------------------------------------
// Returns the name of the material given its ID number
string ContactBuilderFactory::materialName( int num )
{
  string res;
  map<string,int>::const_iterator im;

  for (im=m_materials.begin();im!=m_materials.end();im++)
    if ( im->second == num ) res=im->first;

  return ( res );
}




// ----------------------------------------------------------------------------
// Returns the ID number of the material given its name
int ContactBuilderFactory::materialIDnumber( const string &mat )
{
  return ( m_materials[mat] );
}




// ----------------------------------------------------------------------------
// Deletes all contact force models
void ContactBuilderFactory::eraseAllContactForceModels()
{
  map<int,ContactForceModel*>::iterator icl;
  for (icl=m_contactForceModels.begin();icl!=m_contactForceModels.end();icl++)
    delete icl->second;
}




// ----------------------------------------------------------------------------
// Makes sure that materials for obstacles only are identified after
// a reload as these data are not stored in the reload file
void ContactBuilderFactory::set_materialsForObstaclesOnly_reload(
  	vector<Particle*> const* referenceParticles )
{
  map<string,int>::const_iterator imatA;
  vector<Particle*>::const_iterator particle;

  for (imatA=m_materials.begin();imatA!=m_materials.end();imatA++)
    m_materialsForObstaclesOnly[imatA->first] = true ;

  for (particle=referenceParticles->begin();
  	particle!=referenceParticles->end();particle++)
    m_materialsForObstaclesOnly[(*particle)->getMaterial()] = false ;
}
