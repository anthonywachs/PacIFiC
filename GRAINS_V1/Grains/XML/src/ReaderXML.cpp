#include "ReaderXML.hh"
#include <iostream>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/parsers/AbstractDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
using namespace std;


DOMBuilder* ReaderXML::m_parser = NULL;


// ----------------------------------------------------------------------------
// Initializes the reader
void ReaderXML::initialize()
{
  // Description du parser
  XMLPlatformUtils::Initialize();

  static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
  DOMImplementation *impl   = 
    DOMImplementationRegistry::getDOMImplementation(gLS);
  m_parser = 
    impl->createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);

  bool doNamespaces=false;
  bool doSchema=false;
  bool schemaFullChecking=false;
  m_parser->setFeature(XMLUni::fgDOMNamespaces,            doNamespaces);
  m_parser->setFeature(XMLUni::fgXercesSchema,             doSchema);
  m_parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);

  m_parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
}




// ----------------------------------------------------------------------------
// Frees the reader
void ReaderXML::terminate()
{
  XMLPlatformUtils::Terminate();
}




// ----------------------------------------------------------------------------
// Returns a node from its name and the root node
DOMNode* ReaderXML::getNode( DOMElement* root, string const& name )
{
  return ( root->getElementsByTagName(XMLString::transcode(name.c_str()))
	->item(0) );
}




// ----------------------------------------------------------------------------
// Returns a node from its name and the root node 
DOMNode* ReaderXML::getNode( DOMNode* root, string const& name )
{
  DOMNode* node = NULL;
  DOMNodeList* nodes = ReaderXML::getNodes(root);
  for (XMLSize_t i=0; i<nodes->getLength() && node==NULL; i++) {
    if (name == ReaderXML::getNodeName(nodes->item(i)))
	node = nodes->item(i);
  }  
  return ( node );
}




// ----------------------------------------------------------------------------
// Returns the scalar value xxx of the attribute of a node "<root name="xxx">"
double ReaderXML::getNodeAttr_Double( DOMNode* root, string const& name )
{
  DOMNamedNodeMap *nodeValues = root->getAttributes();
  DOMNode *value = nodeValues->getNamedItem(XMLString::transcode(name.c_str()));
  return ( atof( XMLString::transcode(value->getNodeValue()) ) );
}





// ----------------------------------------------------------------------------
// Returns the integer value xxx of the attribute of a node "<root name="xxx">"
int ReaderXML::getNodeAttr_Int(DOMNode* root, const string &name)
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XMLCh*           attrName   = XMLString::transcode(name.c_str());
  DOMNode* value = nodeValues->getNamedItem(attrName);
  return ( atoi( XMLString::transcode(value->getNodeValue()) ) );
}




// ----------------------------------------------------------------------------
// Returns the string value xxx of the attribute of a node "<root name="xxx">"
string ReaderXML::getNodeAttr_String( DOMNode* root, string const& name )
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XMLCh*           attrName   = XMLString::transcode(name.c_str());
  DOMNode* node = nodeValues->getNamedItem(attrName);
  const XMLCh* value = node->getNodeValue();
  return ( XMLString::transcode(value) );
}



// ----------------------------------------------------------------------------
// Returns whether the node root has an attribute named name "<node name="xxx">"
bool ReaderXML::hasNodeAttr( DOMNode* root, string const& name )
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XMLCh*           attrName   = XMLString::transcode(name.c_str());
  DOMNode* node = nodeValues->getNamedItem(attrName);
  return ( node != NULL );
}



// ----------------------------------------------------------------------------
// Returns the name of a node 
string ReaderXML::getNodeName( DOMNode const* root )
{
  const XMLCh* name = root->getNodeName();
  char* nodeName = XMLString::transcode(name);
  return ( string(nodeName) );
}




// ----------------------------------------------------------------------------
// Returns the next node wrt to a node
DOMNode* ReaderXML::getNodeNext( DOMNode* root )
{
  DOMNodeList* nodes = root->getChildNodes();
  return ( nodes->item(1) );
}




// ----------------------------------------------------------------------------
// Returns the scalar value xxx of the node "<root>xxx</root>"
double ReaderXML::getNodeValue_Double( DOMNode* root )
{
  return ( 
  	atof( XMLString::transcode(root->getFirstChild()->getNodeValue()) ) );
}




// ----------------------------------------------------------------------------
// Returns the string value xxx of the node "<root>xxx</root>"
string ReaderXML::getNodeValue_String( DOMNode* root )
{
  return ( XMLString::transcode(root->getFirstChild()->getNodeValue()) );
}




// ----------------------------------------------------------------------------
// Returns the list of nodes in a node
DOMNodeList* ReaderXML::getNodes( DOMElement* root, string const& name )
{
  return ( root->getElementsByTagName(XMLString::transcode(name.c_str())) );
}




// ----------------------------------------------------------------------------
// Returns the list of nodes in a node
DOMNodeList* ReaderXML::getNodes( DOMNode* root )
{
  DOMNode* allNodes  = root->cloneNode(true);
  DOMNodeList* nodes = allNodes->getChildNodes();
  for (XMLSize_t i=0; i<nodes->getLength(); i++) {
    DOMNode* node =  nodes->item(i);
    if (node->getNodeType() != 1)
      allNodes->removeChild(node);
  }
  return ( allNodes->getChildNodes() );
}




// ----------------------------------------------------------------------------
// Noeud principal du document (root)
DOMElement* ReaderXML::getRoot( string const& xmlFile )
{
  DOMDocument *doc=NULL;
  DOMElement  *root=NULL;
  try {
    doc  = m_parser->parseURI(xmlFile.c_str());
    root = doc->getDocumentElement();
  } catch (const DOMException& e) {
    XERCES_STD_QUALIFIER cerr << "XML exception " << e.code 
			      << XERCES_STD_QUALIFIER endl
			      << XMLString::transcode(e.msg)
			      << XERCES_STD_QUALIFIER endl;
  }
  return ( root );
}
