#include "ReaderXML.hh"

#include <iostream>
using namespace std;

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


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DOMBuilder* ReaderXML::m_parser = NULL;




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialisation du Reader
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




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Liberation du Reader
void ReaderXML::terminate()
{
  XMLPlatformUtils::Terminate();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recherche du noeud "name" dans le document
DOMNode* ReaderXML::getNode(DOMElement* root, const string &name)
{
    return root->getElementsByTagName(XMLString::transcode(name.c_str()))
    	->item(0);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recherche du noeud "name" a partir du Noeud
DOMNode* ReaderXML::getNode(DOMNode* root, const string &name)
{
  DOMNode* node = NULL;
  DOMNodeList* nodes = ReaderXML::getNodes(root);
  for (XMLSize_t i=0; i<nodes->getLength() && node==NULL; i++) {
    if (name == ReaderXML::getNodeName(nodes->item(i)))
	node = nodes->item(i);
  }
  return node;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Valeur de l'attribut "name" du noeud
double ReaderXML::getNodeAttr_Double(DOMNode* root, const string &name)
{
  DOMNamedNodeMap *nodeValues = root->getAttributes();
  DOMNode *value = nodeValues->getNamedItem(XMLString::transcode(name.c_str()));
  return atof( XMLString::transcode(value->getNodeValue()) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// D. RAKOTINIRINA - Dec 2013 - Creation
// L'attribut "name" existe-t-il?
bool ReaderXML::hasNodeAttr_Double(DOMNode* root, const string &name)
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XMLCh*           attrName   = XMLString::transcode(name.c_str());
  DOMNode* node = nodeValues->getNamedItem(attrName);
  return ( node != NULL );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Valeur de l'attribut "name" du noeud
int ReaderXML::getNodeAttr_Int(DOMNode* root, const string &name)
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XMLCh*           attrName   = XMLString::transcode(name.c_str());
  DOMNode* value = nodeValues->getNamedItem(attrName);
  return atoi( XMLString::transcode(value->getNodeValue()) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Valeur de l'attribut "name" du noeud
string ReaderXML::getNodeAttr_String(DOMNode* root, const string &name)
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XMLCh*           attrName   = XMLString::transcode(name.c_str());
  DOMNode* node = nodeValues->getNamedItem(attrName);
  const XMLCh* value = node->getNodeValue();
  return XMLString::transcode(value);
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// D. RAKOTINIRINA - Dec 2013 - Creation
// L'attribut "name" existe-t-il?
bool ReaderXML::hasNodeAttr_String(DOMNode* root, const string &name)
{
  DOMNamedNodeMap* nodeValues = root->getAttributes();
  XMLCh*           attrName   = XMLString::transcode(name.c_str());
  DOMNode* node = nodeValues->getNamedItem(attrName);
  return ( node != NULL );
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nom du noeud.
string ReaderXML::getNodeName(const DOMNode* root)
{
  const XMLCh* name = root->getNodeName();
  char* nodeName = XMLString::transcode(name);
  return string(nodeName);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @name Noeud suivant
    @param root Noeud de reference
    @return Le noeud interne (NULL si pas de noeud) */
DOMNode* ReaderXML::getNodeNext(DOMNode* root)
{
  DOMNodeList* nodes = root->getChildNodes();
  return nodes->item(1);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Valeur du noeud
double ReaderXML::getNodeValue_Double(DOMNode* root)
{
  return atof( XMLString::transcode(root->getFirstChild()->getNodeValue()) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Valeur du noeud
string ReaderXML::getNodeValue_String(DOMNode* root)
{
  return XMLString::transcode(root->getFirstChild()->getNodeValue());
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recherche la liste des noeuds "name" dans le document
DOMNodeList* ReaderXML::getNodes(DOMElement* root, const string &name)
{
  return root->getElementsByTagName(XMLString::transcode(name.c_str()));
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recherche de la liste des noeuds sous le noeud root
DOMNodeList* ReaderXML::getNodes(DOMNode* root)
{
  DOMNode* allNodes  = root->cloneNode(true);
  DOMNodeList* nodes = allNodes->getChildNodes();
  for (XMLSize_t i=0; i<nodes->getLength(); i++) {
    DOMNode* node =  nodes->item(i);
    if (node->getNodeType() != 1)
      allNodes->removeChild(node);
  }
  return allNodes->getChildNodes();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Noeud principal du document (root)
DOMElement* ReaderXML::getRoot(const string& xmlFile)
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
  return root;
}
