#include "WriterXML.hh"
#include <assert.h>
#include <sstream>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/AbstractDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/IOException.hpp>
#include <iostream>
using namespace std;


DOMDocument* WriterXML::m_document   = NULL;
DOMWriter*   WriterXML::m_serializer = NULL;


// ----------------------------------------------------------------------------
// Initializes the writer and returns the document root node
DOMElement* WriterXML::initialize( string const& root )
{
  // Precondition : Nom racine defini
  assert(root != "");

  // Description du serializer
  XMLPlatformUtils::Initialize();

  static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
  DOMImplementation *impl  = 
    DOMImplementationRegistry::getDOMImplementation(gLS);
  
  XMLCh* name  = XMLString::transcode(root.c_str());
  m_document   = impl->createDocument(0, name, 0);
  m_serializer = impl->createDOMWriter();

  return ( m_document->getDocumentElement() );
}




// ----------------------------------------------------------------------------
// Flushes to the xml file and frees the writer
void WriterXML::terminate( string const& file )
{
  // Precondition : Nom fichier defini
  assert(file != "");

  try {
    XMLFormatTarget* target = new LocalFileFormatTarget( file.c_str() );
    
    // A. Wachs 26/08/14 pour obtenir des retours a la ligne
//    m_serializer->setNewLine(XMLString::transcode("\n"));
    m_serializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, 
    	true );
	
    m_serializer->writeNode( target, *m_document );
    target->flush();
    delete target;

    delete m_serializer;
    XMLPlatformUtils::Terminate();

  } catch (IOException &e) {
    std::cout << "Impossible writing ouput file " << file << endl;
    std::cout << "Execution continue on error" << endl;
  }

}




// ----------------------------------------------------------------------------
// Creates a node from its name and returns the node    
DOMElement* WriterXML::createNode( string const& name )
{
  XMLCh*      data = XMLString::transcode( name.c_str() );
  DOMElement* node = m_document->createElement( data );
  return ( node );
}




// ----------------------------------------------------------------------------
// Creates a node from its name and its root and returns the node   
DOMElement* WriterXML::createNode( DOMElement* root, string const& name )
{
  XMLCh*      data = XMLString::transcode( name.c_str() );
  DOMElement* node = m_document->createElement( data );
  root->appendChild( node );
  return ( node );
}




// ----------------------------------------------------------------------------
// Creates a node attribute with a string value to be later written
void WriterXML::createNodeAttr( DOMElement* root, string const& attr, 
	 string const& value )
{
  XMLCh* dataAttr  =  XMLString::transcode( attr.c_str() );
  XMLCh* dataValue =  XMLString::transcode( value.c_str() );
  root->setAttribute( dataAttr, dataValue );
}




// ----------------------------------------------------------------------------
// Creates a node attribute with a scalar value to be later written
void WriterXML::createNodeAttr( DOMElement* root, string const& attr, 
	double value )
{
  ostringstream strValue;
  strValue << value;
  WriterXML::createNodeAttr( root, attr, strValue.str() );
}




// ----------------------------------------------------------------------------
// Creates a string node value to be later written
void WriterXML::createNodeValue( DOMElement* root, string const& value )
{
  XMLCh*   data = XMLString::transcode( value.c_str() );
  DOMText* node = m_document->createTextNode( data );
  root->appendChild( node );
}




// ----------------------------------------------------------------------------
// Creates a scalar node value to be later written
void WriterXML::createNodeValue( DOMElement* root, double const& value )
{
  ostringstream strValue;
  strValue << value;
  WriterXML::createNodeValue( root, strValue.str() );
}

