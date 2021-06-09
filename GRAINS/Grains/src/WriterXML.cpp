#include "WriterXML.hh"

#include <assert.h>
#include <sstream>
using namespace std;

#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/AbstractDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/IOException.hpp>

#include <iostream>

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DOMDocument* WriterXML::m_document   = NULL;
DOMWriter*   WriterXML::m_serializer = NULL;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Initialisation du Writer
    @param root Nom du du noeud "root" 
    @return Le DOMElement racine du document */
DOMElement* WriterXML::initialize( const string &root )
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

  return m_document->getDocumentElement();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Ecriture du fichier xml & Liberation du Writer
    @param file Nom du fichier xml */
void WriterXML::terminate( const string &file )
{
  // Precondition : Nom fichier defini
  assert(file != "");

  try {
    XMLFormatTarget* target = new LocalFileFormatTarget( file.c_str() );
    
    // Ajouter par A. Wachs 26/08/14 pour obtenir des retours a la ligne
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




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Creation d'un noeud
    @param name Le nom du noeud construire
    @return Le noeud construit */
DOMElement* WriterXML::createNode( const string& name )
{
  XMLCh*      data = XMLString::transcode( name.c_str() );
  DOMElement* node = m_document->createElement( data );
  return node;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Creation & ajout d'un noeud
    @param root Le noeud racine
    @param name Le nom du noeud construire
    @return Le noeud construit */
DOMElement* WriterXML::createNode( DOMElement* root, const string& name )
{
  XMLCh*      data = XMLString::transcode( name.c_str() );
  DOMElement* node = m_document->createElement( data );
  root->appendChild( node );
  return node;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Ecriture de l'attribut du noeud
    @param root  Le noeud racine
    @param attr  Le nom de l'attribut
    @param value La valeur de l'attribut */
void WriterXML::createNodeAttr( DOMElement* root, const string &attr, 
	const string& value )
{
  XMLCh* dataAttr  =  XMLString::transcode( attr.c_str() );
  XMLCh* dataValue =  XMLString::transcode( value.c_str() );
  root->setAttribute(dataAttr, dataValue);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Ecriture de l'attribut du noeud
    @param root  Le noeud racine
    @param attr  Le nom de l'attribut
    @param value La valeur de l'attribut */
void WriterXML::createNodeAttr( DOMElement* root, const string &attr, 
	double value )
{
  ostringstream strValue;
  strValue << value;
  WriterXML::createNodeAttr( root, attr, strValue.str() );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Ecriture de la valeur du noeud
    @param root  Le noeud racine
    @param value La valeur du noeud */
void WriterXML::createNodeValue( DOMElement* root, const string &value )
{
  XMLCh*   data = XMLString::transcode( value.c_str() );
  DOMText* node = m_document->createTextNode( data );
  root->appendChild( node );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Ecriture de la valeur du noeud
    @param root  Le noeud racine
    @param value La valeur du noeud */
void WriterXML::createNodeValue( DOMElement* root, const double &value )
{
  ostringstream strValue;
  strValue << value;
  WriterXML::createNodeValue( root, strValue.str() );
}

