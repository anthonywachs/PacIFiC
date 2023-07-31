#ifndef _WRITERXML_HH_
#define _WRITERXML_HH_

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMWriter.hpp>
XERCES_CPP_NAMESPACE_USE
#include <Basic.hh>
#include <string>
using namespace std;


/** @brief The class WriterXML.

    Utilities to write a file in an XML format with Xerces.
    
    @author GRAINS Project - IFP - 2007 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class WriterXML
{
  public:
    /** @name Static methods */
    //@{
    /** @brief Initializes the writer and returns the document root node 
    @param root name of the root node */
    static DOMElement* initialize( string const& root );

    /** @brief Flushes to the xml file and frees the writer
    @param file xml file name */
    static void terminate( string const& file );

    /** @brief Creates a node from its name and returns the node  
    @param name node name */
    static DOMElement* createNode( string const& name );

    /** @brief Creates a node from its name and its root and returns the node 
    @param root root node
    @param name node name */
    static DOMElement* createNode( DOMElement* root, string const& name );

    /** @brief Creates a node attribute with a string value to be later written
    @param root root node
    @param attr attribute name
    @param value string value of the attribute */
    static void createNodeAttr( DOMElement* root, string const& attr, 
	 string const& value );

    /** @brief Creates a node attribute with a scalar value to be later written
    @param root root node
    @param attr attribute name
    @param value scalar value of the attribute */
    static void createNodeAttr( DOMElement* root, string const& attr, 
	double value );

    /** @brief Creates a string node value to be later written
    @param root node
    @param value string value */
    static void createNodeValue( DOMElement* root, string const& value );

    /** @brief Creates a scalar node value to be later written
    @param root node
    @param value scalar value */
    static void createNodeValue( DOMElement* root, double const& value );
    //@}


  private:
    /** @name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    WriterXML();
  
    /** @brief Destructor  (forbidden) */
    ~WriterXML();
    //@}  


    /** @name Parameters */
    //@{
    static DOMWriter* m_serializer; /**< Serializer */
    static DOMDocument* m_document; /**< Document root */
    //@}
};

#endif
