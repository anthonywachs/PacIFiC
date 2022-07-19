#ifndef _READERXML_HH_
#define _READERXML_HH_

#include <xercesc/dom/DOMBuilder.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
XERCES_CPP_NAMESPACE_USE
#include <Basic.hh>
#include <string>
using namespace std;


/** @brief The class WriterXML.

    Utilities to write a file in an XML format with Xerces.
    
    @author GRAINS Project - IFP - 2007 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ReaderXML
{
   public:
    /** @name Static Methods */
    //@{
    /** @brief Initializes the reader */
    static void initialize();
  
    /** @brief Frees the reader */
    static void terminate();

    /** @brief Returns a node from its name and the root node
    @param root root node
    @param name node name */
    static DOMNode* getNode( DOMElement* root, string const& name );
  
    /** @brief Returns a node from its name and the root node 
    @param root root node
    @param name node name */
    static DOMNode* getNode( DOMNode* root, string const& name );

    /** @brief Returns the scalar value xxx of the attribute of a node 
    "<root name="xxx">"
    @param root node
    @param name attribute name */
    static double getNodeAttr_Double( DOMNode* root, string const& name );
  
    /** @brief Returns the integer value xxx of the attribute of a node 
    "<root name="xxx">"
    @param root node
    @param name attribute name */ 
    static int getNodeAttr_Int( DOMNode* root, string const& name );

    /** @brief Returns the string value xxx of the attribute of a node 
    "<root name="xxx">"
    @param root node
    @param name attribute name */ 
    static string getNodeAttr_String( DOMNode* root, string const& name );

    /** @brief Returns whether the node root has an attribute named name
    "<node name="xxx">"
    @param root node
    @param name attribute name */    
    static bool hasNodeAttr( DOMNode* root, string const& name );
  
    /** @brief Returns the name of a node 
    @param root node */
    static string getNodeName( DOMNode const* root );

    /** @brief Returns the next node wrt to a node
    @param root previous node */
    static DOMNode* getNodeNext( DOMNode* root );

    /** @brief Returns the scalar value xxx of the node "<root>xxx</root>"
    @param root node */
    static double getNodeValue_Double( DOMNode* root );

    /** @brief Returns the string value xxx of the node "<root>xxx</root>"
    @param root node */
    static string getNodeValue_String( DOMNode* root );

    /** @brief Returns the list of nodes in a node
    @param root root node
    @param name generic name of each node in the list of nodes */
    static DOMNodeList* getNodes( DOMElement* root, string const& name );
  
    /** @brief Returns the list of nodes in a node
    @param root root node */
    static DOMNodeList* getNodes( DOMNode* root );

    /** @brief Returns the main root node of the document
    @param xmlFile xml file name */
    static DOMElement* getRoot( string const& xmlFile );
    //@}

   private:
    /** @name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    ReaderXML();
  
    /** @brief Destructor  (forbidden) */
    ~ReaderXML();
    //@}  

    /** @name Parameters */
    //@{
    static DOMBuilder* m_parser; /**< Parser */
    //@}  
};

#endif 
