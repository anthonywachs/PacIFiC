#ifndef READERXML_H_
#define READERXML_H_

#include <xercesc/dom/DOMBuilder.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
XERCES_CPP_NAMESPACE_USE

#include <string>
using namespace std;


/** @brief Utilities class : Reader XML
  
    @author G. Ferrer - Grains Project - Institut Francais du Petrole - 2005 */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ReaderXML
{
 public:
  /** @name Static Methods */
  //@{
  /** @brief Initialisation du Reader */
  static void initialize();
  
  /** @brief Liberation du Reader */
  static void terminate();

  /** @brief Recherche du noeud "name" dans le document
  @param root Noeud racine du document
  @param name Nom du noeud
  @return Le noeud demande */
  static DOMNode* getNode(DOMElement* root, const string &name);
  
  /** @brief Recherche du noeud "name" a partir du Noeud
  @param root Noeud reference
  @param name Nom du noeud
  @return Le noeud demande */
  static DOMNode* getNode(DOMNode* root, const string &name);

  /** @brief Valeur de l'attribut "name" du noeud : "<root name="xxx">"
  @param root Noeud de reference
  @param name Nom de l'attribut
  @return Valeur demandee */
  static double getNodeAttr_Double(DOMNode* root, const string &name);
  static bool hasNodeAttr_Double(DOMNode* root, const string &name);
  
  /** @brief Valeur de l'attribut "name" du noeud : "<root name="xxx">"
  @param root Noeud de reference
  @param name Nom de l'attribut
  @return Valeur demandee */
  static int getNodeAttr_Int(DOMNode* root, const string &name);

  /** @brief Valeur de l'attribut "name" du noeud : "<root name="xxx">"
  @param root Noeud de reference
  @param name Nom de l'attribut
  @return Valeur demandee */
  static string getNodeAttr_String(DOMNode* root, const string &name);

  /** @brief l'attribut "name" du noeud : "<root name="xxx">" existe-t-il?
  @param root Noeud de reference
  @param name Nom de l'attribut */
  static bool hasNodeAttr_String(DOMNode* root, const string &name);
  
  /** @brief Nom du noeud.
  @param root Noeud desire.
  @return Nom */
  static string getNodeName(const DOMNode* root);

  /** @brief Noeud suivant
  @param root Noeud de reference
  @return Le noeud interne (NULL si pas de noeud) */
  static DOMNode* getNodeNext(DOMNode* root);

  /** @brief Valeur du noeud : "<root>xxx</root>"
  @param root Noeud 
  @return Valeur du noeud */
  static double getNodeValue_Double(DOMNode* root);

  /** @brief Valeur du noeud : "<root>xxx</root>"
  @param root Noeud 
  @return Valeur du noeud */
  static string getNodeValue_String(DOMNode* root);

  /** @brief Recherche la liste des noeuds "name" dans le document
  @param root Noeud racine du document
  @param name Nom du noeud
  @return Les noeuds demandes */
  static DOMNodeList* getNodes(DOMElement* root, const string &name);
  
  /** @brief Recherche de la liste des noeuds sous le noeud root
  @param root Noeud de reference
  @return Les noeuds presents */
  static DOMNodeList* getNodes(DOMNode* root);

  /** @brief Noeud principal du document (root)
  @param xmlFile Le fichier a parser
  @return Le noeud root */
  static DOMElement* getRoot(const string& xmlFile);
  //@}

 private:
  /** @name Parameters */
  //@{
  static DOMBuilder *m_parser; /**< Parser */
  //@}
  
  /** @name Constructors & Destructor */
  //@{
  /** @brief Static class : No constructor authorized */
  ReaderXML();
  
  /** @brief Static class : No destructor */
  ~ReaderXML();
  //@}
};

#endif 
