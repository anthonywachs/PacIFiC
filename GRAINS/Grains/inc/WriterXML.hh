#ifndef WRITERXML_H_
#define WRITERXML_H_

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMWriter.hpp>
XERCES_CPP_NAMESPACE_USE

#include <string>
using namespace std;


/** @brief Utilities calls : Writer XML

@author GRAINS Project - IFP - 2007 - Creation */
class WriterXML
{
public:
  /** @name Static Methods */
  //@{
  /** @brief Initialisation du Writer
  @param root Nom du du noeud "root" 
  @return Le DOMElement racine du document */
  static DOMElement* initialize(const string &root);

  /** @brief Ecriture du fichier xml & Liberation du Writer
  @param file Nom du fichier xml */
  static void terminate(const string &file);

  /** @brief Creation d'un noeud
  @param name Le nom du noeud construire
  @return Le noeud construit */
  static DOMElement* createNode(const string& name);

  /** @brief Creation & ajout d'un noeud
  @param root Le noeud racine
  @param name Le nom du noeud construire
  @return Le noeud construit */
  static DOMElement* createNode(DOMElement* root, const string& name);

  /** @brief Ecriture de l'attribut du noeud
  @param root  Le noeud racine
  @param attr  Le nom de l'attribut
  @param value La valeur de l'attribut */
  static void createNodeAttr(DOMElement* root, const string &attr, 
			     const string& value);

  /** @brief Ecriture de l'attribut du noeud
  @param root  Le noeud racine
  @param attr  Le nom de l'attribut
  @param value La valeur de l'attribut */
  static void createNodeAttr(DOMElement* root, const string &attr, 
			     double value);

  /** @brief Ecriture de la valeur du noeud
  @param root  Le noeud racine
  @param value La valeur du noeud */
  static void createNodeValue(DOMElement* root, const string &value);

  /** @brief Ecriture de la valeur du noeud
  @param root  Le noeud racine
  @param value La valeur du noeud */
  static void createNodeValue(DOMElement* root, const double &value);
  //@}


private:
  /** @name Constructors & Destructor */
  //@{
  /** @brief Static class : No constructor authorized */
  WriterXML();
  
  /** @brief Static class : No destructor */
  ~WriterXML();
  //@}  


  /** @name Parameters */
  //@{
  static DOMWriter* m_serializer; /**< Serializer */
  static DOMDocument* m_document; /**< Racine du document */
  //@}
};

#endif /*WRITERXML_H_ */
