#ifndef _PostProcessingWriter_BuilderFactory
#define _PostProcessingWriter_BuilderFactory

#include "ReaderXML.hh"
#include <string>
using std::string;
class PostProcessingWriter;


/** @brief Fabrique d'un utilitaire de sortie des resultats

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
// ============================================================================
class PostProcessingWriter_BuilderFactory
{
public:
  /** @name Methods Static */
  //@{
  /** @brief Creation de l'integrateur en temps
  @return L'integrateur en temps 
  @param nPPW noeud XML
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  static PostProcessingWriter* create(DOMNode* nPPW,
  	int const& rank_,int const& nbranks_);
  //@}


private:
  /**@name Constructors */
  //@{
  /** @brief Classe Statique : pas de contructeur disponible. */
  PostProcessingWriter_BuilderFactory() {};

  /** @brief Classe Statique : pas de destructeur disponible. */
  ~PostProcessingWriter_BuilderFactory() {};
  //@}
};

#endif
