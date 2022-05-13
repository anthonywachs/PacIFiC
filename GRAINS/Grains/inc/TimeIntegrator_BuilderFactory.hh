#ifndef _TimeIntegrator_BuilderFactory
#define _TimeIntegrator_BuilderFactory

#include <string>
using std::string;
class TimeIntegrator;


/** @brief Fabrique de l'integrateur en temps

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
// ============================================================================
class TimeIntegrator_BuilderFactory
{
public:
  /** @name Methods Static */
  //@{
  /** @brief Creation de l'integrateur en temps
  @return L'integrateur en temps */
  static TimeIntegrator* create();
  //@}


private:
  /**@name Constructors */
  //@{
  /** @brief Classe Statique : pas de contructeur disponible. */
  TimeIntegrator_BuilderFactory() {};

  /** @brief Classe Statique : pas de destructeur disponible. */
  ~TimeIntegrator_BuilderFactory() {};
  //@}
};

#endif
