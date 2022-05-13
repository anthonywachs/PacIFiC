#ifndef _InterfaceFluide_BuilderFactory
#define _InterfaceFluide_BuilderFactory

class InterfaceFluide;


/** @brief Fabrique de l'interface avec l'application Fluide.

    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation */
// ============================================================================
class InterfaceFluide_BuilderFactory
{
public:
  /** @name Methods Static */
  //@{
  /** @brief Creation de l'interface de couplage avec le Fluide
  @return L'interface avec le fluide */
  static InterfaceFluide* create();
  //@}


private:
  /**@name Constructors */
  //@{
  /** @brief Classe Statique : pas de contructeur disponible. */
  InterfaceFluide_BuilderFactory() {};

  /** @brief Classe Statique : pas de destructeur disponible. */
  ~InterfaceFluide_BuilderFactory() {};
  //@}
};

#endif
