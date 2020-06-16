// Gestion d'une table de correspondance des adresses.
//
// Acquisition de l'adresse d'objet a l'instant precedent
// Gestion en parallele de sa nouvelle adresse
// Un objet qui le reference peut alors connaitre sa nouvelle adresse
//
// G.FERRER - Nove.2000 - Creation
// ============================================================================
#ifndef _SaveTable
#define _SaveTable

#include <string>
#include <map>
using namespace std;



/** Gestion de la correspondance des adresses.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class Table
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation 
  @param str_ nom de l'adresse
  @param ptr_ pointeur */
  Table(const string &str_, void* ptr_) :
    str(str_), ptr(ptr_) {};
    
  /** @brief Destructeur */
  ~Table() {};
  //@}

  /**@name Methods */
  //@{
  /** @brief Adresse de l'objet a cette etape
  @return La nouvelle adresse
  @param str_ Nom de l'ancienne adresse. */
  void* getNew(const string &str_) const
    {return ((str_ == str) ? ptr : NULL);};

  /** @brief Nom_Adresse de l'objet a l'etape precedente
  @return Le nom de l'objet
  @param ptr_ Adresse de l'objet courant. */
  string getOld(const void* ptr_) const
    {return ((ptr_ == ptr) ? str : string());};
  //@}

private:
  /**@name Constructors */
  //@{
  /** @brief  Constructeur par defaut interdit */
  Table() {};
  //@}

  /**@name Parameters */
  //@{
  string str; /**< Nom de l'ancienne adresse */
  void *ptr; /**< Nouvelle adresse */
  //@}
};


struct organize
{
  bool operator()(const string &s1, const string &s2) const
  {
    return s1 < s2;
  }
};


/** @brief Gestion de la table de sauvegarde.

    Cette table permet de retrouver un objet a partir du "nom" de sa 
    precedente adresse.

    @author G.FERRER - Institut Francais du Petrole - Creation 1999 - 
            Modified 2005 (utilisation de la map) */
// ============================================================================
class SaveTable
{
public:
  /**@name Constructors */
  //@{
  /** @brief  Constructeur par defaut */
  SaveTable(){};

  /** @brief  Destructeur */
  ~SaveTable() {clear();};
  //@}

  
  /**@name Methods Static */
  //@{
  /** @brief Destruction des correspondances gerees par la table. */
  static void clear()
    {
      table.clear();
    };

  /** @brief Creation d'un item.
  @param str "Nom-adresse" de l'objet
  @param ptr Adresse de l'objet courant */
  static void create(const string &str, void* ptr)
    {
      table[str] = ptr; 
    };
    
  /** @brief Adresse de l'objet.
  @return Adresse courante (NULL si on ne le trouve pas).
  @param str "Nom-Adresse" de l'objet a l'etape precedente */
  static void* getNew(const string &str)
    {
      return table[str];
    };

  /** @brief  Nom de l'objet (correspond a l'ancienne adresse).
  @return Le nom de l'objet (vide si on ne le trouve pas)
  @param ptr_ Adresse de l'objet recherche
  static string getOld(const void* ptr_)
    {
      string str;
      list<Table*>::iterator iter;
      for (iter=tables.begin() ; iter!=tables.end(); iter++) {
	str = (*iter)->getOld(ptr_);
	if (!(str == "")) return (str);
      }
      return (str);
    } */
  //@}
  
private:
  /**@name Parameters */
  //@{
  static map<const string, void*, organize> table; /**< Table des 
  	correspondances */
  //@}
};

#endif

