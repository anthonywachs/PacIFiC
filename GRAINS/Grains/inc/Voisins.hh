#ifndef _Voisins
#define _Voisins

#include <mpi.h>
#include <Cellule.H>
#include <vector>
#include <iostream>
#include <list>
using std::list;
using std::vector;
using std::ostream;
using std::endl;

struct MPICartInfos
{
  int rank;
  int *coords;
};


/** @brief Gestion des voisins (position - coordonnées) 
    
    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation */
//=============================================================================
class Voisins
{
  public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut 
  @param commgrainsMPI_3D Communicateur MPI lie a la topologie cartesienne 
  @param center_coords coordonnes du processeur dans la topologie cartesienne 
  @param nprocsdir nombre de processeurs dans chaque direction 
  @param period periodicite du pattern MPI */
  Voisins( MPI_Comm &commgrainsMPI_3D, int const* center_coords,
  	int const* nprocsdir, int const* period );

  /** @brief Destructeur */
  ~Voisins();
  //@}
  
  /**@name Access */
  //@{
  /** @brief Renvoie le rang d'un processeur voisin, la position attendue est
  relative (-1,0 ou +1) 
  @param i index dans la direction X (-1,0 ou 1) 
  @param j index dans la direction Y (-1,0 ou 1)   
  @param k index dans la direction Z (-1,0 ou 1) */
  int rank( int i, int j, int k ) const;
  
  /** @brief Renvoie les coordonnes MPI d'un processeur voisin, la position 
  attendue est relative (-1,0 ou +1) 
  @param i index dans la direction X (-1,0 ou 1) 
  @param j index dans la direction Y (-1,0 ou 1)   
  @param k index dans la direction Z (-1,0 ou 1) */
  int const* coordinates( int i, int j, int k ) const;
  
  /** @brief Renvoie le nombre de voisins (avec le processeur lui même) */
  int nbVoisins() const { return m_nneighbors; }
  
  /** @brief Renvoie le rang des voisins (avec le processeur lui même) 
  dans un tableau d'entiers */
  int* rangVoisins() const; 
  
  /** @brief Renvoie le rang des voisins (sans le processeur lui même) */
  list<int> const* rangVoisinsSeuls() const;
  
  /** @brief Renvoie la geolocalisation des voisins (sans le processeur 
  lui même) */
  list<MPIGeoLocalisation> const* geolocVoisinsSeuls() const;     
  //@} 
  
  /**@name I/O methods */
  //@{
  /** @brief Display 
  @param os flux de sortie 
  @param M objet Voisins */
  friend ostream& operator <<( ostream &os, Voisins const& M );
  //@}    


private:
  vector<struct MPICartInfos> m_data; /**< rang et coordonnes MPI */
  int m_nneighbors; /**< nombre de voisins */
  list<int> m_rankNeighborsOnly; /**< rang des voisins sans le processeur 
  	lui meme */
  list<MPIGeoLocalisation> m_geolocNeighborsOnly; /**< geolocalisation des 
  	voisins sans le processeur lui meme */	

  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  Voisins();
  //@}
};

#endif
