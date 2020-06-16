#ifndef _InterfaceFluide3D
#define _InterfaceFluide3D

#include "InterfaceFluide.hh"

#include <iostream>
using namespace std;


/** @brief Interface avec un code Navier & Stokes pour la resolution
    de l'hydrodynamique autour des particules

    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class InterfaceFluide3D : public InterfaceFluide
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur */
  InterfaceFluide3D();

  /** @brief Destructeur */
  ~InterfaceFluide3D();
  //@}


  /** @name Methods */
  //@{
  /**
    @brief Mise � jour de la vitesse des particules par le fluide.
    @param particules liste des particules
    @param dt Valeur de l'increment de temps.
    @param b_set_velocity_nm1_and_diff mise ou non a jour des vitesses et des
      differences de vitesse du pas de temps fluide precedent
  */
  void UpdateParticulesVelocities( list<Particule*>& particules,
        Scalar dt, const bool &b_set_velocity_nm1_and_diff );

  /**
    @brief Mise � jour de la vitesse des particules par le fluide.
    @param particules liste des particules
    @param dt Valeur de l'increment de temps
    @param velocities tableau contenant les vitesse de translation & rotation
    @param b_set_velocity_nm1_and_diff mise ou non a jour des vitesses et des
      differences de vitesse du pas de temps flduide precedent
    @param b_MPI si Grains tourne en MPI ou en sequentiel; en sequentiel, les
      particules sont numerotees de 0 � N-1 (N nombre de particules) et elles
      se trouvent toutes sur le meme processeur (le master), en MPI ce n'est
      plus le cas et il faut utliser l'ID (numero) des particules
  */
  void UpdateParticulesVelocities( list<Particule*>& particules,
        Scalar dt, const vector<vector<double> > &velocities,
        const bool &b_set_velocity_nm1_and_diff,
        const bool &b_MPI = false );

  /**
    @brief Ecriture des particules dans un fichier
    @param particules liste des particules
    @param filename nom du fichier
  */
  void WriteParticulesInFluid( list<Particule*> const& particules,
  	const string &filename ) const;

  /**
    @brief Ecriture des particules dans le istringstream
    @param particules liste des particules
    @param obstaclesToFluid liste des particules
    @param is flux d'entr�e
    @par Content
    For particles
      -
    For Obstacles
      - Loads kinematics and informations in stream
      - Calls Composant::writePositionInFluid
  */
  void WriteParticulesInFluid( list<Particule*> const& particules,
  	list<Obstacle*> const& obstaclesToFluid, istringstream &is ) const;

  /**
    @brief Ecriture des particules dans le istringstream
    @param allProcParticules vecteur contenant l'ensemble des particules sur tous
  	  les processeurs
    @param obstaclesToFluid liste des particules
    @param is flux d'entr�e
  */
  void WriteParticulesInFluid( vector<Particule*> const* allProcParticules,
      list<Obstacle*> const& obstaclesToFluid, istringstream &is ) const;

  /** @brief Ecriture de la vitesse et du centre de gravite des particules dans
  un fichier
  @param particules liste des particules
  @param filename nom du fichier */
  void WritePVGCInFluid( list<Particule*> const& particules,
      const string &filename ) const;

  /** @brief Ecriture de la vitesse et du centre de gravite des particules dans
  un flux
  @param particules liste des particules
  @param is flux d'entr�e */
  void WritePVGCInFluid( list<Particule*> const& particules,
      istringstream &is ) const;

  /** @brief Ecriture de la vitesse et du centre de gravite des particules dans
  un flux
  @param allProcParticules vecteur contenant l'ensemble des particules sur tous
  	les processeurs
  @param is flux d'entr�e */
  void WritePVGCInFluid( vector<Particule*> const* allProcParticules,
  	istringstream &is ) const;
  //@}

  /**
   @brief Ecriture des particules dans le BasiliskDataStructure
   @param allProcParticules vecteur contenant l'ensemble des particules sur tous
   les processeurs
   @param obstaclesToFluid liste des particules
   @param is structure Basilisk
*/
void WriteParticulesInFluid(list<Particule*> const& particules,
           list<Obstacle*> const& obstaclesToFluid, BasiliskDataStructure * b) const;



/**
  @brief Mise � jour de la vitesse des particules par Basilisk.
  @param particules liste des particules
  @param dt Valeur de l'increment de temps
  @param Structure de donne Basilisk
  @param b_set_velocity_nm1_and_diff mise ou non a jour des vitesses et des
    differences de vitesse du pas de temps flduide precedent
  @param b_MPI si Grains tourne en MPI ou en sequentiel; en sequentiel, les
    particules sont numerotees de 0 � N-1 (N nombre de particules) et elles
    se trouvent toutes sur le meme processeur (le master), en MPI ce n'est
    plus le cas et il faut utliser l'ID (numero) des particules
*/
void UpdateParticulesVelocities( list<Particule*>& particules,
      Scalar dt, BasiliskDataStructure * b,
      const bool &b_set_velocity_nm1_and_diff,
      const bool &b_MPI = false );






//@}

};

#endif
