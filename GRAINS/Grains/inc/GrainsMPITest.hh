#ifndef _GrainsMPITest
#define _GrainsMPITest

#include <mpi.h>
#include "MPIWrapperGrains.hh"
#include "GrainsMPI.H"

#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "ReaderXML.hh"

class Voisins;

/** @brief Classe de test du probleme en parallele GrainsMPITest.

    Gestion de la simulation par d�composition de domaines (MPI).

@author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
//=============================================================================
class GrainsMPITest : virtual public GrainsMPI
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  GrainsMPITest();

  /** @brief Destructeur */
  virtual ~GrainsMPITest();
  //@}


  /**@name Methods Virtual */
  //@{
  /** @brief Appel a la simulation granulaire.
  @param predict Simulation en etat de prediction
  @param isPredictorCorrector le sch�ma de couplage avec le fluide est il de
  type pr�dicteur-correcteur
  @param explicit_added_mass la masse ajout�e est elle trait�e de mani�re
  explicite */
  virtual void Simulation(bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false);
  //@}
};

#endif
  
