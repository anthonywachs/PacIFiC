#include "Grains_Exec.hh"
#include "InterfaceFluide3D.hh"
#include "Particule.H"
#include "MonObstacle.H"
#include "Convex.H"
#include "FormeVdW.H"
#include "Polygon.H"
#include "Sphere.H"
#include "Transform.H"
#include <fstream>
#include <string>


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
InterfaceFluide3D::InterfaceFluide3D() :
  InterfaceFluide()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// @brief Destructeur
InterfaceFluide3D::~InterfaceFluide3D()
{
}




// ============================================================================
// Version sequentielle avec lecture a partir d'un fichier de nom fixe
// Pour couplage avec GRIFF
void InterfaceFluide3D::UpdateParticulesVelocities(
	list<Particule*>& particules,
	Scalar dt,
	const bool &b_set_velocity_nm1_and_diff )
{
  ifstream velocitiesFile;
  velocitiesFile.open( "Res/particles_velocities.dat" );

  // Nouvelle vitesse des particules
  int id;
  Vecteur translation, rotation;
  list<Particule*>::iterator particule;
  for (particule=particules.begin(); particule!=particules.end();
       particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0
    	&& (*particule)->getTag() != 2 )
    {
      velocitiesFile >> id
		   >> translation[X] >> translation[Y] >> translation[Z]
		   >> rotation[X]    >> rotation[Y]    >> rotation[Z];

      (*particule)->setVitesseTranslation( translation );
      (*particule)->setVitesseRotation( rotation );

      if ( b_set_velocity_nm1_and_diff )
        (*particule)->setVelocityAndVelocityDifferencePreviousTime();
    }
  }

  velocitiesFile.close();
}




// ============================================================================
// Version sequentielle ou parallele avec lecture a partir d'un tableau
// Pour couplage avec PeliGRIFF
void InterfaceFluide3D::UpdateParticulesVelocities(
	list<Particule*>& particules,
	Scalar dt,
	const vector<vector<double> > &velocities,
	const bool &b_set_velocity_nm1_and_diff,
	const bool &b_MPI )
{
  Vecteur translation, rotation;
  list<Particule*>::iterator particule;
  unsigned id = 0;

  if ( b_MPI )
  {
    // Rem:
    // 1. le vecteur velocities contient l'ensemble des particules
    // dans le systeme et la particule numerotee ID correspond bien �
    // velocities[ID]
    // 2. en MPI, chaque processeur poss�de une partie des particules dans
    // le systeme, donc pour la mise a jour, il faut d'abord chercher si la
    // particule de numero (ID) i est situee sur ce processeur: si oui on met �
    // jour sa vitesse, si non on ne fait rien.
    list<Particule*> particules_ = particules;
    bool found = false;
    for (id=0;id<velocities.size();id++)
    {
      found = false;
      for (particule=particules_.begin();
      	particule!=particules_.end() && !found;)
      {
        if ( (*particule)->getID() == int(id) )
	{
          size_t vecSize =  (velocities[id]).size();
          if ( vecSize != 6 )
	    cout << "ERROR: the vector size of velocities in AppFluide2D_FEM "
		<< "is not 6, but " << vecSize << endl;

	  translation[X] = velocities[id][0];
          translation[Y] = velocities[id][1];
          translation[Z] = velocities[id][2];
          rotation[X] = velocities[id][3];
          rotation[Y] = velocities[id][4];
          rotation[Z] = velocities[id][5];

          (*particule)->setVitesseTranslation( translation );
          (*particule)->setVitesseRotation( rotation );

          if ( b_set_velocity_nm1_and_diff )
            (*particule)->setVelocityAndVelocityDifferencePreviousTime();

	  found = true;
	  particule = particules_.erase( particule );
	}
	else particule++;
      }
    }
  }
  else
  {
    // !!! WARNING !!!
    // Prob de numerotation a verifier car lorsque l'insertion des particules
    // s'est faite en "Aleatoire", les numeros cote Grains3D et PeliGRIFF
    // ne correspondent pas

    if ( velocities.size() != particules.size() )
      cout << "PAY ATTENTION: the particles number in AllSolidComponents and "
		<< "GRAINS are different" << endl;

    for (particule=particules.begin(); particule!=particules.end();
  	particule++, id++)
    {
      if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
      {
        size_t vecSize =  (velocities[id]).size();
        if ( vecSize != 6 )
	  cout << "ERROR: the vector size of velocities in AppFluide2D_FEM is "
		<< "not 6, but " << vecSize << endl;

        translation[X] = velocities[id][0];
        translation[Y] = velocities[id][1];
        translation[Z] = velocities[id][2];
        rotation[X] = velocities[id][3];
        rotation[Y] = velocities[id][4];
        rotation[Z] = velocities[id][5];

        (*particule)->setVitesseTranslation( translation );
        (*particule)->setVitesseRotation( rotation );

        if ( b_set_velocity_nm1_and_diff )
          (*particule)->setVelocityAndVelocityDifferencePreviousTime();
      }
    }
  }
}




// ============================================================================
// Version sequentielle avec ecriture dans un fichier
// Pour couplage avec GRIFF
void InterfaceFluide3D::WriteParticulesInFluid(
	list<Particule*> const& particules,
	const string &filename ) const
{
  ofstream particles_features( filename.c_str(), ios::out );
  particles_features.precision( 10 );
  particles_features << particules.size() << endl;

  list<Particule*>::const_iterator particule;
  int id;
  for (particule=particules.begin(), id=0; particule!=particules.end();
       particule++, id++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      // Informations : donnees de la particule
      Vecteur const* vitesseT = (*particule)->getVitesseTranslation();
      Vecteur const* vitesseR = (*particule)->getVitesseRotation();
      Point const* centre     = (*particule)->getPosition();
      Scalar        masseVol  = (*particule)->getMasseVolumique();
      Scalar        masse     = (*particule)->getMasse();
      const double* inertie   = (*particule)->getInertie();
      int ncorners = (*particule)->getForme()->getConvex()->getNbCorners();
      string particuleType = "P";
      if ( (*particule)->getNombreClonesPeriodiques() ) particuleType = "PP";

      particles_features << id <<'\t'<< ncorners <<endl;
      particles_features
	<< particuleType <<'\t'
	<< (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
	<< (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
	<< masseVol    <<'\t'<< masse       <<'\t'
	<< inertie[0]  <<'\t'<< inertie[1]  <<'\t'<< inertie[2]  <<'\t'
	<< inertie[3]  <<'\t'<< inertie[4]  <<'\t'<< inertie[5]  <<'\t'
	<< (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<'\t'
	<< endl;
      (*particule)->writePositionInFluid( particles_features );
    }
    else
    {
      particles_features << id << '\t' << "1" << endl;
      particles_features << "P "
			 << "0. 0. 0. "
			 << "0. 0. 0. "
			 << "1e8 1."
			 << "1  1  1  1  1  1" << endl;
      particles_features << "0  1"     << endl;
      particles_features << "0. 0. 0." << endl;
    }
  }

  particles_features.close();

}




// ============================================================================
// Version sequentielle avec ecriture dans un istringstream
// Pour couplage avec PeliGRIFF
void InterfaceFluide3D::WriteParticulesInFluid(
	list<Particule*> const& particules,
	list<Obstacle*> const& obstaclesToFluid,
	istringstream &is) const
{
  ostringstream particles_features;

  particles_features << particules.size() + obstaclesToFluid.size() << endl;

  list<Particule*>::const_iterator particule,clone;
  int id=0, ncorners;
  size_t nclonesper;
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;
  const double* inertie = NULL;
  string particuleType = "P";
  Scalar masseVol, masse, rayon;

  for (particule=particules.begin(), id=0; particule!=particules.end();
       particule++, id++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      // Informations : donnees de la particule
      vitesseT      = (*particule)->getVitesseTranslation();
      vitesseR      = (*particule)->getVitesseRotation();
      centre        = (*particule)->getPosition();
      masseVol      = (*particule)->getMasseVolumique();
      masse         = (*particule)->getMasse();
      inertie       = (*particule)->getInertie();
      rayon         = (*particule)->getRayon();
      ncorners = (*particule)->getNbCorners();
      //ncorners = (*particule)->getForme()->getConvex()->getNbCorners();
      nclonesper    = (*particule)->getNombreClonesPeriodiques();
      particuleType = "P";
      if (nclonesper) particuleType = "PP";

      particles_features << id <<'\t'<< ncorners <<endl;

//       particles_features
//         << particuleType <<'\t'
//         << (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
//         << (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
//         << masseVol    <<'\t'<< masse       <<'\t'
//         << inertie[0]  <<'\t'<< inertie[1]  <<'\t'<< inertie[2]  <<'\t'
//         << inertie[3]  <<'\t'<< inertie[4]  <<'\t'<< inertie[5]  <<'\t'
//         << (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<'\t'
//         << endl;

      particles_features
	<< particuleType <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*vitesseT)[X] ) <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*vitesseT)[Y] ) <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*vitesseT)[Z] ) <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*vitesseR)[X] ) <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*vitesseR)[Y] ) <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*vitesseR)[Z] ) <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		masseVol )    <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		masse )       <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		inertie[0] )  <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		inertie[1] )  <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		inertie[2] )  <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		inertie[3] )  <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		inertie[4] )  <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		inertie[5] )  <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*centre)[X] )   <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*centre)[Y] )   <<'\t'
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
		(*centre)[Z] )   <<'\t'
	<< endl;

      if ( particuleType == "PP" )
      {
        particles_features << nclonesper << endl;
        list<Particule*> const* clones = (*particule)->getPeriodicClones();
        Vecteur const* periodic_vecteur = NULL ;
        for (clone=clones->begin();clone!=clones->end();clone++)
        {
          periodic_vecteur = (*clone)->getVecteurPeriodique();
          particles_features << (*periodic_vecteur)[X] << " "
                             << (*periodic_vecteur)[Y] << " "
                             << (*periodic_vecteur)[Z] << endl;
        }
      }
      // Attention dorenavant le rayon est deplace de Forme.cpp l.556
      // vers cette sturcture car il peut evoluer en fonction du temps
      particles_features << rayon ;// <<'\t' << endl;
      (*particule)->writePositionInFluid( particles_features );
    }
    else
    {
      particles_features << id << '\t' << "1" << endl;
      particles_features << "P "
			 << "0. 0. 0. "
			 << "0. 0. 0. "
			 << "1e8 1."
			 << "1  1  1  1  1  1" << endl;
      particles_features << "0  1"     << endl;
      particles_features << "0. 0. 0." << endl;
    }
  }

  // Obstacles a cinematique imposee a transmettre au fluide
  list<Obstacle*>::const_iterator obst;
  string obstacleType = "O";
  for (obst=obstaclesToFluid.begin();obst!=obstaclesToFluid.end();obst++,id++)
  {
    // Informations : donnees de l'obstacle
    vitesseT   = (*obst)->getVitesseTranslation();
    vitesseR   = (*obst)->getVitesseRotation();
    centre     = (*obst)->getPosition();
    rayon      = (*obst)->getRayon();
    ncorners   = (*obst)->getForme()->getConvex()->getNbCorners();

    particles_features << id << '\t' << ncorners << endl;
    particles_features
      << obstacleType <<'\t'
      << (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
      << (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
      << "1000.\t 0.\t 0.\t 0.\t 0.\t 0.\t 0.\t 0.\t"
      << (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<'\t'
      << endl;
    particles_features << rayon ;
    (*obst)->writePositionInFluid( particles_features );
  }

  is.str( particles_features.rdbuf()->str() );

}




// ============================================================================
// Version parallele avec ecriture dans un istringstream
// Le vecteur allProcParticules contient l'ensemble des particules sur tous les
// proc regroupees sur le master
// Pour couplage avec PeliGRIFF
void InterfaceFluide3D::WriteParticulesInFluid(
	vector<Particule*> const* allProcParticules,
	list<Obstacle*> const& obstaclesToFluid,
	istringstream &is ) const
{
  ostringstream particles_features;

  particles_features << allProcParticules->size() + obstaclesToFluid.size()
  	<< endl;

  vector<Particule*>::const_iterator particule;
  int id=0,ncorners;
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;
  const double* inertie = NULL;
  string particuleType = "P";
  Scalar masseVol, masse, rayon;

  for (particule=allProcParticules->begin();
    particule!=allProcParticules->end(); particule++, id++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0
    	&& (*particule)->getTag() != 2 )
    {
      // Informations : donnees de la particule
      vitesseT      = (*particule)->getVitesseTranslation();
      vitesseR      = (*particule)->getVitesseRotation();
      centre        = (*particule)->getPosition();
      masseVol      = (*particule)->getMasseVolumique();
      masse         = (*particule)->getMasse();
      inertie       = (*particule)->getInertie();
      rayon         = (*particule)->getRayon();
      ncorners = (*particule)->getForme()->getConvex()->getNbCorners();
      particuleType = "P";
      if ( (*particule)->getNombreClonesPeriodiques() ) particuleType = "PP";

      particles_features << id <<'\t'<< ncorners <<endl;
      particles_features
	<< particuleType <<'\t'
	<< (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
	<< (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
	<< masseVol    <<'\t'<< masse       <<'\t'
	<< inertie[0]  <<'\t'<< inertie[1]  <<'\t'<< inertie[2]  <<'\t'
	<< inertie[3]  <<'\t'<< inertie[4]  <<'\t'<< inertie[5]  <<'\t'
	<< (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<'\t'
	<< endl;
      // Attention dorenavant le rayon est deplace de Forme.cpp l.556
      // vers cette sturcture car il peut evoluer en fonction du temps
      particles_features << rayon ;//  <<'\t' << endl;
      (*particule)->writePositionInFluid( particles_features );
    }
    else
    {
      particles_features << id << '\t' << "1" << endl;
      particles_features << "P "
			 << "0. 0. 0. "
			 << "0. 0. 0. "
			 << "1e8 1."
			 << "1  1  1  1  1  1" << endl;
      particles_features << "0  1"     << endl;
      particles_features << "0. 0. 0." << endl;
    }
  }

  // Obstacles a cinematique imposee a transmettre au fluide
  list<Obstacle*>::const_iterator obst;
  string obstacleType = "O";
  for (obst=obstaclesToFluid.begin();obst!=obstaclesToFluid.end();obst++,id++)
  {
    // Informations : donnees de l'obstacle
    vitesseT   = (*obst)->getVitesseTranslation();
    vitesseR   = (*obst)->getVitesseRotation();
    centre     = (*obst)->getPosition();
    rayon      = (*obst)->getRayon();
    ncorners   = (*obst)->getForme()->getConvex()->getNbCorners();

    particles_features << id << '\t' << ncorners << endl;
    particles_features
      << obstacleType <<'\t'
      << (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
      << (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
      << "1000.\t 0.\t 0.\t 0.\t 0.\t 0.\t 0.\t 0.\t"
      << (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<'\t'
      << endl;
    particles_features << rayon ;
    (*obst)->writePositionInFluid( particles_features );
  }

  is.str( particles_features.rdbuf()->str() );

}




// ============================================================================
// Version sequentielle avec ecriture dans un fichier
// Pour couplage avec GRIFF
void InterfaceFluide3D::WritePVGCInFluid(
	list<Particule*> const& particules,
	const string &filename) const
{
  ofstream particles_velpos(filename.c_str(),ios::out);
  particles_velpos.precision(10);
  list<Particule*>::const_iterator particule;
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;

  for (particule=particules.begin(); particule!=particules.end();
       particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      // Informations : donnees de la particule
      vitesseT = (*particule)->getVitesseTranslation();
      vitesseR = (*particule)->getVitesseRotation();
      centre   = (*particule)->getPosition();

      particles_velpos
	<< (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
	<< (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
	<< (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<
	endl;
    }
  }

  particles_velpos.close();

}




// ============================================================================
// Version sequentielle avec ecriture dans un istringstream
// Pour couplage avec PeliGRIFF
void InterfaceFluide3D::WritePVGCInFluid(
	list<Particule*> const& particules,
	istringstream &is) const
{
  ostringstream particles_velpos;
  particles_velpos.precision(10);
  list<Particule*>::const_iterator particule;
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;

  for (particule=particules.begin(); particule!=particules.end();
       particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      // Informations : donnees de la particule
      vitesseT = (*particule)->getVitesseTranslation();
      vitesseR = (*particule)->getVitesseRotation();
      centre   = (*particule)->getPosition();

      particles_velpos
	<< (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
	<< (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
	<< (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<
	endl;
    }
  }

  is.str( particles_velpos.rdbuf()->str() );

}




// ============================================================================
// Version parallele avec ecriture dans un istringstream
// Le vecteur allProcParticules contient l'ensemble des particules sur tous les
// proc regroupees sur le master
// Pour couplage avec PeliGRIFF
void InterfaceFluide3D::WritePVGCInFluid(
	vector<Particule*> const* allProcParticules,
	istringstream &is) const
{
  ostringstream particles_velpos;
  particles_velpos.precision(10);
  vector<Particule*>::const_iterator particule;
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;

  for (particule=allProcParticules->begin();
    particule!=allProcParticules->end(); particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      // Informations : donnees de la particule
      vitesseT = (*particule)->getVitesseTranslation();
      vitesseR = (*particule)->getVitesseRotation();
      centre   = (*particule)->getPosition();

      particles_velpos
	<< (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
	<< (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<'\t'
	<< (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<
	endl;
    }
  }

  is.str( particles_velpos.rdbuf()->str() );

}

// ============================================================================
// Sequential version to read particle data in istringstream from
// Direction Splitting solver on fluid side
// Added by Aashish Goyal, Feb 2021
void InterfaceFluide3D::ReadParticulesFromDSFluid(
	list<Particule*> const& particules,
	istringstream &is ) const
{
  ostringstream particles_features;

  particles_features.precision(10);
  list<Particule*>::const_iterator particule;
  Vecteur bodyF;
  Point const* centre;

  for (particule=particules.begin();
      particule!=particules.end(); particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {

      centre   = (*particule)->getPosition();
      bodyF = Vecteur(100000.,100000.,0.);

//      (*particule)->addBodyForce(bodyF);
//      (*particule)->addForceMoment(1000,1000,0,0,0,0);
      (*particule)->addForce(*centre,bodyF);
/*
      // Informations : donnees de la particule
      vitesseT = (*particule)->getVitesseTranslation();
      vitesseR = (*particule)->getVitesseRotation();
      masseVol = (*particule)->getMasseVolumique();
      masse    = (*particule)->getMasse();
      centre   = (*particule)->getPosition();
      orient   = (*particule)->getForme()->getOrientation();

      particles_features << "P" << '\t'
	<< (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'
	<< masseVol       <<'\t'<< masse          <<'\t'
	<< (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'
	<< (*vitesseR)[Z] <<
	endl;

      particles_features << "O" << '\t'
        << (orient)[0][0]   <<'\t'<< (orient)[0][1]   <<'\t'<< (orient)[0][2]   <<'\t'
        << (orient)[1][0]   <<'\t'<< (orient)[1][1]   <<'\t'<< (orient)[1][2]   <<'\t'
        << (orient)[2][0]   <<'\t'<< (orient)[2][1]   <<'\t'<< (orient)[2][2]   <<'\t'
        << endl;
*/
    }
  }
  is.str( particles_features.rdbuf()->str() );
}

// ============================================================================
// Sequential version to write particle data in istringstream for the
// use in Direction Splitting solver
// Added by Aashish Goyal, Feb 2021
void InterfaceFluide3D::WriteParticulesInDSFluid(
	list<Particule*> const& particules,
	istringstream &is ) const
{
  ostringstream particles_features;

  particles_features.precision(10);
  list<Particule*>::const_iterator particule;
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;
//  Quaternion const* qrot;
  Matrix orient;
  Scalar masseVol, masse;

  for (particule=particules.begin();
      particule!=particules.end(); particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      // Informations : donnees de la particule
      vitesseT = (*particule)->getVitesseTranslation();
      vitesseR = (*particule)->getVitesseRotation();
      masseVol = (*particule)->getMasseVolumique();
      masse    = (*particule)->getMasse();
      centre   = (*particule)->getPosition();
//      qrot     = (*particule)->getRotation();
      orient   = (*particule)->getForme()->getOrientation();

      particles_features << "P" << '\t'
	<< (*centre)[X]   <<'\t'<< (*centre)[Y]   <<'\t'<< (*centre)[Z]   <<'\t'
	<< masseVol       <<'\t'<< masse          <<'\t'
//	<< (*qrot)[0]     <<'\t'<< (*qrot)[1]     <<'\t'<< (*qrot)[2]     <<'\t'<< (*qrot)[3]    <<'\t'
	<< (*vitesseT)[X] <<'\t'<< (*vitesseT)[Y] <<'\t'<< (*vitesseT)[Z] <<'\t'
	<< (*vitesseR)[X] <<'\t'<< (*vitesseR)[Y] <<'\t'<< (*vitesseR)[Z] <<
	endl;

      particles_features << "O" << '\t'
	<< (orient)[0][0]   <<'\t'<< (orient)[0][1]   <<'\t'<< (orient)[0][2]   <<'\t'
	<< (orient)[1][0]   <<'\t'<< (orient)[1][1]   <<'\t'<< (orient)[1][2]   <<'\t'
	<< (orient)[2][0]   <<'\t'<< (orient)[2][1]   <<'\t'<< (orient)[2][2]   <<'\t'
	<< endl;
    }
  }
  is.str( particles_features.rdbuf()->str() );
}

// ============================================================================
// Version sequentielle avec ecriture dans une structure de donne
// Pour couplage avec Basilisk C code
void InterfaceFluide3D::WriteParticulesInFluid(
	list<Particule*> const& particules,
	list<Obstacle*> const& obstaclesToFluid,
	BasiliskDataStructure * basiliskdata) const
{



  list<Particule*>::const_iterator particule,clone;
  int id = 0, ncorners;
  // size_t nclonesper;
  Vecteur const* vitesseT;
  Vecteur const* vitesseR;
  Point const* centre;
  const double* inertie = NULL;
  // string particuleType = "P";
  Scalar masseVol, masse, rayon;
  int ii = 0;
  for (particule=particules.begin(), id=0; particule!=particules.end();
       particule++, id++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      basiliskdata[ii].particlessize = particules.size();

      // Informations : donnees de la particule
      vitesseT      = (*particule)->getVitesseTranslation();
      vitesseR      = (*particule)->getVitesseRotation();
      centre        = (*particule)->getPosition();
      masseVol      = (*particule)->getMasseVolumique();
      masse         = (*particule)->getMasse();
      inertie       = (*particule)->getInertie();
      rayon         = (*particule)->getRayon();
      ncorners = (*particule)->getNbCorners();

      // nclonesper    = (*particule)->getNombreClonesPeriodiques();
      // particuleType = "P";
      // if (nclonesper) particuleType = "PP";

      /* Geometric aspect */
      basiliskdata[ii].ncorners = (int) ncorners;
      basiliskdata[ii].centreX = (double)(*centre)[X];
      basiliskdata[ii].centreY = (double)(*centre)[Y];
      basiliskdata[ii].centreZ = (double) (*centre)[Z];
      basiliskdata[ii].rayon = (double) rayon;

      /* Corners coordinates */
      (*particule)->writePositionInFluid (&basiliskdata[ii]);

      /* Velocities */
      basiliskdata[ii].vitesseTX = (double) (*vitesseT)[X];
      basiliskdata[ii].vitesseTY = (double) (*vitesseT)[Y];
      basiliskdata[ii].vitesseTZ = (double) (*vitesseT)[Z];
      basiliskdata[ii].vitesseRX = (double) (*vitesseR)[X];
      basiliskdata[ii].vitesseRY = (double) (*vitesseR)[Y];
      basiliskdata[ii].vitesseRZ = (double) (*vitesseR)[Z];

      /* Physical properties */
      basiliskdata[ii].masseVol = (double) masseVol;
      basiliskdata[ii].masse = (double) masse;

      // Inertia tensor: Grains store them as
      // inertie[0] = Ixx;
      // inertie[1] = Ixy;
      // inertie[2] = Ixz;
      // inertie[3] = Iyy;
      // inertie[4] = Iyz;
      // inertie[5] = Izz;

	for (int ip = 0; ip < 6; ip++) {
	basiliskdata[ii].inertie[ip] = inertie[ip];
      }
      ii++;
    }
    else
    {
      // particles_features << id << '\t' << "1" << endl;
      // particles_features << "P "
      // 			 << "0. 0. 0. "
      // 			 << "0. 0. 0. "
      // 			 << "1e8 1."
      // 			 << "1  1  1  1  1  1" << endl;
      // particles_features << "0  1"     << endl;
      // particles_features << "0. 0. 0." << endl;
    }
  }
}



// ============================================================================
// Version sequentielle ou parallele avec lecture a partir d'un tableau
// Pour couplage avec Basilisk
void InterfaceFluide3D::UpdateParticulesVelocities(
	list<Particule*>& particules,
	Scalar dt,
	BasiliskDataStructure * b,
	const bool &b_set_velocity_nm1_and_diff,
	const bool &b_MPI )
{
  Vecteur translation, rotation;
  list<Particule*>::iterator particule;
  unsigned id = 0;
  //cout << "Inside UpdateParticulesVelocities for Basilisk" << endl;
  if ( b_MPI )
  {
    // Rem:
    // 1. le vecteur velocities contient l'ensemble des particules
    // dans le systeme et la particule numerotee ID correspond bien �
    // velocities[ID]
    // 2. en MPI, chaque processeur poss�de une partie des particules dans
    // le systeme, donc pour la mise a jour, il faut d'abord chercher si la
    // particule de numero (ID) i est situee sur ce processeur: si oui on met �
    // jour sa vitesse, si non on ne fait rien.

    // list<Particule*> particules_ = particules;
    // bool found = false;
    // for (id=0;id<velocities.size();id++)
    // {
    //   found = false;
    //   for (particule=particules_.begin();
    //   	particule!=particules_.end() && !found;)
    //   {
    //     if ( (*particule)->getID() == int(id) )
    // 	{
    //       // size_t vecSize =  (velocities[id]).size();
    //       // if ( vecSize != 6 )
    // 	  //   cout << "ERROR: the vector size of velocities in AppFluide2D_FEM "
    // 	  // 	<< "is not 6, but " << vecSize << endl;

    // 	  translation[X] = b[id].vitesseTX;
    // 	  translation[Y] = b[id].vitesseTY;
    // 	  translation[Z] = b[id].vitesseTZ;
    // 	  rotation[X] = b[id].vitesseRX;
    // 	  rotation[Y] = b[id].vitesseRY;
    // 	  rotation[Z] = b[id].vitesseRZ;

    //       (*particule)->setVitesseTranslation( translation );
    //       (*particule)->setVitesseRotation( rotation );

    //       if ( b_set_velocity_nm1_and_diff )
    //         (*particule)->setVelocityAndVelocityDifferencePreviousTime();

    // 	  found = true;
    // 	  particule = particules_.erase( particule );
    // 	}
    // 	else particule++;
    //   }
    // }

  }
  else
  {
    for (particule=particules.begin(); particule!=particules.end();
  	particule++, id++)
    {
      if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
      {

        translation[X] = b[id].vitesseTX;
        translation[Y] = b[id].vitesseTY;
        translation[Z] = b[id].vitesseTZ;
        rotation[X] = b[id].vitesseRX;
        rotation[Y] = b[id].vitesseRY;
        rotation[Z] = b[id].vitesseRZ;

        (*particule)->setVitesseTranslation( translation );
        (*particule)->setVitesseRotation( rotation );

        if ( b_set_velocity_nm1_and_diff )
          (*particule)->setVelocityAndVelocityDifferencePreviousTime();
	// cout << "Tx = " << translation[X] << endl;
	// cout << "Ty = " << translation[Y] << endl;
	// cout << "Tz = " << translation[Z] << endl;

	// cout << "Rx = " << rotation[X] << endl;
	// cout << "Ry = " << rotation[Y] << endl;
	// cout << "Rz = " << rotation[Z] << endl;

      }
    }
  }
}
