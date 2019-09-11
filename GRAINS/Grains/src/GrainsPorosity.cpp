#include "GrainsPorosity.H"
#include "Contact_BuilderFactory.hh"
#include "LinkedCell.H"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Grains_BuilderFactory.H"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriter_BuilderFactory.hh"
#include "ObstaclePeriodique.hh"
#include "Paraview_PostProcessingWriter.hh"
#include "CompParticule.hh"
#include "PointC.H"
#include "MPIWrapperGrains.hh"
#include "Composant.H"
#include "Particule.H"
#include "Forme.H"
#include "Sphere.H"
#include "Cylinder.H"
#include "Cellule.H"
#include "stdlib.h"
#include "Convex.H"
#include <cassert>
#include <string>
#include <stdlib.h>
#include <iostream>


//-----------------------------------------------------------------------------
// Constructeur par defaut
GrainsPorosity::GrainsPorosity() :
    Grains(),
  m_wrapper( NULL ),
  m_MPIstrategie( "SRLocalCommOpt" )
{
  Grains_Exec::m_isGrainsPorosity = true;
}




//-----------------------------------------------------------------------------
// Destructeur
GrainsPorosity::~GrainsPorosity()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction de la simulation.
void GrainsPorosity::Chargement( DOMElement* rootElement )
{
  if ( m_processorIsActiv )
  {
    assert(rootElement != NULL);

    DOMNode* root = ReaderXML::getNode( rootElement, "Simulation" );
    DOMNode* nPorosity = ReaderXML::getNode( root, "Porosity");
    DOMNode* nPoints = ReaderXML::getNode( nPorosity, "Points");
    DOMNode* nOut = ReaderXML::getNode( nPorosity, "Output");
    DOMNode* nStruct = ReaderXML::getNode( nPorosity, "Structure");
    m_Structure = ReaderXML::getNodeAttr_String( nStruct, "Type");
    DOMNode* nLayers = ReaderXML::getNode( nStruct, "Layers");
    DOMNode* nFenetre = ReaderXML::getNode( nStruct, "Fenetre");


    DOMNodeList* points = ReaderXML::getNodes( nFenetre );
    DOMNode* pointA = points->item( 0 );
    DOMNode* pointB = points->item( 1 );

    /* Domaine de calcul de la Porosite */
    Fenetre fenetre;
//    fenetre.ftype = FENETRE_BOX;
//    fenetre.radius = fenetre.radius_int = fenetre.hauteur = 0. ;
//    fenetre.axisdir = W ;
    fenetre.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
    fenetre.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
    fenetre.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
    fenetre.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
    fenetre.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
    fenetre.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );
    m_fenetres.insert( m_fenetres.begin(), fenetre );

    /* Points de discretisation par proc */
    m_nptsX = ReaderXML::getNodeAttr_Int( nPoints, "nx");
    m_nptsY = ReaderXML::getNodeAttr_Int( nPoints, "ny");
    m_nptsZ = ReaderXML::getNodeAttr_Int( nPoints, "nz");

    m_output_root = ReaderXML::getNodeAttr_String( nOut, "Root");
    m_output_name = ReaderXML::getNodeAttr_String( nOut, "Name");

    if ( m_Structure == "Cartesian" )
    {
      m_lx = ReaderXML::getNodeAttr_Int( nLayers, "lx");
      m_ly = ReaderXML::getNodeAttr_Int( nLayers, "ly");
      m_lz = ReaderXML::getNodeAttr_Int( nLayers, "lz");

      if ( m_rank == 0 )
      {
        cout << "\n----------------------------------------\n";
        cout << "Grains3D : module Porosity\n";
        cout << "----------------------------------------\n";
        cout << "Structure : " << m_Structure << endl;
        cout << "Number of layers : " << endl;
        cout << "   LX = " << m_lx  << endl;
        cout << "   LY = " << m_ly  << endl;
        cout << "   LZ = " << m_lz  << endl;
        cout << "Number of points on each process\n";
        cout << "nx = " << m_nptsX  << "  ny = " << m_nptsY
	     << "  nz = " << m_nptsZ << endl;
        if ( (m_ly != 1 && m_lz != 1) || (m_lx != 1 && m_lz != 1)
          || (m_lx != 1 && m_ly != 1) )
        {
          cout << "----------------------------------------\n";
          cout << "** Error : Only one direction is allowed" << endl;
          cout << "----------------------------------------\n";
          int error_code = 0;
          MPI_Abort( MPI_COMM_WORLD, error_code );
        }
      }
    }
    else if ( m_Structure == "Cylindrical" )
    {
      DOMNode* nCyl = ReaderXML::getNode( nStruct, "Cylindre");

      m_radiusIn = ReaderXML::getNodeAttr_Double( nCyl, "RadiusIn");
      m_radiusOut = ReaderXML::getNodeAttr_Double( nCyl, "RadiusOut");
      m_lr = ReaderXML::getNodeAttr_Int( nLayers, "lr");
      m_lz = ReaderXML::getNodeAttr_Int( nLayers, "lz");
      if ( m_rank == 0 )
      {
        cout << "\n----------------------------------------\n";
        cout << "Grains3D : module Porosity\n";
        cout << "----------------------------------------\n";
        cout << "Structure : " << m_Structure << endl;
        cout << "Number of layers : " << endl;
        cout << "   LR = " << m_lr  << endl;
        cout << "   LZ = " << m_lz  << endl;
        cout << "Number of points on each process\n";
        cout << "nx = " << m_nptsX  << "  ny = " << m_nptsY
	     << "  nz = " << m_nptsZ << endl;
      }
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction des forces actives
void GrainsPorosity::Forces( DOMElement* rootElement )
{
}



// ----------------------------------------------------------------------------
// Gestion de la simulation
void GrainsPorosity::Simulation( bool predict, bool isPredictorCorrector,
	bool explicit_added_mass )
{
  if ( m_processorIsActiv )
  {
    if ( m_rank == 0 )
    {
      clearResultFiles();
       fOut.open( ( m_output_root + "/" + m_output_name
	+ "_porosity.dat" ).c_str(), ios::app );
    }
    m_composants.Link( *m_sec );
    m_sec->affecteVoisinageComplet();

    // For progress display in shell purposes
    double progress = 0.0;
    int barWidth = 70;

    Scalar LX, LY, LZ;
    Scalar Xref, Yref, Zref;
    Scalar XprocMin, YprocMin, ZprocMin;
    if ( m_Structure == "Cartesian" )
    {
      // Taille des couches
      LX = fabs( m_fenetres[0].ptB[X] - m_fenetres[0].ptA[X] )/double(m_lx);
      LY = fabs( m_fenetres[0].ptB[Y] - m_fenetres[0].ptA[Y] )/double(m_ly);
      LZ = fabs( m_fenetres[0].ptB[Z] - m_fenetres[0].ptA[Z] )/double(m_lz);

      for ( int i0=0; i0<m_lx; ++i0 )
        for ( int i1=0; i1<m_ly; ++i1 )
          for ( int i2=0; i2<m_lz; ++i2 )
      {
	// Limite inferieur des couches
        Xref = min( m_fenetres[0].ptA[X], m_fenetres[0].ptB[X] ) +
	double(i0)*LX;
        Yref = min( m_fenetres[0].ptA[Y], m_fenetres[0].ptB[Y] ) +
	double(i1)*LY;
        Zref = min( m_fenetres[0].ptA[Z], m_fenetres[0].ptB[Z] ) +
	double(i2)*LZ;

        // Decomposition du domaine
        int const *coord = m_wrapper->MPI_coordonnees();
        XprocMin = Xref + coord[0] * LX/m_NX;
        YprocMin = Yref + coord[1] * LY/m_NY;
        ZprocMin = Zref + coord[2] * LZ/m_NZ;

        Scalar x, y, z;
        Scalar dx, dy, dz;
        Convex* ThePoint = new PointC() ;
        Transform transform;
        Forme* FormePoint = new Forme( ThePoint, transform );
        dx = LX / ( m_nptsX * m_NX );
        dy = LY / ( m_nptsY * m_NY );
        dz = LZ / ( m_nptsZ * m_NZ );
        Scalar Vtotal = 0.;
        list<Particule*>* particules = NULL;
        list<Particule*>* particulesVoisines = NULL;
        list<Cellule*> const* cellVoisines;
        bool isFound = false;
        Cellule* myCell = NULL;
        Point MyPoint( 0.,0.,0. );

        for ( int ii=0; ii<m_nptsX; ++ii )
        {
          if ( m_rank==0 ){
              // display progress in shell
              progress=(double(ii + m_nptsX*(i0*m_ly*m_lz+i1*m_lz+i2)))/(m_nptsX*(m_lx*m_ly*m_lz));
              std::cout << "[";
              int pos = int(double(barWidth) * progress);
              for (int i = 0; i < barWidth; ++i) {
                  if (i < pos) std::cout << "=";
                  else if (i == pos) std::cout << ">";
                  else std::cout << " ";
              }
              std::cout << "] " << int(progress * 100.0) << " %\r";
              std::cout.flush();
          }
          for ( int jj=0; jj<m_nptsY; ++jj )
            for ( int kk=0; kk<m_nptsZ; ++kk )
          {
	      isFound = false;
	      x = XprocMin + ( ii + 0.5 ) * dx;
	      y = YprocMin + ( jj + 0.5 ) * dy;
	      z = ZprocMin + ( kk + 0.5 ) * dz;
	      FormePoint->setOrigin( x, y, z );
	      MyPoint = *FormePoint->getCentre();
	      myCell = m_sec->getCellule( MyPoint );
	      particules = myCell->getParticules();
	      for ( list<Particule*>::iterator particule=particules->begin();
	        particule!=particules->end() && !isFound; particule++ )
	        if ( isIn( *particule, *FormePoint ) )
	        {
	          Vtotal += dx*dy*dz;
	          isFound = true ;
	        }
	      if ( !isFound )
	      {
	        cellVoisines = myCell->getVoisinageComplet();
	        for ( list<Cellule*>::const_iterator voisine=cellVoisines->
	          begin(); voisine!=cellVoisines->end() && !isFound;
	          voisine++ )
	        {
	          particulesVoisines = (*voisine)->getParticules();
	          for (list<Particule*>::iterator il=particulesVoisines->
	            begin(); il!=particulesVoisines->end() && !isFound; il++ )
	            if ( isIn( *il, *FormePoint ) )
	            {
		      Vtotal += dx*dy*dz;
		      isFound = true;
	            }
	        }
	      }
          }
	    }// Points
        double VolumeParticules = m_wrapper->sum_DOUBLE( Vtotal );
        if ( m_rank == 0 )
        {
          double VolumeBox = LX*LY*LZ;
	  // Calcul dans un parallelepipede
          double VolumeVideBox = VolumeBox - fabs( VolumeParticules );

	  // Si calcul dans un cylindre
	  if ( fabs( m_fenetres[0].ptB[X] - m_fenetres[0].ptA[X] ) ==
               fabs( m_fenetres[0].ptB[Y] - m_fenetres[0].ptA[Y] ) )
	  {
	    double cylRadius = 0.5 * fabs( m_fenetres[0].ptB[X]
	      - m_fenetres[0].ptA[X] );
	    double cylVolume = PI*cylRadius*cylRadius*LZ;
	    double VolumeVideCyl = cylVolume - fabs( VolumeParticules );
	    fOut << Xref + 0.5 * LX  << "\t" << Yref + 0.5 * LY
		 << "\t" << Zref + 0.5 * LZ << "\t" <<
	      ( VolumeVideBox/VolumeBox ) * 100. << "\t" <<
	      ( VolumeVideCyl/cylVolume ) * 100. << endl;
	  }
	  else
	    fOut << Xref + 0.5 * LX  << "\t" << Yref + 0.5 * LY
		 << "\t" << Zref + 0.5 * LZ << "\t" <<
	      ( VolumeVideBox/VolumeBox ) * 100. << endl;
        }
      }
      if (m_rank==0) std::cout << std::endl;

      if ( m_rank == 0 )
      {
        cout << "----------------------------------------\n";
        cout << "Output : " << m_output_root + "/" + m_output_name
            + "_porosity.dat"<< endl;
        cout << "----------------------------------------\n";
        fOut.close();
      }
    }
    else if ( m_Structure == "Cylindrical" )
    {
      // Taille des couches
      LX = fabs( m_fenetres[0].ptB[X] - m_fenetres[0].ptA[X] );
      LY = fabs( m_fenetres[0].ptB[Y] - m_fenetres[0].ptA[Y] );
      LZ = fabs( m_fenetres[0].ptB[Z] - m_fenetres[0].ptA[Z] )/
	double(m_lz);
      double LZproc = LZ/m_NZ;
      m_center = 0.5 * ( m_fenetres[0].ptB + m_fenetres[0].ptA );
      // Delta radius
      double dr = ( m_radiusOut - m_radiusIn ) / double(m_lr);

      if ( LX != LY )
      {
	if ( m_rank == 0 )
	{
	  cout << "-------------------------------------------" << "\n"
	       << "GrainsPorosity::Simulation() " << "\n"
	       << "!!! Error !!!: Check insertion windows" << "\n"
	       << "Structure is 'Cylindrical' " << "\n"
	       << "Something went wrong on X and Y directions" << "\n"
	       << "-------------------------------------------" << "\n";
	  grainsAbort();
	}
      }
      for ( int i1=0; i1<m_lz; ++i1 )
      {
      	// Limite inferieure des couches
      	Xref = min( m_fenetres[0].ptA[X], m_fenetres[0].ptB[X] );
      	Yref = min( m_fenetres[0].ptA[Y], m_fenetres[0].ptB[Y] );
      	Zref = min( m_fenetres[0].ptA[Z], m_fenetres[0].ptB[Z] ) +
      	double(i1)*LZ;
      	//double(i1)*LZ*m_NZ;
	fOut << Zref + 0.5 * LZ;
	for ( int i0=0; i0<m_lr; ++i0 )
	{
      	  /* Decomposition du domaine */
      	  int const *coord = m_wrapper->MPI_coordonnees();
      	  XprocMin = Xref + coord[0] * LX/m_NX;
      	  YprocMin = Yref + coord[1] * LY/m_NY;
      	  ZprocMin = Zref + coord[2] * LZproc;

      	  Scalar x, y, z;
      	  Scalar dx, dy, dz;
      	  Convex* ThePoint = new PointC() ;
      	  Transform transform;
      	  Forme* FormePoint = new Forme( ThePoint, transform );

	  int nbptsLayerZ = int(m_nptsZ/m_NZ);
      	  dx = LX / ( m_nptsX * m_NX );
      	  dy = LY / ( m_nptsY * m_NY );
      	  dz = LZproc / ( nbptsLayerZ );

      	  Scalar Vtotal = 0.;
      	  list<Particule*>* particules = NULL;
      	  list<Particule*>* particulesVoisines = NULL;
      	  list<Cellule*> const* cellVoisines;
      	  bool isFound = false;
      	  Cellule* myCell = NULL;
      	  Point MyPoint( 0.,0.,0. );
      	  for ( int ii=0; ii<m_nptsX; ++ii )
      	    for ( int jj=0; jj<m_nptsY; ++jj )
      	      for ( int kk=0; kk<nbptsLayerZ; ++kk )
      	  {
      	    isFound = false;
      	    x = XprocMin + ( ii + 0.5 ) * dx;
      	    y = YprocMin + ( jj + 0.5 ) * dy;
      	    z = ZprocMin + ( kk + 0.5 ) * dz;
	//    if (m_rank==3)
	//      cout << x << " \t" << y << " \t" << z << " \t" << "1\n";

      	    if ( isInRadius(x, y, z, m_radiusIn+double(i0)*dr,
      	      m_radiusIn + (double(i0)+1.)*dr) )
      	    {
	   //   cout << x << " \t" << y << " \t" << z << " \t" << "1\n";
      	      FormePoint->setOrigin( x, y, z );
      	      MyPoint = *FormePoint->getCentre();
      	      myCell = m_sec->getCellule( MyPoint );
      	      particules = myCell->getParticules();
      	      for ( list<Particule*>::iterator particule=particules->begin();
      	        particule!=particules->end() && !isFound; particule++ )
      	        if ( isIn( *particule, *FormePoint ) )
      	        {
		//  if ( m_rank == 2 )
		//    cout << x << " \t" << y << " \t" << z << " \t" << "1\n";
      	          Vtotal += dx*dy*dz;
      	          isFound = true ;
      	        }
      	      if ( !isFound )
      	      {
      	        cellVoisines = myCell->getVoisinageComplet();
      	        for ( list<Cellule*>::const_iterator voisine=cellVoisines->
      	          begin(); voisine!=cellVoisines->end() && !isFound;
      	          voisine++ )
      	        {
      	          particulesVoisines = (*voisine)->getParticules();
      	          for (list<Particule*>::iterator il=particulesVoisines->
      	            begin(); il!=particulesVoisines->end() && !isFound; il++ )
      	            if ( isIn( *il, *FormePoint ) )
      	            {
		//  if ( m_rank == 2 )
		 //   cout << x << " \t" << y << " \t" << z << " \t" << "1\n";
	              Vtotal += dx*dy*dz;
	              isFound = true;
      	            }
      	        }
      	      } // isFound
      	    }// isInRadius
      	  }// Points
      	  double VolumeParticules = m_wrapper->sum_DOUBLE( Vtotal );
      	  if ( m_rank == 0 )
      	  {
	    double rIn, rOut;
	    rIn = m_radiusIn + double(i0) * dr;
	    rOut = m_radiusIn + (double(i0)+1.) * dr;
      	    double courVolume = PI * ( rOut*rOut - rIn*rIn ) * LZ;
	    cout << "VolumeParticules = " << VolumeParticules << endl;
      	    double epsilon = 1. - fabs( VolumeParticules/courVolume );
      	    fOut << "\t" << fabs(epsilon) * 100.;
      	  }
      	} // rlayers
      	if ( m_rank == 0 ) fOut << "\n";
      }// zlayers
      if ( m_rank == 0 )
      {
        cout << "----------------------------------------\n";
        cout << "Output : " << m_output_root + "/" + m_output_name
            + "_porosity.dat"<< endl;
        cout << "----------------------------------------\n";
        fOut.close();
      }
    }// if cylindrical
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des donnees pour les simulations MPI
void GrainsPorosity::readDomainDecomposition( DOMNode* root,
  	const Scalar& lx, const Scalar& ly, const Scalar& lz )
{
//  int nprocs[3] = { 1, 1, 1 };
//  int coords[3] = { 0, 0, 0 };
//  App::setDlocale( lx, ly, lz );
//  App::setOriginelocale( nprocs, coords );
//  m_processorIsActiv = true;

  // Decomposition de domaine
  m_NZ = 1;
  DOMNode* decomp = ReaderXML::getNode( root, "DomainDecomposition" );
  m_NX = ReaderXML::getNodeAttr_Int( decomp, "NX" );
  m_NY = ReaderXML::getNodeAttr_Int( decomp, "NY" );
  if ( m_dimension == 3 ) m_NZ = ReaderXML::getNodeAttr_Int( decomp, "NZ" );
  int perx = 0, pery = 0, perz = 0;

  m_MPIstrategie = "SRLocalCommOpt";

  // Construction du wrapper MPI
  m_wrapper = new MPIWrapperGrains( m_NX, m_NY, m_NZ, perx, pery, perz );
  Grains_Exec::setComm( m_wrapper );
  Grains_Exec::m_MPI = true;
  m_processorIsActiv = m_wrapper->isActiv();
  m_rank = m_wrapper->rank_ACTIV();
  m_nprocs = m_wrapper->nombre_total_procs_ACTIV();

  // Display le wrapper & la decomposition de domaine
  if ( m_processorIsActiv ) m_wrapper->display( cout );
}




// ----------------------------------------------------------------------------
// Determine si il y a intersection entre deux formes
bool GrainsPorosity::isIn( const Forme &a, const Forme &b )
{
  return( intersect( a, b ) );
}




// ----------------------------------------------------------------------------
// Determine si il y a intersection entre une Forme et un Composant
bool GrainsPorosity::isIn( Composant* composant, const Forme &b )
{
  return( composant->intersectWithForme( b ) );
}




// ----------------------------------------------------------------------------
// Determine si un point se trouve dans une couronne r_in  < x < r_out
bool GrainsPorosity::isInRadius( const double &x, const double &y, const
  double &z, const double &inRad, const double &outRad )
{
  bool inCouronne = false;
  double distance2; // distance au carre (pythagore)
  distance2 = ( x - m_center[X] )*( x - m_center[X] ) +
    ( y - m_center[Y] )*( y - m_center[Y] );
  if ( (inRad*inRad <= distance2) && (distance2 < outRad*outRad) )
    inCouronne = true;

  return inCouronne;
}




/* Efface les fichiers resultats
--------------------------------*/
void GrainsPorosity::clearResultFiles() const
{
  if ( m_rank == 0 )
  {
    ofstream clear_file("PV_clear_exec",ios::out);
    clear_file << "#!/bin/sh" << endl;
    clear_file << "# .dat files" << endl;
    clearOneFile( clear_file, m_output_root + "/" +
	m_output_name + "_porosity.dat" );
    clear_file.close();

    system("chmod 740 PV_clear_exec");
    system("./PV_clear_exec");
    system("/bin/rm -f PV_clear_exec");
  }
}




/* Script pour effacer un fichier
---------------------------------*/
void GrainsPorosity::clearOneFile( ofstream&clear_file,
	string const& filename ) const
{
  clear_file << "if [ -f " << filename << " ]" << endl;
  clear_file << "  then" << endl;
  clear_file << "    rm -f " << filename << endl;
  clear_file << "  fi" << endl << endl;
}
