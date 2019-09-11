#include "GrainsCompFeatures.H"
#include "Contact_BuilderFactory.hh"
#include "LinkedCell.H"
#include "AppFluide_Drag.H"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Grains_BuilderFactory.H"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriter_BuilderFactory.hh"
#include "ObstaclePeriodique.hh"
#include "Paraview_PostProcessingWriter.hh"
#include "CompParticule.hh"
#include "PointC.H"
#include "Composant.H"
#include "Particule.H"
#include "Forme.H"
#include "Sphere.H"
#include "Cylinder.H"
#include "stdlib.h"
#include "Convex.H"
#include <cassert>
#include <string>
#include <iostream>
#include <time.h>


//-----------------------------------------------------------------------------
// Constructeur par defaut
GrainsCompFeatures::GrainsCompFeatures() :
  Grains()
{
  Grains_Exec::m_isGrainsCompFeatures = true;
}




//-----------------------------------------------------------------------------
// Destructeur
GrainsCompFeatures::~GrainsCompFeatures()
{

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du probleme.
void GrainsCompFeatures::Construction( DOMElement* rootElement )
{
  cout << "\n--------------------------------------" << endl;
  cout << "Grains3D : module Volume Inertie" << endl;
  cout << "--------------------------------------" << endl;
  cout << "Domain discretization" << endl;
  cout << "Evaluation of the volume" << endl;

  assert(rootElement != NULL);
  DOMNode* root = ReaderXML::getNode( rootElement, "Construction" );
  DOMNode* ppw = ReaderXML::getNode( rootElement, "PostProcessingWriters");
  DOMNode* para = ReaderXML::getNode( ppw, "CompFeatures");
  if ( !para )
  {
    cout << "!!! ERROR !!!\n"
    << "PostProcessingWriters\nshould be 'CompFeatures'\n" << endl;
    grainsAbort();
  }

  string output_root = ReaderXML::getNodeAttr_String( para, "Root");


  // Particules 
  DOMNode* particules = ReaderXML::getNode( root, "Particules" );
  if ( particules ) 
  {  
    cout << "!!! ERROR !!!\n"
	 << "GrainsCompFeatures is only for non-convex particles\n"
	 << "Check input particles\n";
    exit(10);
  }

  // Particules Composites
  DOMNode* compParticules = ReaderXML::getNode( root, "CompParticules" );
  if ( compParticules ) 
  {
    int nbPC = int(m_composants.getParticuleClassesReference()->size());
    DOMNodeList* allCompParticules = ReaderXML::getNodes( rootElement, 
    	"CompParticule" );

    // nombre de particules elemtaires 
    size_t nbreCompPart = allCompParticules->getLength();

    for (XMLSize_t i=0; i<nbreCompPart; i++) 
    {
      DOMNode* nCompParticule = allCompParticules->item( i );
      int nbre = ReaderXML::getNodeAttr_Int( nCompParticule, "Nombre" );

      // Remarque: les particules de référence ont un numéro générique -1
      // d'où "false" dans le constructeur pour auto_numbering = false
      Particule* particuleRef = new CompParticule( nCompParticule, false, 
      	nbPC+int(i) );
      m_composants.AjouterClasseParticules( particuleRef );
      pair<Particule*,int> ppp( particuleRef, nbre );
      m_newParticules.push_back( ppp );     
      
      /*------------------------------------------------*/
      /* Nombre des points de discretisation du domaine */
      /* de la particule composite                      */
      /*------------------------------------------------*/
      int nx, ny, nz;
      nx = ReaderXML::getNodeAttr_Int( nCompParticule, "nx"); 
      ny = ReaderXML::getNodeAttr_Int( nCompParticule, "ny"); 
      nz = ReaderXML::getNodeAttr_Int( nCompParticule, "nz"); 

      cout << "nx=" << nx << " |ny=" << ny << " |nz=" << nz << endl;
     
      /* recuperation du domaine pour l'evaluation du volume */
      vector<Scalar> box = particuleRef->getCompFeaturesBox();
      vector<Vecteur> InitialPos = particuleRef->getInitialRelativePositions();

      clock_t start_s=clock();
      Scalar lx, ly, lz;

      lx = box[0] - box[1];
      ly = box[2] - box[3];
      lz = box[4] - box[5];

      double dx, dy, dz, dV;
      dx = lx/(nx-1);
      dy = ly/(ny-1);
      dz = lz/(nz-1);
      dV = dx*dy*dz;
      
      Scalar volume_appr = 0.;
      double x = 0., y = 0., z = 0.;
      double xCG = 0., yCG = 0., zCG = 0.;

      Convex* ThePoint = new PointC() ;
      Transform transform;
      Forme* MyPoint = new Forme( ThePoint, transform );

      Scalar Iox, Ioy, Ioz, Ioyz, Ioxz, Ioxy;
      Iox = Ioy = Ioz = Ioyz = Ioxz = Ioxy = 0.;
      
      /* Output python visualisation */
//      ofstream f( ( output_root + "/CompParticule.dat" ).c_str(), ios::out );
      for(int ii = 0; ii < nx-1; ++ii) 
        for(int jj = 0; jj < ny-1; ++jj)
          for(int kk = 0; kk < nz-1; ++kk)
	  {
            x = box[1] + (double(ii)+0.5)*dx ;
  	    y = box[3] + (double(jj)+0.5)*dy ;
  	    z = box[5] + (double(kk)+0.5)*dz ;
  	    MyPoint->setOrigin( x, y, z );
  	    if ( isIn( particuleRef, *MyPoint ) )
  	    {
//	    f<< x << "\t" << y << "\t" << z << "\t"
//	    << isIn( particuleRef, *MyPoint ) << endl;
	      xCG += x*dV;
	      yCG += y*dV;
	      zCG += z*dV;
	      volume_appr += dV;
	      Iox  +=  ( y*y + z*z ) * dV;
	      Ioy  +=  ( x*x + z*z ) * dV;
	      Ioz  +=  ( x*x + y*y ) * dV;
	      Ioyz +=  y * z * dV;
	      Ioxz +=  x * z * dV;
	      Ioxy +=  x * y * dV;
            }
  	  
	  }
      /* Centre de gravite dans le repere local */
      xCG = (xCG / volume_appr);
      yCG = (yCG / volume_appr);
      zCG = (zCG / volume_appr);
      
//      f.close();

      /*-----------------------------------------------------------*/
      /* La matrice d'inertie est sous la forme                    */
      /*          | A -F -E |                                      */
      /* [I(S)] = |-F  B -D |                                      */
      /* 	  |-E -D  C |                                      */
      /* A = Iox, B = Ioy, C = Ioz, D = Ioyz, E = Ioxz et F = Ioxy */
      /*-----------------------------------------------------------*/

      ofstream fOut( ( output_root + "/inertia.xml" ).c_str(), ios::out );
      fOut << "        <Inertie>" << endl;
      fOut << "           <Volume VOLUME=\"" << scientific << 
      setprecision(15)<< volume_appr << "\"/>" << endl; 
      fOut << "           <GravityCenter xCG=\"" << xCG << "\"" << endl;
      fOut << "                          yCG=\"" << yCG << "\"" << endl;
      fOut << "                          zCG=\"" << zCG << "\"/>"<< endl; 
      fOut << "           <Tenseur Ixx=\"" << Iox << "\"" << endl;
      fOut << "                    Iyy=\"" << Ioy << "\"" << endl;
      fOut << "                    Izz=\"" << Ioz << "\"" << endl;
      fOut << "                    Iyz=\"" << Ioyz << "\"" << endl;
      fOut << "                    Ixz=\"" << Ioxz << "\"" << endl;
      fOut << "                    Ixy=\"" << Ioxy << "\"/>" << endl; 
      fOut << "        </Inertie>" << endl;
      fOut.close();

      clock_t stop_s=clock();
      string var = "...";
      cout << "Evaluation time: " <<
      double(stop_s-start_s)/double(CLOCKS_PER_SEC) << " s" << endl;
      cout << "--------------------------------------" << endl;
      cout << "Output file : "<< output_root << "inertia.xml" << endl;
      cout << "--------------------------------------" << endl;
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction de la simulation.
void GrainsCompFeatures::Chargement( DOMElement* rootElement )
{
  assert(rootElement != NULL);
  DOMNode* root = ReaderXML::getNode( rootElement, "Simulation" );

  // Post-processing writers
  DOMNode* nPostProcessors = ReaderXML::getNode( root, 
  	"PostProcessingWriters" );
  if ( nPostProcessors )
  {
    DOMNodeList* allPPW = ReaderXML::getNodes( nPostProcessors );
    for (XMLSize_t i=0; i<allPPW->getLength(); i++)
    {      
      DOMNode* nPPW = allPPW->item( i );
      PostProcessingWriter* ppw = PostProcessingWriter_BuilderFactory::create(
    		nPPW, m_rank, m_nprocs );
      if ( ppw ) m_composants.addPostProcessingWriter( ppw );
    }
  }	
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction des forces actives
void GrainsCompFeatures::Forces( DOMElement* rootElement )
{
}



// ----------------------------------------------------------------------------
// Gestion de la simulation
void GrainsCompFeatures::Simulation( bool predict, bool isPredictorCorrector,
	bool explicit_added_mass )
{
  // IMPORTANT: quel que soit le systeme etudie, si celui ci contient
  // N particules et M obstacles, les particules sont numerotees de 0 à N-1
  // et les obstacles de N à N+M-1
  // Dans le cas de reload avec insertion supplementaire, cela necessite de
  // renumeroter les obstacles
  
  int numPartMax = numeroMaxParticules();
  if (numPartMax) ++numPartMax;  
  list< pair<Particule*,int> >::iterator ipart;

  // Construction des nouvelles particules
  Composant::setNbComposantsCrees( numPartMax );
  for (ipart=m_newParticules.begin();ipart!=m_newParticules.end();ipart++)
  {
    int nbre = ipart->second;
    for (int ii=0; ii<nbre; ii++) 
    {
      Particule* particule = ipart->first->createCloneCopy();
      m_composants.Ajouter( particule );
    }
  }
  
//  Grains::InsertCreateNewParticules();
  m_composants.PostProcessing_start( m_temps, m_dt, m_sec, m_fenetres );
  m_composants.PostProcessing_end();
}




// ----------------------------------------------------------------------------
// Determine si il y a intersection entre deux formes
bool GrainsCompFeatures::isIn( const Forme &a, const Forme &b )
{
  return( intersect( a, b ) );
}




// ----------------------------------------------------------------------------
// Determine si il y a intersection entre une Forme et un Composant
bool GrainsCompFeatures::isIn( Composant* composant, const Forme &b )
{
  return( composant->intersectWithForme( b ) );
}
