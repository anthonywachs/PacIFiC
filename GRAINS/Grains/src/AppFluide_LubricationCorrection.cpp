#include "Grains_Exec.hh"
#include "AppFluide_LubricationCorrection.H"
#include "Torseur.H"
#include "Forme.H"
#include "Box.H"
#include "LinkedCell.H"
#include "ContactLaw.hh"
#include <math.h>


// ----------------------------------------------------------------------------
// Constructeur
// A.ESTEGHAMATIAN - juin.2014 - Creation
AppFluide_LubricationCorrection::AppFluide_LubricationCorrection(
	const double &gridsize, bool is_DNS, double eps_cut ) :
  App()
{
  m_gridsize = gridsize;
  m_isDNS = is_DNS;
  m_eps_cut = eps_cut;
}




// ----------------------------------------------------------------------------
// Destructeur
AppFluide_LubricationCorrection::~AppFluide_LubricationCorrection()
{
}




// ----------------------------------------------------------------------------
// Affecte la viscosite du fluide
void AppFluide_LubricationCorrection::set_viscosity( double mu )
{
  m_viscosity = mu;
}

// ----------------------------------------------------------------------------
// Calcul la force de lubrication
void AppFluide_LubricationCorrection::computeforce( Composant* p0_,
		     Composant* p1_,
		     LinkedCell *LC,
		     Scalar dt )
{
  Point Pos1, Pos2;
  Scalar eps, gap, eps_init, eps_init_trans,
	r, eps_init_rot, coef, coefRot;
  Vecteur const* Vel1= NULL;
  Vecteur const* Vel2= NULL;
  Vecteur const* Om1 = NULL;
  Vecteur const* Om2 = NULL;

  bool compute = false ;
  const Convex* convexA = p0_->getForme()->getConvex();
  const Convex* convexB = p1_->getForme()->getConvex();
  if ( convexA->getConvexType() == SPHERE
  	&& convexB->getConvexType() == SPHERE )
   {
     compute = false ;
     // For explanations on how periodic tests work, see ERCOntact.cpp
     if ( !ContactLaw::is_ClonePer_ClonePer( p0_, p1_ ) )
       compute = true ;
     else if ( ContactLaw::are_ClonePerDirections_Perp( p0_, p1_ ) )
       compute = true ;

     if ( compute )
     {
       Vecteur lubriFNorm1, lubriFNorm2,lubriFTrans1, lubriFTrans2,
               lubriFPerpend2,lubriFPerpend1, lubriTNorm1, lubriTNorm2,
               lubriTPerpend1, lubriTPerpend2, lubriTTrans1, lubriTTrans2,
               lubriF2, lubriF1,lubriT1, lubriT2, TransUnitVec,
               PerpendUnitVec, NormUnitVec;
       Scalar Utrans1, Utrans2, Uperpend1, UrotNorm1, UrotNorm2,
              UrotPerpend1,UrotPerpend2,UrotTrans1,UrotTrans2,Unorm1,Unorm2;
       struct ResistanceMatrix RM_init, RM_init_trans, RM_init_rot,RM_sat, RM;

       // For the moment only equal-sized spheres are considered
       r = p0_->getRayon();
       Pos1 = *(p0_->getPosition());
       Pos2 = *(p1_->getPosition());
       Vel1 = p0_->getVitesseTranslation();
       Om1 = p0_->getVitesseRotation();
       Vel2 = p1_->getVitesseTranslation();
       Om2 = p1_->getVitesseRotation();
       eps = ( Norm( Pos2 - Pos1 )-( 2. * r ) ) / r;
       coef =  6. * acos(-1.) * m_viscosity * r;
       coefRot = 8. * acos(-1.) * m_viscosity * r * r;

       // Critical normalized gap to activate the model. Correlation based
       //  on several computations with different grid sizes(DNS). For meso
       //  simulation critical gap is half radius
       if (m_isDNS)
       {
         eps_init = m_gridsize *
                  ( 1. + ( 2. * r / m_gridsize -16. ) * 0.1/8. )/r;
         eps_init_trans = eps_init;
         eps_init_rot = 0.05;//eps_init;
       }
       else
       {
         eps_init = 0.5;
         eps_init_trans = 0.;
         eps_init_rot = 0.;
       }

       // if m_isDNS we use ResistanceMatrix structure, if !m_isDNS
       // we only need two elements of matrix (in meso-scale only correction normal
       // direction is considered) so we explicitely write it in the method
       // because it's faster (c.f. else if ( m_eps_cut > eps && eps > 0. && !m_isDNS )
       if (m_isDNS)
       {
       // Fill RM_init where eps = critical gap
         fill_ResistanceMatrix( &RM_init, eps_init, 0 );
         fill_ResistanceMatrix( &RM_init_trans, eps_init_trans, 0 );
         fill_ResistanceMatrix( &RM_init_rot, eps_init_rot, 0 );
       // Fill RM_sat(saturating) where eps = cut-off gap
         fill_ResistanceMatrix( &RM_sat, m_eps_cut, 0 );
       }

       if ( eps < eps_init && m_eps_cut < eps && m_isDNS )
       {
         //Center to center unit vector Direction Sphere2->Sphere1
	 //equivalent to e1 (please refer to Dance and Maxy paper)
	 //We define e1 (Normal) and e2 e3 (transverse) unit vectors based on
	 //center2center vector and velocity vector of the 2nd particle
	 //According to Kim and Karrila and also Jeffrey and Onishi
	 // local frame is always on Sphere 2.
         NormUnitVec = ( Pos1 - Pos2 )/  Norm( Pos1 - Pos2 ) ;

         Unorm1 =  (*Vel1) * NormUnitVec ;
         Unorm2 =  (*Vel2) * NormUnitVec ;
	 // equivalent to e2
	 TransUnitVec = (*Vel2) - Unorm2 * NormUnitVec;
         if ( Norm(TransUnitVec) > 1.e-20 )
         {
           // here TransUnitVec is not still unit vector
           // we first compute Utrans2 and then normalize
           Utrans2 = Norm(TransUnitVec);
           TransUnitVec = TransUnitVec/Norm(TransUnitVec);
         }
         // In case of particle/wall 2 being fix, Utrans2 =0
         // This case is treated differently
         // by choosing TransUnitVec(e2)randomly in the perpendicular
         // plane wrt the Normal vector
         else
         {
           TransUnitVec[0] = NormUnitVec[2];
           TransUnitVec[1] = NormUnitVec[2];
           TransUnitVec[2] = -NormUnitVec[0]-NormUnitVec[1] ;
           TransUnitVec = TransUnitVec/Norm(TransUnitVec);
           Utrans2 = 0.;
         }
	 Utrans1 = (*Vel1)*TransUnitVec;
	 // equivalent to e3
	 PerpendUnitVec = NormUnitVec ^ TransUnitVec;
	 // Note that Uperpend2 is zero since e3 is defind as e1^e2 and
	 // it is perpendicular to e1, e2 and accordingly U2
	 Uperpend1 = (*Vel1) * PerpendUnitVec ;
	 UrotNorm1 =  (*Om1) *	NormUnitVec;
	 UrotNorm2 =  (*Om2) * NormUnitVec;
	 UrotTrans1 =  (*Om1) * TransUnitVec;
	 UrotTrans2 =  (*Om2) * TransUnitVec;
	 UrotPerpend1 =  (*Om1) * PerpendUnitVec;
         UrotPerpend2 =  (*Om2) * PerpendUnitVec;

         fill_ResistanceMatrix( &RM, eps, 0 );

         // Lubrication force on particle 0 in normal direcction e1
	 lubriFNorm1 = coef *
       		( Unorm1 * ( RM.A11_1 - RM_init.A11_1 )
		+ Unorm2 * ( RM.A11_2 - RM_init.A11_2 ) )*NormUnitVec;
      	 // Lubrication force on particle 1 in normal direcction e1
	 lubriFNorm2 = coef *
       		( Unorm2 * ( RM.A11_1 - RM_init.A11_1 )
		+ Unorm1 * ( RM.A11_2 - RM_init.A11_2 ) )*NormUnitVec;
         //cout << " lubriFNorm1 "<<lubriFNorm1<<endl;
	 // Lubrication force on particle 0 in transverse directions e2 & e3
	 if ( eps < eps_init_trans )
	 {
	   lubriFTrans1 = coef *
	        ( Utrans1 * (RM.A22_1 - RM_init_trans.A22_1)
		+ Utrans2 * (RM.A22_2 - RM_init_trans.A22_2)
		+ r * UrotPerpend1 * (RM.B23_1 - RM_init_trans.B23_1)
		- r * UrotPerpend2 * (RM.B23_2 - RM_init_trans.B23_2) )*TransUnitVec;
	   lubriFTrans2 = coef *
	 	( Utrans2 * (RM.A22_1 - RM_init_trans.A22_1)
		+ Utrans1 * (RM.A22_2 - RM_init_trans.A22_2)
		- r * UrotPerpend2 * (RM.B23_1 - RM_init_trans.B23_1)
		+ r * UrotPerpend1 * (RM.B23_2 - RM_init_trans.B23_2) )*TransUnitVec;
	   //Remind that Uperpend2 = 0
	   lubriFPerpend1 = coef *
	 	( Uperpend1 * (RM.A22_1 - RM_init_trans.A22_1)
		+ r * UrotTrans1 * (RM.B32_1 - RM_init_trans.B32_1)
		- r * UrotTrans2 * (RM.B32_2 - RM_init_trans.B32_2) )*PerpendUnitVec;
	   lubriFPerpend2 = coef *
	 	( Uperpend1 * (RM.A22_2 - RM_init_trans.A22_2)
		+ r * UrotTrans1 * (RM.B32_2 - RM_init_trans.B32_2)
		- r * UrotTrans2 * (RM.B32_1 - RM_init_trans.B32_1) )*PerpendUnitVec;
	 }
	 else
	 {
	  lubriFTrans1 = 0. * TransUnitVec;
	  lubriFTrans2 = 0. * TransUnitVec;
	 }
	 if ( eps < eps_init_rot)
	 {
	   lubriTNorm1 = coefRot *
	   	( r * UrotNorm1 * (RM.D11_1 - RM_init_rot.D11_1)
		- r * UrotNorm2 * (RM.D11_2 - RM_init_rot.D11_2) )*NormUnitVec;
	   lubriTNorm2 = coefRot *
	   	( r * UrotNorm2 * (RM.D11_1 - RM_init_rot.D11_1)
		- r * UrotNorm1 * (RM.D11_2 - RM_init_rot.D11_2) )*NormUnitVec;
	   //Remind that Uperpend2 = 0
	   lubriTTrans1 = coefRot *
	   	( Uperpend1 * (RM.C23_1 - RM_init_rot.C23_1)
		+ r * UrotTrans1 * (RM.D22_1 - RM_init_rot.D22_1)
		- r * UrotTrans2 * (RM.D22_2 - RM_init_rot.D22_2) )*TransUnitVec;
	   lubriTTrans2 = coefRot *
	   	( - Uperpend1 * (RM.C23_2 - RM_init_rot.C23_2)
		+ r * UrotTrans2 * (RM.D22_1 - RM_init_rot.D22_1)
		- r * UrotTrans1 * (RM.D22_2 - RM_init_rot.D22_2) )*TransUnitVec;
	   //Remind that D33 = D22
	   lubriTPerpend1 = coefRot *
	   	( Utrans1 * (RM.C32_1 - RM_init_rot.C32_1)
		+ Utrans2 * (RM.C32_2 - RM_init_rot.C32_2)
		+ r * UrotPerpend1 * (RM.D22_1 - RM_init_rot.D22_1)
		- r * UrotPerpend2 * (RM.D22_2 - RM_init_rot.D22_2) )*PerpendUnitVec;
	   lubriTPerpend2 = coefRot *
	   	( - Utrans1 * (RM.C32_2 - RM_init_rot.C32_2)
		  - Utrans2 * (RM.C32_1 - RM_init_rot.C32_1)
		  - r * UrotPerpend1 * (RM.D22_2 - RM_init_rot.D22_2)
		  + r * UrotPerpend2 * (RM.D22_1 - RM_init_rot.D22_1) )*PerpendUnitVec;
	 }
	 else
	 {
          lubriTNorm1 = 0. * NormUnitVec;
	  lubriTNorm2 = 0. * NormUnitVec;
	  lubriTTrans1 = 0.* TransUnitVec;
	  lubriTTrans2 = 0.* TransUnitVec;
	  lubriTPerpend1 = 0.* PerpendUnitVec;
	  lubriTPerpend2 = 0.* PerpendUnitVec;
	 }
	 lubriF1 = lubriFNorm1+ lubriFTrans1 + lubriFPerpend1;
	 lubriF2 = lubriFNorm2+ lubriFTrans2 + lubriFPerpend2;
	 lubriT1 = lubriTNorm1 + lubriTTrans1 + lubriTPerpend1;
	 lubriT2 = lubriTNorm2 + lubriTTrans2 + lubriTPerpend2;

       } // end  eps < eps_init && m_eps_cut < eps
       else if ( m_eps_cut > eps && eps > 0. && m_isDNS )
       {
         NormUnitVec = ( Pos1 - Pos2 )/  Norm( Pos1 - Pos2 ) ;
         Unorm1 =  (*Vel1) * NormUnitVec  ;
         Unorm2 =  (*Vel2 ) * NormUnitVec  ;
         // equivalent to e2
	 TransUnitVec = (*Vel1) - Unorm2 * NormUnitVec;
         if ( Norm(TransUnitVec) > 1.e-20 )
         {
           // here TransUnitVec is not still unit vector
           // we first compute Utrans2 and then normalize
           Utrans2 = Norm(TransUnitVec);
           TransUnitVec = TransUnitVec/Norm(TransUnitVec);
         }
         // In case of particle/wall 2 being fix, Utrans2 =0
         // This case is treated differently
         // by choosing TransUnitVec(e2)randomly in the perpendicular
         // plane wrt the Normal vector
         else
         {
           TransUnitVec[0] = NormUnitVec[2];
           TransUnitVec[1] = NormUnitVec[2];
           TransUnitVec[2] = -NormUnitVec[0]-NormUnitVec[1] ;
           TransUnitVec = TransUnitVec/Norm(TransUnitVec);
           Utrans2 = 0.;
         }
	 Utrans1 = (*Vel1)*TransUnitVec;
	 // equivalent to e3
	 PerpendUnitVec = NormUnitVec ^ TransUnitVec;
	 // Note that Uperpend1 is zero since e3 is defind as e1^e2 and
	 // it is perpendicular to e1 e2 and accordingly U1
	 Uperpend1 = (*Vel1) * PerpendUnitVec ;
	 UrotNorm1 =  (*Om1) *	NormUnitVec;
	 UrotNorm2 =  (*Om2) * NormUnitVec;
	 UrotTrans1 =  (*Om1) * TransUnitVec;
	 UrotTrans2 =  (*Om2) * TransUnitVec;
	 UrotPerpend1 =  (*Om1) * PerpendUnitVec;
         UrotPerpend2 =  (*Om2) * PerpendUnitVec;

         // Lubrication force on particle 0 in normal direcction e1
	 lubriFNorm1 = coef *
       		( Unorm1 * ( RM_sat.A11_1 - RM_init.A11_1 )
		+ Unorm2 * ( RM_sat.A11_2 - RM_init.A11_2 ) )*NormUnitVec;
      	 // Lubrication force on particle 1 in normal direcction e1
	 lubriFNorm2 = coef *
       		( Unorm2 * ( RM_sat.A11_1 - RM_init.A11_1 )
		+ Unorm1 * ( RM_sat.A11_2 - RM_init.A11_2 ) )*NormUnitVec;
	 // Lubrication force on particle 0 in transverse directions e2 & e3
	   lubriFTrans1 = coef *
	 	( Utrans1 * (RM_sat.A22_1 - RM_init_trans.A22_1)
		+ Utrans2 * (RM_sat.A22_2 - RM_init_trans.A22_2)
		+ r * UrotPerpend1 * (RM_sat.B23_1 - RM_init_trans.B23_1)
		- r * UrotPerpend2 * (RM_sat.B23_2 - RM_init_trans.B23_2) )*TransUnitVec;
	   lubriFTrans2 = coef *
	 	( Utrans2 * (RM_sat.A22_1 - RM_init_trans.A22_1)
		+ Utrans1 * (RM_sat.A22_2 - RM_init_trans.A22_2)
		- r * UrotPerpend2 * (RM_sat.B23_1 - RM_init_trans.B23_1)
		+ r * UrotPerpend1 * (RM_sat.B23_2 - RM_init_trans.B23_2) )*TransUnitVec;
	   //Remind that Uperpend2 = 0
	   lubriFPerpend1 = coef *
	 	( Uperpend1 * (RM_sat.A22_1 - RM_init_trans.A22_1)
		+ r * UrotTrans1 * (RM_sat.B32_1 - RM_init_trans.B32_1)
		- r * UrotTrans2 * (RM_sat.B32_2 - RM_init_trans.B32_2) )*PerpendUnitVec;
	   lubriFPerpend2 = coef *
	 	( Uperpend1 * (RM_sat.A22_2 - RM_init_trans.A22_2)
		+ r * UrotTrans1 * (RM_sat.B32_2 - RM_init_trans.B32_2)
		- r * UrotTrans2 * (RM_sat.B32_1 - RM_init_trans.B32_1) )*PerpendUnitVec;
	   lubriTNorm1 = coefRot *
	   	( r * UrotNorm1 * (RM_sat.D11_1 - RM_init_rot.D11_1)
		- r * UrotNorm2 * (RM_sat.D11_2 - RM_init_rot.D11_2) )*NormUnitVec;
	   lubriTNorm2 = coefRot *
	   	( r * UrotNorm2 * (RM_sat.D11_1 - RM_init_rot.D11_1)
		- r * UrotNorm1 * (RM_sat.D11_2 - RM_init_rot.D11_2) )*NormUnitVec;
	   //Remind that Uperpend2 = 0
	   lubriTTrans1 = coefRot *
	   	( Uperpend1 * (RM_sat.C23_1 - RM_init_rot.C23_1)
		+ r * UrotTrans1 * (RM_sat.D22_1 - RM_init_rot.D22_1)
		- r * UrotTrans2 * (RM_sat.D22_2 - RM_init_rot.D22_2) )*TransUnitVec;
	   lubriTTrans2 = coefRot *
	   	( - Uperpend1 * (RM_sat.C23_2 - RM_init_rot.C23_2)
		+ r * UrotTrans2 * (RM_sat.D22_1 - RM_init_rot.D22_1)
		- r * UrotTrans1 * (RM_sat.D22_2 - RM_init_rot.D22_2) )*TransUnitVec;
	   //Remind that D33 = D22
	   lubriTPerpend1 = coefRot *
	   	( Utrans1 * (RM_sat.C32_1 - RM_init_rot.C32_1)
		+ Utrans2 * (RM_sat.C32_2 - RM_init_rot.C32_2)
		+ r * UrotPerpend1 * (RM_sat.D22_1 - RM_init_rot.D22_1)
		- r * UrotPerpend2 * (RM_sat.D22_2 - RM_init_rot.D22_2) )*PerpendUnitVec;
	   lubriTPerpend2 = coefRot *
	   	( - Utrans1 * (RM_sat.C32_2 - RM_init_rot.C32_2)
		  - Utrans2 * (RM_sat.C32_1 - RM_init_rot.C32_1)
		  - r * UrotPerpend1 * (RM_sat.D22_2 - RM_init_rot.D22_2)
		  + r * UrotPerpend2 * (RM_sat.D22_1 - RM_init_rot.D22_1) )*PerpendUnitVec;

       	   lubriF1 = lubriFNorm1+ lubriFTrans1 + lubriFPerpend1;
  	   lubriF2 = lubriFNorm2+ lubriFTrans2 + lubriFPerpend2;
	   lubriT1 = lubriTNorm1 + lubriTTrans1 + lubriTPerpend1;
	   lubriT2 = lubriTNorm2 + lubriTTrans2 + lubriTPerpend2;
       }
       else if ( eps < eps_init && m_eps_cut < eps && !m_isDNS )
       {
          double A11_1 =  -1./(4.*eps) + (9./40.)*log(eps) + (3./112.) * eps * log(eps)
             -  0.995;
          double A11_2 =  1./(4.*eps) - (9./40.)*log(eps) - (3./112.) * eps * log(eps)
             + 0.350;
          NormUnitVec = ( Pos1 - Pos2 )/  Norm( Pos1 - Pos2 ) ;
          Unorm1 =  (*Vel1) * NormUnitVec ;
          Unorm2 =  (*Vel2) * NormUnitVec ;
         // Lubrication force on particle 0 in normal direcction e1
	 lubriFNorm1 = coef *
       		( Unorm1 *  A11_1
		+ Unorm2 *  A11_2  )*NormUnitVec;
      	 // Lubrication force on particle 1 in normal direcction e1
	 lubriFNorm2 = coef *
       		( Unorm2 *  A11_1
		+ Unorm1 *  A11_2 )*NormUnitVec;

          lubriF1 = lubriFNorm1;
          lubriF2 = lubriFNorm2;
	  lubriT1 = 0. * (*Om1);
	  lubriT2 = 0. * (*Om2);
       }
       else if ( m_eps_cut > eps && eps > 0. && !m_isDNS )
       {
          double A11_1 =  -1./(4.*m_eps_cut) + (9./40.)*log(m_eps_cut) + (3./112.) *
		 m_eps_cut * log(m_eps_cut) - 0.995;
          double A11_2 =  1./(4.*m_eps_cut) - (9./40.)*log(m_eps_cut) - (3./112.) *
                  m_eps_cut * log(m_eps_cut) + 0.350;
          NormUnitVec = ( Pos1 - Pos2 )/  Norm( Pos1 - Pos2 ) ;
          Unorm1 =  (*Vel1) * NormUnitVec ;
          Unorm2 =  (*Vel2) * NormUnitVec ;
         // Lubrication force on particle 0 in normal direcction e1
	 lubriFNorm1 = coef *
       		( Unorm1 *  A11_1
		+ Unorm2 *  A11_2  )*NormUnitVec;
      	 // Lubrication force on particle 1 in normal direcction e1
	 lubriFNorm2 = coef *
       		( Unorm2 *  A11_1
		+ Unorm1 *  A11_2 )*NormUnitVec;

          lubriF1 = lubriFNorm1;
          lubriF2 = lubriFNorm2;
	  lubriT1 = 0. * (*Om1);
	  lubriT2 = 0. * (*Om2);
       }
       else
       {
          lubriF1 = 0. * (*Om1);
          lubriF2 = 0. * (*Om2);
	  lubriT1 = 0. * (*Om1);
	  lubriT2 = 0. * (*Om2);
       }

       // Composant 0
       if ( p0_->getID() != -2 )
       {
	  p0_->addBodyForce( lubriF1 ) ;
	  p0_->addLubriForcePP( lubriF1) ;
	  p0_->addMoment( lubriT1 ) ;
       }
       else
         if ( !ContactLaw::is_ClonePer_ParticuleWithClonePerSameDirection(
	 	p0_, p1_, LC ) )
	 {
	   p0_->addBodyForce( lubriF1 ) ;
	   p0_->addLubriForcePP( lubriF1 ) ;
	   p0_->addMoment( lubriT1 ) ;
	 }

       // Composant 1
       if ( p1_->getID() != -2 )
       {
          p1_->addBodyForce( lubriF2 ) ;
	  p1_->addLubriForcePP( lubriF2 ) ;
	  p1_->addMoment( lubriT2 ) ;
       }
       else
         if ( !ContactLaw::is_ClonePer_ParticuleWithClonePerSameDirection(
	  	p1_, p0_, LC ) )
	 {
            p1_->addBodyForce( lubriF2 ) ;
	    p1_->addLubriForcePP( lubriF2 ) ;
	    p1_->addMoment( lubriT2 ) ;
	 }
     } // end compute

   } // end Sphere-Sphere
   else if ( convexA->getConvexType() == SPHERE
  	&& convexB->getConvexType() == BOX )
   {
      // No need to test periodic particles are interactions between periodic
      // clones and walls are already discarded in
      // LinkedCell::CalculerForcesClonesPeriodiques


      Vecteur lubriFNorm1,lubriFTrans1,lubriFPerpend1,lubriTNorm1,
              lubriTPerpend1,lubriTTrans1,lubriF1,lubriT1,TransUnitVec,
              PerpendUnitVec, NormUnitVec ;
      Scalar  Utrans1,Uperpend1,Unorm1,UrotTrans1,UrotPerpend1,UrotNorm1,
              Utrans2,Unorm2;
      struct ResistanceMatrix RM_init, RM_init_rot, RM_init_trans, RM, RM_sat;
      Point ProjectedPoint( 0., 0., 0. ), Projected;
      Pos1 = *(p0_->getPosition());
      r = p0_->getRayon();
      Vel1 = p0_->getVitesseTranslation();
      Om1 = p0_->getVitesseRotation();
      Vel2 = p1_->getVitesseTranslation();

      if (m_isDNS)
      {
        eps_init = m_gridsize *
                 ( 1. + ( 2. * r / m_gridsize -16. ) * 0.1/8. )/r;
        eps_init_trans = eps_init;
        eps_init_rot = 0.05;
      }
      else
      {
        eps_init = 0.5;
        eps_init_trans = 0.5;
        eps_init_rot = 0.5;
      }

      fill_ResistanceMatrix(&RM_init, eps_init, 1);
      // Fill RM_init where eps = critical gap
      fill_ResistanceMatrix( &RM_init, eps_init, 1 );
      fill_ResistanceMatrix( &RM_init_trans, eps_init_trans, 1 );
      fill_ResistanceMatrix( &RM_init_rot, eps_init_rot, 1 );

      // nullify RM_init if we are in meso model.
      if (!m_isDNS)
      {
        nullify_ResistanceMatrix( &RM_init, eps_init, 1 );
        nullify_ResistanceMatrix( &RM_init_trans, eps_init_trans, 1 );
        nullify_ResistanceMatrix( &RM_init_rot, eps_init_rot, 1 );
      }

      Box const* convexBoxB = (Box const*)(convexB);
      const Transform* transfB = (p1_->getForme())->getTransform();
      Transform w2b;
      w2b.invert( *transfB );
      ProjectedPoint = convexBoxB->ProjectedPointSPHERE( w2b(Pos1),
      		r, gap );
      eps = gap / r;
      coef =  6. * acos(-1.) * m_viscosity * r;
      coefRot = 8. * acos(-1.) * m_viscosity * r * r;

      if ( eps < eps_init && m_eps_cut < eps )
      {
        Projected = (*transfB)( ProjectedPoint );
        // Normal vector direction sphere->box e1 direction
        NormUnitVec = ( Pos1 - Projected ) / Norm( Pos1 - Projected ) ;
        Unorm1 =  (*Vel1) * NormUnitVec ;
        Unorm2 =  (*Vel2) * NormUnitVec ;

        // e2 direction
	TransUnitVec = (*Vel2) - Unorm2 * NormUnitVec;
        if ( Norm(TransUnitVec) > 1.e-20 )
        {
          // here TransUnitVec is not still unit vector
          // we first compute Utrans2 and then normalize
          Utrans2 = Norm(TransUnitVec);
    	  TransUnitVec = TransUnitVec/Norm(TransUnitVec);
        }
        // In case of particle/wall 2 being fix, Utrans2 =0
        // This case is treated differently
        // by choosing TransUnitVec(e2)randomly in the perpendicular
        // plane wrt the Normal vector
        else
        {
          TransUnitVec[0] = NormUnitVec[2];
          TransUnitVec[1] = NormUnitVec[2];
          TransUnitVec[2] = -NormUnitVec[0]-NormUnitVec[1] ;
          TransUnitVec = TransUnitVec/Norm(TransUnitVec);
          Utrans2 = 0.;
        }
	Utrans1 = (*Vel1)*TransUnitVec;
        // equivalent to e3
	PerpendUnitVec = NormUnitVec ^ TransUnitVec;
	// Note that Uperpend2 is zero since e3 is defind as e1^e2 and
	// it is perpendicular to e1, e2 and accordingly U2
	Uperpend1 = (*Vel1) * PerpendUnitVec ;
	UrotNorm1 =  (*Om1) *	NormUnitVec;
	UrotTrans1 =  (*Om1) * TransUnitVec;
	UrotPerpend1 =  (*Om1) * PerpendUnitVec;

	fill_ResistanceMatrix(&RM, eps, 1);
        // Lubrication force on particle 0 in normal direction e1
        // Slightly different than particle-particle interaction
        // RM.A11 = RM.A11_1 = - RM.A11_2
        lubriFNorm1 = coef *
               ( Unorm1 * ( RM.A11 - RM_init.A11 )
               - Unorm2 * ( RM.A11 - RM_init.A11 ) )*NormUnitVec;

	 // Lubrication force on particle 0 in transverse directions e2 & e3
	 if ( eps < eps_init_trans && m_isDNS)
	 {
           lubriFTrans1 = coef *
                 ( Utrans1 * (RM.A22 - RM_init_trans.A22)
                 - Utrans2 * (RM.A22 - RM_init_trans.A22)
                 + r * UrotPerpend1 * (RM.B23 - RM_init_trans.B23) )*TransUnitVec;
	   lubriFPerpend1 = coef *
                 ( Uperpend1 * (RM.A22 - RM_init_trans.A22)
                 + r * UrotTrans1 * (RM.B32 - RM_init_trans.B32) )*PerpendUnitVec;
         }
	 else
	   lubriFTrans1 = 0. * TransUnitVec;
         if ( eps < eps_init_rot && m_isDNS)
         {
           lubriTNorm1 = coefRot *
               ( r * UrotNorm1 * (RM.D11 - RM_init_rot.D11) )*NormUnitVec;
           // Uperpend2 = 0
           lubriTTrans1 = coefRot *
               ( Uperpend1 * (RM.C23 - RM_init_rot.C23)
                + r * UrotTrans1 * (RM.D22 - RM_init_rot.D22) )*TransUnitVec;
           // D33 = D22
           lubriTPerpend1 = coefRot *
               ( Utrans1 * (RM.C32 - RM_init_rot.C32)
               - Utrans2 * (RM.C32 - RM_init_rot.C32)
               + r * UrotPerpend1 * (RM.D22 - RM_init_rot.D22) )*PerpendUnitVec;
         }
         else
         {
           lubriTNorm1 = 0. * NormUnitVec;
	   lubriTTrans1 = 0.* TransUnitVec;
	   lubriTPerpend1 = 0.* PerpendUnitVec;
         }
	 lubriF1 = lubriFNorm1+ lubriFTrans1 + lubriFPerpend1;
	 lubriT1 = lubriTNorm1 + lubriTTrans1 + lubriTPerpend1;
      }
      else if ( m_eps_cut > eps && eps > 0. )
      {
        Projected = (*transfB)( ProjectedPoint );

        // Normal vector direction sphere->box e1 direction
        NormUnitVec = ( Pos1 - Projected ) / Norm( Pos1 - Projected ) ;
        Unorm1 =  (*Vel1) * NormUnitVec ;
        Unorm2 =  (*Vel2) * NormUnitVec ;

        // e2 direction
	TransUnitVec = (*Vel2) - Unorm2 * NormUnitVec;
        if ( Norm(TransUnitVec) > 1.e-20 )
        {
          // here TransUnitVec is not still unit vector
          // we first compute Utrans2 and then normalize
          Utrans2 = Norm(TransUnitVec);
    	  TransUnitVec = TransUnitVec/Norm(TransUnitVec);
        }
        // In case of particle/wall 2 being fix, Utrans2 =0
        // This case is treated differently
        // by choosing TransUnitVec(e2)randomly in the perpendicular
        // plane wrt the Normal vector
        else
        {
          TransUnitVec[0] = NormUnitVec[2];
          TransUnitVec[1] = NormUnitVec[2];
          TransUnitVec[2] = -NormUnitVec[0]-NormUnitVec[1] ;
          TransUnitVec = TransUnitVec/Norm(TransUnitVec);
          Utrans2 = 0.;
        }
	Utrans1 = (*Vel1)*TransUnitVec;
        // equivalent to e3
	PerpendUnitVec = NormUnitVec ^ TransUnitVec;
	// Note that Uperpend2 is zero since e3 is defind as e1^e2 and
	// it is perpendicular to e1, e2 and accordingly U2
	Uperpend1 = (*Vel1) * PerpendUnitVec ;
	UrotNorm1 =  (*Om1) *	NormUnitVec;
	UrotTrans1 =  (*Om1) * TransUnitVec;
	UrotPerpend1 =  (*Om1) * PerpendUnitVec;

        // Fill RM_sat(saturating) where eps = cut-off gap
        fill_ResistanceMatrix( &RM_sat, m_eps_cut, 1 );

        // Lubrication force on particle 0 in normal direction e1
        // Slightly different than particle-particle interaction
        // RM.A11 = RM.A11_1 = - RM.A11_2
        lubriFNorm1 = coef *
               ( Unorm1 * ( RM_sat.A11 - RM_init.A11 )
               - Unorm2 * ( RM_sat.A11 - RM_init.A11 ) )*NormUnitVec;

	 // Lubrication force on particle 0 in transverse directions e2 & e3
         lubriFTrans1 = coef *
               ( Utrans1 * (RM_sat.A22 - RM_init_trans.A22)
               - Utrans2 * (RM_sat.A22 - RM_init_trans.A22)
               + r * UrotPerpend1 * (RM_sat.B23 - RM_init_trans.B23) )*TransUnitVec;
	 lubriFPerpend1 = coef *
                 ( Uperpend1 * (RM_sat.A22 - RM_init_trans.A22)
                 + r * UrotTrans1 * (RM_sat.B32 - RM_init_trans.B32) )*PerpendUnitVec;

         lubriTNorm1 = coefRot *
             ( r * UrotNorm1 * (RM_sat.D11 - RM_init_rot.D11) )*NormUnitVec;
         // Uperpend2 = 0
         lubriTTrans1 = coefRot *
             ( Uperpend1 * (RM_sat.C23 - RM_init_rot.C23)
              + r * UrotTrans1 * (RM_sat.D22 - RM_init_rot.D22) )*TransUnitVec;
         // D33 = D22
         lubriTPerpend1 = coefRot *
             ( Utrans1 * (RM_sat.C32 - RM_init_rot.C32)
             - Utrans2 * (RM_sat.C32 - RM_init_rot.C32)
             + r * UrotPerpend1 * (RM_sat.D22 - RM_init_rot.D22) )*PerpendUnitVec;

         if (m_isDNS)
         {
	   lubriF1 = lubriFNorm1+ lubriFTrans1 + lubriFPerpend1;
	   lubriT1 = lubriTNorm1 + lubriTTrans1 + lubriTPerpend1;
         }
         else
         {
           lubriF1 = lubriFNorm1;
	   lubriT1 = 0. * (*Vel1);
         }

      } //end els if m_eps_cut > eps && eps > 0.
      else
      {
	lubriF1 = 0. * (*Vel1);
	lubriT1 = 0. * (*Vel1);
      }

      p0_->addBodyForce( lubriF1 ) ;
      p0_->addLubriForcePP( lubriF1 ) ;
      p0_->addMoment( lubriT1 ) ;

      p1_->addBodyForce( - lubriF1 ) ;
      p1_->addLubriForcePP( - lubriF1 ) ;
   } // end Sphere-Box
   else
   {
      // According to Particule::InterAction and MonObstacle::InterAction,
      // convex A should be always Particle and Convex B can be either a
      // Particle or Box
      cout << "PAY ATTENTION : the lubrication method is only designed for "
    	"sphere/sphere or sphere/wall interactions" << endl;
      cout << "convexA = " << convexA->getConvexType() << endl;
      cout << "convexB = " << convexB->getConvexType() << endl;
   }
}




// ----------------------------------------------------------------------------
// Fill the resistance matrix for Particle-Particle or Particle-Wall
// cases (refer to Dance and Maxy paper)
void AppFluide_LubricationCorrection::fill_ResistanceMatrix
	( struct ResistanceMatrix *RM,	Scalar eps, int PPorPW )
{
  switch (PPorPW)
  {
    case 0:
     // sqeezing effect normal direction e1
     RM->A11_1 =  -1./(4.*eps) + (9./40.)*log(eps) + (3./112.) * eps * log(eps)
     	 -  0.995;
     RM->A11_2 =  1./(4.*eps) - (9./40.)*log(eps) - (3./112.) * eps * log(eps)
         + 0.350;
     // translational shearing effect transverse direction e2 & e3
     // Note that A22 = A33, so we do not define A33
     RM->A22_1 = log(eps)/6. +3.*eps*log(eps)/100. - 0.998;
     RM->A22_2 = -log(eps)/6. -3.*eps*log(eps)/100. + 0.274;
     // B23: rotational shearing effect transverse directions e2 & e3
     RM->B23_1 = -log(eps)/6. - eps * log(eps)/12. - 0.159 ;
     RM->B23_2 = log(eps)/6. + eps * log(eps)/12. +0.001;
     RM->B32_1 = - RM->B23_1;
     RM->B32_2 = - RM->B23_2;
     // Rotation correction
     RM->C23_1 = RM->B32_1;
     RM->C23_2 = RM->B32_2;
     RM->C32_1 = RM->B23_1;
     RM->C32_2 = RM->B23_2;
     RM->D11_1 = (1./8.) * eps * log(eps) - 1.052;
     RM->D11_2 = -(1./8.) * eps * log(eps) -0.150;
     // Note that D33 = D22, so we do not define D33
     RM->D22_1 = (1./5.) * log(eps) + (47./250.) * eps * log(eps) - 0.703;
     RM->D22_2 = -(1./20.) * log(eps) - (31./500.) * eps * log(eps) - 0.027;
     break;

    case 1:
     // sqeezing effect normal direction e1
     RM->A11 = -1./eps + log(eps)/5. + (1./21.) * eps * log(eps)- 0.848;
     // Previous constant : TO BE CHECKED!!!
     // RM->A11 = -1./eps + log(eps)/5. + (1./21.) * eps * log(eps) - 0.9713;
     RM->A22 = (8./15.)*log(eps) + (64./375.)*eps*log(eps) - 0.952;
     // B32 = -B23
     RM->B23 = -(2./15.)*log(eps) - (86./375.)*eps*log(eps) - 0.257;
     RM->B32 = - RM->B23;
     // Rotation correction
     RM->C23 = RM->B32;
     RM->C32 = RM->B23;
     RM->D11 = (1./2.)*eps*log(eps) - 1.202;
     RM->D22 = (2./5.)*log(eps) + (66./125.)*eps*log(eps) - 0.371;
     break;
  }

}





// ----------------------------------------------------------------------------
// Nullify the resistance matrix
void AppFluide_LubricationCorrection::nullify_ResistanceMatrix
	( struct ResistanceMatrix *RM,	Scalar eps, int PPorPW )
{
  switch (PPorPW)
  {
    case 0:
     // sqeezing effect normal direction e1
     RM->A11_1 =  0.;
     RM->A11_2 =  0.;
     // translational shearing effect transverse direction e2 & e3
     // Note that A22 = A33, so we do not define A33
     RM->A22_1 = 0.;
     RM->A22_2 = 0.;
     // B23: rotational shearing effect transverse directions e2 & e3
     RM->B23_1 = 0.;
     RM->B23_2 = 0.;
     RM->B32_1 = 0.;
     RM->B32_2 = 0.;
     // Rotation correction
     RM->C23_1 = 0.;
     RM->C23_2 = 0.;
     RM->C32_1 = 0.;
     RM->C32_2 = 0.;
     RM->D11_1 = 0.;
     RM->D11_2 = 0.;
     // Note that D33 = D22, so we do not define D33
     RM->D22_1 = 0.;
     RM->D22_2 = 0.;
     break;

    case 1:
     // sqeezing effect normal direction e1
     RM->A11 = 0.;
     // Previous constant : TO BE CHECKED!!!
     // RM->A11 = -1./eps + log(eps)/5. + (1./21.) * eps * log(eps) - 0.9713;
     RM->A22 = 0.;
     RM->B23 = 0.;
     RM->D11 = 0.;
     RM->D22 = 0.;
     break;
	}
}




// ----------------------------------------------------------------------------
// Calcul la force de lubrication
void AppFluide_LubricationCorrection::CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules )
{ cout << "CalculerForces should not be called in "<<
        "AppFluide_LubricationCorrection" << endl;

}
