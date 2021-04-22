// F.PRADEL - Janvier 2000 - Creation
// ============================================================================
#include "Torseur.H"

#include "Point.H"
#include "Vecteur.H"


//-----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur par defaut          
Torseur::Torseur()
{
}




//-----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur a partir d'un point et de deux vecteurs.
Torseur::Torseur(const Point &pt, const Vecteur &f, const Vecteur &m)
{
  ptReduction = pt;
  force       = f;
  moment      = m;
}




//-----------------------------------------------------------------------------
// Constructeur a partir d'un point et du vecteur de force
// G.FERRER - Aout.2000 - Creation
Torseur::Torseur(const Point &pt, const Vecteur &f) :
  ptReduction(pt), force(f)
{
}




//-----------------------------------------------------------------------------
// Constructeur a partir d'un point
// G.FERRER - Aout.2000 - Creation
Torseur::Torseur(const Point &pt)
{
  ptReduction = pt;
}




//-----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur a partir de pointeurs.
Torseur::Torseur(double *point, double *load, double *momentum)
{
  ptReduction = Point(point[X], point[Y], point[Z]);
  force       = Vecteur(load[X], load[Y], load[Z]);
  moment      = Vecteur(momentum[X], momentum[Y], momentum[Z]);
}




//-----------------------------------------------------------------------------
// F.PRADEL - Fevrier 2000 - Creation
// Constructeur de copie.
Torseur::Torseur(const Torseur &tau)
{
  ptReduction = tau.ptReduction;
  force       = tau.force;
  moment      = tau.moment;
}




//-----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation 
// Destructeur
Torseur::~Torseur()
{
}




// ----------------------------------------------------------------------------
// F.PRADEL - Fevrier 2000 - Creation
// Changement du point de reduction du Torseur.
void Torseur::Changer_Point(const Point &newpoint)
{
  Vecteur vecteurA = ptReduction - newpoint;
  moment += vecteurA ^ force;
  ptReduction = newpoint;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Donne le point de reduction actuel du torseur.
Point const* Torseur::getPoint() const
{
  return &ptReduction;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Donne la force du torseur.
Vecteur const* Torseur::getForce() const
{
  return &force;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Donne le moment du torseur.
Vecteur const* Torseur::getMoment() const
{
 return &moment;
}




// ----------------------------------------------------------------------------
// Affectation de la force du torseur
// G.FERRER - Janv.2002 - Creation
void Torseur::setForce(const Vecteur &force_)
{
  force = force_;
}





// ----------------------------------------------------------------------------
// Affectation du point de reduction
// G.FERRER - Aout.2000 - Creation
void Torseur::setPoint(const Point &point)
{
  ptReduction = point;
}




// ----------------------------------------------------------------------------
// Affectation de la force, du point de reduction et initialisation
// du moment a zero (utile pour initialiser le torseur avec la gravite)
void Torseur::setToBodyForce(const Point &point,const Vecteur &force_)
{
  ptReduction = point;  
  force = force_;
  moment[X]=moment[Y]=moment[Z]=0.;
}




// ----------------------------------------------------------------------------
// Ajouter une force qui agit au point de reduction du torseur
void Torseur::addForce(const Vecteur &force_)
{
  force += force_;
}




// ----------------------------------------------------------------------------
// Ajouter une force qui agit à un point different du point de 
// reduction du torseur
void Torseur::addForce(const Point &point,const Vecteur &force_)
{
  force += force_;
  Vecteur vecteurA = point - ptReduction;
  moment += (vecteurA ^ force_);
} 




// ----------------------------------------------------------------------------
// Ajouter une force qui agit au point de reduction du torseur
void Torseur::addForce(const double &fx,const double &fy,const double &fz)
{
  force[X] += fx;
  force[Y] += fy;  
  force[Z] += fz;  
} 




// ----------------------------------------------------------------------------
// Ajouter un moment
void Torseur::addMoment(const Vecteur &moment_)
{
  moment += moment_;
} 




// ----------------------------------------------------------------------------
// Ajouter un moment
void Torseur::addMoment(const double &mx,const double &my,const double &mz)
{
  moment[X] += mx;
  moment[Y] += my;  
  moment[Z] += mz;  
} 

  


// ----------------------------------------------------------------------------
// F.PRADEL - Fevrier 2000 - Creation
// Affectation de Torseur. 
Torseur& Torseur::operator = (const Torseur &rhs)
{
  if ( &rhs != this )
  {  
    ptReduction = rhs.ptReduction;
    force       = rhs.force;
    moment      = rhs.moment;
  }
  return *this;
}




// ----------------------------------------------------------------------------
// Operateur d'addition des torseurs.
// F.PRADEL - Janv.2000 - Creation
Torseur Torseur::operator + (Torseur &k2)
{
  Torseur somme;
  somme.ptReduction = ptReduction;
  somme.force = force + k2.force;
  if (ptReduction == k2.ptReduction) {
    somme.moment = moment + k2.moment;
  } else {
    Vecteur vecteurA = k2.ptReduction - ptReduction;
    somme.moment = moment + k2.moment + (vecteurA ^ k2.force);
  }   
  return somme ;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Operateur de soustraction des torseurs.
Torseur Torseur::operator - (Torseur& k2)
{
  Torseur somme;
  somme.ptReduction = ptReduction;
  somme.force= force - k2.force;
  if (ptReduction==k2.ptReduction) {
    somme.moment= moment - k2.moment;
  } else {
    Vecteur vecteurA = k2.ptReduction - ptReduction;
    somme.moment = moment - k2.moment - (vecteurA ^ k2.force);
  }   
  return somme;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Fevrier 2000 - Creation
// operateur unaire d'inversion.
Torseur Torseur::operator-()
{
  Vecteur tmpforce, tmpmom;
  tmpforce = -1.0*force;
  tmpmom   = -1.0*moment;
  
  return Torseur(ptReduction, tmpforce, tmpmom);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Operateur de soustraction des torseurs et affectation.
Torseur& Torseur::operator -= (Torseur & k2)
{
  force -= k2.force;
  if (ptReduction==k2.ptReduction) {
    moment -= k2.moment;
  } else {
    Vecteur vecteurA = k2.ptReduction - ptReduction;
    moment -= k2.moment - (vecteurA ^ k2.force);
  }   
  return *this ;
}




// ----------------------------------------------------------------------------
// Operateur d'addition des torseurs et affectation.
// F.PRADEL - Janvier 2000 - Creation
// G.FERRER - Fevr.2004 - Arrondi du vecteur (coherent MyContact::InterAction)
Torseur& Torseur::operator += (const Torseur& k2)
{
  force += k2.force;
  if (ptReduction==k2.ptReduction) {
    moment += k2.moment;
  } else {
    Vecteur vecteurA = k2.ptReduction - ptReduction;
    vecteurA.round();
    moment += k2.moment + (vecteurA ^ k2.force);
  }   
  return *this ;
}




// ----------------------------------------------------------------------------
// Incrementation du torseur courant avec un torseur dont on specifie
// le point de reduction pour le calcul de la contribution au moment provenant 
// du bras de levier ^ force appliquee. Utile pour les particules periodiques.
// A.WACHS - Sept. 2009 - Creation
void Torseur::addWithReductionPoint(const Torseur& rhs,const Point& rp_rhs)
{
  force += rhs.force;
  Vecteur vecteurA = rhs.ptReduction - rp_rhs;
  vecteurA.round();
  moment += rhs.moment + (vecteurA ^ rhs.force);    
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Operateur de comparaison.
bool Torseur::operator == (Torseur& top2)
{
  Vecteur vecteurA = top2.ptReduction - ptReduction;  
  return (moment == 
	  (top2.moment+(vecteurA^top2.force)) && force == top2.force);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Operateur de difference
bool Torseur::operator != (Torseur& top2)
{
 return !(*this==top2);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Operateur de multiplication des torseurs par un scalaire.
Torseur Torseur::operator * (double d)
{
  Torseur multi;
  multi.ptReduction=ptReduction;
  multi.force= force*d;
  multi.moment= moment*d;
  return multi ;
}




// ----------------------------------------------------------------------------
// Ecriture du torseur sous forme Force
// G.FERRER - Janv.2001 - Creation
void Torseur::write(ostream &fileOut)
{
  fileOut << "*PtContact\t" << ptReduction;
  fileOut << "*Fn+Ft\t"     << force;
  if (Norm(moment) < 1.e20)
    fileOut << "*Moment\t"    << moment;
  else
    fileOut << "*Moment   0.   0.   0.\n";
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janv.2000 - Creation
// Ecriture de l'objet
ostream &operator << (ostream &fileOut, const Torseur &objet)
{
  fileOut << "Point reduction = " << objet.ptReduction
	  << "Force = " << objet.force
	  << "Moment = " << objet.moment;
  return (fileOut);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janv.2000 - Creation
// Lecture de l'objet
istream &operator >> (istream &fileIn, Torseur &objet)
{
  fileIn >> objet.ptReduction 
	 >> objet.force 
	 >> objet.moment;
  return (fileIn);
}




// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Copie force & moment exercés sur la particule dans le vecteur fm 
// en débutant à la position i 
void Torseur::copyForceMoment(double *fm,int i) const
{
  for (int j=0 ;j<3; j++) fm[i+j] = force[j];
  for (int j=0 ;j<3; j++) fm[i+3+j] = moment[j];    
} 
