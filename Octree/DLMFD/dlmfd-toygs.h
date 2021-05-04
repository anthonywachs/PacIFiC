/* The ToyGS plugin */ 
# include "dlmfd-plugin.h"


#if DLM_Moving_particle
#if TRANSLATION
void compute_contact_distance (particle * p, const coord gci, double * delta) 
{
  GeomParameter gc = p->g;
  double radius = gc.radius;
  *delta = gci.y - radius;
}




double compute_fwo (const double en, const double v, const double wo) 
{

  double k = log(en)/sqrt(sq(pi) + sq(log(en)));
  double g = atan(-sqrt(1 - sq(k))/k);
  double val;
  return val = v*exp(k*g/sqrt(1- sq(k)))*sin(g)/(wo*sqrt(1 - sq(k)));
}




double derivative_fwo (const double en, const double v, const double wo) 
{
  double k = log(en)/sqrt(sq(pi) + sq(log(en)));
  double g = atan(-sqrt(1 - sq(k))/k);
  double val;
  return val = -v*exp(k*g/sqrt(1 - sq(k)))*sin(g)/(sqrt(1 - sq(k))*sq(wo)); 
}




void compute_wo (particle * p) 
{  
  /* compute an estimate of (kn,gamman) derived from (deltamax,en) */

  /* Guess value for wo */
  double wo = 100.;
  double epsilo = 1e-2;
  int maxiter = 200;
  double aa = 0.;
  double bb = 0.;
  GeomParameter gc = p->g;
  double r = gc.radius;
  double deltamax = (p->toygsp->wished_ratio)*r;
  double en = p->toygsp->en;
  double v = p->toygsp->vzero;
    
  /* Newton iteration to find the root of delta(w0) */
  if (en < 1) {
    for (int pp = 1; pp <= maxiter; pp++) {

      aa = compute_fwo (wo, v, en) - deltamax;
      bb = derivative_fwo (wo, v, en);

      if (abs(aa-deltamax)/deltamax < epsilo) {
	break;
	  }
      wo += -aa/bb;
    }
    p->toygsp->kn = (p->M)*sq(wo);
  }
  else {
    /* if en = 1 the formula simplifies as such */
    p->toygsp->kn = (p->M)*sq(p->toygsp->vzero)/(sq(deltamax));
  }  
}




void compute_Fontact (coord * Fc, particle * p, coord * gci, coord * U, 
	const double gamman) 
{
  double delta_colision = 0.;

  foreach_dimension() {
    (*Fc).x = 0.;
  }
  
  compute_contact_distance(p, *gci, &delta_colision);
    
  if (delta_colision < 0.) {
    double kn = p->toygsp->kn;
    double M = p->M;
    coord vrel = *U;
    coord Fel = {0., 0., 0.};
    coord Fdm = {0., 0., 0.};
    coord normalvec = p->toygsp->normalvector;

    foreach_dimension() {
      /* compute Hookean elastic restoring force */
      Fel.x = -kn*delta_colision*normalvec.x;
    
      /* compute viscous dynamic force */
      Fdm.x = -2.*gamman*M*fabs(vrel.x);
    
      /* Fc is the sum of these two forces */
      (*Fc).x = Fel.x + Fdm.x;
    }
  }
}




void granular_subproblem (particle * p, const int gravity_flag, 
	const double dt, const double rho_f) 
{
  // Mini Granular solver, which solves 
  //   (1 - rho_f/rho_s)MdU/dt = (1-rho_f/rho_s)Mg + sum_all_particles F_c
  //   (1 - rho_f/rho_s) Ip dw/dt = sum_all_particles (r x F_c) 
  //					- (1 - rho_f/rho_s)w x Ip w
  // where F_c is the contact force, g the gravity acceleration and 
  // M the particle's mass and Ip the inertia tensor

  double wo = 0. ;
  double gamman = 0.;
  
  double T_c = 0., dtg = 0., M = 0., kn = 0., en = 0.;
  
  int miter, gi;
  coord Uold, Xold;
  coord k1, k2, k3, k4, Xtemp, Utemp;
  double fsf;

  /* particle's structure pointers  */
  GeomParameter * gci;
  GeomParameter * gcinm1;
  coord * U;
  coord * Unm1;
 
  coord * gravity;
  coord decal = {X0, Y0, Z0};
  
  for (int k = 0; k < NPARTICLES; k++) {
    fsf = (1. - (rho_f)/(p[k].rho_s));
   
    miter = 0;
    M = p[k].M;
    kn = p[k].toygsp->kn;
    en = p[k].toygsp->en;
    
    
    /* compute wo and gamman */
    wo = sqrt(kn/M); // wo = sqrt(2kn/M) for sphere-sphere contact and 
    		     // wo = sqrt(kn/M) for sphere-wall contact 
		     // (Powder tech. Rakotonirina 2018)
    gamman = -wo*log(en)/sqrt(sq(pi) + sq(log(en)));

    /* compute the contact time Tc */
    T_c = pi/(sqrt(sq(wo) - sq(gamman)));

    /* set granular timestep */
    dtg = T_c/20.;
    miter = ceil(dt/dtg);
    dtg = dt/miter;
    
    if (k == 0) fprintf (stderr,"Tc = %g, dt = %20.18f, dtg = %20.18f, "
    	"Tc/dtg = %g, miter = %d, particle = %d\n",T_c, dt, dtg, T_c/dtg, 
	miter, k);
    
    /* get particle structure pointers */

    /* position */
    gcinm1 = &(p[k].toygsp->gnm1);
    gci = &(p[k].g);

    /* translational velocity */
    U = &(p[k].U);
    Unm1 = &(p[k].Unm1);

    gravity = &(p[k].gravity);
    
    /* fprintf (stderr,"gravity = (%f,%f,%f)\n",(*gravity).x,(*gravity).y,
    	(*gravity).z); */
    
    if (gravity_flag) {
      /* Before solving the granular problem: save previous particle's
	 position and velocity (predictor step
	 with gravity, first subproblem after N-S) */
      
      foreach_dimension() {
	(*Unm1).x = (*U).x;

	
	/* Check if the domain is periodic, if yes shift the particle's
	   position when the end of the domain is reached */
	if (Period.x) {
	  if ((*gci).center.x > (L0 + decal.x)) {
	    (*gci).center.x -= L0;
	  }
	  if ((*gci).center.x < (0. + decal.x)) {
	    (*gci).center.x += L0;
	  }
	}

	
	(*gcinm1).center.x = (*gci).center.x;
      }
    }
    
    else {
      /* corrector step without gravity: fourth subproblem after the
	 fictitious domain problem (correction step) */
      foreach_dimension() { 
	(*gci).center.x = (*gcinm1).center.x;
      }
    }
    

    /*  integrating in time with rk4 */
    for (gi = 1; gi <= miter; gi++) {
      Uold = (*U);
      Xold = (*gci).center;

      /* compute k1 */
      compute_Fontact (&k1, &p[k], &Xold, &Uold, gamman);
      foreach_dimension() {
	k1.x /= (fsf*M);
	k1.x += gravity_flag*(*gravity).x;
	Xtemp.x = Xold.x + 0.5*dtg*Uold.x;
	Utemp.x = Uold.x + 0.5*dtg*k1.x;
      }
	
      /* compute k2 */
      compute_Fontact (&k2, &p[k], &Xtemp, &Utemp, gamman);
      foreach_dimension() {
	k2.x /= (fsf*M);
	k2.x += gravity_flag*(*gravity).x;
	Xtemp.x = Xold.x + 0.5*dtg*Uold.x + 0.25*sq(dtg)*k1.x;
	Utemp.x = Uold.x + 0.5*dtg*k2.x;
      }
	
      /* compute k3 */
      compute_Fontact (&k3, &p[k], &Xtemp, &Utemp, gamman);
      foreach_dimension() {
	k3.x /= (fsf*M);
	k3.x += gravity_flag*(*gravity).x;
	Xtemp.x = Xold.x + dtg*Uold.x + 0.5*sq(dtg)*k2.x;
	Utemp.x = Uold.x + dtg*k3.x;
      }
	
      /* compute k4 */
      compute_Fontact (&k4, &p[k], &Xtemp, &Utemp, gamman);
      foreach_dimension() {
	k4.x /= (fsf*M);
	k4.x += gravity_flag*(*gravity).x;
      }

      
      foreach_dimension () {
	(*gci).center.x = Xold.x + dtg*Uold.x + sq(dtg)*(k1.x + k2.x + k3.x)/6.;
	
	/* Check if the domain is periodic, if yes shift the particle's
	   position when the end of the domain is reached */
	if (Period.x) {
	  if ((*gci).center.x > (L0 + decal.x)) {
	    (*gci).center.x -= L0;
	  }
	  if ((*gci).center.x < (0. + decal.x)) {
	    (*gci).center.x += L0;
	  }
	}



	(*U).x = Uold.x + dtg*(k1.x + 2*k2.x + 2*k3.x + k4.x)/6.;
      }
    }
    if (k == 0) {
      if (gravity_flag) {
	fprintf (stderr,"Prediction: particle-0's velocity on thread %d "
		"is (%20.18f, %20.18f, %20.18f)\n",pid(), (*U).x, (*U).y, 
		(*U).z);
	fprintf (stderr,"Prediction: particle-0's position on thread %d "
		"is (%20.18f, %20.18f, %20.18f)\n",pid(), (*gci).center.x, 
		(*gci).center.y, (*gci).center.z);
      
      }
      else {
	fprintf(stderr,"Correction: particle-0's velocity on thread %d "
		"is (%20.18f, %20.18f, %20.18f)\n",pid(), (*U).x, (*U).y, 
		(*U).z);
	fprintf(stderr,"Correction: particle-0's position on thread %d "
		"is (%20.18f, %20.18f, %20.18f)\n",pid(), (*gci).center.x, 
		(*gci).center.y, (*gci).center.z);
      }
    }
  }
}
#endif
#endif





/* Here we overload the generic events defined in the global DLMFD plugin
   dlmfd-plugin.h such that it uses my toy granular solver */

/** Overloading of the granular solver init event */
// -------------------------------------------------
event GranularSolver_init (t < -1.) {
  particle * pp = particles;
  for (int k = 0; k < NPARTICLES; k++) {
    /* Contact model parameters needed to setup the granular time-step */
    pp[k].toygsp->wished_ratio = 0.1;
    pp[k].toygsp->en = 1;
    pp[k].toygsp->vzero = 1;
    compute_wo (&pp[k]);

    /* add this term to make sure that the gravity is added only once in 
    the granular subproblem (it is already present in the granular solver) */
    foreach_dimension()
      pp[k].addforce.x = -pp[k].gravity.x;
  }
} 

/** Overloading of the granular solver predictor event */
// ------------------------------------------------------
event GranularSolver_predictor (t < -1.) {
  if ( pid() == 0 ) printf("Prediction: my rk-4 toy granular solver\n");
  particle * pp = particles;
  granular_subproblem (pp, 1, dt, rhoval);
}

/** Overloading of the granular solver velocity update event */
// ------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.) {
  if ( pid() == 0 ) printf("Correction: my rk-4 toy granular solver\n");
  particle * pp = particles;
  granular_subproblem (pp, 0, dt, rhoval);
}

