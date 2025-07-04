/**
# Flow past a fixed prolate spheroid, rotated or not, at different Reynolds number

We solve here the Navier-Stokes equations and add the prolate spheroid
using an [embedded boundary](/src/embed.h). */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myspheroid.h"
#include "../myperfs.h"
#include "view.h"

/**
We also use the $\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) for vortex detection. */

#include "lambda2.h"

/**
## Reference solution

We first define the geometry of the spheroid, with $r$ the equatorial
radius and $\gamma$ the aspect ratio ($a=r\gamma$ is the polar
radius). */

#if GAMMA // 1, 2, 3, 4, 5, 6
#define gamma  ((double) (GAMMA))
#else // gamma = 5/2
#define gamma  (5./2.) // aspect ratio, a/r
#endif // GAMMA
#define radius (0.5)     // Equatorial radius
#define d      (spheroid_diameter_sphere ((radius), (gamma))) // Diameter of the sphere of equivalent volume

/**
We then define the incidience angle $\phi$ (rotation along the
$z$-axis). */

#if PHI // 0, 22, 45, 66, 90
#define phi   (((double) (PHI))/180.*M_PI)
#else // phi = 0
#define phi   (0.) // incidence angle
#endif // PHI

/**
We then define the Reynolds number. */

#if RE // 1, 10, 100
#define Re   ((double) (RE))
#else // Re = 0.
#define Re   (0.) // Particle Reynolds number Re = ud/nu
#endif // RE

#define uref (1.) // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u

/**
We also define the shape of the domain. Note here that *p_p* is the
position of the center of the spheroid. */

#define spheroid(x,y,z) (sq ((x)/((radius)*(gamma))) + sq ((y)/(radius)) + sq ((z)/(radius)) - 1.)
coord p_p;

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar psi[];
  foreach_vertex() {
    
    // Rotated coordinates around the z-axis at the position p
    double xrot = p.x + (x - p.x)*cos ((phi)) - (y - p.y)*sin ((phi));
    double yrot = p.y + (x - p.x)*sin ((phi)) + (y - p.y)*cos ((phi));
    double zrot = z;
    
    // Distance function
    psi[] = (spheroid ((xrot - p.x), (yrot - p.y), (zrot - p.z)));
  }
  boundary ({psi});
  fractions (psi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We also define a reference velocity field. */

scalar un[];

/**
Finally, we define the mesh adaptation parameters. */

#if DLENGTH == 32
#define lmin (6) // Min mesh refinement level (l=6 is 2pt/d for L0=32)
#elif DLENGTH == 64
#define lmin (7) // Min mesh refinement level (l=7 is 2pt/d for L0=64)
#elif DLENGTH == 128
#define lmin (8) // Min mesh refinement level (l=8 is 2pt/d for L0=128)
#elif DLENGTH == 256
#define lmin (9) // Min mesh refinement level (l=9 is 2pt/d for L0=256)
#else // L0 = 128
#define lmin (8) // Min mesh refinement level (l=8 is 2pt/d for L0=128)
#endif // DLENGTH

#if LMAX // 10, 11, 12, 13, 14
#define lmax ((int) (LMAX))
#else // 12
#define lmax (12) // Max mesh refinement level (l=12 is 32pt/d for L0=128)
#endif // LMAX

#if CMAX
#define cmax (((double) (CMAX))*1.e-3*(uref))
#else // cmax = 1.e-2
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field
#endif // CMAX

int main ()
{
  /**
  The domain is $128\times 128 \times 128$. */

#if DLENGTH // 32, 64, 128, 256
  L0 = ((double) (DLENGTH))*(d);
#else // L0 = 128
  L0 = 128.*(d);
#endif // DLENGTH
  size (L0);
  origin (-(L0)/2., -(L0)/2., -(L0)/2.);
  
  /**
  We set the maximum timestep. */

#if DTMAX
  DT = ((double) (DTMAX))*1.e-3*(tref);
#else // DT = 1.e-2  
  DT = 1.e-2*(tref);
#endif // DTMAX
 
  /**
  We set the tolerance of the Poisson solver. */
  
#if RE == 0
  stokes = true;
#endif // RE
  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  
  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We use inlet boundary conditions. */

u.n[left] = dirichlet ((uref));
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
u.r[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = (uref);
uf.n[bottom] = 0;
uf.n[top]    = 0;
uf.n[back]   = 0;
uf.n[front]  = 0;

/**
## Properties */

event properties (i++)
{
  foreach_face() {
#if RE
    muv.x[] = (uref)*(d)/(Re)*fm.x[];
#else // Re = 0
    muv.x[] = fm.x[];
#endif // RE
  }
  boundary ((scalar *) {muv});
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p, pf})
    s.third = false;
#else
  for (scalar s in {u, p, pf})
    s.third = true;
#endif // ORDER2

  /**
  We use a slope-limiter to reduce the errors made in small-cells. */

#if SLOPELIMITER
  for (scalar s in {u}) {
    s.gradient = minmod2;
  }
#endif // SLOPELIMITER
  
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the full cell center of
  cut-cells. */
#endif // TREE

  /**
  We first define the particle's position. */

  foreach_dimension()
    p_p.x = 0.;

  /**
  If the simulation is not restarted, we define the initial mesh and
  the initial velocity. */

  if (!restore (file = "restart")) {

#if TREE
    /**
    When using *TREE*, we refine the mesh around the embedded
    boundary. */
  
    astats ss;
    int ic = 0;
    do {
      ic++;
      p_shape (cs, fs, p_p);
      ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			  maxlevel = (lmax), minlevel = (1));
    } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE

    p_shape (cs, fs, p_p);

    /**
    We initialize the velocity. */

    foreach() 
      u.x[] = cs[]*(uref);
    boundary ((scalar *) {u});
  }
  else { // restart

    /**
    When restarting a simulation, only the volume fraction *cs* is
    dumped and we therefore need to re-initialize the face fraction
    *fs* after a restart. */
  
    p_shape (cs, fs, p_p);
  }

  /**
  Whether restarting or not, we define the volume fraction at the
  previous timestep *csm1=cs*. */
    
  csm1 = cs;
    
  /**
  We define the boundary conditions for the velocity. */
   
  u.n[embed] = dirichlet (0);
  u.t[embed] = dirichlet (0);
  u.r[embed] = dirichlet (0);
  p[embed]   = neumann (0);

  uf.n[embed] = dirichlet (0);
  uf.t[embed] = dirichlet (0);
  uf.r[embed] = dirichlet (0);
  pf[embed]   = neumann (0);
  
  /**
  We finally initialize, even when restarting, the reference velocity
  field. */
  
  foreach()
    un[] = u.x[];
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We also reset the embedded fractions to avoid interpolation errors
  on the geometry. This would affect in particular the computation of
  the pressure contribution to the hydrodynamic forces. */

  p_shape (cs, fs, p_p);
}
#endif // TREE

/**
## Restarts and dumps

Every 100 time steps, we dump the fluid data for restarting
purposes. We use a relatively small time step interval as the
simulations with a large value of *lmax* are relatively slow. */

#if DUMP
event dump_data (i += 100)
{
  // Dump fluid
  dump ("dump");
}
#endif // DUMP

event dump_end (t = end)
{
  // Dump fluid
  dump ("dump-final");  
}

/**
## Profiling */

#if TRACE > 1
event profiling (i += 20)
{
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

/**
## Outputs */

double CDm1 = 0., Tzm1 = 0.;

event logfile (i++; t < 200.*(tref))
{
  double du = change (u.x, un);

  coord Fp, Fmu, Tp, Tmu;
  embed_force  (p, u, mu, &Fp, &Fmu);
  embed_torque (p, u, mu, (p_p), &Tp, &Tmu);

  // Drag and lifts
  double CD  = (Fp.x + Fmu.x)/(0.5*sq ((uref))*(M_PI)*sq ((d))/4.);
  double CLy = (Fp.y + Fmu.y)/(0.5*sq ((uref))*(M_PI)*sq ((d))/4.);
  double CLz = (Fp.z + Fmu.z)/(0.5*sq ((uref))*(M_PI)*sq ((d))/4.);

  // Torques and pitching torque
  double Tx = (Tp.x + Tmu.x)/(0.5*sq ((uref))*(M_PI)*sq ((d))/8.);
  double Ty = (Tp.y + Tmu.y)/(0.5*sq ((uref))*(M_PI)*sq ((d))/8.);
  double Tz = (Tp.z + Tmu.z)/(0.5*sq ((uref))*(M_PI)*sq ((d))/8.);

  double dF = fabs (CDm1 - CD);
  double dT = fabs (Tzm1 - Tz);
  CDm1 = CD, Tzm1 = Tz;  

  fprintf (stderr, "%g %g %g %g %g %g %d %d %d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   gamma, phi, Re,
	   spheroid_eccentricity ((gamma)), spheroid_sphericity ((gamma)),
	   (L0)/(d), (lmin), (lmax),
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD, CLy, CLz, dF,
	   Tx, Ty, Tz, dT,
	   du);
  fflush (stderr);
}

/**
## Snapshots */

event snapshot (t += 5.*(tref))
{
  scalar n2u[];
  foreach() {
    if (cs[] <= 0.)
      n2u[] = nodata;
    else
      n2u[] = sqrt (sq (u.x[]) + sq (u.y[]) + sq (u.z[]));
  }
  boundary ({n2u});
  stats sn2u = statsf (n2u);
  
  char name2[80];
  
  /**
  We first plot the entire domain. */
 
  view (fov = 20, camera = "front",
	tx = 1.e-12, ty = 1.e-12,
	bg = {1,1,1},
	width = 800, height = 800);

  squares ("cs", n = {0,0,1}, alpha = 1.e-12, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 1.e-12);
  sprintf (name2, "vof-t-%.0f.png", t/(tref));
  save (name2);

  squares ("n2u", n = {0,0,1}, alpha = 1.e-12, min = 0, max = sn2u.max, map = jet);
  sprintf (name2, "nu-t-%.0f.png", t/(tref));
  save (name2);
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", n = {0,0,1}, alpha = 1.e-12, map = cool_warm);
  sprintf (name2, "ux-t-%.0f.png", t/(tref));
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", n = {0,0,1}, alpha = 1.e-12, map = cool_warm);
  sprintf (name2, "uy-t-%.0f.png", t/(tref));
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", n = {0,0,1}, alpha = 1.e-12, map = cool_warm);
  sprintf (name2, "p-t-%.0f.png", t/(tref));
  save (name2);

  /**
  We then zoom on the particle. */

  view (fov = 5.*(128./((L0)/(d))), camera = "front",
	tx = -(p_p.x + 10.*(d))/(L0), ty = 1.e-12,
	bg = {1,1,1},
	width = 800, height = 800);

  squares ("cs", n = {0,0,1}, alpha = 1.e-12, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 1.e-12);
  sprintf (name2, "vof-xy-t-%.0f.png", t/(tref));
  save (name2);
  
  squares ("n2u", n = {0,0,1}, alpha = 1.e-12, min = 0, max = sn2u.max, map = jet);
  sprintf (name2, "nu-xy-t-%.0f.png", t/(tref));
  save (name2);
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", n = {0,0,1}, alpha = 1.e-12, map = cool_warm);
  sprintf (name2, "ux-xy-t-%.0f.png", t/(tref));
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", n = {0,0,1}, alpha = 1.e-12, map = cool_warm);
  sprintf (name2, "uy-xy-t-%.0f.png", t/(tref));
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", n = {0,0,1}, alpha = 1.e-12, map = cool_warm);
  sprintf (name2, "p-xy-t-%.0f.png", t/(tref));
  save (name2);

  /**
  Here we compute two new fields, $\lambda_2$ and the vorticity
  component in the $y-z$ plane. */
  
  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  boundary ({vyz});
  lambda2 (u, l2);

  view (fov = 12.*(128/((L0)/(d))),
	quat = {0.072072,0.245086,0.303106,0.918076},
	tx = -0.307321, ty = 0.22653, bg = {1,1,1},
	width = 802, height = 634);
  draw_vof ("cs", "fs", lw = 5, lc = {0,0,0}, fc = {0,0,0});
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  sprintf (name2, "l2-t-%.0f.png", t/(tref));
  save (name2);
}

/**
## Results for a sphere in a Stokes flow

#### Domain convergence

~~~gnuplot Time evolution of the drag coefficient $C_D$
reset
set terminal svg font ",16"
set key top right font ",12" spacing 0.9
set grid ytics
set xtics 0,20,200
set ytics 0,0.1,100
set xlabel 't/(d_{eq}/u_{inf})'
set ylabel 'kappa'
set xrange [0:200]
set yrange [0.9:1.25]
set title "gamma=1, dtmax=1.e-2, cmax=1e-2"

# Stokes drag
Re = 1;
CDStokes = 24/Re;
kStokes = 1.;

# Asymptotic drag

f11(x) = a11 + b11/(x**(c11)); # L0 = 32, lmax = 9
f12(x) = a12 + b12/(x**(c12)); # L0 = 32, lmax = 10
f13(x) = a13 + b13/(x**(c13)); # L0 = 32, lmax = 11

f21(x) = a21 + b21/(x**(c21)); # L0 = 64, lmax = 10
f22(x) = a22 + b22/(x**(c22)); # L0 = 64, lmax = 11
f23(x) = a23 + b23/(x**(c23)); # L0 = 64, lmax = 12

f31(x) = a31 + b31/(x**(c31)); # L0 = 128, lmax = 11
f32(x) = a32 + b32/(x**(c32)); # L0 = 128, lmax = 12
f33(x) = a33 + b33/(x**(c33)); # L0 = 128, lmax = 13

f41(x) = a41 + b41/(x**(c41)); # L0 = 256, lmax = 12
f42(x) = a42 + b42/(x**(c42)); # L0 = 256, lmax = 13
f43(x) = a43 + b43/(x**(c43)); # L0 = 256, lmax = 14

set fit errorvariables

a11 = kStokes; a12 = kStokes; a13 = kStokes;
fit [25:*][0.5:*] f11(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-32/lmax-9/cmax-10e-3/log'   u 10:($22/CDStokes) via a11,b11,c11; # L0 = 32, lmax = 9
fit [25:*][0.5:*] f12(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-32/lmax-10/cmax-10e-3/log'  u 10:($22/CDStokes) via a12,b12,c12; # L0 = 32, lmax = 10
fit [25:*][0.5:*] f13(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-32/lmax-11/cmax-10e-3/log'  u 10:($22/CDStokes) via a13,b13,c13; # L0 = 32, lmax = 11

a21 = kStokes; a22 = kStokes; a23 = kStokes;
fit [25:*][0.5:*] f21(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-64/lmax-10/cmax-10e-3/log'  u 10:($22/CDStokes) via a21,b21,c21; # L0 = 64, lmax = 10
fit [25:*][0.5:*] f22(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-64/lmax-11/cmax-10e-3/log'  u 10:($22/CDStokes) via a22,b22,c22; # L0 = 64, lmax = 11
fit [25:*][0.5:*] f23(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-64/lmax-12/cmax-10e-3/log'  u 10:($22/CDStokes) via a23,b23,c23; # L0 = 64, lmax = 12

a31 = kStokes; a32 = kStokes; a33 = kStokes;
fit [25:*][0.5:*] f31(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-128/lmax-11/cmax-10e-3/log' u 10:($22/CDStokes) via a31,b31,c31; # L0 = 128, lmax = 11
fit [25:*][0.5:*] f32(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-128/lmax-12/cmax-10e-3/log' u 10:($22/CDStokes) via a32,b32,c32; # L0 = 128, lmax = 12
fit [25:*][0.5:*] f33(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-128/lmax-13/cmax-10e-3/log' u 10:($22/CDStokes) via a33,b33,c33; # L0 = 128, lmax = 13

a41 = kStokes; a42 = kStokes; a43 = kStokes;
fit [25:*][0.5:*] f41(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log' u 10:($22/CDStokes) via a41,b41,c41; # L0 = 256, lmax = 12
fit [25:*][0.5:*] f42(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-256/lmax-13/cmax-10e-3/log' u 10:($22/CDStokes) via a42,b42,c42; # L0 = 256, lmax = 13
fit [25:*][0.5:*] f43(x) 'gamma-1/phi-0/Re-0/dtmax-10e-3/L0-256/lmax-14/cmax-10e-3/log' u 10:($22/CDStokes) via a43,b43,c43; # L0 = 256, lmax = 14

set print "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-32-cmax-10e-3-fit.dat"
print 1./32*2**9  , a11 , a11_err , 100.*abs (a11 - kStokes)/kStokes
print 1./32*2**10 , a12 , a12_err , 100.*abs (a12 - kStokes)/kStokes
print 1./32*2**11 , a13 , a13_err , 100.*abs (a13 - kStokes)/kStokes
set print

set print "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-64-cmax-10e-3-fit.dat"
print 1./64*2**10 , a21 , a21_err , 100.*abs (a21 - kStokes)/kStokes
print 1./64*2**11 , a22 , a22_err , 100.*abs (a22 - kStokes)/kStokes
print 1./64*2**12 , a23 , a23_err , 100.*abs (a23 - kStokes)/kStokes
set print

set print "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-128-cmax-10e-3-fit.dat"
print 1./128*2**11 , a31 , a31_err , 100.*abs (a31 - kStokes)/kStokes
print 1./128*2**12 , a32 , a32_err , 100.*abs (a32 - kStokes)/kStokes
print 1./128*2**13 , a33 , a33_err , 100.*abs (a33 - kStokes)/kStokes
set print

set print "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-256-cmax-10e-3-fit.dat"
print 1./256*2**12 , a41 , a41_err , 100.*abs (a41 - kStokes)/kStokes
print 1./256*2**13 , a42 , a42_err , 100.*abs (a42 - kStokes)/kStokes
print 1./256*2**14 , a43 , a43_err , 100.*abs (a43 - kStokes)/kStokes
set print

plot kStokes w l lw 2 lc rgb "brown"  t "Stokes, 1851",		\
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-32/lmax-9/cmax-10e-3/log"   u 10:($22/CDStokes) w lp pt 6 ps 1.25 pi 800 lw 2 lc rgb "black"      t "L0/d=32, lmax=9", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-32/lmax-10/cmax-10e-3/log"  u 10:($22/CDStokes) w lp pt 3 ps 1.25 pi 800 lw 2 lc rgb "black"      t "lmax=10", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-32/lmax-11/cmax-10e-3/log"  u 10:($22/CDStokes) w lp pt 8 ps 1.25 pi 800 lw 2 lc rgb "black"      t "lmax=11", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-64/lmax-10/cmax-10e-3/log"  u 10:($22/CDStokes) w lp pt 6 ps 1.25 pi 800 lw 2 lc rgb "blue"       t "L0/d=64, lmax=10", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-64/lmax-11/cmax-10e-3/log"  u 10:($22/CDStokes) w lp pt 3 ps 1.25 pi 800 lw 2 lc rgb "blue"       t "lmax=11", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-64/lmax-12/cmax-10e-3/log"  u 10:($22/CDStokes) w lp pt 8 ps 1.25 pi 800 lw 2 lc rgb "blue"       t "lmax=12", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-128/lmax-11/cmax-10e-3/log" u 10:($22/CDStokes) w lp pt 6 ps 1.25 pi 800 lw 2 lc rgb "red"        t "L0/d=128, lmax=11", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-128/lmax-12/cmax-10e-3/log" u 10:($22/CDStokes) w lp pt 3 ps 1.25 pi 800 lw 2 lc rgb "red"        t "lmax=12", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-128/lmax-13/cmax-10e-3/log" u 10:($22/CDStokes) w lp pt 8 ps 1.25 pi 800 lw 2 lc rgb "red"        t "lmax=13", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log" u 10:($22/CDStokes) w lp pt 6 ps 1.25 pi 800 lw 2 lc rgb "sea-green"  t "L0/d=256, lmax=12", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-256/lmax-13/cmax-10e-3/log" u 10:($22/CDStokes) w lp pt 3 ps 1.25 pi 800 lw 2 lc rgb "sea-green"  t "lmax=13", \
     "gamma-1/phi-0/Re-0/dtmax-10e-3/L0-256/lmax-14/cmax-10e-3/log" u 10:($22/CDStokes) w lp pt 8 ps 1.25 pi 800 lw 2 lc rgb "sea-green"  t "lmax=14"
~~~

~~~gnuplot Drag coefficient $C_D$
set xtics 2,2,512
set xlabel 'd/Delta_{max}'
set xrange [8:128]
set yrange [0.9:1.15]
set title "gamma=1, dtmax=1.e-2, cmax=1e-2"
set logscale x
plot kStokes w l lw 2 lc rgb "brown" t "Stokes, 1851", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-32-cmax-10e-3.dat"      u 6:($8/CDStokes) w lp pt 6 ps 1.25 lw 2 lc rgb "black"     t "L_0/d=32", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-32-cmax-10e-3-fit.dat"  u 1:2 w lp pt 3 ps 1.25 lw 2 lc rgb "black"                 t "Asymptotic fit", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-64-cmax-10e-3.dat"      u 6:($8/CDStokes) w lp pt 6 ps 1.25 lw 2 lc rgb "blue"      t "L_0/d=64", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-64-cmax-10e-3-fit.dat"  u 1:2 w lp pt 3 ps 1.25 lw 2 lc rgb "blue"                  t "Asymptotic fit", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-128-cmax-10e-3.dat"     u 6:($8/CDStokes) w lp pt 6 ps 1.25 lw 2 lc rgb "red"       t "L_0/d=128", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-128-cmax-10e-3-fit.dat" u 1:2 w lp pt 3 ps 1.25 lw 2 lc rgb "red"                   t "Asymptotic fit", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-256-cmax-10e-3.dat"     u 6:($8/CDStokes) w lp pt 6 ps 1.25 lw 2 lc rgb "sea-green" t "L_0/d=256", \
     "gamma-1-phi-0-Re-0-dtmax-10e-3-L0-256-cmax-10e-3-fit.dat" u 1:2 w lp pt 3 ps 1.25 lw 2 lc rgb "sea-green"             t "Asymptotic fit"
~~~
*/
