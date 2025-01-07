/**
# Eulerian representation of capsules

We wish to use the VOF framework to represent capsules,
i.e. thin mecapsranes enclosing an inner fluid. Our approach is greatly inspired
from the work of [Ii et al., 2012](https://www.cambridge.org/core/journals/communications-in-computational-physics/article/abs/full-eulerian-fluidmembrane-coupling-method-with-a-smoothed-volumeoffluid-approach/D86755460C14F598A46869E6BDAF8EC5)
and [Ii et al., 2018](https://www.sciencedirect.com/science/article/pii/S0021999118305229).

In this approach, and in contrast to the continuum surface force from Brackbill
et al., the membrane stresses are transmitted to the fluid over a zone of
influence $\Gamma$ of width greater than one mesh size (typically a width of
about 7\Delta in the direction normal to the mecapsrane is a good choice). This
is because some Eulerian quantities (needed to capture the elastic
behaviour of the membrane) are to be defined and advected in the
membrane region only.

The scalar $\phi \equiv caps$ is defined as a smoothed version of the fraction scalar $f$.
Since we have access to a 5x5 stencil but wished to have a 7x7 stencil, we apply
a smoothing twice - hence the definition of $ngcaps$. We also define the
gradient of $caps$, which is non zero in $\Gamma$ only.

*/
typedef struct { double x, y, z;}   pseudo_v;
typedef struct { pseudo_v x, y, z;} pseudo_t;

#include "two-phase.h"

#define IS_INTERFACE_CELL(point, f) (interfacial(point, f))

#ifndef GRAD_THRESHOLD
  #define GRAD_THRESHOLD 1.e-1
#endif
#ifndef SHARP_DIRAC
  #define SHARP_DIRAC 0
#endif
#ifndef RECOMPUTE_NGCAPS
  #define RECOMPUTE_NGCAPS 0
#endif
#ifndef COMPUTE_NORMALS
  #define COMPUTE_NORMALS 1
#endif
#ifndef USE_MYC
  #define USE_MYC 0
#endif
#ifndef USE_YOUNGS
  #define USE_YOUNGS 0
#endif
#if (USE_MYC) || (USE_YOUNGS)
  #define USE_HEIGHT_FUNCTIONS 0
#endif
#ifndef USE_HEIGHT_FUNCTIONS
  #define USE_HEIGHT_FUNCTIONS 1
#endif
#ifndef RENORMALIZE_EXTENDED_NORMALS
  #define RENORMALIZE_EXTENDED_NORMALS 0
#endif


#if (USE_MYC)
  #include "myc.h"
#endif

#define BGHOSTS 2
#define EPS 1.e-14
#define GAMMA (ngcaps[] > GRAD_THRESHOLD)
#define NBG_GAMMA(X,Y,Z) ((ngcaps[X,Y,Z] > GRAD_THRESHOLD))
/** Below we implement the vertex-averaged and the double-vertex-averaged macros. We use them to define the membrane region by smoothing the color function.*/
#if dimension<=2
  #define VAVG(F,X,Y) ((4.*F[X,Y] + 2.*(F[X+1,Y,0] + F[X,Y+1] + F[X-1,Y] + \
  F[X,Y-1]) + F[X+1,Y+1] + F[X-1,Y+1] + F[X+1,Y-1] + F[X-1,Y-1])/16.)
#else
  #define VAVG(F,X,Y,Z) ((8*F[X,Y,Z] + 4*(F[X+1,Y,Z] + F[X,Y+1,Z] + F[X,Y,Z+1] \
  + F[X-1,Y,Z] + F[X,Y-1,Z] + F[X,Y,Z-1]) + 2*(F[X+1,Y+1,Z] + F[X+1,Y,Z+1] + \
  F[X,Y+1,Z+1] + F[X-1,Y+1,Z] + F[X+1,Y-1,Z] + F[X-1,Y,Z+1] + F[X+1,Y,Z-1] + \
  F[X,Y+1,Z-1] + F[X,Y-1,Z+1] + F[X-1,Y-1,Z] + F[X-1,Y,Z-1] + F[X,Y-1,Z-1]) + \
  F[X+1,Y+1,Z+1] + F[X-1,Y+1,Z+1] + F[X+1,Y-1,Z+1] + F[X+1,Y+1,Z-1] + \
  F[X-1,Y-1,Z+1] + F[X-1,Y+1,Z-1] + F[X+1,Y-1,Z-1] + F[X-1,Y-1,Z-1])/64.)
#endif
#if dimension<=2
  #define DVAVG(F) ((4.*VAVG(F,0,0) + 2.*(VAVG(F,1,0) + VAVG(F,0,1) + \
  VAVG(F,-1,0) + VAVG(F,0,-1)) + VAVG(F,1,1) + VAVG(F,-1,1) + VAVG(F,1,-1) + \
  VAVG(F,-1,-1))/16.)
#else
  #define DVAVG(F) ((8*VAVG(F,0,0,0) + 4*(VAVG(F,1,0,0) + VAVG(F,0,1,0) + \
  VAVG(F,0,0,1) + VAVG(F,-1,0,0) + VAVG(F,0,-1,0) + VAVG(F,0,0,-1)) + \
  2*(VAVG(F,1,1,0) + VAVG(F,1,0,1) + VAVG(F,0,1,1) + VAVG(F,-1,1,0) + \
  VAVG(F,1,-1,0) + VAVG(F,-1,0,1) + VAVG(F,1,0,-1) + VAVG(F,0,1,-1) + \
  VAVG(F,0,-1,1) + VAVG(F,-1,-1,0) + VAVG(F,-1,0,-1) + VAVG(F,0,-1,-1)) + \
  VAVG(F,1,1,1) + VAVG(F,-1,1,1) + VAVG(F,1,-1,1) + VAVG(F,1,1,-1) + \
  VAVG(F,-1,-1,1) + VAVG(F,-1,1,-1) + VAVG(F,1,-1,-1) + VAVG(F,-1,-1,-1))/64.)
#endif

scalar caps[];
vector grad_caps[];
scalar ngcaps[];
int nb_iter_extension = 100;

symmetric tensor Ts[];
vector extended_n[];
scalar aen[]; //for debug only
vector centered_ae[]; //for debug only
#include "normal_extension.h"

event defaults (i = 0) {
  for (scalar s in {extended_n, Ts}) {
      s.v.x.i = -1;
      foreach_dimension() {
        s[left] = neumann(0.);
        s[right] = neumann(0.);
      }
  }
  foreach() {
    foreach_dimension() {
      Ts.x.x[] = 0.;
      extended_n.x[] = 0.;
    }
    Ts.x.y[] = 0.;
    Ts.x.z[] = 0.;
    Ts.y.z[] = 0.;
  }
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

event init (i = 0) {
  boundary((scalar*){Ts, extended_n});
}

event pp_vof (i++) {
  if (i>2) nb_iter_extension = MAX_ITER_EXTENSION;
  foreach() {
    caps[] = DVAVG(f);
  }
  boundary({caps});

  foreach() {
    foreach_dimension() {
      grad_caps.x[] = (caps[1,0] - caps[-1,0])/(2.*Delta);
    }
    ngcaps[] = norm(grad_caps);
  }
  boundary({ngcaps, grad_caps});

  #if (COMPUTE_NORMALS)
    #if (USE_HEIGHT_FUNCTIONS)
      if ((i == 0) && (normf(extended_n.x).max + normf(extended_n.y).max +
        normf(extended_n.z).max < .5)) {
    #endif
      foreach() {
        if (GAMMA) {
          #if (USE_MYC) || (USE_YOUNGS)
            #if (USE_MYC)
              coord my_normals = mycs(point, caps);
            #else
              coord my_normals = youngs_normals(point, caps);
            #endif
            foreach_dimension() {
              extended_n.x[] = my_normals.x;
          #else
            foreach_dimension() {
              extended_n.x[] = -grad_caps.x[]/ngcaps[];
          #endif
          }
        }
        else {
          foreach_dimension() {
            extended_n.x[] = 0.;
          }
        }
      }
      boundary((scalar*){extended_n});
    #if (USE_HEIGHT_FUNCTIONS)
        }
      else if (i > 0) {
        initialize_normals_for_extension(extended_n);
        normal_vector_extension(extended_n);
      }
    #endif
  #endif
}


/**
We are now ready to add the body force to the face acceleration field.
*/
event acceleration (i++) {
  face vector ae = a;
  foreach() f[] = clamp (f[], 0., 1.);
  foreach_face() {
    if (NBG_GAMMA(0,0,0) && NBG_GAMMA(-1,0,0)) {
      double dxtxx, dytxx, dztxx, dxtxy, dytxy, dztxy, dytxz, dztxz, dxtxz;
      double nxf, nyf, nzf, ngcapsf;
      dxtxx = Ts.x.x[] - Ts.x.x[-1];
      dxtxy = Ts.x.y[] - Ts.x.y[-1];
      dxtxz = Ts.x.z[] - Ts.x.z[-1];
      dytxx = (Ts.x.x[-1,1,0] + Ts.x.x[0,1,0] - Ts.x.x[-1,-1,0]
          - Ts.x.x[0,-1,0])/4.;
      dytxy = (Ts.x.y[-1,1,0] + Ts.x.y[0,1,0] - Ts.x.y[-1,-1,0]
          - Ts.x.y[0,-1,0])/4.;
      dytxz = (Ts.x.z[-1,1,0] + Ts.x.z[0,1,0] - Ts.x.z[-1,-1,0]
          - Ts.x.z[0,-1,0])/4.;
      dztxx = (Ts.x.x[-1,0,1] + Ts.x.x[0,0,1] - Ts.x.x[-1,0,-1]
          - Ts.x.x[0,0,-1])/4.;
      dztxy = (Ts.x.y[-1,0,1] + Ts.x.y[0,0,1] - Ts.x.y[-1,0,-1]
          - Ts.x.y[0,0,-1])/4.;
      dztxz = (Ts.x.z[-1,0,1] + Ts.x.z[0,0,1] - Ts.x.z[-1,0,-1]
          - Ts.x.z[0,0,-1])/4.;
      nxf = .5*(extended_n.x[] + extended_n.x[-1]);
      nyf = .5*(extended_n.y[] + extended_n.y[-1]);
      nzf = .5*(extended_n.z[] + extended_n.z[-1]);
      #if SHARP_DIRAC
        double dxf, dyf, dzf;
        dxf = f[] - f[-1];
        dyf = (f[-1,1,0] + f[0,1,0] - f[-1,-1,0] - f[0,-1,0])/4.;
        dzf = (f[-1,0,1] + f[0,0,1] - f[-1,0,-1] - f[0,0,-1])/4.;
        ngcapsf = sqrt(sq(dxf) + sq(dyf) + sq(dzf))/Delta;
      #else
        #if RECOMPUTE_NGCAPS
          double dxcaps, dycaps, dzcaps;
          dxcaps = caps[] - caps[-1];
          dycaps = (caps[-1,1,0] + caps[0,1,0] - caps[-1,-1,0]
            - caps[0,-1,0])/4.;
          dzcaps = (caps[-1,0,1] + caps[0,0,1] - caps[-1,0,-1]
            - caps[0,0,-1])/4.;
          ngcapsf = sqrt(sq(dxcaps) + sq(dycaps) + sq(dzcaps))/Delta;
        #else
          ngcapsf = .5*(ngcaps[] + ngcaps[-1]);
        #endif
      #endif
      ae.x[] += ( ((1 - sq(nxf))*dxtxx - nxf*nyf*dytxx - nxf*nzf*dztxx)
        + (-nxf*nyf*dxtxy + (1 - sq(nyf))*dytxy - nyf*nzf*dztxy)
        + (-nxf*nzf*dxtxz - nyf*nzf*dytxz + (1 - sq(nzf))*dztxz)
        )*(ngcapsf*alpha.x[]/Delta);
    }
  }
  boundary((scalar *){ae});
  // --- DEBUG ---
  foreach() {
    if (GAMMA) {
      foreach_dimension()
        centered_ae.x[] = (ae.x[1,0] + ae.x[])/2.;
      aen[] = sqrt( sq(centered_ae.x[]) + sq(centered_ae.y[]) + sq(centered_ae.z[]) );
    }
    else {
      foreach_dimension() centered_ae.x[] = 0.;
      aen[] = 0.;
    }
  }
  // --- END DEBUG ---
}

/**
# Test cases
* [Moving capsule without elasticity](../test-cases/moving_capsule_no_el.c)
*/
