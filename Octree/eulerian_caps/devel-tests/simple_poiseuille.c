#ifndef LEVEL
  #define LEVEL 6
#endif
#ifndef SIZE
  #define SIZE 10
#endif
#ifndef RADIUS
  #define RADIUS 1.
#endif
#ifndef RHO
  #define RHO 1.
#endif
#ifndef MU
  #define MU 1.
#endif
#ifndef MU_RATIO
  #define MU_RATIO 5.
#endif
#ifndef U_MAX
  #define U_MAX 1.
#endif
#ifndef GRAD_P
  #define GRAD_P (4*MU*U_MAX/(SIZE*SIZE))
#endif
#ifndef CA
  #define CA 0.2
#endif
#ifndef G_SHEAR
  #define G_SHEAR (MU*U_MAX/CA)
#endif
#ifndef ND_EB
  #define ND_EB 0.0005
#endif
#ifndef E_BEND
  #define E_BEND (ND_EB*RADIUS*RADIUS*G_SHEAR)
#endif
#ifndef REF_CURV
  #define REF_CURV 1
#endif
#ifndef T_END
  #define T_END 100.
#endif
#ifndef DT_MAX
  #define DT_MAX (0.1)
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-10
#endif
#ifndef U_TOL
  #define U_TOL (1.e-2)
#endif
#ifndef OUTPUT_FREQ
  // #define OUTPUT_FREQ (t += 0.01)
  #define OUTPUT_FREQ (t += 0.5)
#endif
#ifndef AMR
  #define AMR 0
#endif

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "view.h"

FILE * fp = NULL;
scalar noise_gamma[], sdivu[];

int main() {
  periodic(right);
  L0 = SIZE;
  origin(-.5*L0, -.5*L0);
  stokes = true;
  // rho = RHO;
  const face vector muc[] = {MU,MU,MU};
  mu = muc;
  N = 1 << LEVEL;
  TOLERANCE = MY_TOLERANCE;
  DT = DT_MAX;
  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.t[top] = dirichlet(0);
// p[left] = dirichlet(-GRAD_P*SIZE/2);
// p[right] = dirichlet(GRAD_P*SIZE/2);

event init (i=0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

event acceleration (i++) {
  face vector ap = a;
  foreach_face(x) {
    ap.x[] += GRAD_P;
  }
}

event movie (OUTPUT_FREQ) {
  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("u.x", linear = false);
  box (notics = true);
  save ("ux_cells.mp4");
}

// event end (t = T_END) {
event end (t = 100) {
  return 1.;
}
