/**
# Bending force for front-tracking membranes


*/
#define _BENDING_FT 1
#ifndef E_B
  #define E_B 1.
#endif
#ifndef REF_CURV
  #define REF_CURV 1
#endif
#ifndef GLOBAL_REF_CURV
  #if REF_CURV
    #define GLOBAL_REF_CURV 1
  #else
    #define GLOBAL_REF_CURV 0
  #endif
#endif
#if GLOBAL_REF_CURV
  #ifndef C0
  /** $c_0/a = -2.09$ is the typical reference curvature of a red blood cell. */
    // #define C0 (-2.09/RADIUS)
  #endif
#endif
#ifndef LINEAR_BENDING
  #define LINEAR_BENDING 0
#endif

#include "curvature-ft.h"


event acceleration (i++) {
/*Compute borders of the current proc*/
compute_proc_borders(&proc_max, &proc_min);

  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) 
    {
      lagMesh* mesh = &(CAPS(i));
      if(is_capsule_in_boundingbox(proc_max, proc_min, &CAPS(i))) 
      {
        comp_curvature(mesh);
        for(int j=0; j<mesh->nln; j++) {
          #if (!LINEAR_BENDING)
            double curv = mesh->nodes[j].curv;
            double rcurv = mesh->nodes[i].ref_curv;
            double gcurv = mesh->nodes[j].gcurv;
          #endif
          double lbcurv = laplace_beltrami(mesh, j, true);
          // E_B=(ND_EB*E_S*sq(RADIUS))
          double Eb = ND_EB * mesh->cap_es * sq(mesh->cap_radius);
          // printf("ND_EB: %g, Eb: %lf, es: %g, radius: %g\n", ND_EB, Eb, mesh->cap_es, mesh->cap_radius);
          double bending_surface_force = 2*Eb*(lbcurv
            #if (!LINEAR_BENDING)
              + 2*(curv - rcurv)*(sq(curv) - gcurv + rcurv*curv)
            #endif
            );
          /** We now have to compute the area associated with each node */
          double area = compute_node_area(mesh, j);
          /** The bending force is ready to be added to the Lagrangian force of
          the considered node. */
          foreach_dimension()
            mesh->nodes[j].lagForce.x +=
              mesh->nodes[j].normal.x*bending_surface_force*area;
        }
      }
    }
  }
}
