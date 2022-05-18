/** In this file we compute the neo-Hookean force on each Lagrangian node of
the membrane(s).*/

#ifndef E_S
  #define E_S 1.
#endif

#define DWDL1(L1, L2) (E_S*L1*(1. - 1./(sq(L1*L2)))/3.)
#define DWDL2(L1, L2) (E_S*L2*(1. - 1./(sq(L1*L2)))/3.)

#include "elasticity-ft.h"

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) comp_elastic_stress(&mbs.mb[i]);
}
