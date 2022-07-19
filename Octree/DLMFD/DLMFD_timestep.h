// note: u is weighted by fm
double dlmfd_timestep (const face vector u, double dtmax)
{
  static double previous = 0.;
  
  static unsigned int counter = 0;
  if ( !counter ) previous = dtrestart;
  ++counter;  
  
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = Delta/fabs(u.x[]);
#if EMBED
      assert (fm.x[]);
      dt *= fm.x[];
#else
      dt *= cm[];
#endif
      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
