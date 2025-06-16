/**
# Bell-Collela-Glaz advection scheme
*/

foreach_dimension()
double get_sf_x( Point point, scalar f, int i, int j ) {
	double result = 0;
	if (BIGAMMA(i,j)) {
		result = f[i,j];
	}
	else {
		// accessed_cells[i,j] = 1.;
		double grad_f = 0.;
		int sn = (i < 0) ? 1 : -1;
		int st = (j < 0) ? 1 : -1;
		int ext;
		bool normal = (i > 0) ? BIGAMMA(i + sn*3, j) : BIGAMMA(i + sn*2, j);
		bool tangential = (j > 0) ? BIGAMMA(i, j + st*3) : BIGAMMA(i, j + st*2);
		if (normal) {
			ext = (BIGAMMA(i + sn, j)) ? 0 : 1;
			if ((ext) && !(BIGAMMA(i + 2*sn, j))) ext++;
			if (!(is_boundary(point) && ext>0)) {
				grad_f = (f[i + sn*(1 + ext), j] - f[i + sn*(2 + ext), j])/Delta;
				result = f[i + sn*(1 + ext), j] + grad_f * (1 + ext)*Delta;
				// accessed_cells[i + sn*(1 + ext), j] = 1.;
				// accessed_cells[i + sn*(2 + ext), j] = 1.;
			}
			else normal = 0.;
		}
		if (tangential) {
			ext = (BIGAMMA(i, j + st)) ? 0 : 1;
			if ((ext) && !(BIGAMMA(i, j + 2*st))) ext++;
			if (!(is_boundary(point) && ext>0)) {
				grad_f = (f[i, j + st*(1 + ext)] - f[i, j + st*(2 + ext)])/Delta;
				result = f[i, j + st*(1 + ext)] + grad_f * (1 + ext)*Delta;
				// accessed_cells[i, j + st*(1 + ext)] = 1.;
				// accessed_cells[i, j + st*(2 + ext)] = 1.;
			}
			else tangential = 0.;
		}
		if (normal && tangential) result /= 2.;
		else if (!(normal) && !(tangential)) result = f[];

	}
	return result;
}

/**
The function below implements the 2nd-order, unsplit, upwind scheme of
[Bell-Collela-Glaz, 1989](references.bib#bell89). Given a centered
scalar field *f*, a face vector field *uf* (possibly weighted by a
face metric), a timestep *dt* and a source term field *src*, it fills
the face vector field *flux* with the components of the advection
fluxes of *f*. */

void my_tracer_fluxes (scalar f,
		    face vector uf,
		    face vector flux,
		    double dt,
		    (const) scalar src)
{
  /**
  For each face, the flux is composed of two parts... */

  foreach_face() {
		if (BIGAMMA(0,0) || BIGAMMA(-1,0)) {
	    /**
	    A normal component... (Note that we cheat a bit here, `un` should
	    strictly be `dt*(uf.x[i] + uf.x[i+1])/((fm.x[] +
	    fm.x[i+1])*Delta)` but this causes trouble with boundary
	    conditions (when using narrow '1 ghost cell' stencils)). */

	    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
	    int i = -(s + 1.)/2.;
			double grad_f = (un > 0.) ? (get_sf_x(point, f, 0, 0)
																	- get_sf_x(point,f, -2, 0))/(2.*Delta)
															 : (get_sf_x(point, f, 1, 0)
															   - get_sf_x(point, f, -1, 0))/(2.*Delta);
			double f2 = f[i] + (src[] + src[-1])*dt/4.
									+ s*(1. - s*un)*grad_f*Delta/2.;

	    /**
	    and tangential components... */

	    #if dimension > 1
	    if (fm.y[i] && fm.y[i,1]) {
	      double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
				double fyy = vn < 0. ? get_sf_x(point, f, i, 1)
															 - get_sf_x(point, f, i, 0)
														 : get_sf_x(point, f, i, 0)
														 	 - get_sf_x(point, f, i, -1);
	      f2 -= dt*vn*fyy/(2.*Delta);
	    }
	    #endif
	    #if dimension > 2
	    if (fm.z[i] && fm.z[i,0,1]) {
	      double wn = (uf.z[i] + uf.z[i,0,1])/(fm.z[i] + fm.z[i,0,1]);
	      double fzz = wn < 0. ? f[i,0,1] - f[i] : f[i] - f[i,0,-1];
	      f2 -= dt*wn*fzz/(2.*Delta);
	    }
	    #endif

	    flux.x[] = f2*uf.x[];
		}
		else {
			flux.x[] = 0.;
		}
  }

  /**
  Boundary conditions ensure the consistency of fluxes across
  variable-resolution boundaries (on adaptive meshes). */

  boundary_flux ({flux});
}

/**
The function below uses the *tracer_fluxes* function to integrate the
advection equation, using an explicit scheme with timestep *dt*, for
each tracer in the list. */

void local_advection (struct Advection p)
{

  /**
  If *src* is not provided we set all the source terms to zero. */

  scalar * lsrc = p.src;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  for (f,src in p.tracers,lsrc) {
    face vector flux[];
    my_tracer_fluxes (f, p.u, flux, p.dt, src);
#if !EMBED
    foreach() {
			if (GAMMA) {
	      foreach_dimension()
	        f[] += p.dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
				}
	#else // EMBED
	    update_tracer (f, p.u, flux, p.dt);
	#endif // EMBED
		}
  }
  boundary (p.tracers);

  if (!p.src)
    free (lsrc);
}
