#if dimension == 1
#elif dimension == 2
#define PNT_MAX(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp(1.0, (lvl)) / L0;                                   \
    (pnt).i = (int)ceil(((crd).x - X0) * _invdlt + GHOSTS);                    \
    (pnt).j = (int)ceil(((crd).y - Y0) * _invdlt + GHOSTS);                    \
  }
#define PNT_MIN(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp(1.0, (lvl)) / L0;                                   \
    (pnt).i = (int)floor(((crd).x - X0) * _invdlt + GHOSTS);                   \
    (pnt).j = (int)floor(((crd).y - Y0) * _invdlt + GHOSTS);                   \
  }
#define CELL_CEN(pnt, crd)                                                     \
  {                                                                            \
    double _dlt = L0 / ldexp(1.0, (pnt.level));                                \
    crd = (coord){X0 + (pnt.i - GHOSTS + 0.5) * _dlt,                          \
                  Y0 + (pnt.j - GHOSTS + 0.5) * _dlt, 0};                      \
  }
#define INT_CONT_CRD(crd, min, max)                                            \
  ((min.x <= crd.x && crd.x <= max.x) && (min.y <= crd.y && crd.y <= max.y))

#else // dimension == 3
#define PNT_MAX(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp(1.0, (lvl)) / L0;                                   \
    (pnt).i = (int)ceil(((crd).x - X0) * _invdlt + GHOSTS);                    \
    (pnt).j = (int)ceil(((crd).y - Y0) * _invdlt + GHOSTS);                    \
    (pnt).k = (int)ceil(((crd).z - Z0) * _invdlt + GHOSTS);                    \
  }
#define PNT_MIN(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp(1.0, (lvl)) / L0;                                   \
    (pnt).i = (int)floor(((crd).x - X0) * _invdlt + GHOSTS);                   \
    (pnt).j = (int)floor(((crd).y - Y0) * _invdlt + GHOSTS);                   \
    (pnt).k = (int)floor(((crd).z - Z0) * _invdlt + GHOSTS);                   \
  }

#define CELL_CEN(pnt, crd)                                                     \
  {                                                                            \
    double _dlt = L0 / ldexp(1.0, (pnt.level));                                \
    crd = (coord){X0 + (pnt.i - GHOSTS + 0.5) * _dlt,                          \
                  Y0 + (pnt.j - GHOSTS + 0.5) * _dlt,                          \
                  Z0 + (pnt.k - GHOSTS + 0.5) * _dlt};                         \
  }
#define INT_CONT_CRD(crd, min, max)                                            \
  ((min.x <= crd.x && crd.x <= max.x) && (min.y <= crd.y && crd.y <= max.y) && \
   (min.z <= crd.z && crd.z <= max.z))
#endif


macro2 foreach_region_plus_plus(coord rmin, coord rmax) {
  {
    int ig = 0;
    int jg = 0;
    int kg = 0;

    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    NOT_UNUSED(kg);

    Point point = {0};
    Point pmin = {0}, pmax = {0};

#   if TREE // Tree version

    for (point.level = depth(); point.level >= 0; point.level--) {
      PNT_MIN(pmin, rmin, point.level);
      PNT_MAX(pmax, rmax, point.level);

      // bool no_leaves = true;

#if dimension == 1
      for (point.i = pmin.i; point.i < pmax.i; point.i++) {
        if (allocated(0, 0, 0) && is_local(cell)) {
          coord centre = {0};
          CELL_CEN(point, centre);
          if (INT_CONT_CRD(centre, rmin, rmax)) {
            // clang-format off
            {...}
            // clang-format on
            // if (is_leaf(cell))
            //   no_leaves = false;
          }
        }
      }
#elif dimension == 2
      for (point.i = pmin.i; point.i < pmax.i; point.i++) {
        for (point.j = pmin.j; point.j < pmax.j; point.j++) {
          if (allocated(0, 0, 0) && is_local(cell)) {
            coord centre = {0};
            CELL_CEN(point, centre);
            if (INT_CONT_CRD(centre, rmin, rmax)) {
              // clang-format off
              {...}
              // clang-format on
              // if (is_leaf(cell))
              //   no_leaves = false;
            }
          }
        }
      }
#else
      for (point.i = pmin.i; point.i < pmax.i; point.i++) {
        for (point.j = pmin.j; point.j < pmax.j; point.j++) {
          for (point.k = pmin.k; point.k < pmax.k; point.k++) {
            if (allocated(0, 0, 0) && is_local(cell)) {
              coord centre = {0};
              CELL_CEN(point, centre);
              if (INT_CONT_CRD(centre, rmin, rmax)) {
                // clang-format off
                {...}
                // clang-format on
                // if (is_leaf(cell))
                //   no_leaves = false;
              }
            }
          }
        }
      }
#endif
      // if (no_leaves)
      //   break;
    }

# else // Multigrid version   
    coord gridmin = {0}, gridmax = {0}, intersectmin = {0}, intersectmax = {0};

#   if  _MPI    
      // Local grid min/max coordinates
      gridmin.x = ( (double)mpi_coords[0] / (double)Dimensions.x ) * L0 + X0;
      gridmax.x = ( (double)( mpi_coords[0] + 1 ) / (double)Dimensions.x ) * L0 
    		+ X0;      
#     if dimension >= 2    
        gridmin.y = ( (double)mpi_coords[1] / (double)Dimensions.x ) * L0 + Y0;
        gridmax.y = ( (double)( mpi_coords[1] + 1 ) / (double)Dimensions.x ) 
		* L0 + Y0;
#     endif
#     if dimension >= 3  		  	    
        gridmin.z = ( (double)mpi_coords[2] / (double)Dimensions.x ) * L0 + Z0;
        gridmax.z = ( (double)( mpi_coords[2] + 1 ) / (double)Dimensions.x ) 
		* L0 + Z0;
#     endif		 
#   else
      gridmin.x = X0;
      gridmax.x = X0 + L0;
#     if dimension >= 2    
        gridmin.y = Y0;
        gridmax.y = Y0 + L0 * (double)Dimensions.y / (double)Dimensions.x;
#     endif
#     if dimension >= 3  		  	    
        gridmin.z = Z0;      
        gridmax.z = Z0 + L0 * (double)Dimensions.z / (double)Dimensions.x;
#     endif	      
#   endif

    // Do the local grid and region AABB intersect ?
    // If they do, compute the intersection box min/max coordinates
#   if dimension == 1
      if ( rmin.x <= gridmax.x && rmax.x >= gridmin.x ) 
#   elif dimension == 2
      if ( ( rmin.x <= gridmax.x && rmax.x >= gridmin.x ) &&
	( rmin.y <= gridmax.y && rmax.y >= gridmin.y ) ) 
#   else // dimension == 3 
      if ( ( rmin.x <= gridmax.x && rmax.x >= gridmin.x ) &&
	( rmin.y <= gridmax.y && rmax.y >= gridmin.y ) &&
	( rmin.z <= gridmax.z && rmax.z >= gridmin.z ) )
#   endif 	 
      {
        foreach_dimension() 
        {
          // Note the +/- 1.e-8 * L0 avoids round off errors
	  intersectmin.x = max( rmin.x, gridmin.x ) + 1.e-8 * L0;
          intersectmax.x = min( rmax.x, gridmax.x ) - 1.e-8 * L0;        
        }
      }              

    pmin = locate( intersectmin.x, intersectmin.y, intersectmin.z );
    pmax = locate( intersectmax.x, intersectmax.y, intersectmax.z );
    point.level = depth();
    foreach_dimension() point.n.x = pmin.n.x; // is this really needed ?
    
#if dimension == 1
      for (point.i = pmin.i; point.i <= pmax.i; point.i++) {
        // clang-format off
        {...}
        // clang-format on
        // if (is_leaf(cell))
        //   no_leaves = false;
      }
#elif dimension == 2
      for (point.i = pmin.i; point.i <= pmax.i; point.i++) {
        for (point.j = pmin.j; point.j <= pmax.j; point.j++) {
          // clang-format off
          {...}
          // clang-format on
          // if (is_leaf(cell))
          //   no_leaves = false;
        }
      }
#else
      for (point.i = pmin.i; point.i <= pmax.i; point.i++) {
        for (point.j = pmin.j; point.j <= pmax.j; point.j++) {
          for (point.k = pmin.k; point.k <= pmax.k; point.k++) {
                {...}
                // clang-format on
                // if (is_leaf(cell))
                //   no_leaves = false;
          }
        }
      }
#endif

#endif // Tree or multigrid    

  }
}
