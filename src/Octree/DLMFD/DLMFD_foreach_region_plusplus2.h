#if dimension == 1
#elif dimension == 2
#define INT_CONT_CRD(crd, min, max)                                            \
  ((min.x <= crd.x && crd.x <= max.x) && (min.y <= crd.y && crd.y <= max.y))
#else // dimension == 3
#define INT_CONT_CRD(crd, min, max)                                            \
  ((min.x <= crd.x && crd.x <= max.x) && (min.y <= crd.y && crd.y <= max.y) && \
   (min.z <= crd.z && crd.z <= max.z))
#endif

#if TREE
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
#endif
#else // !TREE
#if _MPI
#if dimension == 1
#elif dimension == 2
#define PNT_MAX(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)ceil(((crd).x - X0) * _invdlt + GHOSTS -                    \
                        mpi_coords[0] * (1 << lvl));                           \
    (pnt).j = (int)ceil(((crd).y - Y0) * _invdlt + GHOSTS -                    \
                        mpi_coords[1] * (1 << lvl));                           \
    (pnt).n.x = (1 << lvl);                                                    \
    (pnt).n.y = (1 << lvl);                                                    \
  }
#define PNT_MIN(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)floor(((crd).x - X0) * _invdlt + GHOSTS -                   \
                         mpi_coords[0] * (1 << lvl));                          \
    (pnt).j = (int)floor(((crd).y - Y0) * _invdlt + GHOSTS -                   \
                         mpi_coords[1] * (1 << lvl));                          \
    (pnt).n.x = (1 << lvl);                                                    \
    (pnt).n.y = (1 << lvl);                                                    \
  }
#define CELL_CEN(pnt, crd)                                                     \
  {                                                                            \
    double _dlt = L0 / ldexp((double)Dimensions_scale, (pnt.level));           \
    crd = (coord){                                                             \
        X0 + (pnt.i - GHOSTS + mpi_coords[0] * (1 << pnt.level) + 0.5) * _dlt, \
        Y0 + (pnt.j - GHOSTS + mpi_coords[1] * (1 << pnt.level) + 0.5) * _dlt, \
        0};                                                                    \
  }
#else // dimension == 3
#define PNT_MAX(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)ceil(((crd).x - X0) * _invdlt + GHOSTS -                    \
                        mpi_coords[0] * (1 << lvl));                           \
    (pnt).j = (int)ceil(((crd).y - Y0) * _invdlt + GHOSTS -                    \
                        mpi_coords[1] * (1 << lvl));                           \
    (pnt).k = (int)ceil(((crd).z - Z0) * _invdlt + GHOSTS -                    \
                        mpi_coords[2] * (1 << lvl));                           \
    (pnt).n.x = (1 << lvl);                                                    \
    (pnt).n.y = (1 << lvl);                                                    \
    (pnt).n.z = (1 << lvl);                                                    \
  }
#define PNT_MIN(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)floor(((crd).x - X0) * _invdlt + GHOSTS -                   \
                         mpi_coords[0] * (1 << lvl));                          \
    (pnt).j = (int)floor(((crd).y - Y0) * _invdlt + GHOSTS -                   \
                         mpi_coords[1] * (1 << lvl));                          \
    (pnt).k = (int)floor(((crd).z - Z0) * _invdlt + GHOSTS -                   \
                         mpi_coords[2] * (1 << lvl));                          \
    (pnt).n.x = (1 << lvl);                                                    \
    (pnt).n.y = (1 << lvl);                                                    \
    (pnt).n.z = (1 << lvl);                                                    \
  }

#define CELL_CEN(pnt, crd)                                                     \
  {                                                                            \
    double _dlt = L0 / ldexp((double)Dimensions_scale, (pnt.level));           \
    crd = (coord){                                                             \
        X0 + (pnt.i - GHOSTS + mpi_coords[0] * (1 << pnt.level) + 0.5) * _dlt, \
        Y0 + (pnt.j - GHOSTS + mpi_coords[1] * (1 << pnt.level) + 0.5) * _dlt, \
        Z0 +                                                                   \
            (pnt.k - GHOSTS + mpi_coords[2] * (1 << pnt.level) + 0.5) * _dlt}; \
  }
#endif
#else // !_MPI
#if dimension == 1
#elif dimension == 2
#define PNT_MAX(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)ceil(((crd).x - X0) * _invdlt + GHOSTS);                    \
    (pnt).j = (int)ceil(((crd).y - Y0) * _invdlt + GHOSTS);                    \
    (pnt).n.x = (1 << lvl) * Dimensions.x;                                     \
    (pnt).n.y = (1 << lvl) * Dimensions.y;                                     \
  }
#define PNT_MIN(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)floor(((crd).x - X0) * _invdlt + GHOSTS);                   \
    (pnt).j = (int)floor(((crd).y - Y0) * _invdlt + GHOSTS);                   \
    (pnt).n.x = (1 << lvl) * Dimensions.x;                                     \
    (pnt).n.y = (1 << lvl) * Dimensions.y;                                     \
  }
#define CELL_CEN(pnt, crd)                                                     \
  {                                                                            \
    double _dlt = L0 / ldexp((double)Dimensions_scale, (pnt.level));           \
    crd = (coord){X0 + (pnt.i - GHOSTS + 0.5) * _dlt,                          \
                  Y0 + (pnt.j - GHOSTS + 0.5) * _dlt, 0};                      \
  }
#else // dimension == 3
#define PNT_MAX(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)ceil(((crd).x - X0) * _invdlt + GHOSTS);                    \
    (pnt).j = (int)ceil(((crd).y - Y0) * _invdlt + GHOSTS);                    \
    (pnt).k = (int)ceil(((crd).z - Z0) * _invdlt + GHOSTS);                    \
    (pnt).n.x = (1 << lvl) * Dimensions.x;                                     \
    (pnt).n.y = (1 << lvl) * Dimensions.y;                                     \
    (pnt).n.z = (1 << lvl) * Dimensions.z;                                     \
  }
#define PNT_MIN(pnt, crd, lvl)                                                 \
  {                                                                            \
    double _invdlt = ldexp((double)Dimensions_scale, (lvl)) / L0;              \
    (pnt).i = (int)floor(((crd).x - X0) * _invdlt + GHOSTS);                   \
    (pnt).j = (int)floor(((crd).y - Y0) * _invdlt + GHOSTS);                   \
    (pnt).k = (int)floor(((crd).z - Z0) * _invdlt + GHOSTS);                   \
    (pnt).n.x = (1 << lvl) * Dimensions.x;                                     \
    (pnt).n.y = (1 << lvl) * Dimensions.y;                                     \
    (pnt).n.z = (1 << lvl) * Dimensions.z;                                     \
  }
#define CELL_CEN(pnt, crd)                                                     \
  {                                                                            \
    double _dlt = L0 / ldexp((double)Dimensions_scale, (pnt.level));           \
    crd = (coord){X0 + (pnt.i - GHOSTS + 0.5) * _dlt,                          \
                  Y0 + (pnt.j - GHOSTS + 0.5) * _dlt,                          \
                  Z0 + (pnt.k - GHOSTS + 0.5) * _dlt};                         \
  }
#endif
#endif
#endif

macro2 foreach_region_plus_plus(coord rmin, coord rmax) 
{
  {
    int ig = 0;
    int jg = 0;
    int kg = 0;

    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    NOT_UNUSED(kg);

    Point point = {0};
    Point pmin = {0}, pmax = {0};

    for (point.level = depth(); point.level >= 0; point.level--) {

      PNT_MIN(pmin, rmin, point.level);
      PNT_MAX(pmax, rmax, point.level);

#if !TREE
      foreach_dimension() { point.n.x = pnt.n.x; }
#endif

      for (point.i = pmin.i; point.i < pmax.i; point.i++)
#if dimension >= 2
        for (point.j = pmin.j; point.j < pmax.j; point.j++)
#endif
#if dimension >= 3
          for (point.k = pmin.k; point.k < pmax.k; point.k++)
#endif

#if TREE
            if (allocated(0, 0, 0) && is_local(cell))
#else // !TREE
        if (point.i >= GHOSTS && point.i < point.n.x + GHOSTS)
#if dimension >= 2
          if (point.j >= GHOSTS && point.j < point.n.y + GHOSTS)
#endif
#if dimension >= 3
            if (point.k >= GHOSTS && point.k < point.n.z + GHOSTS)
#endif
#endif
            {
              coord centre = {0};
              CELL_CEN(point, centre);
              if (INT_CONT_CRD(centre, rmin, rmax)) {
                // clang-format off
                {...}
                // clang-format on
              }
            }
    }
  }
}
