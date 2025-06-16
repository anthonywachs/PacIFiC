#include "range.h"
#include "vtkType.h"

typedef struct {
  // MPI-Only
  size_t n_ranges;
  range_t *ranges;
  size_t *offsets;

  // Global variables that are similar across all procs
  double *x, *y, *z;
  size_t n_x, n_y, n_z;
  int64_t *number_of_vertices_per_depth;
  size_t max_vertices, total_vertices, max_depth;
  double time;

  // Master variables that are only held on pid 0
  Bit_t *descriptors;
  size_t descriptors_size;
  size_t n_descriptors;

  // Distributed variables: These are mostly used if called by a memsource
  // function
} vtkHyperTreeGridData;

void get_coordinates(vtkHyperTreeGridData *vtk_hypertreegrid_data) {

#if dimension >= 1
  vtk_hypertreegrid_data->n_x = 2;
  vtk_hypertreegrid_data->x = malloc(2 * sizeof(double));
  vtk_hypertreegrid_data->x[0] = X0;
  vtk_hypertreegrid_data->x[1] = X0 + L0;
#else
  vtk_hypertreegrid_data->n_x = 1;
  vtk_hypertreegrid_data->x = malloc(1 * sizeof(double));
  vtk_hypertreegrid_data->x[0] = 0;
#endif

#if dimension >= 2
  vtk_hypertreegrid_data->n_y = 2;
  vtk_hypertreegrid_data->y = malloc(2 * sizeof(double));
  vtk_hypertreegrid_data->y[0] = Y0;
  vtk_hypertreegrid_data->y[1] = Y0 + L0;
#else
  vtk_hypertreegrid_data->n_x = 1;
  vtk_hypertreegrid_data->y = malloc(1 * sizeof(double));
  vtk_hypertreegrid_data->y[0] = 0;
#endif

#if dimension >= 3
  vtk_hypertreegrid_data->n_z = 2;
  vtk_hypertreegrid_data->z = malloc(2 * sizeof(double));
  vtk_hypertreegrid_data->z[0] = Z0;
  vtk_hypertreegrid_data->z[1] = Z0 + L0;
#else
  vtk_hypertreegrid_data->n_z = 1;
  vtk_hypertreegrid_data->z = malloc(1 * sizeof(double));
  vtk_hypertreegrid_data->z[0] = 0;
#endif

  // return vtk_hypertreegrid_data;
}

void get_number_of_vertices_per_depth(
    vtkHyperTreeGridData *vtk_hypertreegrid_data) {
  int depth;

#if _MPI
  // use grid->maxdepth (Basilisk’s actual field), not grid->depth
  MPI_Allreduce(&grid->maxdepth, &depth, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
  depth = grid->maxdepth;
#endif

  size_t max_depth_old = vtk_hypertreegrid_data->max_depth;
  vtk_hypertreegrid_data->max_depth = (size_t)(depth + 1);

  if (vtk_hypertreegrid_data->number_of_vertices_per_depth != NULL) {
    if (vtk_hypertreegrid_data->max_depth != max_depth_old) {
      free(vtk_hypertreegrid_data->number_of_vertices_per_depth);
      vtk_hypertreegrid_data->number_of_vertices_per_depth =
          malloc(vtk_hypertreegrid_data->max_depth * sizeof(int64_t));
    }
  } else {
    vtk_hypertreegrid_data->number_of_vertices_per_depth =
        malloc(vtk_hypertreegrid_data->max_depth * sizeof(int64_t));
  }

  memset(vtk_hypertreegrid_data->number_of_vertices_per_depth, 0,
         vtk_hypertreegrid_data->max_depth * sizeof(int64_t));

  // count your local cells per level
  foreach_cell_BFS() {
    if (is_local(cell)) {
      vtk_hypertreegrid_data->number_of_vertices_per_depth[level]++;
    }
  }

#if _MPI
  MPI_Allreduce(
      MPI_IN_PLACE, vtk_hypertreegrid_data->number_of_vertices_per_depth,
      vtk_hypertreegrid_data->max_depth, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif

  // compute the max and total
  vtk_hypertreegrid_data->max_vertices = 0;
  vtk_hypertreegrid_data->total_vertices = 0;
  for (size_t i = 0; i < vtk_hypertreegrid_data->max_depth; i++) {
    int64_t nv = vtk_hypertreegrid_data->number_of_vertices_per_depth[i];
    if (nv > vtk_hypertreegrid_data->max_vertices)
      vtk_hypertreegrid_data->max_vertices = nv;
    vtk_hypertreegrid_data->total_vertices += nv;
  }
}

void get_ranges_and_offsets(vtkHyperTreeGridData *vtk_hypertreegrid_data) {
  if (vtk_hypertreegrid_data->ranges != NULL) {
    free(vtk_hypertreegrid_data->ranges);
    vtk_hypertreegrid_data->ranges = NULL;
  }

  if (vtk_hypertreegrid_data->offsets != NULL) {
    free(vtk_hypertreegrid_data->offsets);
    vtk_hypertreegrid_data->offsets = NULL;
  }

#if _MPI
  int n_ranges;
  vtk_hypertreegrid_data->ranges = build_range_global(&n_ranges);
  vtk_hypertreegrid_data->offsets =
      build_offsets(vtk_hypertreegrid_data->ranges, &n_ranges);
  vtk_hypertreegrid_data->n_ranges = n_ranges;
#else
  vtk_hypertreegrid_data->ranges = NULL;
  vtk_hypertreegrid_data->offsets = NULL;
  vtk_hypertreegrid_data->n_ranges = 0;
#endif
}

void get_descriptors(vtkHyperTreeGridData *vtk_hypertreegrid_data) {

  size_t total_vertices = vtk_hypertreegrid_data->total_vertices;
  size_t max_depth = vtk_hypertreegrid_data->max_depth;

  size_t n_descriptors = total_vertices - vtk_hypertreegrid_data->number_of_vertices_per_depth[max_depth - 1];

  vtk_hypertreegrid_data->n_descriptors = n_descriptors;

  size_t descriptors_size = (size_t)(((n_descriptors + 7) / 8));

  vtk_hypertreegrid_data->descriptors_size = descriptors_size;

  if (pid() == 0) {
    vtk_hypertreegrid_data->descriptors =
        malloc(descriptors_size * sizeof(Bit_t));
    memset(vtk_hypertreegrid_data->descriptors, 0, descriptors_size);
  } else {
    vtk_hypertreegrid_data->descriptors =
        malloc(descriptors_size * sizeof(Bit_t));
    memset(vtk_hypertreegrid_data->descriptors, 0, descriptors_size);
  }

#if _MPI
  for (int ri = 0; ri < vtk_hypertreegrid_data->n_ranges; ri++) {
    range_t r = vtk_hypertreegrid_data->ranges[ri];
    size_t offset = vtk_hypertreegrid_data->offsets[ri];
    if (pid() == r.pid) {
      for (uint64_t c = r.begin; c <= r.end; c++) {
        Point point = morton_to_point(c);
        {
          size_t out_byte = offset / 8;
          size_t out_bit = offset % 8;
          if (!is_leaf(cell)) {
            vtk_hypertreegrid_data->descriptors[out_byte] |=
                (Bit_t)(1 << (7 - out_bit));
          }
          offset++;
        }
      }
    }
  }

  if (pid() == 0) {
    MPI_Reduce(MPI_IN_PLACE, vtk_hypertreegrid_data->descriptors,
               vtk_hypertreegrid_data->descriptors_size, MPI_UNSIGNED_CHAR,
               MPI_BOR, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(vtk_hypertreegrid_data->descriptors, NULL,
               vtk_hypertreegrid_data->descriptors_size, MPI_UNSIGNED_CHAR,
               MPI_BOR, 0, MPI_COMM_WORLD);
  }

#else
  uint8_t bit_count = 0;
  size_t byte_count = 0;

  foreach_cell_BFS() {
    if (!is_leaf(cell))
      vtk_hypertreegrid_data->descriptors[byte_count] |=
          (uint8_t)(1 << (7 - bit_count));

    bit_count++;

    if (bit_count == 8) {
      bit_count = 0;
      byte_count++;
    }
  }
#endif
}

void vtk_hypertreegrid_data_free(vtkHyperTreeGridData *vtk_hypertreegrid_data) {
  if (vtk_hypertreegrid_data) {
    free(vtk_hypertreegrid_data->ranges);
    free(vtk_hypertreegrid_data->offsets);
    free(vtk_hypertreegrid_data->x);
    free(vtk_hypertreegrid_data->y);
    free(vtk_hypertreegrid_data->z);
    free(vtk_hypertreegrid_data->number_of_vertices_per_depth);
    free(vtk_hypertreegrid_data->descriptors);
  }
}

vtkHyperTreeGridData *vtk_hypertreegrid_data_init(void) {
  vtkHyperTreeGridData *vtk_hypertreegrid_data =
      calloc(1, sizeof(vtkHyperTreeGridData));

  if (!vtk_hypertreegrid_data) {
    // handle malloc error
  };

  if (!grid) {
    fprintf(stderr, "Grid has not yet been initialized: Please initialize the "
                    "grid first.\n");
    // abort();
    exit(1);
  }

  get_coordinates(vtk_hypertreegrid_data);
  get_number_of_vertices_per_depth(vtk_hypertreegrid_data);
  get_ranges_and_offsets(vtk_hypertreegrid_data);
  get_descriptors(vtk_hypertreegrid_data);

  return vtk_hypertreegrid_data;
}

void vtk_hypertreegrid_data_update(
    vtkHyperTreeGridData *vtk_hypertreegrid_data) {
  get_number_of_vertices_per_depth(vtk_hypertreegrid_data);
  get_ranges_and_offsets(vtk_hypertreegrid_data);
  get_descriptors(vtk_hypertreegrid_data);
}