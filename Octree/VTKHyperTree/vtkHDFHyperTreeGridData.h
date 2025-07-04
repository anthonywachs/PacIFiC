#include "foreach_cell_bfs.h"
#include "vtkType.h"

typedef struct {
  size_t max_vertices;

  int64_t depth_per_tree;

  Bit_t *descriptors;
  size_t descriptors_size;
  int64_t n_descriptors;

  int64_t number_of_cells;
  int64_t *number_of_cells_per_tree_depth;
  int64_t number_of_depths;
  int64_t number_of_trees;
  int64_t tree_ids;

  bool has_mask;
  Bit_t *mask;
  size_t mask_size;
  size_t n_masks;

  double *x, *y, *z;
  size_t n_x, n_y, n_z;

} vtkHDFHyperTreeGridData;

void hdf_get_coordinates(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

#if dimension >= 1
  vtk_hdf_hypertreegrid_data->n_x = 2;
  vtk_hdf_hypertreegrid_data->x = malloc(2 * sizeof(double));
  vtk_hdf_hypertreegrid_data->x[0] = X0;
  vtk_hdf_hypertreegrid_data->x[1] = X0 + L0;
#else
  vtk_hdf_hypertreegrid_data->n_x = 1;
  vtk_hdf_hypertreegrid_data->x = malloc(1 * sizeof(double));
  vtk_hdf_hypertreegrid_data->x[0] = 0;
#endif

#if dimension >= 2
  vtk_hdf_hypertreegrid_data->n_y = 2;
  vtk_hdf_hypertreegrid_data->y = malloc(2 * sizeof(double));
  vtk_hdf_hypertreegrid_data->y[0] = Y0;
  vtk_hdf_hypertreegrid_data->y[1] = Y0 + L0;
#else
  vtk_hdf_hypertreegrid_data->n_y = 1;
  vtk_hdf_hypertreegrid_data->y = malloc(1 * sizeof(double));
  vtk_hdf_hypertreegrid_data->y[0] = 0;
#endif

#if dimension >= 3
  vtk_hdf_hypertreegrid_data->n_z = 2;
  vtk_hdf_hypertreegrid_data->z = malloc(2 * sizeof(double));
  vtk_hdf_hypertreegrid_data->z[0] = Z0;
  vtk_hdf_hypertreegrid_data->z[1] = Z0 + L0;
#else
  vtk_hdf_hypertreegrid_data->n_z = 1;
  vtk_hdf_hypertreegrid_data->z = malloc(1 * sizeof(double));
  vtk_hdf_hypertreegrid_data->z[0] = 0;
#endif

  // return vtk_hdf_hypertreegrid_data;
}

void hdf_get_number_of_vertices_per_depth(
    vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

  int depth;
  depth = grid->maxdepth;

  int64_t depth_per_tree_old = vtk_hdf_hypertreegrid_data->depth_per_tree;
  vtk_hdf_hypertreegrid_data->depth_per_tree = (int64_t)(depth + 1);

  if (vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth != NULL) {
    if (vtk_hdf_hypertreegrid_data->depth_per_tree != depth_per_tree_old) {
      free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth);
      vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth =
          malloc(vtk_hdf_hypertreegrid_data->depth_per_tree * sizeof(int64_t));
    }
  } else {
    vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth =
        malloc(vtk_hdf_hypertreegrid_data->depth_per_tree * sizeof(int64_t));
  }

  memset(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth, 0,
         vtk_hdf_hypertreegrid_data->depth_per_tree * sizeof(int64_t));

  // count your local cells per levelmax_depth
  foreach_cell_BFS() {
    vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[level]++;
  }

  // compute the max and total
  vtk_hdf_hypertreegrid_data->max_vertices = 0;
  vtk_hdf_hypertreegrid_data->number_of_cells = 0;
  for (size_t i = 0; i < vtk_hdf_hypertreegrid_data->depth_per_tree; i++) {
    int64_t nv = vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[i];
    if (nv > vtk_hdf_hypertreegrid_data->max_vertices)
      vtk_hdf_hypertreegrid_data->max_vertices = nv;
    vtk_hdf_hypertreegrid_data->number_of_cells += nv;
  }
}

void hdf_get_local_mask(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

  size_t number_of_cells = vtk_hdf_hypertreegrid_data->number_of_cells;
  // size_t max_depth = vtk_hdf_hypertreegrid_data->depth_per_tree;

  size_t n_masks = number_of_cells;
  vtk_hdf_hypertreegrid_data->n_masks = n_masks;
  size_t mask_size = (size_t)(((n_masks + 7) / 8));
  vtk_hdf_hypertreegrid_data->mask_size = mask_size;

  vtk_hdf_hypertreegrid_data->mask = malloc(mask_size * sizeof(Bit_t));
  memset(vtk_hdf_hypertreegrid_data->mask, 0, mask_size);

  uint8_t bit_count = 0;
  size_t byte_count = 0;

  foreach_cell_BFS() {
    if (!(is_local(cell)) && is_leaf(cell)) {
      bool val = true;
      foreach_neighbor(1) {
        if (is_local(cell)) {
          val = false;
        }
      }
      vtk_hdf_hypertreegrid_data->mask[byte_count] |=
          (uint8_t)(val << (7 - bit_count));
    }

    // if (!(is_local(cell)) && is_leaf(cell)) {
    //   vtk_hdf_hypertreegrid_data->mask[byte_count] |=
    //       (uint8_t)(1 << (7 - bit_count));
    // }

    bit_count++;

    if (bit_count == 8) {
      bit_count = 0;
      byte_count++;
    }
  }

  vtk_hdf_hypertreegrid_data->has_mask = true;
}

void hdf_get_descriptors(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

  size_t total_vertices = vtk_hdf_hypertreegrid_data->number_of_cells;
  size_t max_depth = vtk_hdf_hypertreegrid_data->depth_per_tree;

  size_t n_descriptors =
      total_vertices -
      vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[max_depth - 1];

  vtk_hdf_hypertreegrid_data->n_descriptors = (int64_t)n_descriptors;

  size_t descriptors_size = (size_t)(((n_descriptors + 7) / 8));

  vtk_hdf_hypertreegrid_data->descriptors_size = descriptors_size;

  vtk_hdf_hypertreegrid_data->descriptors =
      malloc(descriptors_size * sizeof(Bit_t));
  memset(vtk_hdf_hypertreegrid_data->descriptors, 0, descriptors_size);

  uint8_t bit_count = 0;
  size_t byte_count = 0;

  foreach_cell_BFS() {
    if (!is_leaf(cell))
      vtk_hdf_hypertreegrid_data->descriptors[byte_count] |=
          (uint8_t)(1 << (7 - bit_count));

    bit_count++;

    if (bit_count == 8) {
      bit_count = 0;
      byte_count++;
    }
  }
}

void vtk_hdf_hypertreegrid_data_free(
    vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {
  if (vtk_hdf_hypertreegrid_data) {
    free(vtk_hdf_hypertreegrid_data->x);
    free(vtk_hdf_hypertreegrid_data->y);
    free(vtk_hdf_hypertreegrid_data->z);
    free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth);
    free(vtk_hdf_hypertreegrid_data->descriptors);
    free(vtk_hdf_hypertreegrid_data->mask);
  }
}

vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data_init(void) {
  vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data =
      calloc(1, sizeof(vtkHDFHyperTreeGridData));

  if (!vtk_hdf_hypertreegrid_data) {
    // handle malloc error
  };

  if (!grid) {
    fprintf(stderr, "Grid has not yet been initialized: Please initialize the "
                    "grid first.\n");
    // abort();
    exit(1);
  }

  vtk_hdf_hypertreegrid_data->number_of_trees = 1;
  vtk_hdf_hypertreegrid_data->tree_ids = 0;

  hdf_get_coordinates(vtk_hdf_hypertreegrid_data);
  hdf_get_number_of_vertices_per_depth(vtk_hdf_hypertreegrid_data);
  hdf_get_descriptors(vtk_hdf_hypertreegrid_data);
#if MPI_SINGLE_FILE
  hdf_get_local_mask(vtk_hdf_hypertreegrid_data);
#endif

  return vtk_hdf_hypertreegrid_data;
}

void vtk_hdf_hypertreegrid_data_update(
    vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {
  hdf_get_number_of_vertices_per_depth(vtk_hdf_hypertreegrid_data);
  hdf_get_descriptors(vtk_hdf_hypertreegrid_data);
}