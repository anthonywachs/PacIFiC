#include "foreach_cell_bfs.h"
#include "morton.h"

typedef struct {
  int begin, end;
  int pid;
} range_t;

static int cmp_range(const void *a, const void *b) {
  const range_t *ra = a;
  const range_t *rb = b;
  if (ra->begin < rb->begin)
    return -1;
  if (ra->begin > rb->begin)
    return 1;
  return 0;
}

void qsort_range(range_t *range, int n) {
  qsort(range, n, sizeof *range, cmp_range);
}

void print_range(range_t *range, int n) {
  for (int i = 0; i < n; i++) {
    range_t r = range[i];
    printf("[Proc %i]: (%i,%i)\n", r.pid, r.begin, r.end);
  }
}

void print_range_cells(range_t *range, int n) {
  for (int i = 0; i < n; i++) {
    range_t r = range[i];
    if (r.pid == pid()) {
      for (uint64_t c = r.begin; c <= r.end; c++) {
        Point p = morton_to_point(c);
        printf("[Proc %i]: i = %d, j = %d, level = %d\n", pid(), p.i, p.j,
               p.level);
      }
    }
  }
}

range_t *build_range_local(int *size) {
  // 1) count how many local cells there are
  int cell_count = 0;
  foreach_cell() {
    if (is_local(cell))
      cell_count++;
  }

  // 2) allocate enough space for the worst case
  range_t *range_local = malloc(cell_count * sizeof *range_local);
  int rc = 0;     // how many ranges we've actually filled
  int in_run = 0; // are we in the middle of a run?
  uint64_t run_begin = 0, prev = 0;
  int me = pid();

  // 3) scan in BFS order
  foreach_cell_BFS() {
    if (is_local(cell)) {
      uint64_t code = point_to_morton(point);

      if (!in_run) {
        // start the first run
        run_begin = prev = code;
        in_run = 1;
      } else if (code == prev + 1) {
        // extend current run
        prev = code;
      } else {
        // gap!  close out the old run ...
        range_local[rc++] =
            (range_t){.begin = (int)run_begin, .end = (int)prev, .pid = me};
        // … and immediately start a new one
        run_begin = prev = code;
      }
    }
  }

  // 4) flush the final run if any
  if (in_run) {
    range_local[rc++] =
        (range_t){.begin = (int)run_begin, .end = (int)prev, .pid = me};
  }

  // 5) shrink down to exactly rc entries
  range_t *tmp = realloc(range_local, rc * sizeof *range_local);
  if (tmp)
    range_local = tmp;

  // 6) hand back the count & pointer
  *size = rc;
  return range_local;
}

range_t *build_range_global(int *size) {
  // Get the local range
  int rc_local = 0;
  range_t *range_local = build_range_local(&rc_local);

#if _MPI
  // MPI Information
  int npe = npe();

  // Gather the value of rc_local on every proc
  int *rc_all = malloc((int)npe * sizeof *rc_all);
  MPI_Allgather(&rc_local, 1, MPI_INT, rc_all, 1, MPI_INT, MPI_COMM_WORLD);

  // Count the total number of ranges globally
  int rc_global = 0;
  for (int i = 0; i < npe; i++) {
    rc_global += rc_all[i];
  }

  MPI_Datatype MPI_RANGE_T;
  {
    range_t dummy;
    MPI_Aint base, addr_begin, addr_end, addr_pid;
    MPI_Aint offsets[3];
    int blocklens[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};

    /* 1) get absolute addresses of each field */
    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.begin, &addr_begin);
    MPI_Get_address(&dummy.end, &addr_end);
    MPI_Get_address(&dummy.pid, &addr_pid);

    /* 2) convert to relative offsets */
    offsets[0] = addr_begin - base;
    offsets[1] = addr_end - base;
    offsets[2] = addr_pid - base;

    /* 3) build and commit the struct type */
    MPI_Type_create_struct(3, blocklens, offsets, types, &MPI_RANGE_T);
    MPI_Type_commit(&MPI_RANGE_T);
  }

  // Compute the displacements into the global range array
  int *displs = malloc(npe * sizeof *displs);
  displs[0] = 0;

  for (int r = 1; r < npe; r++)
    displs[r] = displs[r - 1] + rc_all[r - 1];

  int total = displs[npe - 1] + rc_all[npe - 1];

  range_t *range_global = malloc(total * sizeof *range_global);

  MPI_Allgatherv(range_local,  // send buffer
                 rc_local,     // send count
                 MPI_RANGE_T,  // send type
                 range_global, // recv buffer
                 rc_all,       // recv counts
                 displs,       // displacements
                 MPI_RANGE_T,  // recv type
                 MPI_COMM_WORLD);

  free(rc_all);
  free(displs);
  MPI_Type_free(&MPI_RANGE_T);

  // Sort the global ranges
  qsort_range(range_global, rc_global);

  *size = rc_global;
  return range_global;
#else
  *size = rc_local;
  return range_local;
#endif
}

size_t *build_offsets(range_t *ranges, int *size) {
  size_t *sizes = malloc((*size) * sizeof(size_t));

  for (int i = 0; i < *size; i++) {
    sizes[i] = (size_t)(ranges[i].end - ranges[i].begin + 1);
  }

  size_t *offsets = malloc((*size) * sizeof(size_t));
  size_t running = 0;

  for (int i = 0; i < *size; i++) {
    offsets[i] = running;
    running += sizes[i];
  }

  free(sizes);
  return offsets;
}