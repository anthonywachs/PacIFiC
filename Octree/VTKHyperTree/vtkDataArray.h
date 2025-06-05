#include "vtkFormat.h"
#include "vtkType.h"

// Writer
#if _MPI
typedef void (*vtkDataArrayWriter)(void *ctx, vtkFormat format, MPI_Offset base,
                                   MPI_File fp);
#else
typedef void (*vtkDataArrayWriter)(void *ctx, vtkFormat format, FILE *fp);
#endif

// Data array

typedef struct {
  vtkType type;
  char *name;
  size_t number_of_tuples;
  size_t number_of_components;
  void *writer_ctx;
  vtkDataArrayWriter writer;
} vtkDataArray;

void vtk_data_array_free(vtkDataArray *vtk_data_array) {
  free(vtk_data_array->name);
}

// Data array list
#define VTK_DATA_ARRAYS_CHUNK_SIZE 8

typedef struct {
  vtkDataArray **vtk_data_arrays;
  size_t count;
  size_t capacity;
  bool *list_owns;
} vtkDataArrayList;

vtkDataArrayList vtk_data_array_list_init() {
  return (vtkDataArrayList){
      .vtk_data_arrays = NULL, .count = 0, .capacity = 0, .list_owns = NULL};
}

void vtk_data_array_list_add(vtkDataArrayList *list,
                             vtkDataArray *vtk_data_array, bool list_owns) {
  if (list->count == list->capacity) {
    size_t new_capacity = list->capacity + VTK_DATA_ARRAYS_CHUNK_SIZE;

    // grow pointer array
    vtkDataArray **tmp_arrays =
        realloc(list->vtk_data_arrays, new_capacity * sizeof *tmp_arrays);
    if (!tmp_arrays) {
      perror("realloc vtk_data_arrays");
      exit(EXIT_FAILURE);
    }
    list->vtk_data_arrays = tmp_arrays;

    // grow ownership flags in lock-step
    bool *tmp_owns = realloc(list->list_owns, new_capacity * sizeof *tmp_owns);
    if (!tmp_owns) {
      perror("realloc list_owns");
      exit(EXIT_FAILURE);
    }
    list->list_owns = tmp_owns;

    list->capacity = new_capacity;
  }

  // take ownership by default
  list->vtk_data_arrays[list->count] = vtk_data_array;
  list->list_owns[list->count] = list_owns;
  list->count++;
}

void vtk_data_array_list_clear(vtkDataArrayList *list) {
  // free each owned entry
  for (size_t i = 0; i < list->count; i++) {
    if (list->list_owns[i]) {
      // free inner name
      vtk_data_array_free(list->vtk_data_arrays[i]);
      // free the struct itself
      free(list->vtk_data_arrays[i]);
    }
    // else: don't touch non-owned pointers
  }

  // tear down arrays
  free(list->vtk_data_arrays);
  free(list->list_owns);

  list->vtk_data_arrays = NULL;
  list->list_owns = NULL;
  list->capacity = 0;
  list->count = 0;
}

void vtk_data_array_list_free(vtkDataArrayList *list) {
  vtk_data_array_list_clear(list);
  // (if you ever heap-alloc vtkDataArrayList itself, free it here)
}