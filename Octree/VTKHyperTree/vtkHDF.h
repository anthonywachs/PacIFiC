#ifndef MPI_SINGLE_FILE
  #if _MPI
    #define MPI_SINGLE_FILE 1
  #else
    #define MPI_SINGLE_FILE 0
  #endif
#else
  #if _MPI == 0
    #define MPI_SINGLE_FILE 0
  #endif
#endif

#include <hdf5.h>

typedef struct {
  hid_t file_id;
  hid_t fapl;
  hid_t grp_vtkhdf_id;
} vtkHDF;

void vtk_HDF_close(vtkHDF *vtk_hdf) {
  if (vtk_hdf->grp_vtkhdf_id >= 0)
    H5Gclose(vtk_hdf->grp_vtkhdf_id);
  if (vtk_hdf->fapl >= 0)
    H5Pclose(vtk_hdf->fapl);
  if (vtk_hdf->file_id >= 0)
    H5Fclose(vtk_hdf->file_id);
}

void vtk_HDF_error(vtkHDF *vtk_hdf) { vtk_HDF_close(vtk_hdf); }

vtkHDF vtk_HDF_init(const char *fname) {
  /* File and group identifiers */
  vtkHDF vtk_hdf = {.file_id = -1, .fapl = -1, .grp_vtkhdf_id = -1};

#if _MPI
  vtk_hdf.fapl = H5Pcreate(H5P_FILE_ACCESS);
  if (vtk_hdf.fapl < 0)
    vtk_HDF_error(&vtk_hdf);

#if MPI_SINGLE_FILE
  if (H5Pset_fapl_mpio(vtk_hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
    vtk_HDF_error(&vtk_hdf);
#else
  if (H5Pset_fapl_mpio(vtk_hdf.fapl, MPI_COMM_SELF, MPI_INFO_NULL) < 0)
    vtk_HDF_error(&vtk_hdf);
#endif

  vtk_hdf.file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, vtk_hdf.fapl);
  if (vtk_hdf.file_id < 0)
    vtk_HDF_error(&vtk_hdf);

  vtk_hdf.grp_vtkhdf_id = H5Gcreate2(vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.grp_vtkhdf_id < 0)
    vtk_HDF_error(&vtk_hdf);

#else
  vtk_hdf.file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.file_id < 0)
    vtk_HDF_error(&vtk_hdf);

  vtk_hdf.grp_vtkhdf_id = H5Gcreate2(vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.grp_vtkhdf_id < 0)
    vtk_HDF_error(&vtk_hdf);
#endif

  return vtk_hdf;
}
