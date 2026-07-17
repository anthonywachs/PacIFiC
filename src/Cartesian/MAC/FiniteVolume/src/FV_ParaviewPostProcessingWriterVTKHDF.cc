// #include <H5Ipublic.h>
// #include <H5Tpublic.h>
// #include <H5public.h>
#include <H5Gpublic.h>
#include <H5Tpublic.h>
#include <H5public.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include <cmath>
#include <cstdlib>
#include <cstring>

namespace {

/**
 * @brief vtkHDF file type
 *
 * @memberof vtkHDF
 */
typedef enum {
  IMAGEDATA,
  POLYDATA,
  UNSTRUCTUREDGRID,
  HYPERTREEGRID,
  OVERLAPPINGAMR,
  PARTITIONEDDATASETCOLLECTION,
  MULTIBLOCKDATASET
} vtkHDFType;

/**
 * @brief Major and minor version of the vtkHDF file format
 *
 * @memberof vtkHDF
 */
typedef struct {
  int major;
  int minor;
} vtkHDFVersion;

/**
 * @brief Parent class with many helper functions for writing VTKHDF files
 */
typedef struct {
  hid_t file_id;       /**< File id */
  hid_t fapl;          /**< File access property list  */
  hid_t grp_vtkhdf_id; /**< VTKHDF group id */

  hid_t attr_id;       /**< Attribute id */
  hid_t attr_space_id; /**< Attribute dataspace id */
  hid_t attr_dtype_id; /**< Attribute datatype id */

  hid_t dset_id;       /**< Dataset id */
  hid_t dset_space_id; /**< Dataset dataspace id */
  hid_t dset_dtype_id; /**< Dataset datatype id */

  hid_t dcpl_id; /**< Dataset property-list id */

  hid_t file_space; /**< MPI-IO file space */
  hid_t mem_space;  /**< MPI-IO memory space */
  hid_t xfer_plist; /**< MPI-IO  */
} vtkHDF;

/**
 * @brief Close a vtkHDF file
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_close(vtkHDF *vtk_hdf) {
  if (vtk_hdf->grp_vtkhdf_id >= 0 && H5Iis_valid(vtk_hdf->grp_vtkhdf_id) > 0) {
    H5Gclose(vtk_hdf->grp_vtkhdf_id);
    vtk_hdf->grp_vtkhdf_id = H5I_INVALID_HID;
  }

  if (vtk_hdf->xfer_plist >= 0 && H5Iis_valid(vtk_hdf->xfer_plist) > 0) {
    H5Pclose(vtk_hdf->xfer_plist);
    vtk_hdf->xfer_plist = H5I_INVALID_HID;
  }

  if (vtk_hdf->mem_space >= 0 && H5Iis_valid(vtk_hdf->mem_space) > 0) {
    H5Sclose(vtk_hdf->mem_space);
    vtk_hdf->mem_space = H5I_INVALID_HID;
  }

  if (vtk_hdf->file_space >= 0 && H5Iis_valid(vtk_hdf->file_space) > 0) {
    H5Sclose(vtk_hdf->file_space);
    vtk_hdf->file_space = H5I_INVALID_HID;
  }

  if (vtk_hdf->fapl >= 0 && H5Iis_valid(vtk_hdf->fapl) > 0) {
    H5Pclose(vtk_hdf->fapl);
    vtk_hdf->fapl = H5I_INVALID_HID;
  }

  if (vtk_hdf->file_id >= 0 && H5Iis_valid(vtk_hdf->file_id) > 0) {
    H5Fclose(vtk_hdf->file_id);
    vtk_hdf->file_id = H5I_INVALID_HID;
  }
}

/**
 * @brief Function called on any error which immediately closes the file to
 * avoid corruptions
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_error(vtkHDF *vtk_hdf) { vtk_HDF_close(vtk_hdf); }

/**
 * @brief Check for error of some operation that returns herr_t type
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 * @param result Result of a HDF5 operation to check
 *
 * @memberof vtkHDF
 */
void vtk_HDF_check_result(vtkHDF *vtk_hdf, herr_t result) {
  if (result < 0) {
    vtk_HDF_error(vtk_hdf);
  }
}

/**
 * @brief Check for error of some operation that returns an object id
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 * @param object_id Object id to check
 *
 * @memberof vtkHDF
 */
void vtk_HDF_check_object(vtkHDF *vtk_hdf, hid_t object_id) {
  if (object_id <= H5I_INVALID_HID) {
    vtk_HDF_error(vtk_hdf);
  }
}

/**
 * @brief Initialize the vtkHDF struct by writing a new file
 *
 * @memberof vtkHDF
 */
vtkHDF vtk_HDF_init(const char *fname, bool overwrite) {
  vtkHDF vtk_hdf = {
      .file_id = H5I_INVALID_HID,       /**< File id */
      .fapl = H5I_INVALID_HID,          /**< File access property list  */
      .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
      .attr_id = H5I_INVALID_HID,       /**< Attribute id */
      .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
      .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
      .dset_id = H5I_INVALID_HID,       /**< Dataset id */
      .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
      .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
      .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
      .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
      .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
      .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  // Create the actual file on disk
  vtk_hdf.file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf, vtk_hdf.file_id);

  // Create the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gcreate2(vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf, vtk_hdf.grp_vtkhdf_id);

  return vtk_hdf;
}

vtkHDF vtk_HDF_init_MPIIO(const char *fname, bool overwrite) {
  vtkHDF vtk_hdf = {
      .file_id = H5I_INVALID_HID,       /**< File id */
      .fapl = H5I_INVALID_HID,          /**< File access property list  */
      .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
      .attr_id = H5I_INVALID_HID,       /**< Attribute id */
      .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
      .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
      .dset_id = H5I_INVALID_HID,       /**< Dataset id */
      .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
      .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
      .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
      .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
      .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
      .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  // Create a file access property list
  vtk_hdf.fapl = H5Pcreate(H5P_FILE_ACCESS);
  if (vtk_hdf.fapl < 0)
    vtk_HDF_error(&vtk_hdf);

  // If we use MPI_SINGLE_FILE, we set all procs to have access and write into
  // the same file
  if (H5Pset_fapl_mpio(vtk_hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
    vtk_HDF_error(&vtk_hdf);

  // Create the acutal file on disk
  vtk_hdf.file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, vtk_hdf.fapl);
  if (vtk_hdf.file_id < 0)
    vtk_HDF_error(&vtk_hdf);

  // Create the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gcreate2(vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.grp_vtkhdf_id < 0)
    vtk_HDF_error(&vtk_hdf);

  return vtk_hdf;
}

/**
 * @brief Initialize the vtkHDF struct by opening an existing file
 *
 * @memberof vtkHDF
 */
vtkHDF vtk_HDF_open(const char *fname) {
  vtkHDF vtk_hdf = {
      .file_id = H5I_INVALID_HID,       /**< File id */
      .fapl = H5I_INVALID_HID,          /**< File access property list  */
      .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
      .attr_id = H5I_INVALID_HID,       /**< Attribute id */
      .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
      .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
      .dset_id = H5I_INVALID_HID,       /**< Dataset id */
      .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
      .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
      .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
      .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
      .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
      .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  // Open the actual file on disk
  vtk_hdf.file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf, vtk_hdf.file_id);

  // Open the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gopen(vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf, vtk_hdf.grp_vtkhdf_id);

  return vtk_hdf;
}

/**
 * @brief Initialize the vtkHDF struct by opening an existing file
 *
 * @memberof vtkHDF
 */
vtkHDF vtk_HDF_open_MPIIO(const char *fname) {
  vtkHDF vtk_hdf = {
      .file_id = H5I_INVALID_HID,       /**< File id */
      .fapl = H5I_INVALID_HID,          /**< File access property list  */
      .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
      .attr_id = H5I_INVALID_HID,       /**< Attribute id */
      .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
      .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
      .dset_id = H5I_INVALID_HID,       /**< Dataset id */
      .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
      .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
      .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
      .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
      .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
      .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  herr_t result = 0;
  // Create a file access property list
  vtk_hdf.fapl = H5Pcreate(H5P_FILE_ACCESS);
  vtk_HDF_check_object(&vtk_hdf, vtk_hdf.fapl);

  // If we use MPI_SINGLE_FILE, we set all procs to have access and write into
  // the same file
  result = H5Pset_fapl_mpio(vtk_hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  vtk_HDF_check_result(&vtk_hdf, result);

  // Open the actual file on disk
  vtk_hdf.file_id = H5Fopen(fname, H5F_ACC_RDWR, vtk_hdf.fapl);
  vtk_HDF_check_object(&vtk_hdf, vtk_hdf.file_id);

  // Open the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gopen(vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf, vtk_hdf.grp_vtkhdf_id);

  return vtk_hdf;
}

/**
 * @brief Write new attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param dims The dimensions of the data
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_attribute(const char *attribute_name, const void *data,
                             const hid_t dtype_id, const hid_t group_id,
                             const hsize_t dims[], vtkHDF *vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  vtk_hdf->attr_space_id = H5Screate_simple(1, dims, dims);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_space_id);

  // Set the datatype
  vtk_hdf->attr_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_dtype_id);

  // Create the attribute
  vtk_hdf->attr_id =
      H5Acreate2(group_id, attribute_name, vtk_hdf->attr_dtype_id,
                 vtk_hdf->attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_id);

  // Write the attribute
  result = H5Awrite(vtk_hdf->attr_id, vtk_hdf->attr_dtype_id, data);
  vtk_HDF_check_object(vtk_hdf, result);

  // Close the open objects
  H5Aclose(vtk_hdf->attr_id);
  H5Tclose(vtk_hdf->attr_dtype_id);
  H5Sclose(vtk_hdf->attr_space_id);
}

/**
 * @brief Write new scalar attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */

void vtk_HDF_write_scalar_attribute(const char *attribute_name,
                                    const void *data, const hid_t dtype_id,
                                    const hid_t group_id, vtkHDF *vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  vtk_hdf->attr_space_id = H5Screate(H5S_SCALAR);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_space_id);

  // Set the datatype
  vtk_hdf->attr_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_dtype_id);

  // Create the attribute
  vtk_hdf->attr_id =
      H5Acreate2(group_id, attribute_name, vtk_hdf->attr_dtype_id,
                 vtk_hdf->attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_id);

  // Write the attribute
  result = H5Awrite(vtk_hdf->attr_id, vtk_hdf->attr_dtype_id, data);
  vtk_HDF_check_object(vtk_hdf, result);

  // Close the open objects
  H5Aclose(vtk_hdf->attr_id);
  H5Tclose(vtk_hdf->attr_dtype_id);
  H5Sclose(vtk_hdf->attr_space_id);
}

/**
 * @brief Write new attribute
 *
 * @param type_name The name of the attrbiute
 * @param group_id The HDF5 group to write the attribute in
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_type_attribute(const char *type_name, const hid_t group_id,
                                  vtkHDF *vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  vtk_hdf->attr_space_id = H5Screate(H5S_SCALAR);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_space_id);

  // Set the datatype
  vtk_hdf->attr_dtype_id = H5Tcopy(H5T_C_S1);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_dtype_id);

  // Set the size
  result = H5Tset_size(vtk_hdf->attr_dtype_id, strlen(type_name));
  vtk_HDF_check_result(vtk_hdf, result);

  // Pad the string
  result = H5Tset_strpad(vtk_hdf->attr_dtype_id, H5T_STR_NULLTERM);
  vtk_HDF_check_result(vtk_hdf, result);

  // Set the string value
  result = H5Tset_cset(vtk_hdf->attr_dtype_id, H5T_CSET_ASCII);
  vtk_HDF_check_result(vtk_hdf, result);

  // Create the attribute
  vtk_hdf->attr_id =
      H5Acreate2(group_id, "Type", vtk_hdf->attr_dtype_id,
                 vtk_hdf->attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->attr_id);

  // Write the attribute
  result = H5Awrite(vtk_hdf->attr_id, vtk_hdf->attr_dtype_id, type_name);
  vtk_HDF_check_object(vtk_hdf, result);

  // Close the open objects
  H5Aclose(vtk_hdf->attr_id);
  H5Tclose(vtk_hdf->attr_dtype_id);
  H5Sclose(vtk_hdf->attr_space_id);
}

/**
 * @brief Write new unchunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_dataset(const char *dataset_name, const void *data,
                           const hid_t dtype_id, const hid_t group_id,
                           const int rank, const hsize_t dims[],
                           vtkHDF *vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple(rank, dims, dims);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dcpl_id);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2(group_id, dataset_name, vtk_hdf->dset_dtype_id,
                                vtk_hdf->dset_space_id, H5P_DEFAULT,
                                vtk_hdf->dcpl_id, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_id);

  result = H5Dwrite(vtk_hdf->dset_id, vtk_hdf->dset_dtype_id, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, data);
  vtk_HDF_check_result(vtk_hdf, result);

  // Close opened objects
  H5Dclose(vtk_hdf->dset_id);
  H5Tclose(vtk_hdf->dset_dtype_id);
  H5Pclose(vtk_hdf->dcpl_id);
  H5Sclose(vtk_hdf->dset_space_id);
}

/**
 * @brief Write new chunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 */
void vtk_HDF_write_chunked_dataset(const char *dataset_name, const void *data,
                                   const hid_t dtype_id, const hid_t group_id,
                                   const int rank, const hsize_t dims[],
                                   const hsize_t max_dims[],
                                   const hsize_t chunk_dims[],
                                   vtkHDF *vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result(vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2(group_id, dataset_name, vtk_hdf->dset_dtype_id,
                                vtk_hdf->dset_space_id, H5P_DEFAULT,
                                vtk_hdf->dcpl_id, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_id);

  result = H5Dwrite(vtk_hdf->dset_id, vtk_hdf->dset_dtype_id, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, data);
  vtk_HDF_check_result(vtk_hdf, result);

  // Close opened objects
  H5Dclose(vtk_hdf->dset_id);
  H5Tclose(vtk_hdf->dset_dtype_id);
  H5Pclose(vtk_hdf->dcpl_id);
  H5Sclose(vtk_hdf->dset_space_id);
}

/**
 * @brief Append existing chunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset to append
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 */

void vtk_HDF_append_chunked_dataset(const char *dataset_name, const void *data,
                                    const hid_t dtype_id, const hid_t group_id,
                                    const int rank, const hsize_t dims[],
                                    vtkHDF *vtk_hdf) {
  herr_t result = 0;

  // Open existing dataset
  vtk_hdf->dset_id = H5Dopen2(group_id, dataset_name, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_id);

  // Get current dataspace and its dimensions
  vtk_hdf->dset_space_id = H5Dget_space(vtk_hdf->dset_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  hsize_t cur_dims[H5S_MAX_RANK];
  hsize_t cur_max_dims[H5S_MAX_RANK];
  int file_rank =
      H5Sget_simple_extent_dims(vtk_hdf->dset_space_id, cur_dims, cur_max_dims);
  vtk_HDF_check_result(vtk_hdf, file_rank);

  // Check that the rank should match
  if (file_rank != rank) {
    perror("Rank mismatch when appending chunked dataset.\n");
    // exit(1);
  }

  // We assume we append along dimension 0, and that all other dimensions stay
  // the same.
  hsize_t new_dims[H5S_MAX_RANK];
  for (int i = 0; i < rank; ++i) {
    new_dims[i] = cur_dims[i];
  }
  new_dims[0] = cur_dims[0] + dims[0];

  // Ensure other dims are unchanged
  for (int i = 1; i < rank; ++i) {
    if (dims[i] != cur_dims[i]) {
      perror("Append dims must match existing dims in non-unlimited axes.\n");
      // exit(1);
    }
  }

  // Extend the dataset to the new size
  result = H5Dset_extent(vtk_hdf->dset_id, new_dims);
  vtk_HDF_check_result(vtk_hdf, result);

  // After extending, the old dataspace is invalid; close and re-open
  H5Sclose(vtk_hdf->dset_space_id);
  vtk_hdf->dset_space_id = H5Dget_space(vtk_hdf->dset_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  // 4. Select hyperslab corresponding to the newly appended region
  hsize_t start[H5S_MAX_RANK];
  hsize_t count[H5S_MAX_RANK];

  start[0] = cur_dims[0]; // start right after old end
  count[0] = dims[0];     // append this many along dim 0

  for (int i = 1; i < rank; ++i) {
    start[i] = 0;
    count[i] = dims[i]; // whole extent in other dims
  }

  result = H5Sselect_hyperslab(vtk_hdf->dset_space_id, H5S_SELECT_SET, start,
                               NULL, // stride
                               count,
                               NULL); // block

  vtk_HDF_check_result(vtk_hdf, result);

  // Create memspace for the data block we’re appending
  vtk_hdf->mem_space = H5Screate_simple(rank, dims, NULL);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->mem_space);

  // Write the new chunk into the extended portion
  result = H5Dwrite(vtk_hdf->dset_id,
                    dtype_id,               // in-memory type
                    vtk_hdf->mem_space,     // mem dataspace selection
                    vtk_hdf->dset_space_id, // file dataspace selection
                    H5P_DEFAULT,            // xfer plist
                    data);
  vtk_HDF_check_result(vtk_hdf, result);

  // Close everything
  H5Sclose(vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;

  H5Sclose(vtk_hdf->dset_space_id);
  vtk_hdf->dset_space_id = H5I_INVALID_HID;

  H5Dclose(vtk_hdf->dset_id);
  vtk_hdf->dset_id = H5I_INVALID_HID;
}

/**
 * @brief Write new chunked and compressed dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param compression_level The level of compression
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_compressed_dataset(
    const char *dataset_name, const void *data, const hid_t dtype_id,
    const hid_t group_id, const int rank, const hsize_t dims[],
    const hsize_t max_dims[], const hsize_t chunk_dims[],
    unsigned int compression_level, vtkHDF *vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result(vtk_hdf, result);

  // Set compression on the dataset
  result = H5Pset_deflate(vtk_hdf->dcpl_id, compression_level);
  vtk_HDF_check_result(vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2(group_id, dataset_name, vtk_hdf->dset_dtype_id,
                                vtk_hdf->dset_space_id, H5P_DEFAULT,
                                vtk_hdf->dcpl_id, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_id);

  // Write the dataset
  result = H5Dwrite(vtk_hdf->dset_id, vtk_hdf->dset_dtype_id, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, data);
  vtk_HDF_check_result(vtk_hdf, result);

  // Close opened objects
  H5Dclose(vtk_hdf->dset_id);
  H5Tclose(vtk_hdf->dset_dtype_id);
  H5Pclose(vtk_hdf->dcpl_id);
  H5Sclose(vtk_hdf->dset_space_id);
}

/**
 * @brief Write a collectively written dataset with MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_collective_write_dataset(
    const char *dataset_name, const void *data, const hid_t dtype_id,
    const hid_t group_id, const int rank, const hsize_t dims[],
    const hsize_t local_size[], const hsize_t local_offset[], vtkHDF *vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple(rank, dims, dims);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id =
      H5Dcreate2(group_id, dataset_name, vtk_hdf->dset_dtype_id,
                 vtk_hdf->dset_space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_id);

  // Create a property list for MPI-IO transfer
  vtk_hdf->xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio(vtk_hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  vtk_HDF_check_result(vtk_hdf, result);

  // Get the dataset space
  vtk_hdf->file_space = H5Dget_space(vtk_hdf->dset_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab(vtk_hdf->file_space, H5S_SELECT_SET,
                               local_offset, NULL, local_size, NULL);
  vtk_HDF_check_result(vtk_hdf, result);

  // Create a memory dataspace for our process
  vtk_hdf->mem_space = H5Screate_simple(rank, local_size, NULL);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->mem_space);

  // Actually do the collective write to the file
  result =
      H5Dwrite(vtk_hdf->dset_id,       /* dataset handle */
               vtk_hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
               vtk_hdf->mem_space,     /* memory dataspace [local_nx] */
               vtk_hdf->file_space, /* file dataspace with hyperslab selected */
               vtk_hdf->xfer_plist, /* collective MPI‐IO transfer property */
               data                 /* pointer to local data */
      );
  vtk_HDF_check_result(vtk_hdf, result);

  // Close opened objects
  H5Pclose(vtk_hdf->xfer_plist);
  vtk_hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose(vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;
  H5Sclose(vtk_hdf->file_space);
  vtk_hdf->file_space = H5I_INVALID_HID;
  H5Tclose(vtk_hdf->dset_dtype_id);
  H5Sclose(vtk_hdf->dset_space_id);
  H5Dclose(vtk_hdf->dset_id);
}

/**
 * @brief Write a chunked dataset collectively using MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param dims The dimensions of the dataset
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_collective_write_chunked_dataset(
    const char *dataset_name, const void *data, const hid_t dtype_id,
    const hid_t group_id, const int rank, const hsize_t dims[],
    const hsize_t max_dims[], const hsize_t chunk_dims[],
    const hsize_t local_size[], const hsize_t local_offset[], vtkHDF *vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result(vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2(group_id, dataset_name, vtk_hdf->dset_dtype_id,
                                vtk_hdf->dset_space_id, H5P_DEFAULT,
                                vtk_hdf->dcpl_id, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_id);

  // Create a property list for MPI-IO transfer
  vtk_hdf->xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio(vtk_hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  vtk_HDF_check_result(vtk_hdf, result);

  // Get the dataset space
  vtk_hdf->file_space = H5Dget_space(vtk_hdf->dset_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab(vtk_hdf->file_space, H5S_SELECT_SET,
                               local_offset, NULL, local_size, NULL);
  vtk_HDF_check_result(vtk_hdf, result);

  // Create a memory dataspace for our process
  vtk_hdf->mem_space = H5Screate_simple(rank, local_size, NULL);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->mem_space);

  // Actually do the collective write to the file
  result =
      H5Dwrite(vtk_hdf->dset_id,       /* dataset handle */
               vtk_hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
               vtk_hdf->mem_space,     /* memory dataspace */
               vtk_hdf->file_space, /* file dataspace with hyperslab selected */
               vtk_hdf->xfer_plist, /* collective MPI‐IO transfer property */
               data                 /* pointer to local data */
      );
  vtk_HDF_check_result(vtk_hdf, result);

  // Close opened objects
  H5Pclose(vtk_hdf->xfer_plist);
  vtk_hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose(vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;
  H5Sclose(vtk_hdf->file_space);
  vtk_hdf->file_space = H5I_INVALID_HID;
  H5Tclose(vtk_hdf->dset_dtype_id);
  H5Pclose(vtk_hdf->dcpl_id);
  H5Sclose(vtk_hdf->dset_space_id);
  H5Dclose(vtk_hdf->dset_id);
}

/**
 * @brief Write a collectively written dataset using compression with MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param compression_level The level of compression
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */

void vtk_HDF_collective_write_compressed_dataset(
    const char *dataset_name, const void *data, const hid_t dtype_id,
    const hid_t group_id, const int rank, const hsize_t dims[],
    const hsize_t max_dims[], const hsize_t chunk_dims[],
    const hsize_t local_size[], const hsize_t local_offset[],
    unsigned int compression_level, vtkHDF *vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result(vtk_hdf, result);

  // Set compression on the dataset
  result = H5Pset_deflate(vtk_hdf->dcpl_id, compression_level);
  vtk_HDF_check_result(vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy(dtype_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2(group_id, dataset_name, vtk_hdf->dset_dtype_id,
                                vtk_hdf->dset_space_id, H5P_DEFAULT,
                                vtk_hdf->dcpl_id, H5P_DEFAULT);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->dset_id);

  // Create a property list for MPI-IO transfer
  vtk_hdf->xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio(vtk_hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  vtk_HDF_check_result(vtk_hdf, result);

  // Get the dataset space
  vtk_hdf->file_space = H5Dget_space(vtk_hdf->dset_id);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab(vtk_hdf->file_space, H5S_SELECT_SET,
                               local_offset, NULL, local_size, NULL);
  vtk_HDF_check_result(vtk_hdf, result);

  // Create a memory dataspace for our process
  vtk_hdf->mem_space = H5Screate_simple(rank, local_size, NULL);
  vtk_HDF_check_object(vtk_hdf, vtk_hdf->mem_space);

  // Actually do the collective write to the file
  result =
      H5Dwrite(vtk_hdf->dset_id,       /* dataset handle */
               vtk_hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
               vtk_hdf->mem_space,     /* memory dataspace [local_nx] */
               vtk_hdf->file_space, /* file dataspace with hyperslab selected */
               vtk_hdf->xfer_plist, /* collective MPI‐IO transfer property */
               data                 /* pointer to local data */
      );
  vtk_HDF_check_result(vtk_hdf, result);

  // Close opened objects
  H5Pclose(vtk_hdf->xfer_plist);
  vtk_hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose(vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;
  H5Sclose(vtk_hdf->file_space);
  vtk_hdf->file_space = H5I_INVALID_HID;
  H5Tclose(vtk_hdf->dset_dtype_id);
  H5Pclose(vtk_hdf->dcpl_id);
  H5Sclose(vtk_hdf->dset_space_id);
  H5Dclose(vtk_hdf->dset_id);
}

}; // namespace

#include <FV_ParaviewPostProcessingWriterVTKHDF.hh>

#include <FV.hh>
#include <FV_DomainAndFields.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>

#include <MAC_Communicator.hh>
#include <MAC_Data.hh>
#include <MAC_KeywordDataIterator.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_assertions.hh>

#include <doubleArray2D.hh>

#include <fstream>
#include <sstream>
#include <vector>

using std::ostringstream;

namespace {

struct VTKHDFPackedField {
  std::string dataset_name;
  size_t ncomp = 0;
  size_t ntuples = 0;
  std::vector<float> values; // tuple-major, component-minor
};

bool pack_cell_centered_field_data(FV_DiscreteField const *f,
                                   VTKHDFPackedField &out) {
  if (f == 0 || f->paraview_location() != "at_cell_centers") {
    return false;
  }

  MAC_Module *point_data = MAC_Module::create(0, "PointData");
  MAC_Module *cell_data = MAC_Module::create(0, "CellData");
  f->write_field(point_data, cell_data);

  if (!cell_data->has_entry()) {
    point_data->destroy();
    cell_data->destroy();
    return false;
  }

  MAC_KeywordDataIterator *it = cell_data->create_entry_iterator(0);
  it->start();
  bool found_array = false;
  for (; it->is_valid(); it->go_next()) {
    MAC_KeywordDataPair *pair = it->item();
    MAC_Data const *data = pair->data();
    if (data->data_type() == MAC_Data::DoubleArray2D) {
      doubleArray2D const &arr = data->to_double_array2D();
      out.dataset_name = pair->keyword();
      out.ncomp = arr.index_bound(0);
      out.ntuples = arr.index_bound(1);
      out.values.assign(out.ncomp * out.ntuples, 0.0f);

      for (size_t t = 0; t < out.ntuples; ++t) {
        for (size_t c = 0; c < out.ncomp; ++c) {
          out.values[t * out.ncomp + c] = static_cast<float>(arr(c, t));
        }
      }
      found_array = true;
      break;
    }
  }

  it->destroy();
  point_data->destroy();
  cell_data->destroy();
  return found_array;
}

} // namespace

FV_ParaviewPostProcessingWriterVTKHDF const
    *FV_ParaviewPostProcessingWriterVTKHDF::PROTOTYPE =
        new FV_ParaviewPostProcessingWriterVTKHDF("paraview_vtkhdf");

//----------------------------------------------------------------------
FV_ParaviewPostProcessingWriterVTKHDF::FV_ParaviewPostProcessingWriterVTKHDF(
    std::string const &a_name)
    : FV_PostProcessingWriter(a_name), EXP(0), COM(0), FVFIELDS(),
      PRIMARY_GRID(0), RES_DIRECTORY("Res"), BASE_FILENAME("save"),
      VTKHDF_SERIES_FILENAME("Res/save.vtkhdf.series"),
      VTKHDF_SERIES_IDX_FILENAME("Res/save.vtkhdf.series.idx"),
      TIME_SERIES_RECORDS(), TIME_SERIES_LOADED(false),
      TIME_SERIES_RESTART_MODE(false), TIME_SERIES_PREVIOUS_CYCLE(0),
      CYCLE_NUMBER(0), BINARY(true)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: "
            "FV_ParaviewPostProcessingWriterVTKHDF");
}

//----------------------------------------------------------------------
FV_PostProcessingWriter *FV_ParaviewPostProcessingWriterVTKHDF::create_replica(
    MAC_Object *a_owner, MAC_ModuleExplorer const *exp,
    MAC_Communicator const *com, list<FV_DiscreteField const *> a_fields,
    FV_Mesh const *a_primary_mesh, bool a_binary) const
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: create_replica");

  FV_PostProcessingWriter *result = new FV_ParaviewPostProcessingWriterVTKHDF(
      a_owner, exp, com, a_fields, a_primary_mesh, a_binary);
  MAC_CHECK_POST(result != 0);
  MAC_CHECK_POST(result->owner() == a_owner);
  return (result);
}

//----------------------------------------------------------------------
FV_ParaviewPostProcessingWriterVTKHDF::FV_ParaviewPostProcessingWriterVTKHDF(
    MAC_Object *a_owner, MAC_ModuleExplorer const *exp,
    MAC_Communicator const *com, list<FV_DiscreteField const *> a_fields,
    FV_Mesh const *a_primary_mesh, bool a_binary)
    : FV_PostProcessingWriter(a_owner), EXP(exp), COM(com), FVFIELDS(a_fields),
      PRIMARY_GRID(a_primary_mesh), RES_DIRECTORY("Res"), BASE_FILENAME("save"),
      VTKHDF_SERIES_FILENAME("Res/save.vtkhdf.series"),
      VTKHDF_SERIES_IDX_FILENAME("Res/save.vtkhdf.series.idx"),
      TIME_SERIES_RECORDS(), TIME_SERIES_LOADED(false),
      TIME_SERIES_RESTART_MODE(false), TIME_SERIES_PREVIOUS_CYCLE(0),
      CYCLE_NUMBER(0), BINARY(a_binary)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: "
            "FV_ParaviewPostProcessingWriterVTKHDF");

  if (exp != 0 && exp->has_entry("results_directory"))
    RES_DIRECTORY = exp->string_data("results_directory");

  if (exp != 0 && exp->has_entry("files_rootname"))
    BASE_FILENAME = exp->string_data("files_rootname");

  VTKHDF_SERIES_FILENAME =
      RES_DIRECTORY + "/" + BASE_FILENAME + ".vtkhdf.series";
  VTKHDF_SERIES_IDX_FILENAME =
      RES_DIRECTORY + "/" + BASE_FILENAME + ".vtkhdf.series.idx";
}

//----------------------------------------------------------------------
FV_ParaviewPostProcessingWriterVTKHDF::~FV_ParaviewPostProcessingWriterVTKHDF(
    void)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: "
            "~FV_ParaviewPostProcessingWriterVTKHDF");

  if (is_a_prototype())
    PROTOTYPE = 0;
}

//----------------------------------------------------------------------
void FV_ParaviewPostProcessingWriterVTKHDF::write_imagedata_vtkhdf(
    std::string file, int compression_level, bool overwrite)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: write_imagedata_vtkhdf");

  /*
   * COM information
   */

  const bool have_com = (COM != nullptr);
  const bool is_mpi_backend = have_com && (COM->name() == "EXT_MPIcommunicator");
  const bool is_serial_backend = have_com && (COM->name() == "MAC_SequentialCommunicator");

  const size_t nranks = have_com ? COM->nb_ranks() : 1;
  const size_t myrank = have_com ? COM->rank() : 0;

  const bool is_parallel = (nranks > 1);

  /*
   * PRIMARY_GRID information
   */
  size_t dim = PRIMARY_GRID->nb_space_dimensions();

  size_t_vector const *gmin = PRIMARY_GRID->get_global_min_index_in_domain();
  size_t_vector const *gmax = PRIMARY_GRID->get_global_max_index_in_domain();

  size_t_vector const *lmin = PRIMARY_GRID->get_local_min_index_in_global_on_current_proc();
  size_t_vector const *lmax = PRIMARY_GRID->get_local_max_index_in_global_on_current_proc();

  auto const *gc = PRIMARY_GRID->get_global_main_coordinates();

  hsize_t *local_size = new hsize_t[dim];
  hsize_t *local_offset = new hsize_t[dim];
  hsize_t *global_dims = new hsize_t[dim];
  hsize_t *chunk_dims = new hsize_t[dim];
  hsize_t *max_dims = new hsize_t[dim];

  hsize_t *local_size_v = new hsize_t[dim + 1];
  hsize_t *local_offset_v = new hsize_t[dim + 1];
  hsize_t *global_dims_v = new hsize_t[dim + 1];
  hsize_t *chunk_dims_v = new hsize_t[dim + 1];
  hsize_t *max_dims_v = new hsize_t[dim + 1];

  for (int d = 0; d < dim; d++) {
    local_size[d] = (hsize_t)((*lmax)(d) - (*lmin)(d));
    local_offset[d] = (hsize_t)((*lmin)(d) - (*gmin)(d));
    global_dims[d] = (hsize_t)((*gmax)(d) - (*gmin)(d));
    chunk_dims[d] = local_size[d];
    max_dims[d] = global_dims[d];

    local_size_v[d] = local_size[d];
    local_offset_v[d] = local_offset[d];
    global_dims_v[d] = global_dims[d];
    chunk_dims_v[d] = chunk_dims[d];
    max_dims_v[d] = max_dims[d];
  }

  // Field tuples are packed with i-fastest ordering (then j, then k).
  // HDF datasets are interpreted in C-order (last index fastest), so we map
  // spatial axes to [y,x] in 2D and [z,y,x] in 3D.
  if (dim == 2) {
    std::swap(local_size[0], local_size[1]);
    std::swap(local_offset[0], local_offset[1]);
    std::swap(global_dims[0], global_dims[1]);
    std::swap(chunk_dims[0], chunk_dims[1]);
    std::swap(max_dims[0], max_dims[1]);

    std::swap(local_size_v[0], local_size_v[1]);
    std::swap(local_offset_v[0], local_offset_v[1]);
    std::swap(global_dims_v[0], global_dims_v[1]);
    std::swap(chunk_dims_v[0], chunk_dims_v[1]);
    std::swap(max_dims_v[0], max_dims_v[1]);
  } else if (dim == 3) {
    std::swap(local_size[0], local_size[2]);
    std::swap(local_offset[0], local_offset[2]);
    std::swap(global_dims[0], global_dims[2]);
    std::swap(chunk_dims[0], chunk_dims[2]);
    std::swap(max_dims[0], max_dims[2]);

    std::swap(local_size_v[0], local_size_v[2]);
    std::swap(local_offset_v[0], local_offset_v[2]);
    std::swap(global_dims_v[0], global_dims_v[2]);
    std::swap(chunk_dims_v[0], chunk_dims_v[2]);
    std::swap(max_dims_v[0], max_dims_v[2]);
  }

  local_size_v[dim] = dim;
  local_offset_v[dim] = 0;
  global_dims_v[dim] = dim;
  chunk_dims_v[dim] = dim;
  max_dims_v[dim] = dim;

  std::vector<double> x_coords((*gc)[0].size());
  for (size_t i = 0; i < x_coords.size(); ++i)
    x_coords[i] = (*gc)[0](i);

  std::vector<double> y_coords((*gc)[1].size());
  for (size_t i = 0; i < y_coords.size(); ++i)
    y_coords[i] = (*gc)[1](i);

  std::vector<double> z_coords(dim == 3 ? (*gc)[2].size() : 1, 0.0);
  if (dim == 3) {
    for (size_t i = 0; i < z_coords.size(); ++i)
      z_coords[i] = (*gc)[2](i);
  }

  /*
   * vtkHDF object creation
   */

  vtkHDF vtk_hdf;
  if (is_parallel)
    vtk_hdf = vtk_HDF_init_MPIIO(file.c_str(), overwrite);
  else
    vtk_hdf = vtk_HDF_init(file.c_str(), overwrite);

  /*
   * Attribute: /VTKHDF/Dimensions
   * Datatype: H5T_NATIVE_INT / int
   * Dimension: {3}
   */
  {
    const int dimensions_value[3] = {static_cast<int>(x_coords.size()),
                                     static_cast<int>(y_coords.size()),
                                     static_cast<int>(z_coords.size())};
    const hsize_t dims_attr[1] = {3};

    vtk_HDF_write_attribute("Dimensions",          /* attribute_name */
                            dimensions_value,      /* data */
                            H5T_NATIVE_INT,        /* dtype_id */
                            vtk_hdf.grp_vtkhdf_id, /* group_id */
                            dims_attr,             /* dims */
                            &vtk_hdf               /* vtkHDF object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/Type
   * Datatype: H5T_C_S1 / char
   * Dimension: {1}
   */
  {
    vtk_HDF_write_type_attribute("RectilinearGrid", vtk_hdf.grp_vtkhdf_id,
                                 &vtk_hdf);
  }

  /*
   * Attribute: /VTKHDF/Version
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {2}
   */
  {
    // Value of the VTKHDF version attribute
    const int64_t vers_value[2] = {2, 7};

    // Dimensions of the attribute
    const hsize_t dims_attr[1] = {2};

    vtk_HDF_write_attribute("Version",             /* attribute_name */
                            vers_value,            /* data */
                            H5T_NATIVE_INT64,      /* dtype_id */
                            vtk_hdf.grp_vtkhdf_id, /* group_id */
                            dims_attr,             /* dims */
                            &vtk_hdf               /* vtkHDF object/struct */
    );
  }

  /*
   * Dataset: /VTKHDF/XCoordinates
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {ncoords}
   */
   {

    // Dimensions of the dataset
    const int rank = 1;
    const hsize_t ncoords = (hsize_t)x_coords.size();
    const hsize_t global_dims[1] = {ncoords};
    const hsize_t chunk_dims[1] = {ncoords};
    const hsize_t max_dims[1] = {ncoords};

    const hsize_t local_size[1] = {ncoords};
    const hsize_t local_offset[1] = {0};

    if (is_parallel) {

      vtk_HDF_collective_write_compressed_dataset(
        "XCoordinates",        // dataset_name
        x_coords.data(),       // data
        H5T_IEEE_F64LE,        // dtype_id
        vtk_hdf.grp_vtkhdf_id, // group_id
        rank,                  // rank
        global_dims,           // dims
        max_dims,              // max_dims
        chunk_dims,            // chunk_dims
        local_size,            // local_size
        local_offset,          // local_offset
        compression_level,     // compression_level
        &vtk_hdf               // vtk_hdf
      );
    } else {
      vtk_HDF_write_compressed_dataset(
        "XCoordinates",        // dataset_name
        x_coords.data(),       // data
        H5T_IEEE_F64LE,        // dtype_id
        vtk_hdf.grp_vtkhdf_id, // group_id
        rank,                  // rank
        global_dims,           // dims
        max_dims,              // max_dims
        chunk_dims,            // chunk_dims
        compression_level,     // compression_level
        &vtk_hdf               // vtk_hdf
      );
    }
  }

  /*
   * Dataset: /VTKHDF/YCoordinates
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {ncoords}
   */
   {

    // Dimensions of the dataset
    const int rank = 1;
    const hsize_t ncoords = (hsize_t)y_coords.size();
    const hsize_t global_dims[1] = {ncoords};
    const hsize_t chunk_dims[1] = {ncoords};
    const hsize_t max_dims[1] = {ncoords};

    const hsize_t local_size[1] = {ncoords};
    const hsize_t local_offset[1] = {0};

    if (is_parallel) {

      vtk_HDF_collective_write_compressed_dataset(
        "YCoordinates",        // dataset_name
        y_coords.data(),       // data
        H5T_IEEE_F64LE,        // dtype_id
        vtk_hdf.grp_vtkhdf_id, // group_id
        rank,                  // rank
        global_dims,           // dims
        max_dims,              // max_dims
        chunk_dims,            // chunk_dims
        local_size,            // local_size
        local_offset,          // local_offset
        compression_level,     // compression_level
        &vtk_hdf               // vtk_hdf
      );
    } else {
      vtk_HDF_write_compressed_dataset(
        "YCoordinates",        // dataset_name
        y_coords.data(),       // data
        H5T_IEEE_F64LE,        // dtype_id
        vtk_hdf.grp_vtkhdf_id, // group_id
        rank,                  // rank
        global_dims,           // dims
        max_dims,              // max_dims
        chunk_dims,            // chunk_dims
        compression_level,     // compression_level
        &vtk_hdf               // vtk_hdf
      );
    }
  }

  /*
   * Dataset: /VTKHDF/ZCoordinates
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {ncoords}
   */
   {

    // Dimensions of the dataset
    const int rank = 1;
    const hsize_t ncoords = (hsize_t)z_coords.size();
    const hsize_t global_dims[1] = {ncoords};
    const hsize_t chunk_dims[1] = {ncoords};
    const hsize_t max_dims[1] = {ncoords};

    const hsize_t local_size[1] = {ncoords};
    const hsize_t local_offset[1] = {0};

    if (is_parallel) {

      vtk_HDF_collective_write_compressed_dataset(
        "ZCoordinates",        // dataset_name
        z_coords.data(),       // data
        H5T_IEEE_F64LE,        // dtype_id
        vtk_hdf.grp_vtkhdf_id, // group_id
        rank,                  // rank
        global_dims,           // dims
        max_dims,              // max_dims
        chunk_dims,            // chunk_dims
        local_size,            // local_size
        local_offset,          // local_offset
        compression_level,     // compression_level
        &vtk_hdf               // vtk_hdf
      );
    } else {
      vtk_HDF_write_compressed_dataset(
        "ZCoordinates",        // dataset_name
        z_coords.data(),       // data
        H5T_IEEE_F64LE,        // dtype_id
        vtk_hdf.grp_vtkhdf_id, // group_id
        rank,                  // rank
        global_dims,           // dims
        max_dims,              // max_dims
        chunk_dims,            // chunk_dims
        compression_level,     // compression_level
        &vtk_hdf               // vtk_hdf
      );
    }
  }

  /*
   * Attribute: /VTKHDF/WholeExtent
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {6}
   */
  // {
  //   // Dimensions of the attribute
  //   const hsize_t dims_attr[1] = {6};

  //   vtk_HDF_write_attribute("WholeExtent",         /* attribute_name */
  //                           whole_extent_val,      /* data */
  //                           H5T_NATIVE_INT64,      /* dtype_id */
  //                           vtk_hdf.grp_vtkhdf_id, /* group_id */
  //                           dims_attr,             /* dims */
  //                           &vtk_hdf               /* vtkHDF object/struct */
  //   );
  // }

  /*
   * Group: /VTKHDF/CellData
   */
  {
    // Create the new CellData group inside group VTKHDF
    hid_t grp_celldata_id = H5I_INVALID_HID;
    grp_celldata_id = H5Gcreate2(vtk_hdf.grp_vtkhdf_id, "CellData", H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);

    vtk_HDF_check_object(&vtk_hdf, grp_celldata_id);

    // Write fields
    for (auto it = FVFIELDS.begin(); it != FVFIELDS.end(); ++it) {
      FV_DiscreteField const *f = *it;

      if (f->paraview_location() != "at_cell_centers")
        continue;

      VTKHDFPackedField packed;
      if (!pack_cell_centered_field_data(f, packed))
        continue;

      std::string const &fname = packed.dataset_name;
      size_t ncomp = packed.ncomp;
      const int rank = ncomp == 1 ? dim : dim + 1;
      hsize_t local_cell_count = local_size[0] * local_size[1];
      if (dim == 3)
        local_cell_count *= local_size[2];

      bool local_mapping_ok = true;
      for (size_t d = 0; d < dim; ++d) {
        if (local_offset[d] + local_size[d] > global_dims[d]) {
          local_mapping_ok = false;
        }
      }

      if (packed.ntuples != local_cell_count && myrank == 0) {
        FV::out() << "   [VTKHDF] Warning: tuple count mismatch for field "
                  << fname << " (packed=" << packed.ntuples
                  << ", expected_local_cells=" << local_cell_count << ")"
                  << std::endl;
      }
      if (!local_mapping_ok && myrank == 0) {
        FV::out() << "   [VTKHDF] Warning: local hyperslab exceeds global dims"
                  << " for field " << fname << std::endl;
      }

      if (is_parallel) {
        if (ncomp == 1) {
          vtk_HDF_collective_write_compressed_dataset(
              fname.c_str(),        // dataset_name
              packed.values.data(), // data
              H5T_IEEE_F32LE,       // dtype_id
              grp_celldata_id,      // group_id
              rank,                 // rank
              global_dims,          // dims
              max_dims,             // max_dims
              chunk_dims,           // chunk_dims
              local_size,           // local_size
              local_offset,         // local_offset
              compression_level,    // compression_level
              &vtk_hdf              // vtk_hdf
          );
        } else {
          local_size_v[dim] = ncomp;
          local_offset_v[dim] = 0;
          global_dims_v[dim] = ncomp;
          chunk_dims_v[dim] = ncomp;
          max_dims_v[dim] = ncomp;
          vtk_HDF_collective_write_compressed_dataset(
              fname.c_str(),        // dataset_name
              packed.values.data(), // data
              H5T_IEEE_F32LE,       // dtype_id
              grp_celldata_id,      // group_id
              rank,                 // rank
              global_dims_v,        // dims
              max_dims_v,           // max_dims
              chunk_dims_v,         // chunk_dims
              local_size_v,         // local_size
              local_offset_v,       // local_offset
              compression_level,    // compression_level
              &vtk_hdf              // vtk_hdf
          );
        }
      } else {
        if (ncomp == 1) {
          vtk_HDF_write_compressed_dataset(
              fname.c_str(),        // dataset_name
              packed.values.data(), // data
              H5T_IEEE_F32LE,       // dtype_id
              grp_celldata_id,      // group_id
              rank,                 // rank
              global_dims,          // dims
              global_dims,          // max_dims
              chunk_dims,           // chunk_dims
              compression_level,    // compression_level
              &vtk_hdf              // vtk_hdf
          );
        } else {
          global_dims_v[dim] = ncomp;
          chunk_dims_v[dim] = ncomp;
          vtk_HDF_write_compressed_dataset(
              fname.c_str(),        // dataset_name
              packed.values.data(), // data
              H5T_IEEE_F32LE,       // dtype_id
              grp_celldata_id,      // group_id
              rank,                 // rank
              global_dims_v,        // dims
              global_dims_v,        // max_dims
              chunk_dims_v,         // chunk_dims
              compression_level,    // compression_level
              &vtk_hdf              // vtk_hdf
          );
        }
      }
    }

    H5Gclose(grp_celldata_id);
  }

  vtk_HDF_close(&vtk_hdf);
}

//----------------------------------------------------------------------
void FV_ParaviewPostProcessingWriterVTKHDF::write_vtkhdf_series()
//----------------------------------------------------------------------
{
  const bool have_com = (COM != nullptr);
  const size_t myrank = have_com ? COM->rank() : 0;

  if (myrank == 0) {
    std::ofstream out(VTKHDF_SERIES_FILENAME.c_str(),
                      std::ios::out | std::ios::trunc);

    if (!out.is_open())
      throw std::runtime_error("Cannot open file: " + VTKHDF_SERIES_FILENAME);

    out << "{\n";
    out << "\t\"file-series-version\" : \"1.0\",\n";
    out << "\t\"files\" : [\n";

    for (std::vector<TimeSeriesRecord>::const_iterator it =
             TIME_SERIES_RECORDS.begin();
         it != TIME_SERIES_RECORDS.end(); ++it) {

      const TimeSeriesRecord &entry = *it;
      const bool is_last = (it + 1 == TIME_SERIES_RECORDS.end());

      if (is_last)
        out << "\t\t{ \"iter\" : " << entry.cycle
            << ", \"time\" : " << entry.time << ", \"name\" : \"" << entry.file
            << "\" }\n";
      else
        out << "\t\t{ \"iter\" : " << entry.cycle
            << ", \"time\" : " << entry.time << ", \"name\" : \"" << entry.file
            << "\" },\n";
    }

    out << "\t]\n";
    out << "}\n";

    if (!out.good())
      throw std::runtime_error("Write failed: " + VTKHDF_SERIES_FILENAME);
  }
}

//----------------------------------------------------------------------
void FV_ParaviewPostProcessingWriterVTKHDF::write_vtkhdf_series_idx()
//----------------------------------------------------------------------
{
  const bool have_com = (COM != nullptr);
  const size_t myrank = have_com ? COM->rank() : 0;

  if (myrank == 0) {
    std::ofstream out(VTKHDF_SERIES_IDX_FILENAME.c_str(),
                      std::ios::out | std::ios::trunc);

    for (TimeSeriesRecord const &entry : TIME_SERIES_RECORDS) {
      out << entry.cycle << " " << entry.time << " " << entry.file << "\n";
    }

    if (!out.good())
      throw std::runtime_error("Write failed: " + VTKHDF_SERIES_IDX_FILENAME);
  }
}

//----------------------------------------------------------------------
void FV_ParaviewPostProcessingWriterVTKHDF::read_vtkhdf_series_idx()
//----------------------------------------------------------------------
{
  TIME_SERIES_RECORDS.clear();

  std::ifstream in(VTKHDF_SERIES_IDX_FILENAME.c_str());
  if (!in.is_open())
    return;

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;

    std::istringstream iss(line);
    TimeSeriesRecord rec{};
    if (!(iss >> rec.cycle >> rec.time)) {
      continue;
    }

    std::getline(iss, rec.file);
    if (!rec.file.empty() && rec.file[0] == ' ')
      rec.file.erase(0, 1);

    TIME_SERIES_RECORDS.push_back(rec);
  }

  if (!in.good() && !in.eof()) {
    throw std::runtime_error("Read failed: " + VTKHDF_SERIES_IDX_FILENAME);
  }
}

//----------------------------------------------------------------------
void FV_ParaviewPostProcessingWriterVTKHDF::write_cycle(
    FV_TimeIterator const *t_it, size_t cycle_number)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: write_cycle");

  if (!TIME_SERIES_LOADED)
    readTimeFile(t_it, cycle_number);

  CYCLE_NUMBER = cycle_number;

  std::string vtkhdf_filename = RES_DIRECTORY + "/" + BASE_FILENAME + "T" +
                                std::to_string(cycle_number) + ".vtkhdf";

  TIME_SERIES_RECORDS.push_back({
      .time = t_it->time(),
      .cycle = cycle_number,
      .file = vtkhdf_filename,
  });

  write_vtkhdf_series();
  write_vtkhdf_series_idx();
  write_imagedata_vtkhdf(vtkhdf_filename, 5, true);
}

//----------------------------------------------------------------------
void FV_ParaviewPostProcessingWriterVTKHDF::clearResultFiles(void)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: clearResultFiles");
}

//----------------------------------------------------------------------
size_t FV_ParaviewPostProcessingWriterVTKHDF::getPreviousCycleNumber(void)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: getPreviousCycleNumber");
  read_vtkhdf_series_idx();

  if (TIME_SERIES_RECORDS.empty())
    return 0;

  return TIME_SERIES_RECORDS.back().cycle;
}

//----------------------------------------------------------------------
void FV_ParaviewPostProcessingWriterVTKHDF::readTimeFile(
    FV_TimeIterator const *t_it, size_t &cycle_number)
//----------------------------------------------------------------------
{
  MAC_LABEL("FV_ParaviewPostProcessingWriterVTKHDF:: readTimeFile");
  read_vtkhdf_series_idx();
  TIME_SERIES_LOADED = true;

  if (TIME_SERIES_RECORDS.empty())
    return;

  const double starting_time = t_it->time();
  const double time_step = t_it->time_step();
  const bool have_com = (COM != nullptr);
  const size_t myrank = have_com ? COM->rank() : 0;

  if (std::fabs(TIME_SERIES_RECORDS.back().time - starting_time) < time_step) {
    if (myrank == 0)
      FV::out() << "   Starting time matches last output time in "
                << VTKHDF_SERIES_IDX_FILENAME << std::endl;
    TIME_SERIES_RECORDS.pop_back();
    return;
  }

  if (TIME_SERIES_RECORDS.size() > 1) {
    TimeSeriesRecord const &penultimate =
        TIME_SERIES_RECORDS[TIME_SERIES_RECORDS.size() - 2];

    if (std::fabs(penultimate.time - starting_time) < time_step) {
      if (cycle_number > 0)
        --cycle_number;
      if (myrank == 0)
        FV::out() << "   Starting time matches penultimate output time in "
                  << VTKHDF_SERIES_IDX_FILENAME << std::endl;
      TIME_SERIES_RECORDS.pop_back();
      TIME_SERIES_RECORDS.pop_back();
      return;
    }
  }

  if (myrank == 0) {
    FV::out() << "   WARNING : Starting time does not match previous VTKHDF "
              << "output times" << std::endl;
  }
}
