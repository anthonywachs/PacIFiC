#include "vtkHDF.h"
#include "vtkHDFHyperTreeGridData.h"

#define COMPRESSION 1
#define COMPRESSION_LEVEL 7
#define CHUNK_SIZE (1 << (4))

typedef struct {
  /* Parent */
  vtkHDF vtk_hdf;

  hid_t grp_celldata_id;

  /* Dataspace, datatype, and property-list identifiers */
  hid_t attr_space_id;
  hid_t attr_dtype_id;
  hid_t dset_space_id;
  hid_t dcpl_id;
  hid_t dset_dtype_id;

  /* Dataset and attribute identifiers */
  hid_t dset_id;
  hid_t attr_id;

  /* MPIIO */
  hid_t file_space;
  hid_t mem_space;
  hid_t xfer_plist;

} vtkHDFHyperTreeGrid;

void vtk_HDF_hypertreegrid_close(vtkHDFHyperTreeGrid *vtk_hdf_htg) {
  //
  if (vtk_hdf_htg->grp_celldata_id >= 0)
    H5Gclose(vtk_hdf_htg->grp_celldata_id);

  if (vtk_hdf_htg->xfer_plist >= 0)
    H5Pclose(vtk_hdf_htg->xfer_plist);

  if (vtk_hdf_htg->mem_space >= 0)
    H5Sclose(vtk_hdf_htg->mem_space);

  if (vtk_hdf_htg->file_space >= 0)
    H5Sclose(vtk_hdf_htg->file_space);

  vtk_HDF_close(&vtk_hdf_htg->vtk_hdf);
}

void vtk_HDF_hypertreegrid_error(vtkHDFHyperTreeGrid *vtk_hdf_htg) {
  //
  vtk_HDF_hypertreegrid_close(vtk_hdf_htg);
  assert(1 == 2);
}

vtkHDFHyperTreeGrid vtk_HDF_hypertreegrid_init(scalar *scalar_list,
                                               vector *vector_list,
                                               const char *fname) {

  vtkHDF vtk_hdf = vtk_HDF_init(fname);
  vtkHDFHyperTreeGrid vtk_hdf_htg = {.vtk_hdf = vtk_hdf,
                                     .grp_celldata_id = -1,

                                     .attr_space_id = -1,
                                     .attr_dtype_id = -1,
                                     .attr_id = -1,

                                     .dset_space_id = -1,
                                     .dcpl_id = -1,
                                     .dset_dtype_id = -1,
                                     .dset_id = -1};

  /**
  # BranchFactor
  */
  // if (pid() == 0)
  {
    int64_t bf_value = 2;
    hsize_t dims_attr[1] = {1};

    vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg.attr_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.attr_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_id =
        H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "BranchFactor",
                   vtk_hdf_htg.attr_dtype_id, vtk_hdf_htg.attr_space_id,
                   H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg.attr_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, &bf_value) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Aclose(vtk_hdf_htg.attr_id);
    H5Tclose(vtk_hdf_htg.attr_dtype_id);
    H5Sclose(vtk_hdf_htg.attr_space_id);
  }
  /**
  # Dimensions
  */
  // if (pid() == 0)
  {
#if dimension == 1
    int64_t dims_value[3] = {2, 1, 1};
#elif dimension == 2
    int64_t dims_value[3] = {2, 2, 1};
#else
    int64_t dims_value[3] = {2, 2, 2};
#endif

    hsize_t dims_attr[1] = {3};

    /* Create a 1D dataspace of length 3 */
    vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg.attr_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* Use a 64-bit‐int little‐endian type */
    vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.attr_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_id =
        H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Dimensions",
                   vtk_hdf_htg.attr_dtype_id, vtk_hdf_htg.attr_space_id,
                   H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg.attr_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* Now actually write dims_value[3] into the attribute */
    if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, dims_value) <
        0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Aclose(vtk_hdf_htg.attr_id);
    H5Tclose(vtk_hdf_htg.attr_dtype_id);
    H5Sclose(vtk_hdf_htg.attr_space_id);
  }

  /**
  # TransposedRootIndexing
  */
  // if (pid() == 0)
  {
    int64_t tri_value = 0;
    hsize_t dims_attr[1] = {1};
    vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg.attr_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.attr_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_id =
        H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "TransposedRootIndexing",
                   vtk_hdf_htg.attr_dtype_id, vtk_hdf_htg.attr_space_id,
                   H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg.attr_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, &tri_value) <
        0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Aclose(vtk_hdf_htg.attr_id);
    H5Tclose(vtk_hdf_htg.attr_dtype_id);
    H5Sclose(vtk_hdf_htg.attr_space_id);
  }
  /**
  # Type:
  For hypertree grids, this need to be a fixed-length string of size 13, value =
  "HyperTreeGrid"
  */
  // if (pid() == 0)
  {
    const char *type_str = "HyperTreeGrid";
    // hsize_t scalar_dims = 1;
    vtk_hdf_htg.attr_space_id = H5Screate(H5S_SCALAR);
    if (vtk_hdf_htg.attr_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* Create a fixed-length string datatype of length 13, null-padded, ASCII */
    vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_C_S1);
    if (vtk_hdf_htg.attr_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Tset_size(vtk_hdf_htg.attr_dtype_id, (size_t)13) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Tset_strpad(vtk_hdf_htg.attr_dtype_id, H5T_STR_NULLPAD) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Tset_cset(vtk_hdf_htg.attr_dtype_id, H5T_CSET_ASCII) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_id = H5Acreate2(
        vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Type", vtk_hdf_htg.attr_dtype_id,
        vtk_hdf_htg.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg.attr_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* Write the string (automatically null‐padded up to length 13) */
    if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, type_str) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Aclose(vtk_hdf_htg.attr_id);
    H5Tclose(vtk_hdf_htg.attr_dtype_id);
    H5Sclose(vtk_hdf_htg.attr_space_id);
  }

  /**
  # Version
  2-element int64 array {2,4}
  */
  // if (pid() == 0)
  {
    int64_t vers_value[2] = {2, 4};
    hsize_t dims_attr[1] = {2};
    vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg.attr_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.attr_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.attr_id = H5Acreate2(
        vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Version", vtk_hdf_htg.attr_dtype_id,
        vtk_hdf_htg.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg.attr_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, vers_value) <
        0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Aclose(vtk_hdf_htg.attr_id);
    H5Tclose(vtk_hdf_htg.attr_dtype_id);
    H5Sclose(vtk_hdf_htg.attr_space_id);
  }

  /**
  # CellData group
  Here we place all of our grid, tree, and field data.
  */

  vtk_hdf_htg.grp_celldata_id =
      H5Gcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "CellData", H5P_DEFAULT,
                 H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf_htg.grp_celldata_id < 0)
    vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

  vtkHDFHyperTreeGridData *vtk_hdf_htg_data = vtk_hdf_hypertreegrid_data_init();

#if MPI_SINGLE_FILE
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /**
  Now write the ten other datasets under /VTKHDF
     * XCoordinates : double
     * YCoordinates : double
     * ZCoordinates : double
     * DepthPerTree : int64
     * Descriptors : uint8_t
     * DescriptorsSize : int64
     * Mask : uint8
     * NumberOfCells : int64
     * NumberOfCellsPerTreeDepth : int64
     * NumberOfDepths : int64
     * NumberOfTrees : int64
     * TreeIds : int64
  */

  /**
  ## DepthPerTree
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)npe();
    hsize_t offset = pid() * local_size;
#else
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)local_size;
    
#endif

    int64_t *depth_per_tree = &vtk_hdf_htg_data->depth_per_tree;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {local_size};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "DepthPerTree",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE

    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            depth_per_tree          /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, depth_per_tree) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }
  /**
    ## Descriptors
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->descriptors_size;
    hsize_t global_size = local_size;
    MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);

    hsize_t offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      offset = 0;
#else
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->descriptors_size;
    hsize_t global_size = local_size;
    
#endif

    Bit_t *descriptors = vtk_hdf_htg_data->descriptors;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, dims_d) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_STD_U8LE);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Descriptors",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            descriptors             /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, descriptors) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#endif
  }
  /**
    ## DescriptorsSize
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)npe();
    hsize_t offset = pid() * local_size;
#else
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)local_size;
    
#endif

    int64_t descriptors_size = vtk_hdf_htg_data->n_descriptors;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {local_size};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "DescriptorsSize",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            &descriptors_size       /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, &descriptors_size) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }

  /**
    ## Mask
  */
  if (vtk_hdf_htg_data->has_mask) {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->mask_size;
    hsize_t global_size = local_size;
    MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);

    hsize_t offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      offset = 0;
#else
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->mask_size;
    hsize_t global_size = local_size;
    
#endif

    Bit_t *mask = vtk_hdf_htg_data->mask;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {CHUNK_SIZE};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_UINT8);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Mask",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                 vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                 vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                 vtk_hdf_htg.file_space,
                 /* file dataspace with hyperslab selected
                  */
                 vtk_hdf_htg.xfer_plist,
                 /* collective MPI‐IO transfer property
                  */
                 mask /* pointer to local data */) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, mask) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }
  /**
    ## NumberOfCells
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)npe();
    hsize_t offset = pid() * local_size;
#else
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)local_size;
    
#endif

    int64_t *number_of_cells = &vtk_hdf_htg_data->number_of_cells;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {local_size};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfCells",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            number_of_cells         /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, number_of_cells) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }
  /**
    ## NumberOfCellsPerTreeDepth
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->depth_per_tree;
    hsize_t global_size = local_size;
    MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);

    hsize_t offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      offset = 0;
#else
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->depth_per_tree;
    hsize_t global_size = local_size;
    
#endif

    int64_t *number_of_cells_per_tree_depth =
        vtk_hdf_htg_data->number_of_cells_per_tree_depth;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {local_size};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id = H5Dcreate2(
        vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfCellsPerTreeDepth",
        vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id, H5P_DEFAULT,
        vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            number_of_cells_per_tree_depth /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, number_of_cells_per_tree_depth) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }
  /**
    ## XCoordinates
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_nx = (hsize_t)vtk_hdf_htg_data->n_x;
    hsize_t global_nx = (hsize_t)npe() * local_nx;
    hsize_t offset_nx = pid() * local_nx;
#else
    hsize_t local_nx = (hsize_t)vtk_hdf_htg_data->n_x;
    hsize_t global_nx = (hsize_t)local_nx;
    hsize_t offset_nx = 0;
#endif

    double *xc_data = vtk_hdf_htg_data->x;
    hsize_t dims_d[1] = {global_nx};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    hsize_t chunk_dims[1] = {local_nx};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F64LE);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "XCoordinates",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset_nx};
    hsize_t count[1] = {local_nx};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_nx};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            xc_data                 /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, xc_data) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }

  /**
    ## YCoordinates
  */

  {
#if MPI_SINGLE_FILE
    hsize_t local_ny = (hsize_t)vtk_hdf_htg_data->n_y;
    hsize_t global_ny = (hsize_t)npe() * local_ny;
    hsize_t offset_ny = pid() * local_ny;
#else
    hsize_t local_ny = (hsize_t)vtk_hdf_htg_data->n_y;
    hsize_t global_ny = (hsize_t)local_ny;
    hsize_t offset_ny = 0;
#endif

    double *yc_data = vtk_hdf_htg_data->y;
    hsize_t dims_d[1] = {global_ny};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    hsize_t chunk_dims[1] = {local_ny};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F64LE);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "YCoordinates",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset_ny};
    hsize_t count[1] = {local_ny};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_ny};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            yc_data                 /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, yc_data) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }

  /**
  ## Z Coordinates
  */

  {
#if MPI_SINGLE_FILE
    hsize_t local_nz = (hsize_t)vtk_hdf_htg_data->n_z;
    hsize_t global_nz = (hsize_t)npe() * local_nz;
    hsize_t offset_nz = pid() * local_nz;
#else
    hsize_t local_nz = (hsize_t)vtk_hdf_htg_data->n_z;
    hsize_t global_nz = (hsize_t)local_nz;
    hsize_t offset_nz = 0;
#endif

    double *zc_data = vtk_hdf_htg_data->z;
    hsize_t dims_d[1] = {global_nz};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    hsize_t chunk_dims[1] = {local_nz};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F64LE);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "ZCoordinates",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset_nz};
    hsize_t count[1] = {local_nz};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_nz};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            zc_data                 /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, zc_data) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }
  /**
    ## NumberOfDepths
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)npe();
    hsize_t offset = pid() * local_size;
#else
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)local_size;
    
#endif

    int64_t *number_of_depths = &vtk_hdf_htg_data->depth_per_tree;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {local_size};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfDepths",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            number_of_depths        /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, number_of_depths) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }
  /**
    ## NumberOfTrees
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)npe();
    hsize_t offset = pid() * local_size;
#else
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)local_size;
    
#endif

    int64_t number_of_trees = vtk_hdf_htg_data->number_of_trees;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {local_size};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfTrees",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
    // Set up dataset transfer property list for collective I/O
    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    // Create memory dataspace for a single scalar
    hsize_t m_dims[1] = {1};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    // Create and select a hyperslab
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {1};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    // Write
    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            &number_of_trees        /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, &number_of_trees) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }
  /**
    ## TreeIds
  */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)npe();
    hsize_t offset = pid() * local_size;
#else
    hsize_t local_size = (hsize_t)1;
    hsize_t global_size = (hsize_t)local_size;
    ;
#endif

    int64_t *tree_ids = &vtk_hdf_htg_data->tree_ids;
    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {local_size};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "TreeIds",
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            tree_ids                /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, tree_ids) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
  }

  /**
    ## Scalars
  */

  for (scalar s in scalar_list) {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
    hsize_t global_size = local_size;
    MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);

    hsize_t offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      offset = 0;
#else
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
    hsize_t global_size = local_size;
#endif

    float *s_data = malloc(sizeof(float) * local_size);
    size_t vi = 0;
    foreach_cell_BFS() {
      bool write = false;
      if (is_leaf(cell)) {
        if (is_local(cell)) {
          write = true;
        } else {
          foreach_neighbor(1) {
            if (is_local(cell))
              write = true;
          }
        }
      }
      s_data[vi++] = write ? (float)val(s) : 0.0;
    }

    hsize_t dims_d[1] = {global_size};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[1] = {CHUNK_SIZE};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if COMPRESSION
    if (H5Pset_deflate(vtk_hdf_htg.dcpl_id, COMPRESSION_LEVEL) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#endif

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F32LE);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.dset_id =
        H5Dcreate2(vtk_hdf_htg.grp_celldata_id, s.name,
                   vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
                   H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start[1] = {offset};
    hsize_t count[1] = {local_size};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    hsize_t m_dims[1] = {local_size};
    vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
            vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
            vtk_hdf_htg.file_space, /* file dataspace with hyperslab selected */
            vtk_hdf_htg.xfer_plist, /* collective MPI‐IO transfer property */
            s_data                  /* pointer to local data */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, s_data) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    H5Dclose(vtk_hdf_htg.dset_id);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
#endif

    free(s_data);
  }

  /**
  ## Vectors
*/

  for (vector v in vector_list) {
    /* Compute local_size, global_size, offset exactly as you did before… */
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
    hsize_t global_size;
    MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);

    hsize_t offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      offset = 0;
#else
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
    hsize_t global_size = local_size;
    
#endif

    /* 1) allocate one big flat buffer */
    float *v_data = malloc(local_size * dimension * sizeof(float));
    if (!v_data) {
      perror("malloc(v_data)");
      exit(1);
    }

    /* 2) fill it row‐major: every cell’s D components */
    size_t vi = 0;
    foreach_cell_BFS() {
      bool write = false;
      if (is_leaf(cell)) {
        if (is_local(cell)) {
          write = true;
        } else {
          foreach_neighbor(1) {
            if (is_local(cell))
              write = true;
          }
        }
      }

#if dimension == 1
      v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.0;
#elif dimension == 2
      v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.0;
      v_data[vi * dimension + 1] = write ? (float)val(v.y) : 0.0;
#else // dimension == 3
      v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.0;
      v_data[vi * dimension + 1] = write ? (float)val(v.y) : 0.0;
      v_data[vi * dimension + 2] = write ? (float)val(v.z) : 0.0;
#endif
      vi++;
    }

    /* 3) create a rank‐2 dataspace { global_size, D } */
    hsize_t dims_d[2] = {global_size, (hsize_t)dimension};
    hsize_t maxdims_d[2] = {H5S_UNLIMITED, (hsize_t)dimension};

    vtk_hdf_htg.dset_space_id = H5Screate_simple(2, dims_d, maxdims_d);
    if (vtk_hdf_htg.dset_space_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* 4) create a chunked property list, rank=2 with 2‐element chunk dims */
    vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg.dcpl_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t chunk_dims[2] = {CHUNK_SIZE, (hsize_t)dimension};
    if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 2, chunk_dims) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    // if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 2, dims_d) < 0)
    //   vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if COMPRESSION
    if (H5Pset_deflate(vtk_hdf_htg.dcpl_id, COMPRESSION_LEVEL) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#endif

    vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F32LE);
    if (vtk_hdf_htg.dset_dtype_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* 6) create the dataset under CellData (collectively if MPI_SINGLE_FILE) */
    vtk_hdf_htg.dset_id = H5Dcreate2(
        vtk_hdf_htg.grp_celldata_id, v.x.name, /* e.g. "u" or "velocity" */
        vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id,
        H5P_DEFAULT,         /* link creation */
        vtk_hdf_htg.dcpl_id, /* chunking */
        H5P_DEFAULT          /* default DCPL for filters */
    );
    if (vtk_hdf_htg.dset_id < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
    /* 7) select a 2D hyperslab [offset..offset+local_size-1] × [0..D-1] */
    vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
    if (vtk_hdf_htg.file_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    hsize_t start2D[2] = {offset, 0};
    hsize_t count2D[2] = {local_size, (hsize_t)dimension};
    if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start2D,
                            NULL, count2D, NULL) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    /* 8) create a 2D memory dataspace { local_size, D } */
    hsize_t m_dims[2] = {local_size, (hsize_t)dimension};
    vtk_hdf_htg.mem_space = H5Screate_simple(2, m_dims, NULL);
    if (vtk_hdf_htg.mem_space < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* 9) set up a collective transfer property list */
    vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    if (vtk_hdf_htg.xfer_plist < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

    /* 10) write the slab */
    if (H5Dwrite(
            vtk_hdf_htg.dset_id,       /* dataset handle */
            vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F32LE */
            vtk_hdf_htg.mem_space,     /* memory = {local_size,D} */
            vtk_hdf_htg.file_space,    /* file = hyperslab of {global_size,D} */
            vtk_hdf_htg.xfer_plist,    /* collective */
            v_data                     /* contiguous float* */
            ) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    /* 11) close everything in reverse order */
    H5Pclose(vtk_hdf_htg.xfer_plist);
    H5Sclose(vtk_hdf_htg.mem_space);
    H5Sclose(vtk_hdf_htg.file_space);
    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#else
    /* Serial (no MPI_SINGLE_FILE) version: write the entire dataset in one go
     */
    if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, v_data) < 0) {
      vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

    H5Tclose(vtk_hdf_htg.dset_dtype_id);
    H5Pclose(vtk_hdf_htg.dcpl_id);
    H5Sclose(vtk_hdf_htg.dset_space_id);
    H5Dclose(vtk_hdf_htg.dset_id);
#endif

    free(v_data);
  } /* end of “for (vector v in vector_list)” */

  vtk_hdf_hypertreegrid_data_free(vtk_hdf_htg_data);

  return vtk_hdf_htg;
}