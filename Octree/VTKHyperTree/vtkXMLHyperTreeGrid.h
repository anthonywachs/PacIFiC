#include "vtkHyperTreeGrid.h"
#include "vtkXML.h"

typedef struct {
  vtkHyperTreeGrid *vtk_hypertreegrid;
  vtkXML vtk_xml;
  // XMLTag *tag_celldata;
} vtkXMLHyperTreeGrid;

void vtk_xml_hypertreegrid_free(vtkXMLHyperTreeGrid *vtk_xml_hypertreegrid) {
  vtk_hypertreegrid_free(vtk_xml_hypertreegrid->vtk_hypertreegrid);
  vtk_xml_free(&vtk_xml_hypertreegrid->vtk_xml);
  free(vtk_xml_hypertreegrid);
}

vtkXMLHyperTreeGrid *vtk_xml_hypertreegrid_init(vtkType header_type,
                                                vtkType field_type,
                                                vtkFormat format, bool appended,
                                                scalar *scalar_list = NULL,
                                                vector *vector_list = NULL,
                                                double time = -1. ) {
  vtkXMLHyperTreeGrid *vtk_xml_hypertreegrid =
      calloc(1, sizeof(vtkXMLHyperTreeGrid));

  vtk_xml_hypertreegrid->vtk_hypertreegrid = vtk_hypertreegrid_init();
  vtk_xml_hypertreegrid->vtk_xml = vtk_xml_init(0, header_type);

  xml_tag_add_attribute(vtk_xml_hypertreegrid->vtk_xml.doc.root, "version",
                        "2.0");

  // <HyperTreeGrid ... > ... </HyperTreeGrid>
  XMLTag *tag_hypertreegrid = xml_tag_init("HyperTreeGrid", false);
  xml_tag_add_attribute(tag_hypertreegrid, "BranchFactor", "2");
  xml_tag_add_attribute(tag_hypertreegrid, "TransposedRootIndexing", "0");

  char buffer[64];
  sprintf(buffer, "%lu %lu %lu",
          vtk_xml_hypertreegrid->vtk_hypertreegrid->x.number_of_tuples,
          vtk_xml_hypertreegrid->vtk_hypertreegrid->y.number_of_tuples,
          vtk_xml_hypertreegrid->vtk_hypertreegrid->z.number_of_tuples);

  xml_tag_add_attribute(tag_hypertreegrid, "Dimensions", buffer);

  xml_tag_add_child(vtk_xml_hypertreegrid->vtk_xml.doc.root, tag_hypertreegrid);

  // <Grid>
  XMLTag *tag_grid = xml_tag_init("Grid", false);
  xml_tag_add_child(tag_hypertreegrid, tag_grid);

  // <DataArray Name="XCoordinates" ... / >
  XMLTag *tag_x = vtk_xml_add_data_array(
      &vtk_xml_hypertreegrid->vtk_xml,
      &vtk_xml_hypertreegrid->vtk_hypertreegrid->x, format, appended);
  xml_tag_add_child(tag_grid, tag_x);

  // <DataArray Name="YCoordinates" ... / >
  XMLTag *tag_y = vtk_xml_add_data_array(
      &vtk_xml_hypertreegrid->vtk_xml,
      &vtk_xml_hypertreegrid->vtk_hypertreegrid->y, format, appended);
  xml_tag_add_child(tag_grid, tag_y);

  // <DataArray Name="ZCoordinates" ... / >
  XMLTag *tag_z = vtk_xml_add_data_array(
      &vtk_xml_hypertreegrid->vtk_xml,
      &vtk_xml_hypertreegrid->vtk_hypertreegrid->z, format, appended);
  xml_tag_add_child(tag_grid, tag_z);

  // <Trees> </Trees>
  XMLTag *tag_trees = xml_tag_init("Trees", false);
  xml_tag_add_child(tag_hypertreegrid, tag_trees);

  // <DataArray Name="Descriptors" ... / >
  XMLTag *tag_descriptors = vtk_xml_add_data_array(
      &vtk_xml_hypertreegrid->vtk_xml,
      &vtk_xml_hypertreegrid->vtk_hypertreegrid->descriptors, format, appended);
  xml_tag_add_child(tag_trees, tag_descriptors);

  // <DataArray Name="NumberOfVerticesPerDepth" ... />
  XMLTag *tag_number_of_vertices_per_depth = vtk_xml_add_data_array(
      &vtk_xml_hypertreegrid->vtk_xml,
      &vtk_xml_hypertreegrid->vtk_hypertreegrid->number_of_vertices_per_depth,
      format, appended);
  xml_tag_add_child(tag_trees, tag_number_of_vertices_per_depth);

  // <DataArray Name="TreeIds" ... />
  XMLTag *tag_tree_ids = vtk_xml_add_data_array(
      &vtk_xml_hypertreegrid->vtk_xml,
      &vtk_xml_hypertreegrid->vtk_hypertreegrid->tree_ids, format, appended);
  xml_tag_add_child(tag_trees, tag_tree_ids);

  // <DataArray Name="DepthPerTree" ... />
  XMLTag *tag_depth_per_tree = vtk_xml_add_data_array(
      &vtk_xml_hypertreegrid->vtk_xml,
      &vtk_xml_hypertreegrid->vtk_hypertreegrid->depth_per_tree, format,
      appended);
  xml_tag_add_child(tag_trees, tag_depth_per_tree);

  // <CellData></CellData>
  XMLTag *tag_celldata = xml_tag_init("CellData", false);
  xml_tag_add_child(tag_hypertreegrid, tag_celldata);

  for (scalar s in scalar_list) {
    vtkDataArray *vtk_data_array_scalar = vtk_hypertreegrid_add_scalar(
        vtk_xml_hypertreegrid->vtk_hypertreegrid, s, field_type);
    vtk_data_array_scalar->number_of_components = 0;
    vtk_data_array_scalar->number_of_tuples =
        vtk_xml_hypertreegrid->vtk_hypertreegrid->data->total_vertices;

    XMLTag *tag_scalar =
        vtk_xml_add_data_array(&vtk_xml_hypertreegrid->vtk_xml,
                               vtk_data_array_scalar, format, appended);
    xml_tag_add_child(tag_celldata, tag_scalar);
  }

  for (vector v in vector_list) {
    vtkDataArray *vtk_data_array_vector = vtk_hypertreegrid_add_vector(
        vtk_xml_hypertreegrid->vtk_hypertreegrid, v, field_type);
    vtk_data_array_vector->number_of_components = dimension;
    vtk_data_array_vector->number_of_tuples =
        vtk_xml_hypertreegrid->vtk_hypertreegrid->data->total_vertices;

    XMLTag *tag_vector =
        vtk_xml_add_data_array(&vtk_xml_hypertreegrid->vtk_xml,
                               vtk_data_array_vector, format, appended);
    xml_tag_add_child(tag_celldata, tag_vector);
  }

  // <FieldData></FieldData>
  XMLTag *tag_fielddata = xml_tag_init("FieldData", false);
  xml_tag_add_child(tag_hypertreegrid, tag_fielddata);

  if (time >= 0) {
//    printf("[vtk_xml_hypertreegrid_init]: t=%f\n", time);
    vtk_xml_hypertreegrid->vtk_hypertreegrid->data->time = time;
//     printf("%f %f \n", time,
//            vtk_xml_hypertreegrid->vtk_hypertreegrid->data->time);
    XMLTag *tag_timevalue = vtk_xml_add_data_array(
        &vtk_xml_hypertreegrid->vtk_xml,
        &vtk_xml_hypertreegrid->vtk_hypertreegrid->time_value, format,
        appended);
    xml_tag_add_child(tag_fielddata, tag_timevalue);
  }

  return vtk_xml_hypertreegrid;
}

void vtk_xml_hypertreegrid_to_file(vtkXMLHyperTreeGrid *vtk_xml_hypertreegrid,
#if _MPI
                                   MPI_File fp
#else
                                   FILE *fp
#endif
) {
  vtk_xml_to_file(&vtk_xml_hypertreegrid->vtk_xml, fp);
}
