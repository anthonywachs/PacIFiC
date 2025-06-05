#include "XMLDocument.h"
#include "vtkByteOrder.h"
#include "vtkDataArray.h"
#include "vtkFileType.h"
#include "vtkFormat.h"
#include "vtkType.h"

typedef struct {
  vtkFileType file_type;
  vtkByteOrder byte_order;
  vtkType header_type;
  vtkDataArrayList inline_da;
  vtkDataArrayList appended_da;
  vtkFormat appended_fmt;
  size_t appended_offset;
  XMLDocument doc;
} vtkXML;

vtkXML vtk_xml_init(vtkFileType file_type, vtkType header_type) {
  vtkXML vtk_xml = {.file_type = file_type,
                    .byte_order = BYTE_ORDER,
                    .header_type = header_type,
                    .inline_da = vtk_data_array_list_init(),
                    .appended_da = vtk_data_array_list_init(),
                    .doc = doc_init()};

  vtk_xml.doc.root = xml_tag_init("VTKFile", false);

  xml_tag_add_attribute(vtk_xml.doc.root, "type",
                        vtkFileType_to_string(file_type));
  xml_tag_add_attribute(vtk_xml.doc.root, "byte_order",
                        vtkByteOrder_to_string(BYTE_ORDER));
  xml_tag_add_attribute(vtk_xml.doc.root, "header_type",
                        vtkType_to_string(header_type));

  return vtk_xml;
}

void vtk_xml_free(vtkXML *vtk_xml) {
  vtk_data_array_list_clear(&vtk_xml->inline_da);
  vtk_data_array_list_clear(&vtk_xml->appended_da);
  doc_free(&vtk_xml->doc);
  // free(vtk_xml);
}

size_t _vtk_data_array_compute_header(vtkDataArray *vtk_data_array) {
  size_t number_of_tuples = vtk_data_array->number_of_tuples;
  size_t number_of_components = vtk_data_array->number_of_components;
  number_of_components = number_of_components > 1 ? number_of_components : 1;
  size_t number_of_bytes = number_of_components * number_of_tuples;

  if (vtk_data_array->type == 10) {
    number_of_bytes = (number_of_bytes + 7) / 8;
  } else {
    number_of_bytes *= vtkType_sizeof(vtk_data_array->type);
  }

  return number_of_bytes;
}

/**
 *
 */

typedef struct {
  vtkType header_type;
  vtkFormat format;
  vtkDataArray *vtk_data_array;
} vtkDataArrayWriterCtx;

static void _xmlcontentwriter_to_vtkdataarraywriter_thunk(void *v,
#if _MPI
                                                          MPI_Offset base,
                                                          MPI_File fp
#else
                                                          FILE *fp
#endif
) {
  vtkDataArrayWriterCtx *c = v;
#if _MPI
  c->vtk_data_array->writer(c->vtk_data_array->writer_ctx, c->format, base, fp);
#else
  if (c->format == 1) {
    fputc('_', fp);
    size_t header = _vtk_data_array_compute_header(c->vtk_data_array);

    switch (c->header_type) {
    case 2: {
      UInt32_t header_cast = (UInt32_t)header;
      fwrite(&header_cast, sizeof(UInt32_t), 1, fp);
      break;
    }
    case 3: {
      UInt64_t header_cast = (UInt64_t)header;
      fwrite(&header_cast, sizeof(UInt64_t), 1, fp);
      break;
    }
    default: {
      UInt32_t header_cast = (UInt32_t)header;
      fwrite(&header_cast, sizeof(UInt32_t), 1, fp);
      break;
    }
    }
  }
  c->vtk_data_array->writer(c->vtk_data_array->writer_ctx, c->format, fp);
#endif
}

/**
 */

XMLTag *_vtk_data_array_to_tag(vtkDataArray *vtk_data_array,
                               vtkType header_type, bool appended) {
  XMLTag *tag = xml_tag_init("DataArray", appended);

  xml_tag_add_attribute(tag, "type", vtkType_to_string(vtk_data_array->type));
  xml_tag_add_attribute(tag, "Name", vtk_data_array->name);

  char buffer_val[32];
  sprintf(buffer_val, "%lu", vtk_data_array->number_of_tuples);
  xml_tag_add_attribute(tag, "NumberOfTuples", buffer_val);

  vtkDataArrayWriterCtx *c = malloc(sizeof *c);
  c->vtk_data_array = vtk_data_array;
  c->header_type = header_type;

  if (appended) {
    xml_tag_add_attribute(tag, "format", "appended");
    c->format = raw; // raw
  } else {
    xml_tag_add_attribute(tag, "format", "ascii");
    c->format = ascii; // ascii
  }

  tag->content_writer = _xmlcontentwriter_to_vtkdataarraywriter_thunk;
  tag->content_writer_ctx = c;

  return tag;
}

/*

*/

typedef struct {
  vtkDataArrayList *list;
  vtkType header_type;
} vtkXMLAppendedCtx;

static void _vtk_xml_appended_content_writer(void *v,
#if _MPI
                                             MPI_Offset base, MPI_File fp
#else
                                             FILE *fp
#endif
) {
  vtkXMLAppendedCtx *c = v;

#if _MPI
  size_t *data_offsets = malloc(c->list->count * sizeof(size_t));
  size_t *header_offsets = malloc(c->list->count * sizeof(size_t));
  size_t *header = malloc(c->list->count * sizeof(size_t));
  size_t running_offset = ((size_t)base) + 1;

  for (size_t i = 0; i < c->list->count; i++) {
    header_offsets[i] = running_offset;
    running_offset += vtkType_sizeof(c->header_type);
    header[i] = _vtk_data_array_compute_header(c->list->vtk_data_arrays[i]);
    data_offsets[i] = running_offset;
    running_offset += header[i];
  }

  if (pid() == 0) {
    MPI_File_write_at(fp, base, "_", strlen("_"), MPI_BYTE, MPI_STATUS_IGNORE);
  }

  for (size_t i = 0; i < c->list->count; i++) {
    vtkDataArray *vtk_data_array = c->list->vtk_data_arrays[i];

    switch (c->header_type) {
    case 2: { // UInt32_e
      UInt32_t header_cast = (UInt32_t)header[i];
      MPI_File_write_at(fp, header_offsets[i], &header_cast, sizeof(UInt32_t),
                        MPI_BYTE, MPI_STATUS_IGNORE);
      break;
    }
    case 3: { // UInt64_e
      UInt64_t header_cast = (UInt64_t)header[i];
      MPI_File_write_at(fp, header_offsets[i], &header_cast, sizeof(UInt64_t),
                        MPI_BYTE, MPI_STATUS_IGNORE);
      break;
    }
    default: {
      UInt32_t header_cast = (UInt32_t)header[i];
      MPI_File_write_at(fp, header_offsets[i], &header_cast, sizeof(UInt32_t),
                        MPI_BYTE, MPI_STATUS_IGNORE);
      break;
    }
    }

    vtk_data_array->writer(vtk_data_array->writer_ctx, 1, data_offsets[i],
                           fp);
  }

  free(data_offsets);
  free(header_offsets);
  free(header);

#else
  fputc('_', fp);

  for (size_t i = 0; i < c->list->count; i++) {
    vtkDataArray *vtk_data_array = c->list->vtk_data_arrays[i];
    size_t header = _vtk_data_array_compute_header(vtk_data_array);

    switch (c->header_type) {
    case 2: {
      UInt32_t header_cast = (UInt32_t)header;
      fwrite(&header_cast, sizeof(UInt32_t), 1, fp);
      break;
    }
    case 3: {
      UInt64_t header_cast = (UInt64_t)header;
      fwrite(&header_cast, sizeof(UInt64_t), 1, fp);
      break;
    }
    default: {
      UInt32_t header_cast = (UInt32_t)header;
      fwrite(&header_cast, sizeof(UInt32_t), 1, fp);
      break;
    }
    }

    vtk_data_array->writer(vtk_data_array->writer_ctx, 1, fp);
  }
#endif
}

XMLTag *vtk_xml_get_appended_tag(vtkXML *vtk_xml) {
  XMLTag *xml_tag = xml_tag_init("AppendedData", false);

  xml_tag_add_attribute(xml_tag, "encoding", "raw");

  vtkXMLAppendedCtx *c = malloc(sizeof *c);
  c->header_type = vtk_xml->header_type;
  c->list = &vtk_xml->appended_da;

  xml_tag->content_writer_ctx = c;
  xml_tag->content_writer = _vtk_xml_appended_content_writer;

  return xml_tag;
}

XMLTag *vtk_xml_add_data_array(vtkXML *vtk_xml, vtkDataArray *vtk_data_array,
                               vtkFormat format, bool appended) {

  vtkFormat actual_format;
  bool actual_appended;

  if (appended) {
    actual_appended = appended;
    if (format == 1) {
      actual_format = 1;
    } else {
      actual_format = format;
    }
  } else {
#if _MPI
    actual_format = raw;
    actual_appended = true;
#else
    actual_format = ascii;
    actual_appended = false;
#endif
  }

  // Create a new XMLTag
  XMLTag *tag = xml_tag_init("DataArray", appended);

  xml_tag_add_attribute(tag, "type", vtkType_to_string(vtk_data_array->type));
  xml_tag_add_attribute(tag, "Name", vtk_data_array->name);

  char buffer_val[32];
  sprintf(buffer_val, "%lu", vtk_data_array->number_of_tuples);
  xml_tag_add_attribute(tag, "NumberOfTuples", buffer_val);

  if (vtk_data_array->number_of_components != 0) {
    sprintf(buffer_val, "%lu", vtk_data_array->number_of_components);
    xml_tag_add_attribute(tag, "NumberOfComponents", buffer_val);
  }

  if (actual_appended) {
    // Set the format
    xml_tag_add_attribute(tag, "format", "appended");

    // Add the <Appended></Appended> tags if not already at the root
    if (vtk_xml->appended_da.count == 0) {
      XMLTag *xml_tag_appended = vtk_xml_get_appended_tag(vtk_xml);
      xml_tag_add_child(vtk_xml->doc.root, xml_tag_appended);
    }

    // Store our vtkDataArray
    vtk_data_array_list_add(&vtk_xml->appended_da, vtk_data_array, false);

    // Set the offset
    sprintf(buffer_val, "%lu", vtk_xml->appended_offset);
    xml_tag_add_attribute(tag, "offset", buffer_val);

    // Compute the offset for the next appended data array
    vtk_xml->appended_offset += vtkType_sizeof(vtk_xml->header_type) +
                                _vtk_data_array_compute_header(vtk_data_array);

    return tag;
  } else {
    // Set the format
    switch (actual_format) {
    case 0: // ascii
      xml_tag_add_attribute(tag, "format", "ascii");
      break;
    case 1: // raw
      xml_tag_add_attribute(tag, "format", "raw");
      break;
    default:
      assert(false);
      break;
    }

    // Add the vtkDataArrayWriter as an XMLContentWriter to the tag
    vtkDataArrayWriterCtx *c = malloc(sizeof *c);
    c->format = actual_format;
    c->vtk_data_array = vtk_data_array;
    c->header_type = vtk_xml->header_type;

    tag->content_writer = _xmlcontentwriter_to_vtkdataarraywriter_thunk;
    tag->content_writer_ctx = c;

    // Store our vtkDataArray
    vtk_data_array_list_add(&vtk_xml->inline_da, vtk_data_array, false);
  }

  return tag;
}

void vtk_xml_to_file(vtkXML *vtk_xml,
#if _MPI
                     MPI_File fp
#else
                     FILE *fp
#endif
) {
  doc_to_file(&vtk_xml->doc, fp);
}
