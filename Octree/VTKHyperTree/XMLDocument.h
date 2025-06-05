#include "XMLTag.h"

typedef struct {
  char *version;
  XMLTag *root;
} XMLDocument;

XMLDocument doc_init(void) {
  XMLDocument doc = {
      .version = strdup("1.0"),
      .root = NULL,
  };

  return doc;
}

void doc_free(XMLDocument *doc) {
  if (!doc)
    return;
  free(doc->version);
  if (doc->root) {
    xml_tag_free(doc->root);
  }
}

void doc_to_file(XMLDocument *doc,
#if _MPI
                 MPI_File fp
#else
                 FILE *fp
#endif
) {
#if _MPI
  if (doc->version != NULL) {
    char buffer[64];
    if (pid() == 0) {
      size_t len = sprintf(buffer, "<?xml version=\"%s\"?>\n", doc->version);
      MPI_File_write_shared(fp, buffer, len, MPI_BYTE, MPI_STATUS_IGNORE);
    }
  }
#else
  if (doc->version != NULL) {
    fprintf(fp, "<?xml version=\"%s\"?>\n", doc->version);
  }
#endif
  xml_tag_to_file(doc->root, 0, fp);
}