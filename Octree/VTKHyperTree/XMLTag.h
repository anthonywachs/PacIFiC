#define CHILDREN_CHUNK_SIZE 16
#define INDENT_STRING "  "

/* a single key/value pair */

typedef void (*XMLContentWriter)(void *ctx,
#if _MPI
                                 MPI_Offset base, MPI_File fp
#else
                                 FILE *fp
#endif
);

#include "XMLAttribute.h"

typedef struct XMLTag XMLTag;

struct XMLTag {
  // Tag <Name>
  char *name;

  // Tag <Name attribute="value">
  XMLAttributeList attribute_list;
  XMLAttributeListFactory attribute_list_factory;

  // Inner-tag <Name>content</Name>, either a function pointer or a string
  char *content;

  // Children
  XMLTag **children;
  size_t n_children;
  size_t cap_children;

  XMLContentWriter content_writer;
  void *content_writer_ctx;

  // If true, <Name/> instead of <Name></Name>
  bool self_closing;
};

/* allocate & initialize a tag */
XMLTag *xml_tag_init(const char *name, bool self_closing) {
  XMLTag *tag = calloc(1, sizeof(XMLTag));
  if (!tag)
    return NULL;

  tag->name = name ? strdup(name) : NULL;
  tag->self_closing = self_closing;
  tag->n_children = 0;
  tag->cap_children = self_closing ? 0 : CHILDREN_CHUNK_SIZE;
  tag->children =
      self_closing ? NULL : malloc(tag->cap_children * sizeof *tag->children);
  tag->attribute_list = xml_attribute_list_init();
  tag->attribute_list_factory = NULL;
  tag->content = NULL;
  tag->content_writer = NULL;
  tag->content_writer_ctx = NULL;

  return tag;
}

/* free a tag and its entire subtree */
void xml_tag_free(XMLTag *tag) {
  if (!tag)
    return;
  for (size_t i = 0; i < tag->n_children; i++) {
    xml_tag_free(tag->children[i]);
  }
  free(tag->children);
  free(tag->content);
  free(tag->name);
  xml_attribute_list_clear(&tag->attribute_list);
}

/* append a child, growing in CHUNK-sized steps */
void xml_tag_add_child(XMLTag *tag, XMLTag *child) {
  if (tag->self_closing)
    return;
  if (tag->n_children == tag->cap_children) {
    size_t new_cap = tag->cap_children + CHILDREN_CHUNK_SIZE;
    XMLTag **tmp = realloc(tag->children, new_cap * sizeof *tmp);
    if (!tmp) {
      perror("realloc children");
      exit(EXIT_FAILURE);
    }
    tag->children = tmp;
    tag->cap_children = new_cap;
  }
  tag->children[tag->n_children++] = child;
}

/* append an attribute, growing in CHUNK-sized steps */
void xml_tag_add_attribute(XMLTag *tag, const char *key, const char *value) {
  xml_attribute_list_add(&tag->attribute_list, key, value);
}

void xml_tag_set_attribute_list_factory(
    XMLTag *tag, XMLAttributeListFactory attribute_list_factory) {
  tag->attribute_list_factory = attribute_list_factory;
}

void xml_tag_to_file(XMLTag *tag, size_t indent,
#if _MPI
                     MPI_File fp
#else
                     FILE *fp
#endif
) {
#if _MPI
  MPI_Offset pos_shared = 0;
  MPI_Offset pos = 0;
  if (pid() == 0) {
    // Indent the opening tag
    for (size_t i = 0; i < indent; i++) {
      MPI_File_write_shared(fp, INDENT_STRING, strlen(INDENT_STRING), MPI_BYTE,
                            MPI_STATUS_IGNORE);
    }

    // Open the beginning tag
    MPI_File_write_shared(fp, "<", strlen("<"), MPI_BYTE, MPI_STATUS_IGNORE);
    MPI_File_write_shared(fp, tag->name, strlen(tag->name), MPI_BYTE,
                          MPI_STATUS_IGNORE);

    // Write attributes for the tag
    xml_attribute_list_to_file(&tag->attribute_list, fp);

    // Write attribute list factory attributes for the tag
    if (tag->attribute_list_factory != NULL) {
      XMLAttributeList attribute_list = tag->attribute_list_factory();
      xml_attribute_list_to_file(&attribute_list, fp);
      xml_attribute_list_clear(&attribute_list);
    }
  }

  if (tag->self_closing) {
    if (pid() == 0) {
      MPI_File_write_shared(fp, "/>\n", strlen("/>\n"), MPI_BYTE,
                            MPI_STATUS_IGNORE);
    }
  } else {
    if (pid() == 0) {
      MPI_File_write_shared(fp, ">\n", strlen(">\n"), MPI_BYTE,
                            MPI_STATUS_IGNORE);
    }

    for (size_t i = 0; i < tag->n_children; i++)
      xml_tag_to_file(tag->children[i], indent + 1, fp);

    // char *content, *content_writer;

    if (tag->content) {
      //
      if (pid() == 0) {
        for (size_t i = 0; i < indent + 1; i++) {
          MPI_File_write_shared(fp, INDENT_STRING, strlen(INDENT_STRING),
                                MPI_BYTE, MPI_STATUS_IGNORE);
        }
        // capture the shared pointer’s new position
        MPI_File_get_position_shared(fp, &pos_shared);
        MPI_File_get_position(fp, &pos);
        pos = pos > pos_shared ? pos : pos_shared;
      }

      // MPI_Barrier(MPI_COMM_WORLD);
      // MPI_Bcast(&pos, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);

      // MPI_File_seek(fp,pos,MPI_SEEK_SET);
      // MPI_File_seek_shared(fp,pos,MPI_SEEK_SET);

      // MPI_File_write_shared(fp,tag->content,strlen(tag->content),MPI_BYTE,MPI_STATUS_IGNORE);
      // MPI_Barrier(MPI_COMM_WORLD);

      if (pid() == 0) {
        MPI_File_write_shared(fp, tag->content, strlen(tag->content), MPI_BYTE,
                              MPI_STATUS_IGNORE);
        MPI_File_write_shared(fp, "\n", strlen("\n"), MPI_BYTE,
                              MPI_STATUS_IGNORE);
      }
    } else if (tag->content_writer) {
      if (pid() == 0) {
        for (size_t i = 0; i < indent + 1; i++) {
          // use shared‐pointer write so the shared file pointer advances
          MPI_File_write_shared(fp, INDENT_STRING, strlen(INDENT_STRING),
                                MPI_BYTE, MPI_STATUS_IGNORE);
        }
        // capture the shared pointer’s new position
        MPI_File_get_position_shared(fp, &pos_shared);
        MPI_File_get_position(fp, &pos);
        pos = pos > pos_shared ? pos : pos_shared;
      }

      // broadcast that base offset to everyone
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&pos, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);

      // seek both pointers to pos and flush
      MPI_File_seek(fp, pos, MPI_SEEK_SET);
      MPI_File_seek_shared(fp, pos, MPI_SEEK_SET);

      tag->content_writer(tag->content_writer_ctx, pos, fp);

      // MPI_File_get_position_shared(fp, &pos_shared);
      // MPI_File_get_position(fp, &pos);
      // pos = pos > pos_shared ? pos : pos_shared;

      // if (pid() == 0)
      //     MPI_Reduce(MPI_IN_PLACE,&pos,1,MPI_OFFSET,MPI_MAX,0,MPI_COMM_WORLD);
      // else
      //     MPI_Reduce(&pos,NULL,1,MPI_OFFSET,MPI_MAX,0,MPI_COMM_WORLD);

      MPI_File_sync(fp);
      MPI_File_get_size(fp, &pos);

      MPI_File_seek(fp, pos, MPI_SEEK_SET);
      MPI_File_seek_shared(fp, pos, MPI_SEEK_SET);

      // advance rank 0’s pointers past all appended data
      if (pid() == 0) {
        MPI_File_write_shared(fp, "\n", 1, MPI_BYTE, MPI_STATUS_IGNORE);
      }
    }

    if (pid() == 0) {
      for (size_t i = 0; i < indent; i++) {
        MPI_File_write_shared(fp, INDENT_STRING, strlen(INDENT_STRING),
                              MPI_BYTE, MPI_STATUS_IGNORE);
      }
      MPI_File_write_shared(fp, "</", strlen("</"), MPI_BYTE,
                            MPI_STATUS_IGNORE);
      MPI_File_write_shared(fp, tag->name, strlen(tag->name), MPI_BYTE,
                            MPI_STATUS_IGNORE);
      MPI_File_write_shared(fp, ">\n", strlen(">\n"), MPI_BYTE,
                            MPI_STATUS_IGNORE);
    }
  }

#else
  for (size_t i = 0; i < indent; i++)
    fputs(INDENT_STRING, fp);

  fputc('<', fp);
  fputs(tag->name, fp);

  xml_attribute_list_to_file(&tag->attribute_list, fp);
  if (tag->attribute_list_factory != NULL) {
    XMLAttributeList attribute_list = tag->attribute_list_factory();
    xml_attribute_list_to_file(&attribute_list, fp);
    xml_attribute_list_clear(&attribute_list);
  }

  if (tag->self_closing) {
    fputs("/>\n", fp);
    return;
  } else {
    fputs(">\n", fp);
    for (size_t i = 0; i < tag->n_children; i++) {
      xml_tag_to_file(tag->children[i], indent + 1, fp);
    }
    if (tag->content) {
      for (size_t i = 0; i < indent + 1; i++)
        fputs(INDENT_STRING, fp);
      fputs(tag->content, fp);
      fputc('\n', fp);
    } else if (tag->content_writer) {
      for (size_t i = 0; i < indent + 1; i++)
        fputs(INDENT_STRING, fp);
      tag->content_writer(tag->content_writer_ctx, fp);
      fputc('\n', fp);
    }
    for (size_t i = 0; i < indent; i++)
      fputs(INDENT_STRING, fp);
    fprintf(fp, "</%s>\n", tag->name);
  }
#endif
}