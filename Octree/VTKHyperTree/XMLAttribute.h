#define ATTRIBUTES_CHUNK_SIZE 8

typedef struct {
  char *key;
  char *value;
} XMLAttribute;

XMLAttribute *xml_attribute_init(const char *key, const char *value) {
  XMLAttribute *attribute = malloc(sizeof(XMLAttribute));
  if (!attribute) {
    perror("malloc attribute");
    exit(EXIT_FAILURE);
  }
  attribute->key = strdup(key);
  attribute->value = strdup(value);
  return attribute;
}

void xml_attribute_free(XMLAttribute *attribute) {
  free(attribute->key);
  free(attribute->value);
  free(attribute);
}

typedef struct {
  XMLAttribute **attributes;
  size_t count;
  size_t capacity;
} XMLAttributeList;

XMLAttributeList xml_attribute_list_init() {
  XMLAttributeList list = {.attributes = NULL, .count = 0, .capacity = 0};

  return list;
}

void xml_attribute_list_add(XMLAttributeList *list, const char *key,
                            const char *value) {
  if (list->count == list->capacity) {
    size_t new_capacity = list->capacity + ATTRIBUTES_CHUNK_SIZE;
    XMLAttribute **tmp = realloc(list->attributes, new_capacity * sizeof *tmp);
    if (!tmp) {
      perror("realloc attributes");
      exit(EXIT_FAILURE);
    }
    list->attributes = tmp;
    list->capacity = new_capacity;
  }

  XMLAttribute *attribute = xml_attribute_init(key, value);
  list->attributes[list->count++] = attribute;
}

void xml_attribute_list_clear(XMLAttributeList *list) {
  for (size_t i = 0; i < list->count; i++) {
    xml_attribute_free(list->attributes[i]);
    list->attributes[i] = NULL;
  }

  free(list->attributes);

  list->attributes = NULL;
  list->capacity = 0;
  list->count = 0;
}

void xml_attribute_list_to_file(XMLAttributeList *xml_attribute_list,
#if _MPI
                                MPI_File fp
#else
                                FILE *fp
#endif
) {
#if _MPI
  if (pid() == 0) {
    char buffer[256];
    for (size_t i = 0; i < xml_attribute_list->count; i++) {
      XMLAttribute *attribute = xml_attribute_list->attributes[i];
      size_t len =
          sprintf(buffer, " %s=\"%s\"", attribute->key, attribute->value);
      MPI_File_write_shared(fp, buffer, len, MPI_BYTE, MPI_STATUS_IGNORE);
    }
  }
#else
  for (size_t i = 0; i < xml_attribute_list->count; i++) {
    XMLAttribute *attribute = xml_attribute_list->attributes[i];
    fprintf(fp, " %s=\"%s\"", attribute->key, attribute->value);
  }
#endif
}

typedef XMLAttributeList (*XMLAttributeListFactory)();
