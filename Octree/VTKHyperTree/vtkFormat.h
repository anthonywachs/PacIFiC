typedef enum {
  ascii = 0,
  raw = 1,
  base64 = 2,
} vtkFormat;

char *vtkFormat_to_string(vtkFormat ft) {
  switch (ft) {
  case 0:
    return "ascii";
  case 1:
    return "raw";
  case 2:
    return "base64";
  default:
    return "None";
  }
}
