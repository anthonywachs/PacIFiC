typedef enum {
  HyperTreeGrid = 0,
  PHyperTreeGrid = 1,
} vtkFileType;

char *vtkFileType_to_string(vtkFileType ft) {
  switch (ft) {
  case 0:
    return "HyperTreeGrid";
  case 1:
    return "PHyperTreeGrid";
  default:
    return "None";
  }
}
