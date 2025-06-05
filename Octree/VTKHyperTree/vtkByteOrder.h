typedef enum { BigEndian = 0, LittleEndian = 1 } vtkByteOrder;

#if defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) &&                \
    (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
// BigEndian
#define BYTE_ORDER 0
#else
// LittleEndian
#define BYTE_ORDER 1
#endif

char *vtkByteOrder_to_string(vtkByteOrder bo) {
  switch (bo) {
  case 0:
    return "BigEndian";
  case 1:
    return "LittleEndian";
  default:
    return "None";
  }
}
