typedef uint8_t UInt8_t;
typedef uint16_t UInt16_t;
typedef uint32_t UInt32_t;
typedef uint64_t UInt64_t;
typedef int8_t Int8_t;
typedef int16_t Int16_t;
typedef int32_t Int32_t;
typedef int64_t Int64_t;
typedef float Float32_t;
typedef double Float64_t;
typedef UInt8_t Bit_t;

typedef enum {
  UInt8_e = 0,
  UInt16_e = 1,
  UInt32_e = 2,
  UInt64_e = 3,
  Int8_e = 4,
  Int16_e = 5,
  Int32_e = 6,
  Int64_e = 7,
  Float32_e = 8,
  Float64_e = 9,
  Bit_e = 10,
} vtkType;

// typedef enum vtkType vtkType;

size_t vtkType_sizeof(vtkType t) {
  switch (t) {
  case 0: // UInt8
    return sizeof(UInt8_t);
    break;
  case 1: // UInt16
    return sizeof(UInt16_t);
    break;
  case 2: // UInt32
    return sizeof(UInt32_t);
    break;
  case 3: // UInt64
    return sizeof(UInt64_t);
    break;
  case 4: // Int8
    return sizeof(Int8_t);
    break;
  case 5: // Int16
    return sizeof(Int16_t);
    break;
  case 6: // Int32
    return sizeof(Int32_t);
    break;
  case 7: // Int64
    return sizeof(Int64_t);
    break;
  case 8: // Float32
    return sizeof(Float32_t);
    break;
  case 9: // Float64
    return sizeof(Float64_t);
    break;
  case 10: // Bit
    return sizeof(Bit_t);
    break;
  default:
    return 0;
    break;
  }
}

char *vtkType_to_string(vtkType t) {
  switch (t) {
  case 0:
    return "UInt8";
    break;
  case 1:
    return "UInt16";
    break;
  case 2:
    return "UInt32";
    break;
  case 3:
    return "UInt64";
    break;
  case 4:
    return "Int8";
    break;
  case 5:
    return "Int16";
    break;
  case 6:
    return "Int32";
    break;
  case 7:
    return "Int64";
    break;
  case 8:
    return "Float32";
    break;
  case 9:
    return "Float64";
    break;
  case 10:
    return "Bit";
    break;
  default:
    return "None";
    break;
  }
}

// #define VTK_VOID 0
// #define VTK_BIT 1
// #define VTK_CHAR 2
// #define VTK_SIGNED_CHAR 15
// #define VTK_UNSIGNED_CHAR 3
// #define VTK_SHORT 4
// #define VTK_UNSIGNED_SHORT 5
// #define VTK_INT 6
// #define VTK_UNSIGNED_INT 7
// #define VTK_LONG 8
// #define VTK_UNSIGNED_LONG 9
// #define VTK_FLOAT 10
// #define VTK_DOUBLE 11
// #define VTK_ID_TYPE 12
// #define VTK_STRING 13
// #define VTK_OPAQUE 14
// #define VTK_LONG_LONG 16
// #define VTK_UNSIGNED_LONG_LONG 17
// #define VTK_VARIANT 20
// #define VTK_OBJECT 21