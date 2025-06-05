static inline uint64_t point_to_morton(Point p) {
  uint64_t code = 0;
  int L = p.level;

#if dimension == 2
  // remove ghost‐cell bias
  uint32_t i = p.i - GHOSTS;
  uint32_t j = p.j - GHOSTS;

  // interleave X then Y → SW,SE,NW,NE
  for (int b = 0; b < L; b++) {
    code |= (uint64_t)((i >> b) & 1) << (2 * b);     // X bit
    code |= (uint64_t)((j >> b) & 1) << (2 * b + 1); // Y bit
  }
  // sentinel at bit 2*L
  code |= 1ULL << (2 * L);

#elif dimension == 3
  uint32_t i = p.i - GHOSTS;
  uint32_t j = p.j - GHOSTS;
  uint32_t k = p.k - GHOSTS;
  for (int b = 0; b < L; b++) {
    code |= (uint64_t)((i >> b) & 1) << (3 * b);
    code |= (uint64_t)((j >> b) & 1) << (3 * b + 1);
    code |= (uint64_t)((k >> b) & 1) << (3 * b + 2);
  }
  code |= 1ULL << (3 * L);

#else // dimension == 1
  uint32_t i = p.i - GHOSTS;
  for (int b = 0; b < L; b++) {
    code |= (uint64_t)((i >> b) & 1) << b;
  }
  code |= 1ULL << L;
#endif

  return code;
}

static inline Point morton_to_point(uint64_t code) {
  // locate the sentinel
  int h = 63 - __builtin_clzll(code);
  Point p = {0};

#if dimension == 2
  p.level = h >> 1; // sentinel at 2*level
  uint64_t base = code ^ (1ULL << h);
  uint32_t i = 0, j = 0;
  for (int b = 0; b < p.level; b++) {
    i |= (uint32_t)((base >> (2 * b)) & 1) << b;
    j |= (uint32_t)((base >> (2 * b + 1)) & 1) << b;
  }
  // restore ghost offset
  p.i = i + GHOSTS;
  p.j = j + GHOSTS;

#elif dimension == 3
  p.level = h / 3; // sentinel at 3*level
  uint64_t base = code ^ (1ULL << h);
  uint32_t i = 0, j = 0, k = 0;
  for (int b = 0; b < p.level; b++) {
    i |= (uint32_t)((base >> (3 * b)) & 1) << b;
    j |= (uint32_t)((base >> (3 * b + 1)) & 1) << b;
    k |= (uint32_t)((base >> (3 * b + 2)) & 1) << b;
  }
  p.i = i + GHOSTS;
  p.j = j + GHOSTS;
  p.k = k + GHOSTS;

#else // dimension == 1
  p.level = h; // sentinel at bit level
  uint64_t base = code ^ (1ULL << h);
  uint32_t i = 0;
  for (int b = 0; b < p.level; b++) {
    i |= (uint32_t)((base >> b) & 1) << b;
  }
  p.i = i + GHOSTS;
#endif

  return p;
}
