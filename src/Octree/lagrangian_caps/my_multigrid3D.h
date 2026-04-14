/*We add the cache structures to the original multigrid solver*/

#include "grid/multigrid3D.h"

#define MULT_GRID 1

/* Caches structures */

typedef struct {
  int i;
#if dimension >= 2
  int j;
#endif
#if dimension >= 3
  int k;
#endif  
} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;
#if dimension >= 2
  int j;
#endif
#if dimension >= 3
  int k;
#endif  
  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;

#define BSIZE 128

// static void cache_level_append (CacheLevel * c, Point p)
// {
//   if (c->n >= c->nm) {
//     c->nm += BSIZE;
//     qrealloc (c->p, c->nm, IndexLevel);
//   }
//   c->p[c->n].i = p.i;
// #if dimension >= 2
//   c->p[c->n].j = p.j;
// #endif
// #if dimension >= 3
//   c->p[c->n].k = p.k;
// #endif
//   c->n++;
// }

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/BSIZE + 1)*BSIZE) {
    c->nm = (c->n/BSIZE + 1)*BSIZE;
    assert (c->nm > c->n);
    c->p = (IndexLevel *) realloc (c->p, sizeof (Index)*c->nm);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += BSIZE;
    qrealloc (c->p, c->nm, Index);
  }
  c->p[c->n].i = p.i;
#if dimension >= 2
  c->p[c->n].j = p.j;
#endif
#if dimension >= 3
  c->p[c->n].k = p.k;
#endif  
  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}

#undef BSIZE


#define update_cache() { if (tree->dirty) update_cache_f(); }

macro2 foreach_cache (Cache cache, Reduce reductions = None)
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    point.i = GHOSTS;
#if dimension > 1
    point.j = GHOSTS;
#endif
#if dimension > 2
    point.k = GHOSTS;
#endif
    int _k; unsigned short _flags; NOT_UNUSED(_flags);
    OMP(omp for schedule(static))
      for (_k = 0; _k < cache.n; _k++) {
	point.i = cache.p[_k].i;
#if dimension >= 2
	point.j = cache.p[_k].j;
#endif
#if dimension >= 3
	point.k = cache.p[_k].k;
#endif
	point.level = cache.p[_k].level;
	_flags = cache.p[_k].flags;
	{...}
      }
  }
}

macro2 foreach_cache_level (Cache cache, int _l, Reduce reductions = None)
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0}; NOT_UNUSED (point);
    point.i = GHOSTS;
#if dimension > 1
    point.j = GHOSTS;
#endif
#if dimension > 2
    point.k = GHOSTS;
#endif
    point.level = _l;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 0; _k < cache.n; _k++) {
	point.i = cache.p[_k].i;
#if dimension >= 2
	point.j = cache.p[_k].j;
#endif
#if dimension >= 3
	point.k = cache.p[_k].k;
#endif
	{...}
      }
  }
}


// static void free_cache (CacheLevel * c)
// {
//   for (int l = 0; l <= depth(); l++)
//     free (c[l].p);
//   free (c);
// }


