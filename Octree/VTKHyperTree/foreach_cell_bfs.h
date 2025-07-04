macro2 foreach_cell_root_BFS(Point root)
{
  {

#if dimension == 1
    Point point = {0};
    int ig = 0; 
    NOT_UNUSED(ig);
    typedef struct { int l, i; } queue_t;
#elif dimension == 2
    Point point = {0};
    int ig = 0, jg = 0; 
    NOT_UNUSED(ig); NOT_UNUSED(jg);
    typedef struct { int l, i, j; } queue_t;
#else // dimension == 3
    Point point = {0};
    int ig = 0, jg = 0;	int kg = 0; 
    NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    typedef struct { int l, i, j, k; } queue_t;
#endif

    size_t QMAX = 0;
    foreach_cell() { QMAX++; }

    queue_t queue[QMAX];
    size_t front = 0, rear = 0;

#if dimension == 1
        queue[rear++] = (queue_t){root.level, root.i};
#elif dimension == 2
        queue[rear++] = (queue_t){root.level, root.i, root.j};        
#else // dimension == 3
        queue[rear++] = (queue_t){root.level, root.i, root.j, root.k};
#endif

    while (front < rear) {
      queue_t current = queue[front++];
      point.level = current.l;
#if dimension == 1                
          point.i = current.i;
#elif dimension == 2
          point.i = current.i;
          point.j = current.j;
#else // dimension == 3
          point.i = current.i;
          point.j = current.j;
          point.k = current.k;
#endif

      /* Check if the cell is allocated (macro 'allocated' expects 'point') */
      if (!allocated(0,0,0))
          continue;

      {...} // User defined code block

#if dimension == 1
      if (point.level < grid->depth && !is_leaf(cell)) {
          if (rear + 2 > QMAX) abort();
          queue[rear++] = (queue_t){point.level+1, _LEFT };
          queue[rear++] = (queue_t){point.level+1, _RIGHT};
      }
#elif dimension == 2
      if (point.level < grid->depth && !is_leaf(cell)) {
          if (rear + 4 > QMAX) abort();
          queue[rear++] = (queue_t){point.level+1,_LEFT ,_BOTTOM};
          queue[rear++] = (queue_t){point.level+1,_RIGHT,_BOTTOM};
          queue[rear++] = (queue_t){point.level+1,_LEFT ,_TOP   };
          queue[rear++] = (queue_t){point.level+1,_RIGHT,_TOP   };
      }
#else // dimension == 3
      if (point.level < grid->depth && !is_leaf(cell)) {
          if (rear + 8 > QMAX) abort();
          queue[rear++] = (queue_t){point.level+1,_LEFT ,_BOTTOM,_BACK };
          queue[rear++] = (queue_t){point.level+1,_RIGHT,_BOTTOM,_BACK };
          queue[rear++] = (queue_t){point.level+1,_LEFT ,_TOP   ,_BACK };
          queue[rear++] = (queue_t){point.level+1,_RIGHT,_TOP   ,_BACK };
          queue[rear++] = (queue_t){point.level+1,_LEFT ,_BOTTOM,_FRONT};
          queue[rear++] = (queue_t){point.level+1,_RIGHT,_BOTTOM,_FRONT};
          queue[rear++] = (queue_t){point.level+1,_LEFT ,_TOP   ,_FRONT};
          queue[rear++] = (queue_t){point.level+1,_RIGHT,_TOP   ,_FRONT};
      }
#endif
    } // End of loop over the current level
  }
}

macro2 foreach_cell_BFS()
{
  {
#if dimension == 1
      Point root = {GHOSTS,0};
#elif dimension == 2
      Point root = {GHOSTS,GHOSTS,0};
#else // dimension == 3
      Point root = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
    foreach_cell_root_BFS (root)
      {...}
  }
}
