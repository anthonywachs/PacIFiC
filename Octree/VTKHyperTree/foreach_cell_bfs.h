macro2 foreach_cell_root_BFS(Point root)
{
    {
        int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
        Point point = {0};

        #if dimension == 1
            typedef struct { int l, i; } queue_t;
        #elif dimension == 2
            typedef struct { int l, i, j; } queue_t;
        #else // dimension == 3
            int kg = 0; NOT_UNUSED(kg);
            typedef struct { int l, i, j, k; } queue_t;
        #endif

        int current_count = 0;
        int next_count = 0;    

        int current_size_max = 1;
        #if dimension == 1
            int next_size_max = 2;
        #elif dimension == 2
            int next_size_max = 4;
        #else // dimension == 3
            int next_size_max = 8;
        #endif

        queue_t *queue_current = (queue_t *) malloc(current_size_max * sizeof(queue_t));
        queue_t *queue_next    = (queue_t *) malloc(next_size_max    * sizeof(queue_t));

        if (!queue_current || !queue_next) {
            exit(EXIT_FAILURE);
        }

        #if dimension == 1
            queue_current[0] = (queue_t){root.level, root.i};
        #elif dimension == 2
            queue_current[0] = (queue_t){root.level, root.i, root.j};        
        #else // dimension == 3
            queue_current[0] = (queue_t){root.level, root.i, root.j, root.k};
        #endif

        current_count = 1;

        while (current_count > 0) {
            for (int i = 0; i < current_count; i++) {
                point.level = queue_current[i].l;
                #if dimension == 1                
                    point.i = queue_current[i].i;
                #elif dimension == 2
                    point.i = queue_current[i].i;
                    point.j = queue_current[i].j;
                #else // dimension == 3
                    point.i = queue_current[i].i;
                    point.j = queue_current[i].j;
                    point.k = queue_current[i].k;
                #endif

                /* Check if the cell is allocated (macro 'allocated' expects 'point') */
                if (!allocated(0,0,0))
                    continue;

                {...} // User defined code block

                #if dimension == 1
                    if (point.level < grid->depth && !is_leaf(cell)) {
                        queue_next[next_count++] = (queue_t){point.level+1, _LEFT };
                        queue_next[next_count++] = (queue_t){point.level+1, _RIGHT};
                    }
                #elif dimension == 2
                    if (point.level < grid->depth && !is_leaf(cell)) {
                        queue_next[next_count++] = (queue_t){point.level+1,_LEFT ,_BOTTOM};
                        queue_next[next_count++] = (queue_t){point.level+1,_RIGHT,_BOTTOM};
                        queue_next[next_count++] = (queue_t){point.level+1,_LEFT ,_TOP   };
                        queue_next[next_count++] = (queue_t){point.level+1,_RIGHT,_TOP   };
                    }
                #elif dimension == 3
                    if (point.level < grid->depth && !is_leaf(cell)) {
                        queue_next[next_count++] = (queue_t){point.level+1,_LEFT ,_BOTTOM,_BACK };
                        queue_next[next_count++] = (queue_t){point.level+1,_RIGHT,_BOTTOM,_BACK };
                        queue_next[next_count++] = (queue_t){point.level+1,_LEFT ,_TOP   ,_BACK };
                        queue_next[next_count++] = (queue_t){point.level+1,_RIGHT,_TOP   ,_BACK };
                        queue_next[next_count++] = (queue_t){point.level+1,_LEFT ,_BOTTOM,_FRONT};
                        queue_next[next_count++] = (queue_t){point.level+1,_RIGHT,_BOTTOM,_FRONT};
                        queue_next[next_count++] = (queue_t){point.level+1,_LEFT ,_TOP   ,_FRONT};
                        queue_next[next_count++] = (queue_t){point.level+1,_RIGHT,_TOP   ,_FRONT};
                    }
                #endif
            } // End of loop over the current level

           current_count = next_count;
           next_count = 0;
   
           free(queue_current);
           queue_current = malloc(next_size_max * sizeof(queue_t));
           if (!queue_current) {
               exit(EXIT_FAILURE);
           }
           memcpy(queue_current, queue_next, current_count * sizeof(queue_t));
           current_size_max = next_size_max;
           
           #if dimension == 1
               next_size_max *= 2;
           #elif dimension == 2
               next_size_max *= 4;
           #elif dimension == 3
               next_size_max *= 8;
           #endif
           
           free(queue_next);
           queue_next = malloc(next_size_max * sizeof(queue_t));
           if (!queue_next) {
               exit(EXIT_FAILURE);
           }

        } // End of while loop
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
