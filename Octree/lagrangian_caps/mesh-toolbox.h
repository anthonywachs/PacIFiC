/**
# Toolbox to perform operations on a Lagrangian mesh

From defining useful macros, subdividing triangles, or dumping the mesh in
restart files: below is a collection of helpful functions to deal with
triangular meshes.
*/


/**
## Useful macros

The macros below are useful to define an icosahedron
*/
#define GET_LD(NODE) ((fabs(fabs(NODE.pos.x) - ll) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - ll) < 1.e-8 ? 1 : 2)))
#define GET_LD_SIGN(NODE) ((GET_LD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_LD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_SD(NODE) ((fabs(fabs(NODE.pos.x) - sl) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - sl) < 1.e-8 ? 1 : 2)))
#define GET_SD_SIGN(NODE) ((GET_SD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_SD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_ZD(NODE) ((fabs(NODE.pos.x) < 1.e-8) ? 0 : \
  ((fabs(NODE.pos.y) < 1.e-8 ? 1 : 2)))

/**
## Operations on edges

The function below returns true if the two nodes $i$ and $j$ are neighbors
*/
bool is_neighbor(lagMesh* mesh, int i, int j) {
  for(int k=0; k<mesh->nodes[i].nb_neighbors; k++) {
    if (mesh->nodes[i].neighbor_ids[k] == j) return true;
  }
  return false;
}

/** The function below returns true if there is an edge connecting nodes i and
j */
bool edge_exists(lagMesh* mesh, int j, int k) {
  for(int i=0; i<mesh->nle; i++) {
    if ((mesh->edges[i].node_ids[0] == j && mesh->edges[i].node_ids[1] == k)
      || (mesh->edges[i].node_ids[0] == k
      && mesh->edges[i].node_ids[1] == j)) return true;
  }
  return false;
}

/** The function below returns true if the edge $i$ is across a priodic
boundary. */
bool is_edge_across_periodic(lagMesh* mesh, int i) {
  int n[2];
  for(int k=0; k<2; k++) n[k] = mesh->edges[i].node_ids[k];
  if (ACROSS_PERIODIC(mesh->nodes[n[0]].pos.x, mesh->nodes[n[1]].pos.x)
    || ACROSS_PERIODIC(mesh->nodes[n[0]].pos.y, mesh->nodes[n[1]].pos.y)
    || ACROSS_PERIODIC(mesh->nodes[n[0]].pos.z, mesh->nodes[n[1]].pos.z))
      return true;
  return false;
}


struct _write_edge {
  lagMesh* mesh;
  int i;
  int j;
  int k;
  bool new_mesh;
  bool overwrite;
};

/** The function below writes edge i, connecting nodes j and k. If the edge
exists, the function returns false (no edge creation), true otherwise (edge
creation) */
bool write_edge(struct _write_edge p) {
  lagMesh* mesh = p.mesh;
  int i = p.i;
  int j = p.j;
  int k = p.k;
  bool new_mesh = (p.new_mesh) ? p.new_mesh : false;
  bool overwrite = (p.overwrite) ? p.overwrite : false;
  if (!overwrite && edge_exists(mesh, j, k)) return false;
  else {
    mesh->edges[i].node_ids[0] = j;
    mesh->edges[i].node_ids[1] = k;
    for(int ii=0; ii<2; ii++) mesh->edges[i].triangle_ids[ii] = -1;
    if (new_mesh) {
      mesh->nodes[j].neighbor_ids[mesh->nodes[j].nb_neighbors] = k;
      mesh->nodes[k].neighbor_ids[mesh->nodes[k].nb_neighbors] = j;
      mesh->nodes[j].edge_ids[mesh->nodes[j].nb_neighbors] = i;
      mesh->nodes[k].edge_ids[mesh->nodes[k].nb_neighbors] = i;
      mesh->nodes[j].nb_neighbors++;
      mesh->nodes[k].nb_neighbors++;
    }
    return true;
  }
}

/** The function below creates a new edge between nodes i and j, and updates the
connectivity information of its nodes (but not its triangles, since they
don't exist yet). */
void new_edge(lagMesh* mesh, int i, int j) {
  int eid = mesh->nle; // id of the new edge
  int nodes[2];
  nodes[0] = i; nodes[1] = j;
  for(int k=0; k<2; k++) {
    mesh->edges[eid].node_ids[k] = nodes[k];

    /** Add the edge id to the newly connected nodes */
    for(int l=0; l<mesh->nodes[nodes[k]].nb_neighbors; l++) {
      if (mesh->nodes[nodes[k]].edge_ids[l] == -1) {
        mesh->nodes[nodes[k]].edge_ids[l] = eid;
        break;
      }
    }

    /** Update the neighbors' list of the newly connected nodes */
    for(int l=0; l<mesh->nodes[nodes[k]].nb_neighbors; l++) {
      if (mesh->nodes[nodes[k]].neighbor_ids[l] == -1) {
        mesh->nodes[nodes[k]].neighbor_ids[l] = nodes[(k+1)%2];
        break;
      }
    }

    /** The newly created edge is not yet surrounded by any triangle */
    mesh->edges[eid].triangle_ids[k] = -1;
  }
  mesh->nle++;
}

/** The function below splits an edge in two smaller edges, creating a node
at its midpoint. */
void split_edge(lagMesh* mesh, int i) {
  int nid[2];
  for(int j=0; j<2; j++) nid[j] = mesh->edges[i].node_ids[j];

  /** Create new node */
  foreach_dimension()
    mesh->nodes[mesh->nlp].pos.x =
      .5*(mesh->nodes[nid[0]].pos.x + mesh->nodes[nid[1]].pos.x);
  mesh->nodes[mesh->nlp].nb_neighbors = 6;
  mesh->nodes[mesh->nlp].nb_triangles = 6;
  mesh->nodes[mesh->nlp].neighbor_ids[0] = nid[0];
  mesh->nodes[mesh->nlp].neighbor_ids[1] = nid[1];
  mesh->nodes[mesh->nlp].edge_ids[0] = i;
  mesh->nodes[mesh->nlp].edge_ids[1] = mesh->nle;
  for(int j=0; j<6; j++) {
    mesh->nodes[mesh->nlp].triangle_ids[j] = -1;
    if (j>1) {
      mesh->nodes[mesh->nlp].neighbor_ids[j] = -1;
      mesh->nodes[mesh->nlp].edge_ids[j] = -1;
    }
  }

  /** Create new edge and update current one */
  write_edge(mesh, i, nid[0], mesh->nlp, overwrite = true);
  write_edge(mesh, mesh->nle, nid[1], mesh->nlp);
  for (int j=0; j<2; j++) {
    mesh->edges[i].triangle_ids[j] = -1;
    mesh->edges[mesh->nle].triangle_ids[j] = -1;
  }

  /** Update node information: neighboring nodes, connecting edges */
  for(int j=0; j<mesh->nodes[nid[0]].nb_neighbors; j++)
    if (mesh->nodes[nid[0]].neighbor_ids[j] == nid[1])
      mesh->nodes[nid[0]].neighbor_ids[j] = mesh->nlp;
  for(int j=0; j<mesh->nodes[nid[1]].nb_neighbors; j++)
    if (mesh->nodes[nid[1]].neighbor_ids[j] == nid[0])
      mesh->nodes[nid[1]].neighbor_ids[j] = mesh->nlp;
  for(int j=0; j<mesh->nodes[nid[1]].nb_neighbors; j++)
    if (mesh->nodes[nid[1]].edge_ids[j] == i)
      mesh->nodes[nid[1]].edge_ids[j] = mesh->nle;

  mesh->nlp++;
  mesh->nle++;
}

/**
## Operations on triangles

The function below returns true if the triangle connecting nodes i,j and k
already exists in the mesh. */
bool triangle_exists(lagMesh* mesh, int i, int j, int k) {
  for(int t=0; t<mesh->nlt; t++) {
    for(int a=0; a<3; a++) {
      if (mesh->triangles[t].node_ids[a] == i) {
        for(int b=0; b<3; b++) {
          if (b != a && mesh->triangles[t].node_ids[b] == j) {
            for(int c=0; c<3; c++) {
              if (c != a && c != b && mesh->triangles[t].node_ids[c] == k) {
                return true;
              }
            }
          }
        }
      }
    }
  }
  return false;
}

/** The function below returns true if the triangle_ids $i$ is across a periodic
boundary. */
bool is_triangle_across_periodic(lagMesh* mesh, int i) {
  int e[3];
  for(int k=0; k<3; k++) e[k] = mesh->triangles[i].edge_ids[k];
  if (is_edge_across_periodic(mesh, e[0])
    || is_edge_across_periodic(mesh, e[1])
    || is_edge_across_periodic(mesh, e[2]))
      return true;
  return false;
}

struct _write_triangle {
  lagMesh* mesh;
  int tid;
  int i;
  int j;
  int k;
  bool overwrite;
};

/** The function below writes a triangle at index location tid, connecting nodes
i, j and k. It updates the connectivity information of its nodes and edges.*/
bool write_triangle(struct _write_triangle p) {
  lagMesh* mesh = p.mesh;
  int tid = p.tid;
  int i = p.i;
  int j = p.j;
  int k = p.k;
  bool overwrite = (p.overwrite) ? p.overwrite : false;
  if (!overwrite && triangle_exists(mesh, i, j, k)) return false;
  else {
    mesh->triangles[tid].node_ids[0] = i;
    mesh->triangles[tid].node_ids[1] = j;
    mesh->triangles[tid].node_ids[2] = k;
    int c = 0;
    for(int a=0; a<3; a++) {
      int va = mesh->triangles[tid].node_ids[a];
      int b=(a+1)%3;
      int vb = mesh->triangles[tid].node_ids[b];
      for(int m=0; m<mesh->nodes[va].nb_neighbors; m++) {
        for(int n=0; n<mesh->nodes[vb].nb_neighbors; n++) {
          if (mesh->nodes[va].edge_ids[m] == mesh->nodes[vb].edge_ids[n]) {
            mesh->triangles[tid].edge_ids[c] = mesh->nodes[va].edge_ids[m];
            c++;
            int p = (mesh->edges[mesh->nodes[va].edge_ids[m]].triangle_ids[0]
              == -1) ? 0 : 1;
            mesh->edges[mesh->nodes[va].edge_ids[m]].triangle_ids[p] = tid;
          }
        }
      }
    }
    mesh->nodes[i].triangle_ids[mesh->nodes[i].nb_triangles] = tid;
    mesh->nodes[j].triangle_ids[mesh->nodes[j].nb_triangles] = tid;
    mesh->nodes[k].triangle_ids[mesh->nodes[k].nb_triangles] = tid;
    mesh->nodes[i].nb_triangles++;
    mesh->nodes[j].nb_triangles++;
    mesh->nodes[k].nb_triangles++;
    return true;
  }
}

struct _new_triangle {
  lagMesh* mesh;
  int i;
  int j;
  int k;
  int prev_tid;
};

/** The function below also writes a new triangle connecting nodes i, j and k,
but in the specific context of the refinement process. The variable $prev\_tid$
is used to update the connectivity information of the nodes and edges. */
void new_triangle(lagMesh* mesh, int i, int j, int k, int prev_tid) {
  assert((i < mesh->nlp) && (j < mesh->nlp) && (k < mesh->nlp));
  int nodes[3];
  nodes[0] = i; nodes[1] = j; nodes[2] = k;

  int tid = mesh->nlt;
  for(int k=0; k<3; k++) {
    mesh->triangles[tid].node_ids[k] = nodes[k];
    /** Specify the id of the new triangle for node k */
    bool replaced_tid = false;
    for(int l=0; l<mesh->nodes[nodes[k]].nb_triangles; l++)
      if (mesh->nodes[nodes[k]].triangle_ids[l] == prev_tid) {
        mesh->nodes[nodes[k]].triangle_ids[l] = tid;
        replaced_tid = true;
        break;
      }
    if (!replaced_tid) {
      for(int l=0; l<mesh->nodes[nodes[k]].nb_triangles; l++)
        if (mesh->nodes[nodes[k]].triangle_ids[l] == -1) {
          mesh->nodes[nodes[k]].triangle_ids[l] = tid;
          replaced_tid = true;
          break;
        }
    }
    /** Find the edges and: (i) add them to the list of edges of the new
    triangle; (ii) specify the id of the new triangle for the edges */
    /** First, we find the edge that connects node i to node i+1 */
    int ce = -1; // ce for "current edge
    for(int l=0; l<mesh->nodes[nodes[k]].nb_neighbors; l++) {
      ce = mesh->nodes[nodes[k]].edge_ids[l];
      int cn = (mesh->edges[ce].node_ids[0] == nodes[k]) ?
        mesh->edges[ce].node_ids[1] : mesh->edges[ce].node_ids[0];
      if (cn == nodes[(k+1)%3]) break;
    }
    /** Then, we add this edge to the list of edges of our new triangle */
    mesh->triangles[tid].edge_ids[k] = ce;
    /** And we have to update the triangle id of the edge */
    replaced_tid = false;
    for(int l=0; l<2; l++)
      if (mesh->edges[ce].triangle_ids[l] == prev_tid) {
        mesh->edges[ce].triangle_ids[l] = tid;
        replaced_tid = true;
        break;
      }
    if (!replaced_tid && mesh->edges[ce].triangle_ids[0] != tid) {
      for(int l=0; l<2; l++)
        if (mesh->edges[ce].triangle_ids[l] == -1) {
          mesh->edges[ce].triangle_ids[l] = tid;
          replaced_tid = true;
          break;
        }
    }
  }
  mesh->nlt++;
}

void overwrite_triangle(lagMesh* mesh, int tid, int i, int j, int k) {
  int nodes[3];
  nodes[0] = i; nodes[1] = j; nodes[2] = k;
  for(int i=0; i<3; i++) {
    /** 1. Take care of the nodes of the triangle */
    mesh->triangles[tid].node_ids[i] = nodes[i];
    /** Update the list of triangles for the node $i$ */
    bool already_there = false;
    for(int j=0; j<mesh->nodes[nodes[i]].nb_triangles; j++) {
      if (mesh->nodes[nodes[i]].triangle_ids[j] == tid) already_there = true;
      else if (mesh->nodes[nodes[i]].triangle_ids[j] == -1 && !already_there) {
        mesh->nodes[nodes[i]].triangle_ids[j] = tid;
        break;
      }
    }

    /** 2. Take care of the egdes of the triangle */
    /** 2.1. Identify the edge id connecting nodes $i$, $i+1$ */
    int eid = -1; // eid for "edge id"
    for(int j=0; j<mesh->nodes[nodes[i]].nb_neighbors; j++) {
      eid = mesh->nodes[nodes[i]].edge_ids[j];
      int cn = (mesh->edges[eid].node_ids[0] == nodes[i]) ?
        mesh->edges[eid].node_ids[1] : mesh->edges[eid].node_ids[0];
      if (cn == nodes[(i+1)%3]) break;
    }
    /** 2.2 Add this edge to the triangle's list of edges */
    mesh->triangles[tid].edge_ids[i] = eid;
    /** 2.3 Add the triangle id to the edge's list of triangles */
    int index = (mesh->edges[eid].triangle_ids[0] > -1 &&
      mesh->edges[eid].triangle_ids[0] != tid) ? 1 : 0;
    mesh->edges[eid].triangle_ids[index] = tid;
  }
}

/** The function below returns true if node j is a vertex of triangle i */
bool is_triangle_vertex(lagMesh* mesh, int i, int j) {
  for(int k=0; k<3; k++) {
    if (mesh->triangles[i].node_ids[k] == j) return true;
  }
  return false;
}

/**
The function below returns the (positive) angle between the two vectors formed
by the nodes [*n1, *n2] and [*n1, *n3].
*/
double comp_angle(lagNode* n1, lagNode* n2, lagNode* n3) {
  double theta = 0.;
  foreach_dimension() theta += (n1->pos.x - n2->pos.x)*(n1->pos.x - n3->pos.x);
  double norm1, norm2;
  norm1 = 0.; norm2 = 0.;
  foreach_dimension() {
    norm1 += sq(n1->pos.x - n2->pos.x);
    norm2 += sq(n1->pos.x - n3->pos.x);
  }
  theta /= sqrt(norm1)*sqrt(norm2);
  theta = acos(fabs(theta));
  return theta;
}

/**
The function below returns true if the angle of node $i$ ($0 < i < 2$) of
triangle $tid$ is greater than $\pi$ radians.
*/
bool is_obtuse_node(lagMesh* mesh, int tid, int i) {
  lagNode* n[3];
  for(int j=0; j<3; j++)
    n[j] = &(mesh->nodes[mesh->triangles[tid].node_ids[j]]);
  if (comp_angle(n[i], n[(i+1)%3], n[(i+2)%3]) > pi) return true;
  else return false;
}

/**
The function below returns true if the triangle $tid$ is obtuse at any of its
three angles.
*/
bool is_obtuse_triangle(lagMesh* mesh, int tid) {
  for(int i=0; i<3; i++) if (is_obtuse_node(mesh, tid, i)) return true;
  return false;
}

/**
## Uniform refinement of a mesh by subdividing its triangles

The function below loops through all triangles in the mesh and divide them
in four smaller ones. Keeping the correct structure of the mesh, i.e. updating
the relationship between nodes, edges, triangles and their neighbors results
in the somewhat complicated implementation below. */
void refine_mesh(lagMesh* mesh) {
  /** Perform the loop subdivision algorithm: until we reach the desired number
  of nodes, we split each triangles into four smaller ones */
  int cnt = mesh->nlt; // current number of triangles
  for(int i=0; i<cnt; i++) {
    int mid_ids[3];
    /** If not done yet, we split each edge into two */
    for(int j=0; j<3; j++) {
      int edge_id = mesh->triangles[i].edge_ids[j];
      if (mesh->edges[edge_id].triangle_ids[1] > -1) {
        split_edge(mesh, edge_id);
        mid_ids[j] = mesh->nlp-1;
      }
      else {
        mid_ids[j] = is_triangle_vertex(mesh, i,
          mesh->edges[edge_id].node_ids[0]) ?
          mesh->edges[edge_id].node_ids[1] :
          mesh->edges[edge_id].node_ids[0];
      }
    }

    /** Connect the three midpoints with edges, and create corner triangles */
    for(int j=0; j<3; j++) {
      /** create edges between midpoints */
      new_edge(mesh, mid_ids[j], mid_ids[(j+1)%3]);

      /** Create the corner triangle with the new edge */
      int corner_id = -1;
      for(int k=0; k<3; k++) { // Loop over the current triangle vertices
        corner_id = mesh->triangles[i].node_ids[k];
        if (is_neighbor(mesh, corner_id, mid_ids[j]) &&
          is_neighbor(mesh, corner_id, mid_ids[(j+1)%3])) {
          break;
        }
      }
      new_triangle(mesh, mid_ids[j], mid_ids[(j+1)%3], corner_id, i);
    }
    /** Shrink the original (big) triangle into the center smaller one */
    overwrite_triangle(mesh, i, mid_ids[0], mid_ids[1], mid_ids[2]);
  }
}

/**
## Adding capsules after restart

In case we add membranes to a simulation after restart, we need to call
the function below. */
void initialize_membranes() {
  mbs.nbmb = NCAPS;
  for(int i=0; i<mbs.nbmb; i++) {
    initialize_empty_mb(&mbs.mb[i]);
  }
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

void initialize_membranes_stencils() {
  for(int i=0; i<mbs.nbmb; i++) {
    for(int j=0; j<MB(i).nlp; j++) {
      MB(i).nodes[j].stencil.n = STENCIL_SIZE;
      MB(i).nodes[j].stencil.nm = STENCIL_SIZE;
      MB(i).nodes[j].stencil.p = (Index*) malloc(STENCIL_SIZE*sizeof(Index));
    }
  }
}

/**
## Output
The functions below write/read the lagrangian mesh to/from a file in order
to restart simulations.

At the moment these functions are only valid for three-dimensional simulations.
*/

void dump_lagmesh(FILE* fp, lagMesh* mesh) {
  /** Step 1: The array of Lagrangian nodes is dumped */
  fwrite(&(mesh->nlp), sizeof(int), 1, fp);
  for(int i=0; i<mesh->nlp; i++) {
    foreach_dimension() fwrite(&(mesh->nodes[i].pos.x), sizeof(double), 1, fp);
    fwrite(&(mesh->nodes[i].ref_curv), sizeof(double), 1, fp);
    fwrite(&(mesh->nodes[i].nb_neighbors), sizeof(int), 1, fp);
    for(int j=0; j<6; j++)
      fwrite(&(mesh->nodes[i].neighbor_ids[j]), sizeof(int), 1, fp);
    fwrite(&(mesh->nodes[i].nb_neighbors), sizeof(int), 1, fp);
    for(int j=0; j<6; j++)
      fwrite(&(mesh->nodes[i].edge_ids[j]), sizeof(int), 1, fp);
    fwrite(&(mesh->nodes[i].nb_triangles), sizeof(int), 1, fp);
    for(int j=0; j<6; j++)
      fwrite(&(mesh->nodes[i].triangle_ids[j]), sizeof(int), 1, fp);
  }
  /** Step 2: The array of Lagrangian edges is dumped */
  fwrite(&(mesh->nle), sizeof(int), 1, fp);
  for(int i=0; i<mesh->nle; i++) {
    fwrite(&(mesh->edges[i].l0), sizeof(double), 1, fp);
    for(int j=0; j<2; j++) {
      fwrite(&(mesh->edges[i].node_ids[j]), sizeof(int), 1, fp);
      fwrite(&(mesh->edges[i].triangle_ids[j]), sizeof(int), 1, fp);
    }
  }
  /** Step 3: The array of Lagrangian triangles is dumped. If the finite element
  functions are not used, some attributes of the Triangle structure are not
  initialized and we don't have to dump them. */
  fwrite(&(mesh->nlt), sizeof(int), 1, fp);
  for(int i=0; i<mesh->nlt; i++) {
    #if _ELASTICITY_FT
    for(int j=0; j<2; j++) {
      foreach_dimension()
        fwrite(&(mesh->triangles[i].refShape[j].x), sizeof(double), 1, fp);
    }
    for(int j=0; j<3; j++) {
      for(int k=0; k<2; k++) {
        fwrite(&(mesh->triangles[i].sfc[j][k]), sizeof(double), 1, fp);
      }
    }
    #endif
    for(int j=0; j<3; j++)
      fwrite(&(mesh->triangles[i].node_ids[j]), sizeof(int), 1, fp);
    for(int j=0; j<3; j++)
      fwrite(&(mesh->triangles[i].edge_ids[j]), sizeof(int), 1, fp);
  }
}

void restore_lagmesh(FILE* fp, lagMesh* mesh) {
  /** Step 1: The array of Lagrangian nodes is read */
  fread(&(mesh->nlp), sizeof(int), 1, fp);
  mesh->nodes = malloc(mesh->nlp*sizeof(lagNode));
  for(int i=0; i<mesh->nlp; i++) {
    foreach_dimension() fread(&(mesh->nodes[i].pos.x), sizeof(double), 1, fp);
    fread(&(mesh->nodes[i].ref_curv), sizeof(double), 1, fp);
    fread(&(mesh->nodes[i].nb_neighbors), sizeof(int), 1, fp);
    for(int j=0; j<6; j++)
      fread(&(mesh->nodes[i].neighbor_ids[j]), sizeof(int), 1, fp);
    fread(&(mesh->nodes[i].nb_neighbors), sizeof(int), 1, fp);
    for(int j=0; j<6; j++)
      fread(&(mesh->nodes[i].edge_ids[j]), sizeof(int), 1, fp);
    fread(&(mesh->nodes[i].nb_triangles), sizeof(int), 1, fp);
    for(int j=0; j<6; j++)
      fread(&(mesh->nodes[i].triangle_ids[j]), sizeof(int), 1, fp);
  }
  /** Step 1: The array of Lagrangian edges is read */
  fread(&(mesh->nle), sizeof(int), 1, fp);
  mesh->edges = malloc(mesh->nle*sizeof(Edge));
  for(int i=0; i<mesh->nle; i++) {
    fread(&(mesh->edges[i].l0), sizeof(double), 1, fp);
    for(int j=0; j<2; j++) {
      fread(&(mesh->edges[i].node_ids[j]), sizeof(int), 1, fp);
      fread(&(mesh->edges[i].triangle_ids[j]), sizeof(int), 1, fp);
    }
  }
  /** Step 3: The array of Lagrangian triangles is read. If the finite element
  functions are not used, some attributes of the Triangle structure are not
  dumped and we should not read them. */
  fread(&(mesh->nlt), sizeof(int), 1, fp);
  mesh->triangles = malloc(mesh->nlt*sizeof(Triangle));
  for(int i=0; i<mesh->nlt; i++) {
    #if _ELASTICITY_FT
    for(int j=0; j<2; j++) {
      foreach_dimension()
        fread(&(mesh->triangles[i].refShape[j].x), sizeof(double), 1, fp);
    }
    for(int j=0; j<3; j++) {
      for(int k=0; k<2; k++) {
        fread(&(mesh->triangles[i].sfc[j][k]), sizeof(double), 1, fp);
      }
    }
    #endif
    for(int j=0; j<3; j++)
      fread(&(mesh->triangles[i].node_ids[j]), sizeof(int), 1, fp);
    for(int j=0; j<3; j++)
      fread(&(mesh->triangles[i].edge_ids[j]), sizeof(int), 1, fp);
  }
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
}

/** If the simulation contains several membranes, we dump and read one mesh at
a time, but all meshes are stored in the same file. */
void dump_membranes(char* filename) {
  FILE* file = fopen(filename, "w");
  assert(file);
  for(int i=0; i<mbs.nbmb; i++) {
    dump_lagmesh(file, &MB(i));
  }
  fclose(file);
}

void restore_membranes(char* filename) {
  FILE* file = fopen(filename, "r");
  assert(file);
  for(int i=0; i<mbs.nbmb; i++) {
    restore_lagmesh(file, &MB(i));
  }
  fclose(file);
}

/** ### Output membrane position for external post-processing */
void dump_plain_nodes_pos(lagMesh* mesh, char* filename) {
  if (pid() == 0) {
    FILE* file = fopen(filename, "a");
    assert(file);
    fprintf(file, "%g", t);
    for(int i=0; i<mesh->nlp; i++) {
      fprintf(file, ",%g %g %g", mesh->nodes[i].pos.x, mesh->nodes[i].pos.y,
        mesh->nodes[i].pos.z);
    }
    fprintf(file, "\n");
    fclose(file);
  }
}

/** ### Visualization in paraview */
#ifndef PARAVIEW_CAPSULE
  #define PARAVIEW_CAPSULE 0
#endif

#if PARAVIEW_CAPSULE
int pv_timestep = 0;

void pv_output_ascii() {
  char filename[128];
  FILE* file;
  sprintf(filename, "caps_T%d.vtk", pv_timestep);
  file = fopen(filename, "w");

  /** Populate the header and other non-data fields */
  fprintf(file, "# vtk DataFile Version 4.2\n");
  fprintf(file, "Capsules at time %g\n", t);
  fprintf(file, "ASCII\n");
  fprintf(file, "DATASET POLYDATA\n");

  /** Populate the coordinates of all the Lagrangian nodes */
  int nbpts_tot = 0;
  for(int j=0; j<NCAPS; j++) nbpts_tot += MB(j).nlp;
  fprintf(file, "POINTS %d double\n", nbpts_tot);
  for(int j=0; j<NCAPS; j++) {
    for(int k=0; k<MB(j).nlp; k++) {
      fprintf(file, "%g %g %g\n", MB(j).nodes[k].pos.x, MB(j).nodes[k].pos.y,
        MB(j).nodes[k].pos.z);
    }
  }

  /** Populate the connectivity of the triangles */
  int nbtri_tot = 0;
  for(int j=0; j<NCAPS; j++) nbtri_tot += MB(j).nlt;
  fprintf(file, "TRIANGLE_STRIPS %d %d\n", nbtri_tot, 4*nbtri_tot);
  for(int j=0; j<NCAPS; j++) {
    int offset = 0;
    for(int k=0; k<MB(j).nlt; k++) {
      fprintf(file, "%d %d %d %d\n", 3,
        offset + MB(j).triangles[k].node_ids[0],
        offset + MB(j).triangles[k].node_ids[1],
        offset + MB(j).triangles[k].node_ids[2]);
    }
    offset += MB(j).nlp;
  }

  fclose(file);
  pv_timestep++;
}
#endif
