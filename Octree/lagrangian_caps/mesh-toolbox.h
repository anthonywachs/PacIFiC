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

/** The function below returns true if the two nodes $i$ and $j$ are neighbors
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
      mesh->nodes[j].edge_ids[mesh->nodes[j].nb_edges] = i;
      mesh->nodes[k].edge_ids[mesh->nodes[k].nb_edges] = i;
      mesh->nodes[j].nb_edges++;
      mesh->nodes[j].nb_neighbors++;
      mesh->nodes[k].nb_edges++;
      mesh->nodes[k].nb_neighbors++;
    }
    return true;
  }
}

/** The function below returns true if the triangle connecting nodes i,j and k
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

struct _write_triangle {
  lagMesh* mesh;
  int tid;
  int i;
  int j;
  int k;
  bool new_mesh;
  bool overwrite;
  int tid_to_replace;
};

/** The function below writes a triangle at index location tid, connecting nodes
i, j and k. */
bool write_triangle(struct _write_triangle p) {
  lagMesh* mesh = p.mesh;
  int tid = p.tid;
  int i = p.i;
  int j = p.j;
  int k = p.k;
  bool new_mesh = (p.new_mesh) ? p.new_mesh : false;
  bool overwrite = (p.overwrite) ? p.overwrite : false;
  int tid_to_replace = (p.tid_to_replace) ? p.tid_to_replace : -1;
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
      for(int m=0; m<mesh->nodes[va].nb_edges; m++) {
        for(int n=0; n<mesh->nodes[vb].nb_edges; n++) {
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
    if (new_mesh) {
      mesh->nodes[i].triangle_ids[mesh->nodes[i].nb_triangles] = tid;
      mesh->nodes[j].triangle_ids[mesh->nodes[j].nb_triangles] = tid;
      mesh->nodes[k].triangle_ids[mesh->nodes[k].nb_triangles] = tid;
      mesh->nodes[i].nb_triangles++;
      mesh->nodes[j].nb_triangles++;
      mesh->nodes[k].nb_triangles++;
    }
    else {
      for(int m=0; m<3; m++) {
        bool updated_tid = false;
        for(int n=0; n<p.mesh->nodes[
          p.mesh->triangles[tid].node_ids[m]].nb_triangles; n++) {
          if (p.mesh->nodes[p.mesh->triangles[tid].node_ids[m]].triangle_ids[n]
            == tid_to_replace) {
              p.mesh->nodes[p.mesh->triangles[tid].node_ids[m]].triangle_ids[n] = tid;
              updated_tid = true;
            }
        }
        if (!updated_tid) {
          for(int n=0; n<p.mesh->nodes[
            p.mesh->triangles[tid].node_ids[m]].nb_triangles; n++) {
            if (p.mesh->nodes[
              p.mesh->triangles[tid].node_ids[m]].triangle_ids[n] == -1) {
                p.mesh->nodes[
                  p.mesh->triangles[tid].node_ids[m]].triangle_ids[n] = tid;
                updated_tid = true;
              }
          }
        }
      }
    }
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
    int ce; // ce for "current edge
    for(int l=0; l<mesh->nodes[nodes[k]].nb_edges; l++) {
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
    int eid;
    for(int j=0; j<mesh->nodes[nodes[i]].nb_edges; j++) {
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
    for(int l=0; l<mesh->nodes[nodes[k]].nb_edges; l++) {
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
  mesh->nodes[mesh->nlp].nb_edges = 6;
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
  for(int j=0; j<mesh->nodes[nid[1]].nb_edges; j++)
    if (mesh->nodes[nid[1]].edge_ids[j] == i)
      mesh->nodes[nid[1]].edge_ids[j] = mesh->nle;

  mesh->nlp++;
  mesh->nle++;
}


/** The function below returns true if node j is a vertex of triangle i */
bool is_triangle_vertex(lagMesh* mesh, int i, int j) {
  for(int k=0; k<3; k++) {
    if (mesh->triangles[i].node_ids[k] == j) return true;
  }
  return false;
}

/** The function below loops through all triangles in the mesh and divide them
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
    /** Shrink the original, big, triangle into the center smaller one */
    overwrite_triangle(mesh, i, mid_ids[0], mid_ids[1], mid_ids[2]);
  }
}
