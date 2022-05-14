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
