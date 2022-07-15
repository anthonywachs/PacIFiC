/**
# View front-tracking

This file provides functions to display the Lagrangian mesh.
*/

#include "view.h"

static void begin_draw_vertices (bview * view, float color[3], float ps)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  glColor3f (color[0], color[1], color[2]);
  glEnable(GL_POINT_SMOOTH);
  glPointSize (view->samples*(ps > 0. ? ps : 8.));
  _reversed = view->reversed;
  view->reversed = false;
}

static void end_draw_vertices()
{
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  bview * view = draw();
  view->reversed = _reversed;
}

struct _draw_lag {
  lagMesh* mesh; // Compulsory
  bool nodes;
  bool edges;
  bool facets;
  float fc[3], lc[3], nc[3], lw, ns;
};

void draw_lag(struct _draw_lag p) {
  if (pid() == 0) {
    bool edges = p.edges;
    if (p.lw > 0 || (p.lc[0] > 0 || p.lc[1] > 0 || p.lc[2] > 0)) edges = true;
    bool nodes = p.nodes;
    if (p.ns > 0) nodes = true;
    float my_color[3];
    bview * view = draw();
    if (edges) {
      p.lw = (p.lw) ? p.lw : 2.;
      if (p.lc) {my_color[0] = p.lc[0]; my_color[1] = p.lc[1]; my_color[2] = p.lc[2];}
      else {my_color[0] = 0.; my_color[1] = 0.; my_color[2] = 0.;}
      draw_lines(view, my_color, p.lw) {
        for (int i=0; i<p.mesh->nle; i++) {
          bool across_periodic_bc = false;
          int v1, v2;
          v1 = p.mesh->edges[i].node_ids[0];
          v2 = p.mesh->edges[i].node_ids[1];
          foreach_dimension() {
            if (fabs(p.mesh->nodes[v1].pos.x
              - p.mesh->nodes[v2].pos.x) > L0/2.) across_periodic_bc = true;
          }
          /** If the edge crosses a perdiodic boundary (i.e. the edge length
          is larger than L0/2), we simply don't draw it */
          if (!across_periodic_bc) {
            glBegin(GL_LINES);
              #if dimension < 3
                glvertex2d(view, p.mesh->nodes[v1].pos.x,
                  p.mesh->nodes[v1].pos.y);
                glvertex2d(view, p.mesh->nodes[v2].pos.x,
                  p.mesh->nodes[v2].pos.y);
              #else
                glvertex3d(view, p.mesh->nodes[v1].pos.x,
                  p.mesh->nodes[v1].pos.y, p.mesh->nodes[v1].pos.z);
                glvertex3d(view, p.mesh->nodes[v2].pos.x,
                  p.mesh->nodes[v2].pos.y, p.mesh->nodes[v2].pos.z);
              #endif
            glEnd();
            view->ni++;
          }
        }
      }
    }
    if (nodes) {
      p.ns = (p.ns) ? p.ns : 8.;
      if (p.nc) {my_color[0] = p.nc[0]; my_color[1] = p.nc[1]; my_color[2] = p.nc[2];}
      else {my_color[0] = 0.; my_color[1] = 0.; my_color[2] = 0.;}
      draw_vertices(view, my_color, p.ns) {
        glBegin(GL_POINTS);
          for (int i=0; i<p.mesh->nlp; i++)
            #if dimension < 3
              glvertex2d(view, p.mesh->nodes[i].pos.x,
                p.mesh->nodes[i].pos.y);
            #else
            glvertex3d(view, p.mesh->nodes[i].pos.x,
              p.mesh->nodes[i].pos.y, p.mesh->nodes[i].pos.z);
            #endif
        glEnd();
        view->ni++;
      }
    }
    #if dimension > 2
      bool facets = p.facets;
      if (facets) {
        if (p.fc[0]) {my_color[0] = p.fc[0]; my_color[1] = p.fc[1]; my_color[2] = p.fc[2];}
        else {my_color[0] = 1.; my_color[1] = 1.; my_color[2] = 1.;}
        for (int i=0; i<p.mesh->nlt; i++) {
          int nodes[3];
          for(int j=0; j<3; j++) nodes[j] = p.mesh->triangles[i].node_ids[j];
          if (!(is_triangle_across_periodic(p.mesh, i))) {
            glBegin (GL_POLYGON);
            glColor3f (my_color[0], my_color[1], my_color[2]);
              for(int j=0; j<3; j++) {
                glVertex3d(
                  p.mesh->nodes[nodes[j]].pos.x,
                  p.mesh->nodes[nodes[j]].pos.y,
                  p.mesh->nodes[nodes[j]].pos.z);
                }
            glEnd ();
            view->ni++;
          }
        }
      }
    #endif
  }
}
