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
  coord* lagMesh; // Compulsory
  bool edges;
  bool vertices;
  float fc[3], lc[3], vc[3], lw, vs;
};

void draw_lag(struct _draw_lag p) {
  bool edges = true;
  bool vertices = false;
  if (p.edges) edges = true;
  if (p.vertices) vertices = true;
  float color[3];
  bview * view = draw();
  if (edges) {
    p.lw = (p.lw) ? p.lw : 2.;
    if (p.lc) {color[0] = p.lc[0]; color[1] = p.lc[1]; color[2] = p.lc[2];}
    else {color[0] = 0.; color[1] = 0.; color[2] = 0.;}
    draw_lines(view, color, p.lw) {
      foreach_visible(view) {
        glBegin(GL_LINE_LOOP);
          for (int i=0; i<NLP-1; i++)
            glvertex2d(view, p.lagMesh[i].x, p.lagMesh[i].y);
        glEnd();
        view->ni++;
      }
    }
  }
  if (vertices || p.vc || p.vs) {
    p.vs = (p.vs) ? p.vs : 8.;
    if (p.vc) {color[0] = p.vc[0]; color[1] = p.vc[1]; color[2] = p.vc[2];}
    else {color[0] = 0.; color[1] = 0.; color[2] = 0.;}
    draw_vertices(view, color, p.vs) {
      glBegin(GL_POINTS);
        for (int i=0; i<NLP-1; i++)
          glvertex2d(view, p.lagMesh[i].x, p.lagMesh[i].y);
      glEnd();
      view->ni++;
    }
  }
}
