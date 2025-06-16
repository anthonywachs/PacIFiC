#define N_NEWTON 10

#define MC 1.3858189
#define MA (RADIUS/MC)
#define MB 0.207
#define MD 2.003
#define ME (-1.123)

double dDdx(double x, double x0, double y0) {
  return (2*(.5*MA*MC*cos(x)*(MD + MB*sq(cos(x)) + ME*sq(sq(cos(x)))) + .5*MA*MC*sin(x)*(-2*MB*sin(x)*cos(x) - 4*ME*sin(x)*cube(cos(x))))*(.5*MA*MC*sin(x)*(MB*sq(cos(x)) + ME*sq(sq(cos(x))) + MD) - y0) - 2*MA*MC*sin(x)*(MA*MC*cos(x) - x0));
}

double d2Ddx2(double x, double x0, double y0) {
  return (2*sq(MA)*sq(MC)*sq(sin(x)) - 2*MA*MC*cos(x)*(MA*MC*cos(x) - x0) + 2*(.5*MA*MC*sin(x)*(MB*sq(cos(x)) + MD + ME*sq(sq(cos(x)))) + MA*MC*cos(x)*(-2*MB*sin(x)*cos(x) - 4*ME*sin(x)*cube(cos(x))) + .5*MA*MC*sin(x)*(2*MB*sq(sin(x)) - 2*MB*sq(cos(x)) - 4*ME*sq(sq(cos(x))) + 12*ME*sq(sin(x))*sq(cos(x)))) * (.5*MA*MC*sin(x)*(MB*sq(cos(x)) + MD + ME*sq(sq(cos(x)))) - y0) + 2*sq((.5*MA*MC*cos(x)*(MB*sq(cos(x)) + MD + ME*sq(sq(cos(x)))) + .5*MA*MC*sin(x)*(-2*MB*sin(x)*cos(x) - 4*ME*sin(x)*cube(cos(x))))));
}

double get_phi(double x0, double y0) {
  double xn = (y0>=0) ? pi/2 : 3*pi/2;
  for(int k=0; k<N_NEWTON; ++k) {
    xn = xn - dDdx(xn, x0, y0)/d2Ddx2(xn, x0, y0);
  }
  return xn;
}

double dist_biconcave(double x0, double y0) {
  double phi;
  phi = get_phi(x0, y0);
  // int msign = abs(x0/(MA*MC))>=1 ? -1 : sign((- sq(y0/MA) + sq(.5*MC*sin(acos(x0/(MA*MC)))*(.207 + 2.003*sq(cos(acos(x0/(MA*MC)))) - 1.123*sq(sq(cos(acos(x0/(MA*MC)))))))));
  return (sqrt(sq(MA*MC*cos(phi) - x0) + sq(.5*MA*MC*sin(phi)*(MB + MD*sq(cos(phi)) + ME*sq(sq(cos(phi)))) - y0)));
}
