//                           _____      _ _  __                                   //
//                          / ____|    | | |/ /                                   //
//    ___  _ __   ___ _ __ | |  __     | | ' /                                    //
//   / _ \| '_ \ / _ \ '_ \| | |_ |_   | |  <                                     //
//  | (_) | |_) |  __/ | | | |__| | |__| | . \                                    //
//   \___/| .__/ \___|_| |_|\_____|\____/|_|\_\                                   //
//        | |                                                                     //
//        |_|                                                                     //
//                                                                                //
// Copyright 2022 Mattia Montanari, University of Oxford                          //
//                                                                                //
// This program is free software: you can redistribute it and/or modify it under  //
// the terms of the GNU General Public License as published by the Free Software  //
// Foundation, either version 3 of the License. You should have received a copy   //
// of the GNU General Public License along with this program. If not, visit       //
//                                                                                //
//     https://www.gnu.org/licenses/                                              //
//                                                                                //
// This program is distributed in the hope that it will be useful, but WITHOUT    //
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  //
// FOR A PARTICULAR PURPOSE. See GNU General Public License for details.          //

/**
 * @file openGJK.c
 * @author Mattia Montanari
 * @date 1 Jan 2022
 * @brief Source of OpenGJK and its fast sub-algorithm.
 *
 * @see https://www.mattiamontanari.com/opengjk/
 */

#include "Matrix.hh"
#include "openGJK.hh"
#include "GrainsExec.hh"


// We will use the given tolerance in the XML input file instead
// #define eps_rel22        (gkFloat) gkEpsilon * 1e4f
// #define eps_tot22        (gkFloat) gkEpsilon * 1e2f

#define norm2(a)         (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
#define dotProduct(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

#define S3Dregion1234()                                                                                                \
  v[0] = 0;                                                                                                            \
  v[1] = 0;                                                                                                            \
  v[2] = 0;                                                                                                            \
  s->nvrtx = 4;

#define select_1ik()                                                                                                   \
  s->nvrtx = 3;                                                                                                        \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[2][t] = s->vrtx[3][t];                                                                                     \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[1][t] = si[t];                                                                                             \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[0][t] = sk[t];

#define select_1ij()                                                                                                   \
  s->nvrtx = 3;                                                                                                        \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[2][t] = s->vrtx[3][t];                                                                                     \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[1][t] = si[t];                                                                                             \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[0][t] = sj[t];

#define select_1jk()                                                                                                   \
  s->nvrtx = 3;                                                                                                        \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[2][t] = s->vrtx[3][t];                                                                                     \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[1][t] = sj[t];                                                                                             \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[0][t] = sk[t];

#define select_1i()                                                                                                    \
  s->nvrtx = 2;                                                                                                        \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[1][t] = s->vrtx[3][t];                                                                                     \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[0][t] = si[t];

#define select_1j()                                                                                                    \
  s->nvrtx = 2;                                                                                                        \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[1][t] = s->vrtx[3][t];                                                                                     \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[0][t] = sj[t];

#define select_1k()                                                                                                    \
  s->nvrtx = 2;                                                                                                        \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[1][t] = s->vrtx[3][t];                                                                                     \
  for (t = 0; t < 3; t++)                                                                                              \
    s->vrtx[0][t] = sk[t];

#define getvrtx(point, location)                                                                                       \
  point[0] = s->vrtx[location][0];                                                                                     \
  point[1] = s->vrtx[location][1];                                                                                     \
  point[2] = s->vrtx[location][2];

#define calculateEdgeVector(p1p2, p2)                                                                                  \
  p1p2[0] = p2[0] - s->vrtx[3][0];                                                                                     \
  p1p2[1] = p2[1] - s->vrtx[3][1];                                                                                     \
  p1p2[2] = p2[2] - s->vrtx[3][2];

#define S1Dregion1()                                                                                                   \
  v[0] = s->vrtx[1][0];                                                                                                \
  v[1] = s->vrtx[1][1];                                                                                                \
  v[2] = s->vrtx[1][2];                                                                                                \
  s->nvrtx = 1;                                                                                                        \
  s->vrtx[0][0] = s->vrtx[1][0];                                                                                       \
  s->vrtx[0][1] = s->vrtx[1][1];                                                                                       \
  s->vrtx[0][2] = s->vrtx[1][2];

#define S2Dregion1()                                                                                                   \
  v[0] = s->vrtx[2][0];                                                                                                \
  v[1] = s->vrtx[2][1];                                                                                                \
  v[2] = s->vrtx[2][2];                                                                                                \
  s->nvrtx = 1;                                                                                                        \
  s->vrtx[0][0] = s->vrtx[2][0];                                                                                       \
  s->vrtx[0][1] = s->vrtx[2][1];                                                                                       \
  s->vrtx[0][2] = s->vrtx[2][2];

#define S2Dregion12()                                                                                                  \
  s->nvrtx = 2;                                                                                                        \
  s->vrtx[0][0] = s->vrtx[2][0];                                                                                       \
  s->vrtx[0][1] = s->vrtx[2][1];                                                                                       \
  s->vrtx[0][2] = s->vrtx[2][2];

#define S2Dregion13()                                                                                                  \
  s->nvrtx = 2;                                                                                                        \
  s->vrtx[1][0] = s->vrtx[2][0];                                                                                       \
  s->vrtx[1][1] = s->vrtx[2][1];                                                                                       \
  s->vrtx[1][2] = s->vrtx[2][2];

#define S3Dregion1()                                                                                                   \
  v[0] = s1[0];                                                                                                        \
  v[1] = s1[1];                                                                                                        \
  v[2] = s1[2];                                                                                                        \
  s->nvrtx = 1;                                                                                                        \
  s->vrtx[0][0] = s1[0];                                                                                               \
  s->vrtx[0][1] = s1[1];                                                                                               \
  s->vrtx[0][2] = s1[2];

inline gkFloat
determinant(const gkFloat* __restrict__ p, const gkFloat* __restrict__ q, const gkFloat* __restrict__ r) {
  return ( p[0] * ((q[1] * r[2]) - (r[1] * q[2])) - p[1] * (q[0] * r[2] - r[0] * q[2])
         + p[2] * (q[0] * r[1] - r[0] * q[1]) );
}

inline void
crossProduct(const gkFloat* __restrict__ a, const gkFloat* __restrict__ b, gkFloat* __restrict__ c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

inline void
projectOnLine(const gkFloat* __restrict__ p, const gkFloat* __restrict__ q, gkFloat* __restrict__ v) {
  gkFloat pq[3];
  pq[0] = p[0] - q[0];
  pq[1] = p[1] - q[1];
  pq[2] = p[2] - q[2];

  const gkFloat tmp = dotProduct(p, pq) / dotProduct(pq, pq);

  for (int i = 0; i < 3; i++) {
    v[i] = p[i] - pq[i] * tmp;
  }
}

inline void
projectOnPlane(const gkFloat* __restrict__ p, const gkFloat* __restrict__ q, const gkFloat* __restrict__ r, gkFloat* __restrict__ v) {
  gkFloat n[3], pq[3], pr[3];

  for (int i = 0; i < 3; i++) {
    pq[i] = p[i] - q[i];
  }
  for (int i = 0; i < 3; i++) {
    pr[i] = p[i] - r[i];
  }

  crossProduct(pq, pr, n);
  const gkFloat tmp = dotProduct(n, p) / dotProduct(n, n);

  for (int i = 0; i < 3; i++) {
    v[i] = n[i] * tmp;
  }
}

inline  int
hff1(const gkFloat* __restrict__ p, const gkFloat* __restrict__ q) {
  gkFloat tmp = 0;

  for (int i = 0; i < 3; i++) {
    tmp += (p[i] * p[i] - p[i] * q[i]);
  }

  if (tmp > 0) {
    return 1; // keep q
  }
  return 0;
}

inline  int
hff2(const gkFloat* __restrict__ p, const gkFloat* __restrict__ q, const gkFloat* __restrict__ r) {
  gkFloat ntmp[3];
  gkFloat n[3], pq[3], pr[3];

  for (int i = 0; i < 3; i++) {
    pq[i] = q[i] - p[i];
  }
  for (int i = 0; i < 3; i++) {
    pr[i] = r[i] - p[i];
  }

  crossProduct(pq, pr, ntmp);
  crossProduct(pq, ntmp, n);

  return dotProduct(p, n) < 0; // Discard r if true
}

inline int
hff3(const gkFloat* __restrict__ p, const gkFloat* __restrict__ q, const gkFloat* __restrict__ r) {
  gkFloat n[3], pq[3], pr[3];

  for (int i = 0; i < 3; i++) {
    pq[i] = q[i] - p[i];
  }
  for (int i = 0; i < 3; i++) {
    pr[i] = r[i] - p[i];
  }

  crossProduct(pq, pr, n);
  return dotProduct(p, n) <= 0; // discard s if true
}

inline void
S1D(gkSimplex* s, gkFloat* v) {
  const gkFloat* __restrict__ s1p = s->vrtx[1];
  const gkFloat* __restrict__ s2p = s->vrtx[0];

  if (hff1(s1p, s2p)) {
    projectOnLine(s1p, s2p, v); // Update v, no need to update s
    return;                     // Return V{1,2}
  } else {
    S1Dregion1(); // Update v and s
    return;       // Return V{1}
  }
}

inline void
S2D(gkSimplex* s, gkFloat* v) {
  const gkFloat* s1p = s->vrtx[2];
  const gkFloat* s2p = s->vrtx[1];
  const gkFloat* s3p = s->vrtx[0];
  const int hff1f_s12 = hff1(s1p, s2p);
  const int hff1f_s13 = hff1(s1p, s3p);

  if (hff1f_s12) {
    const int hff2f_23 = !hff2(s1p, s2p, s3p);
    if (hff2f_23) {
      if (hff1f_s13) {
        const int hff2f_32 = !hff2(s1p, s3p, s2p);
        if (hff2f_32) {
          projectOnPlane(s1p, s2p, s3p, v); // Update s, no need to update c
          return;                           // Return V{1,2,3}
        } else {
          projectOnLine(s1p, s3p, v); // Update v
          S2Dregion13();              // Update s
          return;                     // Return V{1,3}
        }
      } else {
        projectOnPlane(s1p, s2p, s3p, v); // Update s, no need to update c
        return;                           // Return V{1,2,3}
      }
    } else {
      projectOnLine(s1p, s2p, v); // Update v
      S2Dregion12();              // Update s
      return;                     // Return V{1,2}
    }
  } else if (hff1f_s13) {
    const int hff2f_32 = !hff2(s1p, s3p, s2p);
    if (hff2f_32) {
      projectOnPlane(s1p, s2p, s3p, v); // Update s, no need to update v
      return;                           // Return V{1,2,3}
    } else {
      projectOnLine(s1p, s3p, v); // Update v
      S2Dregion13();              // Update s
      return;                     // Return V{1,3}
    }
  } else {
    S2Dregion1(); // Update s and v
    return;       // Return V{1}
  }
}

inline void
S3D(gkSimplex* s, gkFloat* v) {
  gkFloat s1[3], s2[3], s3[3], s4[3], s1s2[3], s1s3[3], s1s4[3];
  gkFloat si[3], sj[3], sk[3];
  int testLineThree, testLineFour, testPlaneTwo, testPlaneThree, testPlaneFour, dotTotal;
  int i, j, k, t;

  getvrtx(s1, 3);
  getvrtx(s2, 2);
  getvrtx(s3, 1);
  getvrtx(s4, 0);
  calculateEdgeVector(s1s2, s2);
  calculateEdgeVector(s1s3, s3);
  calculateEdgeVector(s1s4, s4);

  int hff1_tests[3];
  hff1_tests[2] = hff1(s1, s2);
  hff1_tests[1] = hff1(s1, s3);
  hff1_tests[0] = hff1(s1, s4);
  testLineThree = hff1(s1, s3);
  testLineFour = hff1(s1, s4);

  dotTotal = hff1(s1, s2) + testLineThree + testLineFour;
  if (dotTotal == 0) { /* case 0.0 -------------------------------------- */
    S3Dregion1();
    return;
  }

  const gkFloat det134 = determinant(s1s3, s1s4, s1s2);
  const int sss = (det134 <= 0);

  testPlaneTwo = hff3(s1, s3, s4) - sss;
  testPlaneTwo = testPlaneTwo * testPlaneTwo;
  testPlaneThree = hff3(s1, s4, s2) - sss;
  testPlaneThree = testPlaneThree * testPlaneThree;
  testPlaneFour = hff3(s1, s2, s3) - sss;
  testPlaneFour = testPlaneFour * testPlaneFour;

  switch (testPlaneTwo + testPlaneThree + testPlaneFour) {
    case 3:
      S3Dregion1234();
      break;

    case 2:
      // Only one facing the oring
      // 1,i,j, are the indices of the points on the triangle and remove k from
      // simplex
      s->nvrtx = 3;
      if (!testPlaneTwo) { // k = 2;   removes s2
        for (i = 0; i < 3; i++) {
          s->vrtx[2][i] = s->vrtx[3][i];
        }
      } else if (!testPlaneThree) { // k = 1; // removes s3
        for (i = 0; i < 3; i++) {
          s->vrtx[1][i] = s2[i];
        }
        for (i = 0; i < 3; i++) {
          s->vrtx[2][i] = s->vrtx[3][i];
        }
      } else if (!testPlaneFour) { // k = 0; // removes s4  and no need to reorder
        for (i = 0; i < 3; i++) {
          s->vrtx[0][i] = s3[i];
        }
        for (i = 0; i < 3; i++) {
          s->vrtx[1][i] = s2[i];
        }
        for (i = 0; i < 3; i++) {
          s->vrtx[2][i] = s->vrtx[3][i];
        }
      }
      // Call S2D
      S2D(s, v);
      break;
    case 1:
      // Two triangles face the origins:
      //    The only positive hff3 is for triangle 1,i,j, therefore k must be in
      //    the solution as it supports the the point of minimum norm.

      // 1,i,j, are the indices of the points on the triangle and remove k from
      // simplex
      s->nvrtx = 3;
      if (testPlaneTwo) {
        k = 2; // s2
        i = 1;
        j = 0;
      } else if (testPlaneThree) {
        k = 1; // s3
        i = 0;
        j = 2;
      } else {
        k = 0; // s4
        i = 2;
        j = 1;
      }

      getvrtx(si, i);
      getvrtx(sj, j);
      getvrtx(sk, k);

      if (dotTotal == 1) {
        if (hff1_tests[k]) {
          if (!hff2(s1, sk, si)) {
            select_1ik();
            projectOnPlane(s1, si, sk, v);
          } else if (!hff2(s1, sk, sj)) {
            select_1jk();
            projectOnPlane(s1, sj, sk, v);
          } else {
            select_1k(); // select region 1i
            projectOnLine(s1, sk, v);
          }
        } else if (hff1_tests[i]) {
          if (!hff2(s1, si, sk)) {
            select_1ik();
            projectOnPlane(s1, si, sk, v);
          } else {
            select_1i(); // select region 1i
            projectOnLine(s1, si, v);
          }
        } else {
          if (!hff2(s1, sj, sk)) {
            select_1jk();
            projectOnPlane(s1, sj, sk, v);
          } else {
            select_1j(); // select region 1i
            projectOnLine(s1, sj, v);
          }
        }
      } else if (dotTotal == 2) {
        // Two edges have positive hff1, meaning that for two edges the origin's
        // project fall on the segement.
        //  Certainly the edge 1,k supports the the point of minimum norm, and so
        //  hff1_1k is positive

        if (hff1_tests[i]) {
          if (!hff2(s1, sk, si)) {
            if (!hff2(s1, si, sk)) {
              select_1ik(); // select region 1ik
              projectOnPlane(s1, si, sk, v);
            } else {
              select_1k(); // select region 1k
              projectOnLine(s1, sk, v);
            }
          } else {
            if (!hff2(s1, sk, sj)) {
              select_1jk(); // select region 1jk
              projectOnPlane(s1, sj, sk, v);
            } else {
              select_1k(); // select region 1k
              projectOnLine(s1, sk, v);
            }
          }
        } else if (hff1_tests[j]) { //  there is no other choice
          if (!hff2(s1, sk, sj)) {
            if (!hff2(s1, sj, sk)) {
              select_1jk(); // select region 1jk
              projectOnPlane(s1, sj, sk, v);
            } else {
              select_1j(); // select region 1j
              projectOnLine(s1, sj, v);
            }
          } else {
            if (!hff2(s1, sk, si)) {
              select_1ik(); // select region 1ik
              projectOnPlane(s1, si, sk, v);
            } else {
              select_1k(); // select region 1k
              projectOnLine(s1, sk, v);
            }
          }
        } else {
          // ERROR;
        }

      } else if (dotTotal == 3) {
        // MM : ALL THIS HYPHOTESIS IS FALSE
        // sk is s.t. hff3 for sk < 0. So, sk must support the origin because
        // there are 2 triangles facing the origin.

        int hff2_ik = hff2(s1, si, sk);
        int hff2_jk = hff2(s1, sj, sk);
        int hff2_ki = hff2(s1, sk, si);
        int hff2_kj = hff2(s1, sk, sj);

        if (hff2_ki == 0 && hff2_kj == 0) {
          // mexPrintf("\n\n UNEXPECTED VALUES!!! \n\n");
        }
        if (hff2_ki == 1 && hff2_kj == 1) {
          select_1k();
          projectOnLine(s1, sk, v);
        } else if (hff2_ki) {
          // discard i
          if (hff2_jk) {
            // discard k
            select_1j();
            projectOnLine(s1, sj, v);
          } else {
            select_1jk();
            projectOnPlane(s1, sk, sj, v);
          }
        } else {
          // discard j
          if (hff2_ik) {
            // discard k
            select_1i();
            projectOnLine(s1, si, v);
          } else {
            select_1ik();
            projectOnPlane(s1, sk, si, v);
          }
        }
      }
      break;

    case 0:
      // The origin is outside all 3 triangles
      if (dotTotal == 1) {
        // Here si is set such that hff(s1,si) > 0
        if (testLineThree) {
          k = 2;
          i = 1; // s3
          j = 0;
        } else if (testLineFour) {
          k = 1; // s3
          i = 0;
          j = 2;
        } else {
          k = 0;
          i = 2; // s2
          j = 1;
        }
        getvrtx(si, i);
        getvrtx(sj, j);
        getvrtx(sk, k);

        if (!hff2(s1, si, sj)) {
          select_1ij();
          projectOnPlane(s1, si, sj, v);
        } else if (!hff2(s1, si, sk)) {
          select_1ik();
          projectOnPlane(s1, si, sk, v);
        } else {
          select_1i();
          projectOnLine(s1, si, v);
        }
      } else if (dotTotal == 2) {
        // Here si is set such that hff(s1,si) < 0
        s->nvrtx = 3;
        if (!testLineThree) {
          k = 2;
          i = 1; // s3
          j = 0;
        } else if (!testLineFour) {
          k = 1;
          i = 0; // s4
          j = 2;
        } else {
          k = 0;
          i = 2; // s2
          j = 1;
        }
        getvrtx(si, i);
        getvrtx(sj, j);
        getvrtx(sk, k);

        if (!hff2(s1, sj, sk)) {
          if (!hff2(s1, sk, sj)) {
            select_1jk(); // select region 1jk
            projectOnPlane(s1, sj, sk, v);
          } else if (!hff2(s1, sk, si)) {
            select_1ik();
            projectOnPlane(s1, sk, si, v);
          } else {
            select_1k();
            projectOnLine(s1, sk, v);
          }
        } else if (!hff2(s1, sj, si)) {
          select_1ij();
          projectOnPlane(s1, si, sj, v);
        } else {
          select_1j();
          projectOnLine(s1, sj, v);
        }
      }
      break;
    default:
      printf("\nERROR:\tunhandled");
  }
}

inline static void
subalgorithm(gkSimplex* s, gkFloat* v) {
  switch (s->nvrtx) {
    case 4:
      S3D(s, v);
      break;
    case 3:
      S2D(s, v);
      break;
    case 2:
      S1D(s, v);
      break;
    default:
      printf("\nERROR:\t invalid simplex\n");
  }
}



double closest_points_GJK_SV2( Convex const& a, 
                       Convex const& b, 
                       Transform const& a2w,
                       Transform const& b2w, 
                       Point3& pa,
                       Point3& pb,
                       int& nbIter )
{
  double relError = GrainsExec::m_colDetTolerance;
  double absError = 1.e-4 * relError;
  bool acceleration = GrainsExec::m_colDetAcceleration;
  double momentum = 0., oneMinusMomentum = 1.;

  unsigned int numIterations = 0;   /**< Iteration counter                 */
  const int maxNumIterations = 30;  /**< Maximum number of GJK iterations  */
  unsigned int i;
  double v[3];
  double mu = 0.;
  

  // Vectors and initial search direction
  Vector3 wVec = a2w( a.support( Vector3Nul ) ) - 
                 b2w( b.support( Vector3Nul ) );
  Vector3 vVec = wVec;
  Vector3 dVec = wVec;
  
  
  /* Initialise simplex */
  double dist = Norm( vVec );
  gkSimplex s = { 1, { 0. } };
  for (int t = 0; t < 3; ++t)
    s.vrtx[0][t] = wVec[t];


  /* Begin GJK iteration */
  do {
    numIterations++;

    // Finding the suitable direction using either Nesterov or original
    // The number 8 is hard-coded. Emprically, it shows the best convergence for
    // superquadrics. For the rest of shapes, we really do not need to use
    // Nesterov as the improvemenet is marginal.
    if ( acceleration && numIterations % 8 != 0 )
    {
      momentum = numIterations / ( numIterations + 2. );
      oneMinusMomentum = 1. - momentum;
      dVec = momentum * dVec + 
             momentum * oneMinusMomentum * vVec +
             oneMinusMomentum * oneMinusMomentum * wVec;
    }
    else
      dVec = vVec;
    
    wVec = a2w( a.support( ( -dVec ) * a2w.getBasis() ) ) - 
           b2w( b.support( (  dVec ) * b2w.getBasis() ) );
    
    // termination criteria
    if ( mu < dist * relError ||
         mu < absError )
    {
      if ( acceleration )
      {
        if ( Norm( dVec - vVec ) < EPSILON )
          break;
        acceleration = false;
        wVec = a2w( a.support( ( -vVec ) * a2w.getBasis() ) ) - 
               b2w( b.support( (  vVec ) * b2w.getBasis() ) );
      }
      else
        break;
    }

    /* Add new vertex to simplex */
    i = s.nvrtx;
    for (int t = 0; t < 3; ++t) {
      s.vrtx[i][t] = wVec[t];
    }
    s.nvrtx++;

    /* Invoke distance sub-algorithm */
    subalgorithm(&s, v);
    dist = sqrt( norm2( v ) );
    vVec.setValue( v[0], v[1], v[2] );

  } while ( ( s.nvrtx != 4 ) && 
            ( numIterations != maxNumIterations ) && 
            ( dist > EPSILON2 ) );

  // It is not the best way to return the witness points.
  // It seems to be inconsistent with the way that it has been done in the
  // original GJK.
  // We can use the latest update on wVec to get a better estimate on witness
  // points.
  pa = a.support( ( -vVec ) * a2w.getBasis() );
  pb = a.support( (  vVec ) * b2w.getBasis() );
  nbIter = numIterations;
  return ( dist );
}