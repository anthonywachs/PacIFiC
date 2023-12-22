#include "run.h"
#include "embed.h"
#include "view.h"
#if p_n > 0
#include "ibm.h"
#endif

/*
==========================================================================================
We use D2Q9 lattice model for 2D simulations and D3Q19 model for 3D simulations. 
Parameter list:
q:             discrete index
opposite[q]:   opposite direction
e[q]:          discrete direction
w[q]:          weight
nu:            kinetic viscosity
omega_0:       relaxation parameter
rho_0:         initial density, set to 1 in general cases
uc:            characteristic velocity
step_max:      total time steps
gravity:       driven-force in fluid
==========================================================================================
*/
#if dimension == 2
#define q (9) 
#else // dimension == 3
#define q (19) 
#endif

//physical properties
#define rho_0   (1.)

// #if uc
// #define uc  (1.) 
// #endif
double uc;
double lref;
double tref;
double mref;

double flow_rate;
double flow_rate_pre = 100;
double flow_rate_set;

bool exceed = false;
bool below  = false;
bool turn   = false;

int cell_tn = 0;
double fvf  = 0.;
coord utot  = {0, 0, 0};
coord utot_pre  = {0, 0, 0};
coord uftot = {0, 0, 0};

double vf;
double x1, x2;
double Dref;

coord gravity;

// lattice model
//	6	2	5
//	3	0	1
//	7	4	8
#if dimension == 2 
const coord e[q] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
const int opposite[q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
const double w[q] = {4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.};
#else // dimension == 3
const coord e[q] = {{0,0,0}, {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}, {1,1,0}, {-1,1,0}, {1,-1,0}, {-1,-1,0}, {1,0,1}, {-1,0,1}, {1,0,-1}, {-1,0,-1}, {0,1,1}, {0,-1,1}, {0,1,-1}, {0,-1,-1}};
const int opposite[q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
const double w[q] = {1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36, 1./36., 1./36, 1./36, 1./36, 1./36, 1./36, 1./36, 1./36, 1./36, 1./36, 1./36};
#endif

//viscosity
double nu;

//relaxation parameter
#define omega_0 (1. / (0.5 + nu*3.))

//maximum iteration number
int step_max;


/*
==========================================================================================
Basilsik-type variables:
f:       distribution function
rho:     density
omega:   viscosity
u:       velocity
force:   force experienced by each fluid cell
cc & ff: VOF functions for post-processing
==========================================================================================
*/
#if dimension == 2
scalar f0[], f1[], f2[], f3[], f4[], f5[], f6[], f7[], f8[];
scalar * f = {f0, f1, f2, f3, f4, f5, f6, f7, f8};
#else // dimension == 3
scalar f0[], f1[], f2[], f3[], f4[], f5[], f6[], f7[], f8[], f9[], f10[], f11[], f12[], f13[], f14[], f15[], f16[], f17[], f18[];
scalar * f = {f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18};
#endif
scalar rho[];
scalar omega[];
vector u[];
#if p_n == 0
vector force[];
#endif

scalar cc[];
face vector ff[];
scalar cc_act[];
face vector ff_act[];

typedef struct
{
    int type;
    coord value;
}direction;

struct
{
    direction left;
    direction right;
    direction top;
    direction bottom;
    direction front;
    direction back;
}bd = {{0,{0,0,0}},{0,{0,0,0}},{0,{0,0,0}},{0,{0,0,0}},{0,{0,0,0}},{0,{0,0,0}}};
/*
Boundary type 0: Periodic
Boundary type 1: Dirichlet
Boundary type 2: Neumann
Boundary type 3: In_Out
*/

//apply boundary conditions
#define is_neighbor(...) (allocated(__VA_ARGS__))

#ifndef right_bd
#define right_bd  (L0 / 2. - Delta / 2.)
#endif

#ifndef left_bd
#define left_bd   (-L0 / 2. + Delta / 2.)
#endif

#ifndef top_bd
#define top_bd    (L0 / 2. - Delta / 2.)
#endif

#ifndef bottom_bd
#define bottom_bd (-L0 / 2. + Delta / 2.)
#endif

#if dimension == 3
#ifndef front_bd
#define front_bd  (L0 / 2. - Delta / 2.)
#endif
#ifndef back_bd
#define back_bd   (-L0 / 2. + Delta / 2.)
#endif
#endif // dimension == 3

#ifndef my_periodic
#define my_periodic 0
#endif

void Calculate_feq(int k, scalar rho, vector u, scalar feq)
{
    foreach()
    {
        double eu = e[k].x * u.x[] + e[k].y * u.y[];
        double u_square = u.x[] * u.x[] + u.y[] * u.y[];
        feq[] = w[k] * rho[] * (1. + (3. * eu + 4.5 * eu * eu - 1.5 * u_square)); 
    }
}

event MRT_Collision_Force(i++)
{
    if (i > 0)
    {
        double S[q];

        for (int i = 0; i < q; i++)
        {
            S[i] = 1.;
        }
        #if dimension == 2
        S[1] = 1.4;
        S[2] = 1.4;
        S[4] = 1.2;
        S[6] = 1.2;
        S[7] = omega_0;
        S[8] = omega_0;
        #else // dimension == 3
        S[1] = 1.19;
        S[2] = 1.4;
        S[4] = 1.2;
        S[6] = S[4];
        S[8] = S[4];
        S[9] = omega_0;
        S[10] = 1.4;
        S[11] = S[9];
        S[12] = S[10];
        S[13] = omega_0;
        S[14] = S[13];
        S[15] = S[13];
        S[16] = 1.98;
        S[17] = S[16];
        S[18] = S[16];
        #endif

        foreach()
        {
            double m_temp[q];
            double meq_temp[q];
            double force_la[q];
            {
                coord f_vel, f_g, f_j;
                foreach_dimension()
                {
                    f_vel.x = u.x[] + 0.5 * (force.x[] + gravity.x) / rho[];
                    f_g.x   = (force.x[] + gravity.x) * rho[];
                    f_j.x   = f_vel.x * rho[];
                }
                #if dimension == 2
                m_temp[0] = f0[] + f1[] + f2[] + f3[] + f4[] + f5[] + f6[] + f7[] + f8[];
                m_temp[1] = -4. * f0[] - f1[] - f2[] - f3[] - f4[] + 2. * (f5[] + f6[] + f7[] + f8[]);
                m_temp[2] = 4. * f0[] - 2. * (f1[] + f2[] + f3[] + f4[]) + f5[] + f6[] + f7[] + f8[];
                m_temp[3] = f1[] - f3[] + f5[] - f6[] - f7[] + f8[];
                m_temp[4] = -2. * (f1[] - f3[]) + f5[] - f6[] - f7[] + f8[];
                m_temp[5] = f2[] - f4[] + f5[] + f6[] - f7[] - f8[];
                m_temp[6] = -2. * (f2[] - f4[]) + f5[] + f6[] - f7[] - f8[];
                m_temp[7] = f1[] - f2[] + f3[] - f4[];
                m_temp[8] = f5[] - f6[] + f7[] - f8[];

                meq_temp[0] = rho[];
                meq_temp[1] = -2. * rho[] + 3. * (sq(f_j.x) + sq(f_j.y));
                meq_temp[2] = rho[] - 3. * (sq(f_j.x) + sq(f_j.y));
                meq_temp[3] = f_j.x;
                meq_temp[4] = -f_j.x;
                meq_temp[5] = f_j.y;
                meq_temp[6] = -f_j.y;
                meq_temp[7] = sq(f_j.x) - sq(f_j.y);
                meq_temp[8] = f_j.x * f_j.y;

                force_la[0] = 0.;
                force_la[1] = 6.0 * (f_vel.x * f_g.x + f_vel.y * f_g.y);
                force_la[2] = -6.0 * (f_vel.x * f_g.x + f_vel.y * f_g.y);
                force_la[3] = f_g.x;
                force_la[4] = -f_g.x;
                force_la[5] = f_g.y;
                force_la[6] = -f_g.y;
                force_la[7] = 2.0 * (f_vel.x * f_g.x - f_vel.y * f_g.y);
                force_la[8] = f_vel.x * f_g.y + f_vel.y * f_g.x;

                for (int k = 0; k < q; ++k)
                {
                    m_temp[k] = m_temp[k] - (m_temp[k] - meq_temp[k]) * S[k] + (1.0 - 0.5 * S[k]) * force_la[k];
                }

                f0[] = (4. * (m_temp[0] - m_temp[1] + m_temp[2])) / 36.;
                f1[] = (4. * m_temp[0] - m_temp[1] - 2. * m_temp[2] + 6. * (m_temp[3] - m_temp[4]) + 9. * m_temp[7]) / 36.;
                f2[] = (4. * m_temp[0] - m_temp[1] - 2. * m_temp[2] + 6. * (m_temp[5] - m_temp[6]) - 9. * m_temp[7]) / 36.;
                f3[] = (4. * m_temp[0] - m_temp[1] - 2. * m_temp[2] - 6. * (m_temp[3] - m_temp[4]) + 9. * m_temp[7]) / 36.;
                f4[] = (4. * m_temp[0] - m_temp[1] - 2. * m_temp[2] - 6. * (m_temp[5] - m_temp[6]) - 9. * m_temp[7]) / 36.;
                f5[] = (4. * m_temp[0] + 2. * m_temp[1] + m_temp[2] + 6. * m_temp[3] + 3. * m_temp[4] + 6. * m_temp[5] + 3. * m_temp[6] + 9. * m_temp[8]) / 36.;
                f6[] = (4. * m_temp[0] + 2. * m_temp[1] + m_temp[2] - 6. * m_temp[3] - 3. * m_temp[4] + 6. * m_temp[5] + 3. * m_temp[6] - 9. * m_temp[8]) / 36.;
                f7[] = (4. * m_temp[0] + 2. * m_temp[1] + m_temp[2] - 6. * m_temp[3] - 3. * m_temp[4] - 6. * m_temp[5] - 3. * m_temp[6] + 9. * m_temp[8]) / 36.;
                f8[] = (4. * m_temp[0] + 2. * m_temp[1] + m_temp[2] + 6. * m_temp[3] + 3. * m_temp[4] - 6. * m_temp[5] - 3. * m_temp[6] - 9. * m_temp[8]) / 36.;
            
                #else // dimension == 3
                m_temp[0] = f0[] + f1[] + f2[] + f3[] + f4[] + f5[] + f6[] + f7[] + f8[] + f9[] + f10[] + f11[] + f12[] + f13[] + f14[] + f15[] + f16[] + f17[] + f18[];
                m_temp[1] = -30. * f0[] - 11. * (f1[] + f2[] + f3[] + f4[] + f5[] + f6[]) + 8. * (f7[] + f8[] + f9[] + f10[] + f11[] + f12[] + f13[] + f14[] + f15[] + f16[] + f17[] + f18[]);
                m_temp[2] = 12. * f0[] - 4. * (f1[] + f2[] + f3[] + f4[] + f5[] + f6[]) + (f7[] + f8[] + f9[] + f10[] + f11[] + f12[] + f13[] + f14[] + f15[] + f16[] + f17[] + f18[]);
                m_temp[3] = f1[] - f2[] + f7[] - f8[] + f9[] - f10[] + f11[] - f12[] + f13[] - f14[];
                m_temp[4] = -4. * f1[] + 4. * f2[] + f7[] - f8[] + f9[] - f10[] + f11[] - f12[] + f13[] - f14[];
                m_temp[5] = f3[] - f4[] + f7[] + f8[] - f9[] - f10[] + f15[] - f16[] + f17[] - f18[];
                m_temp[6] = -4. * f3[] + 4. * f4[] + f7[] + f8[] - f9[] - f10[] + f15[] - f16[] + f17[] - f18[];
                m_temp[7] = f5[] - f6[] + f11[] + f12[] - f13[] - f14[] + f15[] + f16[] - f17[] - f18[];
                m_temp[8] = -4. * f5[] + 4. * f6[] + f11[] + f12[] - f13[] - f14[] + f15[] + f16[] - f17[] - f18[];
                m_temp[9] = 2. * (f1[] + f2[]) - (f3[] + f4[] + f5[] + f6[]) + f7[] + f8[] + f9[] + f10[] + f11[] + f12[] + f13[] + f14[] - 2. * (f15[] + f16[] + f17[] + f18[]);
                m_temp[10] = -4. * (f1[] + f2[]) + 2. * (f3[] + f4[] + f5[] + f6[]) + f7[] + f8[] + f9[] + f10[] + f11[] + f12[] + f13[] + f14[] - 2. * (f15[] + f16[] + f17[] + f18[]);
                m_temp[11] = (f3[] + f4[]) - (f5[] + f6[]) + f7[] + f8[] + f9[] + f10[] - (f11[] + f12[] + f13[] + f14[]);
                m_temp[12] = -2 * (f3[] + f4[]) + 2. * (f5[] + f6[]) + f7[] + f8[] + f9[] + f10[] - (f11[] + f12[] + f13[] + f14[]);
                m_temp[13] = f7[] - f8[] - f9[] + f10[];
                m_temp[14] = f15[] - f16[] - f17[] + f18[];
                m_temp[15] = f11[] - f12[] - f13[] + f14[];
                m_temp[16] = f7[] - f8[] + f9[] - f10[] - f11[] + f12[] - f13[] + f14[];
                m_temp[17] = -f7[] - f8[] + f9[] + f10[] + f15[] - f16[] + f17[] - f18[];
                m_temp[18] = f11[] + f12[] - f13[] - f14[] - f15[] - f16[] + f17[] + f18[];

                meq_temp[0] = rho[];
                meq_temp[1] = -11. * rho[] + 19. * (sq(f_j.x) + sq(f_j.y) + sq(f_j.z));
                meq_temp[2] = 3. * rho[] - 11 * (sq(f_j.x) + sq(f_j.y) + sq(f_j.z)) / 2.;
                meq_temp[3] = f_j.x;
                meq_temp[4] = -2. * f_j.x / 3.;
                meq_temp[5] = f_j.y;
                meq_temp[6] = -2. * f_j.y / 3.;
                meq_temp[7] = f_j.z;
                meq_temp[8] = -2. * f_j.z / 3.;
                meq_temp[9] = 2. * sq(f_j.x) - (sq(f_j.y) + sq(f_j.z));
                meq_temp[10] = -0.5 * (2. * sq(f_j.x) - (sq(f_j.y) + sq(f_j.z)));
                meq_temp[11] = sq(f_j.y) - sq(f_j.z);
                meq_temp[12] = -0.5 * (sq(f_j.y) - sq(f_j.z));
                meq_temp[13] = f_j.x * f_j.y;
                meq_temp[14] = f_j.y * f_j.z;
                meq_temp[15] = f_j.x * f_j.z;
                meq_temp[16] = 0.;
                meq_temp[17] = 0.;
                meq_temp[18] = 0.;

                force_la[0] = 0.;
                force_la[1] = 38.0 * (f_vel.x * f_g.x + f_vel.y * f_g.y + f_vel.z * f_g.z);
                force_la[2] = -11.0 * (f_vel.x * f_g.x + f_vel.y * f_g.y + f_vel.z * f_g.z);
                force_la[3] = f_g.x;
                force_la[4] = -2. * f_g.x / 3.;
                force_la[5] = f_g.y;
                force_la[6] = -2. * f_g.y / 3.;
                force_la[7] = f_g.z;
                force_la[8] = -2. * f_g.z / 3.;
                force_la[9] = 2. * (2. * f_vel.x * f_g.x - f_vel.y * f_g.y - f_vel.z * f_g.z);
                force_la[10] = -(2. * f_vel.x * f_g.x - f_vel.y * f_g.y - f_vel.z * f_g.z);
                force_la[11] = 2. * (f_vel.y * f_g.y - f_vel.z * f_g.z);
                force_la[12] = -(f_vel.y * f_g.y - f_vel.z * f_g.z);
                force_la[13] = f_vel.x * f_g.y + f_vel.y * f_g.x;
                force_la[14] = f_vel.y * f_g.z + f_vel.z * f_g.y;
                force_la[15] = f_vel.z * f_g.x + f_vel.x * f_g.z;
                force_la[16] = 0.;
                force_la[17] = 0.;
                force_la[18] = 0.;

                for (int k = 0; k < q; ++k)
                {
                    m_temp[k] = m_temp[k] - (m_temp[k] - meq_temp[k]) * S[k] + (1.0 - 0.5 * S[k]) * force_la[k];
                }

                f0[] = m_temp[0] / 19. - 5. * m_temp[1] / 399. + m_temp[2] / 21.;
                f1[] = m_temp[0] / 19. - 11. * m_temp[1] / 2394. - m_temp[2] / 63. + (m_temp[3] - m_temp[4]) / 10. + (m_temp[9] - m_temp[10]) / 18.;
                f2[] = m_temp[0] / 19. - 11. * m_temp[1] / 2394. - m_temp[2] / 63. - (m_temp[3] - m_temp[4]) / 10. + (m_temp[9] - m_temp[10]) / 18.;
                f3[] = m_temp[0] / 19. - 11. * m_temp[1] / 2394. - m_temp[2] / 63. + (m_temp[5] - m_temp[6]) / 10. - (m_temp[9] - m_temp[10]) / 36. + (m_temp[11] - m_temp[12]) / 12.;
                f4[] = m_temp[0] / 19. - 11. * m_temp[1] / 2394. - m_temp[2] / 63. - (m_temp[5] - m_temp[6]) / 10. - (m_temp[9] - m_temp[10]) / 36. + (m_temp[11] - m_temp[12]) / 12.;
                f5[] = m_temp[0] / 19. - 11. * m_temp[1] / 2394. - m_temp[2] / 63. + (m_temp[7] - m_temp[8]) / 10. - (m_temp[9] - m_temp[10]) / 36. - (m_temp[11] - m_temp[12]) / 12.;
                f6[] = m_temp[0] / 19. - 11. * m_temp[1] / 2394. - m_temp[2] / 63. - (m_temp[7] - m_temp[8]) / 10. - (m_temp[9] - m_temp[10]) / 36. - (m_temp[11] - m_temp[12]) / 12.;

                f7[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (m_temp[3] + m_temp[5]) / 10. + (m_temp[4] + m_temp[6]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. + m_temp[11] / 12. + m_temp[12] / 24. + m_temp[13] / 4. + (m_temp[16] - m_temp[17]) / 8.;
                f8[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (-m_temp[3] + m_temp[5]) / 10. + (-m_temp[4] + m_temp[6]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. + m_temp[11] / 12. + m_temp[12] / 24. - m_temp[13] / 4. + (-m_temp[16] - m_temp[17]) / 8.;
                f9[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (m_temp[3] - m_temp[5]) / 10. + (m_temp[4] - m_temp[6]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. + m_temp[11] / 12. + m_temp[12] / 24. - m_temp[13] / 4. + (m_temp[16] + m_temp[17]) / 8.;
                f10[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. - (m_temp[3] + m_temp[5]) / 10. - (m_temp[4] + m_temp[6]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. + m_temp[11] / 12. + m_temp[12] / 24. + m_temp[13] / 4. + (-m_temp[16] + m_temp[17]) / 8.;

                f11[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (m_temp[3] + m_temp[7]) / 10. + (m_temp[4] + m_temp[8]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. - m_temp[11] / 12. - m_temp[12] / 24. + m_temp[15] / 4. - (m_temp[16] - m_temp[18]) / 8.;
                f12[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (-m_temp[3] + m_temp[7]) / 10. + (-m_temp[4] + m_temp[8]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. - m_temp[11] / 12. - m_temp[12] / 24. - m_temp[15] / 4. + (m_temp[16] + m_temp[18]) / 8.;
                f13[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (m_temp[3] - m_temp[7]) / 10. + (m_temp[4] - m_temp[8]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. - m_temp[11] / 12. - m_temp[12] / 24. - m_temp[15] / 4. - (m_temp[16] + m_temp[18]) / 8.;
                f14[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. - (m_temp[3] + m_temp[7]) / 10. - (m_temp[4] + m_temp[8]) / 40. + m_temp[9] / 36. + m_temp[10] / 72. - m_temp[11] / 12. - m_temp[12] / 24. + m_temp[15] / 4. + (m_temp[16] - m_temp[18]) / 8.;

                f15[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (m_temp[5] + m_temp[7]) / 10. + (m_temp[6] + m_temp[8]) / 40. - m_temp[9] / 18. - m_temp[10] / 36. + m_temp[14] / 4 + (m_temp[17] - m_temp[18]) / 8.;
                f16[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (-m_temp[5] + m_temp[7]) / 10. + (-m_temp[6] + m_temp[8]) / 40. - m_temp[9] / 18. - m_temp[10] / 36. - m_temp[14] / 4 + (-m_temp[17] - m_temp[18]) / 8.;
                f17[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. + (m_temp[5] - m_temp[7]) / 10. + (m_temp[6] - m_temp[8]) / 40. - m_temp[9] / 18. - m_temp[10] / 36. - m_temp[14] / 4 + (m_temp[17] + m_temp[18]) / 8.;
                f18[] = m_temp[0] / 19. + 4. * m_temp[1] / 1197. + m_temp[2] / 252. - (m_temp[5] + m_temp[7]) / 10. - (m_temp[6] + m_temp[8]) / 40. - m_temp[9] / 18. - m_temp[10] / 36. + m_temp[14] / 4 + (-m_temp[17] + m_temp[18]) / 8.;
                #endif
            }
        }
    }
}

#if motion_type == 3
event Domain_Boundary(i++)
{
    if (i > 0)
    {
        foreach()
        {
/*
==========================================================================================
Inlet & Outlet: Dirichlet velocity BC at the inlet and Neumann pressure(density in LBM) at
the outlet. Here we only consider in/outlet along the same axis, and flow is always headi-
ng to the positive direction. The equilibrium scheme is adopted for the in/outlet.
==========================================================================================
*/
            //X-axis
            //Inlet on the left.
            if (bd.left.type == 3)
            {
                if (x == left_bd)
                {
                    //Velocity field is given by the Dirichlet BCs.
                    u.x[] = 4.0 * uc * ((y - bottom_bd) * (top_bd - bottom_bd) - sq(y - bottom_bd)) / sq(top_bd - bottom_bd);
                    u.y[] = 0.;
                    #if dimension == 3
                    u.z[] = 0.;
                    #endif // dimension == 3
                    //Periodic boundary condition is considered for the fully developed flow.
                    rho[] = rho[-2,0];
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_o += e[i].x * u.x[];
                        }
                        //The distributions are estimated by their equilibrium counterparts.
                        f_temp[] = w[i] * rho[-2,0] * (1 + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        i++;
                    }
                }

                // if (x == left_bd)
                // {
                //     //Velocity field is given by the Dirichlet BCs.
                //     u.x[] = 4.0 * uc * ((y - bottom_bd) * (top_bd - bottom_bd) - sq(y - bottom_bd)) / sq(top_bd - bottom_bd);
                //     u.y[] = 0.;
                //     #if dimension == 3
                //     u.z[] = 0.;
                //     #endif // dimension == 3
                //     //Periodic boundary condition is considered for the fully developed flow.
                //     rho[] = rho[-2,0];
                //     double u_square_o = 0.;
                //     double u_square_i = 0.;
                //     foreach_dimension()
                //     {
                //         u_square_i += u.x[-2,0] * u.x[-2,0];
                //         u_square_o += u.x[] * u.x[];
                //     }
                //     int i = 0;
                //     for (scalar f_temp in f)
                //     {
                //         double eu_i = 0.;
                //         double eu_o = 0.;
                //         foreach_dimension()
                //         {
                //             eu_i += e[i].x * u.x[-2,0];
                //             eu_o += e[i].x * u.x[];
                //         }
                //         //The distributions are estimated by their equilibrium counterparts.
                //         double feq_i = w[i] * rho[-2,0] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                //         double feq_o = w[i] * rho[-2,0] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                //         f_temp[] = f_temp[-2,0] - feq_i + feq_o;
                //         i++;
                //     }
                // }
            } 

            //Outlet on the right.
            if (bd.right.type == 3)
            {
                if (x == right_bd)
                {
                    foreach_dimension()
                        u.x[] = u.x[-1,0];
                    rho[] = rho_0;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_o += e[i].x * u.x[];
                        }
                        f_temp[] = w[i] * rho_0 * (1 + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        i++;
                    }
                }

                // if (x == right_bd)
                // {
                //     //Velocity field is given by the Dirichlet BCs.
                //     foreach_dimension()
                //         u.x[] = u.x[-1,0];
                //     //Periodic boundary condition is considered for the fully developed flow.
                //     rho[] = rho_0;
                //     double u_square_o = 0.;
                //     double u_square_i = 0.;
                //     foreach_dimension()
                //     {
                //         u_square_i += u.x[-1,0] * u.x[-1,0];
                //         u_square_o += u.x[] * u.x[];
                //     }
                //     int i = 0;
                //     for (scalar f_temp in f)
                //     {
                //         double eu_i = 0.;
                //         double eu_o = 0.;
                //         foreach_dimension()
                //         {
                //             eu_i += e[i].x * u.x[-1,0];
                //             eu_o += e[i].x * u.x[];
                //         }
                //         //The distributions are estimated by their equilibrium counterparts.
                //         double feq_i = w[i] * rho[-1,0] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                //         double feq_o = w[i] * rho[] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                //         f_temp[] = f_temp[-1,0] - feq_i + feq_o;
                //         i++;
                //     }
                // }
            }

            //Y-axis
            //Inlet on the bottom.
            if (bd.bottom.type == 3)
            {
                if (y == bottom_bd)
                {
                    u.x[] = 0.;
                    u.y[] = 4.0 * uc * ((x - left_bd) * (right_bd - left_bd) - sq(x - left_bd)) / sq(right_bd - left_bd);
                    #if dimension == 3
                    u.z[] = 0.;
                    #endif
                    rho[] = rho[0,-2];
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_o += e[i].x * u.x[];
                        }
                        f_temp[] = w[i] * rho[0,-2] * (1 + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        i++;
                    }
                }
            }
            
            //Outlet on the top.
            if (bd.top.type == 3) 
            {
                if (y == top_bd)
                {
                    u.x[] = u.x[0,-1];
                    u.y[] = u.y[0,-1];
                    #if dimension == 3
                    u.z[] = u.z[0,-1];
                    #endif
                    rho[] = rho_0;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_o += e[i].x * u.x[];
                        }
                        f_temp[] = w[i] * rho_0 * (1 + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        i++;
                    }
                }
            }

            //Z-axis
            #if dimension == 3
            //Inlet on the back.
            if (bd.back.type == 3)
            {
                if (z == back_bd)
                {
                    u.x[] = 0.;
                    u.y[] = 0.;
                    #if dimension == 3
                    u.z[] = uc;
                    #endif
                    rho[] = rho[0,0,-2];
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_o += e[i].x * u.x[];
                        }
                        f_temp[] = w[i] * rho[0,0,-2] * (1 + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        i++;
                    }
                }
            }
            
            //Outlet on the front.
            if (bd.front.type == 3) 
            {
                if (y == top_bd)
                {
                    u.x[] = u.x[0,0,-1];
                    u.y[] = u.y[0,0,-1];
                    #if dimension == 3
                    u.z[] = u.z[0,0,-1];
                    #endif
                    rho[] = rho_0;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_o += e[i].x * u.x[];
                        }
                        f_temp[] = w[i] * rho_0 * (1 + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        i++;
                    }
                }
            }
            #endif

/*
==========================================================================================
The non-equilibrium extrapolation scheme is adopted for the Dirichlet BCs (for velocity).
==========================================================================================
*/
            //X-axis
            if(bd.left.type == 1)
            {
                if (x == left_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.left.value.x;
                    rho[] = rho[1];
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[1] * u.x[1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho[1] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho[1] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            if(bd.right.type == 1)
            {
                if (x == right_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.right.value.x;
                    rho[] = rho[-1];
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[-1] * u.x[-1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[-1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho[-1] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho[-1] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[-1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            //Y-axis
            if(bd.top.type == 1)
            {
                if (y == top_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.top.value.x;
                    rho[] = rho[0,-1];
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,-1] * u.x[0,-1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,-1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho[0,-1] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho[0,-1] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,-1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            if(bd.bottom.type == 1)
            {
                if (y == bottom_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.bottom.value.x;
                    rho[] = rho[0,1];
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,1] * u.x[0,1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho[0,1] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho[0,1] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            //Z-axis
            #if dimension == 3
            if(bd.front.type == 1)
            {
                if (z == front_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.front.value.x;
                    rho[] = rho[0,0,-1];
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,0,-1] * u.x[0,0,-1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,0,-1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho[0,0,-1] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho[0,0,-1] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,0,-1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            if(bd.back.type == 1)
            {
                if (z == back_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.back.value.x;
                    rho[] = rho[0,0,1];
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,0,1] * u.x[0,0,1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,0,1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho[0,0,1] * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho[0,0,1] * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,0,1] - feq_i + feq_o;
                        i++;
                    }
                }
            }
            #endif
/*
==========================================================================================
The non-equilibrium extrapolation scheme is adopted for the Neumann BCs (for velocity).
==========================================================================================
*/
            //X-axis
            if(bd.left.type == 2)
            {
                if (x == left_bd)
                {
                    foreach_dimension()
                        u.x[] = u.x[1] - bd.left.value.x*Delta;
                    rho[] = rho_0;
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[1] * u.x[1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho_0 * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho_0 * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            if(bd.right.type == 2)
            {
                if (x == right_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.right.value.x*Delta + u.x[-1];
                    rho[] = rho_0;
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[-1] * u.x[-1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[-1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho_0 * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho_0 * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[-1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            //Y-axis
            if(bd.top.type == 2)
            {
                if (y == top_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.bottom.value.x*Delta + u.x[0,-1];
                    rho[] = rho_0;
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,-1] * u.x[0,-1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,-1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho_0 * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho_0 * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,-1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            if(bd.bottom.type == 2)
            {
                if (y == bottom_bd)
                {
                    foreach_dimension()
                        u.x[] = u.x[0,1] - bd.top.value.x*Delta;
                    rho[] = rho_0;
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,1] * u.x[0,1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho_0 * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho_0 * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            //Z-axis
            #if dimension == 3
            if(bd.front.type == 2)
            {
                if (z == front_bd)
                {
                    foreach_dimension()
                        u.x[] = bd.back.value.x*Delta + u.x[0,0,-1];
                    rho[] = rho_0;
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,0,-1] * u.x[0,0,-1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,0,-1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho_0 * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho_0 * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,0,-1] - feq_i + feq_o;
                        i++;
                    }
                }
            }

            if(bd.back.type == 2)
            {
                if (z == back_bd)
                {
                    foreach_dimension()
                        u.x[] = u.x[0,0,1] - bd.front.value.x*Delta;
                    rho[] = rho_0;
                    double u_square_i = 0.;
                    double u_square_o = 0.;
                    foreach_dimension()
                    {
                        u_square_i += u.x[0,0,1] * u.x[0,0,1];
                        u_square_o += u.x[] * u.x[];
                    }
                    int i = 0;
                    for (scalar f_temp in f)
                    {
                        double eu_i = 0.;
                        double eu_o = 0.;
                        foreach_dimension()
                        {
                            eu_i += e[i].x * u.x[0,0,1];
                            eu_o += e[i].x * u.x[];
                        }
                        double feq_i = w[i] * rho_0 * (1. + (3. * eu_i + 4.5 * eu_i * eu_i - 1.5 * u_square_i));
                        double feq_o = w[i] * rho_0 * (1. + (3. * eu_o + 4.5 * eu_o * eu_o - 1.5 * u_square_o));
                        f_temp[] = f_temp[0,0,1] - feq_i + feq_o;
                        i++;
                    }
                }
            }
            #endif
        }
    }
}
#endif // motion_type

scalar f_new[];
//streaming operator (Lax-Wendroff scheme)
event Streaming_LaxWendroff(i++)
{
    if (i > 0)
    {
        int i = 0;
        for (scalar f_old in f)
        {
            foreach()
            {
                #if my_periodic == 0
                #if dimension == 2
                if(x != left_bd && x != right_bd && y != top_bd && y != bottom_bd)
                #elif dimension == 3 
                if(x != left_bd && x != right_bd && y != top_bd && y != bottom_bd && z != front_bd && z != back_bd)
                #endif // dimension
                #elif my_periodic == 1
                #if dimension == 2
                if(y != top_bd && y != bottom_bd)
                #elif dimension == 3 
                if(y != top_bd && y != bottom_bd && z != front_bd && z != back_bd)
                #endif // dimension
                #endif // my_periodic
                {
                    double cfl = 1. / Delta;
                    if (cfl <= 0.99)
                    {
                        #if dimension == 2
                        f_new[] = cfl / 2. * (1. + cfl) * f_old[-(int)e[i].x, -(int)e[i].y] + (1. - cfl * cfl) * f_old[] - cfl / 2. * (1. - cfl) * f_old[(int)e[i].x, (int)e[i].y];
                        #else // dimension == 3
                        f_new[] = cfl / 2. * (1. + cfl) * f_old[-(int)e[i].x, -(int)e[i].y, -(int)e[i].z] + (1. - cfl * cfl) * f_old[] - cfl / 2. * (1. - cfl) * f_old[(int)e[i].x, (int)e[i].y, (int)e[i].z];
                        #endif
                        // f_new[] = cfl / 2. * (1. + cfl) * f_old[-ex[i], -ey[i], -ez[i]] + (1. - cfl * cfl) * f_old[] - cfl / 2. * (1. - cfl) * f_old[ex[i], ey[i], ez[i]];
                    }
                    else
                    {
                        #if dimension == 2
                        f_new[] = f_old[-(int)e[i].x, -(int)e[i].y];
                        #else // dimension == 3
                        f_new[] = f_old[-(int)e[i].x, -(int)e[i].y, -(int)e[i].z];
                        #endif
                    }
                }
            }

            foreach()
            {
                #if my_periodic == 0
                #if dimension == 2
                if(x != left_bd && x != right_bd && y != top_bd && y != bottom_bd)
                #elif dimension == 3 
                if(x != left_bd && x != right_bd && y != top_bd && y != bottom_bd && z != front_bd && z != back_bd)
                #endif // dimension
                #elif my_periodic == 1
                #if dimension == 2
                if(y != top_bd && y != bottom_bd)
                #elif dimension == 3 
                if(y != top_bd && y != bottom_bd && z != front_bd && z != back_bd)
                #endif // dimension
                #endif // my_periodic
                {
                    f_old[] = f_new[];
                }
            }

            i++;
        }
    }
}


//evaluations of the macro-scopic properties
event Macro_Properties(i++)
{
    if (i > 0)
    {
        foreach()
        {
            #if my_periodic == 0
            #if dimension == 2
            if(x != left_bd && x != right_bd && y != top_bd && y != bottom_bd)
            #elif dimension == 3 
            if(x != left_bd && x != right_bd && y != top_bd && y != bottom_bd && z != front_bd && z != back_bd)
            #endif // dimension
            #elif my_periodic == 1
                #if dimension == 2
                if(y != top_bd && y != bottom_bd)
                #elif dimension == 3 
                if(y != top_bd && y != bottom_bd && z != front_bd && z != back_bd)
                #endif // dimension
            #endif // my_periodic
            {
                rho[] = 0.;
                foreach_dimension()
                    u.x[] = 0.;
                int i = 0;
                for (scalar f_temp in f)
                {
                    rho[] += f_temp[];
                    foreach_dimension()
                        u.x[] += e[i].x * f_temp[];
                    i++;
                }
                foreach_dimension()
                    u.x[] /= rho[];
            }
        }
    }
}

