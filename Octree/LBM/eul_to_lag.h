#include "ibm.h"

#if _MPI
double vel_x[p_n-1][tn_node];
double vel_y[p_n-1][tn_node];
double vel_x_sum[p_n-1][tn_node];
double vel_y_sum[p_n-1][tn_node];
#if dimension == 3
double vel_z[p_n-1][tn_node];
double vel_z_sum[p_n-1][tn_node];
#endif
double lvel_x[(int)(r_ratio*r_ratio)*tn_node];
double lvel_y[(int)(r_ratio*r_ratio)*tn_node];
double lvel_x_sum[(int)(r_ratio*r_ratio)*tn_node];
double lvel_y_sum[(int)(r_ratio*r_ratio)*tn_node];
#if dimension == 3
double lvel_z[(int)(r_ratio*r_ratio)*tn_node];
double lvel_z_sum[(int)(r_ratio*r_ratio)*tn_node];
#endif
#endif

event interpolation(i++)
{
    if (i > 0)
    {
        for (int k = 0; k < p_n; k++)
        {
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                foreach_dimension()
                    particle[k].node[n].vel.x = 0;
                #if _MPI
                if (k == (p_n-1))
                {
                    lvel_x[n] = 0.;
                    lvel_y[n] = 0.;
                    lvel_x_sum[n] = 0.;
                    lvel_y_sum[n] = 0.;
                    #if dimension == 3
                    lvel_z[n] = 0.;
                    lvel_z_sum[n] = 0.;
                    #endif
                }
                else
                {
                    vel_x[k][n] = 0.;
                    vel_y[k][n] = 0.;
                    vel_x_sum[k][n] = 0.;
                    vel_y_sum[k][n] = 0.;
                    #if dimension == 3
                    vel_z[k][n] = 0.;
                    vel_z_sum[k][n] = 0.;
                    #endif
                }
                #endif
            }
        }

        for (int k = 0; k < p_n; k++)
        {
            // reset node velocity
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                foreach_cache(particle[k].node[n].my_stencil)
                {
                    // compute distance between object node and fluid lattice node.
                    const double dist_x = particle[k].node[n].pos.x - x;
                    const double dist_y = particle[k].node[n].pos.y - y;
                    #if dimension == 3
                    const double dist_z = particle[k].node[n].pos.z - z;
                    #endif

                    #if dimension == 2
                    #if (IBM_stencil == 1 || IBM_stencil == 11)
                    if (fabs(dist_x) <= 1.5 && fabs(dist_y) <= 1.5)
                    #else // (IBM_stencil == 4 || IBM_stencil == 14)
                    if (fabs(dist_x) <= 2.5 && fabs(dist_y) <= 2.5)
                    #endif
                    #elif dimension == 3
                    #if (IBM_stencil == 1 || IBM_stencil == 11)
                    if (fabs(dist_x) <= 1.5 && fabs(dist_y) <= 1.5 && fabs(dist_z) <= 1.5)
                    #else // (IBM_stencil == 4 || IBM_stencil == 14)
                    if (fabs(dist_x) <= 2.5 && fabs(dist_y) <= 2.5 && fabs(dist_z) <= 2.5)
                    #endif
                    #endif
                    {
                        // compute interpolation weights for x- and y-direction based on the distance.
                        const double weight_x = stencil(dist_x);
                        const double weight_y = stencil(dist_y);
                        #if dimension == 3
                        const double weight_z = stencil(dist_z);
                        #endif

                        // compute node velocity.
                        #if _MPI
                        if (k == (p_n-1))
                        {
                            #if dimension == 2
                            lvel_x[n] += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y;
                            lvel_y[n] += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y;
                            #elif dimension == 3
                            lvel_x[n] += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y * weight_z;
                            lvel_y[n] += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y * weight_z;
                            lvel_z[n] += (u.z[] + 0.5 * force.z[] / rho[]) * weight_x * weight_y * weight_z;
                            #endif
                        }
                        else
                        {
                            #if dimension == 2
                            vel_x[k][n] += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y;
                            vel_y[k][n] += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y;
                            #elif dimension == 3
                            vel_x[k][n] += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y * weight_z;
                            vel_y[k][n] += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y * weight_z;
                            vel_z[k][n] += (u.z[] + 0.5 * force.z[] / rho[]) * weight_x * weight_y * weight_z;
                            #endif
                        }
                        #else
                        #if dimension == 2
                        particle[k].node[n].vel.x += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y;
                        particle[k].node[n].vel.y += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y;
                        #elif dimension == 3
                        particle[k].node[n].vel.x += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y * weight_z;
                        particle[k].node[n].vel.y += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y * weight_z;
                        particle[k].node[n].vel.z += (u.z[] + 0.5 * force.z[] / rho[]) * weight_x * weight_y * weight_z;
                        #endif
                        #endif
                    }
                }
            }
        }
    }
}

event update_particle_position(i++)
{
    if (i > 0)
    {
        for (int k = 0; k < p_n; k++)
        {
            // foreach_dimension()
            //     particle[k].center.pos.x = 0.;

            #if _MPI
            if (k == (p_n-1))
            {
                MPI_Reduce(lvel_x, lvel_x_sum, particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Bcast(lvel_x_sum, particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Reduce(lvel_y, lvel_y_sum, particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Bcast(lvel_y_sum, particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Reduce(vel_x[k], vel_x_sum[k], particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Bcast(vel_x_sum[k], particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Reduce(vel_y[k], vel_y_sum[k], particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Bcast(vel_y_sum[k], particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            #if dimension == 3
            if (k == (p_n-1))
            {
                MPI_Reduce(lvel_z, lvel_z_sum, particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Bcast(lvel_z_sum, particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Reduce(vel_z[k], vel_z_sum[k], particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Bcast(vel_z_sum[k], particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            #endif
            #endif

            // update node and center positions
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                #if _MPI
                if (k == (p_n-1))
                {
                    particle[k].node[n].pos.x += lvel_x_sum[n];
                    particle[k].node[n].pos.y += lvel_y_sum[n];
                    #if dimension == 3
                    particle[k].node[n].pos.z += lvel_z_sum[n];
                    #endif
                    // particle[k].center.pos.x += lvel_x_sum[n] / particle[k].num_nodes;
                    // particle[k].center.pos.y += lvel_y_sum[n] / particle[k].num_nodes;
                    // particle[k].center.pos.z += lvel_z_sum[n] / particle[k].num_nodes;
                }
                else
                {
                    particle[k].node[n].pos.x += vel_x_sum[k][n];
                    particle[k].node[n].pos.y += vel_y_sum[k][n];
                    #if dimension == 3
                    particle[k].node[n].pos.z += vel_z_sum[k][n];
                    #endif
                    particle[k].center.pos.x += vel_x_sum[k][n] / particle[k].num_nodes;
                    particle[k].center.pos.y += vel_y_sum[k][n] / particle[k].num_nodes;
                    particle[k].center.pos.z += vel_z_sum[k][n] / particle[k].num_nodes;
                }
                #else
                foreach_dimension()
                    particle[k].node[n].pos.x += particle[k].node[n].vel.x;
                #endif
                foreach_dimension()
                {
                    if (particle[k].node[n].pos.x < -L0 / 2)
                        particle[k].node[n].pos.x += L0;
                    else if (particle[k].node[n].pos.x > L0 / 2)
                        particle[k].node[n].pos.x -= L0;
                }
            }

            foreach_dimension()
            {
                if (particle[k].center.pos.x < -L0 / 2)
                    particle[k].center.pos.x += L0;
                else if (particle[k].center.pos.x > L0 / 2)
                    particle[k].center.pos.x -= L0;
            }

            if (k < (p_n-1))
            {
                foreach_dimension()
                {
                    particle[k].F_tot.x = 0.;
                    particle[k].T_tot.z = 0.;
                }

                for (int n = 0; n < particle[k].num_nodes; n++)
                {
                    foreach_dimension()
                    {
                        particle[k].F_tot.x += particle[k].node[n].lag_F.x;
                        particle[k].T_tot.z += particle[k].node[n].lag_T.z;
                    }
                }

                foreach_dimension()
                {
                    // particle[k].center.vel.x = particle[k].center.vel_pre.x + particle[k].center.acl.x / rho_ratio - particle[k].F_tot.x / (particle_volume * rho_ratio) + (1. - 1. / rho_ratio) * fa.x;
                    particle[k].center.vel.x = particle[k].center.vel_pre.x - particle[k].F_tot.x / (particle_volume * rho_ratio) + (1. - 1. / rho_ratio) * fa.x;
		            particle[k].center.acl.x = particle[k].center.vel.x - particle[k].center.vel_pre.x;
                    particle[k].center.pos_ref.x += particle[k].center.vel_pre.x + particle[k].center.acl.x / 2.;
                    particle[k].center.vel_pre.x = particle[k].center.vel.x;

                    if (particle[k].center.pos_ref.x < -L0 / 2)
                        particle[k].center.pos_ref.x += L0;
                    else if (particle[k].center.pos_ref.x > L0 / 2)
                        particle[k].center.pos_ref.x -= L0;

                    // particle[k].center.agl_vel.z = particle[k].center.agl_vel_pre.z + particle[k].center.agl_acl.z / rho_ratio 
                    //                           - particle[k].T_tot.z / (moment_inertia * rho_ratio);
		            particle[k].center.agl_vel.z = particle[k].center.agl_vel_pre.z - particle[k].T_tot.z / (moment_inertia * rho_ratio);
                    particle[k].center.agl_acl.z = particle[k].center.agl_vel.z - particle[k].center.agl_vel_pre.z;
                    particle[k].angle.z += particle[k].center.agl_vel_pre.z + particle[k].center.agl_acl.z / 2.;
                    particle[k].center.agl_vel_pre.z = particle[k].center.agl_vel.z;
                }

                #if dimension == 2
                for (int n = 0; n < particle[k].num_nodes; n++)
                {
                    double xref = particle[k].radius * cos(2. * M_PI * (double)n / particle[k].num_nodes);
                    double yref = particle[k].radius * sin(2. * M_PI * (double)n / particle[k].num_nodes);
                    particle[k].node[n].pos_ref.x = particle[k].center.pos_ref.x + xref * cos(particle[k].angle.z) - yref * sin(particle[k].angle.z);
                    particle[k].node[n].pos_ref.y = particle[k].center.pos_ref.y + xref * sin(particle[k].angle.z) + yref * cos(particle[k].angle.z);
                }
                #elif dimension == 3
                for (int n = 0; n < particle[k].num_nodes; n++)
                {
                    double indice = (double)n + 0.5;
                    double phi = acos(1. - 2. * indice / particle[k].num_nodes);
                    double theta = M_PI * (1. + sqrt(5.)) * indice;
                    double xref = particle[k].radius * cos(theta) * sin(phi);
                    double yref = particle[k].radius * sin(theta) * sin(phi);
                    double zref = particle[k].radius * cos(phi);
                    particle[k].node[n].pos_ref.x = particle[k].center.pos_ref.x + xref * cos(particle[k].angle.z) * cos(particle[k].angle.y) + yref * (cos(particle[k].angle.z) * sin(particle[k].angle.y) * sin(particle[k].angle.x) - sin(particle[k].angle.z) * cos(particle[k].angle.x)) + zref * (cos(particle[k].angle.z) * sin(particle[k].angle.y) * cos(particle[k].angle.x) + sin(particle[k].angle.z) * sin(particle[k].angle.x));
                    particle[k].node[n].pos_ref.y = particle[k].center.pos_ref.y + xref * sin(particle[k].angle.z) * cos(particle[k].angle.y) + yref * (sin(particle[k].angle.z) * sin(particle[k].angle.y) * sin(particle[k].angle.x) + cos(particle[k].angle.z) * cos(particle[k].angle.x)) + zref * (sin(particle[k].angle.z) * sin(particle[k].angle.y) * cos(particle[k].angle.x) - cos(particle[k].angle.z) * sin(particle[k].angle.x));
                    particle[k].node[n].pos_ref.z = particle[k].center.pos_ref.z - xref * sin(particle[k].angle.y) + yref * cos(particle[k].angle.y) * sin(particle[k].angle.x) + zref * cos(particle[k].angle.y) * cos(particle[k].angle.x);
                    foreach_dimension()
                    {
                        if (particle[k].node[n].pos_ref.x < -L0 / 2)
                            particle[k].node[n].pos_ref.x += L0;
                        else if (particle[k].node[n].pos_ref.x > L0 / 2)
                            particle[k].node[n].pos_ref.x -= L0;

                        // if (fabs(particle[k].node[n].pos_ref.x - particle[k].node[n].pos.x) > L0 / 2)
                        //     fprintf(stderr,"Wrong:%10d\n",n);
                    }
                }
                #endif // dimension
            }
        }
    }
}
