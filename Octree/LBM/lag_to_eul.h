#include "ibm.h"

event compute_particle_forces(i++)
{
    if (i > 0)
    {
        for (int k = 0; k < p_n; k++)
        {
            #if dimension == 2
            const double area = 2. * M_PI * particle[k].radius / particle[k].num_nodes;
            #else // dimension == 3
            const double area = 4. * M_PI * particle[k].radius * particle[k].radius / particle[k].num_nodes;
            #endif
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                #if motion_type == 2
                foreach_dimension()
                    particle[k].center.pos_ref.x = pc[k].x;
                #if dimension == 2
                particle[k].node[n].pos_ref.x = particle[k].center.pos_ref.x + particle[k].radius * cos(2. * M_PI * (double)n / particle[k].num_nodes);
                particle[k].node[n].pos_ref.y = particle[k].center.pos_ref.y + particle[k].radius * sin(2. * M_PI * (double)n / particle[k].num_nodes);
                #elif dimension == 3
                double indice = (double)n + 0.5;
                double phi = acos(1. - 2. * indice / particle[k].num_nodes);
                double theta = M_PI * (1. + sqrt(5.)) * indice;
                particle[k].node[n].pos.x = particle[k].center.pos.x + particle[k].radius * cos(theta) * sin(phi);
                if (particle[k].node[n].pos.x < -L0 / 2)
                    particle[k].node[n].pos.x += L0;
                else if (particle[k].node[n].pos.x > L0 / 2)
                    particle[k].node[n].pos.x -= L0;
                particle[k].node[n].pos.y = particle[k].center.pos.y + particle[k].radius * sin(theta) * sin(phi);
                if (particle[k].node[n].pos.y < -L0 / 2)
                    particle[k].node[n].pos.y += L0;
                else if (particle[k].node[n].pos.y > L0 / 2)
                    particle[k].node[n].pos.y -= L0;
                particle[k].node[n].pos.z = particle[k].center.pos.z + particle[k].radius * cos(phi);
                if (particle[k].node[n].pos.z < -L0 / 2)
                    particle[k].node[n].pos.z += L0;
                else if (particle[k].node[n].pos.z > L0 / 2)
                    particle[k].node[n].pos.z -= L0;
                #endif // dimension
                #endif // motion_type
            
                coord add_period = {0,0,0};
                coord mov_period = {0,0,0};
                foreach_dimension()
                {
                    if ((particle[k].center.pos_ref.x + L0/2) < particle[k].radius && particle[k].node[n].pos.x > 0.)
                        add_period.x = -L0;
                    else if ((L0/2 - particle[k].center.pos_ref.x) < particle[k].radius && particle[k].node[n].pos.x < 0.)
                        add_period.x = L0;
                    if (particle[k].node[n].pos_ref.x - particle[k].node[n].pos.x < -L0/2)
                    {
                        // fprintf(stderr,"%12.8f\t%12.8f\n",particle[k].node[n].pos_ref.x,particle[k].node[n].pos.x);
                        mov_period.x = -L0;
                    }
                    else if (particle[k].node[n].pos_ref.x - particle[k].node[n].pos.x > L0/2)
                        mov_period.x = L0;
                }

                foreach_dimension()
                    particle[k].node[n].lag_F.x = -particle[k].stiffness * (particle[k].node[n].pos.x - particle[k].node[n].pos_ref.x + mov_period.x) * area;

                #if dimension == 2
                particle[k].node[n].lag_T.z = (particle[k].node[n].pos.x - particle[k].center.pos_ref.x)*particle[k].node[n].lag_F.y
                                              - (particle[k].node[n].pos.y - particle[k].center.pos_ref.y)*particle[k].node[n].lag_F.x;
                #elif dimension == 3
                foreach_dimension()
                    particle[k].node[n].lag_T.z = (particle[k].node[n].pos.x - particle[k].center.pos_ref.x + add_period.x)*particle[k].node[n].lag_F.y 
                                              - (particle[k].node[n].pos.y - particle[k].center.pos_ref.y + add_period.y)*particle[k].node[n].lag_F.x;
                #endif // dimension
            }
        }
    }
}

event spread(i++)
{
    if (i > 0)
    {
        Rearrange_indices(particle);
        // reset forces
        reset({force}, 0.);

        for (int k = 0; k < p_n; k++)
        {

            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                foreach_cache(particle[k].node[n].my_stencil)
                {
                    if (level > -1)
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
                    // fprintf(stderr,"HHHHHHHHHHHHHHHH\n");
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

                        // compute lattice force.
                        foreach_dimension()
                            #if dimension == 2
                            force.x[] += (particle[k].node[n].lag_F.x * weight_x * weight_y);
                            #else // dimension == 3
                            force.x[] += (particle[k].node[n].lag_F.x * weight_x * weight_y * weight_z);
                            #endif
                    }
                    }
                }
            }
        }
    }
}
