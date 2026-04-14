/**
# MPI synchronization of the Lagrangian mesh

In the current implementation, each processor has its own copy of the Lagrangian
mesh. In order to update the position of all the local copies of the Lagrangian
mesh, the processors need to communicate to each other the contributions of
their flow field to the velocities of each Lagrangian nodes.
*/

bool on_face2(double p, int n, double l0) {
  if ((fabs(p/(l0/n)) - ((int)fabs(p/(l0/n)))) < 1.e-10) return true;
  else return false;
}

double correct_point_pos(double ptx, double a, double L) {
   double origin = a + L/2; // L is the length of the domain in the current dimension
    if (on_face2(ptx, N, L))
      ptx += 1.e-10;
    //FIXME: the nodes should not be sent to the other side of the domain if the boundary is not periodic...
    if (ptx > origin + L/2)
      ptx -= L;
    else if (ptx < origin - L/2)
      ptx += L;
  return ptx;
}

/*Two-sides periodicity*/
#define POS_PBC_X(X) ((u.x.boundary[left] != periodic_bc) ? (X) : correct_point_pos(X, X0, L0))
#define POS_PBC_Y(Y) ((u.x.boundary[top] != periodic_bc) ? (Y) : correct_point_pos(Y, Y0, L0*Dimensions.y/Dimensions.x))
#define POS_PBC_Z(Z) ((u.x.boundary[front] != periodic_bc) ? (Z) : correct_point_pos(Z, Z0, L0*Dimensions.z/Dimensions.x))

/*Onside periodicity, might be useful for visualization*/
#define vPOS_PBC_X(X) ((u.x.boundary[left] != periodic_bc) ? (X) : (((X - (X0 +\
  L0/2)) > L0/2.) ? (X) - L0 : (X)))
#define vPOS_PBC_Y(Y) ((u.x.boundary[top] != periodic_bc) ? (Y) : (((Y - (Y0 +\
  L0*Dimensions.y/Dimensions.x/2)) > L0*Dimensions.y/Dimensions.x/2.) ? (Y) - L0*Dimensions.y/Dimensions.x : (Y)))
#define vPOS_PBC_Z(Z) ((u.x.boundary[front] != periodic_bc) ? (Z) : (((Z - (Z0 +\
  L0*Dimensions.z/Dimensions.x/2)) > L0*Dimensions.z/Dimensions.x/2.) ? (Z) - L0*Dimensions.z/Dimensions.x : (Z)))


void reduce_lagVel(lagMesh* mesh) {
  int li = dimension; // length of one item
  double* send_data = (double*)malloc(mesh->nln*li*sizeof(double));
  double* recv_data = (double*)malloc(mpi_npe*mesh->nln*li*sizeof(double));
  for(int i=0; i<mesh->nln; i++) {
    send_data[li*i] = mesh->nodes[i].lagVel.x;
    send_data[li*i+1] = mesh->nodes[i].lagVel.y;
    #if dimension > 2
    send_data[li*i+2] = mesh->nodes[i].lagVel.z;
    #endif
  }
  MPI_Allgather(send_data, mesh->nln*li, MPI_DOUBLE, recv_data,
    mesh->nln*li, MPI_DOUBLE, MPI_COMM_WORLD);

  for(int i=0; i<mesh->nln; i++) {
    mesh->nodes[i].lagVel.x = 0.;
    mesh->nodes[i].lagVel.y = 0.;
    #if dimension > 2
    mesh->nodes[i].lagVel.z = 0.;
    #endif
  }

  for(int k=0; k<mpi_npe; k++) {
    for(int i=0; i<mesh->nln; i++) {
      mesh->nodes[i].lagVel.x += recv_data[k*li*mesh->nln+li*i];
      mesh->nodes[i].lagVel.y += recv_data[k*li*mesh->nln+li*i+1];
      #if dimension > 2
      mesh->nodes[i].lagVel.z += recv_data[k*li*mesh->nln+li*i+2];
      #endif
    }
  }
  free(send_data);
  free(recv_data);
}

void reduce_alllagVel() 
{ 
  int total_nln = 0;
  int* pos_in_pack = (int*)malloc(NCAPS*sizeof(int));
  for(int i=0; i<NCAPS; i++) 
  {
      if (CAPS(i).isactive) 
      {
        // eul2lag(&CAPS(i));
        pos_in_pack[i] = total_nln;
        total_nln += CAPS(i).nln;
      }
  }   

  int li = dimension; // length of one item
  double* send_data_pack = (double*)malloc(total_nln*li*sizeof(double));
  double* recv_data_pack = (double*)malloc(total_nln*li*sizeof(double));

  for(int i=0; i<NCAPS; i++) 
    {
      if (CAPS(i).isactive) 
      { 
        for(int j=0; j<CAPS(i).nln; j++) 
        {
            send_data_pack[pos_in_pack[i]*li + li*j] = CAPS(i).nodes[j].lagVel.x;
            send_data_pack[pos_in_pack[i]*li + li*j+1] = CAPS(i).nodes[j].lagVel.y;
            #if dimension > 2
            send_data_pack[pos_in_pack[i]*li + li*j+2] = CAPS(i).nodes[j].lagVel.z;
            #endif
        }
      }
    }   

  MPI_Allreduce(send_data_pack, recv_data_pack, total_nln*li, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(int i=0; i<NCAPS; i++) 
    {
      if (CAPS(i).isactive) 
      { 
        for(int j=0; j<CAPS(i).nln; j++) 
        {
            CAPS(i).nodes[j].lagVel.x = recv_data_pack[pos_in_pack[i]*li + li*j];
            CAPS(i).nodes[j].lagVel.y = recv_data_pack[pos_in_pack[i]*li + li*j+1];
            #if dimension > 2
            CAPS(i).nodes[j].lagVel.z = recv_data_pack[pos_in_pack[i]*li + li*j+2];
            #endif
        }
      }
    }   

  free(pos_in_pack);
  free(send_data_pack);
  free(recv_data_pack);
}


bool is_capsule_in_proc(lagMesh* mesh) {

  Cache c = {0};
  foreach()
  {
    coord checkpt={0};
    checkpt.x = x;
    checkpt.y = y;
    checkpt.z = z;
    double tetative_dist = sqrt(GENERAL_SQNORM(checkpt, mesh->centroid));
    if(tetative_dist <= mesh->circum_radius)
    {
	    cache_append( &c, point, 0 );
    }
  }

  foreach_cache(c)
  {
   if(point.level>-1)
   {
     free(c.p); //free cache
     return true;
   }
  }
  free(c.p); //free cache
  return false;
}


void compute_proc_borders(coord* proc_max, coord* proc_min)
{
  Cache c = {0};
  foreach() cache_append( &c, point, 0 );
  foreach_cache(c)
  {
   coord checkpt={x, y, z}; 
   if(point.level>-1)
   {
    foreach_dimension() if (proc_max->x < checkpt.x) proc_max->x = checkpt.x;
    foreach_dimension() if (proc_min->x > checkpt.x) proc_min->x = checkpt.x;
   }
  }
  free(c.p); //free cache
}

bool is_capsule_in_boundingbox(coord proc_max, coord proc_min, lagMesh* mesh) 
{
  /* Check only if the point is in the AABB (Axed-Aligned-Bounding-Box) */
  coord cap_max = {mesh->centroid.x + mesh->circum_radius, mesh->centroid.y + mesh->circum_radius, mesh->centroid.z + mesh->circum_radius}; 
  coord cap_min = {mesh->centroid.x - mesh->circum_radius, mesh->centroid.y - mesh->circum_radius, mesh->centroid.z - mesh->circum_radius}; 

  /*Check if the capsule box overlating with the current proc*/
  if ( ( cap_max.x >= proc_min.x ) && ( cap_min.x <= proc_max.x ) ) 
    if ( ( cap_max.y >= proc_min.y ) && ( cap_min.y <= proc_max.y ) )
      if ( ( cap_max.z >= proc_min.z ) && ( cap_min.z <= proc_max.z ) )
        return true;      

  /*In case of periodicity, check only if the point is in the AABB (Axed-Aligned-Bounding-Box)*/
  foreach_dimension()
  {
      if(POS_PBC_X(cap_max.x) != cap_max.x) cap_min.x -= L0*L0_ratio.x;
      if(POS_PBC_X(cap_min.x) != cap_min.x) cap_max.x += L0*L0_ratio.x;
  }

  /*Check if the copy capsule box overlating with the current proc*/
  if ( ( cap_max.x >= proc_min.x ) && ( cap_min.x <= proc_max.x ) ) 
    if ( ( cap_max.y >= proc_min.y ) && ( cap_min.y <= proc_max.y ) )
      if ( ( cap_max.z >= proc_min.z ) && ( cap_min.z <= proc_max.z ) )
        return true;    

  /*If not return false*/
  return false;
}


