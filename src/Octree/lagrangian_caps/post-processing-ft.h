/**
# Post-processing functions for triangulated surfaces
*/

/** File names definition and global variables */
# ifndef fluid_dump_filename
#   define fluid_dump_filename "Savings/dump"
# endif
# ifndef dump_dir
#   define dump_dir "Savings"
# endif
# ifndef result_dir
#   define result_dir "Res"
# endif
# ifndef result_fluid_rootfilename
#   define result_fluid_rootfilename "fluid_basilisk"
# endif
# ifndef figs_dir
#   define figs_dir "Figs"
# endif

/** ## Output membrane position in plain text for external post-processing */
void dump_plain_nodes_pos(lagMesh* mesh, char* filename) {
  if (pid() == 0) {
    FILE* file = fopen(filename, "a");
    assert(file);
    fprintf(file, "%g", t);
    for(int i=0; i<mesh->nln; i++) {
      fprintf(file, ",%g %g %g", mesh->nodes[i].pos.x, mesh->nodes[i].pos.y,
      mesh->nodes[i].pos.z);
    }
    fprintf(file, "\n");
    fclose(file);
    }
}

void dump_plain_triangles(lagMesh* mesh, char* filename) {
  if (pid() == 0) {
    FILE* file = fopen(filename, "a");
    assert(file);
    fprintf(file, "%g", t);
    for(int i=0; i<mesh->nlt; i++) {
      fprintf(file, ",%d %d %d %g", mesh->triangles[i].node_ids[0],
      mesh->triangles[i].node_ids[1], mesh->triangles[i].node_ids[2],
      mesh->triangles[i].area);
    }
    fprintf(file, "\n");
    fclose(file);
  }
}

/** ## Visualization in paraview */
struct _pv_output_ascii {
  char* name;
  FILE* fp;
};

// int pv_timestep = 0;
// void pv_output_ascii(struct _pv_output_ascii p) {
void pv_output_ascii(int pv_timestep) {
  if (pid() == 0) {
    char name[128];
    char default_name[15];
    sprintf(default_name, "%s", result_dir );
    strcat(default_name, "/" );
    strcat(default_name, "caps\0" );
    // char* prefix = p.name ? p.name : default_name;
    char* prefix = default_name;
    char suffix[64];
    sprintf(suffix, "_T%d.vtk", pv_timestep);
    sprintf(name, "%s%s", prefix, suffix);
    // FILE* file = p.fp ? p.fp : fopen(name, "w");
    FILE* file = fopen(name, "w");
    assert(file);

    /* Populate the header and other non-data fields */
    fprintf(file, "# vtk DataFile Version 4.2\n");
    fprintf(file, "Capsules at time %g\n", t);
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    
    /**
    First, we find the number of nodes and triangles to display
    */
    int nbpts_tot = 0;
    int nbtri_tot = 0;
    for(int j=0; j<NCAPS; j++) {
      if (CAPS(j).isactive) {
        nbpts_tot += CAPS(j).nln;
        nbtri_tot += CAPS(j).nlt;
      }
    }

    /* Populate the coordinates of all the Lagrangian nodes */
    fprintf(file, "POINTS %d double\n", nbpts_tot);
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nln; k++) {
        coord node_pos = correct_periodic_node_pos(
          CAPS(j).nodes[k].pos, CAPS(j).centroid);
        fprintf(file, "%g %g %g\n", node_pos.x, node_pos.y, node_pos.z);
      }
    }

    /* Populate the connectivity of the triangles */
    fprintf(file, "TRIANGLE_STRIPS %d %d\n", nbtri_tot, 
    4*nbtri_tot);
    int node_offset = 0;
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++) {
      fprintf(file, "%d %d %d %d\n", 3,
        CAPS(j).triangles[k].node_ids[0] + node_offset,
        CAPS(j).triangles[k].node_ids[1] + node_offset,
        CAPS(j).triangles[k].node_ids[2] + node_offset);
      }
      node_offset += CAPS(j).nln;
    }
    fprintf(file, "CELL_DATA %d\n", nbtri_tot);
    fprintf(file, "SCALARS T1 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {  
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].tension[0]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS T2 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].tension[1]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Tmax double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", max(CAPS(j).triangles[k].tension[0], 
          CAPS(j).triangles[k].tension[1]));
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Tavg double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", .5*(CAPS(j).triangles[k].tension[0] + 
          CAPS(j).triangles[k].tension[1]));
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Lambda_1 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].stretch[0]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Lambda_2 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].stretch[1]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Radius double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).cap_radius);
    }
    fprintf(file, "\n");
    fclose(file);
  }
  pv_timestep++;
}

// return true if the file specified by the filename exists
bool file_exists(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    bool is_exist = false;
    if (fp != NULL)
    {
        is_exist = true;
        fclose(fp); // close the file
    }
    return is_exist;
}

//----------------------------------------------------------------------------
int reinitialize_restart( void )
//----------------------------------------------------------------------------
{
  int init_cycle_number = -1;

  // Get the last cycle cumber from previous simulation
  char filename_lcn[80] = "";
  sprintf( filename_lcn, "%s", result_dir );
  strcat( filename_lcn, "/" );
  strcat( filename_lcn, result_fluid_rootfilename );
  strcat( filename_lcn, "_lcn_vtk.txt" );

  if(file_exists(filename_lcn))
  {
  FILE * fpvtk = fopen( filename_lcn, "r" );
  fscanf ( fpvtk, "%d", &init_cycle_number );
  fclose( fpvtk );

  char buffer_name[64];
  char dump_name[64];
  sprintf(dump_name, "%s", dump_dir );
  strcat(dump_name, "/" );
  sprintf(buffer_name, "flow_%d.dump", init_cycle_number);
  strcat(dump_name, buffer_name);
  restore(dump_name);

  sprintf(dump_name, "%s", dump_dir );
  strcat(dump_name, "/" );
  sprintf(buffer_name, "caps_%d.dump", init_cycle_number);
  strcat(dump_name, buffer_name);
  restore_capsules(dump_name);
  // return init_cycle_number + 1;
  }
  // else
  // {
    return init_cycle_number;
  // }
}


//----------------------------------------------------------------------------
void dump_cycle_number( int nb_dump )
//----------------------------------------------------------------------------
{
  FILE * fpvtk;
      // Write the last cycle number in a file for restart
  if ( pid() == 0 )
  {
    char filename_lcn[256] = "";
    sprintf( filename_lcn, "%s", result_dir );
    strcat( filename_lcn, "/" );
    strcat( filename_lcn, "fluid_basilisk");
    strcat( filename_lcn, "_lcn_vtk.txt" );

    fpvtk = fopen( filename_lcn, "w" );

    fprintf( fpvtk, "%d\n", nb_dump );

    fclose( fpvtk );
  }

}

//----------------------------------------------------------------------------
void output_caps_node_tri()
//----------------------------------------------------------------------------
{
    for(int k=0; k<NCAPS; k++) {
      char fposname[64];
      char ftriname[64];
      sprintf(fposname, "%s/mb_%d_pos.csv", result_dir, k);
      sprintf(ftriname, "%s/mb_%d_tri.csv", result_dir, k);
      dump_plain_nodes_pos(&CAPS(k), fposname);
      dump_plain_triangles(&CAPS(k), ftriname);
    }
}




void output_physics(int cap_number, int iter) {

  double top_visc_stress = 0;
  int top_nb_cells = 0;
  foreach_boundary(top, reduction(+:top_visc_stress) reduction(+:top_nb_cells)) {
    top_nb_cells++;
    top_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  top_visc_stress *= MU/(sq(L0));

  double bottom_visc_stress = 0;
  int bottom_nb_cells = 0;
  foreach_boundary(bottom, reduction(+:bottom_visc_stress) reduction(+:bottom_nb_cells)) {
    bottom_nb_cells++;
    bottom_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  bottom_visc_stress *= MU/(sq(L0));

  double fluid_visc_stress = (top_visc_stress + bottom_visc_stress) / 2.;

 

////////////////////////////////////// Particle Stresslet

  double* send_stress_pack = (double*)calloc(1*4, sizeof(double));
  double* recv_stress_pack = (double*)calloc(1*4, sizeof(double));

  double pN1 = 0.;
  double pN2 = 0.;
  double pmu = 0.;
  double ppres = 0.;

  int k = cap_number;

    double sigmaxy = 0.;
    double sigmaxx = 0.;
    double sigmayy = 0.;
    double sigmazz = 0.;

    Point point = locate(CAPS(k).centroid.x, CAPS(k).centroid.y, CAPS(k).centroid.z);

    if(point.level > -1)
    {

	for(int i=0; i<CAPS(k).nln; i++) 
	{
	double rx, ry, rz;
        rx = CAPS(k).centroid.x + GENERAL_1DIST(CAPS(k).nodes[i].pos.x, CAPS(k).centroid.x, L0*L0_ratio.x);
        ry = CAPS(k).centroid.y + GENERAL_1DIST(CAPS(k).nodes[i].pos.y, CAPS(k).centroid.y, L0*L0_ratio.y);
        rz = CAPS(k).centroid.z + GENERAL_1DIST(CAPS(k).nodes[i].pos.z, CAPS(k).centroid.z, L0*L0_ratio.z);

	#ifndef CAPS_VISCOSITY
	sigmaxx += - CAPS(k).nodes[i].lagForce.x * rx;
	sigmayy += - CAPS(k).nodes[i].lagForce.y * ry;
	sigmazz += - CAPS(k).nodes[i].lagForce.z * rz;
	sigmaxy += - (CAPS(k).nodes[i].lagForce.x * ry + CAPS(k).nodes[i].lagForce.y * rx) / 2.;
	#else 
	
        /** We now have to compute the area associated with each node */
        double visc_ratio = 1./MUP * MUC;
        double nodal_area = compute_node_area(&(CAPS(k)), i);
        
	sigmaxx += - CAPS(k).nodes[i].lagForce.x * rx + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.x*CAPS(k).nodes[i].normal.x)*nodal_area;
	sigmayy += - CAPS(k).nodes[i].lagForce.y * ry + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.y*CAPS(k).nodes[i].normal.y)*nodal_area;
	sigmazz += - CAPS(k).nodes[i].lagForce.z * rz + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.z*CAPS(k).nodes[i].normal.z)*nodal_area;
	sigmaxy += - (CAPS(k).nodes[i].lagForce.x * ry + CAPS(k).nodes[i].lagForce.y * rx) / 2. 
	           + MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.x*CAPS(k).nodes[i].normal.y + CAPS(k).nodes[i].lagVel.y*CAPS(k).nodes[i].normal.x)*nodal_area;
	#endif	

	}

    pN1 = (sigmaxx - sigmayy);
    pN2 = (sigmayy - sigmazz);
    pmu = sigmaxy;
    ppres = -(sigmaxx + sigmayy + sigmazz)/3.;

    send_stress_pack[0*4] = pN1;
    send_stress_pack[0*4 + 1] = pN2;
    send_stress_pack[0*4 + 2] = pmu;
    send_stress_pack[0*4 + 3] = ppres;
    }

  pN1 = 0.;
  pN2 = 0.;
  pmu = 0.;
  ppres = 0.;
 

MPI_Reduce(send_stress_pack, recv_stress_pack, 4*1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

if (pid() == 0) 
{
    pN1 = 0.;
    pN2 = 0.;
    pmu = 0.;
    ppres = 0.;

  //  for(int k = 0; k < NCAPS; k++)
  //  {
	pN1 += recv_stress_pack[0*4];
	pN2 += recv_stress_pack[0*4 + 1];
	pmu += recv_stress_pack[0*4 + 2];
	ppres += recv_stress_pack[0*4 + 3];
  //  }

}
free(send_stress_pack);
free(recv_stress_pack);

//////////////////////////////////////////////////////////////////////////////////////


  if (pid() == 0) 
  {
    double cap_area = 0;
    double cap_volume = 0;
    cap_volume += CAPS(cap_number).volume/CAPS(cap_number).initial_volume;
    for(int i=0; i<CAPS(cap_number).nlt; i++) cap_area += CAPS(cap_number).triangles[i].area/(4*pi*sq(CAPS(cap_number).cap_radius)); 

  /*Compute average Taylor deformation and angular velocity*/
    double TDmaxmin = 0;
    double TDang = 0;
    double taylor_deform = 0;
    double inclin_angle = 0;
    double ang_vel = 0;
    coord rs = {0., 0., 0.};
    coord centers = {0., 0., 0.};
    compute_taylor_factor(&CAPS(cap_number), &taylor_deform, &inclin_angle, &rs, &TDmaxmin, &TDang);
    ang_vel = CAPS(cap_number).ang_vel.z;
    foreach_dimension() 
    {
      rs.x = rs.x/CAPS(cap_number).cap_radius;
      centers.x = CAPS(cap_number).centroid.x;
    } 

  
    char name[128];
    char default_name[30];
    sprintf(default_name, "%s", result_dir );
    strcat(default_name, "/" );
    strcat(default_name, "Individual_outputs\0" );
    // char* prefix = p.name ? p.name : default_name;
    char* prefix = default_name;
    char suffix[64];
    sprintf(suffix, "_cap%d.txt", cap_number);
    sprintf(name, "%s%s", prefix, suffix);
    FILE* foutput_physics = fopen(name, "a+");
    assert(foutput_physics);

    fprintf(foutput_physics, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", iter, t, 
      fluid_visc_stress, ppres, pmu, pN1, pN2, taylor_deform, inclin_angle, TDmaxmin, TDang,
      rs.x, rs.y, rs.z, ang_vel, cap_area, cap_volume, centers.x, centers.y, centers.z);
    fflush(foutput_physics);

    fclose(foutput_physics);
  }
   
}

void output_bidispse_physics(int N_pops, int iter) {

  double top_visc_stress = 0;
  int top_nb_cells = 0;
  foreach_boundary(top, reduction(+:top_visc_stress) reduction(+:top_nb_cells)) {
    top_nb_cells++;
    top_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  top_visc_stress *= MU/(sq(L0));

  double bottom_visc_stress = 0;
  int bottom_nb_cells = 0;
  foreach_boundary(bottom, reduction(+:bottom_visc_stress) reduction(+:bottom_nb_cells)) {
    bottom_nb_cells++;
    bottom_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  bottom_visc_stress *= MU/(sq(L0));

  double fluid_visc_stress = (top_visc_stress + bottom_visc_stress) / 2.;


////////////////////////////////////// Particle Stresslet

  double* send_stress_pack = (double*)calloc(NCAPS*4, sizeof(double));
  double* recv_stress_pack = (double*)calloc(NCAPS*4, sizeof(double));

  double pN1 = 0.;
  double pN2 = 0.;
  double pmu = 0.;
  double ppres = 0.;

  for(int k = 0; k < NCAPS; k++)
  {

    double sigmaxy = 0.;
    double sigmaxx = 0.;
    double sigmayy = 0.;
    double sigmazz = 0.;

    Point point = locate(CAPS(k).centroid.x, CAPS(k).centroid.y, CAPS(k).centroid.z);

    if(point.level > -1)
    {

	for(int i=0; i<CAPS(k).nln; i++) 
	{
	/** The post-processing is only carried out if we are in the shear plane */ 
	double rx, ry, rz;

        rx = CAPS(k).centroid.x + GENERAL_1DIST(CAPS(k).nodes[i].pos.x, CAPS(k).centroid.x, L0*L0_ratio.x);
        ry = CAPS(k).centroid.y + GENERAL_1DIST(CAPS(k).nodes[i].pos.y, CAPS(k).centroid.y, L0*L0_ratio.y);
        rz = CAPS(k).centroid.z + GENERAL_1DIST(CAPS(k).nodes[i].pos.z, CAPS(k).centroid.z, L0*L0_ratio.z);

	#ifndef CAPS_VISCOSITY
	sigmaxx += - CAPS(k).nodes[i].lagForce.x * rx;
	sigmayy += - CAPS(k).nodes[i].lagForce.y * ry;
	sigmazz += - CAPS(k).nodes[i].lagForce.z * rz;
	sigmaxy += - (CAPS(k).nodes[i].lagForce.x * ry + CAPS(k).nodes[i].lagForce.y * rx) / 2.;
	#else 
	double visc_ratio = 1./MUP * MUC;
          /** We now have to compute the area associated with each node */
        double nodal_area = compute_node_area(&(CAPS(k)), i);
        
	sigmaxx += - CAPS(k).nodes[i].lagForce.x * rx + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.x*CAPS(k).nodes[i].normal.x)*nodal_area;
	sigmayy += - CAPS(k).nodes[i].lagForce.y * ry + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.y*CAPS(k).nodes[i].normal.y)*nodal_area;
	sigmazz += - CAPS(k).nodes[i].lagForce.z * rz + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.z*CAPS(k).nodes[i].normal.z)*nodal_area;
	sigmaxy += - (CAPS(k).nodes[i].lagForce.x * ry + CAPS(k).nodes[i].lagForce.y * rx) / 2. 
	           + MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.x*CAPS(k).nodes[i].normal.y + CAPS(k).nodes[i].lagVel.y*CAPS(k).nodes[i].normal.x)*nodal_area;
	#endif	

	}

    pN1 = (sigmaxx - sigmayy);
    pN2 = (sigmayy - sigmazz);
    pmu = sigmaxy;
    ppres = -(sigmaxx + sigmayy + sigmazz)/3.;

    send_stress_pack[k*4] = pN1;
    send_stress_pack[k*4 + 1] = pN2;
    send_stress_pack[k*4 + 2] = pmu;
    send_stress_pack[k*4 + 3] = ppres;
    }

  pN1 = 0.;
  pN2 = 0.;
  pmu = 0.;
  ppres = 0.;
 }

 MPI_Reduce(send_stress_pack, recv_stress_pack, 4*NCAPS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


int* pop_count = (int*)calloc(N_pops, sizeof(int));
int tot_count = 0;
double* pN1_p = (double*)calloc(N_pops, sizeof(double));
double* pN2_p = (double*)calloc(N_pops, sizeof(double));
double* pmu_p = (double*)calloc(N_pops, sizeof(double));
double* ppres_p = (double*)calloc(N_pops, sizeof(double));

if (pid() == 0) 
{
tot_count = 0;
pN1 = 0.;
pN2 = 0.;
pmu = 0.;
ppres = 0.;

   for(int k = 0; k < NCAPS; k++)
   {
    for (int j = 0; j < N_pops; j++)
    { 
      if(CAPS(k).cap_type == j)
      {
        pN1_p[j] += recv_stress_pack[k*4];
        pN2_p[j] += recv_stress_pack[k*4 + 1];
        pmu_p[j] += recv_stress_pack[k*4 + 2];
        ppres_p[j] += recv_stress_pack[k*4 + 3];
        pop_count[j] ++;
        tot_count ++;
      }
    }
   }

   assert(tot_count == NCAPS && "Error pop count is different from NCAPS!\n");
 
}
free(send_stress_pack);
free(recv_stress_pack);

// We compute the averaged quantities by population
//////////////////////////////////////////////////////////////////////////////////////
  for(int pop_type = 0; pop_type < N_pops; pop_type++)
  {
    if (pid() == 0) 
    {
      double avg_ncaps_area = 0;
      double avg_ncaps_volume = 0;
      for(int k = 0; k < NCAPS; k++) {
        if(CAPS(k).cap_type == pop_type)
        {
          avg_ncaps_volume += CAPS(k).volume/CAPS(k).initial_volume;
          for(int i=0; i<CAPS(k).nlt; i++) avg_ncaps_area += CAPS(k).triangles[i].area/(4*pi*sq(CAPS(k).cap_radius));
        }
      }
      
      avg_ncaps_area /= pop_count[pop_type];
      avg_ncaps_volume /= pop_count[pop_type];

    /*Compute average Taylor deformation and angular velocity*/
      double avg_TDmaxmin = 0;
      double avg_TDang = 0;
      double TDmaxmin = 0;
      double TDang = 0;
      double avg_taylor_deform = 0;
      double avg_inclin_angle = 0;
      double taylor_deform = 0;
      double inclin_angle = 0;
      double avg_ang_vel = 0;
      coord avg_rs = {0., 0., 0.};
      coord rs = {0., 0., 0.};
      for(int k=0; k<NCAPS; k++) {
        if(CAPS(k).cap_type == pop_type)
        {
          compute_taylor_factor(&CAPS(k), &taylor_deform, &inclin_angle, &rs, &TDmaxmin, &TDang);
          foreach_dimension() avg_rs.x += rs.x/CAPS(k).cap_radius; 
          avg_TDmaxmin += TDmaxmin;
          avg_TDang += TDang;
          avg_taylor_deform += taylor_deform;
          avg_inclin_angle += inclin_angle;
          avg_ang_vel += CAPS(k).ang_vel.z;
        }
      }
      
      foreach_dimension() avg_rs.x /= pop_count[pop_type]; 
      avg_TDmaxmin /= pop_count[pop_type];
      avg_TDang /= pop_count[pop_type];
      avg_inclin_angle /= pop_count[pop_type];
      avg_taylor_deform /= pop_count[pop_type];
      avg_ang_vel /= pop_count[pop_type];

      char name[128];
      char default_name[30];
      sprintf(default_name, "%s", result_dir );
      strcat(default_name, "/" );
      strcat(default_name, "Individual_outputs\0" );
      // char* prefix = p.name ? p.name : default_name;
      char* prefix = default_name;
      char suffix[64];
      sprintf(suffix, "_cap%d.txt", pop_type);
      sprintf(name, "%s%s", prefix, suffix);
      FILE* foutput_physics = fopen(name, "a+");
      assert(foutput_physics);

      fprintf(foutput_physics, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", iter, t, 
        ppres_p[pop_type], pmu_p[pop_type], pN1_p[pop_type], pN2_p[pop_type], 
        avg_taylor_deform, avg_inclin_angle, avg_TDmaxmin, avg_TDang,
        avg_rs.x, avg_rs.y, avg_rs.z, avg_ang_vel, avg_ncaps_area, avg_ncaps_volume);
      fflush(foutput_physics);

      fclose(foutput_physics);
    }
  }


// We compute the ensemble average quantities
//////////////////////////////////////////////////////////////////////////////////////
  if (pid() == 0) 
    {
      double avg_ncaps_area = 0;
      double avg_ncaps_volume = 0;
      for(int k = 0; k < NCAPS; k++) {
          avg_ncaps_volume += CAPS(k).volume/CAPS(k).initial_volume;
          for(int i=0; i<CAPS(k).nlt; i++) avg_ncaps_area += CAPS(k).triangles[i].area/(4*pi*sq(CAPS(k).cap_radius));
      }
      
      avg_ncaps_area /= NCAPS;
      avg_ncaps_volume /= NCAPS;

    /*Compute average Taylor deformation and angular velocity*/
      double avg_TDmaxmin = 0;
      double avg_TDang = 0;
      double TDmaxmin = 0;
      double TDang = 0;
      double avg_taylor_deform = 0;
      double avg_inclin_angle = 0;
      double taylor_deform = 0;
      double inclin_angle = 0;
      double avg_ang_vel = 0;
      coord avg_rs = {0., 0., 0.};
      coord rs = {0., 0., 0.};
      for(int k=0; k<NCAPS; k++) {
          compute_taylor_factor(&CAPS(k), &taylor_deform, &inclin_angle, &rs, &TDmaxmin, &TDang);
          foreach_dimension() avg_rs.x += rs.x/CAPS(k).cap_radius; 
          avg_TDmaxmin += TDmaxmin;
          avg_TDang += TDang;
          avg_taylor_deform += taylor_deform;
          avg_inclin_angle += inclin_angle;
          avg_ang_vel += CAPS(k).ang_vel.z;
      }
      
      foreach_dimension() avg_rs.x /= NCAPS; 
      avg_TDmaxmin /= NCAPS;
      avg_TDang /= NCAPS;
      avg_inclin_angle /= NCAPS;
      avg_taylor_deform /= NCAPS;
      avg_ang_vel /= NCAPS;

      char name[128];
      // char default_name[30];
      // sprintf(default_name, "%s", result_dir );
      // strcat(default_name, "/" );
      // strcat(default_name, "Individual_outputs\0" );
      // char* prefix = p.name ? p.name : default_name;
      // char* prefix = default_name;
      // char suffix[64];
      // sprintf(suffix, "_cap%d.txt", pop_type);
      sprintf(name, "output.txt");
      FILE* fouttot_physics = fopen(name, "a+");
      assert(fouttot_physics);

      fprintf(fouttot_physics, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", iter, t, 
        fluid_visc_stress, ppres_p[0]+ppres_p[1], pmu_p[0]+pmu_p[1], pN1_p[0]+pN1_p[1], pN2_p[0]+pN2_p[1], 
        avg_taylor_deform, avg_inclin_angle, avg_TDmaxmin, avg_TDang,
        avg_rs.x, avg_rs.y, avg_rs.z, avg_ang_vel, avg_ncaps_area, avg_ncaps_volume);
      fflush(fouttot_physics);

      fclose(fouttot_physics);
    }
 


free(pop_count);
free(pN1_p);
free(pN2_p);
free(pmu_p);
free(ppres_p);

}
