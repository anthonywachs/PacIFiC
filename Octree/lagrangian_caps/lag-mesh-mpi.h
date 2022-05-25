/**
# MPI synchronization of the Lagrangian mesh

In the current implementation, each processor has its own copy of the Lagrangian
mesh. In order to update the position of all the local copies of the Lagrangian
mesh, the processors need to communicate to each other the contributions of
their flow field to the velocities of each Lagrangian nodes.
*/

void synchronize_mesh(lagMesh* mesh) {
  /** We start by constructing the array of information to send. In this array,
  we store the id of the node, the rank of the processor which manages the
  Eulerian cell the node is located in, the position of the node and its
  velocity */
  int n_localNodes = 0;
  for(int i=0; i<mesh->nlp; i++)
    if (mesh->nodes[i].pid >= 0) n_localNodes++;
  double* send_data = (double*)malloc(n_localNodes*(6*sizeof(double)));
  assert(send_data != NULL);
  int c = 0;
  for(int i=0; i<mesh->nlp; i++) {
    if (mesh->nodes[i].pid >= 0) {
      send_data[6*c] = ((double) i);
      send_data[6*c+1] = ((double) mesh->nodes[i].pid);
      send_data[6*c+2] = mesh->nodes[i].pos.x;
      send_data[6*c+3] = mesh->nodes[i].pos.y;
      send_data[6*c+4] = mesh->nodes[i].lagVel.x;
      send_data[6*c+5] = mesh->nodes[i].lagVel.y;
      c++;
    }
  }

  /** We then broadcast this array to all processors */
  int* nodesPerProc = (int*)malloc(mpi_npe*sizeof(int));
  assert(nodesPerProc != NULL);
  MPI_Allgather(&n_localNodes, 1, MPI_INT, nodesPerProc, 1, MPI_INT,
    MPI_COMM_WORLD);
  int sumNodes = 0;
  int* displacements = (int*)malloc(mpi_npe*sizeof(int));
  for(int i=0; i<mpi_npe; i++) {
    sumNodes += nodesPerProc[i];
    nodesPerProc[i] *= 6;
  }
  assert(sumNodes == mesh->nlp);
  for(int i=0; i<mpi_npe; i++) {
    if (i > 0) displacements[i] = displacements[i-1] + nodesPerProc[i-1];
    else displacements[i] = 0.;
  }
  double* recv_data = (double*)malloc(mesh->nlp*(6*sizeof(double)));
  assert(recv_data != NULL);
  n_localNodes *= 6;
  MPI_Allgatherv(send_data, n_localNodes, MPI_DOUBLE, recv_data, nodesPerProc,
    displacements, MPI_DOUBLE, MPI_COMM_WORLD);
  free(nodesPerProc);
  free(displacements);
  free(send_data);

  /** Finally we read the arrays and update the Lagrangian points */
  for(int i=0; i<mesh->nlp; i++) {
    int k = ((int) recv_data[6*i]);
    mesh->nodes[k].pid = ((int) recv_data[6*i+1]);
    mesh->nodes[k].pos.x = recv_data[6*i+2];
    mesh->nodes[k].pos.y = recv_data[6*i+3];
    mesh->nodes[k].lagVel.x = recv_data[6*i+4];
    mesh->nodes[k].lagVel.y = recv_data[6*i+5];
  }
  free(recv_data);
}

void reduce_lagVel(lagMesh* mesh) {
  int li = dimension; // length of one item
  double* send_data = (double*)malloc(mesh->nlp*li*sizeof(double));
  double* recv_data = (double*)malloc(mpi_npe*mesh->nlp*li*sizeof(double));
  for(int i=0; i<mesh->nlp; i++) {
    send_data[li*i] = mesh->nodes[i].lagVel.x;
    send_data[li*i+1] = mesh->nodes[i].lagVel.y;
    #if dimension > 2
    send_data[li*i+2] = mesh->nodes[i].lagVel.z;
    #endif
  }
  MPI_Allgather(send_data, mesh->nlp*li, MPI_DOUBLE, recv_data,
    mesh->nlp*li, MPI_DOUBLE, MPI_COMM_WORLD);

  for(int i=0; i<mesh->nlp; i++) {
    mesh->nodes[i].lagVel.x = 0.;
    mesh->nodes[i].lagVel.y = 0.;
    #if dimension > 2
    mesh->nodes[i].lagVel.z = 0.;
    #endif
  }

  for(int k=0; k<mpi_npe; k++) {
    for(int i=0; i<mesh->nlp; i++) {
      mesh->nodes[i].lagVel.x += recv_data[k*li*mesh->nlp+li*i];
      mesh->nodes[i].lagVel.y += recv_data[k*li*mesh->nlp+li*i+1];
      #if dimension > 2
      mesh->nodes[i].lagVel.z += recv_data[k*li*mesh->nlp+li*i+2];
      #endif
    }
  }
  free(send_data);
  free(recv_data);
}
