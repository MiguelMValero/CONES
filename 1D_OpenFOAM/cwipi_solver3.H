if (rank == 0) printf("        Delete coupling\n");

//******************************************************************** To fill
// Delete coupling
cwipi_delete_coupling(cl_coupling_name);

//******************************************************************** End To fill

// Freeing memory
free(coords);
free(connecindex);
free(connec);
free(sendvalues);
free(receivedvalues);

//******************************************************************** To fill
// Finalize
cwipi_finalize();
MPI_Finalize();

//******************************************************************** End To fill