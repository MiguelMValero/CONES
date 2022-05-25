int strideCWIPI = 1;
int it;
int nNotLocatedPoints = 0;

scalar t = mesh.time().value();

//******************************************************************** To fill
// Define mesh
if (t == runTime.deltaTValue()){
    cwipi_define_mesh(cl_coupling_name,
                  nPoints,
                  nCells,
                  coords,
                  connecindex,
                  connec);
}

else {
    cwipi_update_location(cl_coupling_name);
}

//******************************************************************** End To fill

// Send receive
if (rank == 0) printf("        Exchange code 1-2\n");

sendvalues = (double*) malloc(sizeof(double) * nPoints);
receivedvalues = (double*) malloc(sizeof(double) * nPoints);

for (int nbpoint(0); nbpoint < nPoints; ++nbpoint){
  sendvalues[nbpoint] = coords[3 * nbpoint + coord_id];
}
  
sprintf(cl_exchange_name, "exch1");
if (coord_id == 1){
  sprintf(cl_sending_field_name, "cooy");
}
else{
  sprintf(cl_sending_field_name, "coox");
}
sprintf(cl_receiving_field_name, "recv");

// We define the iterations

if (t == runTime.deltaTValue()){
  it = 1;
}
else {
  ++it;
}

std::cout << "cp 1\n";
cwipi_exchange_status_t status = cwipi_exchange(cl_coupling_name,
                                cl_exchange_name,
                                strideCWIPI,
                                it,     // n_step
                                t,   // physical_time
                                cl_sending_field_name,
                                sendvalues,
                                cl_receiving_field_name,
                                receivedvalues,
                                &nNotLocatedPoints);

//******************************************************************** End To fill