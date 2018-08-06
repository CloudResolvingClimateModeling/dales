# dales-lib
This is a dedicated fork of Dales, the Dutch Atmospheric Large-Eddy Simulation model. 
The purpose of dales-lib is to be embedded within a column of an atmospheric model, 
having its own MPI communicator and extensive fortran API that allows it to be called from other (Fortran) codes. The API is to some extent tailored to be called from the OMUSE wrapper of DALES: https://bitbucket.org/omuse/omuse/branch/meteo.

### API functions

The API source code is located in `daleslib.f90`, it contains the public functions

* `initialize(path,mpi_comm,date,time)`: initializes a dales model from a fortran namelist file in `path`. Optionally, one can pass an MPI (sub-)communicator `mpi_comm` for the DALES internal messaging, and integers `date` and `time` to set the stard date and time (used for zenith angle and time axes in netCDF output). 
* `step`: performs a dales time step; it is identical to the original DALES time step, except that it adds user-defined large scale forcings to the total tendencies.
* `run_to(tstop)`: iteratively calls stepping until tyhe model time reached `tstop`.
* `finalize`: clears model, deallocates all of its state arrays.
* `get_model_time`: returns the model time in seconds since start time.
* `allocate_3d(a)`: helper function to allocate the array`a` to the size of the model grid.
* `gatherlayeravg(Al,Ag)`: helper function that gathers the array `Al` on all MPI tasks and writes the slab averages (vertical profile) to the array `Ag` on the rank-0 MPI task.
* `gathervol(Al,Ag)`: helper function that gathers 3D arrays `Al` and glues them to the global 3D array `Ag` on the rank-0 MPI task. 

Furthermore the API contains the global profile arrays `u_tend`, `v_tend`, `thl_tend`, `qt_tend` for large-scale tendencies of resp. x- and y-velocities, potential temperature and total humidity. These arrays can be modified at any moment and their values will be automatically picked up in the next time step.
