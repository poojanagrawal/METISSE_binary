
Package contents
======================


All the programs and subroutines are written in Modern Fortran. METISSE has been developed as an alternative to SSE (Hurley et al. 2000) and so, the subprogram files have similar names and functionality. Module files on the other hand, contain more general data structures and subroutines specific to METISSE that can be accessed by other subprogram units as required. 

The package contains three types of files in the source (*src*) directory: 

METISSE Specific Subprograms/Subroutines
-----------------------------------------

- *METISSE_zcnsts.f90*: Controls metallcity (Z) related part of the package. Reads EEP files and sets Z parameters.
- *METISSE_star.f90*: Finds relevant tracks from the set and interpolates in mass. 
- *METISSE_hrdiag.f90*: Computes stellar parameters by interpolating in age for nuclear burning phases. Also deteremines evolutionary phases of the star including the remnant phases and their properties. 
- *METISSE_deltat.f90*: Determines evolution timestep.
- *METISSE_mlwind.f90*: Calculates the mass loss rate (calls mlwind.f for SSE related mass loss). 

Modules
-------

- *track_support.f90*: Contains general data structures and functions needed throughout the program.
- *interp_support.f90*: Contains routines required for interpolation.
- *remnant_support.f90*: Contains functions needed to calculate type and properties of remnant phases.
- *z_support.f90*: Contains routines to read EEP files, load relevant data and EEPS and other metallicity based parameters.
- *sse_support.f90*: Contains routines to calculate SSE related quantities.

Other Files
-----------

A combination of these files are used depending on how METISSE is being used.

In the standalone mode:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- *main_metisse.f90*: Main program for running METISSE. Can only evolve single stars. Reads the input files, sets up relevant parameters and data structures before evolving stars of given masses. 
- *evolv_metisse.f90*: Controls the evolution of each star and writes output to files.
- *assign_commons_main.f90*: Assigns input values from SSE input files to METISSE variables.

As part of other codes:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- *METISSE_miscellaneous*: This file contains miscellaneous subroutines needed by METISSE to work in standalone or otherwise. Ideally, these should be packed in a module but Fortran 77 does not know how to use them. 

Contains:  

  - *alloc_track.f90*: Allocates the track object.
  - *dealloc_track.f90*: Deallocates the track object and arrays within.
  - *initialize_front_end*: Informs METISSE what code is using it.
  - *assign_commons.f90*: Assigns input values from SSE input files to METISSE variables.
  - *set_star_type*: Sets star type to rejuvenated before calling star.


   
