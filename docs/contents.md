  
All the programs and subroutines are written in Modern Fortran. METISSE has been developed as an alternative to SSE (Hurley et al. 2000) and so, the subprogram files have similar names and functionality. Module files on the other hand, contain more general data structures and subroutines specific to METISSE that can be accessed by other subprogram units as required. 

The package contains three types of files in the source (*src*) directory: 

### 1. METISSE specific subprogram/subroutines - 

*METISSE_zcnsts.f90*  -   controls metallcity (Z) related part of the package. Reads EEP files and sets Z parameters.

*METISSE_star.f90*     - find relevant tracks from the set and interpolate in mass. 

*METISSE_hrdiag.f90*   - compute stellar parameters by interpolating in age for nuclear burning phases. Also, compute evolutionary phases of the star including the remnant phases and their properties.

*METISSE_deltat.f90*   - determines evolution timestep.

*METISSE_mlwind.f90*   - calculate the mass loss rate (calls mlwind.f for SSE related mass loss). 


### 2. Modules -

*track_support.f90* -- contains general data structures and functions needed throughtout the program.

*interp_support.f90* -- contains routines required for interpolation 

*remnant_support.f90* -- contains functions needed to calculate type and properties of remnant phases

*z_support.f90* -- contains routines to read EEP files, load relevant data and EEPS and other metallicity based parameters

*sse_support.f90* - contains routines to assign stellar phases to the mass interpolated track, to evolve He stars, calculate SSE related quantities.

### 3. Other files - 

A combination of these files are used depending on how METISSE is being used.

**In the standalone mode**

*main_metisse.f90*      - Main program for running metisse. Can only evolve single stars. Reads the input files, sets up relevant parameters and data structures before evolving stars of given masses. 
*evolv_metisse.f90*    - controls the evolution of each star and write output to files.


*assign_commons_main.f90* - assign input values from SSE input files to metisse variables.

**As part of other codes**

*METISSE_miscellaneous* - This file contains miscellaneous subroutines needed by METISSE to work in standalone or otherwise. Ideally these should be packed in a module but Fortran 77 does not know how to use them. 
    Contains:  
*alloc_track.f90* - allocate the track object. 

*dealloc_track.f90* - deallocate the track object and arrays within.

*initialize_front_end* - Inform METISSE what code is using it.

*assign_commons.f90* - assign input values from SSE input files to metisse variables.

*set_star_type* - set star type to rejuvenated before calling star
