# Installation and usage

**Prequisite:** METISSE requires gcc/6.4.0 or above

## Download METISSE
METISSE is available for download at [here](https://github.com/TeamMETISSE).

*makefile* in folder *make* contains all necessary instructions to compile METISSE in standalone mode. 
Just do `./mk` in the main directory (not in the make directory) to compile the package.

## Using METISSE

Currently, METISSE can be used in the following modes:

1. Standalone mode
2. with BSE (Hurley et al. 2002)
3. with COSMIC (Brievik et al. 2020) (requires additional code from <link> to work)


## Standalone mode

### Supplying input 

METISSE needs to know the values of mass, metallicity and time up to which a star should be evolved. 

Input parameters can be supplied using *evolve_metisse.in* file. It contains two Fortran namelists:

*SSE_input_controls* - For providing parameters such as the mass, metallicity and time up to which a star should evolve.
Only used in standalone mode, ignored otherwise. 

*METISSE_input_controls* - It contains input parameters specific to metisse. 

Both are Fortran namelist, so comments (!) and blank lines can be used freely. Characters are **case-insensitive**. 
Note: Make sure to leave a blank line at the end of the file (after `/` symbol)

Use the format specified in src/defaults/evolve_metisse_defaults.inc for variable names. **Do not modify any file inside the defaults folder**.

### Running METISSE 

To run METISSE in the standalone or main mode, simply do:

`./metisse`

METISSE will evolve a 1 M$_\odot$ star, of input metallicity upto 10 Gyr for you. :blush:

Check *output* directory for output data files.

### Output 
METISSE can produce two types of output files:

**1. files ending with .dat :**

SSE-like output files, controlled by `write_track_to_file` in SSE_input_controls.
These contain the following stellar parameters until max_age. Time and age at hydrogen ZAMS are assumed to be zero.

*time* -- Physical time [Myr]

*age*  -- age of star [Myr]

*mass* -- current mass of the star [M$_\odot$]

*core_mass*        -- mass of dominant core [M$_\odot$]

*He_core*       -- mass of helium core [M$_\odot$]

*CO_core*          -- mass of carbon-oxygen core [M$_\odot$]

*log_L*       -- log of surface luminosity [L$_\odot$]

*log_Teff*     -- log of effective tempertaure [K]

*log_radius*     -- log of radius [R$_\odot$]

*phase*     --   SSE phase number

*e*  -- extra information (not of general relevance- to be removed in future)

**2. files ending with .eep :**

For debugging purposes, METISSE can write mass interpolated files with the same columns as input files (plus phase column) and MIST-style file structure. 
Only contains data from ZAMS to the end of AGB or carbon burning phase.
Can be activated by `write_eep_file` in METISSE_input_controls.

