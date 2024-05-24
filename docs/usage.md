
## Compiling METISSE

    It can be
It has been integrated with the binary evolution code BSE (Hurley et al. 2002) and can be used to evolve populations of isolated binaries.


To run METISSE, it first needs to be compiled (**requires gcc/6.4.0 or above**).
*makefile* in folder *make* contains all necessary instructions to compile METISSE.
Just do `./mk` in main directory (not in the make directory) to compile the package.

## Supplying input parameters 

**By default, METISSE knows how to read MIST (Choi et al. 2016) like files.**

Input values are supplied using *evolve_metisse.in* file. It is a fortran namelist file, so comments (!) and blank lines can be used freely. Characters are case-insensitive. Use the format specified in defaults/evolve_metisse_defaults.in for variable names. **DO NOT modify any file inside defaults folder**.

Similar to SSE, METISSE needs to know the values of mass, metallicity and max age upto which star should be evolved. 

There are two input lists in evolve_metisse.in - 

### 1. SSE_input_controls
These are read only when METISSE is called directly though its main unit. Otherwise, values provided by the overlying code (e.g., SSE, BSE) are used

### 2. Extra_controls
These are input specific to metisse and are always read. Values from default files are used if not supplied by user.


**In any mode INPUT_FILES_DIR (pointing to location of MIST-like files) should be specified in the extra_controls.**

`INPUT_FILES_DIR` -- Path to the folder from where tracks for interpolation are to be read from.


**METALLICITY**:

`Initial_Z` -- Initial metallicity (Z), required for calculations involving white dwarfs.

**MASS**: 

There are two ways mass values can be supplied to METISSE.

1. `(DEFAULT) read_mass_from_file = .false.`

Use `min_mass` to specify the lower limit, `max_mass` for upper limit and `number_of_tracks` to be evolved uniformly distributed in mass between the two limits.

To evolve just one star, specify `number_of_tracks = 1` and input mass value in `min_mass`. Ignore max_mass.

OR 

2. `read_mass_from_file = .true.`

Provide input masses in a text file (one mass value in a line) and specify the location of that file in `input_mass_file`.
Use `number_of_tracks` to tell how many mass values are to be read from the file.

**AGE**:

Use `max_age` to supply maximum age of the star in Myr.

## Running METISSE 
METISSE can be used as a standalone code or in conjncetion with codes like SSE and BSE for evolving a population of single and binary stars.

### Running METISSE in standalone mode

In main mode(not being called from another code), once the code is compiled and at least `INPUT_FILES_DIR` specified, do

`./metisse`

METISSE will evolve a 1 M$_\odot$ star, of input metallicity upto 10 Gyr for you. :blush:

Check *output* directory for output data files.

### Output (Standalone mode)
METISSE can produce two types of output files. 

1. files ending with .eep :

These are mass interpolated files having same columns as input files (plus phase column) and MIST like file strcuture. Only contains data from ZAMS to the end of AGB or Carbon burning phase.

2. files ending with .dat :

These are SSE like output files, containing following stellar parameters until max_age. Time and age at hydrogen ZAMS is assumed to be zero.

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

These files can be be controlled using *advanced_controls* option in *evolve.in* (see below for details). 
*(in progress)*.

## Reading input from other codes

**By default, METISSE knows how to read MIST (Choi et al. 2016) like files.** This means that the input files should have a similar data structure and eep format as in .eep files of MIST package. 

If the input files are in different format, provide details such as location of header for column names, EEP locations etc. using a copy of *format.dat* file. 
Set `read_eep_files = .false.`
The location of format file should also be specified in *evolve.in* using `format_file` option.

## Specifying input columns for interpolation

To reduce memory usage, you can specify selected columns to be used by METISSE by providing the column names in a text file (one name per line) and path to that file using `key_columns_file` option in evolve.in. 
If this file is not provided, all columns from the input files are used.

## Additional options

**writing in progress**

**An alternative way for providing metallicity details:**

If you are using a grid of models, containing sets of masses at different metallicities, you can choose 
`read_files_from_Z = .true.`

and specify a list for getting folder names, these names are appended to INPUT_FILES_DIR

`Z_folder_list = ''`


**Remnant mass prescriptions**

The type and mass of the NS/BH can  be  calculated  from  one  of  the  following  prescriptions:

(1) Deafult SSE

(2) [Belczynski et al. 2002, ApJ, 572, 407](https://iopscience.iop.org/article/10.1086/340304)

(3) [Eldridge J. J., Tout C. A., 2004, MNRAS, 353, 87](https://ui.adsabs.harvard.edu/abs/2004MNRAS.353...87E/abstract) 

(4) [Belczynski et al. 2008, ApJS, 174, 223](https://iopscience.iop.org/article/10.1086/521026)


The corresponding choices for NS/BH are:

`mxns` is the maximum NS mass.

` BHNS_mass_scheme = 'Belczynski2008'`

Options- "original_SSE", "Belczynski2002", "Belczynski2008","Eldridge_Tout2004"


For White dwarfs, two prescriptions are avilable for calculating luminosity:

(1) p. 85 of [Shapiro S. L., Teukolsky S. A., 1983](https://ui.adsabs.harvard.edu/abs/1983bhwd.book.....S/abstract)

(2) [Hurley J. R., Shara M. M., 2003, ApJ, 589, 179](https://iopscience.iop.org/article/10.1086/374637)


`WD_mass_scheme = 'Modified_mestel'`

Options - "Mestel","Modified_mestel"

Ecflag: 

    Electron_capture = .true.       !allow electron capture supernovae


Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800

    Use_Initial_final_mass_relation = .false.       
    construct_wd_track = .true.


**Time stepping options**

Similar to SSE, METISSE determines timesteps as the fractions of the time spent in a phase

pts1 - 95% of MS, HeMS

pts2 - last 5% of MS, cHeBurn, HeHG, HeGB

pts3 - HG, RGB, EAGB, TPAGB


## format files


For Pols et al. 1998 files

`key_columns_file = 'pols_columns.list'

format_file = 'format_pols.dat`

For Bonn files 

`key_columns_file = 'bist_read_cols.txt'

format_file = 'format_bist.dat'`

