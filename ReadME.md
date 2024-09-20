# BSE-METISSE

[Public version](https://astronomy.swin.edu.au/~jhurley/bsedload.html) of BSE with a few updates, modified to be compatible with both SSE and METISSE. METISSE is available at https://github.com/TeamMETISSE/METISSE. The documentation for METISSE can be found at this [link](https://metisse.readthedocs.io). 

For more details see: [Agrawal et al 2023](https://arxiv.org/abs/2303.10187) and [Agrawal et al 2020](https://arxiv.org/abs/2005.13177). 

## Installation: 

To clone BSE together with METISSE: 
`git clone --recurse-submodules https://github.com/poojanagrawal/BSE-METISSE.git`


If you have already cloned BSE-METISSE but METISSE directory is empty (you forgot to add `--recurse-submodules`) then do:

`git submodule update --init --recursive`
## Usage:

To use BSE with METISSE, activate `use_SSE = .false.` in the `bse.f` or `popbin.f`.

METISSE specific inputs are provided through `METISSE_input_controls` inlist in the *evolve_metisse.in* file. 

**Important:** When using BSE, the file *evolve_metisse.in* should be located in the same directory as the bse executable. 


