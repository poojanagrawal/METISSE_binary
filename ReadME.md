# BSE-METISSE

[Public version](https://astronomy.swin.edu.au/~jhurley/bsedload.html) of BSE with a few updates, modified to be compatible with both SSE and METISSE. 

METISSE is available at https://github.com/TeamMETISSE/METISSE.
The documentation for METISSE can be found at this link. 

To use BSE with METISSE, activate `use_SSE = .false.` in the `bse.f` or `popbin.f`.

METISSE specific inputs are provided through `METISSE_input_controls` inlist in the *evolve_metisse.in* file. 

**Important:** When using BSE, the file *evolve_metisse.in* should be located in the same directory as the sse/bse/popbin executable. If you're using BSE, the *evolve_metisse.in* file in the METISSE directory will not be read. Instead, you should copy the file to the BSE directory and make the necessary modifications.

For more details see: [Agrawal et al 2023](https://arxiv.org/abs/2303.10187) and [Agrawal et al 2020](https://arxiv.org/abs/2005.13177). 

