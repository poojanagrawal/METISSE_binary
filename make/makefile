BSE_DIR=..

SRC_DIR=$(BSE_DIR)/src
SSE_DIR=$(BSE_DIR)/src/SSE
METISSE_DIR=$(BSE_DIR)/src

#com
FC=gfortran
#DEBUG=-g3 -fbacktrace -fcheck=all -fbounds-check
#WARN=-Wall -pedantic -Wno-tabs -Wno-maybe-uninitialized
# -Wconversion -Wextra

#Optimization flag
#OPT=-O3

#NOT USED
#Tells gfortran to compile using OpenMP
#OMP=-fopenmp

WARN=-Wall -Wno-tabs -Wno-maybe-uninitialized -fimplicit-none  -Wsurprising  -Waliasing
DEBUG=-g3 -O2 -fbacktrace -fcheck=all -fbounds-check 


#-O2    -Wall  -Wline-truncation  -Wcharacter-truncation   -Wimplicit-interface  -Wunused-parameter   -fcheck=all      -fbacktrace




FFLAGS1= -ffree-form -fimplicit-none $(WARN) $(DEBUG)
FFLAGS2= -ffixed-form -fimplicit-none $(DEBUG)

POPBIN =popbin

SRCE3 = \
popbin.f comenv.f corerd.f dgcore.f evolv2.f gntage.f \
instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f zfuncs.f \
hrdiag.f star.f zcnsts.f deltat.f mlwind.f \
assign_remnant.f hrdiag_remnant.f

OBJT3= $(SRCE3:.f=.o)

SSE_METISSE=sse_metisse

SRCE4 = \
sse.f evolv1.f kick.f ran3.f mlwind.f mrenv.f zfuncs.f \
hrdiag.f star.f zcnsts.f deltat.f mlwind.f\
assign_remnant.f hrdiag_remnant.f

OBJT4 = $(SRCE4:.f=.o)


BSE_METISSE=bse_metisse

SRCE5 = \
bse.f evolv2.f comenv.f corerd.f dgcore.f gntage.f \
instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f zfuncs.f \
hrdiag.f star.f zcnsts.f deltat.f mlwind.f\
assign_remnant.f hrdiag_remnant.f

OBJT5 = $(SRCE5:.f=.o)


SRCE_SSE = \
SSE_hrdiag.f SSE_star.f SSE_zcnsts.f SSE_deltat.f SSE_mlwind.f SSE_gntage.f
	
OBJT_SSE = $(SRCE_SSE:.f=.o)


SRCE_METISSE = \
track_support.f90 z_support.f90 sse_support.f90 \
remnant_support.f90 interp_support.f90 \
assign_commons_BSE.f90 assign_commons_main.f90 \
alloc_track.f90 dealloc_track.f90 \
comenv_lambda.f90 initialize_front_end.f90 \
METISSE_zcnsts.f90 METISSE_star.f90 METISSE_hrdiag.f90 \
METISSE_deltat.f90 METISSE_mlwind.f90 \
METISSE_gntage.f90
 
OBJT_METISSE = $(SRCE_METISSE:.f90=.o)


track_support.o: $(METISSE_DIR)/track_support.f90
z_support.o: $(METISSE_DIR)/z_support.f90 track_support.o
interp_support.o: $(METISSE_DIR)/interp_support.f90 track_support.o
sse_support.o: $(METISSE_DIR)/sse_support.f90 track_support.o
remnant_support.o: $(METISSE_DIR)/remnant_support.f90 track_support.o sse_support.o

initialize_front_end.o: $(METISSE_DIR)/initialize_front_end.f90 track_support.o
alloc_track.o: $(METISSE_DIR)/alloc_track.f90 track_support.o
dealloc_track.o: $(METISSE_DIR)/dealloc_track.f90 track_support.o
assign_commons_main.o: $(METISSE_DIR)/assign_commons_main.f90 track_support.o remnant_support.o
metisse_mlwind.o: $(METISSE_DIR)/metisse_mlwind.f90 track_support.o
metisse_deltat.o: $(METISSE_DIR)/metisse_deltat.f90 track_support.o
metisse_zcnsts.o: $(METISSE_DIR)/metisse_zcnsts.f90 z_support.o track_support.o
metisse_star.o: $(METISSE_DIR)/metisse_star.f90 track_support.o interp_support.o sse_support.o
metisse_hrdiag.o: $(METISSE_DIR)/metisse_hrdiag.f90 track_support.o interp_support.o sse_support.o remnant_support.o
metisse_gntage.o: $(METISSE_DIR)/metisse_gntage.f90 track_support.o interp_support.o

comenv_lambda.o: $(METISSE_DIR)/comenv_lambda.f90 track_support.o

assign_commons.o: $(METISSE_DIR)/assign_commons_bse.f90 track_support.o remnant_support.o



%.o : $(SRC_DIR)/%.f
	$(FC) $(FFLAGS2) -c $<

%.o : $(SSE_DIR)/%.f
	$(FC) $(FFLAGS2) -c $<
	
%.o : $(METISSE_DIR)/%.f90
	$(FC) $(FFLAGS1) -c $<
	

$(SSE_METISSE): $(OBJT4) $(OBJT_SSE) $(OBJT_METISSE)
	$(FC) $(OMP) -o $(BSE_DIR)/$(SSE_METISSE) $^


$(BSE_METISSE): $(OBJT5) $(OBJT_SSE) $(OBJT_METISSE)
	$(FC) $(OMP) -o $(BSE_DIR)/$(BSE_METISSE) $^


$(POPBIN): $(OBJT3) $(OBJT_SSE) $(OBJT_METISSE)
	$(FC) $(OMP) -o $(BSE_DIR)/$(POPBIN) $^


all:$(BSE_METISSE) $(SSE_METISSE) $(POPBIN)


clean:
	rm -f $(BSE_DIR)/$(SSE_METISSE) $(BSE_DIR)/$(BSE_METISSE) $(BSE_DIR)/$(POPBIN)
	rm -f *.mod *.o

