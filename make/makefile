BSE_DIR=..

SRC_DIR=$(BSE_DIR)/src

#com
FC=gfortran
DEBUG=-g3 -fbacktrace -fcheck=all -fbounds-check
WARN=-Wall -pedantic -Wno-tabs -Wextra
# -Wconversion

#Optimization flag
OPT=-O3

#NOT USED
#Tells gfortran to compile using OpenMP
OMP=-fopenmp


FFLAGS1= -ffree-form -fimplicit-none $(OPT) $(WARN) $(DEBUG)
FFLAGS2= -ffixed-form -fimplicit-none 

SRCE3 = \
popbin.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT3 = $(SRCE3:.f=.o)

#popbin: $(OBJT3) $(LFLAGS)
#	$(CMPLR) $(FFLAGS) $(OBJT3) -o popbin


SSE_METISSE=sse_metisse


track_support.o: $(SRC_DIR)/track_support.f90
z_support.o: $(SRC_DIR)/z_support.f90 track_support.o 
sse_support.o: $(SRC_DIR)/sse_support.f90 track_support.o z_support.o 
remnant_support.o: $(SRC_DIR)/remnant_support.f90 track_support.o sse_support.o
interp_support.o: $(SRC_DIR)/interp_support.f90 track_support.o z_support.o sse_support.o
  
alloc_track.o: $(SRC_DIR)/alloc_track.f90 track_support.o
dealloc_track.o: $(SRC_DIR)/dealloc_track.f90 track_support.o
assign_commons.o: $(SRC_DIR)/assign_commons.f90 track_support.o remnant_support.o
metisse_mlwind.o: $(SRC_DIR)/metisse_mlwind.f90 track_support.o
deltat.o:$(SRC_DIR)/deltat.f90 track_support.o
zcnsts.o: $(SRC_DIR)/zcnsts.f90 z_support.o track_support.o
star.o: $(SRC_DIR)/star.f90 track_support.o interp_support.o sse_support.o 
hrdiag.o: $(SRC_DIR)/hrdiag.f90 track_support.o interp_support.o remnant_support.o sse_support.o


SRCE4 = \
sse.f evolv1.f kick.f ran3.f mlwind.f mrenv.f zfuncs.f

OBJT4 = $(SRCE4:.f=.o)


BSE_METISSE=bse_metisse

all:$(BSE_METISSE) $(SSE_METISSE)

SRCE5 = \
bse.f evolv2.f comenv.f corerd.f dgcore.f gntage.f \
instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f zfuncs.f


OBJT5 = $(SRCE5:.f=.o)


%.o : $(SRC_DIR)/%.f
	$(FC) $(OMP) $(FFLAGS2) -c $<

%.o : $(SRC_DIR)/%.f90
	$(FC) $(OMP) $(FFLAGS1) -c $<

$(SSE_METISSE): $(OBJT4) star.o sse_support.o hrdiag.o remnant_support.o z_support.o track_support.o interp_support.o zcnsts.o metisse_mlwind.o deltat.o alloc_track.o dealloc_track.o assign_commons.o
	$(FC) $(OMP) -o $(BSE_DIR)/$(SSE_METISSE) $^


$(BSE_METISSE): $(OBJT5) star.o sse_support.o hrdiag.o remnant_support.o z_support.o track_support.o interp_support.o zcnsts.o metisse_mlwind.o deltat.o alloc_track.o dealloc_track.o assign_commons.o
	$(FC) $(OMP) -o $(BSE_DIR)/$(BSE_METISSE) $^


clean:
	rm -f $(BSE_DIR)/$(SSE_METISSE) $(BSE_DIR)/$(BSE_METISSE)
	rm -f *.mod *.o

