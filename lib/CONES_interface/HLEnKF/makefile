
.PHONY: all

.DEFAULT: all

all: ../../../build/HLEnKF_main.exe

CC      =  mpicc
CCFLAGS = -std=c99 -C
CCFLAGS = -std=gnu99
C++     =  mpic++
C++FLAGS = -std=c++17 -m64 -Dlinux64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wno-unused-parameter -Wno-unused-but-set-variable -Wnon-virtual-dtor -O3 -DNoRepository -ftemplate-depth-100 -fPIC

CCSanitizer = -Wall -Werror -fsanitize=address -fno-omit-frame-pointer -g
CCSanitizer_address = -Wall -fsanitize=address -fno-omit-frame-pointer -g

C++EXTINC = -I$(FOAM_SRC)/MomentumTransportModels/momentumTransportModels/lnInclude -I$(FOAM_SRC)/MomentumTransportModels/incompressible/lnInclude -I$(FOAM_SRC)/transportModels/lnInclude -I$(FOAM_SRC)/finiteVolume/lnInclude -I$(FOAM_SRC)/fvOptions/lnInclude -I$(FOAM_SRC)/sampling/lnInclude -I$(FOAM_SRC)/dynamicFvMesh/lnInclude -I$(FOAM_SRC)/dynamicMesh/lnInclude -I$(FOAM_SRC)/meshTools/lnInclude -I$(FOAM_SRC)/OpenFOAM/lnInclude -I$(FOAM_SRC)/OSspecific/POSIX/lnInclude 
C++EXTINC += -I/usr/include/eigen3
C++EXTLIB = -lutil -lstdc++ -lmpi -lmpi_mpifh -lmpi_cxx -fuse-ld=bfd -Xlinker --add-needed -L$(FOAM_LIBBIN) -lOpenFOAM -lmomentumTransportModels -lincompressibleMomentumTransportModels -lfiniteVolume -lsampling -ldynamicFvMesh -ltopoChangerFvMesh -ldynamicMesh -lmeshTools -ldl -lm

# CWIPI_DIR  = ../cwipi-1.1.0/build
CWIPI_DIR  = /home/tom/Documents/cwipi-1.1.0/build
CWIPI_INC  = $(CWIPI_DIR)/include
CWIPI_LIB  = $(CWIPI_DIR)/lib
CWIPI_LINK =  -Wl,-rpath=$(CWIPI_LIB) -L$(CWIPI_LIB) -lcwp

DAFUNCDIR = src/DAFunctions
OBSDIR = src/DAClasses/observations
PROBEDIR = src/DAClasses/probe
REGIONDIR = src/DAClasses/region

region.o: $(REGIONDIR)/region.C
	$(C++) $(C++FLAGS) -c -o $@ $< -I$(CWIPI_INC) $(C++EXTINC) $(C++EXTLIB) $(CWIPI_LINK) -lm 

probe.o: $(PROBEDIR)/probe.C
	$(C++) $(C++FLAGS) -c -o $@ $< -I$(CWIPI_INC) $(C++EXTINC) $(C++EXTLIB) $(CWIPI_LINK) -lm 

observations.o: $(OBSDIR)/observations.C
	$(C++) $(C++FLAGS) -c -o $@ $< -I$(CWIPI_INC) $(C++EXTINC) $(C++EXTLIB) $(CWIPI_LINK) -lm 

HLEnKF_functions.o: $(DAFUNCDIR)/HLEnKF_functions.C
	$(C++) $(C++FLAGS) -c -o $@ $< -I$(CWIPI_INC) $(C++EXTINC) $(C++EXTLIB) $(CWIPI_LINK) -lm 

../../../build/HLEnKF_main.exe: HLEnKF_functions.o region.o probe.o observations.o HLEnKF_main.C | ../../../build
	$(C++) $(C++FLAGS) -o $@ $^ -I$(CWIPI_INC) $(C++EXTINC) $(C++EXTLIB) $(CWIPI_LINK)

../../../build:
	mkdir -p ../../../build

allclean:
	rm -f  ../../../build/HLEnKF_*
	rm -f *.o
