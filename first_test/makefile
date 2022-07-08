
.PHONY: all

.DEFAULT: all

all: KF_coupling.exe

CC      =  mpicc
CCFLAGS = -std=c99 -C
CCFLAGS = -std=gnu99
C++		=  mpic++
C++FLAGS = -m64 -Dlinux64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-unused-but-set-variable -Wnon-virtual-dtor -O3 -DNoRepository -ftemplate-depth-100 -fPIC

CCSanitizer = -Wall -Werror -fsanitize=address -fno-omit-frame-pointer -g
CCSanitizer_address = -Wall -fsanitize=address -fno-omit-frame-pointer -g

C++EXTINC = -I/home/villanul/OpenFOAM/villanul-8/run/cwipi_tests/first_test -I$(FOAM_SRC)/MomentumTransportModels/momentumTransportModels/lnInclude -I$(FOAM_SRC)/MomentumTransportModels/incompressible/lnInclude -I$(FOAM_SRC)/transportModels/lnInclude -I$(FOAM_SRC)/finiteVolume/lnInclude -I$(FOAM_SRC)/fvOptions/lnInclude -I$(FOAM_SRC)/sampling/lnInclude -I$(FOAM_SRC)/dynamicFvMesh/lnInclude -I$(FOAM_SRC)/dynamicMesh/lnInclude -I$(FOAM_SRC)/meshTools/lnInclude -I$(FOAM_SRC)/OpenFOAM/lnInclude -I$(FOAM_SRC)/OSspecific/POSIX/lnInclude -I/home/villanul/Documents/Travail/tools/eigen-3.4.0/
C++EXTLIB =  /usr/lib64/libstdc++.so.6 -lc -lgcc_s
C++EXTLIB = -lutil -lstdc++ -lgfortran -lmpi -lmpi_mpifh -lmpi_cxx -fuse-ld=bfd -Xlinker --add-needed -L$(FOAM_LIBBIN) -lOpenFOAM -lmomentumTransportModels -lincompressibleMomentumTransportModels -lincompressibleTransportModels -lfiniteVolume -lfvOptions -lsampling -ldynamicFvMesh -ltopoChangerFvMesh -ldynamicMesh -lmeshTools -ldl -lm

CWIPI_DIR  = /home/villanul/Documents/Travail/tools/cwipi/build
CWIPI_INC  = $(CWIPI_DIR)/include
CWIPI_LIB  = $(CWIPI_DIR)/lib
CWIPI_LINK =  -Wl,-rpath=$(CWIPI_LIB) -L$(CWIPI_LIB) -lcwpf -lcwp

KF_coupling.exe: KF_coupling_functions.o KF_coupling.C
	$(C++) $(C++FLAGS) -o KF_coupling.exe KF_coupling.C KF_coupling_functions.o -I$(CWIPI_INC) $(C++EXTINC) $(C++EXTLIB) $(CWIPI_LINK)

KF_coupling_functions.o:
	$(C++) $(C++FLAGS) -c KF_coupling_functions.C -I$(CWIPI_INC) $(C++EXTINC) $(C++EXTLIB) $(CWIPI_LINK) -lm 

clean:
	rm -f -r cwipi U*

allclean:
	rm -f -r KF_coupling.exe KF_coupling_functions.o cwipi U*
