CC      = mpicc
CCFLAGS = -std=c99
CCFLAGS = -std=gnu99

CCSanitizer = -Wall -Werror -fsanitize=address -fno-omit-frame-pointer -g
CCSanitizer_address = -Wall -fsanitize=address -fno-omit-frame-pointer -g

CWIPI_DIR  = /home/miguel/2022_CWIPI_Training/cwipi-0.12.0/build
CWIPI_INC  = $(CWIPI_DIR)/include
CWIPI_LIB  = $(CWIPI_DIR)/lib
CWIPI_LINK =  -Wl,-rpath,$(CWIPI_LIB) -L$(CWIPI_LIB) -lcwpf -lcwp

MPI_DIR = /usr/lib/x86_64-linux-gnu/openmpi
MPI_INC = $(MPI_DIR)/include

all: c_cavity_coupling c_cavity_coupling2 MyIcoFoam

c_cavity_coupling: c_cavity_coupling.o
	$(CC) $(CC_FLAGS) c_cavity_coupling.o -o c_cavity_coupling -I$(CWIPI_INC) -I$(MPI_INC) $(CWIPI_LINK) -lm 

c_cavity_coupling.o: c_cavity_coupling.c
	$(CC) $(CC_FLAGS) -c c_cavity_coupling.c -I$(CWIPI_INC) -I$(MPI_INC)

c_cavity_coupling2: c_cavity_coupling2.o
	$(CC) $(CC_FLAGS) c_cavity_coupling2.o -o c_cavity_coupling2 $(CCSanitizer_address) -I$(CWIPI_INC) -I$(MPI_INC) $(CWIPI_LINK) -lm 

c_cavity_coupling2.o: c_cavity_coupling2.c
	$(CC) $(CC_FLAGS) -c c_cavity_coupling2.c $(CCSanitizer_address) -I$(CWIPI_INC) -I$(MPI_INC)

MyIcoFoam: MyIcoFoam.o
	g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam8/src/finiteVolume/lnInclude -I/opt/openfoam8/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam8/src/OpenFOAM/lnInclude -I/opt/openfoam8/src/OSspecific/POSIX/lnInclude -fPIC -fuse-ld=bfd -Xlinker --add-needed -Xlinker --no-as-needed MyIcoFoam.o -L/opt/openfoam8/platforms/linux64GccDPInt32Opt/lib \
    -lfiniteVolume -lmeshTools -lOpenFOAM -ldl  \
     -I$(CWIPI_INC) -I$(MPI_INC) $(CWIPI_LINK) -lm -o MyIcoFoam

MyIcoFoam.o: MyIcoFoam.C createFields.H
	g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam8/src/finiteVolume/lnInclude -I/opt/openfoam8/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam8/src/OpenFOAM/lnInclude -I/opt/openfoam8/src/OSspecific/POSIX/lnInclude -fPIC -c MyIcoFoam.C -o MyIcoFoam.o -I$(CWIPI_INC) -I$(MPI_INC)

clean:
	rm -f c_cavity_coupling MyIcoFoam c_cavity_coupling2 *.o -r cwipi
