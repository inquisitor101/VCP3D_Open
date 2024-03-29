# ******************************************************************************
# *                 Makefile to create the executable.                         *
# ******************************************************************************

# ******************************************************************************
# * Define the compiler and linker and their flags.
# ******************************************************************************

#MAKE = make -j 16
MAKE = make -j 8
DEPENDENCIES = ./Dependencies

CXX = mpicxx -std=c++11 -fopenmp 
#CXX = g++ -std=c++11 -fopenmp

LINKER = $(CXX)

MPI_FLAGS  = -DHAVE_MPI
OMP_FLAGS  = -DHAVE_OPENMP
#PREC_FLAGS = -DUSE_SINGLE_PRECISION

CXX_OPTFLAGS = -msse4.2 -O2 -funroll-loops -ftree-vectorize -march=native

#CXX_DEBUG_FLAGS = -g -Wall -Wextra -Wno-unused-parameter -Wno-empty-body -DENABLE_NAN_CHECK

CXX_FLAGS = $(CXX_OPTFLAGS) $(CXX_DEBUG_FLAGS) $(MPI_FLAGS) $(OMP_FLAGS) \
	    $(PREC_FLAGS) -I../include

CXX_FLAGS_1 = $(CXX_OPTFLAGS) $(CXX_DEBUG_FLAGS) $(MPI_FLAGS) $(OMP_FLAGS) \
	      $(PREC_FLAGS) -I../include

CXX_FLAGS_2 = $(CXX_OPTFLAGS) $(CXX_DEBUG_FLAGS) $(MPI_FLAGS) \
	      $(PREC_FLAGS) -I../include


LAPACK_LINK_FLAGS = -llapack -lblas
LINKER_FLAGS      = $(CXX_OPTFLAGS) $(LAPACK_LINK_FLAGS) $(OMP_FLAGS)

# ******************************************************************************
# *                      Name of the executable                                *
# ******************************************************************************

EXEC = ../bin/VCP3D_CurviLinear

# ******************************************************************************
# * Rules to make the objects. Redefine .SUFFIXES to be sure all the desired   *
# * ones are included.                                                         *
# ******************************************************************************

.SUFFIXES: .o .f .F .f90 .F90 .c .C .cc .cpp .h .hh .H

../obj/%.o : %.cpp Makefile
	$(CXX) $(CXX_FLAGS) -c $< -o $@
	@echo
	@echo "        --- Compiled $*.cpp successfully ---"
	@echo

# ******************************************************************************
# *                      Names of the source files.                            *
# ******************************************************************************

CXX_SOURCES_1 = InputParamClass.cpp\
	        SolverClass.cpp\
	        StandAloneFunctions.cpp

CXX_SOURCES_2 = adt_structure.cpp\
	        BoundaryState.cpp\
					SpongeLayer.cpp\
	        CG_MassMatrix.cpp\
	        ElementClass.cpp\
	        ExchangeDataWallModelClass.cpp\
	        FluxesFace.cpp\
	        GaussJacobiQuadrature.cpp\
	        SGSModel.cpp\
	        StandardElementClass.cpp\
	        SubfaceClasses.cpp\
	        TensorProductResIFace.cpp\
	        TensorProductResIFaceK2.cpp\
	        TensorProductResIFaceK3.cpp\
	        TensorProductResIFaceK4.cpp\
	        TensorProductResIFaceK5.cpp\
	        TensorProductResIFaceK6.cpp\
	        TensorProductResIFaceK7.cpp\
	        TensorProductResIFaceK8.cpp\
	        TensorProductResIFaceK9.cpp\
	        TensorProductResIFaceK10.cpp\
	        TensorProductResJFace.cpp\
	        TensorProductResJFaceK2.cpp\
	        TensorProductResJFaceK3.cpp\
	        TensorProductResJFaceK4.cpp\
	        TensorProductResJFaceK5.cpp\
	        TensorProductResJFaceK6.cpp\
	        TensorProductResJFaceK7.cpp\
	        TensorProductResJFaceK8.cpp\
	        TensorProductResJFaceK9.cpp\
	        TensorProductResJFaceK10.cpp\
	        TensorProductResKFace.cpp\
	        TensorProductResKFaceK2.cpp\
	        TensorProductResKFaceK3.cpp\
	        TensorProductResKFaceK4.cpp\
	        TensorProductResKFaceK5.cpp\
	        TensorProductResKFaceK6.cpp\
	        TensorProductResKFaceK7.cpp\
	        TensorProductResKFaceK8.cpp\
	        TensorProductResKFaceK9.cpp\
	        TensorProductResKFaceK10.cpp\
	        TensorProductSolAndGradIFace.cpp\
	        TensorProductSolAndGradIFaceK2.cpp\
	        TensorProductSolAndGradIFaceK3.cpp\
	        TensorProductSolAndGradIFaceK4.cpp\
	        TensorProductSolAndGradIFaceK5.cpp\
	        TensorProductSolAndGradIFaceK6.cpp\
	        TensorProductSolAndGradIFaceK7.cpp\
	        TensorProductSolAndGradIFaceK8.cpp\
	        TensorProductSolAndGradIFaceK9.cpp\
	        TensorProductSolAndGradIFaceK10.cpp\
	        TensorProductSolAndGradJFace.cpp\
	        TensorProductSolAndGradJFaceK2.cpp\
	        TensorProductSolAndGradJFaceK3.cpp\
	        TensorProductSolAndGradJFaceK4.cpp\
	        TensorProductSolAndGradJFaceK5.cpp\
	        TensorProductSolAndGradJFaceK6.cpp\
	        TensorProductSolAndGradJFaceK7.cpp\
	        TensorProductSolAndGradJFaceK8.cpp\
	        TensorProductSolAndGradJFaceK9.cpp\
	        TensorProductSolAndGradJFaceK10.cpp\
	        TensorProductSolAndGradKFace.cpp\
	        TensorProductSolAndGradKFaceK2.cpp\
	        TensorProductSolAndGradKFaceK3.cpp\
	        TensorProductSolAndGradKFaceK4.cpp\
	        TensorProductSolAndGradKFaceK5.cpp\
	        TensorProductSolAndGradKFaceK6.cpp\
	        TensorProductSolAndGradKFaceK7.cpp\
	        TensorProductSolAndGradKFaceK8.cpp\
	        TensorProductSolAndGradKFaceK9.cpp\
	        TensorProductSolAndGradKFaceK10.cpp\
	        TensorProductSolAndGradVolume.cpp\
	        TensorProductSolAndGradVolumeK2.cpp\
	        TensorProductSolAndGradVolumeK3.cpp\
	        TensorProductSolAndGradVolumeK4.cpp\
	        TensorProductSolAndGradVolumeK5.cpp\
	        TensorProductSolAndGradVolumeK6.cpp\
	        TensorProductSolAndGradVolumeK7.cpp\
	        TensorProductSolAndGradVolumeK8.cpp\
	        TensorProductSolAndGradVolumeK9.cpp\
	        TensorProductSolAndGradVolumeK10.cpp\
	        TensorProductVolumeResidual.cpp\
	        TensorProductVolumeResidualK2.cpp\
	        TensorProductVolumeResidualK3.cpp\
	        TensorProductVolumeResidualK4.cpp\
          TensorProductVolumeResidualK5.cpp\
	        TensorProductVolumeResidualK6.cpp\
          TensorProductVolumeResidualK7.cpp\
	        TensorProductVolumeResidualK8.cpp\
          TensorProductVolumeResidualK9.cpp\
	        TensorProductVolumeResidualK10.cpp\
	        WallModel.cpp\
	        VCP3D_CurviLinear.cpp

# ******************************************************************************
# * Names of the object files to be generated.                                 *
# ******************************************************************************

CXX_OBJECTS_1 = $(patsubst %.cpp,../obj/%.o,$(CXX_SOURCES_1))
CXX_OBJECTS_2 = $(patsubst %.cpp,../obj/%.o,$(CXX_SOURCES_2))

CXX_OBJECTS := $(CXX_OBJECTS_1) $(CXX_OBJECTS_2)

# ******************************************************************************
# * Object dependant compiler flags.                                            *
# ******************************************************************************

$(CXX_OBJECTS_1): CXX_FLAGS := $(CXX_FLAGS_1)
$(CXX_OBJECTS_2): CXX_FLAGS := $(CXX_FLAGS_2)

# ******************************************************************************
# *                                  Targets.                                  *
# ******************************************************************************

default: all

all:
	$(MAKE) CXX_SRC
	$(MAKE) $(EXEC)

depend:
	touch $(DEPENDENCIES)
	makedepend -p../obj/ -o.o -Y -- $(CXX_FLAGS_1) -- $(CXX_SOURCES_1) \
	$(CXX_SOURCES_2) -f$(DEPENDENCIES)
	/bin/rm -f $(DEPENDENCIES).bak

clean:
	@echo "       Making clean in src"
	/bin/rm -f ../obj/* *~ ../include/*~ ../bin/VCP3D_CurviLinear

CXX_SRC: $(CXX_OBJECTS)

# ******************************************************************************
# *                      Rule to make the executable.                          *
# ******************************************************************************

$(EXEC):$(CXX_OBJECTS)
	  @echo
	  @echo "        creating  $(EXEC)..."
	  @echo
	  $(LINKER) -o $(EXEC) $(CXX_OBJECTS) $(LINKER_FLAGS)
	  @echo
	  @echo "       $(EXEC) created."
	  @echo

# ****************************************************************************
# *                                                                          *
# *                              Dependencies.                               *
# * These are generated by makedepend and must be executed once before the   *
# * actual make command in order to have the dependencies work correctly.    *
# *                                                                          *
# ****************************************************************************

-include ${DEPENDENCIES}
