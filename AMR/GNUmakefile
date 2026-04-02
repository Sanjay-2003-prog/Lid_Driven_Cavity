# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= ../../ImportCodes/amrex

COMP    = gnu
DIM     = 2

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_HIP   = FALSE

# Enable exceptions and RTTI
CXXFLAGS_FROM_MAKE += -DAMREX_USE_EXCEPTION=1

# Include necessary packages
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Include your source files
include ./Make.package

# Include AMReX base components
include $(AMREX_HOME)/Src/Base/Make.package

# Include MLMG components (CRITICAL - this was missing)
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
